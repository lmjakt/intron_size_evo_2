require('entropy')
## a simple wrapper.. 
ensembl.connect <- function(db.name){
    ## dbConnect(MySQL(), user='anonymous', dbname=db.name, port=5306,
    ##           host='ensembldb.ensembl.org', password='')
    dbConnect(MySQL(), user='lmj', dbname=db.name,
                    host='localhost')
}

## Gets exons from protein coding genes only
## which is what we will use because these will be the most well annotated.
get.ens.exon.info <- function(db.name){
    query <-paste("select a.gene_id, b.name as 'chr', a.seq_region_strand as 'strand', e.transcript_id as 'transcript',",
                  "d.exon_id, d.seq_region_start as 'start', d.seq_region_end as 'end', c.rank",
                  "from gene a",
                  "inner join seq_region b on a.seq_region_id=b.seq_region_id",
                  "inner join transcript e on e.gene_id=a.gene_id",
                  "inner join exon_transcript c on e.transcript_id=c.transcript_id",
                  "inner join exon d on c.exon_id=d.exon_id",
                  "where a.biotype='protein_coding'",
                  "order by a.gene_id, e.transcript_id, c.rank;")
    db <- ensembl.connect(db.name)
    exons <- dbGetQuery(db, query);
    exons <- exons[ !grepl("ALT|MT", exons[ ,'chr']), ]
    dbDisconnect(db)
    exons
}

get.genome.stats <- function(db.name){
    query <- "select * from genome_statistics;"
    db <- ensembl.connect(db.name)
    stats <- dbGetQuery(db, query)
    dbDisconnect(db)
    stats
}

get.ens.canonical.translations <- function(db.name){
    query <- paste("select a.stable_id as 'gene', b.stable_id as 'transcript', c.stable_id as 'translation'",
                   "from gene a",
                   "inner join transcript b on a.canonical_transcript_id=b.transcript_id",
                   "inner join translation c on c.transcript_id=b.transcript_id",
                   "where a.biotype='protein_coding';");
    db <- ensembl.connect(db.name)
    data <- dbGetQuery(db, query);
    dbDisconnect(db)
    data
}

get.ensembl.taxonomy <- function(db.name='ensembl_compara_98', root.name='Ensembl'){
    query <- paste("select a.*, c.rank from species_tree_node a",
                   sprintf("inner join species_tree_root b on a.root_id=b.root_id and b.label='%s'", root.name),
                   "inner join ncbi_taxa_node c on a.taxon_id=c.taxon_id",
                   "order by a.left_index;")
    db <- ensembl.connect(db.name)
    data <- dbGetQuery(db, query)
    dbDisconnect(db)
    rownames(data) <- data[,'node_id']
    data
}

## tree is a table as obtained by get.ensembl.taxonomy
## sp is species name
get.ensembl.lineage <- function( tree, sp ){
    i <- which(tree$node_name == sp)
    if(length(i) == 0)
        stop(paste(sp, "not found"))
    if(length(i) > 1)
        warning(paste("Found", length(i), "entries for", sp))
    i <- i[1]  ## we have no other real information as to which one to choose
    id <- as.character( tree[i, 'node_id'] )
    lineage <- tree[id, , drop=FALSE]
    while( !is.na( tree[id, 'parent_id'] )){
        id <- as.character( tree[ as.character( tree[ id, 'parent_id' ] ), 'node_id' ] )
        lineage <- rbind(lineage, tree[id, ])
    }
    lineage
}

get.ensembl.dist <- function(tree, sp1, sp2){
    lin1 <- get.ensembl.lineage(tree, sp1)
    lin2 <- get.ensembl.lineage(tree, sp2)
    d1 <- sum( lin1[ !(lin1$node_name %in% lin2$node_name), 'distance_to_parent' ] )
    d2 <- sum( lin2[ !(lin2$node_name %in% lin1$node_name), 'distance_to_parent' ] )
    c(d1, d2, d1 + d2)
}

## returns the intron length, rank and inverse rank
extract.intron.lengths <- function(df, id.col='transcript'){
    get.i <- function(i){
        if(is.unsorted( df[i,'rank'])){
            print(df[i,])
            stop("Exons not ordered")
        }
        if(length(i) == 0)
            return(c())
        if(df[i[1], 'strand'] == 1){
            j <- i[1:length(i)]
            rank <- df[i[-length(i)], 'rank']
        }else{
            j <- i[length(i):1]
            rank <- df[ j[-1], 'rank' ]
        }
        sz <- df[ j[-1], 'start'] - df[ j[-length(j)], 'end']
        cbind(sz, rank, 1 + length(rank) - rank)
    }
    tmp <- tapply( 1:nrow(df), df[,id.col], get.i)
    cbind( 'size'=unlist(sapply(tmp, function(x){ x[,1] })),
          'rank'=unlist(sapply(tmp, function(x){ x[,2] })),
          'irank'=unlist(sapply(tmp, function(x){ x[,3] })) )
}

## get family members from compara
## this is a bit wasteful as it also gets the family descriptions and so on. But I don't
## think we are running short of memory at the moment.
get.family.members <- function(sp.name, db, db.is.name=FALSE){
    family.query.fmt <-
     "select a.stable_id as 'gene.id', a.description, b.stable_id as 'protein.id',
     d.stable_id as 'family.id', d.description as 'family', e.name as 'species'
     from gene_member a
     inner join seq_member b on a.gene_member_id=b.gene_member_id
     inner join family_member c on b.seq_member_id=c.seq_member_id
     inner join family d on c.family_id=d.family_id
     inner join genome_db e on a.genome_db_id=e.genome_db_id
     and e.name='%s';"
    query <- sprintf(family.query.fmt, sp.name)
    if(db.is.name)
        db <- ensembl.connect(db)
    data <- dbGetQuery(db, query)
    if(db.is.name)
        dbDisconnect(db)
    data
}

## get genes from the reference assembly only
## This uses 16 to indicate non-ref. That is a bad thing, I should really make a join to attrib_type and
## specify that we are lookin gfor 'Non Ref'
get.all.ref.genes <- function(db.name, max.coord.rank=1){
    query.1 <- paste("select a.stable_id, b.name as 'chr', a.seq_region_start as 'start', a.seq_region_end as 'end', a.seq_region_strand as 'strand',",
                     "a.biotype",
                     "from gene a",
                     "inner join seq_region b on a.seq_region_id=b.seq_region_id",
                     "inner join coord_system c on b.coord_system_id=c.coord_system_id and c.rank <= %d")
    query.2 <- paste("(select a.stable_id",
                     "from gene a",
                     "inner join seq_region b on a.seq_region_id=b.seq_region_id",
                     "inner join coord_system c on b.coord_system_id=c.coord_system_id and c.rank <= %d",
                     "inner join seq_region_attrib d on a.seq_region_id=d.seq_region_id and d.attrib_type_id=16 and d.value=1)")
    query.1 <- sprintf(query.1, max.coord.rank)
    query.2 <- sprintf(query.2, max.coord.rank)
    query <- paste(query.1, "WHERE a.stable_id NOT IN", query.2, ";")
    db <- ensembl.connect(db.name)
    data <- dbGetQuery(db, query)
    dbDisconnect(db)
    data
}


## h.all = intron.lengths.2.h, h.1 -intron.lengths.2.
plot.intron.hist <- function(h.all, h.1, h.2, oma=c(2.1, 4.1, 4.1, 4.1), mar=c(0.1, 0.1, 0.1, 0.1), lab.cex=1.5, lwds=c(2,1,1),
                             col.labels=c('species', 'density', 'counts', 'all', 'rank=1', 'rank>1'),
                             l.cols=c('white', rgb(0.8, 0, 0), rgb(0, 0, 0.9)), pl.bg=rgb(0.8, 0.8, 0.8),
                             v.lines=NULL, v.lin.col=rgb(1,1,0.8), v.lin.lwd=1){
    par('oma'=oma)
    par('mar'=mar)
    layout(matrix( 1:((2 + length(h.all)) * 6), ncol=6, byrow=TRUE ),
           widths=c(0.5, 1, 1, 0.2, 0.2, 0.2) )
    for(i in 1:6){
        plot.new()
        plot.window(xlim=c(0,1), ylim=c(0,1))
        text( 0.5, 0.5, col.labels[i], cex=lab.cex, font=2 )
    }
    plot.row <- function(sp, hists, cols=l.cols, bg.col=pl.bg){
        max.cnt <- max(unlist(sapply(hists, function(x){ x$counts })))
        max.dens <- max(unlist(sapply(hists, function(x){ x$density })))
        plot.new()
        plot.window( xlim=c(0,1), ylim=c(0,1))
        text(0.5, 0.5, sp, cex=lab.cex)
        plot( hists[[1]]$mids, hists[[1]]$density, ylim=c(0, max.dens), axes=FALSE, col=cols[1], type='n', xaxs='i' )
        with(par(), rect(usr[1], usr[3], usr[2], usr[4], col=bg.col, border=NA ))
        if(!is.null(v.lines))
            abline(v=v.lines, lwd=v.lin.lwd, lty=2, col=v.lin.col)
        points( hists[[1]]$mids, hists[[1]]$density, col=cols[1], type='l', lwd=lwds[1])
        points( hists[[1]]$mids, hists[[2]]$density, col=cols[2], type='l', lwd=lwds[2])
        points( hists[[1]]$mids, hists[[3]]$density, col=cols[3], type='l', lwd=lwds[3] )
        plot( hists[[1]]$mids, hists[[1]]$counts, ylim=c(0, max.cnt), axes=FALSE, col=cols[1], type='n', xaxs='i' )
        with(par(), rect(usr[1], usr[3], usr[2], usr[4], col=bg.col, border=NA ))
        if(!is.null(v.lines))
            abline(v=v.lines, lwd=v.lin.lwd, lty=2, col=v.lin.col)
        points( hists[[1]]$mids, hists[[1]]$counts, col=cols[1], type='l', lwd=lwds[1] )
        points( hists[[1]]$mids, hists[[2]]$counts, col=cols[2], type='l', lwd=lwds[2] )
        points( hists[[1]]$mids, hists[[3]]$counts, col=cols[3], type='l', lwd=lwds[3] )
        for(i in 1:length(hists)){
            plot.new()
            plot.window(xlim=c(0,1), ylim=c(0,1))
            text(0.5, 0.5, sum( hists[[i]]$counts ), cex=lab.cex)
        }
    }
    for(i in 1:length(h.all)){
        plot.row( names(h.all)[i], list(h.all[[i]], h.1[[i]], h.2[[i]]) )
    }
    plot.new()
    plot.new()
    plot.new()
    plot.window(xlim=c(0,10), ylim=c(0,1))
    with(par(), rect(usr[1], usr[3], usr[2], usr[4], col=pl.bg, border=NA ))
    segments( seq(2, 8, 3), 0.5, seq(3, 9, 3)-0.3, 0.5, lwd=lwds, col=l.cols )
    text( seq(2, 8, 3)-0.2, 0.5, c('all', 'rank=1', 'rank>1'), adj=c(1, 0.5))
}

mutual.info <- function(x1, x2, tform=log2, numBins=20){
    x1 <- tform(x1)
    x2 <- tform(x2)
    b <- !(is.na(x1) | is.na(x2)) & (is.finite(x1) & is.finite(x2))
    x1 <- x1[b]
    x2 <- x2[b]
    x.2d <- discretize2d( x1, x2, numBins1=numBins, numBins2=numBins )
    x.2d.r <- discretize2d( x1, sample(x2), numBins1=numBins, numBins2=numBins )
    rs <- rowSums( x.2d )
    cs <- colSums( x.2d )
    list('m'=x.2d, 'm.r'=x.2d.r, 'H1'=entropy(rs), 'H2'=entropy(cs), 'H12'=entropy(x.2d),
         'rs'=rs, 'cs'=cs, 'mi'=mi.empirical(x.2d), 'mi.r'=mi.empirical(x.2d.r), 'n'=sum(b))
}

intron.med.mi.sp <- function( int.s=intron.s, gene.st=gene.stats, sp1, sp2, numBins=20, tform=log2){
    med.l <- tapply( 1:length(int.s), gene.st$family, function(i){
        sp.dup <- duplicated( gene.st$sp[i] )
        sp.b1 <- gene.st$sp[i] == sp1 & !sp.dup
        sp.b2 <- gene.st$sp[i] == sp2 & !sp.dup
        if(sum(sp.b1) == 0 || sum(sp.b2) == 0)
            return(c())
        return( c(median( int.s[[ i[sp.b1] ]]),
                  median( int.s[[ i[sp.b2] ]]) ) )
    })
    b <- sapply( med.l, length) == 2
    med <- sapply( med.l[b], eval )
    ## NAs can occur where we have single exon genes..
    med <- med[ ,!is.na( med[1,] ) & !is.na(med[2,]) ]
    med <- tform(med)
    med.2d <- discretize2d( med[1,], med[2,], numBins1=numBins, numBins2=numBins )
    med.2d.r <- discretize2d( med[1,], med[2, sample(1:ncol(med)) ], numBins1=numBins, numBins2=numBins )
    med.2d.i <- discretize2d( med[1,], med[1,], numBins1=numBins, numBins2=numBins )
    rs <- rowSums( med.2d )
    cs <- colSums( med.2d )
    list('m'=med.2d, 'H1'=entropy(rs), 'H2'=entropy(cs), 'H12'=entropy(med.2d), 'rs'=rs, 'cs'=cs,
         'mi'=mi.empirical(med.2d), 'mi.r'=mi.empirical(med.2d.r), 'mi.i'=mi.empirical(med.2d.i))
}

plot.med.mi <- function(x, val.cm=0.3, val.v=0.8, tform=eval, ps.count=0, lab.cex=1){
    width <- ncol( x$m )
    height <- nrow(x$m) ## these should be the same

    plot.new()
    ## let us plot both the observed and the expected.. 
    plot.window(xlim=c(0, 4 + 2 * width + 3), ylim=c(-0.5, 2 + height ))
    x2 <- matrix(1:width, ncol=width, nrow=height, byrow=TRUE )
    x1 <- x2 - 1
    y2 <- 2 + matrix(1:height, ncol=width, nrow=height )
    y1 <- y2 - 1
    m.exp <- tform( ps.count + (matrix(x$rs, ncol=1) %*% matrix(x$cs, nrow=1)) / sum(x$m) )
    m.obs <- tform( ps.count + x$m)
    rs <- tform( ps.count + x$rs )
    cs <- tform( ps.count + x$cs )
    cols <- hsvScale( c(m.obs, m.exp), val=val.v * (val.cm + c(m.obs, m.exp)) / (val.cm + max(c(m.obs, m.exp))))
    rect( x1 + 2, y1, x2 + 2, y2, col=cols[1:length(x1)],
         border='grey', lwd=0.25 )
    rect( x1 + 4 + width, y1, x2 + 4 + width, y2, col=cols[(1+length(x1)):length(cols)],
         border='grey', lwd=0.25 )
    rect( x1[1,] + 2, 0.5, x2[1,] + 2, 1.5, col=hsvScale( cs, val=val.v * (val.cm + cs) / (val.cm + max(cs))),
         border='grey', lwd=0.25)
    rect( x1[1,] + 4 + width, 0.5, x2[1,] + 4 + width, 1.5, col=hsvScale( cs, val=val.v * (val.cm + cs) / (val.cm + max(cs))),
         border='grey', lwd=0.25)
    rect( 0.5, y1[,1], 1.5, y2[,1], col=hsvScale( rs, val=val.v * (val.cm + rs) / (val.cm + max(rs))),
         border='grey', lwd=0.25)
    rect( width + 2.5, y1[,1], width + 3.5, y2[,1], col=hsvScale( rs, val=val.v * (val.cm + rs) / (val.cm + max(rs))),
         border='grey', lwd=0.25)

    text( mean(x1) + 2, -0.5, sub("_", " ",x$sp2), cex=lab.cex )
    text( mean(x1) + width + 4, -0.5, sub("_", " ",x$sp2), cex=lab.cex )
    text( 0, mean(y1), sub("_", " ", x$sp1), srt=90, cex=lab.cex)
    x3 <- 4 + 2 * width + 2
    y <- seq(2, max(y1), length.out=100 )
    yd <- diff(y)[1]
    rect( x3, y, x3+1, y + yd, col=hsvScale( y, val=val.v * (val.cm + y) / (val.cm + max(y))), border=NA )
    arrows(x3-0.5, 2, x3-0.5, max(y), lwd=2, length=0.15 )
}

### for reading in introns
### horrendously inefficient. Absolutely horrid, in order to avoid looping
### through lines..
read.transcripts <- function(fname){
    exons.l <- readLines(fname)
    id.l <- grep(">", exons.l)
    s.e <- c(id.l - 1, length(exons.l))[-1]
    exons <- lapply( 1:length(id.l), function(i){
        exons.l[(id.l[i]+1):s.e[i]]
    })
    names(exons) <- sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][1] })
    exon.dbs <- unname(sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][4] }))

    lapply( 1:length(exons), function(i){
        list('id'=sub("^>", "", names(exons)[i]), 'db'=exon.dbs[i], 'sp'=sub("([^_]+)_([^_]+)_.+$", "\\1 \\2", exon.dbs[i]),
             'l'=unname(sapply(exons[[i]], nchar)), 's'=paste(exons[[i]], collapse='') )
    })
}

## f1.2 => probability of 1 passing given that 2 has passed
## f2.1 => probability of 2 passing givent that 1 has passed
pairwise.phyper <- function(x1, x2, min.v=NULL, cond1=function(a){ a >= min.v }, cond2=cond1){
    b <- !(is.na(x1) | is.na(x2))
    x1 <- x1[b]
    x2 <- x2[b]
    b1 <- cond1(x1)
    b2 <- cond2(x2)
    p <- phyper( sum(b1 & b2), sum(b1), sum(!b1), sum(b2), lower.tail=FALSE )
    list('n'=sum(b), 'f1'=sum(b1)/sum(b), 'f2'=sum(b2)/sum(b), 'f1.2'=sum(b1&b2) / sum(b2),
         'f2.1'=sum(b1&b2) / sum(b1), 'p'=p)
}

## non-parametric; how many observed
## ranks sorted by x1 ranks and number of x2
## within x1 rank returned for each rank.
rank.pairwise.sum <- function(x1, x2){
    b <- !(is.na(x1) | is.na(x2))
    x1 <- x1[b]
    x2 <- x2[b]
    r1 <- rank(x1, ties.method='max')
    r2 <- rank(x2, ties.method='max')
    rr1 <- rank(-x1, ties.method='max')
    rr2 <- rank(-x2, ties.method='max')
    o1 <- order(r1)
    o1r <- order(rr1)
    ## a loop works much better
    m <- matrix(0, nrow=length(r1), ncol=2)
    c1 <- 0
    c2 <- 0
    l <- length(x1)
    for(i in 1:length(r2)){
        m[ i, 1] <- sum( r2[ o1[1:i] ] <= r1[ o1[i] ] )
        m[ i, 2] <- sum( rr2[ o1r[1:i] ] <= rr1[ o1r[i] ] )
    }
    m
}

N50 <- function(regions, size, p=0.5){
    ln <- sort(regions[,'length'], decreasing=TRUE)
    b.i <- which(cumsum(ln) >= (size * p))
    if(!length(b.i))
        return(0)
    return(ln[ b.i[1] ])
}
