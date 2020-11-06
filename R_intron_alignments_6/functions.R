var.par <- function(x){
    c('n'=sum(!is.na(x)), 'mean'=mean(x, na.rm=TRUE), 'sd'=sd(x, na.rm=TRUE),
      'min'=min(x, na.rm=TRUE), 'max'=max(x, na.rm=TRUE), quantile(x, probs=seq(0,1,0.1), na.rm=TRUE))
}

sp.name <- function(x, from=".", to="_", toCap=TRUE, toLow=FALSE){
    x <- sub(from, to, x, fixed=TRUE)
    if(toCap)
        substr(x, 1,1) <- toupper( substr(x, 1,1))
    if(toLow)
        substr(x, 1,1) <- tolower( substr(x, 1,1))
    x
}

to.sp <- function(x){
    x <- sub("[._]", " ", x)
    substr(x,1,1) <- toupper( substr(x,1,1) )
    x
}

read.ortho <- function(path="../R_172_genomes/", conv.sp.names=FALSE){
    files <- paste(path, c('dr_intron_fam.txt', 'dr_intron_orthology_id.txt',
                           'dr_intron_orthology_tr.txt', 'dr_intron_orthology_i.txt',
                           'dr_intron_orthology_l.txt'), sep="")
    orth <- lapply(files, read.table, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    if(conv.sp.names)
        for(i in 2:length(orth))
            colnames(orth[[i]]) <- sp.name( colnames(orth[[i]]) )
    
    names(orth) <- c('fam', 'id', 'tr', 'i', 'l')
    orth
}


## this makes use of a number 
## of global variables; use with care.. 
## max.l is the maximum length of the score tables
## i is the row int fam.id / orth.i / from which
## we extract b.seq
## if j == i then we align b.seq to the orthologous
## intron is al.sp; otherwise as a control we align to
## a control intron chosen to have a simlar size.
align.introns <- function(i, j=i, al.sp='danio.rerio',
                          gap=c(-10,-1), min.width=15, min.score=40,
                          sm=sub.matrix, max.l=max.l, r.c=FALSE,
                          orth.fam=orth$fam, ortho=orth){
    fam.id <- orth.fam[i, 1]
    tr <- ortho$tr[i,]
    int.i <- ortho$i[i,]
    
    int.seq <- read.exons( intron.f[ fam.id ] )
    int.meta <- sapply( int.seq, function(x){ c('id'=x$id, 'tr'=x$tr, 'sp'=x$sp, 'n'=length(x$l) )})
    ## select the sequences that we will use
    b <- int.meta['tr',] %in% tr
    int.seq <- int.seq[b]
    int.meta <- int.meta[,b]
    ## so we can order by species
    names(int.seq) <- sub(" ", ".", int.meta['sp',])
    colnames(int.meta) <- sub(" ", ".", int.meta['sp',])
    
    int.i <- int.i[ names(int.i) %in% names(int.seq) ]
    al.i <- which( names(int.i) == al.sp )
    if(i == j){
        a.seq <- int.seq[[al.sp]]$e[ as.numeric(int.i[al.i]) ]
    }else{
        a.int.seq <- read.exons( intron.f[ orth.fam[j,1] ] )
        a.int.meta <- sapply( a.int.seq, function(x){ c('id'=x$id, 'tr'=x$tr,
                                                        'sp'=x$sp, 'n'=length(x$l) )})
        a.seq.b <- a.int.meta['tr',] == ortho$tr[j,al.sp]
        a.seq <- a.int.seq[[ which(a.seq.b) ]]$e[ ortho$i[j,al.sp] ]
    }
    ##        
    b.seq <- sapply( names(int.i[-al.i]), function(sp){
        int.seq[[sp]]$e[ as.numeric(int.i[sp]) ] })
    b.seq <- b.seq[ !is.na(b.seq) ]
    b.sp <- names(b.seq)
    ## reverse complement if r.c
    if(r.c){
        b.seq <- rev.comp(b.seq)
        names(b.seq) <- b.sp
    }
    ##
    aligns <- vector(mode='list', length=length(b.seq))
    names(aligns) <- names(b.seq)
    for(sp in names(b.seq)){
        if(as.double(nchar(a.seq)) * as.double(nchar(b.seq[sp])) <= max.l )
            aligns[[sp]] <- local.aligns( a.seq, b.seq[sp], sm$offset,
                                        sm$size, sm$sm,
                                        gap=gap, min.width=min.width, min.score=min.score )
    }
    list(sp1=al.sp, sp=names(b.seq), al=aligns, l1=nchar(a.seq), l2=unname(sapply(b.seq, nchar)),
         i=i, j=j, al.i=al.i, seq.meta=int.meta )
}


align.introns.ctl <- function( i, max.l, gap, ... ){
    ## find a suitable j to make the alignment from
    rank <- ga.int.l.r[ i ]
    rank.r <- -1:1 + rank
    rr.i <- ga.int.l.or[ rank.r ]
    rr.i <- rr.i[!is.na(rr.i) & rr.i != i]
    align.introns(i=i, j=rr.i[1], max.l=max.l, gap=gap, ...)
}

### lets see if we can plot some of these alignments and what we can do
### with them.
alignment.figure <- function(al, sp.col=sp.col, class.col=class.col,
                             col.tform=log2, col.f=score.col,
                             id.tbl=orth$id, tr.tbl=orth$tr,
                             i.tbl=orth$i, l.tbl=orth$l){
    id <- id.tbl[ al$i, al$sp1 ]
    tr <- tr.tbl[ al$i, al$sp1 ]
    int.i <- i.tbl[ al$i, al$sp1 ]
    int.l <- l.tbl[ al$i, al$sp1 ]
    layout(matrix(2:1, nrow=1), widths=c(0.075, 0.95))
    par(mar=c(5.1, 0.1, 4.1, 2.1))
    tmp <- plot.alignments(al, col.f=col.f, col.tform=col.tform, sp.col=sp.col)
    ##
    mtext( sprintf("%s : %s intron no: %d length: %d", id, tr, int.i, int.l), cex=1.2 )
    if(nrow(tmp$coord) == 0){
        with(par(), text(usr[1], usr[4], sprintf("%d alignments above max size (%1.1e)",
                                                 sum( al$l1 * al$l2 > max.l ), max.l ),
                         adj=c(0, 1)))
        return()
    }
    par(mar=c(5.1, 1.1, 4.1, 0.1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,100))
    text(0.05, 100 - 1.5 * strheight("A") * 1:length(class.col),
         capitalise(names(class.col)), col=class.col, adj=c(0,0), font=2)
    hm.v <- seq( tmp$cols$range[1], tmp$cols$range[2], length.out=100 )
    hm.y <- seq( 0, 25, length.out=length(hm.v))
    hm.yd <- diff(hm.y)[1]
    rect( 0.85, hm.y, 1.0, hm.y + hm.yd, col=col.f(hm.v)$cols, border=NA )
    lab.i <- seq(1, length(hm.v), length.out=5)
    text( 0.82, hm.y[ lab.i ] + hm.yd/2, sprintf("%.1e", hm.v[lab.i]), adj=c(1,0.5),
         cex=0.85)
}
