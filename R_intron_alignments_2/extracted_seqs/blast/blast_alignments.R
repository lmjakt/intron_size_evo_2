source('../../../R_common_functions/functions.R') ## for ensembl functions
require("RMySQL")

## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

sp.name <- function(x){
    x <- sub(".", "_", x, fixed=TRUE)
    substr(x, 1,1) <- toupper( substr(x, 1,1))
    x
}

## unfortunately we will need the full set of alignments as well
## read in the alignment data and assign to a named list.

aligns <- list(
    'long' = readRDS("../../../R_intron_alignments/al_top_500.rds"),
    'long.ctl' = readRDS("../../../R_intron_alignments/al_top_500_ctl.rds"),
    'sampled' = readRDS("../../../R_intron_alignments/al_ctl_500.rds"),
    'sampled.ctl' = readRDS("../../../R_intron_alignments/al_ctl_500_ctl.rds")
)

aligns.top.cl <- readRDS("../../aligns_top_al_cl.rds")
## I have also created an extened set of control alignments. I will need to read these into
## this analysis.
aligns.ctl.top.cl <- readRDS("../../../R_intron_alignments_5/ctl_aligns_top_cl.rds")
aligns.ctl <- readRDS("../../../R_intron_alignments_5/ctl_aligns.rds")

## set each on to use the rownumber (i) as a name
for(i in 1:length(aligns))
    names( aligns[[i]] ) <- as.character( sapply( aligns[[i]], function(x){ x$i }) )

names(aligns.ctl) <- as.character( sapply( aligns.ctl, function(x){ x$i }) )

## first lets get a set of Ensembl databases from the first line
## of an ouput file.
tmp <- readLines( "../../../family_members/vertebrate_family_members_1_130.txt", n=1)
dbs <- strsplit( tmp, '\t' )[[1]]
rm(tmp)
## lets set the names of the dbs

tmp <- sub("^([^_]+_[^_]+)_.+", "\\1", dbs)
substr(tmp, 1, 1) <- toupper(substr(tmp, 1, 1))
names(dbs) <- tmp

sp.class <- readRDS( "../../../R_trees_distances/sp_class_3.rds")
class.col <-  readRDS( "../../../R_trees_distances/class_col_3.rds")
substr( rownames(sp.class), 1, 1) <- toupper( substr( rownames(sp.class), 1, 1))

bl.files <- c( list.files("./", pattern="\\.bl$"),
               list.files("../../../R_intron_alignments_5/extracted_seqs/blast", pattern="\\.bl$", full.names=TRUE) )

## since the second set have directory names in them we need to remove these
## to parse the class and samples
bl.files.bn <- sub(".*?([^/]+)$", "\\1",  bl.files)
bl.files.cl <- sub("([^_]+)_.+", "\\1", bl.files.bn)
bl.files.sm <- sub("([^_]+)_([^_]+)_.+", "\\2", bl.files.bn)
bl.files.sp <- sub(".+?_([A-Z].+?)\\.bl$", "\\1", bl.files.bn) ## with capitalized genus
names(bl.files) <- paste( bl.files.cl, bl.files.sm, bl.files.sp, sep="_")

bl.columns <- strsplit('qseqid sseqid qlen qstart qend sstart send evalue bitscore score length pident nident qcovs qcovhsp', ' ')[[1]]

## this is incredibly wasteful, but whatever.
bl.data <- lapply( bl.files, function(x){
    data <- read.table( x, sep="\t", quote="", stringsAsFactors=FALSE)
    colnames(data) <- bl.columns
    data
})
names(bl.data) <- names(bl.files)

## and the extracted aligned sequences which were used as blast query sequences
fa.files <- c(list.files("../", pattern="fasta$", full.names=TRUE),
              list.files("../../../R_intron_alignments_5/extracted_seqs", pattern="fasta$", full.names=TRUE ))
fa.files.bn <- sub(".*?([^/]+)$", "\\1", fa.files )
fa.files.cl <- sub("([^_]+)_.+", "\\1", fa.files.bn)
fa.files.sm <- sub("([^_]+)_(.+?)\\.fasta", "\\2", fa.files.bn)
names(fa.files) <- paste(fa.files.cl, fa.files.sm, sep="_")

## the fasta sequences are on a single line here, so we can read in the full set
## using:

read.fa <- function( fn ){
    tmp <- readLines(fn)
    des <- tmp[ seq(1, length(tmp), 2) ]
    ids <- unname(sapply(des, function(x){ sub("^>([^ ]+).+", "\\1", x) }))
    introns <- as.numeric(sapply( strsplit(ids, '_'), function(x){ x[2] }))
    genes <- unname(sapply( strsplit(ids, '_'), function(x){ x[1] }))
    des.terms <- strsplit( des, " " )
    tr <- unname( sapply( des.terms, function(x){ x[2] }))
    sp1 <- unname( sapply( des.terms, function(x){ x[5] }))
    sp2 <- unname( sapply( des.terms, function(x){ x[6] }))
    coords <- t(sapply( des.terms, function(x){ as.numeric(x[7:12]) }))
    colnames(coords) <- c('a_beg', 'a_end', 'b_beg', 'b_end', 'score', 'length')
    seq <- tmp[ seq(2, length(tmp), 2) ]
    list(id=ids, gene=genes, intron=introns, des=des,
         seq=seq, tr=tr, sp1=sp1, sp2=sp2, coords=coords)
}

query.seq <- lapply( fa.files, read.fa )

#### Then we can decide on how to handle this.
#### Since we have a record of the species with the best alignment we can
#### look for that file.

conv.sp.names <- function(x){
    x <- sub(".", "_", x, fixed=TRUE)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

## For more advanced plotting we need to know which intron was looked at in each case.
## In order to know which family and gene we need to include the following information
orth.fam <- read.table("../../../R_172_genomes/dr_intron_fam.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)[1]
orth.id <-  read.table("../../../R_172_genomes/dr_intron_orthology_id.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.tr <-  read.table("../../../R_172_genomes/dr_intron_orthology_tr.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.i <-  read.table("../../../R_172_genomes/dr_intron_orthology_i.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.l <-  read.table("../../../R_172_genomes/dr_intron_orthology_l.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)

colnames(orth.id) <- conv.sp.names(colnames(orth.id))
colnames(orth.tr) <- colnames(orth.id)
colnames(orth.i) <- colnames(orth.id)
colnames(orth.l) <- colnames(orth.id)

log.tf <- function(x){ log2(x + 1) }

db.cred <- read.credentials( 'db.cred' )
## modify password part here.. and remove before saving
## This would seem to be a situation where using an external pointer
## would be the way to go as that doesn't get remembered.
## ex.height refers to the half height
plot.transcript.region <- function(orth.row, species, seed.coords, bl, query.id, alignments, max.nm=1e8,
                                   ex.height=0.4, bl.height=0.3){
    tr <- orth.tr[ orth.row, species ]
    int.i <- orth.i[ orth.row, species ]
    db <- dbs[ species ]
    tr.exons <- get.transcript.exons( db.cred, db, tr )
    tr.description <- get.transcript.description( db.cred, db, tr )
    if(!nrow(tr.exons)){
        plot.new()
        plot.window(xlim=c(0,1), ylim=c(0,1))
        mtext(species)
        plot.new()
        return()
    }
    strand = tr.exons[1,'strand']
    int.i.beg <- ifelse(strand == 1, tr.exons[ int.i, 'ex_end'], -tr.exons[ int.i, 'ex_start'] )
    int.i.end <- ifelse(strand == 1, tr.exons[ int.i+1, 'ex_start'], -tr.exons[ int.i+1, 'ex_end'] )
    xlim <- range( strand * tr.exons[,c('ex_start', 'ex_end')] )
    plot.new()
    plot.window(xlim=xlim, ylim=c(0,4))
    lines( xlim, c(1,1) )
    rect( strand * tr.exons[,'ex_start'], 1-ex.height, strand * tr.exons[,'ex_end'], 1 + ex.height )
    segments( int.i.beg, 1, int.i.end, 1, lwd=4 )
    ## if(strand == 1)
    ##     segments( strand * tr.exons[int.i,'ex_end'], 1, strand * tr.exons[int.i + 1, 'ex_start'], 1, lwd=4 )
    ## else
    ##     segments( strand * tr.exons[int.i,'ex_start'], 1, strand * tr.exons[int.i + 1, 'ex_end'], 1, lwd=4 )
##
    ## lets draw the blast matches that fit within this directory
    bl <- bl[ bl$qseqid == query.id & sub("gb\\|([^|]+).+", "\\1", bl$sseqid) == tr.exons$name[1], ]
    if(nrow(bl))
        rect( strand * bl$sstart, 2-bl.height, strand * bl$send, 2+bl.height, density=20 )
    axis(1)
    ## lets try to draw the alignment that we have from the initial alignments.
    ## I don't really need to have a special case for the seed species here; but
    ## 
    species.al <- tolower(sub("_", ".", species))
    align <- alignments[[ as.character(orth.row) ]]
    sp1 <- sub(".", "_", align$sp1, fixed=TRUE)
    substr(sp1, 1,1) <- toupper(substr(sp1, 1,1))
    if(length(align)){
        if(align$sp1 == species.al){
            segments( int.i.beg + seed.coords[1], 1, int.i.beg + seed.coords[2], 1, col='red', lwd=4 ) 
        }else{
            al.i <- which( align$sp == species.al )
            if(length(al.i) == 1){
                al.pos <- align$al[[al.i]]$pos[1,]
                if( as.double(align$l1) * as.double(align$l2[al.i]) > max.nm )
                    segments( int.i.beg, 1, int.i.end, 1, lwd=4, col=rgb(0.6, 0.6, 0.6) )
                if(!is.null( al.pos ))
                    segments( int.i.beg + al.pos[3], 1, int.i.beg + al.pos[4], 1, col='red', lwd=4 )
            }
        }
    }
    mtext(species, font=2)
    with(par(), text(usr[1], usr[4],
                     paste( tr.description[,c('gene.id', 'gene', 'description')], collapse=" "),
                     adj=c(0,2) ))
    ## we are now repeating ourselves a bit, but, let's have a go at this
    plot.new()
    if(is.na( orth.l[orth.row, species] ) || is.na(orth.l[orth.row, sp1])){
        plot.window(xlim=c(0,1), ylim=c(0,1))
        return()
    }
    plot.window( xlim=c(0, orth.l[orth.row, sp1]), ylim=c(0, orth.l[orth.row, species]), asp=1 )
    al.i <- which( align$sp == species.al )
    segments( 1, 1, c(orth.l[orth.row, sp1], 1), c(1, orth.l[orth.row, species]), lwd=2 )
    lines( seed.coords[1:2], c(1,1), lwd=2, col='red' )
    axis(1)
    axis(2)
    if(length(al.i) == 1){
        al.pos <- align$al[[al.i]]$pos
        if(!is.null( al.pos ))
            segments( al.pos[,'a_beg'], al.pos[,'b_beg'], al.pos[,'a_end'], al.pos[,'b_end'], col='black', lwd=1 )
        ## if( align$l1 * align$l2[al.i] > max.nm )
        ##     segments( int.i.beg, 1, int.i.end, 1, lwd=4, col=rgb(0.6, 0.6, 0.6) )
    }
}

## take a fasta data entry as an argument and then open the appropriate blast files
## and all sorts of global variables
plot.coverages <- function(fa.i, sp.b, alignments=aligns$long, start.from=1, interactive=TRUE){
    ## a function that returns a depth of coverage vector
    al.cov <- function(bl, id, l, max.e=1e-5, min.pid=65, exclude.sub=NULL,
                       tform=log.tf){
        bl <- bl[ bl$qseqid == id & bl$qlen == l, ]
        if(!is.null(exclude.sub))
            bl <- bl[ !grepl(exclude.sub, bl$sseqid), ]
##        l <- bl[1,'qlen']
        cov <- rep(0, l)
        if(!nrow(bl))
            return(cov)
        for(i in 1:nrow(bl)){
            r <- bl[i,'qstart']:bl[i,'qend']
            cov[r] <- cov[r] + 1
        }
        tform(cov)
    }
    plot.cov <- function(cov, main=NULL, poly.col=rgb(0.7, 0.7, 0.7), ...){
        pts.polygon <- function(y, y.min=0){
            x <- c(1, 1:length(y), length(y))
            y <- c(0, y, 0)
            cbind(x=x, y=y)
        }
        ylim=c(0, max(cov)) * 2
        if(ylim[2] == 0)
            ylim[2] <- 1
        if(is.null(nrow(cov))){
            plot(1:length(cov), cov, type='n', ylim=ylim, main=main, xaxs='i', yaxs='i',
                 frame.plot=FALSE, yaxt='n', ...)
            polygon(pts.polygon(cov), col=poly.col)
        }else{
            plot(1:nrow(cov), 1:nrow(cov), type='n', ylim=ylim, main=main, xaxs='i', yaxs='i',
                 frame.plot=FALSE, yaxt='n', ...)
            for(i in 1:ncol(cov))
                polygon( pts.polygon(cov), col=poly.col)
        }
        axis(side=2, at=seq(0, max(cov), ceiling( max(cov) / 4 )))
    }
##    
    fa.d <- query.seq[[fa.i]]
    fa.cl <- fa.files.cl[fa.i]
    fa.sm <- fa.files.sm[fa.i]
    bl.b <- bl.files.cl == fa.cl & bl.files.sm == fa.sm
    bl.i1 <- which(bl.b &  bl.files.sp == "Danio_rerio")
    bl.i2 <- which(bl.b &  bl.files.sp != "Danio_rerio" &
                   sp.b[ bl.files.sp ])
    qnt.probs <- seq(0, 1, 0.1)
    cov.qnts <- matrix(0, nrow=length(fa.d$id), ncol=length(qnt.probs))
    colnames(cov.qnts) <- names( quantile(1:100, qnt.probs))
    for(i in start.from:length( fa.d$id )){
        sp <- fa.d$sp2[i]
        id <- fa.d$id[i]
        tr <- fa.d$tr[i]
        al.i <- i
        while( (alignments[[al.i]]$seq.meta['tr','danio.rerio'] != tr || orth.i[ alignments[[al.i]]$i, 'Danio_rerio'] != fa.d$intron[i]) &&
              al.i <= length(fa.d$id) )
            al.i <- al.i + 1
        if( al.i > length(fa.d$id)){
            plot.new()
            plot.new()
            plot.new()
            next
        }
        ## this is so ugly; and because I did not include this information
        ## when I exported the data.
        orth.row <- alignments[[al.i]]$i
##        orth.row <- which(orth.tr[,'Danio_rerio'] == tr)[fa.d$intron[i]]
        tr.2 <- orth.tr[orth.row, bl.files.sp[ bl.i2 ]]
        tr.2.i <- orth.i[orth.row, bl.files.sp[ bl.i2 ]]
        ## if(length(bl.i2) != 1)
        ##     bl.i2 <- which(bl.b & bl.files.sp == "Oryzias_latipes")
        ## if(length(bl.i2) != 1)
        ##     stop(paste("What?? i is ", i))
        ## rows.1 <-  bl.data[[bl.i1]]$qseqid == id
        ## rows.2 <-  bl.data[[bl.i2]]$qseqid == id
        q.l <- nchar(query.seq[[fa.i]]$seq[i])
        cov.1 <- al.cov( bl.data[[bl.i1]], id, q.l, exclude.sub='ALT', tform=eval)
        cov.qnts[i, ] <- quantile( cov.1, qnt.probs, na.rm=TRUE )
        cov.1 <- log.tf(cov.1)
        cov.2 <- sapply( bl.data[bl.i2], al.cov, id=id, l=q.l )
##        cov.2 <- al.cov( bl.data[[bl.i2]], id, nchar(query.seq[[fa.i]]$seq[i]) )
        par(oma=c(2.1, 2.1, 2.1, 2.1))
        par(mar=c(2.1, 3.1, 2.1, 2.1))
        layout.nr <- 1 + ncol(cov.2)
        layout.m <- cbind( 1:layout.nr, matrix( layout.nr + 1:(2*layout.nr), ncol=2, byrow=TRUE ))
##        par(mfrow=c(1 + ncol(cov.2),1))
        layout(layout.m, widths=c(0.2, 0.8, 0.2))
        plot.cov(cov.1, main='')
        for(j in 1:ncol(cov.2))
            plot.cov( cov.2[,j], main='') ## bl.files.sp[ bl.i2[j] ] )
#        plot.cov(cov.2, main='Others')
##        inpt <- readline(paste(i, ":"))
##        layout( matrix(1:(2 * (1 + length(bl.i2))), ncol=2, byrow=TRUE), widths=c(0.8, 0.2) )
        par(mar=c(2.1, 1.1, 2.1, 2.1))
        cat(i, ": ")
        plot.transcript.region( orth.row, 'Danio_rerio', fa.d$coords[i,], bl.data[[bl.i1]], id, alignments )
        for(j in bl.i2){
            cat(j, ", ")
            plot.transcript.region( orth.row, bl.files.sp[j], fa.d$coords[i,], bl.data[[j]], id, alignments )
        }
        cat("\n")
        if(interactive){
            inpt <- readline(paste(i, ":"))
            if(inpt == 'q')
                break
        }
    }
    invisible( cov.qnts )
}
        


cairo_pdf("tel_long__blast.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
tel.long.dr.bl <- plot.coverages( which( names(query.seq) == 'tel_long'), sp.class[,'teleostei'], aligns$long, interactive=FALSE )
dev.off()

plot.coverages( which( names(query.seq) == 'tel_long'), sp.class[,'teleostei'], aligns$long, start.from=6 )
plot.coverages( which( names(query.seq) == 'tel_long'), sp.class[,'mammalia'], aligns$long )

plot.coverages( which( names(query.seq) == 'mam_long'), sp.class[,'mammalia'], aligns$long, start.from=1 )

cairo_pdf("mam_long_blast.pdf",  width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
mam.long.dr.bl <- plot.coverages( which( names(query.seq) == 'mam_long'), sp.class[,'mammalia'], aligns$long, interactive=FALSE )
dev.off()

plot.coverages( which( names(query.seq) == 'sau_long'), sp.class[,'sauria'], aligns$long )

plot.coverages( which( names(query.seq) == 'tel_sampled'), sp.class[,'teleostei'], aligns$sampled, start.from=159 )

cairo_pdf("tel_sampled_blast.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
tel.sampled.dr.bl <- plot.coverages( which( names(query.seq) == 'tel_sampled'), sp.class[,'teleostei'], aligns$sampled, interactive=FALSE )
dev.off()

cairo_pdf("mam_sampled_blast.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
mam.sampled.dr.bl <- plot.coverages( which( names(query.seq) == 'mam_sampled'), sp.class[,'mammalia'], aligns$sampled, interactive=FALSE )
dev.off()

plot.coverages( which(names(query.seq) == 'tel_ctl'), sp.class[,'teleostei'], aligns.ctl, start.from=400 )

cairo_pdf("tel_ctl_blast.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
tel.ctl.dr.bl <- plot.coverages( which(names(query.seq) == 'tel_ctl'), sp.class[,'teleostei'], aligns.ctl, interactive=FALSE )
dev.off()

cairo_pdf("tel_mam_blast.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
mam.ctl.dr.bl <- plot.coverages( which(names(query.seq) == 'mam_ctl'), sp.class[,'mammalia'], aligns.ctl, interactive=FALSE )
dev.off()


plot.coverages( which( names(query.seq) == 'tel_sampled'), sp.class[,'teleostei'], aligns$sampled )

## we can now have a quick look at..
## These are the scores where we remove those that have no alignments
## which should make them 
tel.long.sc <- unlist(sapply( aligns.top.cl$tel$long, function(x){
    if(!length(x)) return(NULL)
    x[[1]]$coords['score'] }))

## take a section of aligns.top.cl 
get.mn <- function(x){
    ex.func <- function(y){
        if(!length(y)) return(NULL)
        if(is.null(y[[1]]$i)) return(0)
        l <- orth.l[ y[[1]]$i, sp.name( c(y[[1]]$sp1, y[[1]]$sp2) ) ]
        as.double(l[1]) * as.double(l[2])
    }
    unlist( lapply( x, ex.func ) )
}

## it is actually more reasonable to calculate:
## m * sum(n) as that represents the true search space used
get.mn.2 <- function(x){
    ex.func <- function(y){
        if(!length(y)) return(NULL)
        if(is.null(y[[1]]$i)) return(0)
        l1 <- as.double( orth.l[ y[[1]]$i, sp.name( y[[1]]$sp1) ])
        sp2 <- sp.name( unlist( sapply(y, function(z){ z$sp2 })))
        l2 <- as.double( orth.l[ y[[1]]$i, sp2 ] )
        l1 * sum(l2)
    }
    unlist( lapply(x, ex.func))
}

tel.long.mn <- get.mn( aligns.top.cl$tel$long )
tel.sampled.mn <- get.mn( aligns.top.cl$tel$sampled )
tel.ctl.mn <- get.mn( aligns.ctl.top.cl$tel )

mam.long.mn <- get.mn( aligns.top.cl$mam$long )
mam.sampled.mn <- get.mn( aligns.top.cl$mam$sampled )
mam.ctl.mn <- get.mn( aligns.ctl.top.cl$mam )

tel.long.mn.2 <- get.mn.2( aligns.top.cl$tel$long )
tel.sampled.mn.2 <- get.mn.2( aligns.top.cl$tel$sampled )
tel.ctl.mn.2 <- get.mn.2( aligns.ctl.top.cl$tel )

mam.long.mn.2 <- get.mn.2( aligns.top.cl$mam$long )
mam.sampled.mn.2 <- get.mn.2( aligns.top.cl$mam$sampled )
mam.ctl.mn.2 <- get.mn.2( aligns.ctl.top.cl$mam )


## I don't need this as I have the coords table in the various query.seq data structures..
plot( log2(1+tel.long.dr.bl[,'50%']), log2(query.seq$tel_long$coords[,'score']) )
plot( log2(1+tel.sampled.dr.bl[,'50%']), log2(query.seq$tel_sampled$coords[,'score']) )
plot( log2(1+tel.ctl.dr.bl[,'50%']), log2(query.seq$tel_ctl$coords[,'score']) )

plot( log2(1+tel.long.dr.bl[,'80%']), log2(query.seq$tel_long$coords[,'score']) )
plot( log2(1+tel.sampled.dr.bl[,'80%']), log2(query.seq$tel_sampled$coords[,'score']) )
plot( log2(1+tel.ctl.dr.bl[,'80%']), log2(query.seq$tel_ctl$coords[,'score']) )


## The blast alignment depth should be between 1 and 2 (i.e. allow two copies) to remove
## common sequences
tel.long.b <- tel.long.dr.bl[,'50%'] > 0 & tel.long.dr.bl[,'50%'] < 3
tel.sampled.b <- tel.sampled.dr.bl[,'50%'] > 0 & tel.sampled.dr.bl[,'50%'] < 3
tel.ctl.b <- tel.ctl.dr.bl[,'50%'] > 0 & tel.ctl.dr.bl[,'50%'] < 3
mam.long.b <- mam.long.dr.bl[,'50%'] > 0 & mam.long.dr.bl[,'50%'] < 3
mam.sampled.b <- mam.sampled.dr.bl[,'50%'] > 0 & mam.sampled.dr.bl[,'50%'] < 3
mam.ctl.b <- mam.ctl.dr.bl[,'50%'] > 0 & mam.ctl.dr.bl[,'50%'] < 3


hist( log2( query.seq$tel_long$coords[tel.long.b, 'score'] ))
hist( log2( query.seq$tel_sampled$coords[tel.sampled.b, 'score'] ))
hist( log2( query.seq$tel_ctl$coords[tel.ctl.b, 'score'] ))


tel.lm.long <- lm( log2(query.seq$tel_long$coords[tel.long.b, 'score']) ~ log2(tel.long.mn[tel.long.b]))

tel.lm.sampled <- lm( log2(query.seq$tel_sampled$coords[tel.sampled.b, 'score']) ~ log2(tel.sampled.mn[tel.sampled.b]))

tel.lm.ctl <- lm( log2(query.seq$tel_ctl$coords[tel.ctl.b, 'score']) ~ log2(tel.ctl.mn[tel.ctl.b]))

tel.lm.ctl.2 <- lm( log2(c(query.seq$tel_ctl$coords[tel.ctl.b, 'score'],
                           query.seq$tel_sampled$coords[tel.sampled.b, 'score']))
                         ~ log2(c(tel.ctl.mn[tel.ctl.b], tel.sampled.mn[tel.sampled.b])) )


mam.lm.long <- lm( log2(query.seq$mam_long$coords[mam.long.b, 'score']) ~ log2(mam.long.mn[mam.long.b]))

mam.lm.sampled <- lm( log2(query.seq$mam_sampled$coords[mam.sampled.b, 'score']) ~ log2(mam.sampled.mn[mam.sampled.b]))

mam.lm.ctl <- lm( log2(query.seq$mam_ctl$coords[mam.ctl.b, 'score']) ~ log2(mam.ctl.mn[mam.ctl.b]))
mam.lm.ctl.2 <- lm( log2(c(query.seq$mam_ctl$coords[mam.ctl.b, 'score'],
                           query.seq$mam_sampled$coords[mam.sampled.b, 'score']))
                         ~ log2(c(mam.ctl.mn[mam.ctl.b], mam.sampled.mn[mam.sampled.b])) )



## try some more models using mn.2 instead:
tel.lm.long.2 <- lm( log2(query.seq$tel_long$coords[tel.long.b, 'score']) ~
                         log2(tel.long.mn.2[tel.long.b]))

tel.lm.sampled.2 <- lm( log2(query.seq$tel_sampled$coords[tel.sampled.b, 'score']) ~
                          log2(tel.sampled.mn.2[tel.sampled.b]))

#tel.lm.ctl.2 <- lm( log2(query.seq$tel_ctl$coords[tel.ctl.b, 'score']) ~
#                     log2(tel.ctl.mn.2[tel.ctl.b]))

tel.lm.ctl.2.2 <- lm( log2(c(query.seq$tel_ctl$coords[tel.ctl.b, 'score'],
                           query.seq$tel_sampled$coords[tel.sampled.b, 'score']))
                         ~ log2(c(tel.ctl.mn.2[tel.ctl.b], tel.sampled.mn.2[tel.sampled.b])) )


mam.lm.long.2 <- lm( log2(query.seq$mam_long$coords[mam.long.b, 'score']) ~
                         log2(mam.long.mn.2[mam.long.b]))

mam.lm.sampled.2 <- lm( log2(query.seq$mam_sampled$coords[mam.sampled.b, 'score']) ~
                            log2(mam.sampled.mn.2[mam.sampled.b]))

#mam.lm.ctl.2 <- lm( log2(query.seq$mam_ctl$coords[mam.ctl.b, 'score']) ~ log2(mam.ctl.mn.2[mam.ctl.b]))
mam.lm.ctl.2.2 <- lm( log2(c(query.seq$mam_ctl$coords[mam.ctl.b, 'score'],
                           query.seq$mam_sampled$coords[mam.sampled.b, 'score']))
                         ~ log2(c(mam.ctl.mn.2[mam.ctl.b], mam.sampled.mn.2[mam.sampled.b])) )
## these all show much stronger correlations..

## These are nice in that the long ones do not show a relationship:
summary(tel.lm.long)
## Coefficients:
##                               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                    7.58988    1.13565   6.683 1.15e-10 ***
## log2(tel.long.mn[tel.long.b])  0.04766    0.04567   1.044    0.298    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 1.139 on 298 degrees of freedom
## Multiple R-squared:  0.003641,	Adjusted R-squared:  0.0002976 
## F-statistic: 1.089 on 1 and 298 DF,  p-value: 0.2975
summary(tel.lm.sampled)
## Coefficients:
##                                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                           4.8087     0.4291   11.21  < 2e-16 ***
## log2(tel.sampled.mn[tel.sampled.b])   0.1052     0.0202    5.21 5.24e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.4808 on 176 degrees of freedom
## Multiple R-squared:  0.1336,	Adjusted R-squared:  0.1287 
## F-statistic: 27.14 on 1 and 176 DF,  p-value: 5.244e-07

summary(tel.lm.ctl)
## Coefficients:
##                             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                  3.17889    0.24148   13.16   <2e-16 ***
## log2(tel.ctl.mn[tel.ctl.b])  0.18476    0.01135   16.28   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.5528 on 821 degrees of freedom
## Multiple R-squared:  0.2441,	Adjusted R-squared:  0.2432 
## F-statistic: 265.1 on 1 and 821 DF,  p-value: < 2.2e-16
pf( 265.1, 1, 821, lower.tail=FALSE )
## [1] 7.287097e-52

summary(tel.lm.ctl.2)
## similar to the two above, but with better scores

summary(mam.lm.long)
## Coefficients:
##                               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                    8.54370    1.08278   7.891 3.85e-13 ***
## log2(mam.long.mn[mam.long.b]) -0.04001    0.04334  -0.923    0.357    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.7522 on 166 degrees of freedom
## Multiple R-squared:  0.005109,	Adjusted R-squared:  -0.0008846 
## F-statistic: 0.8524 on 1 and 166 DF,  p-value: 0.3572
summary(mam.lm.sampled)
## Coefficients:
##                                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                          4.49356    0.26129  17.197  < 2e-16 ***
## log2(mam.sampled.mn[mam.sampled.b])  0.09121    0.01332   6.849 1.23e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.3529 on 174 degrees of freedom
## Multiple R-squared:  0.2123,	Adjusted R-squared:  0.2078 
## F-statistic: 46.91 on 1 and 174 DF,  p-value: 1.231e-10


## there is probably a function that does this
pt.resids <- function(model, x, y){
    coeff <- model$coefficients
    y.m <- coeff[1] + coeff[2] * x
    data.frame('x'=x, 'y'=y, 'y.m'=y.m, 'resid'=y-y.m)
}

tel.long.res <- pt.resids( tel.lm.sampled,
                          log2(tel.long.mn[tel.long.b]),
                          log2( query.seq$tel_long$coords[tel.long.b, 'score']) )

tel.sampled.res <- pt.resids( tel.lm.sampled,
                          log2(tel.sampled.mn[tel.sampled.b]),
                          log2( query.seq$tel_sampled$coords[tel.sampled.b, 'score']) )

tel.ctl.res <- pt.resids( tel.lm.ctl,
                          log2(tel.ctl.mn[tel.ctl.b]),
                          log2(query.seq$tel_ctl$coords[tel.ctl.b, 'score']) )

tel.ctl.2.res <- pt.resids( tel.lm.ctl.2,
                           log2(c(tel.ctl.mn[tel.ctl.b], tel.sampled.mn[tel.sampled.b])),
                           log2(c(query.seq$tel_ctl$coords[tel.ctl.b, 'score'],
                                  query.seq$tel_sampled$coords[tel.sampled.b, 'score'])) )

mam.long.res <- pt.resids( mam.lm.sampled,
                          log2(mam.long.mn[mam.long.b]),
                          log2( query.seq$mam_long$coords[mam.long.b, 'score']) )

mam.sampled.res <- pt.resids( mam.lm.sampled,
                          log2(mam.sampled.mn[mam.sampled.b]),
                          log2( query.seq$mam_sampled$coords[mam.sampled.b, 'score']) )

mam.ctl.2.res <- pt.resids( mam.lm.ctl.2,
                           log2(c(mam.ctl.mn[mam.ctl.b], mam.sampled.mn[mam.sampled.b])),
                           log2(c(query.seq$mam_ctl$coords[mam.ctl.b, 'score'],
                                  query.seq$mam_sampled$coords[mam.sampled.b, 'score'])) )

#### and do the same for model.2 copying and pasting in an absolutely horrible
#### manner:
tel.long.res.2 <- pt.resids( tel.lm.sampled.2,
                          log2(tel.long.mn.2[tel.long.b]),
                          log2( query.seq$tel_long$coords[tel.long.b, 'score']) )

tel.sampled.res.2 <- pt.resids( tel.lm.sampled.2,
                          log2(tel.sampled.mn.2[tel.sampled.b]),
                          log2( query.seq$tel_sampled$coords[tel.sampled.b, 'score']) )

## tel.ctl.res.2 <- pt.resids( tel.lm.ctl.2,
##                           log2(tel.ctl.mn.2[tel.ctl.b]),
##                           log2(query.seq$tel_ctl$coords[tel.ctl.b, 'score']) )

tel.ctl.2.res.2 <- pt.resids( tel.lm.ctl.2.2,
                           log2(c(tel.ctl.mn.2[tel.ctl.b], tel.sampled.mn.2[tel.sampled.b])),
                           log2(c(query.seq$tel_ctl$coords[tel.ctl.b, 'score'],
                                  query.seq$tel_sampled$coords[tel.sampled.b, 'score'])) )

mam.long.res.2 <- pt.resids( mam.lm.sampled.2,
                          log2(mam.long.mn.2[mam.long.b]),
                          log2( query.seq$mam_long$coords[mam.long.b, 'score']) )

mam.sampled.res.2 <- pt.resids( mam.lm.sampled.2,
                          log2(mam.sampled.mn.2[mam.sampled.b]),
                          log2( query.seq$mam_sampled$coords[mam.sampled.b, 'score']) )

mam.ctl.2.res.2 <- pt.resids( mam.lm.ctl.2.2,
                           log2(c(mam.ctl.mn.2[mam.ctl.b], mam.sampled.mn.2[mam.sampled.b])),
                           log2(c(query.seq$mam_ctl$coords[mam.ctl.b, 'score'],
                                  query.seq$mam_sampled$coords[mam.sampled.b, 'score'])) )


###

q.prob <- seq(0, 1, 0.05)
tel.long.res.q <- quantile( tel.long.res[,'resid'], probs=q.prob )
tel.sampled.res.q <- quantile( tel.sampled.res[,'resid'], probs=q.prob )
tel.ctl.res.q <- quantile( tel.ctl.res[,'resid'], probs=q.prob )
tel.ctl.2.res.q <- quantile( tel.ctl.2.res[,'resid'], probs=q.prob )
mam.long.res.q <- quantile( mam.long.res[,'resid'], probs=q.prob )
mam.sampled.res.q <- quantile( mam.sampled.res[,'resid'], probs=q.prob )
mam.ctl.2.res.q <- quantile( mam.ctl.2.res[,'resid'], probs=q.prob )

cairo_pdf("alignment_scores_filtered.pdf", width=a4.w * pdf.m * 0.6, height=a4.w * pdf.m * 0.6 )
par(mfrow=c(2,2))
xlim <- range( c(tel.long.res[,'x'], tel.sampled.res[,'x'], tel.ctl.res[,'x']) )
ylim <- range( c(tel.long.res[,'y'], tel.sampled.res[,'y'], tel.ctl.res[,'y']) )
with(tel.ctl.2.res, plot(x, y, xlim=xlim, ylim=ylim,
     col=ifelse(resid >= tel.ctl.2.res.q['95%'], 'red', 'black'),
     xlab='log2 mn', ylab='log2 score'))
abline( tel.lm.ctl.2, col='red', lwd=2 )
with(par(), mtext("A", at=usr[1], cex=mt.cex, line=1))
with(tel.long.res, plot(x, y, xlim=xlim, ylim=ylim,
     col=ifelse(resid >= tel.ctl.2.res.q['95%'], 'red', 'black'),
     xlab='log2 mn', ylab='log2 score'))
abline( tel.lm.ctl.2, col='red', lwd=2 )
with(par(), mtext("B", at=usr[1], cex=mt.cex, line=1))
##
xlim <- range( c(mam.long.res[,'x'], mam.sampled.res[,'x']) )
ylim <- range( c(mam.long.res[,'y'], mam.sampled.res[,'y']) )
with(mam.ctl.2.res, plot(x, y, xlim=xlim, ylim=ylim,
     col=ifelse(resid >= mam.ctl.2.res.q['95%'], 'red', 'black'),
     xlab='log2 mn', ylab='log2 score'))
abline( mam.lm.ctl.2, col='red', lwd=2 )
with(par(), mtext("D", at=usr[1], cex=mt.cex, line=1))
with(mam.long.res, plot(x, y, xlim=xlim, ylim=ylim,
     col=ifelse(resid >= mam.ctl.2.res.q['95%'], 'red', 'black'),
     xlab='log2 mn', ylab='log2 score'))
abline( mam.lm.ctl.2, col='red', lwd=2 )
with(par(), mtext("C", at=usr[1], cex=mt.cex, line=1))
dev.off()

### and for the second models:
q.prob <- seq(0, 1, 0.05)
tel.long.res.q.2 <- quantile( tel.long.res.2[,'resid'], probs=q.prob )
tel.sampled.res.q.2 <- quantile( tel.sampled.res.2[,'resid'], probs=q.prob )
##tel.ctl.res.q.2 <- quantile( tel.ctl.res.2[,'resid'], probs=q.prob )
tel.ctl.2.res.q.2 <- quantile( tel.ctl.2.res.2[,'resid'], probs=q.prob )
mam.long.res.q.2 <- quantile( mam.long.res.2[,'resid'], probs=q.prob )
mam.sampled.res.q.2 <- quantile( mam.sampled.res.2[,'resid'], probs=q.prob )
mam.ctl.2.res.q.2 <- quantile( mam.ctl.2.res.2[,'resid'], probs=q.prob )

### This analysis has now become too messy and needs to be refilteres. But it seems quite clear that
### using m * sum(n) is the way to go. The models have much better correlations..
cairo_pdf("alignment_scores_filtered_2.pdf", width=a4.w * pdf.m * 0.6, height=a4.w * pdf.m * 0.6 )
par(mfrow=c(2,2))
xlim <- range( c(tel.long.res.2[,'x'], tel.sampled.res.2[,'x'], tel.ctl.res.2[,'x']) )
ylim <- range( c(tel.long.res.2[,'y'], tel.sampled.res.2[,'y'], tel.ctl.res.2[,'y']) )
with(tel.ctl.2.res.2, plot(x, y, xlim=xlim, ylim=ylim,
                           col=ifelse(resid >= tel.ctl.2.res.q.2['95%'], 'red', 'black'),
                           xlab='log2 mn', ylab='log2 score'))
abline( tel.lm.ctl.2.2, col='red', lwd=2 )
with(par(), mtext("A", at=usr[1], cex=mt.cex, line=1))
with(tel.long.res.2, plot(x, y, xlim=xlim, ylim=ylim,
     col=ifelse(resid >= tel.ctl.2.res.q.2['95%'], 'red', 'black'),
     xlab='log2 mn', ylab='log2 score'))
##
abline( tel.lm.ctl.2.2, col='red', lwd=2 )
with(par(), mtext("B", at=usr[1], cex=mt.cex, line=1))
##
xlim <- range( c(mam.long.res.2[,'x'], mam.ctl.2.res.2[,'x']) )
ylim <- range( c(mam.long.res.2[,'y'], mam.ctl.2.res.2[,'y']) )
## 
with(mam.ctl.2.res.2, plot(x, y, xlim=xlim, ylim=ylim,
     col=ifelse(resid >= mam.ctl.2.res.q.2['95%'], 'red', 'black'),
     xlab='log2 mn', ylab='log2 score'))
abline( mam.lm.ctl.2, col='red', lwd=2 )
with(par(), mtext("C", at=usr[1], cex=mt.cex, line=1))
with(mam.long.res.2, plot(x, y, xlim=xlim, ylim=ylim,
     col=ifelse(resid >= mam.ctl.2.res.q.2['95%'], 'red', 'black'),
     xlab='log2 mn', ylab='log2 score'))
abline( mam.lm.ctl.2.2, col='red', lwd=2 )
with(par(), mtext("D", at=usr[1], cex=mt.cex, line=1))
dev.off()


############ Looking at these plots, it looks like genes with long introns mostly contain long introns
########### This can be checked by doing
###
## sum( log( n_a / n_b ))
## for the introns of each gene where n_a and n_b are the number of introns in the genome
## which are equal to and shorter or longer respectively. A random selection should tend
## to 0, where the deviation from 0 should be possible to estimate for a given number of
## introns.

## ed : the empirical distribution (all values sorted)
## s  : the sampled values
lr.div <- function(ed, s){
    if(is.unsorted(ed))
        ed <- sort(ed)
    ## s.i1 The first index of the value
    ## s.i2 The last index of the value
    ## Hence the number that are equal to or smaller
    ## is s.i2
    ## and the the number that are equal to or larger
    ## is (1 + length - s.i1)
    s.i1 <- match(s, ed)
    s.i2 <- 1 + length(ed) - match(s, rev(ed))
    sum( log2( s.i2 / (1 + length(ed) - s.i1) ) )
}

tmp <- sort(orth.l[, 'Danio_rerio'])

tmp.lr <- sapply(1:10000, function(i){ lr.div(tmp, sample(tmp, 10)) })

hist(tmp.lr)

tmp.s <- sapply(1:10000, function(i){ sample(tmp, 10) })
tmp.lr2 <- apply(tmp.s, 2, function(x){ lr.div(tmp, x) })

par(mfrow=c(1,2))
hist(tmp.lr)
hist(tmp.lr2)

tmp.dr <- tapply( orth.l[,'Danio_rerio'], orth.id[,'Danio_rerio'], function(x){
    c('n'=length(x), 'lr'=lr.div(tmp, x)) })

tmp.dr <- t(sapply( tmp.dr, eval ))

hist(tmp.dr[,'n'])
hist(tmp.dr[,'lr'])

hist(tmp.dr[ tmp.dr[,'n'] == 10, 'lr'], breaks=30)
hist(tmp.dr[ tmp.dr[,'n'] == 8, 'lr'])

image(log2(1 + table(tmp.dr[,'n'], as.integer(tmp.dr[,'lr']))))
## but actually I do not see very much here. Some assymetry in the
## distribution certainly, but nothing obvious. 

