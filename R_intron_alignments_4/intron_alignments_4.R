## perform similar sets of alignments as in intron_alignments, but this
## time use Danio rerio reverse complement sequences.

## run alignments of long introns in parallel

library(parallel)

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
source("../R_intron_alignments_2/functions.R")

max.l <- 1e8
## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2


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
                          sm=sub.matrix, max.l=1e7, r.c=FALSE){
    fam.id <- orth.fam[i]
    tr <- orth.tr[i,]
    int.i <- orth.i[i,]
    
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
        a.int.seq <- read.exons( intron.f[ orth.fam[j] ] )
        a.int.meta <- sapply( a.int.seq, function(x){ c('id'=x$id, 'tr'=x$tr,
                                                        'sp'=x$sp, 'n'=length(x$l) )})
        a.seq.b <- a.int.meta['tr',] == orth.tr[j,al.sp]
        a.seq <- a.int.seq[[ which(a.seq.b) ]]$e[ orth.i[j,al.sp] ]
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

## variance in teleosts and 
tel.var <-  read.table("../R_172_genomes/dr_teleost_var.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)

## the files containing intron sequences
intron.f <-  read.table("../R_172_genomes/intron_files.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
tmp <- rownames(intron.f)
intron.f <- intron.f[,1]
names(intron.f) <- tmp
rm(tmp)

## some of the following are used by the align.introns function
orth.fam <- read.table("../R_172_genomes/dr_intron_fam.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)[,1]
orth.id <-  read.table("../R_172_genomes/dr_intron_orthology_id.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.tr <-  read.table("../R_172_genomes/dr_intron_orthology_tr.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.i <-  read.table("../R_172_genomes/dr_intron_orthology_i.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.l <-  read.table("../R_172_genomes/dr_intron_orthology_l.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)

## numbers in tel.var are log2 transformed. Numbers in orth.l are
## linear..
b <- tel.var[,'min'] > 10 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10 &
    !is.na(orth.l[,'danio.rerio'])
sum(b)
## 1255

## the test indices
test.i <- which(b)
test.o <- order( tel.var[b,'sd'] )
test.i <- test.i[test.o]

## The control
ctl.b <- tel.var[,'X50.'] < 8 & orth.l[,'danio.rerio'] > 1024 & tel.var[,'n'] > 10
sum(ctl.b, na.rm=TRUE)  ## 3130
ctl.i <- which(ctl.b)
ctl.ii <- sample(ctl.i, size=500)

sub.matrix <- make.sub.matrix()

do.rc.alignment <- function(i, max.l=1e8, gap=c(-10,-2), al.sp='danio.rerio', ctl=FALSE, ...){
    if(!ctl)
        align.introns(i, al.sp=al.sp, gap=gap, r.c=TRUE, max.l=max.l)
    else
        align.introns.ctl(i, al.sp=al.sp, gap=gap, max.l=max.l)
}

aligns.rc <- vector(mode='list', length=4)

aligns.rc[[1]] <- mclapply( test.i[1:500], do.rc.alignment, max.l=max.l, mc.cores=10 )
aligns.rc[[2]] <- mclapply( ctl.ii, do.rc.alignment, max.l=max.l, mc.cores=10 )

## the controls rely on the following datat structures:
## we will use ranks to find similar sized control alignments
## use ranks
ga.int.l <- orth.l[,'danio.rerio']
ga.int.l[ is.na(ga.int.l) ] <- 0
ga.int.l.o <- order( ga.int.l )
ga.int.l.r <- rank( ga.int.l, ties.method='random' )
ga.int.l.or <- (1:length(ga.int.l))[ga.int.l.o]


aligns.rc[[3]] <- mclapply( test.i[1:500], do.rc.alignment, max.l=max.l, mc.cores=10, ctl=TRUE )
aligns.rc[[4]] <- mclapply( ctl.ii, do.rc.alignment, max.l=max.l, mc.cores=10, ctl=TRUE )
## These finished with some errors, but think this is in the R-code; for where no alignment
## was returned.

## use the same names for these as we have used before:
names(aligns.rc) <- c('long', 'sampled', 'long.ctl', 'sampled.ctl')

## to Visualise the results:
sp.class <- get.sp.taxonomy()
##

##
aligns.rc.scores <- lapply(aligns.rc, extract.scores.l, max.l=max.l)
##

aligns.rc.scores.cl <- lapply(colnames(sp.class), function(cl){
    sp.b <- sp.class[,cl]
    names(sp.b) <- rownames(sp.class)
    lapply(aligns.rc, extract.scores.l, max.l=max.l, sp.b=sp.b)
})
names(aligns.rc.scores.cl) <- colnames(sp.class)

##
### let us make a different set of quantiles that we can slice
### the distributions up by.
##
qnt.probs <- seq(0, 1, 0.05)
##
all.qnt <- quantile( do.call(c, lapply(aligns.rc.scores, function(x){ x[,'mn'] })), probs=qnt.probs )
long.qnt <- quantile( do.call(c, lapply(aligns.rc.scores[c(1,3)], function(x){ x[,'mn'] })), probs=qnt.probs )
sampled.qnt <- quantile( do.call(c, lapply(aligns.rc.scores[c(2,4)], function(x){ x[,'mn'] })), probs=qnt.probs )
##
qnt <- list('all'=all.qnt, 'long'=long.qnt, 'sampled'=sampled.qnt)
##

## we can then use these to make sliced histograms...
aligns.rc.scores.mnq <- lapply(aligns.rc.scores, function(x){
    lapply(qnt, function(q){ mnq.class(x, q) }) })
##
aligns.rc.scores.cl.mnq <- lapply(aligns.rc.scores.cl, function(al){
    lapply(al, function(x){
            lapply(qnt, function(q){ mnq.class(x, q) }) }) })


##
### can we use the same quantiles for both the long and sampled
### alignments?
##
sapply( aligns.rc.scores.mnq, function(x){
    table(c(2:length(all.qnt) - 1, x$all)) - 1 })

##    long sampled long.ctl sampled.ctl
## 1     4    5774        4        5782
## 2     4    5792        4        5785
## 3     3    5787        3        5787
## 4     7    5782        7        5781
## 5     5    5783        5        5783
## 6     9    5779        9        5780
## 7    14    5758       14        5757
## 8    15    5788       15        5788
## 9    32    5758       32        5758
## 10   88    5701       88        5701
## 11  414    5374      414        5374
## 12 1769    4019     1769        4019
## 13 3771    2017     3773        2017
## 14 4822     966     4821         967
## 15 5277     512     5277         511
## 16 5521     266     5523         266
## 17 5653     134     5655         134
## 18 5692      96     5693          96
## 19 5730      60     5727          60
## 20 5738      48     5743          48

## here we have a really big difference in the numbers of alignments
## for the different classes. We may need to do things differently,
## but for now, let us see what we get.

aligns.rc.scores.h <- hist( log2(do.call( c, lapply( aligns.rc.scores, function(x){ x[,'score'] }))),
                        breaks=20, plot=FALSE )
##
aligns.rc.scores.qs1.h <- lapply( 1:length(aligns.rc.scores), function(i){
    sliced.h( log2(aligns.rc.scores[[i]][,'score']), slices=aligns.rc.scores.mnq[[i]]$all,
             breaks=aligns.rc.scores.h$breaks )})
names(aligns.rc.scores.qs1.h) <- names(aligns.rc.scores)
##
##
aligns.rc.scores.cl.qs1.h <- lapply( 1:length(aligns.rc.scores.cl), function(i){
    lapply( 1:length(aligns.rc.scores.cl[[i]]), function(j){
        tmp <- sliced.h( log2(aligns.rc.scores.cl[[i]][[j]][,'score']),
                        slices=aligns.rc.scores.cl.mnq[[i]][[j]]$all,
                        breaks=aligns.rc.scores.h$breaks )
        names(tmp) <- names(aligns.rc.scores.cl[[i]])
        tmp})
})
names(aligns.rc.scores.cl.qs1.h) <- names(aligns.rc.scores.cl)


plot.dens <- TRUE
par(mfrow=c(2,2))
plot.counts( aligns.rc.scores.qs1.h[[1]], dens=plot.dens )
plot.counts( aligns.rc.scores.qs1.h[[2]], dens=plot.dens )
plot.counts( aligns.rc.scores.qs1.h[[3]], dens=plot.dens )
plot.counts( aligns.rc.scores.qs1.h[[4]], dens=plot.dens )

## looks like very, very little there.
## plot example alignments:
## to do this more elegantly we should have better species information
sp.class.3 <- readRDS("../R_trees_distances/sp_class_3.rds")
class.col.3 <- readRDS("../R_trees_distances/class_col_3.rds")
sp.col.3 <-  readRDS("../R_trees_distances/sp_col_3.rds")
sp.y <-  readRDS("../R_trees_distances/sp_y.rds")

## mangle the names to fit
rownames(sp.class.3) <- sub("_", ".", rownames(sp.class.3), fixed=TRUE )
names(sp.col.3) <- sub("_", ".", names(sp.col.3), fixed=TRUE )
names(sp.y) <- sub("_", ".", names(sp.y), fixed=TRUE )

## In order to know which family and gene we need to include the following information
orth.fam <- read.table("../R_172_genomes/dr_intron_fam.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)[,1]
orth.id <-  read.table("../R_172_genomes/dr_intron_orthology_id.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.tr <-  read.table("../R_172_genomes/dr_intron_orthology_tr.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.i <-  read.table("../R_172_genomes/dr_intron_orthology_i.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
orth.l <-  read.table("../R_172_genomes/dr_intron_orthology_l.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)

### lets see if we can plot some of these alignments and what we can do
### with them.
alignment.figure <- function(al, sp.col=sp.col.3, class.col=class.col.3,
                             col.tform=log2, col.f=score.col,
                             id.tbl=orth.id, tr.tbl=orth.tr,
                             i.tbl=orth.i, l.tbl=orth.l){
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
    text(0.05, 100 - 1.5 * strheight("A") * 1:length(class.col.3),
         capitalise(names(class.col.3)), col=class.col, adj=c(0,0), font=2)
    hm.v <- seq( tmp$cols$range[1], tmp$cols$range[2], length.out=100 )
    hm.y <- seq( 0, 25, length.out=length(hm.v))
    hm.yd <- diff(hm.y)[1]
    rect( 0.85, hm.y, 1.0, hm.y + hm.yd, col=col.f(hm.v)$cols, border=NA )
    lab.i <- seq(1, length(hm.v), length.out=5)
    text( 0.82, hm.y[ lab.i ] + hm.yd/2, sprintf("%.1e", hm.v[lab.i]), adj=c(1,0.5),
         cex=0.85)
}


par(mfrow=c(1,1))

for(i in 1:length(aligns.rc$long)){
    alignment.figure( aligns.rc$long[[i]] )
    inpt <- readline(paste(i, ": "))
    if(inpt == 'q')
        break
}
