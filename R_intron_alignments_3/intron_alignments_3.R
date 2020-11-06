## perform similar sets of alignments as in intron_alignments, but this
## time use gasterosteus aculeates as the alignments species.

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
                          sm=sub.matrix, max.l=1e7){
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
    b.seq <- sapply( names(int.i[-al.i]), function(sp){
        int.seq[[sp]]$e[ as.numeric(int.i[sp]) ] })
    b.seq <- b.seq[ !is.na(b.seq) ]
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
    !is.na(orth.l[,'gasterosteus.aculeatus'])
sum(b)
## 909
test.i <- which(b)
test.o <- order( tel.var[b,'sd'] )

## it would probably be better to use tel.var[,'X30.'] < 8
## which gives 492 rows.
b2 <- tel.var[,'X50.'] < 9 & orth.l[,'gasterosteus.aculeatus'] > 1024 & tel.var[,'n'] > 10
sum(b2, na.rm=TRUE)  ## 509
ctl.i <- which(b2)

sub.matrix <- make.sub.matrix()

test.al <- mclapply( test.i[ test.o[1:500] ], align.introns, max.l=1e8,
                    gap=c(-10, -2), al.sp='gasterosteus.aculeatus',
                    mc.cores=14 )
##
ctl.al <- mclapply( ctl.i[ 1:500 ], align.introns, max.l=1e8,
                    gap=c(-10, -2), al.sp='gasterosteus.aculeatus',
                    mc.cores=14 )
##
##
## we will use ranks to find similar sized control alignments
## use ranks
ga.int.l <- orth.l[,'gasterosteus.aculeatus']
ga.int.l[ is.na(ga.int.l) ] <- 0
ga.int.l.o <- order( ga.int.l )
ga.int.l.r <- rank( ga.int.l, ties.method='random' )
ga.int.l.or <- (1:length(ga.int.l))[ga.int.l.o]
## 
##

## have to redo these and all the subsequent steps..
test.ctl.al <- mclapply( test.i[ test.o[1:500] ], align.introns.ctl, max.l=1e8,
                        gap=c(-10, -2), al.sp='gasterosteus.aculeatus',
                        mc.cores=14 )
ctl.ctl.al <- mclapply( ctl.i[ 1:500 ], align.introns.ctl, max.l=1e8,
                        gap=c(-10, -2), al.sp='gasterosteus.aculeatus',
                        mc.cores=14 )
##
## mimic the steps in ../R_intron_alignments_2
## to visualise the results
##
sp.class <- get.sp.taxonomy()
##
## the names are not that appropriate as I did not sampling for
## the controls, but simply took the first 500 
aligns <- list('long'=test.al, 'long.ctl'=test.ctl.al,
               'sampled'=ctl.al, 'sampled.ctl'=ctl.ctl.al)
##
aligns.scores <- lapply(aligns, extract.scores.l, max.l=max.l)
##

aligns.scores.cl <- lapply(colnames(sp.class), function(cl){
    sp.b <- sp.class[,cl]
    names(sp.b) <- rownames(sp.class)
    lapply(aligns, extract.scores.l, max.l=max.l, sp.b=sp.b)
})
names(aligns.scores.cl) <- colnames(sp.class)

##
### let us make a different set of quantiles that we can slice
### the distributions up by.
##
qnt.probs <- seq(0, 1, 0.05)
##
all.qnt <- quantile( do.call(c, lapply(aligns.scores, function(x){ x[,'mn'] })), probs=qnt.probs )
long.qnt <- quantile( do.call(c, lapply(aligns.scores[1:2], function(x){ x[,'mn'] })), probs=qnt.probs )
sampled.qnt <- quantile( do.call(c, lapply(aligns.scores[3:4], function(x){ x[,'mn'] })), probs=qnt.probs )
##
qnt <- list('all'=all.qnt, 'long'=long.qnt, 'sampled'=sampled.qnt)
##

## we can then use these to make sliced histograms...
aligns.scores.mnq <- lapply(aligns.scores, function(x){
    lapply(qnt, function(q){ mnq.class(x, q) }) })
##
aligns.scores.cl.mnq <- lapply(aligns.scores.cl, function(al){
    lapply(al, function(x){
            lapply(qnt, function(q){ mnq.class(x, q) }) }) })

##
### can we use the same quantiles for both the long and sampled
### alignments?
##
sapply( aligns.scores.mnq, function(x){
    table(c(2:length(all.qnt) - 1, x$all)) - 1 })
##
##
## that suggests that the numbers are large enough that we can use
## all.qnt to divide up the distributions. Smallest ones are actually
## the short alignments from the long teleost introns.. 
##

## to define a good set of breaks let us do:
aligns.scores.h <- hist( log2(do.call( c, lapply( aligns.scores, function(x){ x[,'score'] }))),
                        breaks=20, plot=FALSE )
##
aligns.scores.qs1.h <- lapply( 1:length(aligns.scores), function(i){
    sliced.h( log2(aligns.scores[[i]][,'score']), slices=aligns.scores.mnq[[i]]$all,
             breaks=aligns.scores.h$breaks )})
names(aligns.scores.qs1.h) <- names(aligns.scores)
##
##
aligns.scores.cl.qs1.h <- lapply( 1:length(aligns.scores.cl), function(i){
    lapply( 1:length(aligns.scores.cl[[i]]), function(j){
        tmp <- sliced.h( log2(aligns.scores.cl[[i]][[j]][,'score']),
                        slices=aligns.scores.cl.mnq[[i]][[j]]$all,
                        breaks=aligns.scores.h$breaks )
        names(tmp) <- names(aligns.scores.cl[[i]])
        tmp})
})
names(aligns.scores.cl.qs1.h) <- names(aligns.scores.cl)


plot.dens <- TRUE
par(mfrow=c(2,2))
plot.counts( aligns.scores.qs1.h[[1]], dens=plot.dens )
plot.counts( aligns.scores.qs1.h[[2]], dens=plot.dens )
plot.counts( aligns.scores.qs1.h[[3]], dens=plot.dens )
plot.counts( aligns.scores.qs1.h[[4]], dens=plot.dens )


cl.labels <- c("Teleostei", "Mammalia", "Sauria")
hist.labels <- c("long", "long ctl", "ctl", "ctl ctl")


cairo_pdf('quantile_sliced_score_distributions_supp.pdf', width=a4.w * 0.6 * pdf.m,
          height=a4.w * 0.6 * pdf.m, onefile=TRUE)
par(oma=c(2.1, 2.1, 2.1, 2.1))
par(mfrow=c(2,2))
for(j in 1:length(aligns.scores.qs1.h)){
        plot.counts( aligns.scores.qs1.h[[j]], dens=plot.dens,
                    xlab='log2 alignment score', ylab='log2 count' )
        with(par(), mtext(LETTERS[j], at=usr[1], cex=mt.cex, line=1))
        with(par(), text(usr[2], usr[4], hist.labels[j], adj=c(1.1,1.5), cex=0.75*mt.cex))
}
mtext("All species", side=1, outer=TRUE, at=0.1, cex=0.75 * mt.cex)
for(i in 1:length(aligns.scores.cl.qs1.h) ){
    for(j in 1:length(aligns.scores.cl.qs1.h[[i]])){
        plot.counts( aligns.scores.cl.qs1.h[[i]][[j]], dens=plot.dens,
                    xlab='log2 alignment score', ylab='log2 count')
        with(par(), mtext(LETTERS[j], at=usr[1], cex=mt.cex, line=1))
        with(par(), text(usr[2], usr[4], hist.labels[j], adj=c(1.1,1.5), cex=0.75*mt.cex))
    }
    mtext(cl.labels[i], side=1, outer=TRUE, at=0.1, cex=0.75 * mt.cex)
##    inpt <- readline("next: ")
}
dev.off()

## ok, then we can make a figure out of this with reasonable axis labels and
## a legend explaining the colour.
## It's also the case that we should try this with one more control set as
## the differences seem almost too big to be true.

p.cols <- plot.counts( aligns.scores.qs1.h$long )

## how about..
cairo_pdf('quantile_sliced_score_distributions.pdf', width=a4.w * 0.8 * pdf.m,
          height=a4.w * 0.5 * pdf.m)
layout(matrix(1:6, nrow=2, byrow=TRUE))
y <- qnt$all
x <- seq(0, 100, 5)
plot(x, y, xlab='quantile', ylab='log2 mn', type='n')
rect( x[-length(x)], y[-length(y)], x[-1], y[-1], col=p.cols, border=NA )
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))
plot.counts( aligns.scores.qs1.h$long, xlab='log2 alignment score', ylab='log2 count' )
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'long', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.qs1.h$sampled, xlab='log2 alignment score', ylab='log2 count' )
with(par(), mtext('C', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'control', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.cl.qs1.h$tel[[1]], xlab='log2 alignment score', ylab='log2 count' )
with(par(), mtext('D', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'Teleostei', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.cl.qs1.h$mam[[1]], xlab='log2 alignment score', ylab='log2 count' )
with(par(), mtext('E', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'Mammalia', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.cl.qs1.h$sau[[1]], xlab='log2 alignment score', ylab='log2 count' )
with(par(), mtext('F', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'Sauria', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
dev.off()


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


for(i in 1:length(aligns$long)){
    alignment.figure( aligns$long[[i]] )
    inpt <- readline(paste(i, ": "))
    if(inpt == 'q')
        break
}


for(i in 1:length(aligns$sampled)){
    alignment.figure( aligns$sampled[[i]] )
    inpt <- readline(paste(i, ": "))
    if(inpt == 'q')
        break
}

### And there is a definite difference between these two. There are more
### teleost alignments for the gasterosteus than for danio, but that is
### not that strange since danio is a bit of a weird fish.

### In many ways summarising is better done by looking at the best alignment
### per alignment (i.e. the single best scoring hit). This may miss specific
### regions that are aligned, but for the purpose of statistical testing it is
### easier.
### We can then look at the best alignment per species class as well, and
### that can then be summarised in a simpler manner.

aligns.top.scores.cl <- lapply( colnames(sp.class), function(cl){
    sp.b <- sp.class[,cl]
    names(sp.b) <- rownames(sp.class)
    lapply( aligns, extract.top.scores.l, max.l=max.l, sp.b=sp.b )
})
names(aligns.top.scores.cl) <- colnames(sp.class)

plot.best.scores <- function(al.scores, best.only=FALSE, ...){
    all.s <- do.call(c, lapply( al.scores, function(x){ x[,'score'] }) )
    all.i <- do.call(c, lapply( al.scores, function(x){ x[,'i'] }) )
    ylim <- range(log2(all.s))
    plot(1, 1, xlab='rank', ylab='log2 score', xlim=range(all.i), ylim=ylim)
    all.y <- list()
    for(i in 1:length(al.scores)){
        y <- log2(al.scores[[i]][,'score'])
        x <- al.scores[[i]][,'i']
        if(best.only){
            y <- tapply( y, x, max )
            x <- as.numeric(names(y))
        }
        all.y <- c(all.y, list(y))
        points( x, y, col=i, ...)
    }
    invisible(all.y)
}

plot.best.scores( aligns.top.scores.cl[[1]][c(1,3)], cex=0.4 )
plot.best.scores( aligns.top.scores.cl[[2]][c(1,3)], cex=0.4 )
plot.best.scores( aligns.top.scores.cl[[3]][c(1,3)], cex=0.4 )

plot.best.scores( aligns.top.scores.cl[[1]][c(1,3)], cex=1, best.only=TRUE )
plot.best.scores( aligns.top.scores.cl[[2]][c(1,3)], cex=1, best.only=TRUE )
plot.best.scores( aligns.top.scores.cl[[3]][c(1,3)], cex=1, best.only=TRUE )

plot.best.scores( aligns.top.scores.cl[[1]][c(1,2)], cex=1, best.only=TRUE )
plot.best.scores( aligns.top.scores.cl[[1]][c(3,4)], cex=1, best.only=TRUE )

tel.y1 <- plot.best.scores( aligns.top.scores.cl[[1]][c(1,2)], cex=1, best.only=TRUE )
tel.y2 <- plot.best.scores( aligns.top.scores.cl[[1]][c(3,4)], cex=1, best.only=TRUE )
mam.y1 <- plot.best.scores( aligns.top.scores.cl[[2]][c(1,2)], cex=1, best.only=TRUE )
mam.y2 <- plot.best.scores( aligns.top.scores.cl[[2]][c(3,4)], cex=1, best.only=TRUE )
sau.y1 <- plot.best.scores( aligns.top.scores.cl[[3]][c(1,2)], cex=1, best.only=TRUE )
sau.y2 <- plot.best.scores( aligns.top.scores.cl[[3]][c(3,4)], cex=1, best.only=TRUE )

plot(sort( tel.y1[[1]] - tel.y1[[2]] ))
points(sort( tel.y2[[1]] - tel.y2[[2]]), col='red')
abline(h=0, lty=2)

plot(sort( mam.y1[[1]] - mam.y1[[2]] ))
points(sort( mam.y2[[1]] - mam.y2[[2]]), col='red')
abline(h=0, lty=2)

plot(sort( sau.y1[[1]] - sau.y1[[2]] ))
points(sort( sau.y2[[1]] - sau.y2[[2]]), col='red')

## we can try the same but with best.only=FALSE
## but that becomes more messy, so let us only use the above
