## This is a cleaned up analysis of the alignments produced in
## in ../R_intron_alignments. It is primarily rewritten simply in
## order to be more readable. We may also extend this to include
## more alignments. In particular I really would like to do reverse
## complemented alignments too. But I need to make some figures
## first.

## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2


## I may need these if I want to run more alignments.
library(parallel)
dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
source("~/R/general_functions.R")
source("functions.R")

## a simple species classification
sp.class <- get.sp.taxonomy()

## read in the alignment data and assign to a named list.

al.top.500 <- readRDS("../R_intron_alignments/al_top_500.rds")
al.top.500.ctl <- readRDS("../R_intron_alignments/al_top_500_ctl.rds")
al.ctl.500 <- readRDS("../R_intron_alignments/al_ctl_500.rds")
al.ctl.500.ctl <- readRDS("../R_intron_alignments/al_ctl_500_ctl.rds")

## each of these data structures is a list where each element contains
## a set of orthologous introns aligned to, either the danio rerio
## orthologous intron, or a non-orthologous danio rerio intron of
## a similar size (the .ctl set).
##

## al.top.500 contains alignments for introns whose minimal teleost
## length is 1024 and which are represented by at least 10 species
## ordered by the teleost intron length variance.
##
## al.top.500 contain alignments for the same set of orthologous
## introns, but where the danio rerio intron is non-orthologous

## al.ctl.500 contains alignments of a random selection of introns
..## whose length in danio rerio is longer than 1024 but which have a
## median teleost intron length of less than 256 bp.
##
## al.ctl.500.ctl contains the respective controls for these.

### assign to a named list for further analyses
aligns <- list('long'=al.top.500, 'long.ctl'=al.top.500.ctl,
               'sampled'=al.ctl.500, 'sampled.ctl'=al.ctl.500.ctl)
max.l <- 1e8

aligns.scores <- lapply(aligns, extract.scores.l, max.l=max.l)

aligns.scores.cl <- lapply(colnames(sp.class), function(cl){
    sp.b <- sp.class[,cl]
    names(sp.b) <- rownames(sp.class)
    lapply(aligns, extract.scores.l, max.l=max.l, sp.b=sp.b)
})
names(aligns.scores.cl) <- colnames(sp.class)

### let us make a different set of quantiles that we can slice
### the distributions up by.

qnt.probs <- seq(0, 1, 0.05)

all.qnt <- quantile( do.call(c, lapply(aligns.scores, function(x){ x[,'mn'] })), probs=qnt.probs )
long.qnt <- quantile( do.call(c, lapply(aligns.scores[1:2], function(x){ x[,'mn'] })), probs=qnt.probs )
sampled.qnt <- quantile( do.call(c, lapply(aligns.scores[3:4], function(x){ x[,'mn'] })), probs=qnt.probs )

qnt <- list('all'=all.qnt, 'long'=long.qnt, 'sampled'=sampled.qnt)

## we can then use these to make sliced histograms...
aligns.scores.mnq <- lapply(aligns.scores, function(x){
    lapply(qnt, function(q){ mnq.class(x, q) }) })

aligns.scores.cl.mnq <- lapply(aligns.scores.cl, function(al){
    lapply(al, function(x){
            lapply(qnt, function(q){ mnq.class(x, q) }) }) })

### can we use the same quantiles for both the long and sampled
### alignments?

sapply( aligns.scores.mnq, function(x){
    table(c(2:length(all.qnt) - 1, x$all)) - 1 })

##    long long.ctl sampled sampled.ctl
## 1     4        4    5829        5829
## 2     3        3    5828        5830
## 3     3        3    5831        5829
## 4     5        5    5825        5829
## 5     6        6    5828        5824
## 6     9        9    5824        5824
## 7    16       16    5817        5818
## 8    10       10    5820        5825
## 9    31       31    5800        5798
## 10   74       75    5763        5758
## 11  381      380    5452        5452
## 12 1699     1699    4129        4138
## 13 3730     3732    2079        2124
## 14 4740     4739    1115        1072
## 15 5357     5356     481         470
## 16 5568     5574     262         262
## 17 5694     5691     140         140
## 18 5718     5721     113         113
## 19 5752     5751      81          81
## 20 5768     5770      64          64


## that suggests that the numbers are large enough that we can use
## all.qnt to divide up the distributions. Smallest ones are actually
## the short alignments from the long teleost introns.. 

## to define a good set of breaks let us do:
aligns.scores.h <- hist( log2(do.call( c, lapply( aligns.scores, function(x){ x[,'score'] }))),
                        breaks=20 )

aligns.scores.qs1.h <- lapply( 1:length(aligns.scores), function(i){
    sliced.h( log2(aligns.scores[[i]][,'score']), slices=aligns.scores.mnq[[i]]$all,
             breaks=aligns.scores.h$breaks )})
names(aligns.scores.qs1.h) <- names(aligns.scores)

plot.dens <- TRUE
par(mfrow=c(2,2))
plot.counts( aligns.scores.qs1.h[[1]], dens=plot.dens )
plot.counts( aligns.scores.qs1.h[[2]], dens=plot.dens )
plot.counts( aligns.scores.qs1.h[[3]], dens=plot.dens )
plot.counts( aligns.scores.qs1.h[[4]], dens=plot.dens )

aligns.scores.cl.qs1.h <- lapply( 1:length(aligns.scores.cl), function(i){
    lapply( 1:length(aligns.scores.cl[[i]]), function(j){
        tmp <- sliced.h( log2(aligns.scores.cl[[i]][[j]][,'score']),
                        slices=aligns.scores.cl.mnq[[i]][[j]]$all,
                        breaks=aligns.scores.h$breaks )
        names(tmp) <- names(aligns.scores.cl[[i]])
        tmp})
})
names(aligns.scores.cl.qs1.h) <- names(aligns.scores.cl)

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
plot.counts( aligns.scores.qs1.h$long, xlab='log2 alignment score', ylab='log2 count', dens=plot.dens )
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'long', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.qs1.h$sampled, xlab='log2 alignment score', ylab='log2 count', dens=plot.dens )
with(par(), mtext('C', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'control', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.cl.qs1.h$tel[[1]], xlab='log2 alignment score', ylab='log2 count', dens=plot.dens )
with(par(), mtext('D', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'Teleostei', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.cl.qs1.h$mam[[1]], xlab='log2 alignment score', ylab='log2 count', dens=plot.dens )
with(par(), mtext('E', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'Mammalia', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
plot.counts( aligns.scores.cl.qs1.h$sau[[1]], xlab='log2 alignment score', ylab='log2 count', dens=plot.dens )
with(par(), mtext('F', at=usr[1], cex=mt.cex, line=1))
with(par(), text(usr[2], usr[4], 'Sauria', cex=mt.cex * 0.75, adj=c(1.1,1.5)))
dev.off()



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

alignment.figure( aligns$long[[1]] )

pdf("alignments_long.pdf", title="Alignments against long teleost introns",
    width=a4.h * pdf.m * 0.9, height=a4.w * pdf.m * 0.9 )
for(i in 1:100){
    alignment.figure( aligns$long[[i]] )
}
dev.off()

pdf("alignments_sampled.pdf", title="Alignments against sampled teleost introns",
    width=a4.h * pdf.m * 0.9, height=a4.w * pdf.m * 0.9 )
for(i in 1:100){
    alignment.figure( aligns$sampled[[i]] )
}
dev.off()


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

par(mfrow=c(2,1))
plot.best.scores( aligns.top.scores.cl[[1]][c(1,2)], cex=1, best.only=TRUE )
plot.best.scores( aligns.top.scores.cl[[1]][c(3,4)], cex=1, best.only=TRUE )

plot.best.scores( aligns.top.scores.cl[[1]][c(1,3)], cex=1, best.only=TRUE )
plot.best.scores( aligns.top.scores.cl[[1]][c(2,4)], cex=1, best.only=TRUE )


tel.y1 <- plot.best.scores( aligns.top.scores.cl[[1]][c(1,2)], cex=1, best.only=TRUE )
tel.y2 <- plot.best.scores( aligns.top.scores.cl[[1]][c(3,4)], cex=1, best.only=TRUE )
mam.y1 <- plot.best.scores( aligns.top.scores.cl[[2]][c(1,2)], cex=1, best.only=TRUE )
mam.y2 <- plot.best.scores( aligns.top.scores.cl[[2]][c(3,4)], cex=1, best.only=TRUE )
sau.y1 <- plot.best.scores( aligns.top.scores.cl[[3]][c(1,2)], cex=1, best.only=TRUE )
sau.y2 <- plot.best.scores( aligns.top.scores.cl[[3]][c(3,4)], cex=1, best.only=TRUE )

plot(sort( tel.y1[[1]] - tel.y1[[2]] ), xlim=c(1, length(tel.y2[[1]])))
points(sort( tel.y2[[1]] - tel.y2[[2]]), col='red')
abline(h=0, lty=2)

plot(sort( mam.y1[[1]] - mam.y1[[2]] ), xlim=c(1, length(mam.y2[[1]])))
points(sort( mam.y2[[1]] - mam.y2[[2]]), col='red')
abline(h=0, lty=2)

plot(sort( sau.y1[[1]] - sau.y1[[2]] ), xlim=c(1, length(sau.y2[[1]])))
points(sort( sau.y2[[1]] - sau.y2[[2]]), col='red')

## we can try the same but with best.only=FALSE
## but that becomes more messy, so let us only use the above

### This might be a reasonable idea
par(mfrow=c(1,3))
with( aligns.top.scores.cl[[1]], plot(log2(long[,'l1']), long[,'score'], cex=0.5))
with( aligns.top.scores.cl[[1]], points(log2(sampled[,'l1']), sampled[,'score'], cex=0.5, col='red'))
with( aligns.top.scores.cl[[2]], plot(log2(long[,'l1']), long[,'score'], cex=0.5))
with( aligns.top.scores.cl[[2]], points(log2(sampled[,'l1']), sampled[,'score'], cex=0.5, col='red'))
with( aligns.top.scores.cl[[3]], plot(log2(long[,'l1']), long[,'score'], cex=0.5))
with( aligns.top.scores.cl[[3]], points(log2(sampled[,'l1']), sampled[,'score'], cex=0.5, col='red'))

## or:
par(mfrow=c(1,3))
with( aligns.top.scores.cl[[1]], plot(long[,'mn'], long[,'score'], cex=0.5))
with( aligns.top.scores.cl[[1]], points(sampled[,'mn'], sampled[,'score'], cex=0.5, col='red'))
with( aligns.top.scores.cl[[2]], plot(long[,'mn'], long[,'score'], cex=0.5))
with( aligns.top.scores.cl[[2]], points(sampled[,'mn'], sampled[,'score'], cex=0.5, col='red'))
with( aligns.top.scores.cl[[3]], plot(long[,'mn'], long[,'score'], cex=0.5))
with( aligns.top.scores.cl[[3]], points(sampled[,'mn'], sampled[,'score'], cex=0.5, col='red'))

## but it seems that we actually need the distributions of these for different
## windows of l1 or mn. We can maybe do this using discretize2d as as in figure
## 2, and scale the columns...


## this needs to be plotted with some sort of heatmap..
require(entropy)
hist.2d <- function(x1, x2, numBins=20){
    b <- !is.na(x1) & !is.na(x2) & is.finite(x1) & is.finite(x2)
    h <- discretize2d( x1[b], x2[b], numBins1=numBins, numBins2=numBins )
    ## try to get the breaks
    rn <- sub("[]]", "", rownames(h))
    cn <- sub("[]]", "", colnames(h))
    rn <- sub("[[)(]", "", rn)
    cn <- sub("[[)(]", "", cn)
    row.breaks <- unique(as.numeric(unlist(strsplit(rn, ','))))
    col.breaks <- unique(as.numeric(unlist(strsplit(cn, ','))))
    list('b'=b, 'h'=h, 'rb'=row.breaks, 'cb'=col.breaks)
}

v2col <- function(v, cm=0.3){
    vv= v - min(v, na.rm=TRUE)
    val <- (cm + vv)/(cm + max(vv, na.rm=TRUE))
    hsvScale(v, val=val)
}


## OK. We will make a plot with these types of values:
## (cm is a colour moderator that moderates the colour value)
## x and y are the breaks of the rectangles, (x[-1] == x2)
## (inner) margins are bottom, left, top, right
rect_image <- function(m, x=0:ncol(m), y=0:nrow(m),
                       col.f=v2col, col.args=NULL, margins=rep(0, 4),
                       border=NA, clear=TRUE, axis=TRUE,
                       xlab=NULL, ylab=NULL, lab.cex=1.5, cm=0.3){
    x1 <- x[-length(x)]
    x2 <- x[-1]
    y1 <- y[-length(y)]
    y2 <- y[-1]
    if(clear)
        plot.new()
    xr <- range(x)
    yr <- range(y)
    xr <- xr + diff(xr) * c(-margins[2], margins[4])
    yr <- yr + diff(yr) * c(-margins[1], margins[3])
    plot.window(xlim=xr, ylim=yr, xaxs='i', yaxs='i')
    ## and then make x1, x2, y1, and y2 the appropriate matrices
    x1 <- matrix(rep(x1, nrow(m)), nrow=nrow(m), byrow=TRUE)
    x2 <- matrix(rep(x2, nrow(m)), nrow=nrow(m), byrow=TRUE)
    y1 <- matrix(rep(y1, ncol(m)), nrow=nrow(m))
    y2 <- matrix(rep(y2, ncol(m)), nrow=nrow(m))
    if(is.null(col.args))
        cols <- col.f(m, cm)
    else
        cols <- do.call(col.f, c(list(v=m), col.args))
    rect(x1, y1, x2, y2, col=cols, border=border)
    if(axis){
        axis(1)
        axis(2)
    }
    if(!is.null(xlab))
        mtext(xlab, side=1, line=2.5, cex=lab.cex)
    if(!is.null(ylab))
        mtext(ylab, side=2, line=2.5, cex=lab.cex)
}


top.l1.scores.2d <- lapply( aligns.top.scores.cl, function(x){
    lapply( x, function(y){
        hist.2d( y[,'score'], log2(y[,'l1']) )
    })
})

top.mn.scores.2d <- lapply( aligns.top.scores.cl, function(x){
    lapply( x, function(y){
        hist.2d( y[,'score'], y[,'mn'] )
    })
})

top.mn.scores.lx.2d <- lapply( aligns.top.scores.cl, function(x){
    lapply( x, function(y){
        hist.2d( y[,'score'], log2(y[,'mn']) )
    })
})

top.mn.scores.ll.2d <- lapply( aligns.top.scores.cl, function(x){
    lapply( x, function(y){
        hist.2d( log2(y[,'score']), log2(y[,'mn']) )
    })
})


top.i.scores.2d <- lapply( aligns.top.scores.cl, function(x){
    lapply( x, function(y){
        hist.2d( log2(y[,'score']), y[,'i'] )
    })
})


par(mfrow=c(2,3))
with( top.l1.scores.2d[[1]]$long, rect_image( log2(h+1), x=cb, y=rb ))
with( top.l1.scores.2d[[2]]$long, rect_image( log2(h+1), x=cb, y=rb ))
with( top.l1.scores.2d[[3]]$long, rect_image( log2(h+1), x=cb, y=rb ))

with( top.l1.scores.2d[[1]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.l1.scores.2d[[2]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.l1.scores.2d[[3]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))

par(mfrow=c(2,3))
with( top.mn.scores.2d[[1]]$long, rect_image( log2(h+1), x=cb, y=rb ))
with( top.mn.scores.2d[[2]]$long, rect_image( log2(h+1), x=cb, y=rb ))
with( top.mn.scores.2d[[3]]$long, rect_image( log2(h+1), x=cb, y=rb ))

with( top.mn.scores.2d[[1]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.mn.scores.2d[[2]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.mn.scores.2d[[3]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))

par(mfrow=c(2,3))
with( top.mn.scores.ll.2d[[1]]$long, rect_image( scale(h), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[2]]$long, rect_image( scale(h), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[3]]$long, rect_image( scale(h), x=cb, y=rb ))

with( top.mn.scores.ll.2d[[1]]$sampled, rect_image( scale(h), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[2]]$sampled, rect_image( scale(h), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[3]]$sampled, rect_image( scale(h), x=cb, y=rb ))

par(mfrow=c(2,3))
with( top.mn.scores.ll.2d[[1]]$long, rect_image( log2(h + 1), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[2]]$long, rect_image( log2(h + 1), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[3]]$long, rect_image( log2(h + 1), x=cb, y=rb ))

with( top.mn.scores.ll.2d[[1]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[2]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.mn.scores.ll.2d[[3]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))



par(mfrow=c(2,3))
with( top.i.scores.2d[[1]]$long, rect_image( log2(h+1), x=cb, y=rb ))
with( top.i.scores.2d[[2]]$long, rect_image( log2(h+1), x=cb, y=rb ))
with( top.i.scores.2d[[3]]$long, rect_image( log2(h+1), x=cb, y=rb ))

with( top.i.scores.2d[[1]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.i.scores.2d[[2]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))
with( top.i.scores.2d[[3]]$sampled, rect_image( log2(h+1), x=cb, y=rb ))

par(mfrow=c(2,3))
with( top.i.scores.2d[[1]]$long, rect_image( h, x=cb, y=rb ))
with( top.i.scores.2d[[2]]$long, rect_image( h, x=cb, y=rb ))
with( top.i.scores.2d[[3]]$long, rect_image( h, x=cb, y=rb ))

with( top.i.scores.2d[[1]]$sampled, rect_image( h, x=cb, y=rb ))
with( top.i.scores.2d[[2]]$sampled, rect_image( h, x=cb, y=rb ))
with( top.i.scores.2d[[3]]$sampled, rect_image( h, x=cb, y=rb ))

## take the top scoring alignment for each species for each clade..
top.top <- function(m){
    tt <- tapply( 1:nrow(m), m[,'i'], function(i){
        j <- which.max( m[i,'score'] )
        m[i[j], c('i', 'mn', 'score')]
    })
    do.call(rbind, tt)
}

aligns.top.top.cl <- lapply(aligns.top.scores.cl, function(x){
    lapply(x, top.top)})

aligns.top.top.cl.mod <- lapply(aligns.top.top.cl, function(x){
    lapply(x, function(y){
        list('lin'=lm( y[,'score'] ~ y[,'mn'] ),
             'log'=lm( log2(y[,'score']) ~ y[,'mn'] ))
    })
})


aligns.top.mod.2 <- lapply( aligns.top.scores.cl, function(x){
    lapply(x, function(y){
        pars <- tapply( 1:nrow(y), y[,'i'], function(i){
            c('l1'=as.double(y[i[1],'l1']), 'l2'=as.double(sum( y[i,'l2'] )),
              'score'=as.double(max(y[i,'score'])))
        })
        pars <- sapply(pars, eval)
        mod <- lm( log2(pars['score',]) ~ log2( pars['l1',] * pars['l2',] ) )
        list('p'=pars, 'lm'=mod)
    })
})
                                          

par(mfrow=c(2,2))
with(aligns.top.mod.2$tel$long.ctl, plot(log2(p['l1',] * p['l2',]), log2(p['score',])))
with(aligns.top.mod.2$tel$sampled.ctl, points(log2(p['l1',] * p['l2',]), log2(p['score',]), col='red'))
with(aligns.top.mod.2$tel$long.ctl, abline(lm, col='red', lwd=2))
with(aligns.top.mod.2$tel$sampled.ctl, abline(lm, col='purple', lwd=2))

with(aligns.top.mod.2$tel$long, plot(log2(p['l1',] * p['l2',]), log2(p['score',])))
with(aligns.top.mod.2$tel$sampled, points(log2(p['l1',] * p['l2',]), log2(p['score',]), col='red'))
with(aligns.top.mod.2$tel$long.ctl, abline(lm, col='red', lwd=2))
with(aligns.top.mod.2$tel$sampled.ctl, abline(lm, col='purple', lwd=2))

## now the model for the sampled.ctl and the long.ctl look very similar.
## giving more evidence that these are the models that we should use.

par('mfrow'=c(2,3))
for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(long[,'mn'], log2(long[,'score'])) )
    abline( aligns.top.top.cl.mod[[i]]$long$log, col='red' )
}

for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(sampled[,'mn'], log2(sampled[,'score'])) )
    abline( aligns.top.top.cl.mod[[i]]$sampled$log, col='red' )
}

for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(long.ctl[,'mn'], log2(long.ctl[,'score'])) )
    abline( aligns.top.top.cl.mod[[i]]$long.ctl$log, col='red' )
}

for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(sampled.ctl[,'mn'], log2(sampled.ctl[,'score'])) )
    abline( aligns.top.top.cl.mod[[i]]$sampled.ctl$log, col='red' )
}

par('mfrow'=c(2,3))
for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(long.ctl[,'mn'], (long.ctl[,'score']), col='red') )
    with(aligns.top.top.cl[[i]], points(sampled.ctl[,'mn'], (sampled.ctl[,'score']), col='blue') )
    abline( aligns.top.top.cl.mod[[i]]$long.ctl$lin, col='red', lwd=3 )
    abline( aligns.top.top.cl.mod[[i]]$sampled.ctl$lin, col='blue', lwd=3 )
}

for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(long.ctl[,'mn'], log2(long.ctl[,'score']), col='red') )
    with(aligns.top.top.cl[[i]], points(sampled.ctl[,'mn'], log2(sampled.ctl[,'score']), col='blue') )
    abline( aligns.top.top.cl.mod[[i]]$long.ctl$log, col='red', lwd=3 )
    abline( aligns.top.top.cl.mod[[i]]$sampled.ctl$log, col='blue', lwd=3 )
}

par('mfrow'=c(2,3))
for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(long[,'mn'], (long[,'score']), col='red') )
    with(aligns.top.top.cl[[i]], points(sampled[,'mn'], (sampled[,'score']), col='blue') )
    abline( aligns.top.top.cl.mod[[i]]$long.ctl$lin, col='red', lwd=3 )
    abline( aligns.top.top.cl.mod[[i]]$sampled.ctl$lin, col='blue', lwd=3 )
}

for(i in 1:length(aligns.top.top.cl)){
    with(aligns.top.top.cl[[i]], plot(long[,'mn'], log2(long[,'score']), col='red') )
    with(aligns.top.top.cl[[i]], points(sampled[,'mn'], log2(sampled[,'score']), col='blue') )
    abline( aligns.top.top.cl.mod[[i]]$long.ctl$log, col='red', lwd=3 )
    abline( aligns.top.top.cl.mod[[i]]$sampled.ctl$log, col='blue', lwd=3 )
}

## this is unnecessary; aligns.top.top.cl is a single alignment per sequence already
plot.top <- function(x, pts, mod, mods, cols=c('red', 'blue'), mod.t='lin', tform=eval){
    xy <- lapply(pts, function(smp){
        t(sapply( tapply( 1:nrow( x[[smp]] ), x[[smp]][,'i'], function(i){
            j <- which.max(x[[smp]][i,'score'])
            c(x[[smp]][i[j],c('i', 'mn', 'score')])
        }), eval))
    })
    all.x <- unlist(sapply(xy, function(x){ x[,'mn'] }))
    all.y <- unlist(sapply(xy, function(x){ tform(x[,'score']) }))
    plot(all.x, all.y, type='n', xlab='log2 mn', ylab='score')
    for(i in 1:length(xy)){
        points( xy[[i]][,'mn'], tform(xy[[i]][,'score']), col=cols[1 + (i-1) %% length(cols) ] )
    }
    for(i in 1:length(mods))
        abline( mod[[mods[i]]][[mod.t]], col=cols[1 + (i-1) %% length(cols) ], lwd=2 )
    invisible(xy)
}

xy <- plot.top( aligns.top.top.cl[[1]], pts=c('long', 'sampled'),
               mod=aligns.top.top.cl.mod[[1]], mods=c('long.ctl', 'sampled.ctl'),
               tform=log2, mod.t='log')

xy <- plot.top( aligns.top.top.cl[[1]], pts=c('long', 'sampled'),
               mod=aligns.top.top.cl.mod[[1]], mods=c('long.ctl', 'sampled.ctl'),
               tform=eval, mod.t='lin')

par('mfrow'=c(2,3))
for(i in 1:length(aligns.top.top.cl)){
    plot.top( aligns.top.top.cl[[i]], pts=c('long', 'sampled'),
             mod=aligns.top.top.cl.mod[[i]], mods=c('long.ctl', 'sampled.ctl'),
             tform=eval, mod.t='lin' )
}

for(i in 1:length(aligns.top.top.cl)){
    plot.top( aligns.top.top.cl[[i]], pts=c('long', 'sampled'),
             mod=aligns.top.top.cl.mod[[i]], mods=c('long.ctl', 'sampled.ctl'),
             tform=log2, mod.t='log' )
}

for(i in 1:length(aligns.top.top.cl)){
    plot.top( aligns.top.top.cl[[i]], pts=c('long.ctl', 'sampled.ctl'),
             mod=aligns.top.top.cl.mod[[i]], mods=c('long.ctl', 'sampled.ctl'),
             tform=eval, mod.t='lin' )
}

for(i in 1:length(aligns.top.top.cl)){
    plot.top( aligns.top.top.cl[[i]], pts=c('long.ctl', 'sampled.ctl'),
             mod=aligns.top.top.cl.mod[[i]], mods=c('long.ctl', 'sampled.ctl'),
             tform=log2, mod.t='log' )
}

for(i in 1:length(aligns.top.top.cl)){
    plot.top( aligns.top.top.cl[[i]], pts=c('long', 'sampled'),
             mod=aligns.top.top.cl.mod[[i]], mods=c('long.ctl', 'sampled.ctl'),
             tform=eval, mod.t='lin' )
}

### The above functions are a bit unnecessary as we already have a single alignment
### per seed sequence.

## What are the high scoring control alignments?
o <- order( aligns.top.top.cl[[1]]$long.ctl[,'score'], decreasing=TRUE)
head(aligns.top.top.cl[[1]]$long.ctl[o,] )
##       i       mn score
## 418 418 26.24893  1376
## 387 387 26.23373   846
## 97   97 26.27064   726
## 160 160 25.43801   722
## 91   91 26.21902   714
## 92   92 24.54633   680

with( aligns.top.scores.cl$tel, head( long.ctl[ long.ctl[,'i'] == 418, ] ))
##        i    l1   l2 score       mn
## [1,] 418 22962 3957   114 26.43715
## [2,] 418 22962 3473  1376 26.24893
## [3,] 418 22962 1933   134 25.40359
## [4,] 418 22962 3946   140 26.43314
## [5,] 418 22962 4273   110 26.54799
## [6,] 418 22962 1515    88 25.05206
## the second one of these has the high scores, but that is the only one it
## seems. Very strange indeed. It seems that instead of taking top.top, we might
## be better with top.median. But lets try to have a look at that..

##
aligns$long.ctl[[418]]$l1  ## 22962
## and
which((aligns$long.ctl[[418]]$l2) == 3473)
## [1] 33

## so we can have a look at the alignment
head(aligns$long.ctl[[418]]$al[[3]]$pos)

##      a_beg a_end b_beg b_end score length
## [1,]  3049  3839  1902  2836  1376    963
## [2,]   232   429   974  1163    84    216
## so we have found it and can do..

align.print( aligns$long.ctl[[418]]$al[[33]]$seq[[1]], w=80 )
## and that really looks like a reasonable motif.. What are the species here?
names( aligns$long.ctl[[418]]$al )[33]
## [1] "clupea.harengus"
## that is atlantic herring ??
## This sequence would seem to be pretty common across the teleosts; and similar
## sequenes can be found in many places in both danio rerio and clupea harengus
## This is clearly a repeat. How about the second best sequence that we have
## with a score of 846.

with( aligns.top.scores.cl$tel, head( long.ctl[ long.ctl[,'i'] == 387, ]))
##        i   l1    l2 score       mn
## [1,] 387 4840 16304   846 26.23373
## [2,] 387 4840 12665   102 25.86935
## [3,] 387 4840 15037   758 26.11702

which((aligns$long.ctl[[387]]$l2) == 16304)
## 1
head( aligns$long.ctl[[387]]$al[[1]]$pos ) ## that is correct
align.print( aligns$long.ctl[[387]]$al[[1]]$seq[[1]], w=80 )

## to copy and paste the sequences to ncbi blast
gsub("-", "", aligns$long.ctl[[387]]$al[[1]]$seq[[1]][1])
gsub("-", "", aligns$long.ctl[[387]]$al[[1]]$seq[[1]][2])
## again these sequences similar to these are found in many species and many
## locations of the genome.

## What if we look at the sequences which we have got from the long; what do
## these look like.

o <- order( aligns.top.top.cl[[1]]$long[,'score'], decreasing=TRUE)
head(aligns.top.top.cl[[1]]$long[o,] )
##       i       mn score
## 5     5 25.00018  6568
## 139 139 26.29454  5936
## 1     1 23.40884  4270
## 56   56 24.47072  2550
## 94   94 25.61145  2504
## 180 180 25.33034  2462

with( aligns.top.scores.cl$tel, head( long[ long[,'i'] == 5, ]))
##      i   l1   l2 score       mn
## [1,] 5 5871 4954  4622 24.79377
## [2,] 5 5871 4400  4380 24.62268
## [3,] 5 5871 5104  4742 24.83680
## [4,] 5 5871 5094  4698 24.83397

## many high scores indeed.
which( aligns$long[[5]]$l2 == 4954 ) ## 1 ?
head( aligns$long[[5]]$al[[1]]$pos ) ## that is correct

align.print( aligns$long[[5]]$al[[1]]$seq[[1]], w=80 )

gsub("-", "", aligns$long[[5]]$al[[1]]$seq[[1]][1])
gsub("-", "", aligns$long[[5]]$al[[1]]$seq[[1]][2])
## And that is pretty fantastic; looks like a single alignment in a range of
## species. So we have an ultra-conserved element here.

## lets look at one more..
with( aligns.top.scores.cl$tel, head( long[ long[,'i'] == 139, ]))
with( aligns.top.scores.cl$tel, head( long[ long[,'i'] == 139, ], n=7))
##        i   l1    l2 score       mn
## [1,] 139 8548  7653  3308 25.96318
## [2,] 139 8548  7051  3044 25.84498
## [3,] 139 8548  7309  3272 25.89683
## [4,] 139 8548  7334  3226 25.90176
## [5,] 139 8548  7083  3318 25.85152
## [6,] 139 8548  6856  3014 25.80452
## [7,] 139 8548 11367  5406 26.53394

which( aligns$long[[139]]$l2 == 11367 ) ## 1 ?
head( aligns$long[[139]]$al[[13]]$pos ) ## that is correct
align.print( aligns$long[[139]]$al[[13]]$seq[[1]], 80 )

gsub("-", "", aligns$long[[139]]$al[[13]]$seq[[1]][1])
## Unfortunately, it seems that this is an intron which contains an alternative
## exon (which includes it's own promoter; so this is an alternative intron)
## but it does look like the intronic sequences themselves are also aligned
## across multiple species.
## So here we have an example of this. I will clearly have to extract sequences
## and draw the alignments 

## It seems clear from these that I should actually extract the top alignment
## from each class and then run blast of the sequences to determine the nature
## of the alignments.

aligns.top.al.cl <- lapply( colnames(sp.class), function(cl){
    sp.b <- sp.class[,cl]
    names(sp.b) <- rownames(sp.class)
    lapply( aligns, extract.top.al.l, max.l=max.l, sp.b=sp.b )
})
names(aligns.top.al.cl) <- colnames(sp.class)

saveRDS( aligns.top.al.cl, "aligns_top_al_cl.rds" )

## lets make fasta sequences...
for( cl in names(aligns.top.al.cl) ){
    for(sample in names(aligns.top.al.cl[[cl]])){
        fname <- paste("extracted_seqs/", cl, "_", sample, ".fasta", sep="")
        lines <- do.call(c,
                         lapply(aligns.top.al.cl[[cl]][[sample]], function(x){
                             if(!length(x)) return(NULL)
                             c(paste(c(">", orth.id[ x[[1]]$i, x[[1]]$sp1],
                                       "_", orth.i[ x[[1]]$i, x[[1]]$sp1],
                                       " ", orth.tr[ x[[1]]$i, x[[1]]$sp1],
                                       " ", orth.id[ x[[1]]$i, x[[1]]$sp2],
                                       "_", orth.i[ x[[1]]$i, x[[1]]$sp2],
                                       " ", orth.tr[ x[[1]]$i, x[[1]]$sp2],
                                       " ", x[[1]]$sp1, " ", x[[1]]$sp2, " ",
                                       paste(x[[1]]$coords, collapse=" ")), collapse=""),
                               degap.seq( x[[1]]$seq[1] ))
                         }))
        writeLines(lines, fname)
    }
}

## and lets simply extract out the coordinates for the danio rerio orthologues
## Then just use a perl script to extract the complete intron sequences and
## use these a blast queries against the set of databases we have used for the
## the aligned sections.

for( sample in names(aligns) ){
    fname <- paste("extracted_introns/", sample, ".txt", sep="")
    lines <- do.call(c, lapply(aligns[[sample]], function(x){
        i <- x$i
        sp <- x$sp1
        paste(orth.fam[i], orth.id[i, sp], orth.tr[i, sp], orth.i[i, sp], i,
              sep="\t")
    }))
    writeLines(lines, fname)
}
        

## takes a table with columns 'l1', 'l2', 'score'
## to something approaching a log10( e-value 0
tform.s.1 <- function(x, lambda=0.01){
##    log( x[,'l1'] ) + log(x[,'l2']) - lambda * x[,'score']
    -log10(x[,'l1'] * x[,'l2'] * exp(-lambda * x[,'score']))
}

## tforms to something approximating a p-value
tform.s.2 <- function(x, lambda=0.01){
##    log( x[,'l1'] ) + log(x[,'l2']) - lambda * x[,'score']
    E <- (x[,'l1'] * x[,'l2'] * exp(-lambda * x[,'score']))
    1 - exp(-E)
}


plot.normalised.dists <- function(x, breaks=20, density=FALSE, cols=1:length(x),
                                  s.tform=tform.s.1, lambda=0.01, type='l', ...){
    all.v <- do.call(c, lapply(x, function(y){
        s.tform(y, lambda=lambda) }))
    h.all <- hist(all.v, breaks=breaks, plot=FALSE)
    h.s <- lapply(x, function(y){
        hist( s.tform(y, lambda=lambda), plot=FALSE, breaks=h.all$breaks ) })
    if(density)

    else
        v <- sapply( h.s, function(x){ x$counts })
    plot( h.all$mids, h.all$counts, type='n', ylim=range(v), xlab='log( score/mn )',
         ylab=ifelse(density, 'density', 'counts') )
    for(i in 1:ncol(v))
        lines( h.all$mids, v[,i], col=cols[i], type=type, ... )
    legend('topright', legend=names(x), lty=1, col=cols, ...)
    invisible(v)
}

## a different way of looking at it perhaps..
par(mfrow=c(2,2))
with( aligns.top.scores.cl[[1]], hist(log2(long[,'score'] / long[,'mn']) ))
with( aligns.top.scores.cl[[1]], hist(log2(long.ctl[,'score'] / long.ctl[,'mn']) ))
with( aligns.top.scores.cl[[1]], hist(log2(sampled[,'score'] / sampled[,'mn']) ))
with( aligns.top.scores.cl[[1]], hist(log2(sampled.ctl[,'score'] / sampled.ctl[,'mn']) ))

par(mfrow=c(1,3))
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.1, density=TRUE )
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.01, density=TRUE )
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.001, density=TRUE )

par(mfrow=c(1,3))
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.15,
                      density=TRUE, s.tform=tform.s.2, type='b' )
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.1,
                      density=TRUE, s.tform=tform.s.2 )
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.055,
                      density=TRUE, s.tform=tform.s.2 )

plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.1,
                      density=TRUE, s.tform=tform.s.2 )

## to tune lambda, we can ask what proportion of p-values are below 0.95 for a given
## random sample. It seems that we do have to tune both lambda and K in order for
## this to work. After a bit of scribbling it seems that for alignments of unrelated
## sequences we would expect to see one of the given score. That means that we can rearrange
##
## E = Kmn*exp(-lambda * S)
## where
## E is the expected number of alignments, S, is the score, mn = m*n and K is a constant
## if we take the log of both sides we det
##
## log(E) = log(K) + log(mn) - (lambda * S)
## and since we expect one*, that reduces to
## 0 = log(K) + log(mn) - (lambda * S)
##
## that can be rearrange to:
##
## S = log(K)/lambda + log(mn)/lambda
##
## and since the first part is a constant we can express the above as
## S ~ log(mn)
##
## and use lm to solve the equation using the control alignments.
## *this assumption is probably a bit flawed, but we can possibly get some
## reasonable estimate of K and lambda doing the above.

estimate.params <- function(){
    lapply( aligns.top.scores.cl, function(x){
        lapply( x[c('long.ctl', 'sampled.ctl')], function(y){
            lg.nm <- log( y[,'l1'] * y[,'l2'] )
            scores <- y[,'score']
            model <- lm( scores ~ lg.nm )
            list('model'=model, 'lg.nm'=lg.nm, 'scores'=scores)
        })
    })
}

al.params <- estimate.params()
## these vary quite dramatically. Especially with respect to the intercept.
## to get the base parameters
hist( al.params$tel$long.ctl$model$coefficients[2] / al.params$tel$long.ctl$lg.nm )
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.9, density=TRUE )

plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.95,
                      density=TRUE, s.tform=tform.s.2, type='b' )
## so this really does not seem to make any sense. I really do need a better way to
## modify these so that I get similar distributions..
## This is not making much sense to me. I have to put it down for now..


par(mfrow=c(1,3))
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, lambda=0.9, density=TRUE )

plot.normalised.dists( aligns.top.scores.cl[[2]], breaks=30, lwd=2, lambda=0.1, density=TRUE)

plot.normalised.dists( aligns.top.scores.cl[[3]], breaks=30, lwd=2 )

par(mfrow=c(1,3))
plot.normalised.dists( aligns.top.scores.cl[[1]], breaks=30, lwd=2, density=FALSE )
plot.normalised.dists( aligns.top.scores.cl[[2]], breaks=30, lwd=2, density=FALSE )
plot.normalised.dists( aligns.top.scores.cl[[3]], breaks=30, lwd=2, density=FALSE )
