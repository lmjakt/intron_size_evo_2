## Redo analyses in ../R_intron_alignments_summary/
## with the proper set of alignments. Previous alignments
## had been done with incorrectly selected introns due to a
## difference in row and column orders of tables.

## Much of this is copied and then modified from above:

## Refactor the analyses carried out in :
## ../R_intron_alignments_2/extracted_seqs/blast/blast_alignments.R
##
## These make use of alignment data from:
##
## ../R_intron_alignments/
## ../R_intron_alignments_2/
## ../R_intron_alignments_5/
##
## and blast data from:
## ../R_intron_alignments_2/extracted_seqs/blast/
## ../R_intron_alignments_5/extracted_seqs/blast/
##
## orthology tables in:
## ../R_172_genomes/
##
## sp classification in:
## ../R_trees_distances

source('~/R/general_functions.R')
source('../R_common_functions/functions.R') ## for ensembl functions
source("../R_intron_alignments_2/functions.R")
source('functions.R') ## for local functions
require("RMySQL")
require("parallel")

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")

## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

## the max.nm allowed
max.l <- 1e8

### let us read in the orthology tables first
orth <- read.ortho()

saveRDS(orth, 'orth.rds')

## add supplementary including local alignment scores at intron positions
orth.sc <- readRDS("../R_intron_orthology/intron_orth_lscore.rds")
## we actually only want to use the 'danio rerio' part here:
orth.sc <- orth.sc[['danio rerio']]

saveRDS(orth.sc, 'orth_sc_dr.rds')
## modify the col names to fit with the orth tables
for(i in 1:length(orth.sc))
    colnames(orth.sc[[i]]) <- sub(" ", ".",  colnames(orth.sc[[i]]), fixed=TRUE)

all( orth.sc$i == orth$i, na.rm=TRUE )  ## TRUE..


## and the database credentials: (which have an unset password: set this manually, when needed)
db.cred <- read.credentials("../R_intron_alignments_2/extracted_seqs/blast/db.cred")


## Some species classification:
sp.class <- readRDS( "../R_trees_distances/sp_class_3.rds")
class.col <-  readRDS( "../R_trees_distances/class_col_3.rds")
sp.col <- readRDS("../R_trees_distances/sp_col_3.rds")
names(sp.col) <- sub("_", ".", names(sp.col))
rownames(sp.class) <- sub("_", ".", rownames(sp.class))

## we will at some point want to make use of the tel.var an so on.
var.par <- function(x){
    c('n'=sum(!is.na(x)), 'mean'=mean(x, na.rm=TRUE), 'sd'=sd(x, na.rm=TRUE),
      'min'=min(x, na.rm=TRUE), 'max'=max(x, na.rm=TRUE), quantile(x, probs=seq(0,1,0.1), na.rm=TRUE))
}

int.s.b <- readRDS("../R_trees_distances/int_s_b.rds")
names(int.s.b) <- sub("_", ".", names(int.s.b))
all(names(int.s.b) == rownames(sp.class)) ## TRUE

tel.var <- t(apply( log2(orth$l[, rownames(sp.class)[ sp.class[,'teleostei'] & int.s.b ] ]),
                    1, var.par ))

mam.var <- t(apply( log2(orth$l[, rownames(sp.class)[ sp.class[,'mammalia'] & int.s.b ] ]),
                    1, var.par ))

saveRDS(tel.var, 'tel_var.rds')
saveRDS(mam.var, 'mam_var.rds')

## The full set of alignments:
aligns <- list(
    'long' = readRDS("../R_intron_alignments_6/long_aligns.rds"),
    'med' = readRDS("../R_intron_alignments_6/med_aligns.rds"),
    'ctl' = readRDS("../R_intron_alignments_6/ctl_aligns.rds")
)


## ctl alignments. These should be complemented
## with more later on. But for now we can at least make a model
## from these:

aligns.ctl <- list(
    'long' = readRDS("../R_intron_alignments_6/long_ctl_aligns.rds"),
    'med' = readRDS("../R_intron_alignments_6/med_ctl_aligns.rds"),
    'ctl' = readRDS("../R_intron_alignments_6/ctl_ctl_aligns.rds")
)

## supplementary aligns; these are not control alignments, but more similar
## the aligns$ctl above. We can call them short...

aligns.short <- readRDS("../R_intron_alignments_8/ctl_2_aligns.rds")
aligns.short.2 <- readRDS("../R_intron_alignments_8/ctl_3_aligns.rds")

## lets output plots for all of these alignments:
make.alignment.figures <- function(al.l, suffix){
    for(i in 1:length(al.l)){
        fname <- paste( names(al.l)[i], "_", suffix[i], ".pdf", sep="")
        cairo_pdf( fname, onefile=TRUE, width=a4.w * pdf.m, height=a4.h * pdf.m )
        par(omi=c(0.3, 0.3, 0.3, 0.3))
        for(j in 1:length(al.l[[i]])){
            alignment.figure( al.l[[i]][[j]] )
        }
        dev.off()
    }
}

make.alignment.figures( c(aligns, 'short_1'=list(aligns.short), 'short_2'=list(aligns.short.2), aligns.ctl),
                        c(rep('test', 5), rep('ctl', 3) ))

## set the rownumbers of aligns so that we can access the raw information
## easily
for(sm in names(aligns)){
    names(aligns[[sm]]) <- as.character( sapply( aligns[[sm]], function(x){ x$i } ))
}

for(sm in names(aligns.ctl)){
    names(aligns.ctl[[sm]]) <- as.character( sapply( aligns.ctl[[sm]], function(x){ x$i } ))
}

names(aligns.short) <- as.character( sapply( aligns.short, function(x){ x$i }) )
names(aligns.short.2) <- as.character( sapply( aligns.short.2, function(x){ x$i }) )

## create new summary structures with a simpler NULL values for missing alignments
aligns.top <- lapply( aligns, function(x){ lapply(x, al.summarise, sp.class=sp.class )})
aligns.ctl.top <- lapply( aligns.ctl, function(x){ lapply(x, al.summarise, sp.class=sp.class )})
aligns.short.top <- lapply( aligns.short, al.summarise, sp.class=sp.class )
aligns.short.2.top <- lapply( aligns.short.2, al.summarise, sp.class=sp.class )

aligns.top$short <- aligns.short.top
aligns.top$short.2 <- aligns.short.2.top

saveRDS(aligns.top, 'aligns_top.rds')
saveRDS(aligns.ctl.top, 'aligns_ctl_top.rds')
saveRDS(aligns.short.top, 'aligns_short_top.rds')
saveRDS(aligns.short.2.top, 'aligns_short_2_top.rds')

## let us redo the blast, this time embedding query identifiers that
## allow us to directly get the relevant data out.

## the teleost sequences
long.tel.f <- export.al.seq( aligns.top$long, 'teleostei', 'long', 'extracted_seqs', orth )
med.tel.f <- export.al.seq( aligns.top$med, 'teleostei', 'med', 'extracted_seqs', orth )
ctl.tel.f <- export.al.seq( aligns.top$ctl, 'teleostei', 'ctl', 'extracted_seqs', orth )
short.tel.f <- export.al.seq( aligns.short.top, 'teleostei', 'short', 'extracted_seqs', orth )
short.tel.2.f <- export.al.seq( aligns.short.2.top, 'teleostei', 'short.2', 'extracted_seqs', orth )

## the mammalian sequnces
long.mam.f <- export.al.seq( aligns.top$long, 'mammalia', 'long', 'extracted_seqs', orth )
med.mam.f <- export.al.seq( aligns.top$med, 'mammalia', 'med', 'extracted_seqs', orth )
ctl.mam.f <- export.al.seq( aligns.top$ctl, 'mammalia', 'ctl', 'extracted_seqs', orth )
short.mam.f <- export.al.seq( aligns.short.top, 'mammalia', 'short', 'extracted_seqs', orth )
short.mam.2.f <- export.al.seq( aligns.short.2.top, 'mammalia', 'short.2', 'extracted_seqs', orth )

## and the control sequences
long.ctl.tel.f <- export.al.seq( aligns.ctl.top$long, 'teleostei', 'long.ctl', 'extracted_seqs', orth )
med.ctl.tel.f <- export.al.seq( aligns.ctl.top$med, 'teleostei', 'med.ctl', 'extracted_seqs', orth )
ctl.ctl.tel.f <- export.al.seq( aligns.ctl.top$ctl, 'teleostei', 'ctl.ctl', 'extracted_seqs', orth )

long.ctl.mam.f <- export.al.seq( aligns.ctl.top$long, 'mammalia', 'long.ctl', 'extracted_seqs', orth )
med.ctl.mam.f <- export.al.seq( aligns.ctl.top$med, 'mammalia', 'med.ctl', 'extracted_seqs', orth )
ctl.ctl.mam.f <- export.al.seq( aligns.ctl.top$ctl, 'mammalia', 'ctl.ctl', 'extracted_seqs', orth )

## Having exported and blasted these we now want to read in the blast data
## and calculate self-coverage..
## Since we may want to make figures for the coverage in the different species
## we also read in the other blast data although this is mostly wasteful.

## let us obtain the names of the databases from:
tmp <- readLines("../family_members/vertebrate_family_members_1_130.txt", n=1)
dbs <- strsplit(tmp, "\t")[[1]]
rm(tmp)
names(dbs) <- sub( "([^_]+)_([^_]+)_.+", "\\1.\\2", dbs )

tmp <- list.files( "blast", pattern='bl$', full.names=TRUE )
bl.files <- parse.bl.names( tmp )
rm(tmp)

## The above command would give a different result now as we
## have additional files present;

tmp <- list.files("blast", pattern='short.+bl$', full.names=TRUE)
bl.files.2 <- parse.bl.names(tmp)
rm(tmp)

## And for the second set of short introns:
tmp <- list.files("blast", pattern='short\\.2.+bl$', full.names=TRUE)
bl.files.3 <- parse.bl.names(tmp)
rm(tmp)

## we will then have to merge the blast files in a reasonable manner..
bl.columns <- strsplit('qseqid sseqid qlen qstart qend sstart send evalue bitscore score length pident nident qcovs qcovhsp', ' ')[[1]]

bl.data <- lapply( bl.files$fn, function(x){
    tmp <- read.table(x, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(tmp) <- bl.columns
    tmp
})
## and that takes up about 10GB of memory.

bl.data.2 <- lapply( bl.files.2$fn, function(x){
    tmp <- read.table(x, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(tmp) <- bl.columns
    tmp
})

bl.data.3 <- lapply( bl.files.3$fn, function(x){
    tmp <- read.table(x, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(tmp) <- bl.columns
    tmp
})

saveRDS(bl.data.3, 'bl_data_3.rds')
## calculate blast coverage for the set of data:

bl.cov <- lapply( c('teleostei', 'mammalia'), function(cl){
    tmp <- lapply( c('long', 'med', 'ctl'), function(sm){
        tmp <- mclapply( unique(bl.files$sp), function(sp){
            i <- which( bl.files$sp == sp &
                        bl.files$cl == cl &
                        bl.files$sm == sm )
            if( length(i) != 1 )
                stop(paste(c(cl, sm, sp, i), collapse=" "))
            query.covs( bl.data[[ i ]], aligns.top[[sm]], cl, exclude.sub='ALT' )
        }, mc.cores=22 )
        names(tmp) <- unique(bl.files$sp)
        tmp
    })
    names(tmp) <- c('long', 'med', 'ctl')
    tmp
})
names(bl.cov) <- c('teleostei', 'mammalia')

bl.cov.2 <- lapply( c('teleostei', 'mammalia'), function(cl){
    sm <- 'short'
    tmp <- mclapply( unique(bl.files.2$sp), function(sp){
        i <- which( bl.files.2$sp == sp &
                    bl.files.2$cl == cl &
                    bl.files.2$sm == sm )
        if( length(i) != 1 )
                stop(paste(c(cl, sm, sp, i), collapse=" "))
        query.covs( bl.data.2[[ i ]], aligns.short.top, cl, exclude.sub='ALT' )
    }, mc.cores=22 )
    names(tmp) <- unique(bl.files$sp)
    tmp
})
names(bl.cov.2) <- c('teleostei', 'mammalia')

## that took far more time than it should have done, and did not seem to use mclapply at all:
## did not see more than one process running.
saveRDS(bl.cov.2, "bl_cov_2.rds")

bl.cov.3 <- lapply( c('teleostei', 'mammalia'), function(cl){
    sm <- 'short.2'
    tmp <- mclapply( unique(bl.files.3$sp), function(sp){
        i <- which( bl.files.3$sp == sp &
                    bl.files.3$cl == cl &
                    bl.files.3$sm == sm )
        if( length(i) != 1 )
                stop(paste(c(cl, sm, sp, i), collapse=" "))
        query.covs( bl.data.3[[ i ]], aligns.short.2.top, cl, exclude.sub='ALT' )
    }, mc.cores=22 )
    names(tmp) <- unique(bl.files.3$sp)
    tmp
})
names(bl.cov.3) <- c('teleostei', 'mammalia')

saveRDS(bl.cov.3, "bl_cov_3.rds")

bl.ctl.cov <- lapply( c('teleostei', 'mammalia'), function(cl){
    tmp.1 <- lapply( c('long', 'med', 'ctl'), function(sm){
        bl.sm <- paste(sm, ".ctl", sep="")
        tmp <- mclapply( unique(bl.files$sp), function(sp){
            i <- which( bl.files$sp == sp &
                        bl.files$cl == cl &
                        bl.files$sm == bl.sm )
            if( length(i) != 1 )
                stop(paste(c(cl, sm, sp, i), collapse=" "))
            cat("i :", i, "\n")
            query.covs( bl.data[[ i ]], aligns.ctl.top[[sm]], cl, exclude.sub='ALT' )
        }, mc.cores=22 )
        names(tmp) <- unique(bl.files$sp)
        tmp
    })
    names(tmp.1) <- c('long', 'med', 'ctl')
    tmp.1
})
names(bl.ctl.cov) <- c('teleostei', 'mammalia')

saveRDS(bl.ctl.cov, "bl_ctl_cov.rds")

## the self-blast to detect repetitive regions
q.probs <- seq(0, 1, 0.05)
bl.dr.cov.q <- lapply(bl.cov, function(x){
    lapply(x, function(y){
        t(sapply(y[['danio.rerio']], quantile, probs=q.probs))
    })
})

bl.dr.cov.2.q <- lapply( bl.cov.2, function(x){
    t(sapply(x[['danio.rerio']], quantile, probs=q.probs))
})

bl.dr.cov.3.q <- lapply( bl.cov.3, function(x){
    t(sapply(x[['danio.rerio']], quantile, probs=q.probs))
})

bl.dr.ctl.cov.q <- lapply(bl.ctl.cov, function(x){
    lapply(x, function(y){
        t(sapply(y[['danio.rerio']], quantile, probs=q.probs))
    })
})

## after this we do not use bl.cov, so we can add the bl.dr.cov.2.q to
## bl.dr.cov;
## this will be complicated at some point in the future. But for now:
bl.dr.cov.q$teleostei$short <- bl.dr.cov.2.q$teleostei
bl.dr.cov.q$mammalia$short <- bl.dr.cov.2.q$mammalia

bl.dr.cov.q$teleostei$short.2 <- bl.dr.cov.3.q$teleostei
bl.dr.cov.q$mammalia$short.2 <- bl.dr.cov.3.q$mammalia

## get the scores for a given class and sample
align.top.pars <- lapply(names(aligns.top), function(sm){
    tmp <- lapply( names(aligns.top[[sm]][[1]]), function(cl){
        get.align.pars( aligns.top, sm, cl)
    })
    names(tmp) <- names(aligns.top[[sm]][[1]])
    tmp
})
names(align.top.pars) <- names(aligns.top)

align.short.top.pars <- lapply( names(aligns.short.top[[1]]), function(cl){
    get.align.pars( aligns.short.top, NULL, cl )
})
names(align.short.top.pars) <- names(aligns.short.top[[1]])

align.top.pars$short <- align.short.top.pars

align.short.2.top.pars <- lapply( names(aligns.short.2.top[[1]]), function(cl){
    get.align.pars( aligns.short.2.top, NULL, cl )
})
names(align.short.2.top.pars) <- names(aligns.short.2.top[[1]])

align.top.pars$short <- align.short.top.pars
align.top.pars$short.2 <- align.short.2.top.pars

align.top.ctl.pars <- lapply(names(aligns.ctl.top), function(sm){
    tmp <- lapply( names(aligns.ctl.top[[sm]][[1]]), function(cl){
        get.align.pars( aligns.ctl.top, sm, cl)
    })
    names(tmp) <- names(aligns.ctl.top[[sm]][[1]])
    tmp
})
names(align.top.ctl.pars) <- names(aligns.ctl.top)

par(mfrow=c(2,5))
align.top.lm <- plot.lm( c('teleostei', 'mammalia'), c('ctl', 'short', 'short.2', 'med', 'long'),
                         align.top.pars, bl.dr.cov.q )

## the min.al.score here actually refers to the score for the intron, and is used
## to remove introns whose orthology are not that clear
align.top.lm2 <- plot.lm( c('teleostei', 'mammalia'), c('ctl', 'short', 'short.2', 'med', 'long'),
                         align.top.pars, bl.dr.cov.q, min.al.score=2 )

par(mfrow=c(2,3))
align.ctl.top.lm <- plot.lm( c('teleostei', 'mammalia'), c('ctl', 'med', 'long'),
                            align.top.ctl.pars, bl.dr.ctl.cov.q )

extract.pts <- function(dt, v, use.b=FALSE){
    unlist(lapply(dt, function(x){
        unlist(lapply(x, function(y){
            if(use.b)
                y$pts[[v]][ y$pts$b ]
            else
                y$pts[[v]]}))
    }))
}

xlim <- range( extract.pts( align.top.lm, 'x', use.b=TRUE ), na.rm=TRUE )
ylim <- range( extract.pts( align.top.lm, 'y', use.b=TRUE ), na.rm=TRUE )

par(mfrow=c(2,5))
for(cl in c('teleostei', 'mammalia')){
    for(sm in c('ctl', 'short', 'short.2', 'med', 'long')){
        x <- align.top.lm[[cl]][[sm]]
        plot.points( x$pts, x$lm, xlim=xlim, ylim=ylim, main=paste(cl, sm) )
        abline( align.ctl.top.lm$teleostei$long$lm, col='purple', lwd=2 )
        abline( align.ctl.top.lm$teleostei$ctl$lm, col='purple', lwd=2, lty=2 )
        abline( align.ctl.top.lm$mammalia$long$lm, col='green', lwd=2 )
        abline( align.ctl.top.lm$mammalia$ctl$lm, col='green', lwd=2, lty=2 )
    }
}

## lets look at the control data and models
par(mfrow=c(2,3))
for(cl in c('teleostei', 'mammalia')){
    for(sm in c('ctl', 'med', 'long')){
        x <- align.ctl.top.lm[[cl]][[sm]]
        plot.points( x$pts, x$lm, xlim=xlim, ylim=ylim, main=paste(cl, sm) )
    }
}

tel.ctl.x <- extract.pts(align.ctl.top.lm['teleostei'], 'x', use.b=TRUE )
tel.ctl.y <- extract.pts(align.ctl.top.lm['teleostei'], 'y', use.b=TRUE )

mam.ctl.x <- extract.pts(align.ctl.top.lm['mammalia'], 'x', use.b=TRUE )
mam.ctl.y <- extract.pts(align.ctl.top.lm['mammalia'], 'y', use.b=TRUE )


## Make a model from these and use that for calculating residuals and
## so on..
tel.ctl.lm <- lm( tel.ctl.y ~ tel.ctl.x )
mam.ctl.lm <- lm( mam.ctl.y ~ mam.ctl.x )
## p-values for these models are of course very good
## e-228, and e-194
## R-squared: 0.45 and 0.46

res.qp <- seq(0, 1, 0.05)
tel.ctl.lm.rq <- quantile( tel.ctl.lm$residuals, probs=res.qp )
mam.ctl.lm.rq <- quantile( mam.ctl.lm$residuals, probs=res.qp )



## calculate residuals on the fly as this is not expensive
plot.lm.q <- function(pts, mod, res.qnt, qnt.l, ...){
    pts.res <- lm.res(pts, mod, res.q=res.qp)
    col <- ifelse( pts.res$res > res.qnt[qnt.l], 'red', 'black' )
    with(pts, plot(x[b], y[b], col=col[b], ...))
    abline(mod, lwd=2, col='red')
    invisible( (pts.res$res > res.qnt[qnt.l])[pts$b] )
}

cex.lab=1.25
cex.axis=1.15

## make this a supplementary figure: then it can be bigger:
cairo_pdf("alignment_scores_filtered.pdf", width=0.9 * a4.h * pdf.m, height=0.6 * a4.w * pdf.m)
par(mfrow=c(2,5))
par(mar=c(5.1, 4.6, 4.1, 2.1))
cl <- 'teleostei'
m.labs <- c(ctl='A', 'short'='B', 'short.2'='C', med='D', long='E')
for(sm in c('ctl', 'short', 'short.2', 'med', 'long')){
    b <- plot.lm.q( align.top.lm[[cl]][[sm]]$pts,
                   tel.ctl.lm, tel.ctl.lm.rq, '95%', xlim=xlim, ylim=ylim,
                   main=paste(cl, sm), xlab='log2 search space', ylab='log2 score',
                   cex.lab=cex.lab, cex.axis=cex.axis)
    text(xlim[1], ylim[2],
         sprintf("%d / %d (%.2f%%)", sum(b, na.rm=TRUE), length(b) - sum(is.na(b)),
                 100 * sum(b, na.rm=TRUE) / (length(b) - sum(is.na(b))) ),
         adj=c(0,1), cex=cex.axis)
    with(par(), mtext(m.labs[sm], at=usr[1], cex=mt.cex, line=1))
}
## 
cl <- 'mammalia'
m.labs <- c(ctl='F', 'short'='G', 'short.2'='H', med='I', long='J')
for(sm in c('ctl', 'short', 'short.2', 'med', 'long')){
    b <- plot.lm.q( align.top.lm[[cl]][[sm]]$pts,
                   mam.ctl.lm, mam.ctl.lm.rq, '95%', xlim=xlim, ylim=ylim,
                   main=paste(cl, sm), xlab='log2 search space', ylab='log2 score',
                   cex.lab=cex.lab, cex.axis=cex.axis)
    text(xlim[1], ylim[2],
         sprintf("%d / %d (%.2f%%)", sum(b, na.rm=TRUE), length(b) - sum(is.na(b)),
                 100 * sum(b, na.rm=TRUE) / (length(b) - sum(is.na(b))) ),
         adj=c(0,1), cex=cex.axis)
    with(par(), mtext(m.labs[sm], at=usr[1], cex=mt.cex, line=1))
}
dev.off()

## The following use align.top.lm2 which are filtered for the intron alignment
## score. The med and long do show an increase in conserved sequences, but the
## ctl sequences are fairly static. The improvements though are fairly marginal
## and it is probably not wortwhile to add this step here.
par(mfrow=c(2,4))
par(mar=c(5.1, 4.6, 4.1, 2.1))
cl <- 'teleostei'
m.labs <- c(ctl='A', 'short.2'='B', med='C', long='D')
for(sm in c('ctl', 'short.2', 'med', 'long')){
    b <- plot.lm.q( align.top.lm2[[cl]][[sm]]$pts,
                   tel.ctl.lm, tel.ctl.lm.rq, '95%', xlim=xlim, ylim=ylim,
                   main=paste(cl, sm), xlab='log2 search space', ylab='log2 score',
                   cex.lab=cex.lab, cex.axis=cex.axis)
    text(xlim[1], ylim[2],
         sprintf("%d / %d (%.2f%%)", sum(b, na.rm=TRUE), length(b) - sum(is.na(b)),
                 100 * sum(b, na.rm=TRUE) / (length(b) - sum(is.na(b))) ),
         adj=c(0,1), cex=cex.axis)
    with(par(), mtext(m.labs[sm], at=usr[1], cex=mt.cex, line=1))
}
##
cl <- 'mammalia'
m.labs <- c(ctl='E', short.2='F', med='G', long='H')
for(sm in c('ctl', 'short.2', 'med', 'long')){
    b <- plot.lm.q( align.top.lm2[[cl]][[sm]]$pts,
                   mam.ctl.lm, mam.ctl.lm.rq, '95%', xlim=xlim, ylim=ylim,
                   main=paste(cl, sm), xlab='log2 search space', ylab='log2 score',
                   cex.lab=cex.lab, cex.axis=cex.axis)
    text(xlim[1], ylim[2],
         sprintf("%d / %d (%.2f%%)", sum(b, na.rm=TRUE), length(b) - sum(is.na(b)),
                 100 * sum(b, na.rm=TRUE) / (length(b) - sum(is.na(b))) ),
         adj=c(0,1), cex=cex.axis)
    with(par(), mtext(m.labs[sm], at=usr[1], cex=mt.cex, line=1))
}


## do the same plots but for the control alignments
## for a supplementary figure
cairo_pdf("ctl_alignment_scores_filtered.pdf", width=0.9 * a4.w * pdf.m, height=0.6 * a4.w * pdf.m)
par(mfrow=c(2,3))
par(mar=c(5.1, 4.6, 4.1, 2.1))
cl <- 'teleostei'
m.labs <- c(ctl='A', med='B', long='D')
for(sm in c('ctl', 'med', 'long')){
    b <- plot.lm.q( align.ctl.top.lm[[cl]][[sm]]$pts,
                   tel.ctl.lm, tel.ctl.lm.rq, '95%', xlim=xlim, ylim=ylim,
                   main=paste(cl, sm), xlab='log2 search space', ylab='log2 score',
                   cex.lab=cex.lab, cex.axis=cex.axis)
    text(xlim[1], ylim[2],
         sprintf("%d / %d (%.2f%%)", sum(b, na.rm=TRUE), length(b) - sum(is.na(b)),
                 100 * sum(b, na.rm=TRUE) / (length(b) - sum(is.na(b))) ),
         adj=c(0,1), cex=cex.axis)
    with(par(), mtext(m.labs[sm], at=usr[1], cex=mt.cex, line=1))
}
##
cl <- 'mammalia'
m.labs <- c(ctl='E', med='F', long='G')
for(sm in c('ctl', 'med', 'long')){
    b <- plot.lm.q( align.ctl.top.lm[[cl]][[sm]]$pts,
                   mam.ctl.lm, mam.ctl.lm.rq, '95%', xlim=xlim, ylim=ylim,
                   main=paste(cl, sm), xlab='log2 search space', ylab='log2 score',
                   cex.lab=cex.lab, cex.axis=cex.axis)
    text(xlim[1], ylim[2],
         sprintf("%d / %d (%.2f%%)", sum(b, na.rm=TRUE), length(b) - sum(is.na(b)),
                 100 * sum(b, na.rm=TRUE) / (length(b) - sum(is.na(b))) ),
         adj=c(0,1), cex=cex.axis)
    with(par(), mtext(m.labs[sm], at=usr[1], cex=mt.cex, line=1))
}
dev.off()

### in order to look at the actual alignments that were responsible for the
### identification of the orthologs that were used for the alignments we can
### either import the alignment structures here, or export the aligns.top
### structure.

## in fact we can simply use:
## ex_align_id.rds
## which takes up about 4GB of space on disk.

ex.align.id <- readRDS("../R_172_genomes/ex_align_id.rds")

## we plotted values from align.top.lm
## for samples long, med, ctl
## and classes teleost and mammalia

## for a given class and sample we have the following vector of values:
## x, y, b, na  (as part of the pts object)
## which does not give us the gene, nor the intron.
##
## To get these we have to go back to the function that created
## align.top.lm
##
## plot.lm( c('teleostei', 'mammalia'), c('ctl', 'med', 'long'),
##                         align.top.pars, bl.dr.cov.q )
##
## align.top.pars holds the alignment scores in the 'i.score' column
##
## the alignments are filtered in the plot.al.scores function.
## align.top.pars [[ sample ]] [[ long ]]
##
## align.top.pars includes everything that we need, 
## the species information is given as the orth.sc column expressed
## as a double.

## we can get the plotted or not from align.top.lm

nuc.col <- 'grey'
gap.col <- 'white'
i.col <- 'red'

al.cols <- c(gap.col, rep(nuc.col, 5), i.col)
names(al.cols) <- c('-', 'A', 'C', 'T', 'G', 'N', 'I')

al.cols.2 <- c('white', 'grey', hsvScale(1:4, sat=1, val=0.75, max.v=5))
names(al.cols.2) <- c('-', 'N', 'A', 'C', 'G', 'T')

plot.gene.align <- function(i, sm='long', cl='teleostei'){
    x <- align.top.pars[[sm]][[cl]][i,]
    orth.j <- x['i']
    sp <- colnames(orth.sc$i)[x['sp2']]
    fam <- orth$fam[ orth.j, 1 ]
    y <- ex.align.id[[fam]][[1]]
    y.i <- which( y$sp.2 == sub(".", " ", sp, fixed=TRUE) )
    if(length(y.i) > 1){
        scores <- sapply(y$al, function(z){ z$stats['score'] })
        y.i <- y.i[ which.max( scores[y.i] ) ]
    }
    ## to see the alignment:
    ## align.print( y$al[[y.i]]$seq, w=60 )
    plot.new()
    xlim <- c(0, nchar(y$al[[y.i]]$seq[1]))
    ylim <- c(-4, 5)
    plot.window(xlim=xlim, ylim=ylim, xaxs='i')
    draw.aligns( y$al[[y.i]], 4, 2, 1, al.cols, id.lwd=0)
    with( y$al[[y.i]], segments( a.pos, 4, a.pos, 5, col='red', lend=2))
    with( y$al[[y.i]], segments( b.pos, 2, b.pos, 3, col='red', lend=2))
    intron.i <- orth.sc$i[ orth.j, sp ]
    intron.pos <- y$al[[y.i]]$b.pos[ intron.i ]
    segments(intron.pos, 1.9, xlim, 0.6, lty=3)
    axis(3)
    a.x <- cumsum(c(0, x['a_beg'], x['length'], x['l1'] - x['a_end']))
    b.x <- cumsum(c(0, x['b_beg'], x['length'], x['l2'] - x['b_end']))
    ## ever so ugly.. 
    if(x['a_beg'] > x['b_beg'])
        b.x <- b.x + (x['a_beg'] - x['b_beg'])
    else
        a.x <- a.x + (x['b_beg'] - x['a_beg'])
    plot.window(ylim=ylim, xlim=range(a.x, b.x), xaxs='i')
    rect( a.x[-4], 0.26, a.x[-1], 0.5, col=c('white', 'grey', 'white') )
    rect( b.x[-4], 0, b.x[-1], 0.24, col=c('white', 'grey', 'white') )
    segments( c(a.x[2], a.x[3]), 0, c(1, max(c(a.x[4], b.x[4]))), -0.8, lty=3 )
    plot.window(ylim=ylim, xlim=c(0, x['length']), xaxs='i')
    draw.aligns( aligns.top[[sm]][[i]][[cl]], -2, -4, 1, al.cols.2, id.lwd=0.2 )
    axis(1)
    tr.2 <- orth$tr[ orth.j, sp ]
    tr.1 <- orth$tr[ orth.j, sub(" ", ".", y$sp.1) ]
    mtext(sprintf("%s (%s), %s (%s)", to.sp(y$sp.1), tr.1, to.sp(sp), tr.2), line=2.5)
}

null.models <- list('teleostei'=tel.ctl.lm, 'mammalia'=mam.ctl.lm )

## do the residuals fit to a extreme values distribution (Gumbel distribution)
require(evd)

tmp.h <- with(null.models$teleostei, hist( residuals, freq=F, breaks=20 ))
tmp.g <- dgumbel( tmp.h$mids, -0.15, 0.3 )
points(tmp.h$mids, tmp.g)
## do they fit a gumbel distribution

null.models.evd <- lapply( null.models, function(x){ fgev( x$residuals ) })
plot(tmp.h, freq=FALSE)
points( tmp.h$mids, dgumbel( tmp.h$mids, loc=null.models.evd$teleostei$estimate['loc'], scale=null.models.evd$teleostei$estimate['scale'] ))

for(sm in c('short', 'short.2', 'long', 'med', 'ctl')){
    for(cl in c('teleostei', 'mammalia')){
        cairo_pdf(paste('exon_intron_aligns_',  sm, '_', cl, '.pdf', sep=""),
                  width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE)
##        par(mfrow=c(1,1))
        par(mfrow=c(6,1))
        par(oma=c(3,2,5,3))
        residue <- lm.res( align.top.lm[[cl]][[sm]]$pts, null.models[[cl]], use.b=FALSE )
        o <- order(residue$res, decreasing=TRUE )
        for(i in o){
##        for(i in 1:length( align.top.lm[[cl]][[sm]]$pts$b )){
            with(align.top.lm[[cl]][[sm]]$pts,
                 if(!is.na(b[i]) && b[i])
                     plot.gene.align(i, sm=sm, cl=cl) 
                 )
##            input <- readline(paste(i, ":"))
        }
        dev.off()
    }
}


## To see the correlation between different paramaters from tel.var
## we can use
## align.top.pars; which contains the columns i and j giving us
## the relevant parameters.

plot.residual.corr <- function(par, test.lm=align.top.lm, test.pars=align.top.pars,
                               ctl.lm=null.models, sms=c('ctl', 'short', 'short.2', 'med', 'long'),
                               par.table=tel.var){
    par(mfrow=c( length(test.lm), length(sms)))
    coords <- vector(mode='list', length=length(test.lm))
    names(coords) <- names(test.lm)
    for(cl in names(coords)){
        coords[[cl]] <- lapply(sms, function(sm){
            residue <- lm.res( test.lm[[cl]][[sm]]$pts, ctl.lm[[cl]], use.b=FALSE )
            b <- test.lm[[cl]][[sm]]$pts$b
            i <- test.pars[[sm]][[cl]][ ,'i']
            par.v <- par.table[i, par]
            plot( par.v[b], residue$res[b], main=paste(sm, cl), xlab=par, ylab='residual' )
            list('p'=par.v[b], 'r'=residue$res[b])
        })
        names(coords[[cl]]) <- sms
    }
    invisible(coords)
}

extract <- function(l, n){ unlist( lapply(l, function(x){ x[[n]] })) }
## med.b has some values of min that are longer than 10; this is because
## I have used int.s.b here to remove species that have too many very small
## introns (presumably due to bad annotation)

dyn.load( "~/R/blurR/src/g_blur.so" )
source(  "~/R/blurR/functions.R" )

plot.coords <- function(coords, qr=0.99, ws=100L, sd=1, set.mfrow=TRUE,
                        breaks=seq(6, 13, 0.5), ctl.res=tel.ctl.lm$residuals,
                        cl='teleostei', xlab='', lab.cex=1.25, main='', main.cex=lab.cex,
                        axis.cex=1.25, ylab.1='residual', ylim=NULL, ylim.2=NULL){
    if(set.mfrow)
        par(mfrow=c(1,1))
##    attach(coords)
    x <- extract(coords[[cl]], n='p')
    y <- extract(coords[[cl]], n='r')
    res <- sort(ctl.res)
    b.99 <- rep(FALSE, length(x))
    b.99[ y > res[ length(res) * qr ] ] <- TRUE
    col <- ifelse(b.99, 'red', rgb(0.5, 0.5, 0.5))
    ## stupid function..
    counts <- sapply(2:length(breaks), function(i){
        b <- x >= breaks[i-1] & x <= breaks[i]
        sum(b.99[b], na.rm=TRUE) / sum(b, na.rm=TRUE) })
    plot.new()
    if(is.null(ylim)){
        ylim <- range(y, na.rm=TRUE)
    }
    plot.window(xlim=range(x, na.rm=TRUE), ylim=ylim)
    mtext(xlab, side=1, line=2.5, cex=lab.cex)
    mtext(main, side=3, line=1.5, cex=lab.cex)
    mtext(ylab.1, side=2, cex=lab.cex, line=2.5)
##    mtext(paste('proportion above ', qr * 100, '%th ctl residual'), side=4, cex=lab.cex, line=2.5)
    points( x, y, col=col, cex=0.5 )
    axis(1, cex.axis=axis.cex)
    axis(2, cex.axis=axis.cex)
    if(is.null(ylim.2)) ylim.2 = c(0,max(counts, na.rm=TRUE))
    plot.window(xlim=range(x, na.rm=TRUE), ylim=ylim.2)
    rect(breaks[-length(breaks)], 0, breaks[-1], counts, col=rgb(0.7,0.7,0.3, 0.3))
    ##
    o <- order(x)
    sq.m <- lin.blur( as.numeric(b.99[o]), x[o], sd, ws )
##    plot.window(xlim=range(x, na.rm=TRUE), ylim=c(0,max(counts, na.rm=TRUE)))
##    lines(sq.m[,'pos'], sq.m[,'bl'], col='blue', lwd=0.5)
    axis(4, cex.axis=axis.cex)
    invisible(cbind(x=x, y=y, b=as.numeric(b.99)))
}


min.coords <- plot.residual.corr( 'min' )
plot.coords(min.coords)

min.ctl.coords <- plot.residual.corr( 'min', test.lm=align.ctl.top.lm, test.pars=align.top.ctl.pars,
                                     ctl.lm=null.models, sms=c('ctl', 'med', 'long') )

plot.coords(min.ctl.coords, sd=0.5, ws=200L, qr=0.95)

par(mfrow=c(1,2))
plot.coords(min.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE)
plot.coords(min.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95)

q10.coords <- plot.residual.corr( '10%' )
q10.ctl.coords <- plot.residual.corr( '10%', test.lm=align.ctl.top.lm, test.pars=align.top.ctl.pars,
                                     ctl.lm=null.models, sms=c('ctl', 'med', 'long') )

## obtain the min and max residual values.. 
par(mfrow=c(1,2))
tmp.1 <- plot.coords(q10.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE, breaks=seq(6, 14, 1))
tmp.2 <- plot.coords(q10.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=seq(6, 14, 1))

tmp.3 <- plot.coords(q10.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE, cl='mammalia', breaks=seq(6, 13, 1), ctl.res=mam.ctl.lm$residuals)
tmp.4 <- plot.coords(q10.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=seq(6, 13, 1), ctl.res=mam.ctl.lm$residuals)

cairo_pdf("excess_residuals_by_length.pdf", width=0.5*a4.w * pdf.m, height=1 * a4.w * pdf.m )
par(mfrow=c(2,1))
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot.coords(q10.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=seq(6, 14, 1), xlab='10th percentile log2 teleost intron length',
            main='Teleostei', lab.cex=1.5, axis.cex=1.3)
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))
plot.coords(q10.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=seq(6, 14, 1), ctl.res=mam.ctl.lm$residuals,
            xlab='10th percentile log2 teleost intron length', main='Mammalia', lab.cex=1.5, axis.cex=1.3)
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
dev.off()


cairo_pdf("excess_residuals_by_length_2.pdf", width=0.9*a4.w * pdf.m, height=0.7 * a4.w * pdf.m )
par(mfrow=c(2,2))
par(mar=c(2.6, 2.1, 2.1, 3.1))
par(oma=c(3, 3, 2, 3))
plot.coords(q10.ctl.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=seq(6, 14, 1), xlab='',
            main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.1[,'y'], tmp.2[,'y']), na.rm=TRUE), ylim.2=c(0,0.7))
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))
mtext('non-orthologous', cex=1.5, line=1)
plot.coords(q10.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=seq(6, 14, 1), xlab='',
            main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.1[,'y'], tmp.2[,'y']), na.rm=TRUE), ylim.2=c(0,0.7))
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
with(par(), mtext('orthologous', cex=1.5, line=1))
mtext('orthologous', cex=1.5, line=1)
mtext('Teleostei', side=4, cex=1.5, line=3)
plot.coords(q10.ctl.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=seq(6, 14, 1), ctl.res=mam.ctl.lm$residuals,
            xlab='', main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.3[,'y'], tmp.4[,'y']), na.rm=TRUE), ylim.2=c(0,0.25))
with(par(), mtext('C', at=usr[1], cex=mt.cex, line=1))
plot.coords(q10.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=seq(6, 14, 1), ctl.res=mam.ctl.lm$residuals,
            xlab='', main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.3[,'y'], tmp.4[,'y']), na.rm=TRUE), ylim.2=c(0,0.25))
with(par(), mtext('D', at=usr[1], cex=mt.cex, line=1))
mtext('10th percentile log2 teleost intron length', side=1, cex=1.5, outer=TRUE, line=1)
mtext('Mammalia', side=4, cex=1.5, line=3)
mtext('residual', side=2, cex=1.5, outer=TRUE, line=1)
dev.off()

## let us check if the correlation against danio.rerio length is better or worse:
dr.l.coords <- plot.residual.corr( 'danio.rerio', par.table=log2(orth$l) )
dr.l.ctl.coords <- plot.residual.corr( 'danio.rerio', test.lm=align.ctl.top.lm, test.pars=align.top.ctl.pars,
                                     ctl.lm=null.models, sms=c('ctl', 'med', 'long'), par.table=log2(orth$l) )


par(mfrow=c(2,2))
tmp.1 <- plot.coords(dr.l.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE, breaks=seq(6, 14, 1))
tmp.2 <- plot.coords(dr.l.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=seq(6, 14, 1))

tmp.3 <- plot.coords(dr.l.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE, cl='mammalia', breaks=seq(6, 13, 1), ctl.res=mam.ctl.lm$residuals)
tmp.4 <- plot.coords(dr.l.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=seq(6, 13, 1), ctl.res=mam.ctl.lm$residuals)

## summarise
cairo_pdf("excess_residuals_by_dr_length_2.pdf", width=0.9*a4.w * pdf.m, height=0.7 * a4.w * pdf.m )
par(mfrow=c(2,2))
par(mar=c(2.6, 2.1, 2.1, 3.1))
par(oma=c(3, 3, 2, 3))
breaks <- c(6.5, seq(7, 16, 1))
plot.coords(dr.l.ctl.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=breaks, xlab='',
                     main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.1[,'y'], tmp.2[,'y']), na.rm=TRUE), ylim.2=c(0,0.55))
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))
mtext('non-orthologous', cex=1.5, line=1)
plot.coords(dr.l.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=breaks, xlab='',
                     main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.1[,'y'], tmp.2[,'y']), na.rm=TRUE), ylim.2=c(0,0.55))
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
with(par(), mtext('orthologous', cex=1.5, line=1))
mtext('orthologous', cex=1.5, line=1)
mtext('Teleostei', side=4, cex=1.5, line=3)
plot.coords(dr.l.ctl.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=breaks, ctl.res=mam.ctl.lm$residuals,
                     xlab='', main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.3[,'y'], tmp.4[,'y']), na.rm=TRUE), ylim.2=c(0,0.2))
with(par(), mtext('C', at=usr[1], cex=mt.cex, line=1))
plot.coords(dr.l.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=breaks, ctl.res=mam.ctl.lm$residuals,
                     xlab='', main='', lab.cex=1.5, axis.cex=1.3, ylab.1='', ylim=range(c(tmp.3[,'y'], tmp.4[,'y']), na.rm=TRUE), ylim.2=c(0,0.2))
with(par(), mtext('D', at=usr[1], cex=mt.cex, line=1))
mtext('log2 Danio rerio intron length', side=1, cex=1.5, outer=TRUE, line=1)
mtext('Mammalia', side=4, cex=1.5, line=3)
mtext('residual', side=2, cex=1.5, outer=TRUE, line=1)
dev.off()

### Note that the small number of points for the non-orthologous alignments is because for these the top-scoring
### alignments were usually from repetitive sequences. This is particularly obvious when comparing
### alignment_scores_filtered.pdf and ctl_alignment_scores_filtered.pdf where teleost long has
## 806 and 279 points respectively. This difference is smaller for the medium where the numbers
## are 540 / 1137. For the control points there is no appreciable difference, as might be expected.

q20.coords <- plot.residual.corr( '20%' )
q20.ctl.coords <- plot.residual.corr( '20%', test.lm=align.ctl.top.lm, test.pars=align.top.ctl.pars,
                                     ctl.lm=null.models, sms=c('ctl', 'med', 'long') )

par(mfrow=c(1,2))
plot.coords(q20.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE, breaks=seq(6, 14, 1))
plot.coords(q20.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=seq(6, 14, 1))

plot.coords(q20.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE, cl='mammalia', breaks=seq(6, 13, 1))
plot.coords(q20.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, cl='mammalia', breaks=seq(6, 13, 1))


sd.coords <- plot.residual.corr( 'sd' )
sd.ctl.coords <- plot.residual.corr('sd', test.lm=align.ctl.top.lm, test.pars=align.top.ctl.pars,
                                     ctl.lm=null.models, sms=c('ctl', 'med', 'long') )

par(mfrow=c(1,2))
plot.coords(sd.ctl.coords, sd=0.5, ws=300L, qr=0.95, set.mfrow=FALSE, breaks=seq(0, 2, 0.1))
plot.coords(sd.coords, set.mfrow=FALSE, sd=0.5, ws=300L, qr=0.95, breaks=seq(0, 2, 0.1))

par(mfrow=c(1,2))
q10.pos <- plot.coords(q10.coords, sd=0.5, ws=300L, qr=0.95, breaks=seq(4,13,0.5), set.mfrow=FALSE)
sd.pos <- plot.coords(sd.coords, sd=0.5, ws=300L, qr=0.95, breaks=seq(0, 2, 0.1), set.mfrow=FALSE)

freq.2d <- function(pos.x, pos.y, numBins1=10, numBins2=10,
                    r1=range(pos.x[,'x'], na.rm=TRUE), r2=range(pos.y[,'x'], na.rm=TRUE)){
    all.c <- discretize2d( pos.x[,'x'], pos.y[,'x'], numBins1, numBins2, r1, r2 )
    b <- as.logical(pos.x[,'b'])
    b.c <- discretize2d( pos.x[b,'x'], pos.y[b,'x'], numBins1, numBins2, r1, r2 )
    r.c <- b.c / all.c
    list('all'=all.c, 'pos'=b.c, 'r'=r.c)
}

q10.sd.2d <- freq.2d( q10.pos, sd.pos, numBins1=14, numBins2=14, r1=c(6,13), r2=c(0,2))

par(mfrow=c(2,2))
plot( q10.pos[,'x'], sd.pos[,'x'], col=ifelse(q10.pos[,'b'], 'red', 'black'), cex=0.5 )
image(q10.sd.2d$all)
image(q10.sd.2d$pos)
image(q10.sd.2d$r)

## make a c blurring function so that we can get a somewhat nicer distribution of
## points.
saveRDS( list('q10'=q10.pos, 'sd'=sd.pos, 'h.2d'=q10.sd.2d), 'disc_2d_example_rds')

## but in general this seems to say that the length is more important than the
## the variance, but long introns have less variance. (That might be an effect
## of using minimum selection; rather than anything else.
## 
