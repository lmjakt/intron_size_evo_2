## More control alignments:

## In order to answer the question if there is a sudden increase
## in conservation as min.length goes above 256 bp

## get suitable data structures from various places:

## run alignments of long introns in parallel

library(parallel)

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
source("../R_intron_alignments_2/functions.R")
source("../R_intron_alignments_6/functions.R") ## for align.introns.. 

max.l <- 1e8
## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

## Some species classification:
sp.class <- readRDS( "../R_trees_distances/sp_class_3.rds")
class.col <-  readRDS( "../R_trees_distances/class_col_3.rds")
rownames(sp.class) <- sub("_", ".", rownames(sp.class))

sp.col <- rep(class.col['others'], nrow(sp.class))
for( cl in c('teleostei', 'mammalia', 'eutheria', 'sauria', 'aves'))
    sp.col[ sp.class[,cl] ] <- class.col[cl]
names(sp.col) <- rownames(sp.class)

## orthology...
orth <- readRDS("../R_intron_alignments_summary_2/orth.rds")
orth.sc <- readRDS("../R_intron_alignments_summary_2/orth_sc_dr.rds")

## telost variance
## this one has excluded species with too many too small introns
tel.var <- readRDS("../R_intron_alignments_summary_2/tel_var.rds")

## we will also need to know where the sequences are and other such things
## the files containing intron sequences
intron.f <-  read.table("../R_172_genomes/intron_files.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
tmp <- rownames(intron.f)
intron.f <- intron.f[,1]
names(intron.f) <- tmp
rm(tmp)

## Since the controls are matched by length it would be nice to redo
## the controls. However, it is not absolutely necessary.
## But if we do, then we will need the following:
ga.int.l <- orth$l[,'danio.rerio']
ga.int.l[ is.na(ga.int.l) ] <- 0
ga.int.l.o <- order( ga.int.l )
ga.int.l.r <- rank( ga.int.l, ties.method='random' )
ga.int.l.or <- (1:length(ga.int.l))[ga.int.l.o]

## read in the previously aligned set of introns:
long.i <- readRDS("../R_intron_alignments_6/long_i.rds")
med.i <- readRDS("../R_intron_alignments_6/med_i.rds")
ctl.i <- readRDS("../R_intron_alignments_6/ctl_i.rds")

## And then devise an additional set of controls; We want tel.var[,'min']
## to be between 6.5 and 8 as this is where we are missing data

ctl.2.b <- ctl.b <- tel.var[,'min'] < 8 & tel.var[,'min'] > 6.5 & tel.var[,'n'] > 10
## do not force d.rerio to be long as this may increase the variance.
sum(ctl.2.b) ## 10242

ctl.2.i <- which(ctl.2.b)
ctl.2.i <- setdiff( ctl.2.i, c(ctl.i, med.i, long.i))
length(ctl.2.i)  ## 10182

## ctl.2.i <- sample(ctl.2.i, 2000)
saveRDS(ctl.2.i, 'ctl_2_i.rds')

sub.matrix <- make.sub.matrix()
gap <- c(-10, -2)

## check..
long.i <- long.i[ order(tel.var[long.i, 'sd']) ]

tmp <- align.introns( long.i[1], al.sp='danio.rerio', gap=gap, max.l=max.l )
alignment.figure( tmp, sp.col=sp.col, class.col=class.col )

ctl.2.aligns <- mclapply( ctl.2.i, align.introns, al.sp='danio.rerio',
                          gap=gap, max.l=max.l, mc.cores=24 )
saveRDS(ctl.2.aligns, 'ctl_2_aligns.rds')

## also create a set where we select those with low variance
## so that we can more easily compare..

ctl.3.i <- which(ctl.2.b)
ctl.3.i <- setdiff( ctl.3.i, c(ctl.2.i, ctl.i, med.i, long.i))
ctl.3.i <- ctl.3.i[ order(tel.var[ ctl.3.i, 'sd']) ]

tmp <- align.introns( ctl.3.i[1], al.sp='danio.rerio', gap=gap, max.l=max.l )
alignment.figure( tmp, sp.col=sp.col, class.col=class.col )

ctl.3.aligns <- mclapply( ctl.3.i[ 1:2000 ], align.introns, al.sp='danio.rerio',
                          gap=gap, max.l=max.l, mc.cores=24 )
saveRDS(ctl.3.aligns, 'ctl_3_aligns.rds')
