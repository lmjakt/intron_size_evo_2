## Redo all the intron alignments.
##
## Unfortunately it turns out that I messed up in the creation of the
## very important teleost.var data structure, and hence ended up
## selecting the wrong set of introns to align. This is bad, but,
## hopefully we will get better distinction when using the correct
## introns.

## The prior error arose from the fact that I had rearranged the order
## of the intron.orth tables, but I had not re-arranged the order of the
## teleost.b data structure telling me which species is a teleost.

## run alignments of long introns in parallel

library(parallel)

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
source("../R_intron_alignments_2/functions.R")
source("functions.R")

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

orth <- read.ortho()
orth.sp <- colnames(orth$l)

sp.class <- sp.class[ match( colnames(orth$l), rownames(sp.class) ), ]
all(rownames(sp.class) == colnames(orth$l))  ## TRUE

tel.var <- t(apply( log2( orth$l[ , sp.class[,'teleostei']]), 1, var.par ))

saveRDS(tel.var, "tel_var.rds")

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

### Define the set of introns to align:
long.b <- tel.var[,'min'] > 10 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10
med.b <- tel.var[,'min'] > 8 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10 & !long.b
ctl.b <- tel.var[,'50%'] < 8 & orth$l[,'danio.rerio'] > 1024 & tel.var[,'n'] > 10

sum(long.b) ## 1594
sum(med.b) ## 4631
sum(ctl.b) ## 10911

## We can subset a subset of these:
## but do not rerun these
######
###### Commented out to stop me from misadvertent sampling.
long.i <- which(long.b)
long.i <- long.i[ order(tel.var[long.i,'sd']) ]
med.i <- which(med.b)
med.i <- (med.i[ order( tel.var[med.i,'sd'] ) ])[1:2000]
## ctl.i <- sample(which(ctl.b), 2000)

saveRDS(long.i, 'long_i.rds')
saveRDS(med.i, 'med_i.rds')
saveRDS(ctl.i, 'ctl_i.rds')

sub.matrix <- make.sub.matrix()
gap <- c(-10, -2)

## check..
tmp <- align.introns( long.i[1], al.sp='danio.rerio', gap=gap, max.l=max.l )
alignment.figure( tmp, sp.col=sp.col, class.col=class.col )

## and that does now seem to be ok. Let
long.aligns <- mclapply( long.i, align.introns, al.sp='danio.rerio',
                         gap=gap, max.l=max.l, mc.cores=16 )
saveRDS( long.aligns, "long_aligns.rds" )
##
med.aligns <- mclapply( med.i, align.introns, al.sp='danio.rerio',
                         gap=gap, max.l=max.l, mc.cores=16 )
saveRDS( med.aligns, "med_aligns.rds" )
## 
ctl.aligns <- mclapply( ctl.i, align.introns, al.sp='danio.rerio',
                         gap=gap, max.l=max.l, mc.cores=16 )
saveRDS( ctl.aligns, "ctl_aligns.rds" )

## have a look at some of the alignments:
for(i in 1:length(long.aligns)){
    alignment.figure( long.aligns[[i]], sp.col=sp.col, class.col=class.col )
    inpt <- readline(paste(i, ":"))
}

for(i in 1:length(med.aligns)){
    alignment.figure( med.aligns[[i]], sp.col=sp.col, class.col=class.col )
    inpt <- readline(paste(i, ":"))
}

for(i in 1:length(ctl.aligns)){
    alignment.figure( ctl.aligns[[i]], sp.col=sp.col, class.col=class.col )
    inpt <- readline(paste(i, ":"))
}

### These do look very good. In fact they are probably better than what
### I had before. But we most definitely need to filter for simple repeats
### as before.

### although we do not absolutely need it, we can run the control alignments
## as well. 

long.ctl.aligns <- mclapply( long.i, align.introns.ctl, al.sp='danio.rerio',
                             gap=gap, max.l=max.l, mc.cores=12 )
saveRDS( long.ctl.aligns, 'long_ctl_aligns.rds')
##
med.ctl.aligns <- mclapply( med.i, align.introns.ctl, al.sp='danio.rerio',
                             gap=gap, max.l=max.l, mc.cores=12 )
saveRDS( med.ctl.aligns, 'med_ctl_aligns.rds')
##

## we got a warning message here:
## In mclapply(ctl.i, align.introns.ctl, al.sp = "danio.rerio", gap = gap,  :
##  scheduled cores 2, 12 did not deliver results, all values of the jobs will be affected

## try redoing to see what happens. 
ctl.ctl.aligns <- mclapply( ctl.i, align.introns.ctl, al.sp='danio.rerio',
                             gap=gap, max.l=max.l, mc.cores=12 )

saveRDS( ctl.ctl.aligns, 'ctl_ctl_aligns.rds' )


