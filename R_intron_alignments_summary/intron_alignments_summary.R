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

## add supplementary including local alignment scores at intron positions
orth.sc <- readRDS("../R_intron_orthology/intron_orth_lscore.rds")
## we actually only want to use the 'danio rerio' part here:
orth.sc <- orth.sc[['danio rerio']]

## modify the col names to fit with the orth tables
for(i in 1:length(orth.sc))
    colnames(orth.sc[[i]]) <- sub(" ", ".",  colnames(orth.sc[[i]]), fixed=TRUE)

all( orth.sc$i == orth$i, na.rm=TRUE )  ## TRUE..

## and the database credentials: (which have an unset password: set this manually, when needed)
db.cred <- read.credentials("../R_intron_alignments_2/extracted_seqs/blast/db.cred")

## Some species classification:
sp.class <- readRDS( "../R_trees_distances/sp_class_3.rds")
class.col <-  readRDS( "../R_trees_distances/class_col_3.rds")
rownames(sp.class) <- sub("_", ".", rownames(sp.class))

## The full set of alignments:
## (needed for visualisation, but not much else)
aligns <- list(
    'long' = c( readRDS("../R_intron_alignments/al_top_500.rds"),
               readRDS("../R_intron_alignments_5/long_2_aligns.rds")),
    'med' = readRDS("../R_intron_alignments_5/long_3_aligns.rds"),
    'ctl' = c(readRDS("../R_intron_alignments/al_ctl_500.rds"),
              readRDS("../R_intron_alignments_5/ctl_aligns.rds"))
)

aligns.ctl <- list(
    'long' = readRDS("../R_intron_alignments/al_top_500_ctl.rds"),
    'ctl' =  readRDS("../R_intron_alignments/al_ctl_500_ctl.rds")
)

## set the rownumbers of aligns so that we can access the raw information
## easily
for(sm in names(aligns)){
    names(aligns[[sm]]) <- as.character( sapply( aligns[[sm]], function(x){ x$i } ))
}

for(sm in names(aligns.ctl)){
    names(aligns.ctl[[sm]]) <- as.character( sapply( aligns.ctl[[sm]], function(x){ x$i } ))
}


## create new summary structures with a simpler NULL values for missing alignments
aligns.top <- lapply( aligns, function(x){ lapply(x, al.summarise, sp.class=sp.class )})
aligns.ctl.top <- lapply( aligns.ctl, function(x){ lapply(x, al.summarise, sp.class=sp.class )})

saveRDS(aligns.top, 'aligns_top.rds')
saveRDS(aligns.ctl.top, 'aligns_ctl_top.rds')

## let us redo the blast, this time embedding query identifiers that
## allow us to directly get the relevant data out.

## the teleost sequences
long.tel.f <- export.al.seq( aligns.top$long, 'teleostei', 'long', 'extracted_seqs', orth )
med.tel.f <- export.al.seq( aligns.top$med, 'teleostei', 'med', 'extracted_seqs', orth )
ctl.tel.f <- export.al.seq( aligns.top$ctl, 'teleostei', 'ctl', 'extracted_seqs', orth )

## the mammalian sequnces
long.mam.f <- export.al.seq( aligns.top$long, 'mammalia', 'long', 'extracted_seqs', orth )
med.mam.f <- export.al.seq( aligns.top$med, 'mammalia', 'med', 'extracted_seqs', orth )
ctl.mam.f <- export.al.seq( aligns.top$ctl, 'mammalia', 'ctl', 'extracted_seqs', orth )

## and the control sequences
long.ctl.tel.f <- export.al.seq( aligns.ctl.top$long, 'teleostei', 'long.ctl', 'extracted_seqs', orth )
ctl.ctl.tel.f <- export.al.seq( aligns.ctl.top$ctl, 'teleostei', 'ctl.ctl', 'extracted_seqs', orth )

long.ctl.mam.f <- export.al.seq( aligns.ctl.top$long, 'mammalia', 'long.ctl', 'extracted_seqs', orth )
ctl.ctl.mam.f <- export.al.seq( aligns.ctl.top$ctl, 'mammalia', 'ctl.ctl', 'extracted_seqs', orth )

## let us obtain the names of the databases from:
tmp <- readLines("../family_members/vertebrate_family_members_1_130.txt", n=1)
dbs <- strsplit(tmp, "\t")[[1]]
rm(tmp)
names(dbs) <- sub( "([^_]+)_([^_]+)_.+", "\\1.\\2", dbs )

tmp <- list.files( "blast", pattern='bl$', full.names=TRUE )
bl.files <- parse.bl.names( tmp )
rm(tmp)
## we do not include the long.ctl the sampled.ctl

## we will then have to merge the blast files in a reasonable manner..
bl.columns <- strsplit('qseqid sseqid qlen qstart qend sstart send evalue bitscore score length pident nident qcovs qcovhsp', ' ')[[1]]

## We can read in the full set of blast data. This will take something like
## 10 GB of memory. Or we can read in data when we need it.

bl.data <- lapply( bl.files$fn, function(x){
    tmp <- read.table(x, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(tmp) <- bl.columns
    tmp
})
## and that takes up about 10GB of memory.

## We now want to calculate query.covs for the aligns.top and aligns.ctl.top
## As these should correspond.
## Note that we only have blast data for teleostei and mammalia, so we
## have to limit ourselves to that.

i <- which( bl.files$sp ==  'danio.rerio' &
            bl.files$cl == 'teleostei' &
            bl.files$sm == 'long' )
system.time(
    tmp <- query.covs( bl.data[[ i ]], aligns.top$long, 'teleostei' )
 )
##    user  system elapsed 
## 807.350  50.383 857.800 
## bloody hell that is slow. Would probably be faster to write a C-extension
## to sort the problem. Nevertheless, We do need to go through it for
## all of the blast data that we have..

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
##
saveRDS(bl.cov, "bl_cov.rds")

bl.ctl.cov <- lapply( c('teleostei', 'mammalia'), function(cl){
    cat(cl, "\n")
    tmp <- lapply( c('long', 'ctl'), function(sm){
        bl.sm <- paste(sm, ".ctl", sep="")
        cat("\t", sm, " : ", bl.sm, "\n")
        tmp <- vector(mode='list', length=length(unique(bl.files$sp)) )
        names(tmp) <- unique(bl.files$sp)
        for(sp in names(tmp)){
            cat("\t", sp, "\n")
            i <- which( bl.files$sp == sp &
                        bl.files$cl == cl &
                        bl.files$sm == bl.sm )
            if( length(i) != 1 )
                stop(paste(c(cl, sm, sp, i), collapse=" "))
            tmp[[sp]] <- query.covs( bl.data[[ i ]], aligns.ctl.top[[sm]], cl, exclude.sub='ALT' )
        }
        tmp
    })
    names(tmp) <- c('long', 'ctl')
    tmp
})
names(bl.ctl.cov) <- c('teleostei', 'mammalia')
##

saveRDS(bl.ctl.cov, "bl_ctl_cov.rds")


## we should do the same for the .ctl, but we can wait with that until we have had
## a look at the others.

## the self-blast to detect repetitive regions
q.probs <- seq(0, 1, 0.05)
bl.dr.cov.q <- lapply(bl.cov, function(x){
    lapply(x, function(y){
        t(sapply(y[['danio.rerio']], quantile, probs=q.probs))
    })
})

bl.dr.ctl.cov.q <- lapply(bl.ctl.cov, function(x){
    lapply(x, function(y){
        t(sapply(y[['danio.rerio']], quantile, probs=q.probs))
    })
})

## get the scores for a given class and sample
align.top.pars <- lapply(names(aligns.top), function(sm){
    tmp <- lapply( names(aligns.top[[sm]][[1]]), function(cl){
        get.align.pars( aligns.top, sm, cl)
    })
    names(tmp) <- names(aligns.top[[sm]][[1]])
    tmp
})
names(align.top.pars) <- names(aligns.top)

align.top.ctl.pars <- lapply(names(aligns.ctl.top), function(sm){
    tmp <- lapply( names(aligns.ctl.top[[sm]][[1]]), function(cl){
        get.align.pars( aligns.ctl.top, sm, cl)
    })
    names(tmp) <- names(aligns.ctl.top[[sm]][[1]])
    tmp
})
names(align.top.ctl.pars) <- names(aligns.ctl.top)

par(mfrow=c(2,3))
align.top.lm <- plot.lm( c('teleostei', 'mammalia'), c('ctl', 'med', 'long'),
                         align.top.pars, bl.dr.cov.q )

align.top.lm2 <- plot.lm( c('teleostei', 'mammalia'), c('ctl', 'med', 'long'),
                         align.top.pars, bl.dr.cov.q, min.al.score=2 )


par(mfrow=c(2,2))
align.ctl.top.lm <- plot.lm( c('teleostei', 'mammalia'), c('ctl', 'long'),
                            align.top.ctl.pars, bl.dr.ctl.cov.q )

extract.pts <- function(dt, v, use.b=FALSE){
    unlist(lapply(dt, function(x){
        unlist(lapply(x, function(y){
            if(use.b)
                (y$pts[[v]])[ y$pts$b ]
            else
                y$pts[[v]]}))
    }))
}
    
xlim <- range( extract.pts( align.top.lm, 'x' ), na.rm=TRUE )
ylim <- range( extract.pts( align.top.lm, 'y' ), na.rm=TRUE )

par(mfrow=c(2,3))
for(cl in c('teleostei', 'mammalia')){
    for(sm in c('ctl', 'med', 'long')){
        x <- align.top.lm[[cl]][[sm]]
        plot.points( x$pts, x$lm, xlim=xlim, ylim=ylim, main=paste(cl, sm) )
        abline( align.ctl.top.lm$teleostei$long$lm, col='purple', lwd=2 )
        abline( align.ctl.top.lm$teleostei$ctl$lm, col='purple', lwd=2, lty=2 )
        abline( align.ctl.top.lm$mammalia$long$lm, col='green', lwd=2 )
        abline( align.ctl.top.lm$mammalia$ctl$lm, col='green', lwd=2, lty=2 )
    }
}

## lets look at the control data and models
par(mfrow=c(2,2))
for(cl in c('teleostei', 'mammalia')){
    for(sm in c('ctl', 'long')){
        x <- align.ctl.top.lm[[cl]][[sm]]
        plot.points( x$pts, x$lm, xlim=xlim, ylim=ylim, main=paste(cl, sm) )
    }
}

tel.ctl.x <- extract.pts(align.ctl.top.lm['teleostei'], 'x', use.b=TRUE )
tel.ctl.y <- extract.pts(align.ctl.top.lm['teleostei'], 'y', use.b=TRUE )

mam.ctl.x <- extract.pts(align.ctl.top.lm['mammalia'], 'x', use.b=TRUE )
mam.ctl.y <- extract.pts(align.ctl.top.lm['mammalia'], 'y', use.b=TRUE )

par(mfrow=c(1,2))
with( align.ctl.top.lm$teleostei$long$pts, plot(x[b], y[b], xlim=xlim, ylim=ylim, col='red'))
with( align.ctl.top.lm$teleostei$ctl$pts, points(x[b], y[b], col='blue'))
with( align.ctl.top.lm$teleostei, abline(long$lm, col='red', lwd=2))
with( align.ctl.top.lm$teleostei, abline(ctl$lm, col='blue', lwd=2))
points(tel.ctl.x, tel.ctl.y, col='grey')

with( align.ctl.top.lm$mammalia$long$pts, plot(x[b], y[b], xlim=xlim, ylim=ylim, col='red'))
with( align.ctl.top.lm$mammalia$ctl$pts, points(x[b], y[b], col='blue'))
with( align.ctl.top.lm$mammalia, abline(long$lm, col='red', lwd=2))
with( align.ctl.top.lm$mammalia, abline(ctl$lm, col='blue', lwd=2))
points(mam.ctl.x, mam.ctl.y, col='grey')

## that confirms that tel.ctl. and mam.ctl.y are the correct points. We can
## now make a model from these and use that for calculating residuals and
## so on..
tel.ctl.lm <- lm( tel.ctl.y ~ tel.ctl.x )
mam.ctl.lm <- lm( mam.ctl.y ~ mam.ctl.x )
## p-values for these models are of course very good
## e-71, and e-33

res.qp <- seq(0, 1, 0.05)
tel.ctl.lm.rq <- quantile( tel.ctl.lm$residuals, probs=res.qp )
mam.ctl.lm.rq <- quantile( mam.ctl.lm$residuals, probs=res.qp )

### Since I am having a hard time considering how to generalise this
### I will make a plot here

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

cairo_pdf("alignment_scores_filtered.pdf", width=0.9 * a4.w * pdf.m, height=0.6 * a4.w * pdf.m)

par(mfrow=c(2,3))
par(mar=c(5.1, 4.6, 4.1, 2.1))
cl <- 'teleostei'
m.labs <- c(ctl='A', med='B', long='C')
for(sm in c('ctl', 'med', 'long')){
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
m.labs <- c(ctl='D', med='E', long='F')
for(sm in c('ctl', 'med', 'long')){
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
par(mfrow=c(2,3))
par(mar=c(5.1, 4.6, 4.1, 2.1))
cl <- 'teleostei'
m.labs <- c(ctl='A', med='B', long='C')
for(sm in c('ctl', 'med', 'long')){
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
m.labs <- c(ctl='D', med='E', long='F')
for(sm in c('ctl', 'med', 'long')){
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

sm <- 'long'
cl <- 'teleostei'

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
    segments( c(a.x[2], a.x[3]), 0, c(1, a.x[4]), -0.8, lty=3 )
    plot.window(ylim=ylim, xlim=c(0, x['length']), xaxs='i')
    draw.aligns( aligns.top[[sm]][[i]][[cl]], -2, -4, 1, al.cols.2, id.lwd=0.2 )
    axis(1)
    tr.2 <- orth$tr[ orth.j, sp ]
    tr.1 <- orth$tr[ orth.j, sub(" ", ".", y$sp.1) ]
    mtext(sprintf("%s (%s), %s (%s)", to.sp(y$sp.1), tr.1, to.sp(sp), tr.2), line=2.5)
}

plot.gene.align(13, 'long', 'teleostei')  ## this looks wrong 


for(sm in c('long', 'med', 'ctl')){
    for(cl in c('teleostei', 'mammalia')){
        cairo_pdf(paste('exon_intron_aligns_',  sm, '_', cl, '.pdf', sep=""),
                  width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE)
        par(mfrow=c(6,1))
        par(oma=c(3,2,5,3))
        for(i in 1:length( align.top.lm[[cl]][[sm]]$pts$b )){
            with(align.top.lm[[cl]][[sm]]$pts,
                 if(!is.na(b[i]) && b[i])
                     plot.gene.align(i, sm=sm, cl=cl) 
                 )
            ##    input <- readline(paste(i, ":"))
        }
        dev.off()
    }
}
    
## we have some strangeness in this. In particular the 9th entry for the long intron size
## has an alignment to a Poecilia latipina, which has a very short intron here. Since
## all introns should have been long in teleosts this is problematic. Indeed the intron
## is long in most telosts, but P. latipina, suggesting that the problem may have come
## from the classification used. That was obtained from Ensembl compara 98 on my workstation
## (R_i72_genomes.R) and was used. Let us have a look at it here to see if we can see if that
## is the reason for the anomaly.

sp.lineages <- readRDS("../R_172_genomes/species_lineages.rds")
teleost.b <- read.table("../R_172_genomes/teleost_b.txt", stringsAsFactors=FALSE, header=TRUE, sep='\t')

teleost.b['poecilia latipinna',]  ## TRUE?
