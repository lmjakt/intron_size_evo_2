## Align large numbers of non-orthologous sequences (at least 10,000)
## of a constant length (1000 bp) in order to estimate the parameters lambda
## and K that describe the probabilities of a given alignment given the
## the scoring system used.

library(parallel)

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
source("../R_intron_alignments_2/functions.R")
source("../R_intron_alignments_6/functions.R") ## for align.introns.. 

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
## in order to choose non-orthologous sequeences
orth <- readRDS("../R_intron_alignments_summary_2/orth.rds")

## align danio rerio introns to introns which are generally long
tel.var <- readRDS("../R_intron_alignments_summary_2/tel_var.rds")
mam.var <- readRDS("../R_intron_alignments_summary_2/mam_var.rds")

## where can we find the sequences
intron.f <-  read.table("../R_172_genomes/intron_files.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
tmp <- rownames(intron.f)
intron.f <- intron.f[,1]
names(intron.f) <- tmp
rm(tmp)


align.length <- 1000

## the seed lengths
tel.b1 <- orth$l[,'danio.rerio'] >= align.length & tel.var[,'50%'] >= log2(align.length) & tel.var[,'n'] >= 20
mam.b1 <- orth$l[,'danio.rerio'] >= align.length & mam.var[,'50%'] >= log2(align.length) & mam.var[,'n'] >= 20
mam.b1[ is.na(mam.b1) ] <- FALSE

sum(tel.b1)  ## 10164
sum(mam.b1) ## 20827

## so we have plenty of potential target sequences..

tel.i <- which(tel.b1)
tel.j <- sample(tel.i)
sum(tel.i == tel.j) ## 0 which is nice.

mam.i <- which(mam.b1)
mam.j <- sample(mam.i)
sum(mam.i == mam.j) ## 0.. so don't change.

sub.matrix <- make.sub.matrix()
gap <- c(-10, -2)

substr <- function(str, n){
    l <- nchar(str)
    if(l < n)
        return(NULL)
    beg <- 1 + as.integer( (l-n) / 2 )
    end <- as.integer( beg + n - 1 )
    substring(str, beg, end)
}

## align a subset of sequences in i to a subset of sequences in j
align.introns <- function(i, j, i.cl='teleostei', j.cl='teleostei',
                          al.n=4, al.l=1000, min.width=15, min.score=20,
                          sm=sub.matrix, orth.fam=orth$fam, gap=c(-10, -2)){
    fam.i <- orth.fam[i, 1]
    fam.j <- orth.fam[j, 1]
    seq.i <- read.exons( intron.f[ fam.i ] )
    seq.j <- read.exons( intron.f[ fam.j ] )
    meta.i <- sapply( seq.i, function(x){ c('id'=x$id, 'tr'=x$tr, 'sp'=x$sp, 'n'=length(x$l) )})
    meta.j <- sapply( seq.j, function(x){ c('id'=x$id, 'tr'=x$tr, 'sp'=x$sp, 'n'=length(x$l) )})
    ## do this with restrictions as to the groupings
    i.b <- sp.class[ sub(' ', '.',  meta.i['sp',]), i.cl ]
    j.b <- sp.class[ sub(' ', '.',  meta.j['sp',]), j.cl ]
    ## Extract the potential intron sequences
    int.i <- unlist( lapply( seq.i[i.b], function(x){ x$e }))
    int.j <- unlist( lapply( seq.j[j.b], function(x){ x$e }))
    ## keep only the ones that are sufficiently long
    int.i <- int.i[ !is.na(int.i) & !grepl("^>", int.i) & sapply(int.i, nchar) >= al.l ]
    int.j <- int.j[ !is.na(int.j) & !grepl("^>", int.j) & sapply(int.j, nchar) >= al.l ]
    ## at the moment the sequences will be ordered by gene, intron; this may not matter so much
    ## so we can leave.
    al.n <- min( c(al.n, length(int.i), length(int.j)) )
    ## if no suitable alignments give up.. 
    if(al.n == 0)
        return(NULL)
    ## extract sub strings
    int.i <- sapply( int.i, substr, n=al.l )
    int.j <- sapply( int.j, substr, n=al.l )
    ## and let us simply reorder both of these so that we do not bias ourselves for specific
    ## species pairs.
    int.i <- int.i[ sample(1:length(int.i)) ]
    int.j <- int.j[ sample(1:length(int.j)) ]
    lapply( 1:al.n, function(k){
        local.aligns( int.i[k], int.j[k], sm$offset, sm$size, sm$sm,
                     gap=gap, min.width=min.width, min.score=min.score )
    })
}

## let us do 1000 * 20 for tel.tel and tel.mam alignments to see what these look like
## need to rerun this as I accidentally overwrote it further below
tel.tel.1000 <- mclapply( 1:5000, function(i){
    align.introns( tel.i[i], tel.j[i], al.n=20, al.l=1000,
                   i.cl='teleostei', j.cl='teleostei' )
}, mc.cores=20)

## and that was a bit faster than I expected..
table( sapply(tel.tel.1000, length ))
##  13   20 
##   2 4998 


tel.mam.1000 <- mclapply( 1:5000, function(i){
    align.introns( tel.i[i], mam.j[i], al.n=20, al.l=1000,
                   i.cl='teleostei', j.cl='mammalia' )
}, mc.cores=20)

table( sapply(tel.mam.1000, length ))
##  13   20 
##   1 4999 

tel.tel.1000.sc <- unlist( sapply(tel.tel.1000, function(x){
    sapply(x, function(y){ y$pos[1,'score'] }) }))

length(tel.tel.1000.sc)
## [1] 96452

tel.mam.1000.sc <- unlist( sapply(tel.mam.1000, function(x){
    sapply(x, function(y){ y$pos[1,'score'] }) }))
length(tel.mam.1000.sc)
## [1] 96064

tel.h <- hist( tel.tel.1000.sc, breaks=100, freq=FALSE )
mam.h <- hist( tel.mam.1000.sc, breaks=tel.h$breaks, freq=FALSE )


## these are very similar.. the mammalian one is shifted slightly to the left.
## i.e. has slighlty lower scores. Now can we estimate the relevant probabilties:
## I found the relevant equations at:
## www.cbcb.umd.edu/confcour/Statistics_Talk_UMD.pdf
## which seems to be a presenation by Stephan Altschul

## seems that I can get the paper at Researchgate
## https://www.researchgate.net/publication/251452626_27_Local_alignment_statistics

require('evd')

plot(tel.h$mids, tel.h$density, type='p', lwd=2, col='blue')
lines(mam.h$mids, mam.h$density, type='p', lwd=2, col='red')
lines(tel.h$mids, dgumbel(tel.h$mids, loc=65, scale=10 ))

tel.evd <- fgev( tel.tel.1000.sc )
mam.evd <- fgev( tel.mam.1000.sc )

## these fit the gumbel distribution nicely in that the shape parameter is
## very close to 0

plot(tel.h$mids, (tel.h$density), type='p', lwd=2, col='blue')
lines(mam.h$mids, (mam.h$density), type='p', lwd=2, col='red')
##
lines(tel.h$mids, (with(tel.evd, dgev(tel.h$mids, loc=estimate['loc'], scale=estimate['scale'], shape=estimate['shape']))), col='blue' )
lines(mam.h$mids, (with(mam.evd, dgev(mam.h$mids, loc=estimate['loc'], scale=estimate['scale'], shape=estimate['shape']))), col='red' )


## looks like I can use these to scale the alignments
## to p values directly.
## But let us make some for longer alignments as well.

## the seed lengths
align.length = 2000
tel.b2 <- tel.var[,'50%'] >= log2(align.length) & tel.var[,'n'] >= 20
mam.b2 <- mam.var[,'50%'] >= log2(align.length) & mam.var[,'n'] >= 20
mam.b2[ is.na(mam.b1) ] <- FALSE

sum(tel.b2)  ## 5921
sum(mam.b2) ## 24976

## so we have plenty of potential target sequences..

tel.i2 <- which(tel.b2)
tel.j2 <- sample(tel.i2)
sum(tel.i2 == tel.j2) ## 0 which is nice.

mam.i2 <- which(mam.b2)
mam.j2 <- sample(mam.i2)
sum(mam.i2 == mam.j2) ## 0.. so don't change.

tel.tel.2000 <- mclapply( 1:5000, function(i){
    align.introns( tel.i2[i], tel.j2[i], al.n=20, al.l=2000,
                   i.cl='teleostei', j.cl='teleostei' )
}, mc.cores=20)
##
tel.mam.2000 <- mclapply( 1:5000, function(i){
    align.introns( tel.i2[i], mam.j2[i], al.n=20, al.l=2000,
                   i.cl='teleostei', j.cl='mammalia' )
}, mc.cores=20)

table( sapply(tel.mam.1000, length ))
##  13   20 
##   1 4999 

tel.tel.2000.sc <- unlist( sapply(tel.tel.2000, function(x){
    sapply(x, function(y){ y$pos[1,'score'] }) }))

length(tel.tel.2000.sc)
## [1] 97096

tel.mam.2000.sc <- unlist( sapply(tel.mam.2000, function(x){
    sapply(x, function(y){ y$pos[1,'score'] }) }))

length(tel.mam.2000.sc)
## [1] 96916

tel.h2 <- hist( tel.tel.2000.sc, breaks=100 )
mam.h2 <- hist( tel.mam.2000.sc, breaks=100 )

plot( sort(tel.tel.1000.sc ))
plot( sort(tel.tel.2000.sc ))

tel2.evd <- fgev( tel.tel.2000.sc )
mam2.evd <- fgev( tel.mam.2000.sc )

## Altschul and Gish provide an alternative parameterisation of the Gumbel distribution
## in particular it would seem that they use an inverse scale.
## To confirm this, generate numbers using the equation used by Altschul and Gish
## to check that we get the same numbers
alt.gumbel <- function(x, lambda, mu){
    exp( -lambda * (x - mu) - exp(-lambda * (x-mu)) )
}

## from the example in their presentation
lambda <- 0.27
mu <- 39.75
scores <- 30:80
## for mn = 1097 * 1097

plot( scores, lambda * alt.gumbel( scores, lambda, mu), type='b' )
points( scores, dgumbel(scores, mu, 1/lambda), type='b', col='red')
## And these are identical. So we can simply 1/scale for the additional calculations

## takes an object returned by fgev
## and reverses the equation: mu = ln(Kmn) / lambda
## K = exp( mu * lambda ) / mn
gumbel2K <- function(d, mn){
    lambda = 1 / d$estimate['scale']
    mu = d$estimate['loc']
    k <- exp( mu * lambda ) / mn
    names(k) <- 'K'
    k
}

## test it for the above values
gumbel2K( tel.evd, 1e6 )

## that seems to work fine. Now let us obtain estimates of K for the two different
## alignments.

gumbel2K( tel.evd, 1e6 )  ## 0.0001876942
gumbel2K( tel2.evd, 4e6 )  ## 4.502316e-05
## hmm. These are very different to each other
## But note that K is smaller for the longer sequences; since
## E = Kmn * exp(-lambda * S)
## this suggests that we will get smaller E values for longer
## alignments. This shift is actually similar to that seen in the
## Altschul and Gish, but much larger in scale.

gumbel2K( mam.evd, 1e6 )  ## 0.0002063511
gumbel2K( mam2.evd, 1e6 )  ## 0.0002380221
## these though are rather more similar to each other
## suggesting that they are closer to random sequences.

## this sugests that we have to generate random sequences. This is sort of easy
## but I need to have the counts first. Before we do that I can run more intron
## alignments as above.

align.sets <- function( align.length ){
    tel.b <- orth$l[,'danio.rerio'] >= align.length & tel.var[,'50%'] >= log2(align.length) & tel.var[,'n'] >= 20
    mam.b <- orth$l[,'danio.rerio'] >= align.length & mam.var[,'50%'] >= log2(align.length) & mam.var[,'n'] >= 20
    mam.b[ is.na(mam.b1) ] <- FALSE
    ##
    tel.i <- which(tel.b1)
    tel.j <- sample(tel.i)
    while( sum(tel.i == tel.j) )
        tel.j <- sample(tel.i)
    mam.i <- which(mam.b1)
    mam.j <- sample(mam.i)
    while( sum(mam.i == mam.j) )
        mam.j <- sample(mam.i)
    ##
    ## then we do the alignments and extract the scores...
    tel.al <- mclapply( 1:1000, function(i){
        align.introns(tel.i[i], tel.j[i], al.n=20, al.l=align.length,
                   i.cl='teleostei', j.cl='teleostei' )
    }, mc.cores=20)
    mam.al <- mclapply( 1:1000, function(i){
        align.introns( tel.i[i], mam.j[i], al.n=20, al.l=align.length,
                      i.cl='teleostei', j.cl='mammalia' )
    }, mc.cores=20)
    ##
    tel.scores <- unlist( sapply( tel.al, function(x){
        sapply(x, function(y){ y$pos[1,'score'] })
    }))
    mam.scores <- unlist( sapply( mam.al, function(x){
        sapply(x, function(y){ y$pos[1,'score'] })
    }))
    tel.l <- unlist( sapply( tel.al, function(x){
        sapply(x, function(y){ y$pos[1,'length'] })
    }))
    mam.l <- unlist( sapply( mam.al, function(x){
        sapply(x, function(y){ y$pos[1,'length'] })
    }))
    ## and run fgev on the scores
    tel.gev <- fgev( tel.scores )
    mam.gev <- fgev( mam.scores )
    list('tel.al'=tel.al, 'mam.al'=mam.al, 'tel.s'=tel.scores, 'tel.l'=tel.l, 'mam.s'=mam.scores, 'mam.l'=mam.l,
         'tel.gev'=tel.gev, 'mam.gev'=mam.gev)
}

al.sets <- lapply( 2^(6:12), align.sets )

al.sets.n <- 2^(6:12)
al.sets.mn <-al.sets.n^2

## ec for corrected for edge effects (t for teleost)
al.sets.mn.t.ec <- (al.sets.n - sapply( al.sets, function(x){ mean( x$tel.l ) }))^2
al.sets.mn.m.ec <- (al.sets.n - sapply( al.sets, function(x){ mean( x$mam.l ) }))^2

al.sets.tgev <- sapply( al.sets, function(x){ x$tel.gev$estimate })
al.sets.mgev <- sapply( al.sets, function(x){ x$mam.gev$estimate })
al.sets.tgev.std <- sapply( al.sets, function(x){ x$tel.gev$std.err })
al.sets.mgev.std <- sapply( al.sets, function(x){ x$mam.gev$std.err })


plot( log(al.sets.mn), al.sets.tgev['loc',], type='b', col='blue' )
abline( lm( al.sets.tgev['loc', ] ~ log(al.sets.mn) ), col='blue', lty=2 )
points( log(al.sets.mn), al.sets.mgev['loc',], type='b', col='red' )

plot( log(al.sets.mn.t.ec), al.sets.tgev['loc',], type='b', col='blue' )
abline( lm( al.sets.tgev['loc', ] ~ log(al.sets.mn.t.ec) ), col='blue', lty=2 )
points( log(al.sets.mn.t.ec), al.sets.mgev['loc',], type='b', col='red' )
abline( lm( al.sets.mgev['loc', ] ~ log(al.sets.mn.m.ec) ), col='red', lty=2 )

summary( lm( al.sets.tgev['loc', ] ~ log(al.sets.mn) ))
## -39.7, 7.533. R^2 0.9926, F: 802.3  (lambda estimate: 0.132)
summary( lm( al.sets.tgev['loc', ] ~ log(al.sets.mn.t.ec) ))
## -31.6, 7.197. R^2 0.9924, F: 779.9 (lambda estimate: 0.139)

## Correcting for edge effects increases the estimate for lambda as in
## Altschul & Gish
## But the lambda estimate is still far too small. So something else is going on here.

plot( al.sets.n, 1 / al.sets.tgev['scale',], type='l', col='blue' )
segments( al.sets.n, 1 / (al.sets.tgev['scale',] - al.sets.tgev.std['scale',]), al.sets.n, 1 / (al.sets.tgev['scale',] + al.sets.tgev.std['scale',]), col='blue')
points( al.sets.n, 1 / al.sets.mgev['scale',], type='l', col='red' )
segments( al.sets.n, 1 / (al.sets.mgev['scale',] - al.sets.mgev.std['scale',]), al.sets.n, 1 / (al.sets.mgev['scale',] + al.sets.mgev.std['scale',]), col='red')

## For both sets I have very strong length effects; lambda (1/scale) decreases with length, but seems to be approaching an asymptote with
## longer sequences. This is not unsimilar to the published effects. From the publication:

## alt for Altschul
## nm is already log transformed
alt.nm <- seq(10.5, 16.0, 0.5)
alt.nm.ec <- c(10.25, 10.78, 11.30, 11.83, 12.35, 12.87, 13.39, 13.91, 14.43, 14.94, 15.45, 15.96)
alt.mu <- c(26.45, 28.31, 30.21, 32.04, 33.92, 35.94, 37.84, 39.75, 41.71, 43.54, 45.53, 47.32)
alt.lambda <- c(0.298, 0.286, 0.282, 0.275, 0.279, 0.273, 0.272, 0.275, 0.268, 0.271, 0.267, 0.270)

plot(alt.nm, alt.mu, type='b', col='blue')
points(alt.nm.ec, alt.mu, type='b', col='red')
## these both give much straighter lines than what I get
abline( lm( alt.mu ~ alt.nm ), col='blue', lty=2 )
abline( lm( alt.mu ~ alt.nm.ec ), col='red', lty=2 )

summary( lm( alt.mu ~ alt.nm ))
## -13.67, 3.81, R^2=0.9999
## F: 1.72e05
summary( lm( alt.mu ~ alt.nm.ec ))
## -11.32, 3.67, R^2=0.9999
## F: 9.91e04

## so the non-corrected data gives a better fit..
## that is the same that we have if I am not mistaken.

## We should also do the same for random sequences to determine if we get similar results to Altshcul et al.
## I counted nucleotide frequencies for the full set of

nuc.counts <- read.table("../family_members/introns_nuc_counts.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(nuc.counts) <- c('nuc', 'count')

## remove ambiguity symbols
b <- nuc.counts[,'nuc'] %in% c('A', 'C', 'T', 'G')
nuc.freq <- nuc.counts[ b, 'count'] / sum( nuc.counts[b, 'count'])
names(nuc.freq) <- nuc.counts[b,'nuc']

## then we can simply use sample() to geneerate random sequences
## eg:

paste( sample( names(nuc.freq), size=20, replace=TRUE, prob=nuc.freq ), collapse="" )

align.random <- function( align.length, nuc.freq, n,
                         min.width=15, min.score=20, sm=sub.matrix, gap=c(-10, -2)){
    ## generate pairs of sequences for convenience:
    a.seq <- sapply( 1:n, function(x){ paste( sample( names(nuc.freq), size=20, replace=TRUE, prob=nuc.freq ), collapse="" ) })
    b.seq <- sapply( 1:n, function(x){ paste( sample( names(nuc.freq), size=20, replace=TRUE, prob=nuc.freq ), collapse="" ) })

    aligns <- mclapply( 1:n, function(i){
        local.aligns( a.seq[i], b.seq[i], sm$offset, sm$size, sm$sm,
                      gap=gap, min.width=min.width, min.score=min.score )
    }, mc.cores=20)
    ## and then let us collect some of the information about the aligns.
    al.scores <- sapply(aligns, function(x){ x$pos[1,'score'] })
}
