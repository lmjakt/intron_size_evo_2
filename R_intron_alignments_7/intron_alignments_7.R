## Align all introns from two closely related species and
## consider the overall similarity.

## perform global alignments;

library(parallel)
library(entropy)

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")
source("../R_intron_alignments_2/functions.R")
source("functions.R")
source("~/R/general_functions.R")

## plotting dimensions
## these would be better to put into a single file that is sourced
## by every R file in order to ensure consistency.
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

max.l <- 1e8

## read the intron orthology in.
## this makes of a read.ortho() function that has been copied around
## the different 'functions.R' files. It should be refactored into
## a single source. Better yet would be to export the orthology as
## an RDS file.

orth <- read.ortho()

## add supplementary including local alignment scores at intron positions
orth.sc <- readRDS("../R_intron_orthology/intron_orth_lscore.rds")
## we actually only want to use the 'danio rerio' part here:
orth.sc <- orth.sc[['danio rerio']]

## modify the col names to fit with the orth tables
for(i in 1:length(orth.sc))
    colnames(orth.sc[[i]]) <- sub(" ", ".",  colnames(orth.sc[[i]]), fixed=TRUE)

all( orth.sc$i == orth$i, na.rm=TRUE )  ## TRUE..
all(colnames(orth$i) == colnames(orth.sc$i)) ## TRUE

## Do alignments of introns from two species against each other. For this
## it makes sense for the species to be closely related, and for the introns to
## not be too big.
##
## Looking at the nj-tree we built before, obvious alternatives are:
##
## Poecilia mexicana  vs. Poecilia formosa
## Maylandia zebra vs. Astatotilapia calliptera
##
## Lets have a look at the distributions of these:

sp1.1 <- 'poecilia.mexicana'
sp1.2 <- 'poecilia.formosa'
sp2.1 <- 'maylandia.zebra'
sp2.2 <- 'astatotilapia.calliptera'

par(mfrow=c(2,2))
with(orth, hist( log2(l[,sp1.1]), main=to.sp(sp1.1) ))
with(orth, hist( log2(l[,sp1.2]), main=to.sp(sp1.2) ))
with(orth, hist( log2(l[,sp2.1]), main=to.sp(sp2.1) ))
with(orth, hist( log2(l[,sp2.2]), main=to.sp(sp2.2) ))

## The Poecilias have very nice bimodality; in the cichlids it
## is more difficult to discern.
sum( orth$l[,sp1.1] > 256, na.rm=TRUE ) ## 29971
sum( orth$l[,sp1.2] > 256, na.rm=TRUE ) ## 31681
sum( orth$l[,sp2.1] > 256, na.rm=TRUE ) ## 31286
sum( orth$l[,sp2.2] > 256, na.rm=TRUE ) ## 30720

## so in that respect quite similar:
par(mfrow=c(1,1))
with(orth, plot(log2(l[,sp1.1]), log2(l[,sp1.2]), xlab=to.sp(sp1.1), ylab=to.sp(sp1.2), cex=0.5))
with(orth, plot(log2(l[,sp2.1]), log2(l[,sp2.2]), xlab=to.sp(sp2.1), ylab=to.sp(sp2.2), cex=0.5))

with(orth, plot(l[,sp1.1], l[,sp1.2], xlab=to.sp(sp1.1), ylab=to.sp(sp1.2), cex=0.5))
with(orth, plot(l[,sp2.1], l[,sp2.2], xlab=to.sp(sp2.1), ylab=to.sp(sp2.2), cex=0.5))


## for the Poecilias (sp1) we have a very nice indication of an insertion of
## a repeat element into Poecila mexicana

## to see how good the orthology is we can do:

lm.1 <- with(orth, lm(log2(l[,sp1.1]) ~ log2(l[,sp1.2]))) ## R^2 = 0.8809
lm.2 <- with(orth, lm(log2(l[,sp2.1]) ~ log2(l[,sp2.2]))) ## R^2 = 0.9452

hist( lm.1$residuals, breaks=100 )
hist( lm.2$residuals, breaks=100 )

sum( abs(lm.1$residuals) < 0.5 ) / length(lm.1$residuals)  ## 0.92
sum( abs(lm.1$residuals) < 0.25 ) / length(lm.1$residuals)  ## 0.89
sum( abs(lm.1$residuals) < 0.125 ) / length(lm.1$residuals)  ## 0.64

sum( abs(lm.2$residuals) < 0.5 ) / length(lm.2$residuals)  ## 0.957
sum( abs(lm.2$residuals) < 0.25 ) / length(lm.2$residuals)  ## 0.937
sum( abs(lm.2$residuals) < 0.125 ) / length(lm.2$residuals)  ## 0.894

## so it looks like we may have a better chance with the maylandia zebra
## and the astatotilapia.
## The residuals above may be caused by differences in the quality of the
## intron orthology or by real changes in intron length.
## We can use orth.sc to look at this:

resid.orth.score <- function(sp1, sp2){
    ## assume that the difference in length is due to
    ## NAs
    l1 <- log2( orth$l[,sp1] )
    l2 <- log2( orth$l[,sp2] )
    b <- !is.na(l1) & !is.na(l2)
    l1 <- l1[b]
    l2 <- l2[b]
    mod <- lm( l1 ~ l2 )
    par(mfrow=c(1,3))
    plot( orth.sc$score[b, sp1], abs( mod$residuals ), cex=0.5 )
    plot( orth.sc$score[b, sp2], abs( mod$residuals ), cex=0.5 )
    plot( pmin(orth.sc$score[b, sp1], orth.sc$score[b, sp2]),
         abs( mod$residuals ), cex=0.5 )
}

resid.orth.score(sp1.1, sp1.2)
resid.orth.score(sp2.1, sp2.2)
## there is very little correlation here; however if we do models of these
## then there is a very significant relationship.

## Let us do alignments for both sets; we can then afterwards load in the
## the exon alignments (as in R_intron_alignments_summary) and see whether
## the introns that diverge in size appear to be non-orthologous
## (could be because of paralogy and so on, but..)

intron.f <-  read.table("../R_172_genomes/intron_files.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
tmp <- rownames(intron.f)
intron.f <- intron.f[,1]
names(intron.f) <- tmp
rm(tmp)

## collect sequences for a number of species
collect.intron.s <- function(sp){
    tapply(1:nrow(orth$fam), orth$fam, function(i){
    introns <- read.exons( intron.f[ orth$fam[i[1],1] ] )
    seqs.id <- sapply( introns, function(x){ x$id })
    seqs <- matrix(nrow=length(i), ncol=length(sp))
    colnames(seqs) <- sp
    for(j in 1:length(sp)){
            gid <- orth$id[ i[1], sp[j] ]
            if(is.na(gid)){
                seqs[,j] <- NA
                next
            }
            k <- which(seqs.id == gid)
            if(!length(k)){
                seqs[,j] <- NA
                next
            }
            if(length(k) != 1)
                stop(paste(c("k is longer than one for", gid, i), collapse=" "))
            seqs[,j] <- introns[[k]]$e[ orth$i[ i, sp[j] ] ]
    }
    data.frame( i, seqs, stringsAsFactors=FALSE )
    })
}

sp.l <- c(sp1.1, sp1.2, sp2.1, sp2.2)
introns.l <- collect.intron.s( sp.l )

## that seems to be ok.. We can now try to do some alignments with it.

sub.matrix <- make.sub.matrix()
gap <- c(-10L, -1L)

########### Run a few tests and check that things are working as we expect
###########
tmp <- align.seqs.mt( introns.l[[1]][1,sp1.1], introns.l[[1]][1,sp1.2],
                      al.offset=sub.matrix$offset, al.size=sub.matrix$size,
                      sub.matrix=sub.matrix$sm, gap=gap, tgaps.free=FALSE,
                      thread.n=1 )

al.cols <- c('white', 'grey', hsvScale(1:4, sat=1, val=0.75, max.v=5))
names(al.cols) <- c('-', 'N', 'A', 'C', 'G', 'T')

plot.new()
plot.window(xlim=c(0, nchar(tmp[[1]]$seq[1])), ylim=c(0,5))
draw.aligns( tmp[[1]], 4, 2, 1, cols=al.cols )


## I get segmentation errors with sequences that are too long
## this may be related to not using size_t where I should have
## or simply malloc failing.
## hence I am limiting the alignment space to 1e8
align.introns <- function( x, sp1, sp2, max.l=1e8 ){
    al.int <- function(row){
        s1 <- x[row, sp1]
        s2 <- x[row, sp2]
        if(is.na(s1) || is.na(s2))
            return(NA)
        if(as.double(nchar(s1)) * as.double(nchar(s2)) > max.l)
            return(NA)
        align.seqs.mt( s1, s2, al.offset=sub.matrix$offset, al.size=sub.matrix$size,
                       sub.matrix=sub.matrix$sm, gap=gap,  tgaps.free=FALSE, thread.n=1 )[[1]]
    }
    cat( x[,'i'], '\n' )
    if( !all(c(sp1, sp2) %in% colnames(x)))
        return(NA)
    list(i=x[,'i'], al=lapply(1:nrow(x), al.int))
}

tmp <- align.introns( introns.l[[1]], sp2.1, sp2.2 )

plot.new()
plot.window(xlim=c(0, nchar(tmp$al[[3]]$seq[1])), ylim=c(0,5))
draw.aligns( tmp$al[[3]], 4, 2, 1, cols=al.cols )
### That all looks OK. We should now be able to do:

sp1.al <- mclapply( introns.l, align.introns, sp1=sp1.1, sp2=sp1.2, mc.cores=10 )
sp2.al <- mclapply( introns.l, align.introns, sp1=sp2.1, sp2=sp2.2, mc.cores=10 )

saveRDS(sp1.al, 'sp1_al.rds')
saveRDS(sp2.al, 'sp2_al.rds')


### let us have a look at these.

for(i in 1:length(sp1.al)){
    for(j in 1:length(sp1.al[[i]]$al)){
        if(is.na(sp1.al[i]$al[[j]]))
            next
        n <- nchar( sp1.al[[i]]$al[[j]]$seq[1] )
        plot.new()
        plot.window(xlim=c(0, n), ylim=c(0,5))
        draw.aligns( sp1.al[[i]]$al[[j]], 4, 2, 1, cols=al.cols )
        input <- readline(paste(i, j, sp1.al[[i]]$i[j], ":"))
    }
}

for(i in 1:length(sp2.al)){
    for(j in 1:length(sp2.al[[i]]$al)){
        if(is.na(sp2.al[[i]]$al[j]))
            next
        n <- nchar( sp2.al[[i]]$al[[j]]$seq[1] )
        plot.new()
        plot.window(xlim=c(0, n), ylim=c(0,5))
        draw.aligns( sp2.al[[i]]$al[[j]], 4, 2, 1, cols=al.cols )
        axis(1)
        input <- readline(paste(i, j, sp2.al[[i]]$i[j], ":"))
    }
}

### these data structures look very good indeed.
### We now want to extract the alignment stats for these
### into tables that have the same structure as the orth
### tables.

sp1.al.stats <- matrix( nrow=nrow(orth$l), ncol=2 + length(sp1.al[[1]]$al[[1]]$stats) )
colnames(sp1.al.stats) <- c('i', 'j', names( sp1.al[[1]]$al[[1]]$stats ))
sp2.al.stats <- sp1.al.stats

for(i in 1:length(sp1.al)){
    for(j in 1:length( sp1.al[[i]]$i )){
        if(!is.na(sp1.al[[i]]$al[j]))
            sp1.al.stats[ sp1.al[[i]]$i[j], ] <-  c(i, j, sp1.al[[i]]$al[[j]]$stats)
        if(!is.na(sp2.al[[i]]$al[j]))
            sp2.al.stats[ sp2.al[[i]]$i[j], ] <-  c(i, j, sp2.al[[i]]$al[[j]]$stats)
    }
}

## plot percentage matching for non-gapped positions

plot( log2( orth$l[,sp1.1] ), sp1.al.stats[ ,'match'] / sp1.al.stats[ ,'al_n'], cex=0.5 )

plot( log2( orth$l[,sp1.2] ), sp1.al.stats[ ,'match'] / sp1.al.stats[ ,'al_n'], cex=0.5,
     ylim=c(0.92, 1))

plot( log2( orth$l[,sp2.1] ), sp2.al.stats[ ,'match'] / sp2.al.stats[ ,'al_n'], cex=0.5,
     ylim=c(0.90, 1))

b1 <- !is.na( sp2.al.stats[,'al_n'] ) & (sp2.al.stats[,'match'] / sp2.al.stats[,'al_n']) > 0.9 

b2 <- !is.na( sp2.al.stats[,'al_n'] ) &
             rowSums(sp2.al.stats[,c('a_gap', 'b_gap')]) / sp2.al.stats[,'length'] < 0.05 

plot( log2(orth$l[b1, sp2.1]), sp2.al.stats[b1 ,'match'] / sp2.al.stats[b1 ,'al_n'], cex=0.5 )
plot( log2(orth$l[b2, sp2.1]), sp2.al.stats[b2 ,'match'] / sp2.al.stats[b2 ,'al_n'], cex=0.5 )


for(i in unique( sp2.al.stats[b2 & !b1, 'i'])){
    for(j in 1:length(sp2.al[[i]]$al)){
        if(is.na(sp2.al[[i]]$al[j]))
            next
        n <- nchar( sp2.al[[i]]$al[[j]]$seq[1] )
        plot.new()
        plot.window(xlim=c(0, n), ylim=c(0,5))
        draw.aligns( sp2.al[[i]]$al[[j]], 4, 2, 1, cols=al.cols )
        axis(1)
        input <- readline(paste(i, j, sp2.al[[i]]$i[j], ":"))
    }
}

tmp <- discretize2d( log2( orth$l[b,sp2.1] ), sp2.al.stats[b ,'match'] / sp2.al.stats[b ,'al_n'],
                    numBins1=100, numBins2=100 )

image(tmp)
image(log2(tmp))

image(t(scale(t(tmp))))
image(log2(1 + t(scale(t(tmp)))))

## difficult to make any conclusions from these simple plots. Need to consider grouping the
## sizes and to plot against tel.var[,'min'] instead.

## whatever we do we also want to have some teleost variance parameters..
sp.class <- readRDS( "../R_trees_distances/sp_class_3.rds")
class.col <-  readRDS( "../R_trees_distances/class_col_3.rds")
rownames(sp.class) <- sub("_", ".", rownames(sp.class))

sp.col <- rep(class.col['others'], nrow(sp.class))
for( cl in c('teleostei', 'mammalia', 'eutheria', 'sauria', 'aves'))
    sp.col[ sp.class[,cl] ] <- class.col[cl]
names(sp.col) <- rownames(sp.class)

sp.class <- sp.class[ match( colnames(orth$l), rownames(sp.class) ), ]
all(rownames(sp.class) == colnames(orth$l))  ## TRUE

var.par <- function(x){
    c('n'=sum(!is.na(x)), 'mean'=mean(x, na.rm=TRUE), 'sd'=sd(x, na.rm=TRUE),
      'min'=min(x, na.rm=TRUE), 'max'=max(x, na.rm=TRUE), quantile(x, probs=seq(0,1,0.1), na.rm=TRUE))
}

tel.var <- t(apply( log2( orth$l[ , sp.class[,'teleostei']]), 1, var.par ))

b1 <- !is.na(sp1.al.stats[,'al_n']) &
    rowSums(sp1.al.stats[,c('a_gap', 'b_gap')]) / sp1.al.stats[,'length'] < 0.05 

sp1.pid.1 <- tapply( which(b1), as.integer(log2(orth$l[b1,sp1.1])), function(i){
    sp1.al.stats[i,'match'] / sp1.al.stats[i,'al_n'] })

sp1.pid.2 <- tapply( which(b1), as.integer(tel.var[b1,'min']), function(i){
    sp1.al.stats[i,'match'] / sp1.al.stats[i,'al_n'] })

sp1.pid.3 <- tapply( which(b1), as.integer(tel.var[b1,'50%']), function(i){
    sp1.al.stats[i,'match'] / sp1.al.stats[i,'al_n'] })

sp1.pid.4 <- tapply( which(b1), as.integer(tel.var[b1,'20%']), function(i){
    sp1.al.stats[i,'match'] / sp1.al.stats[i,'al_n'] })

sp1.pid.5 <- tapply( which(b1), as.integer(tel.var[b1,'sd'] * 5), function(i){
    sp1.al.stats[i,'match'] / sp1.al.stats[i,'al_n'] })


plot( as.integer(names(sp1.pid.1)), sapply(sp1.pid.1, mean) )
plot( as.integer(names(sp1.pid.2)), sapply(sp1.pid.2, mean) )
plot( as.integer(names(sp1.pid.3)), sapply(sp1.pid.3, mean) )
plot( as.integer(names(sp1.pid.4)), sapply(sp1.pid.4, mean) )
plot( as.integer(names(sp1.pid.5)), sapply(sp1.pid.5, mean) )


quantile.l <- function(v, l){
    o <- order(v)
    o.l <- seq(1, length(o), as.integer(length(o) / l))
    if(max(o.l) != length(o))
        o.l <- c(o.l, length(o))
    q.i <- rep(0, length(o))
    for(i in 2:length(o.l))
        q.i[ o[ o.l[i-1]:o.l[i] ] ] <- i-1
    q.i
}

sd.associate <- function(al.stats, b, l, par='sd'){
    sd <- tel.var[b, par]
    sim <- al.stats[b,'match'] / al.stats[b,'al_n']
    sd.l <- quantile.l( sd, l )
    tapply( sim, sd.l, eval )
}

l.associate <- function(ind.v, dep.v, l, b=rep(TRUE, length(ind.v))){
    ind.v <- ind.v[b]
    dep.v <- dep.v[b]
    ind.l <- quantile.l( ind.v, l )
    d <- tapply( dep.v, ind.l, eval )
    i <- tapply( ind.v, ind.l, eval )
    list(d=d, i=i)
}

tmp <- sd.associate( sp1.al.stats, b1, 20 )
plot( as.numeric(names(tmp)), sapply(tmp, mean))

b2 <- !is.na( sp2.al.stats[,'al_n'] ) &
    rowSums(sp2.al.stats[,c('a_gap', 'b_gap')]) / sp2.al.stats[,'length'] < 0.05 

tmp <- sd.associate( sp2.al.stats, b2, 20 )
plot( as.numeric(names(tmp)), sapply(tmp, mean))

tmp <- sd.associate( sp2.al.stats, b2, 50, par='sd' )
plot( as.numeric(names(tmp)), sapply(tmp, mean))

summary(lm( as.numeric(names(tmp[-length(tmp)])) ~ sapply(tmp[-length(tmp)], mean) ))

tmp <- sd.associate( sp2.al.stats, b2, 50, par='50%' )
plot( as.numeric(names(tmp)), sapply(tmp, mean))

summary(lm( as.numeric(names(tmp[-length(tmp)])) ~ sapply(tmp[-length(tmp)], mean) ))


tmp <- l.associate( pmin( orth$l[b1, sp1.1], orth$l[b1, sp1.2]),
                   sp1.al.stats[b1,'match'] / sp1.al.stats[b1,'al_n'],
                   20 )
plot( log2(sapply(tmp$i, mean)), sapply(tmp$d, mean))

tmp <- l.associate( pmax( orth$l[b1, sp1.1], orth$l[b1, sp1.2]),
                   sp1.al.stats[b1,'match'] / sp1.al.stats[b1,'al_n'],
                   20 )
plot( log2(sapply(tmp$i, mean)), sapply(tmp$d, mean))

tmp <- l.associate( tel.var[b1, 'sd'],
                   sp1.al.stats[b1,'match'] / sp1.al.stats[b1,'al_n'],
                   20 )
with(tmp, plot( sapply(i[-length(i)], mean), sapply(d[-length(d)], mean)) )


tmp <- l.associate( pmin( orth$l[b2, sp2.1], orth$l[b2, sp2.2]),
                   sp2.al.stats[b2,'match'] / sp2.al.stats[b2,'al_n'],
                   20 )
plot( log2(sapply(tmp$i, mean)), sapply(tmp$d, mean))

tmp <- l.associate( tel.var[b2, 'sd'],
                   sp2.al.stats[b2,'match'] / sp2.al.stats[b2,'al_n'],
                   20 )
with(tmp, plot( sapply(i[-length(i)], mean), sapply(d[-length(d)], mean)) )

tmp <- l.associate( tel.var[b2, 'min'],
                   sp2.al.stats[b2,'match'] / sp2.al.stats[b2,'al_n'],
                   20 )
with(tmp, plot( sapply(i[-length(i)], mean), sapply(d[-length(d)], mean)) )

tmp <- l.associate( tel.var[b2, '50%'],
                   sp2.al.stats[b2,'match'] / sp2.al.stats[b2,'al_n'],
                   20 )
with(tmp, plot( sapply(i[-length(i)], mean), sapply(d[-length(d)], mean)) )



## hmm
plot( tel.var[,'50%'], tel.var[,'sd'], cex=0.5, col=rgb(0, 0, 0, 0.2) )
plot( tel.var[b2,'50%'], tel.var[b2,'sd'], cex=0.5, col=rgb(0, 0, 0, 0.2) )
plot( tel.var[b1,'50%'], tel.var[b1,'sd'], cex=0.25, col=rgb(0, 0, 0, 0.2) )
