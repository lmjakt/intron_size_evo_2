## run alignments of long introns in parallel

library(parallel)

dyn.load("~/R/exon_aligneR/src/exon_aligneR.so")
source("~/R/exon_aligneR/functions.R")


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

## take a list of sequence pairs
p.id <- function(x){
    y <- lapply( x, strsplit, split='' )
    p.id <- sapply( y, function(z){
        sum( z[[1]] == z[[2]] ) / length(z[[1]])
    })
    p.id
}

max.al.scores <- function(al){
    m.scores <- sapply(al$al, function(x){
        ifelse( is.null( nrow(x$pos) ), 0, x$pos[1,'score'] )
    })
    names(m.scores) <- al$sp
    m.scores
}    

plot.alignments <- function(al, col.pid=FALSE){
    ## find out how many of the species have alignments
    ## and order by the highest scoring alignment
    m.scores <- sapply(al$al, function(x){
        ifelse( is.null( nrow(x$pos) ), 0, x$pos[1,'score'] )
    })
    o <- order(m.scores, decreasing=FALSE)
    coords <- sapply(1:length(o), function(i){
        j <- o[i]
        m <- matrix(nrow=0, ncol=4)
        colnames(m) <- c('a_beg', 'a_end', 'score', 'y')
        if(!is.null(al$al[[j]]$pos))
            m <- cbind( al$al[[j]]$pos[,c('a_beg', 'a_end', 'score'), drop=FALSE],
                       'y'=rep(i, nrow(al$al[[j]]$pos)) )
        m
    })
    coords <- do.call(rbind, coords)
    ## and then we can simply call rect with the appropriate parameters
    plot.new()
    y.min <- sum( m.scores == 0 ) + 1
    plot.window(xlim=c(0, al$l1), ylim=c(y.min, 2 + length(m.scores)) )
    cv <- coords[,'score'] / max(m.scores)
    if(col.pid)
        cv <- unlist(sapply( o, function(i){ p.id( al$al[[i]]$seq ) }))
    cols <- hsv(1, s=cv, v=0.8)
    rect( coords[,'a_beg'], coords[,'y'], coords[,'a_end'], coords[,'y']+0.8,
         col=cols )
    b <- m.scores[o] > 0
    axis(2, at=y.min:length(al$sp) + 0.5, labels=al$sp[o[y.min:length(al$sp)]], las=2)
    axis(1)
    invisible( list( 'coord'=coords, 'm.scores'=m.scores[o], 'sp'=al$sp[o], 'l2'=al$l2[o], 'o'=o ))
}

sub.matrix <- make.sub.matrix()

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

tel.var <-  read.table("../R_172_genomes/dr_teleost_var.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)

intron.f <-  read.table("../R_172_genomes/intron_files.txt", header=TRUE, sep="\t",
                       stringsAsFactors=FALSE)
tmp <- rownames(intron.f)
intron.f <- intron.f[,1]
names(intron.f) <- tmp
rm(tmp)

## we can test this with o[1088]
## --> 55047 (COUPTF)

al.test <- align.introns( 55047 )
plot.alignments(al.test)

## ok, that seems to work..
b <- tel.var[,'min'] > 10 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10  ## 1255 introns;
## a bit too many, but, let us take the top 100 with the least
## amount of variance
b.i <- which(b)
o <- order( tel.var[b,'sd'] )

system.time(
    tmp <- align.introns( b.i[o[1]], max.l=1e8 )
)
##    user  system elapsed 
## 108.000   3.186 111.196 

par(mar=c(4.1, 12.1, 4.1, 2.1))
plot.alignments(tmp)

tmp.ms <- max.al.scores( tmp )

## default of max.l=1e7 is too small; can't handle 3000x3000... 
al.top.400 <- mclapply( b.i[ o[1:400] ], align.introns, max.l=1e8,
                       mc.cores=14 )

al.top.28 <- mclapply( b.i[ o[1:28] ], align.introns, max.l=1e8,
                      gap=c(-10, -2),
                      mc.cores=14 )

## lets have a look at the alignments:
par(mar=c(4.1, 12.1, 4.1, 2.1))
for(i in 1:length(al.top.28)){
    plot.alignments( al.top.28[[i]], col.pid=TRUE )
    inpt <- readline(paste(orth.fam[ b.i[i] ], "next: "))
}


al.ctl.28 <- mclapply( b2.i[ 1:28 ], align.introns, max.l=1e8,
                      gap=c(-10, -2),
                      mc.cores=14 )

par(mar=c(4.1, 12.1, 4.1, 2.1))
for(i in 1:length(al.ctl.28)){
    plot.alignments( al.ctl.28[[i]], col.pid=TRUE )
    inpt <- readline(paste(orth.fam[ b2.i[i] ], "next: "))
}


## do
max.al <- 1e8
## numbers in tel.var are log2 transformed. Numbers in orth.l are
## linear..
b <- tel.var[,'min'] > 10 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10  ## 1255 introns;
sum(b)
## 1255
b.i <- which(b)
o <- order( tel.var[b,'sd'] )
b2 <- tel.var[,'X50.'] < 8 & orth.l[,'danio.rerio'] > 1024 & tel.var[,'n'] > 10
sum(b2)
## 3031
## don't redoe the sampling by accidnent
## b2.i <- sample(which(b2))

al.top.500 <- mclapply( b.i[ o[1:500] ], align.introns, max.l=1e8,
                      gap=c(-10, -2),
                      mc.cores=14 )
al.ctl.500 <- mclapply( b2.i[ 1:500 ], align.introns, max.l=1e8,
                      gap=c(-10, -2),
                      mc.cores=14 )

saveRDS(al.top.500, "al_top_500.rds")
saveRDS(al.ctl.500, "al_ctl_500.rds")

plot.al <- function(aligns, b.i, ...){
    par(mar=c(4.1, 12.1, 4.1, 2.1))
    for(i in 1:length(aligns)){
        plot.alignments( aligns[[i]], ... )
        inpt <- readline(paste(orth.fam[ b.i[i] ], "next: "))
        if(inpt == 'q')
            break
    }
}
    
plot.al( al.top.500, b.i, col.pid=TRUE )
plot.al( al.ctl.500, b.i, col.pid=TRUE )



al.top.500.ms <- unlist( lapply( al.top.500, max.al.scores ))
al.ctl.500.ms <- unlist( lapply( al.ctl.500, max.al.scores ))

hist(log2(1 + al.top.500.ms))
hist(log2(1 + al.ctl.500.ms))

al.top.500.l1 <- sapply( al.top.500, function(x){ x$l1 })
al.ctl.500.l1 <- sapply( al.ctl.500, function(x){ x$l1 })

hist(log2(al.top.500.l1))
hist(log2(al.ctl.500.l1))

al.top.500.l2 <- unlist(lapply( al.top.500, function(x){ x$l2 }))
al.ctl.500.l2 <- unlist(lapply( al.ctl.500, function(x){ x$l2 }))

hist(log2(al.top.500.l2))
hist(log2(al.ctl.500.l2))

## compare the scores:
tmp.h <- hist( log2(1 + c(al.top.500.ms, al.ctl.500.ms)), breaks=30 )
top.ms.h <- hist(log2(1 + al.top.500.ms), breaks=tmp.h$breaks)
ctl.ms.h <- hist(log2(1 + al.ctl.500.ms), breaks=tmp.h$breaks)

par(mfrow=c(1,1))
plot( top.ms.h$mids, top.ms.h$count, type='l' )
lines( ctl.ms.h$mids, ctl.ms.h$count, type='l', col='red')
## but this is problematic because we definitely have longer alignments
## for the top.500 than the ctl, so the increased scores may simply be a
## reflection of this. I will need some way of controlling for this;
## the simplest may be to align introns to non-orthologous introns
## in the same set; but this would need to be done for both
## sets to determine whether what is observed indicates orthology


## we need some sort of control; where we look at the statistics of
## obtaining alignments of a given score from non-orthologous introns;

dr.int.l <- orth.l[,'danio.rerio']
dr.int.l[ is.na(dr.int.l) ] <- 0
dr.int.l.o <- order( dr.int.l )
dr.int.l.r <- rank( dr.int.l, ties.method='random' )
dr.int.l.or <- (1:length(dr.int.l))[dr.int.l.o]

## then if we look at 55047
dr.int.l[ 55047 ]
## [1] 2547
dr.int.l.r[ 55047 ]
## [1] 46582
dr.int.l.or[ 46582 ]
## [1] 57964
dr.int.l[ dr.int.l.or[ 46582 ] ]
## [1] 2547
sum( dr.int.l == 2547 )
## [1] 7
dr.int.l.or[ -10:10 + 46582 ]
## this contains 55047, 
##  [1] 46104 55450 58680 59216  4511  5261 14939 29021 44222 55047 57964  5803
## [13]  8889 15993 23664 23771 33390 40695 41216   774  7623

align.introns.ctl <- function( i, max.l, gap ){
    ## find a suitable j to make the alignment from
    rank <- dr.int.l.r[ i ]
    rank.r <- -1:1 + rank
    rr.i <- dr.int.l.or[ rank.r ]
    rr.i <- rr.i[!is.na(rr.i) & rr.i != i]
    align.introns(i=i, j=rr.i[1], max.l=max.l, gap=gap)
}

max.al <- 1e8
## numbers in tel.var are log2 transformed. Numbers in orth.l are
## linear..
b <- tel.var[,'min'] > 10 & !is.na(tel.var[,'sd']) & tel.var[,'n'] > 10  ## 1255 introns;
sum(b)
## 1255
b.i <- which(b)
o <- order( tel.var[b,'sd'] )

tmp <- align.introns.ctl( b.i[o[1]], max.l=1e8, gap=c(-10,-2) )
plot.alignments( tmp, col.pid=TRUE )
plot.alignments( al.top.500[[1]], col.pid=TRUE )

al.top.500.ctl <- mclapply( b.i[ o[1:500] ], align.introns.ctl, max.l=1e8,
                      gap=c(-10, -2),
                      mc.cores=14 )
al.ctl.500.ctl <- mclapply( b2.i[ 1:500 ], align.introns.ctl, max.l=1e8,
                      gap=c(-10, -2),
                      mc.cores=14 )

saveRDS(al.top.500.ctl, "al_top_500_ctl.rds")
saveRDS(al.ctl.500.ctl, "al_ctl_500_ctl.rds")

for(i in 1:28){
    dev.set(2)
    plot.alignments( al.top.500[[i]], col.pid=TRUE )
    mtext("top 500")
    dev.set(3)
    plot.alignments( al.top.500.ctl[[i]], col.pid=TRUE )
    mtext("controls")
    input <- readline("next: ")
}

plot.h2 <- function(v1, v2, breaks=50, log.y=TRUE, cols=c('red', 'blue'),
                    lwd=2, lty=1, ...){
    if(!length(v1) && !length(v2))
        return()
    h.all <- hist(c(v1, v2), breaks=breaks, plot=FALSE)
    h1 <- hist(v1, breaks=h.all$breaks, plot=FALSE )
    h2 <- hist(v2, breaks=h.all$breaks, plot=FALSE )
    c1 <- if(log.y) log2(1+h1$counts) else h1$counts
    c2 <- if(log.y) log2(1+h2$counts) else h2$counts
    plot( h1$mids, c1, ylim=range(c(c1,c2)), col=cols[1], lty=lty, lwd=lwd, type='l', ...)
    points( h1$mids, c2, col=cols[2], lty=lty, lwd=lwd, type='l')
}

for(i in 1:500){
    dev.set(2)
    t5 <- do.call( c, lapply( al.top.500[[i]]$al, function(x){ x$pos[,'score'] }) )
    tc <- do.call( c,lapply(al.top.500.ctl[[i]]$al, function(x){ x$pos[,'score'] } ))
    c5 <- do.call( c, lapply( al.ctl.500[[i]]$al, function(x){ x$pos[,'score'] }) )
    cc <- do.call( c,lapply(al.ctl.500.ctl[[i]]$al, function(x){ x$pos[,'score'] } ))
    par(mfrow=c(1,2))
    ks1 <- if(length(t5) && length(tc)) ks.test( t5, tc ) else 1
    ks2 <- if(length(c5) && length(cc)) ks.test( c5, cc ) else 1
    plot.h2( log2(1+t5), log2(1 + tc), main=sprintf("%.2f: %.2e", ks1$statistic, ks1$p.value) )
    plot.h2( log2(1+c5), log2(1 + cc), main=sprintf("%.2f: %.2e", ks2$statistic, ks2$p.value) )
    input <- readline(paste(i, "next: "))
}

## From https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
## We get the following equation:
##
## E = Kmn * exp(-lambda * S)
##
## where
## E is the expected number of alignments for the given score S
## m and n are the length of the sequences
## K and lambda are constants that relate to the scoring system
## and which are not trivial to estimate. Nevertheless, it
## does suggest that
## e = mn * exp(-S)
##
## should give us a a scaled e-value that we can use to compare between
## different alignments. Or we can consider:
##

## these are not real expect values, as we do not know the value of lambda
## and K. But they should allow us to compare the values more directly
## We still need a lambda here, or we will get underflow issues as
## exp(-big number) becomes very small very rapidly.

## Unfortunately this does not seem to work terribly well;
## looking at the differences in the distributions seems to be the best way forward
## we can also consider some sort of q-values for paired distributions
## maybe qq plots of the scores are enough.. 

get.e <- function(al, lambda=0.01){
    m <- al$l1
    n <- al$l2
    mn <- m * n
    if(length(mn) != length(al$al))
        stop("mismatched lengths")
    lapply( 1:length(al$al), function(i){
        if(!is.null( al$al[[i]]$pos))
            cbind(i=i, m=m, n=n[i], score=al$al[[i]]$pos[,'score'],
                  e=mn * exp(-lambda * al$al[[i]]$pos[,'score']))
    })
}

for(i in 1:500){
    t5 <- do.call( rbind, get.e(al.top.500[[i]]) )
    tc <- do.call( rbind, get.e(al.top.500.ctl[[i]]) )
    c5 <- do.call( rbind, get.e(al.ctl.500[[i]]) )
    cc <- do.call( rbind, get.e(al.ctl.500.ctl[[i]]) )
    ##
    dev.set(2)
    par(mfrow=c(1,2))
#    ks1 <- if(length(t5) && length(tc)) ks.test( t5, tc ) else 1
#    ks2 <- if(length(c5) && length(cc)) ks.test( c5, cc ) else 1
    if(!is.null(t5) && !is.null(tc))
        plot.h2( log10(t5[,'e']), log10(tc[,'e']) )
    if(!is.null(c5) && !is.null(cc))
        plot.h2( log10(c5[,'e']), log10(cc[,'e']) )
    input <- readline(paste(i, "next: "))
}

## try q-q plots instead and see what they look like:
for(i in 1:500){
    t5 <- do.call( rbind, get.e(al.top.500[[i]]) )
    tc <- do.call( rbind, get.e(al.top.500.ctl[[i]]) )
    c5 <- do.call( rbind, get.e(al.ctl.500[[i]]) )
    cc <- do.call( rbind, get.e(al.ctl.500.ctl[[i]]) )
    ##
    dev.set(2)
    par(mfrow=c(1,2))
    if(!is.null(t5) && !is.null(tc)){
        l <- min( c(nrow(t5), nrow(tc)) )
        plot( sort(log2(t5[,'score']), decreasing=TRUE)[1:l],
              sort(log2(tc[,'score']), decreasing=TRUE)[1:l],
             main='top 500', cex=0.5)
        abline(0,1,col='red')
    }
    if(!is.null(c5) && !is.null(cc)){
        l <- min( c(nrow(c5), nrow(cc)) )
        plot( sort(log2(c5[1:l,'score']), decreasing=TRUE)[1:l],
              sort(log2(cc[1:l,'score']), decreasing=TRUE)[1:l],
             main='control', cex=0.5)
        abline(0,1,col='red')
    }
    input <- readline(paste(i, "next: "))
}

## normalise values in al1 by al2;
## al2 should be something that approximates the NULL hypothesis
## even if not exaclty that.
qnt.norm <- function(al1, al2){
    ## to get all the scores and associated values we can use
    s1 <- do.call( rbind, get.e( al1 ) )
    s2 <- do.call( rbind, get.e( al2 ) )
    if(is.null(s1) || is.null(s2))
        return(NULL)
    s1.l <- log2(s1[,'score'])
    s2.l <- sort(log2(s2[,'score']))
    s1.r <- rank(s1.l, ties.method='random')
    s1.q <- s1.l - s2.l[ round( nrow(s2) * s1.r / nrow(s1) ) ]
    cbind('qn'=s1.q, s1)
}

## let us have a look at that...
for(i in 1:length(al.top.500)){
    q1 <- qnt.norm( al.top.500[[i]], al.top.500.ctl[[i]] )
    q2 <- qnt.norm( al.ctl.500[[i]], al.ctl.500.ctl[[i]] )
    if(is.null(q1) || is.null(q2))
        next
    plot(1:nrow(q1), sort(q1[,'qn']), xlim=c(1, max(nrow(q1), nrow(q2))),
         ylim=range(c(q1[,'qn'], q2[,'qn'])), col='red' )
    points(1:nrow(q2), sort(q2[,'qn']), col='blue')
    inpt <- readline('next')
}

harvest.qn <- function(al1, al2){
    qn <- lapply(1:length(al1), function(i){
        qnt.norm( al1[[i]], al2[[i]] )
    })
    qn.r <- sapply( qn, function(x){ range(x[,'qn']) })
    list('qn'=qn, 'r'=qn.r)
}

al.top.500.qn <- harvest.qn( al.top.500, al.top.500.ctl )
al.ctl.500.qn <- harvest.qn( al.ctl.500, al.ctl.500.ctl )

plot( al.top.500.qn$r[2,] )
plot( al.ctl.500.qn$r[2,] )

## We can now try to plot all of these on a single page.
collect.qn.points <- function(qn, al){
    points <- lapply(1:length(qn$qn), function(i){
        if(!is.null( qn$qn[[i]] ))
            suppressWarnings( cbind(i, qn$qn[[i]][,'qn']) )
    })
    points <- do.call(rbind, points[ !is.null(points) ])
    ## we then want to get the species for these points
    sp <- unlist(lapply( 1:length(qn$qn), function(i){
        al[[i]]$sp[ qn$qn[[i]][,'i'] ]
    }))
    list('pts'=points, 'sp'=sp)
}

plot.qn.points <- function(qn, min.score=1, cols=sp.col, ...){
    b <- qn$pts[,2] >= min.score
    plot(qn$pts[b,1], qn$pts[b,2], col=cols[ qn$sp[b] ], ...)
}

al.top.500.qp <- collect.qn.points( al.top.500.qn, al.top.500 )
al.ctl.500.qp <- collect.qn.points( al.ctl.500.qn, al.ctl.500 )


## let us get species information..
teleost.b <- read.table( "../R_172_genomes/teleost_b.txt", header=TRUE, sep="\t" )
mammal.b <-  read.table( "../R_172_genomes/mammal_b.txt", header=TRUE, sep="\t" )
sauria.b <-  read.table( "../R_172_genomes/sauria_b.txt", header=TRUE, sep="\t" )

rownames(teleost.b) <- sub(" ", ".", rownames(teleost.b))
rownames(mammal.b) <- sub(" ", ".", rownames(mammal.b))
rownames(sauria.b) <- sub(" ", ".", rownames(sauria.b))

sp.col <- rgb( teleost.b[,1], mammal.b[,1], sauria.b[,1] )
names(sp.col) <- sub(" ", ".", rownames(teleost.b))

dev.set(3)
plot.qn.points( al.top.500.qp, main='top 500', lwd=2 )
dev.set(2)
plot.qn.points( al.ctl.500.qp, main='ctl 500', lwd=2 )
## And that does suggest that there is a real difference between these
## which probably is not due to simply having more sequence to align.
## We could also restrict ourselves to the top alignment per species.
## But, let us leave that for the time being.

## lets look a the individual positions;
## takes an entry for a single intron
harvest.conservation <- function( al, min.score=60, max.l=1e8 ){
    id.count <- matrix(0, nrow=4, ncol=al$l1)
    rownames(id.count) <- c('all', 'teleost', 'mammal', 'sauria')
    sp <- al$sp
    sp.b <- (al$l1 * al$l2) <= max.l
    tel.b <- teleost.b[ sp, 1 ]
    mam.b <- mammal.b[ sp, 1 ]
    sau.b <- sauria.b[ sp, 1 ]
    for(i in 1:length(al$al)){
        ## i is per species
        if(is.null(al$al[[i]]$pos))
            next
        for(j in 1:nrow(al$al[[i]]$pos)){
            if(al$al[[i]]$pos[j,'score'] < min.score )
                break
            s <- strsplit(al$al[[i]]$seq[[j]], '')
            ## the positions in the alignment
            b <- s[[1]] != '-'
            p <- cumsum( b ) + al$al[[i]]$pos[j,'a_beg'] - 1
            b2 <- s[[1]] == s[[2]]
            id.count['all', p[b2] ] <- id.count['all', p[b2] ] + 1
            if(tel.b[i])
                id.count['teleost', p[b2] ] <- id.count['teleost', p[b2] ] + 1
            if(mam.b[i])
                id.count['mammal', p[b2] ] <- id.count['mammal', p[b2] ] + 1
            if(sau.b[i])
                id.count['sauria', p[b2] ] <- id.count['sauria', p[b2] ] + 1
        }
    }
    list('count'=id.count, 'sp'=sp, 'sp.b'=sp.b,
         'tel.b'=tel.b, 'mam.b'=mam.b, 'sau.b'=sau.b)
}

tmp <- harvest.conservation( al.top.500[[1]] )

ms <- 200
for(i in 1:500){
    tmp1 <- harvest.conservation( al.top.500[[i]], min.score=ms )
    tmp1.c <- harvest.conservation( al.top.500.ctl[[i]], min.score=ms )
    tmp2 <- harvest.conservation( al.ctl.500[[i]], min.score=ms )
    tmp2.c <- harvest.conservation( al.ctl.500.ctl[[i]], min.score=ms )
    par(mfrow=c(4,1))
    par(mar=c(2.1, 4.1, 2.1, 2.1))
    if(!sum(tmp1$sp.b))
        next
    plot( 1:ncol(tmp1$count), tmp1$count['all',] / sum(tmp1$sp.b), type='l', col='grey' )
    points( 1:ncol(tmp1$count), tmp1$count['teleost',] / sum(tmp1$sp.b & tmp1$tel.b),
           type='l', col='red' )
    points( 1:ncol(tmp1$count), tmp1$count['mammal',] / sum(tmp1$sp.b & tmp1$mam.b),
           type='l', col='blue' )
    points( 1:ncol(tmp1$count), tmp1$count['sauria',] / sum(tmp1$sp.b & tmp1$sau.b),
           type='l', col='green' )
##
    plot( 1:ncol(tmp1.c$count), tmp1.c$count['all',] / sum(tmp1.c$sp.b), type='l', col='grey' )
    points( 1:ncol(tmp1.c$count), tmp1.c$count['teleost',] / sum(tmp1.c$sp.b & tmp1.c$tel.b),
           type='l', col='red' )
    points( 1:ncol(tmp1.c$count), tmp1.c$count['mammal',] / sum(tmp1.c$sp.b & tmp1.c$mam.b),
           type='l', col='blue' )
    points( 1:ncol(tmp1.c$count), tmp1.c$count['sauria',] / sum(tmp1.c$sp.b & tmp1.c$sau.b),
           type='l', col='green' )
##
    plot( 1:ncol(tmp2$count), tmp2$count['all',] / sum(tmp2$sp.b), type='l', col='grey' )
    points( 1:ncol(tmp2$count), tmp2$count['teleost',] / sum(tmp2$sp.b & tmp2$tel.b),
           type='l', col='red' )
    points( 1:ncol(tmp2$count), tmp2$count['mammal',] / sum(tmp2$sp.b & tmp2$mam.b),
           type='l', col='blue' )
    points( 1:ncol(tmp2$count), tmp2$count['sauria',] / sum(tmp2$sp.b & tmp2$sau.b),
           type='l', col='green' )
##
    plot( 1:ncol(tmp2.c$count), tmp2.c$count['all',] / sum(tmp2.c$sp.b), type='l', col='grey' )
    points( 1:ncol(tmp2.c$count), tmp2.c$count['teleost',] / sum(tmp2.c$sp.b & tmp2.c$tel.b),
           type='l', col='red' )
    points( 1:ncol(tmp2.c$count), tmp2.c$count['mammal',] / sum(tmp2.c$sp.b & tmp2.c$mam.b),
           type='l', col='blue' )
    points( 1:ncol(tmp2.c$count), tmp2.c$count['sauria',] / sum(tmp2.c$sp.b & tmp2.c$sau.b),
           type='l', col='green' )
    inpt <- readline("next: ")
}

### How do we make sure that we are not simply dealing with
### simple repeats here? We do have longer sequences and these
### are more likely to have simple repeat content.
#### think about it tomorrow. Hmm, I can 

extract.scores <- function(x, max.l=max.l){
    b.i <- which(as.double(x$l1) * as.double(x$l2) <= max.l)
    do.call(rbind, lapply(b.i, function(i){
        if(is.null( x$al[[i]]$pos )) return(NULL)
        suppressWarnings( cbind( 'l1'=x$l1, 'l2'=x$l2, x$al[[i]]$pos[,'score'] ))
    }))
}

extract.scores.2 <- function(x, max.l=max.l, sp.b=teleost.b){
    b.i <- which(as.double(x$l1) * as.double(x$l2) <= max.l & sp.b[ x$sp, 1 ])
    do.call(rbind, lapply(b.i, function(i){
        if(is.null( x$al[[i]]$pos )) return(NULL)
        suppressWarnings( cbind( 'l1'=x$l1, 'l2'=x$l2, x$al[[i]]$pos[,'score'] ))
    }))
}

df2double <- function(x){
    cn <- colnames(x)
    x <- matrix( as.double(x), ncol=ncol(x))
    colnames(x) <- cn
    x
}

max.l <- 1e8
top.scores <- do.call(rbind, lapply( al.top.500, extract.scores, max.l=max.l ))
top.ctl.scores <- do.call(rbind, lapply( al.top.500.ctl, extract.scores, max.l=max.l ))

top.scores <- matrix( as.double(top.scores), ncol=ncol(top.scores) )
top.ctl.scores <- matrix( as.double(top.ctl.scores), ncol=ncol(top.ctl.scores) )

## lets do the same for the ctl.500
ctl.scores <- do.call(rbind, lapply( al.ctl.500, extract.scores, max.l=max.l ))
ctl.ctl.scores <- do.call(rbind, lapply( al.ctl.500.ctl, extract.scores, max.l=max.l))
ctl.scores <- matrix( as.double(ctl.scores), ncol=ncol(ctl.scores))
ctl.ctl.scores <- matrix( as.double(ctl.ctl.scores), ncol=ncol(ctl.ctl.scores))

colnames(top.scores) <- c('l1', 'l2', 'score')
colnames(top.ctl.scores) <- colnames(top.scores)
colnames(ctl.scores) <- c('l1', 'l2', 'score')
colnames(ctl.ctl.scores) <- colnames(ctl.scores)

par(mfrow=c(2,2))
hist( log2(top.ctl.scores[,3]) )
hist( log2(top.scores[,3]) )

library(MASS)

top.scores.2dh <- kde2d( log2(top.scores[,'l1'] * top.scores[,'l2']), log2(top.scores[,'score']) )

top.ctl.scores.2dh <- kde2d( log2(top.ctl.scores[,'l1'] * top.ctl.scores[,'l2']),
                            log2(top.ctl.scores[,'score']) )

image(top.scores.2dh$x, top.scores.2dh$y, log2(top.scores.2dh$z))
image(top.ctl.scores.2dh$x, top.ctl.scores.2dh$y, log2(top.ctl.scores.2dh$z))

image( top.scores.2dh$x, top.scores.2dh$y, top.scores.2dh$z - top.ctl.scores.2dh$z)

## the above two-dimensional distributions hint at something usable; almost
## all of the alignments are likely to be non-informative in both distributions, but
## that is what we expect. What we should do instead is to do divided histograms
## make life a little bit simpler
top.scores <- cbind(top.scores, 'mn'=log2( top.scores[,'l1'] * top.scores[,'l2'] ) )
top.ctl.scores <- cbind(top.ctl.scores, 'mn'=log2( top.ctl.scores[,'l1'] * top.ctl.scores[,'l2'] ) )

ctl.scores <- cbind(ctl.scores, 'mn'=log2( ctl.scores[,'l1'] * ctl.scores[,'l2'] ))
ctl.ctl.scores <- cbind(ctl.ctl.scores, 'mn'=log2( ctl.ctl.scores[,'l1'] * ctl.ctl.scores[,'l2'] ))

top.all.mn.h <- hist( c(top.scores[,'mn'], top.ctl.scores[,'mn']) )
ctl.all.mn.h <- hist( c(ctl.scores[,'mn'], ctl.ctl.scores[,'mn']) )

qnt.probs <- seq(0, 1, 0.05 )
top.all.qnt <- quantile( c(top.scores[,'mn'], top.ctl.scores[,'mn']),
                        probs=qnt.probs)
ctl.all.qnt <- quantile( c(ctl.scores[,'mn'], ctl.ctl.scores[,'mn']),
                        probs=qnt.probs)

top.all.r <- range( c(top.scores[,'mn'], top.ctl.scores[,'mn']) )

mn.breaks <- seq(top.all.r[1], top.all.r[2], length.out=50)

top.mn.slices <- colSums( sapply( top.scores[,'mn'], function(x){ x > mn.breaks } ))
top.ctl.mn.slices <- colSums( sapply( top.ctl.scores[,'mn'], function(x){ x > mn.breaks } ))

top.mnq.slices <- colSums( sapply( top.scores[,'mn'],
                                  function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))

top.ctl.mnq.slices <- colSums( sapply( top.ctl.scores[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))

ctl.mnq.slices <- colSums( sapply( ctl.scores[,'mn'], function(x){
    x >= ctl.all.qnt[-length(ctl.all.qnt)] }))
ctl.ctl.mnq.slices <- colSums( sapply( ctl.ctl.scores[,'mn'], function(x){
    x >= ctl.all.qnt[-length(ctl.all.qnt)] }))

## given the small numbers, for extremes it might be better to use the quantiles breaks for this
## for most of the data they are pretty well spaced.
top.all.h <- hist( log2( c(top.scores[,'score'], top.ctl.scores[,'score'])), breaks=50)
top.h <- hist(log2(top.scores[,'score']), breaks=top.all.h$breaks)
ctl.h <- hist(log2(top.ctl.scores[,'score']), breaks=top.all.h$breaks)

ctl.all.h <- hist(log2(c(ctl.scores[,'score'], ctl.ctl.scores[,'score'])), breaks=50)

top.scores.sliced.h <- tapply( log2(top.scores[,'score']), top.mn.slices,
                              hist, breaks=top.all.h$breaks, plot=FALSE)

top.ctl.scores.sliced.h <- tapply( log2(top.ctl.scores[,'score']), top.ctl.mn.slices,
                              hist, breaks=top.all.h$breaks, plot=FALSE)

top.scores.qsliced.h <- tapply( log2(top.scores[,'score']), top.mnq.slices,
                               hist, breaks=top.all.h$breaks, plot=FALSE)

top.ctl.scores.qsliced.h <- tapply( log2(top.ctl.scores[,'score']), top.ctl.mnq.slices,
                                   hist, breaks=top.all.h$breaks, plot=FALSE )

ctl.scores.qsliced.h <- tapply( log2(ctl.scores[,'score']), ctl.mnq.slices,
                               hist, breaks=ctl.all.h$breaks )
ctl.ctl.scores.qsliced.h <- tapply( log2(ctl.ctl.scores[,'score']), ctl.ctl.mnq.slices,
                               hist, breaks=ctl.all.h$breaks )

## should really make a function; extract named.. 
plot.counts <- function(h, dens=FALSE, log=TRUE,
                        cols=hsv(1, 0.8, 1:length(h) / length(h)), ...){
    if(dens)
        counts <- sapply(h, function(x){ x$counts })
    else
        counts <- sapply(h, function(x){ x$counts })
    if(log)
        counts <- log2( min(counts) / 2 + counts )
    mids <- h[[1]]$mids
    plot(mids, counts[,1], type='l', lwd=2, ylim=range(counts[is.finite(counts)], na.rm=TRUE),
         col=cols[1], ...)
    invisible( sapply( 2:ncol(counts[,-1]), function(i){
        lines(mids, counts[,i], lwd=2, col=cols[i])}))
}
        
plot.counts( top.scores.sliced.h, dens=TRUE)                         
plot.counts( top.ctl.scores.sliced.h )                         

pdf("qsliced_score_distributions.pdf", width=7, height=7, title='Score distributions by quantile sequence alignment space')
par(mfrow=c(2,2))
plot.counts( top.scores.qsliced.h, main='Long teleost',
            xlab='log2 score', ylab='log2 n')
plot.counts( top.ctl.scores.qsliced.h, main='Long teleost, control',
            xlab='log2 score', ylab='log2 n')
plot.counts( ctl.scores.qsliced.h, main='Sampled teleost',
            xlab='log2 score', ylab='log2 n')
plot.counts( ctl.ctl.scores.qsliced.h, main='Sampled teleost control',
            xlab='log2 score', ylab='log2 n')
dev.off()


### let us make the same plot but restricting the alignments to different species
### clades. We can actually reuse the quantiles that we have used before. We just
### have to check if three is any particular problem.

top.scores.tel <- do.call(rbind, lapply( al.top.500, extract.scores.2, max.l=max.l, sp.b=teleost.b))
top.scores.mam <- do.call(rbind, lapply( al.top.500, extract.scores.2, max.l=max.l, sp.b=mammal.b))
top.scores.sau <- do.call(rbind, lapply( al.top.500, extract.scores.2, max.l=max.l, sp.b=sauria.b))

top.ctl.scores.tel <- do.call(rbind, lapply( al.top.500.ctl,
                                            extract.scores.2, max.l=max.l, sp.b=teleost.b))
top.ctl.scores.mam <- do.call(rbind, lapply( al.top.500.ctl,
                                            extract.scores.2, max.l=max.l, sp.b=mammal.b))
top.ctl.scores.sau <- do.call(rbind, lapply( al.top.500.ctl,
                                            extract.scores.2, max.l=max.l, sp.b=sauria.b))

ctl.scores.tel <- do.call(rbind, lapply( al.ctl.500, extract.scores.2, max.l=max.l, sp.b=teleost.b))
ctl.scores.mam <- do.call(rbind, lapply( al.ctl.500, extract.scores.2, max.l=max.l, sp.b=mammal.b))
ctl.scores.sau <- do.call(rbind, lapply( al.ctl.500, extract.scores.2, max.l=max.l, sp.b=sauria.b))

ctl.ctl.scores.tel <- do.call(rbind, lapply( al.ctl.500.ctl, extract.scores.2, max.l=max.l, sp.b=teleost.b))
ctl.ctl.scores.mam <- do.call(rbind, lapply( al.ctl.500.ctl, extract.scores.2, max.l=max.l, sp.b=mammal.b))
ctl.ctl.scores.sau <- do.call(rbind, lapply( al.ctl.500.ctl, extract.scores.2, max.l=max.l, sp.b=sauria.b))

## we need to convert these to double matrices...
top.scores.tel <- df2double( top.scores.tel )
top.scores.mam <- df2double( top.scores.mam )
top.scores.sau <- df2double( top.scores.sau )
top.ctl.scores.tel <- df2double( top.ctl.scores.tel )
top.ctl.scores.mam <- df2double( top.ctl.scores.mam )
top.ctl.scores.sau <- df2double( top.ctl.scores.sau )

ctl.scores.tel <- df2double( ctl.scores.tel )
ctl.scores.mam <- df2double( ctl.scores.mam )
ctl.scores.sau <- df2double( ctl.scores.sau )
ctl.ctl.scores.tel <- df2double( ctl.ctl.scores.tel )
ctl.ctl.scores.mam <- df2double( ctl.ctl.scores.mam )
ctl.ctl.scores.sau <- df2double( ctl.ctl.scores.sau )

## this is utterly horrible, I should have made some reasonable lists.. I really should
## redo everything in here.. but lets do something quickly.

colnames(top.scores.tel) <- c('l1', 'l2', 'score')
colnames(top.scores.mam) <- c('l1', 'l2', 'score')
colnames(top.scores.sau) <- c('l1', 'l2', 'score')
colnames(top.ctl.scores.tel) <- c('l1', 'l2', 'score')
colnames(top.ctl.scores.mam) <- c('l1', 'l2', 'score')
colnames(top.ctl.scores.sau) <- c('l1', 'l2', 'score')
colnames(ctl.scores.tel) <- c('l1', 'l2', 'score')
colnames(ctl.scores.mam) <- c('l1', 'l2', 'score')
colnames(ctl.scores.sau) <- c('l1', 'l2', 'score')
colnames(ctl.ctl.scores.tel) <- c('l1', 'l2', 'score')
colnames(ctl.ctl.scores.mam) <- c('l1', 'l2', 'score')
colnames(ctl.ctl.scores.sau) <- c('l1', 'l2', 'score')

add.mn <- function(x){
    cbind(x, 'mn'=log2( x[,'l1'] * x[,'l2'] ))
}

top.scores.tel <- add.mn( top.scores.tel )
top.scores.mam <- add.mn( top.scores.mam )
top.scores.sau <- add.mn( top.scores.sau )
top.ctl.scores.tel <- add.mn( top.ctl.scores.tel )
top.ctl.scores.mam <- add.mn( top.ctl.scores.mam )
top.ctl.scores.sau <- add.mn( top.ctl.scores.sau )
ctl.scores.tel <- add.mn( ctl.scores.tel )
ctl.scores.mam <- add.mn( ctl.scores.mam )
ctl.scores.sau <- add.mn( ctl.scores.sau )
ctl.ctl.scores.tel <- add.mn( ctl.ctl.scores.tel )
ctl.ctl.scores.mam <- add.mn( ctl.ctl.scores.mam )
ctl.ctl.scores.sau <- add.mn( ctl.ctl.scores.sau )

top.mnq.tel.slices <- colSums( sapply(top.scores.tel[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
top.mnq.mam.slices <- colSums( sapply(top.scores.mam[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
top.mnq.sau.slices <- colSums( sapply(top.scores.sau[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))

top.ctl.mnq.tel.slices <- colSums( sapply(top.ctl.scores.tel[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
top.ctl.mnq.mam.slices <- colSums( sapply(top.ctl.scores.mam[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
top.ctl.mnq.sau.slices <- colSums( sapply(top.ctl.scores.sau[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))

ctl.mnq.tel.slices <- colSums( sapply(ctl.scores.tel[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
ctl.mnq.mam.slices <- colSums( sapply(ctl.scores.mam[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
ctl.mnq.sau.slices <- colSums( sapply(ctl.scores.sau[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))

ctl.ctl.mnq.tel.slices <- colSums( sapply(ctl.ctl.scores.tel[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
ctl.ctl.mnq.mam.slices <- colSums( sapply(ctl.ctl.scores.mam[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))
ctl.ctl.mnq.sau.slices <- colSums( sapply(ctl.ctl.scores.sau[,'mn'],
                                      function(x){ x >= top.all.qnt[-length(top.all.qnt)] }))

#### Argh, this is horrid..
top.scores.tel.qsliced.h <- tapply( log2(top.scores.tel[,'score']), top.mnq.tel.slices,
                                   hist, breaks=top.all.h$breaks, plot=FALSE)
top.scores.mam.qsliced.h <- tapply( log2(top.scores.mam[,'score']), top.mnq.mam.slices,
                                   hist, breaks=top.all.h$breaks, plot=FALSE)
top.scores.sau.qsliced.h <- tapply( log2(top.scores.sau[,'score']), top.mnq.sau.slices,
                                   hist, breaks=top.all.h$breaks, plot=FALSE)

top.ctl.scores.tel.qsliced.h <- tapply( log2(top.ctl.scores.tel[,'score']), top.ctl.mnq.tel.slices,
                                   hist, breaks=top.all.h$breaks, plot=FALSE)
top.ctl.scores.mam.qsliced.h <- tapply( log2(top.ctl.scores.mam[,'score']), top.ctl.mnq.mam.slices,
                                   hist, breaks=top.all.h$breaks, plot=FALSE)
top.ctl.scores.sau.qsliced.h <- tapply( log2(top.ctl.scores.sau[,'score']), top.ctl.mnq.sau.slices,
                                   hist, breaks=top.all.h$breaks, plot=FALSE)

ctl.scores.tel.qsliced.h <- tapply( log2(ctl.scores.tel[,'score']), ctl.mnq.tel.slices,
                                   hist, breaks=ctl.all.h$breaks, plot=FALSE)
ctl.scores.mam.qsliced.h <- tapply( log2(ctl.scores.mam[,'score']), ctl.mnq.mam.slices,
                                   hist, breaks=ctl.all.h$breaks, plot=FALSE)
ctl.scores.sau.qsliced.h <- tapply( log2(ctl.scores.sau[,'score']), ctl.mnq.sau.slices,
                                   hist, breaks=ctl.all.h$breaks, plot=FALSE)

ctl.ctl.scores.tel.qsliced.h <- tapply( log2(ctl.ctl.scores.tel[,'score']), ctl.ctl.mnq.tel.slices,
                                   hist, breaks=ctl.all.h$breaks, plot=FALSE)
ctl.ctl.scores.mam.qsliced.h <- tapply( log2(ctl.ctl.scores.mam[,'score']), ctl.ctl.mnq.mam.slices,
                                   hist, breaks=ctl.all.h$breaks, plot=FALSE)
ctl.ctl.scores.sau.qsliced.h <- tapply( log2(ctl.ctl.scores.sau[,'score']), ctl.ctl.mnq.sau.slices,
                                   hist, breaks=ctl.all.h$breaks, plot=FALSE)


## finally we can have a look at these..
par(mfrow=c(2,2))
plot.counts( top.scores.qsliced.h, main='long teleost all')
plot.counts( top.scores.tel.qsliced.h, main='long teleost tel')
plot.counts( top.scores.mam.qsliced.h, main='long teleost mam')
plot.counts( top.scores.sau.qsliced.h, main='long teleost sau')

plot.counts( top.scores.mam.qsliced.h, main='long introns')
plot.counts( top.ctl.scores.mam.qsliced.h, main='long introns ctl')
plot.counts( ctl.scores.mam.qsliced.h, main='ctl introns')
plot.counts( ctl.ctl.scores.mam.qsliced.h, main='ctl introns ctl')

### That looks not so bad. I better consider to rewrite all of this though.. 

#### This is probably good enough for now
#### We could certainly try to extract some examples of aligned regions, this
#### and it would be nice to do lots of other stuff as well. But I think I should
#### leave it here and start writing.

## We don't actually see a very strong effect of the alignment space
## size. But it is very clear that there are more high scoring alignments
## for the sampled teleosts. What we would like to do here is also to
## break these down by teleost / mammal and so on, as we would expect
## to see bigger differences there.

##
## but this looks completely different to the kde2d data above?
## what is going on here? 

## we can harvest details from a set of alignments;
## for each one we want to know how many alignments were attempted
## and the statistics of those that were succsessful.

## experiment a bit to see what has happened

plot( al.top.500[[2]]$l2 * al.top.500[[2]]$l1 )
sum( al.top.500[[1]]$l2 * al.top.500[[1]]$l1 > 1e8 )
## all of these should be possible for me to align;

sum( al.top.500[[1]]$l2 * al.top.500[[1]]$l1 > 0 ) ## 160
range( al.top.500[[1]]$l2 * al.top.500[[1]]$l1 ) ##   312816 21794976
range( al.top.500[[1]]$l2 )  ## 98 6828

par(mar=c(4.1, 12.1, 4.1, 2.1))
tmp <- plot.alignments( al.top.500[[2]], col.pid=TRUE )

