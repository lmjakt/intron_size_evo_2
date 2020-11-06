## Plots and figures showing the distribution of vertebrate intron
## size distributions.

source("../R/functions.R")
source("~/R/general_functions.R")


gene.stats <- read.table("../family_members/orthologue_transcripts/exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)
colnames(gene.stats) <- c('sp', 'db', 'family', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

### for all genes:
all.gene.stats <- read.table("all_genes_exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)

colnames(all.gene.stats) <- c('sp', 'db', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

dim(gene.stats)

## [1] 1063202      10
## about a million rows

dim(all.gene.stats)
## [1] 4,368,399       9
## only about 4 times more? Hmm, 4x6 = 24..
## quite reasonble


barplot(tapply( gene.stats[,'gene'], gene.stats[,'sp'], length ), las=2)

par(mar=c(4.1, 12.1, 4.1, 2.1))
barplot(tapply( all.gene.stats[,'gene'], all.gene.stats[,'sp'], length ), las=2,
        cex.names=0.5, horiz=T)

sort(tapply( all.gene.stats[,'gene'], all.gene.stats[,'sp'], length ))
## too many for mus_musculus and homo_sapiens perhaps...
## these may not be all coding genes. 
## lets have a look ath chromosome names
hum.b <- all.gene.stats[,'sp'] == 'homo_sapiens'
mus.b <- all.gene.stats[,'sp'] == 'mus_musculus'

non.chr.b <- grepl( "\\.[0-9]+$", all.gene.stats[,'chr'] )
## that doesn't really help us at very few are both human and non-chromosomal
## presumably these numbers represent the type of gene. I should have included the type
## of gene in the data. But we can anyway have a look at these first.

## Unfortunately this is only for the orthologous transcript set; we actually want
## numbers for the full set of genes. I have those at work, but do not have time to go
## in today. Will need to go tomorrow to copy the relevant data.

## Nevertheless we need to show that the distributions here are the same as for others.
exon.s <- lapply( strsplit( gene.stats$ex.s, "," ), as.numeric )
intron.s <- lapply( strsplit( gene.stats$in.s, ","), as.numeric )

names(exon.s) <- gene.stats$gene
names(intron.s) <- gene.stats$gene

all.intron.s <- lapply( strsplit( all.gene.stats$in.s, ","), as.numeric )
all.exon.s <- lapply( strsplit( all.gene.stats$ex.s, ","), as.numeric )

## divide into first and other introns for each gene

all.intron.h <- hist( log2( 1 + unlist(all.intron.s) ), breaks=50)
abline(v=8, col='red')

intron.s.h <- tapply( all.intron.s, all.gene.stats$sp, function(x){
    hist( log2(1 + unlist(x)), breaks=all.intron.h$breaks, plot=FALSE ) })

first.intron.s.h <- tapply( all.intron.s, all.gene.stats$sp, function(x){
    hist( log2(1 + unlist( sapply(x, function(y){ y[1] } ))), breaks=all.intron.h$breaks,
         plot=FALSE ) })

other.intron.s.h <- tapply( all.intron.s, all.gene.stats$sp, function(x){
    hist( log2(1 + unlist( sapply(x, function(y){ y[-1] } ))), breaks=all.intron.h$breaks,
         plot=FALSE ) })

## I got this function by rerranging the equeation for the two-sample Kolmogorov-Smirnov
## test found on Wikipedia. It is unlikely to be absolutely correct (as it is not used
## by ks.test().
## However, it seems to pretty much give the same results as the ks.test for lower p-values
## though it gives a maximum p value of 2 for d = 0. (2 * exp(0) )
## For our purpose however, this should be fine.
d2p <- function(d, n1, n2){
    n1 <- as.double(n1)
    n2 <- as.double(n2)
    2 * exp( -2*(d^2) / ((n1+n2)/(n1*n2)) )
}

## we can also do a ks test for the difference in distributions for these.
intron.s.ks <- tapply( all.intron.s, all.gene.stats$sp, function(x){
    all <- log2(1 + unlist(x))
    first <- log2(1 + unlist( sapply(x, function(y){ y[1] } )))
    other <- log2(1 + unlist( sapply(x, function(y){ y[-1] } )))
    ctl <- log2(1 + sample(unlist(x), length(x)/2))
##    
    first.ks=ks.test(all, first)
    first.2.ks=ks.test(other, first)
    other.ks=ks.test(all, other)
    other.s.ks=lapply(1:10, function(x){
        ks.test(all, sample(other, size=length(first))) })
    ctl.ks=ks.test(all, ctl)
    list(first=first.ks, first.2=first.2.ks, other=other.ks, other.s=other.s.ks, ctl=ctl.ks)
})

sp.db2sp <- function(x){
    x <- sub( "([^_]+)_([^_]+)_?.*", "\\1 \\2", x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

uc.1 <- function(x){
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

names(intron.s.h) <- sp.db2sp( names(intron.s.h) )
names(first.intron.s.h) <- sp.db2sp( names(first.intron.s.h) )
names(other.intron.s.h) <- sp.db2sp( names(other.intron.s.h) )
names(intron.s.ks) <- sp.db2sp( names(intron.s.ks))

teleost.b <- read.table("../R_172_genomes/teleost_b.txt", sep="\t" )
mammal.b <- read.table("../R_172_genomes/mammal_b.txt", sep="\t" )
sauria.b <- read.table("../R_172_genomes/sauria_b.txt", sep="\t" )

sp.class <- cbind( teleost.b, mammal.b, sauria.b )
colnames(sp.class) <- c('teleost', 'mammal', 'sauria')
rownames(sp.class) <- uc.1( rownames(teleost.b) )

sp.domain <- apply( sp.class, 1, function(x){
    if(!sum(x))
        return('others')
    colnames(sp.class)[ which(x) ] })

class.col <- rgb(c(1, 0.5, 0.5, 0.5), c(0.5, 1, 0.5, 0.5), c(0.5, 0.5, 1, 0.5))
names(class.col) <- c('teleost', 'mammal', 'sauria', 'other')
sp.col <- rgb(sp.class[,1] * 0.5 + 0.5, sp.class[,2] * 0.5 + 0.5, sp.class[,3] * 0.5 + 0.5)
names(sp.col) <- rownames(sp.class)

sp.col.2 <- rgb(sp.class[,1] * 0.5 + 0.2, sp.class[,2] * 0.5 + 0.2, sp.class[,3] * 0.5 + 0.2)
names(sp.col.2) <- rownames(sp.class)

## take colours from the trees_distances table.
sp.col.3 <- readRDS( "../R_trees_distances/sp_col_3.rds" )  ## names are in lowercase _ separated
## positions from taxonomic tree
sp.y <- readRDS("../R_trees_distances/sp_y.rds")
class.col.3 <- readRDS("../R_trees_distances/class_col_3.rds")

tmp <- read.table("../R_172_genomes/genome_sizes.txt")
genome.sizes <- tmp[,1]
names(genome.sizes) <- uc.1( sp.db2sp( rownames(tmp)) )
rm(tmp)

genome.sizes <- sort(genome.sizes)

## first we will make a full figure of all of the distributions
## with a divider at 8 (256 bp)
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

pdf("species_intron_size_distributions.pdf", width=a4.w, height=a4.h,
    title='Intron size distributions')
par(omi=c(0.5, 0.5, 0.5, 0.5))
par(mfrow=c(5, 4))
par(mar=c(5.1, 2.1, 4.1, 2.1))
for(sp in names(genome.sizes)){
    plot( intron.s.h[[sp]], main=sp, xlab="log2 intron size", ylab='',
         cex.main=0.9, cex.lab=0.8, cex.axis=0.8, col=sp.col[sp])
    abline(v=8, col='black')
}
dev.off()

pdf("species_intron_size_distributions_2.pdf", width=a4.w, height=a4.h,
    title='Intron size distributions')
par(omi=c(0.5, 0.5, 0.5, 0.5))
par(mfrow=c(5, 4))
par(mar=c(2.1, 2.1, 4.1, 2.1))
for(sp in names(genome.sizes)){
    d <- cbind( intron.s.h[[sp]]$density, first.intron.s.h[[sp]]$density,
               other.intron.s.h[[sp]]$density )
    m <- intron.s.h[[sp]]$mids
    plot( m, d[,1], ylim=range(d), col='black', lwd=2, main='', xlab='',
         ylab='', type='l' )
    mtext(sp, cex=0.75)
    with(par(), text(usr[2], usr[4], sp.domain[sp],                                                          adj=c(1.05,1.2)), cex=0.4)
    lines(m, d[,2], lwd=2, col='red')
    lines(m, d[,3], lwd=2, col='blue')
    abline(v=log2(c(105, 150, 256)),lty=2, col='black')
}
dev.off()

## let us plot the d-statistics
plot.d <- function(ks=intron.s.ks, col=sp.col, gs=genome.sizes, ...){
    par(mar=c(5.1, 4.1, 4.1, 4.1))
    d <- sapply( ks, function(x){ c('f'=x$first$statistic,
                                    'o'=x$first.2$statistic)})
    colnames(d) <- sub("_", " ", names(ks))
    substr(colnames(d), 1, 1) <- toupper( substr(colnames(d), 1, 1) )
    sp.n <- names(sp.col)
    plot(1:ncol(d), log2(gs[sp.n]), type='n', xlab='species', ylab='',
         xlim=c(0.5, 0.5+ncol(d)), xaxs='i', yaxt='n')
    axis(side=4)
    mtext('log2 genome size', side=4, line=2.5)
    rect(1:ncol(d)-0.5, 0, 1:ncol(d)+0.5, log2(gs[sp.n]), col=rgb(0.9, 0.9, 0.9),
         border=NA)
    plot.window( xlim=c(0.5,0.5+ncol(d)), ylim=c(0, max(d)), xaxs='i')
    axis(side=2)
    mtext('KS D statistic', side=2, line=2.5)
##    points(1:ncol(d), d[1,sp.n], pch=19, col=sp.col)
    points(1:ncol(d), d[2,sp.n], pch=19, col=sp.col)
    legend('topright', legend=names(class.col), col=class.col, pch=19, bty='n')
    invisible(list('x'=1:ncol(d), 'y'=d[2,sp.n], 'sp'=sp.n))
}

ks.coords <- plot.d()
##plot.window(xlim=c(0.5, 0.5+length(ks.coords$x)), ylim=c(0, max(all.intron.s.n[,'5'] / all.intron.s.n[,'n'], na.rm=TRUE)), xaxs='i')
##with(ks.coords, points( x, all.intron.s.n[ sp, '5'] / all.intron.s.n[ sp, 'n'], pch=1, col='black', lwd=3))
ks.pts <- identify( ks.coords$x, ks.coords$y, ks.coords$sp, pos=TRUE )

plt <- recordPlot( )

cairo_pdf("first_intron_ks_d.pdf", width=a4.w*0.8*pdf.m, height=a4.w*0.4*pdf.m)
plot.d()
with(ks.coords, text( x[ks.pts$ind], y[ks.pts$ind], sp[ks.pts$ind], pos=ks.pts$pos, cex=0.65 ))
dev.off()



plot.d.gs <- function(ks=intron.s.ks, col=sp.col, gs=genome.sizes, ...){
    sp.n <- names(gs)
    tel.sp <- sp.n[ sp.class[,'teleost'] ]
    sau.sp <- sp.n[ sp.class[,'sauria'] ]
    mam.sp <- sp.n[ sp.class[,'mammal'] ]
    ## this doesn't really work as the x-axis is species index, not genome size
    ## we would need a different plot for that.
    d <- sapply( ks, function(x){ c('f'=x$first$statistic,
                                    'o'=x$other$statistic)})
    plot(log2(gs), d[1,sp.n], pch=19, col=sp.col[sp.n])
    abline( lm( d[1,tel.sp] ~ log2(gs[tel.sp]) ), col=sp.col[tel.sp])
    abline( lm( d[1,sau.sp] ~ log2(gs[sau.sp]) ), col=sp.col[sau.sp])
    abline( lm( d[1,mam.sp] ~ log2(gs[mam.sp]) ), col=sp.col[mam.sp])
}

plot.d.gs()
## but those relationships are actually not that nice..


## the number of introns with size below min.s,
## 
n.below <- function(x, t.s=c(5,8)){
    x <- x[!is.na(x)]
    y <- c( sapply(t.s, function(t){
        sum(x < t) }), length(x) )
    names(y) <- c(as.character(t.s), 'n')
    y
}

intron.s.n <- tapply( intron.s, gene.stats[,'sp'], function(x){
    n.below( log2(1+unlist(x)) ) })

all.intron.s.n <- tapply( all.intron.s, all.gene.stats[,'sp'], function(x){
    n.below( log2(1+unlist(x))) })

all.first.intron.s.n <- tapply( all.intron.s, all.gene.stats[,'sp'], function(x){
    n.below( log2(1 + sapply( x, function(y){ y[1] })) ) })

all.not.first.intron.s.n <- tapply( all.intron.s, all.gene.stats[,'sp'], function(x){
    n.below( log2(1 +  unlist(lapply( x, function(y){ y[-1] }))) ) })


names(intron.s.n) <- uc.1( sp.db2sp( names(intron.s.n)))
intron.s.n <- t(sapply( intron.s.n, eval ))

names(all.intron.s.n) <- uc.1( sp.db2sp( names(all.intron.s.n)))
all.intron.s.n <- t(sapply( all.intron.s.n, eval ))

names(all.first.intron.s.n) <- uc.1( sp.db2sp( names(all.first.intron.s.n)))
all.first.intron.s.n <- t(sapply( all.first.intron.s.n, eval ))

names(all.not.first.intron.s.n) <- uc.1( sp.db2sp( names(all.not.first.intron.s.n)))
all.not.first.intron.s.n <- t(sapply( all.not.first.intron.s.n, eval ))


par(mfrow=c(1,2))
sp <- names(genome.sizes)

cairo_pdf("unreasonable_intron_proportion.pdf", width=a4.w * 0.5 * pdf.m, height=a4.h * 0.3 * pdf.m)
plot( all.intron.s.n[sp,1] / all.intron.s.n[sp,'n'], col=sp.col.2[ sp ], lwd=2,
     xlab='species', ylab='intron size < 33 bp')
legend('topleft', legend=c('Teleost', 'Mammal', 'Sauria', 'other'),
       pch=1, col=sort(unique( sp.col.2 ), decreasing=TRUE), lwd=2, lty=0 )
dev.off()

sort( intron.s.n[,1] / intron.s.n[,'n'] )
plot( all.intron.s.n[,1] / all.intron.s.n[,'n'] )
## so we do have a few outliers; these should perhaps be excluced from
## additional analyses.

## We should probably include this as a specific layout.
plot.small.r <- function(int.s.r=intron.s.n, pch=1, cols=sp.col, ...){
    sp <- names(genome.sizes)
    plot( log2(genome.sizes[sp]), (int.s.r[sp,2] - int.s.r[sp,1]) /
                                  (int.s.r[sp,'n'] - int.s.r[sp,1]),
         col=cols[sp], xlab='log2 genome size', ylab='proportion < 256 bp',
         pch=pch, ...)
    legend('topright', legend=c('Teleosts', 'Mammals', 'Sauria', 'others'),
           pch=pch, col=sort(unique(cols), decreasing=TRUE))
}

plot.small.r.2 <- function(int.s.r.1=all.not.first.intron.s.n,
                           int.s.r.2=all.not.first.intron.s.n, pch1=1, pch2=2, cols=sp.col,
                           pt.cex=1, ...){
    sp <- names(genome.sizes)
    y1 <- (int.s.r.1[sp,2] - int.s.r.1[sp,1]) / (int.s.r.1[sp,'n'] - int.s.r.1[sp,1])
    y2 <- (int.s.r.2[sp,2] - int.s.r.2[sp,1]) / (int.s.r.2[sp,'n'] - int.s.r.2[sp,1])
    x <- log2(genome.sizes[sp])
    plot( x, y1, ylim=range(c(y1,y2)),
         col=cols[sp], xlab='log2 genome size', ylab='proportion < 256 bp',
         pch=pch1, cex=pt.cex,  ...)
    points( x, y2, col=cols[sp], pch=pch2, cex=pt.cex)
    segments(x, y1, x, y2) 
    legend('topright', legend=c('Teleosts', 'Mammals', 'Sauria', 'others'),
           pch=pch2, col=sort(unique(cols), decreasing=TRUE))
    legend('topright', legend=c('rank > 1','rank = 1' ), pch=c(pch1, pch2), inset=c(0,0.25),
           bty='n')
}


cairo_pdf("genome_size_short_intron_prop.pdf", width=a4.w * 0.4 * pdf.m, height=a4.h * 0.2 * pdf.m )
##plot.small.r(int.s.r=intron.s.n, pch=19)
par('mar'=c(5.1, 4.1, 2.1, 2.1))
plot.small.r(int.s.r=all.intron.s.n, pch=19, cex=1, lwd=3, cols=sp.col)
dev.off()

plot.small.r(int.s.r=all.first.intron.s.n, pch=19, cex=1, lwd=3, cols=sp.col)
plot.small.r(int.s.r=all.not.first.intron.s.n, pch=19, cex=1, lwd=3, cols=sp.col)

cairo_pdf("genome_size_short_intron_prop_2.pdf", width=a4.w * 0.4 * pdf.m, height=a4.h * 0.2 * pdf.m )
par('mar'=c(5.1, 4.1, 2.1, 2.1))
plot.small.r.2( all.not.first.intron.s.n, all.first.intron.s.n, pch1=1, pch2=19, lwd=3, cols=sp.col)
dev.off()


### How to summarise the distributions?
hist.summary.plot <- function(sp, xlim=NULL, cols=sp.col, lwd=1,
                              draw.means=FALSE, main=NULL, ...){
    int.s <- intron.s.h[sp]
    dens <- sapply( int.s, function(x){ x$density } )
    mids <- sapply( int.s, function(x){ x$mids } )
    if(is.null(xlim))
        xlim <- range(mids)
    plot(1,1, type='n', xlab='log2 intron size', ylab='density',
         xlim=xlim, ylim=range(dens), main=main )
    if(draw.means)
        lines( mids[,1], rowMeans( dens ), col='grey', lwd=lwd, ... )
    invisible(sapply(sp, function(id){
        lines( mids[,id], dens[,id], col=cols[id], lwd=lwd, ... )}))
}



sp <- names(genome.sizes)
b <- (all.intron.s.n[sp,1] / all.intron.s.n[sp,'n']) < 0.02
## curiously the number of stupidly short introns is largest in mammals;
## notably it is the smallest in well annotated mammals like mouse and human..
## 0.02, 2% is a bit much, but whatevaer
hist.summary.plot(sp[b], xlim=c(4.5, 20), draw.means=FALSE, lwd=1.5, cols=sp.col.2)
legend('topright', legend=c('teleost', 'mammal', 'sauria', 'others'),
       lwd=1.5, col=sort(unique(sp.col.2), decreasing=TRUE) )

## sp.col.3 <- hsvScale( genome.sizes )
## names(sp.col.3) <- names(genome.sizes)

tel.col <- hsvScale( log2(genome.sizes[ sp.class[,'teleost'] ]))
names(tel.col) <- names(genome.sizes[ sp.class[,'teleost'] ])

mammal.col <- hsvScale( log2(genome.sizes[ sp.class[,'mammal'] ]))
names(mammal.col) <- names(genome.sizes[ sp.class[,'mammal'] ])

sauria.col <- hsvScale( log2(genome.sizes[ sp.class[,'sauria'] ]))
names(sauria.col) <- names(genome.sizes[ sp.class[,'sauria'] ])

cairo_pdf("intron_size_distribution_summary.pdf", width=a4.w * 0.6 * pdf.m, height=a4.h * 0.6 * pdf.m )
##par(mfrow=c(2,2))
layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE))
hist.summary.plot(sp[b & sp.class[,'teleost']], xlim=c(4.5, 20), draw.means=FALSE, cols=tel.col,
                  lwd=2, main='Teleosts')
abline(v=log2(c(105, 150, 256)))
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))
hist.summary.plot(sp[b & sp.class[,'mammal']], xlim=c(4.5, 20), draw.means=FALSE, cols=mammal.col,
                  lwd=2, main='Mammals')
abline(v=log2(c(105, 150, 256)))
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
hist.summary.plot(sp[b & sp.class[,'sauria']], xlim=c(4.5, 20), draw.means=FALSE, cols=sauria.col,
                  main='Sauria')
abline(v=log2(c(105, 150, 256)))
with(par(), mtext('C', at=usr[1], cex=mt.cex, line=1))
plot.small.r(int.s=all.intron.s.n, pch=19, cex=1, lwd=3, cols=sp.col )
with(par(), mtext('D', at=usr[1], cex=mt.cex, line=1))
plot.small.r.2( all.not.first.intron.s.n, all.first.intron.s.n, pch1=1, pch2=19, lwd=1, cols=sp.col,
               pt.cex=1.5)
with(par(), mtext('E', at=usr[1], cex=mt.cex, line=1))
dev.off()

#


tel.sel <- unlist(lapply(c('takifugu', 'betta', 'gambus',
                           ##'clupea',
                           'esox', 'wildenow', 'danio'),
                         function(x){ grep(x, sp, value=TRUE,
                                           ignore.case=TRUE) }))


cols <- hsvScale( 1:length(tel.sel))
names(cols) <- tel.sel
hist.summary.plot(tel.sel, xlim=c(4.5, 20), cols=cols, lwd=2)                                   


## For the supplementary I would like to have a figure that shows exon / intron ratios
exon.intron.r <- tapply( 1:nrow(gene.stats), gene.stats$sp, function(i){
    c(sum(unlist(exon.s[i]), na.rm=T), sum(unlist(intron.s[i]), na.rm=T)) })
exon.intron.r <- t(sapply(exon.intron.r, eval))
colnames(exon.intron.r) <- c('exon', 'intron')

sp.n <- rownames(exon.intron.r)

cairo_pdf("exon_intron_ratios.pdf", width=a4.w * 0.9 * pdf.m, height=a4.w * 0.75 * pdf.m )
par(mfrow=c(2,2))
plot( genome.sizes[uc.1(sp.db2sp(sp.n))], exon.intron.r[,'exon'],
     col=sp.col.3[sp.n], pch=19, xlab='Genome size', ylab='Sum exon lengths')
with(par(), mtext("A", at=usr[1], line=1, cex=mt.cex))
##
plot( genome.sizes[uc.1(sp.db2sp(sp.n))], exon.intron.r[,'intron'],
     col=sp.col.3[sp.n], pch=19, xlab='Genome size', ylab='Sum intron lengths')
with(par(), mtext("B", at=usr[1], line=1, cex=mt.cex))
##
plot( sp.y[sp.n], exon.intron.r[,'exon'] / (exon.intron.r[,'exon'] + exon.intron.r[,'intron']),
     col=sp.col.3[sp.n], pch=19, xlab='Species', ylab='Sum exon length / Sum transcript length')
with(par(), mtext("C", at=usr[1], line=1, cex=mt.cex))
##
plot( genome.sizes[uc.1(sp.db2sp(sp.n))], exon.intron.r[,'exon'] / (exon.intron.r[,'exon'] + exon.intron.r[,'intron']),
     col=sp.col.3[sp.n], pch=19, xlab='Genome size', ylab='Sum exon length / Sum transcript length')
with(par(), mtext("D", at=usr[1], line=1, cex=mt.cex))
legend('topright', uc.1(names(class.col.3)), pch=19, col=class.col.3)
dev.off()


