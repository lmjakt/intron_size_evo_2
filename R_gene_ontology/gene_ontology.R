## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

source("~/R/general_functions.R")


read.data <- function(f){
    read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
}

intron.orth.id <- read.data("../R_172_genomes/dr_intron_orthology_id.txt")
intron.orth.l <- read.data("../R_172_genomes/dr_intron_orthology_l.txt")

## we also want the telost variance, though it seems that we may have made
## some mistakes with that before.

tel.b <- read.data("../R_172_genomes/teleost_b.txt")
mam.b <- read.data("../R_172_genomes/mammal_b.txt")
sau.b <- read.data("../R_172_genomes/sauria_b.txt")

sp.class <- cbind('teleost'=tel.b[,1], 'mammal'=mam.b[,1], 'sauria'=sau.b[,1])
rownames(sp.class) <- sub(" ", ".", rownames(tel.b), fixed=TRUE )

## the rownames in the different tables are now the same so we can do
intron.orth.id <- intron.orth.id[ , rownames(sp.class)]
intron.orth.l <- intron.orth.l[, rownames(sp.class)]

## for convenience we also make this one
int.s.l <- log2(intron.orth.l)

## we will recalculate the variance measuremens as there seemed to
## be some issue with the saved ones.

var.par <- function(x){
    c('n'=sum(!is.na(x)), 'mean'=mean(x, na.rm=TRUE), 'sd'=sd(x, na.rm=TRUE),
      'min'=min(x, na.rm=TRUE), 'max'=max(x, na.rm=TRUE), quantile(x, probs=seq(0,1,0.1), na.rm=TRUE))
}

## I need also the int.s.b...
int.s.b <- readRDS("../R_trees_distances/int_s_b.rds")
names(int.s.b) <- sub("_", ".", names(int.s.b), fixed=TRUE)
all(names(int.s.b) == rownames(sp.class))  ## TRUE

## I only need tel.var, but might as well read the others in for now in case I
## wish to use them.
tel.var <- t(apply( int.s.l[ ,sp.class[,'teleost'] & int.s.b], 1, var.par ))
mam.var <- t(apply( int.s.l[ ,sp.class[,'mammal'] & int.s.b], 1, var.par ))
sau.var <- t(apply( int.s.l[ ,sp.class[,'sauria'] & int.s.b], 1, var.par ))
                             
### packages for gene ontology and human annotation
require("GOstats")
require("org.Hs.eg.db")

ens2eg <- as.list(org.Hs.egENSEMBL2EG)
table( sapply( ens2eg, length ))

##     1     2     3     4     5     7 
## 28813   364    22     6     5     1 

eg2ens <- as.list(org.Hs.egENSEMBL)

## I have human ids in
## intron.orth.t[[1]]$id[,'homo sapiens']
## but note that not all of the rows are present

hum.g.na <- is.na(intron.orth.id[,'homo.sapiens'])
sum(hum.g.na)
## [1] 303

length(unique( intron.orth.id[!hum.g.na,'homo.sapiens']))
## 5708

hum.g <- unique( intron.orth.id[!hum.g.na,'homo.sapiens'])
sum( hum.g %in% names( ens2eg ) )
## 5695  ## this will be our universe..

hum.universe <- unlist(ens2eg[ hum.g ])

do.hp.test <- function(universe, sel.ens, ont='BP'){
    sel.eg <- unique(unlist(ens2eg[ sel.ens ]))
    params <- new("GOHyperGParams",
                  geneIds=sel.eg,
                  universeGeneIds=universe,
                  annotation="org.Hs.eg.db",
                  ontology=ont,
                  conditional=FALSE,
                  testDirection="over")
    hyperGTest(params)
}


tel.min.bp <- lapply(5:15, function(x){
    sel.ens <- unique( intron.orth.id[ tel.var[,'min'] > x, 'homo.sapiens'] )
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens,
         'bp'=do.hp.test( hum.universe, sel.ens, ont='BP'),
         'mf'=do.hp.test( hum.universe, sel.ens, ont='MF'),
         'cc'=do.hp.test( hum.universe, sel.ens, ont='CC'))
})
                         
for(i in 1:length(tel.min.bp)){
    print(head(summary(tel.min.bp[[i]]$bp)))
    print(head(summary(tel.min.bp[[i]]$mf)))
    print(head(summary(tel.min.bp[[i]]$cc)))
    inpt <- readline(paste((5:15)[i], "next: "))
}

### we have some very nice set of data of enrichments there..
### First we want to identify the set of terms that were strongly enriched

tel.min.sum.bp <- do.call(rbind, lapply(1:length(tel.min.bp), function(i){
    tmp <- summary(tel.min.bp[[i]]$bp)
    if(!nrow(tmp))
        return(NULL)
    cbind('length'=i+4, tmp)}))

tel.min.sum.mf <- do.call(rbind, lapply(1:length(tel.min.bp), function(i){
    tmp <- summary(tel.min.bp[[i]]$mf)
    if(!nrow(tmp))
        return(NULL)
    cbind('length'=i+4, tmp)}))

tel.min.sum.cc <- do.call(rbind, lapply(1:length(tel.min.bp), function(i){
    tmp <- summary(tel.min.bp[[i]]$cc)
    if(!nrow(tmp))
        return(NULL)
    cbind('length'=i+4, tmp)}))


### I do not quite understand this, but we hvae a very definite enrichment
### for the cell periperhy. I worry that this may be related to such genes
### in general having larger number of introns. To control for this we
### can do random sampling of introns.
### ideally I should do this more than once, but.. it's quite time consuming.

tel.min.ctl <- lapply(5:15, function(x){
    n <- sum( tel.var[,'min'] > x )
    i <- sample(1:nrow(tel.var), size=n)
    sel.ens <- unique( intron.orth.id[ i, 'homo.sapiens'] )
    sel.ens <- sel.ens[ !is.na(sel.ens) ]
    list('sel'=sel.ens,
         'bp'=do.hp.test( hum.universe, sel.ens, ont='BP'),
         'mf'=do.hp.test( hum.universe, sel.ens, ont='MF'),
         'cc'=do.hp.test( hum.universe, sel.ens, ont='CC'))
})

tel.ctl.sum.bp <- do.call(rbind, lapply(1:length(tel.min.ctl), function(i){
    tmp <- summary(tel.min.ctl[[i]]$bp)
    if(!nrow(tmp))
        return(NULL)
    cbind('length'=i+4, tmp)}))

tel.ctl.sum.mf <- do.call(rbind, lapply(1:length(tel.min.ctl), function(i){
    tmp <- summary(tel.min.ctl[[i]]$mf)
    if(!nrow(tmp))
        return(NULL)
    cbind('length'=i+4, tmp)}))

tel.ctl.sum.cc <- do.call(rbind, lapply(1:length(tel.min.ctl), function(i){
    tmp <- summary(tel.min.ctl[[i]]$cc)
    if(!nrow(tmp))
        return(NULL)
    cbind('length'=i+4, tmp)}))

### For MF, we actually get lower p-values for the random selection of introns.
### This indicates that we such genes are likely to have many introns. The selection
### is not random.
### Nevertheless we need to see what the ctl values for terms identifed by long introns
### are.
### We should probbaly do more than one set of random selections to confirm this.

## An extended set of random analyses, but for a smaller range of sizes;
tel.min.ctl.l <- lapply(1:50, function(z){
    lapply(5:15, function(x){
        n <- sum( tel.var[,'min'] > x )
        i <- sample(1:nrow(tel.var), size=n)
        sel.ens <- unique( intron.orth.id[ i, 'homo.sapiens'] )
        sel.ens <- sel.ens[ !is.na(sel.ens) ]
        list('sel'=sel.ens,
             'bp'=do.hp.test( hum.universe, sel.ens, ont='BP'),
             'mf'=do.hp.test( hum.universe, sel.ens, ont='MF'),
             'cc'=do.hp.test( hum.universe, sel.ens, ont='CC'))
    })
})

tel.ctl.l.sum.bp <- lapply(tel.min.ctl.l, function(x){
    do.call(rbind, lapply(1:length(x), function(i){
        tmp <- summary(x[[i]]$bp)
        if(!nrow(tmp))
            return(NULL)
        cbind('length'=i+4, tmp)
    }))
})

tel.ctl.l.sum.mf <- lapply(tel.min.ctl.l, function(x){
    do.call(rbind, lapply(1:length(x), function(i){
        tmp <- summary(x[[i]]$mf)
        if(!nrow(tmp))
            return(NULL)
        cbind('length'=i+4, tmp)
    }))
})

tel.ctl.l.sum.cc <- lapply(tel.min.ctl.l, function(x){
    do.call(rbind, lapply(1:length(x), function(i){
        tmp <- summary(x[[i]]$cc)
        if(!nrow(tmp))
            return(NULL)
        cbind('length'=i+4, tmp)
    }))
})

collapse.ctl <- function(ctls){
    lapply(1:length(ctls), function(h){
        x <- ctls[[h]]
        tmp <- tapply( 1:nrow(x), x[,2], function(i){
            cbind('i'=h, x[ i[ which.min(x[i,'Pvalue']) ], ]) })
        tmp <- do.call(rbind, tmp)
        tmp <- tmp[ order(tmp$Pvalue), ]
        tmp
    })
}

collapse.collapsed <- function(col){
    col <- do.call(rbind, col)
    tmp <- do.call( rbind, tapply( 1:nrow(col), col[,3], function(i){
        col[ i[ which.min(col[i,'Pvalue'])], ] }))
    tmp[ order(tmp$Pvalue), ]
}

tel.ctl.bp <- collapse.ctl( tel.ctl.l.sum.bp )
tel.ctl.bp.sum <- collapse.collapsed( tel.ctl.bp )

tel.ctl.mf <- collapse.ctl( tel.ctl.l.sum.mf )
tel.ctl.mf.sum <- collapse.collapsed( tel.ctl.mf )

tel.ctl.cc <- collapse.ctl( tel.ctl.l.sum.cc )
tel.ctl.cc.sum <- collapse.collapsed( tel.ctl.cc )

### there are some interesting observations there that most likely relate
### to intron number; that could be tested separately, but there are better
### ways of doing that.
### and that could be done for individual species much more efficiently.
### I can come back to this later to check if I need to.

head( tel.min.sum.bp[ order(tel.min.sum.bp$Pvalue), ] )
##      length     GOBPID       Pvalue OddsRatio ExpCount Count Size
## 1299      9 GO:0048731 2.253415e-23  1.882262 496.5474   651 1457
## 1300      9 GO:0032501 6.436190e-23  1.786480 710.2297   876 2084
## 1301      9 GO:0048856 1.943490e-22  1.804094 602.1957   761 1767
## 1302      9 GO:0007275 5.119830e-22  1.812506 555.8468   710 1631
## 1303      9 GO:0023052 5.492281e-22  1.802916 575.6133   731 1689
## 601       8 GO:0007154 6.236639e-22  1.770049 897.0906  1059 1707
##                                    Term
## 1299                 system development
## 1300   multicellular organismal process
## 1301   anatomical structure development
## 1302 multicellular organism development
## 1303                          signaling
## 601                  cell communication

tel.ctl.sum.bp[ tel.ctl.sum.bp$GOBPID == 'GO:0048731', ]
## <0 rows> (or 0-length row.names)
tel.ctl.sum.bp[ tel.ctl.sum.bp$GOBPID == 'GO:0032501', ]
##     length     GOBPID       Pvalue OddsRatio ExpCount Count Size
## 475     10 GO:0032501 0.0033668154  1.182759 624.3691   669 2084
## 542     11 GO:0032501 0.0004741639  1.286438 341.8599   386 2084
##                                 Term
## 475 multicellular organismal process
## 542 multicellular organismal process

head( tel.min.sum.mf[ order(tel.min.sum.mf$Pvalue), ] )
##     length     GOMFID       Pvalue OddsRatio  ExpCount Count Size
## 183      9 GO:0038023 1.467437e-11  2.736011  63.74386   108  188
## 184      9 GO:0060089 3.165988e-11  2.553616  71.20325   117  210
## 79       8 GO:0060089 6.377209e-11  2.660838 109.87911   155  210
## 278     10 GO:0038023 7.510781e-11  2.826436  37.14545    75  188
## 279     10 GO:0004888 1.724910e-10  3.238056  26.87117    59  136
## 80       8 GO:0038023 1.936587e-10  2.745881  98.36796   140  188
##                                          Term
## 183               signaling receptor activity
## 184             molecular transducer activity
## 79              molecular transducer activity
## 278               signaling receptor activity
## 279 transmembrane signaling receptor activity
## 80                signaling receptor activity

tel.ctl.sum.mf[ tel.ctl.sum.mf$GOMFID == 'GO:0038023', ]
## no rows
tel.ctl.sum.mf[ tel.ctl.sum.mf$GOMFID == 'GO:0060089', ]
## no rows

### we can control for accidental enrichment of the specific groups
### later in the proces when we consider these.
## Let us reduce these tables by removing duplicate entries

reduce.sum.table <- function(x){
    ## second column is the GOxxID
    x <- tapply( 1:nrow(x), x[,2], function(i){
        x[ i[which.min(x[i,'Pvalue'])], ] })
    x <- do.call(rbind, x)
    x[ order(x$Pvalue), ]
}

bp.sum <- reduce.sum.table( tel.min.sum.bp )
mf.sum <- reduce.sum.table( tel.min.sum.mf )
cc.sum <- reduce.sum.table( tel.min.sum.cc )

### Look at the patterns of enrichment for specific groups
## to get all genes annotated with a specific group:
go2eg <- as.list(org.Hs.egGO2ALLEGS)

eg.from.enrichment.table <- function(x, max.p){
    x <- x[ x$Pvalue <= max.p, 2]
    eg <- lapply(x, function(y){ unlist(go2eg[ y ]) })
    names(eg) <- x
    eg
}

## I think the permutation test used here is not valid, as it does not
## take into consideration the effect of intron number in the selection
## But I think that I have to consider that in a separate function
## where we can look at the distributions of the number of introns
## of the set
gene.set.association <- function(eg, intron.stat=tel.var[,'min'],
                                 ids=intron.orth.id,
                                 min.l=5, int.criteria=which.max,
                                 perm.n=10, decreasing=FALSE,
                                 entrez=TRUE){
    if(entrez){
        ens.id <- unlist( eg2ens[eg] )
        ens.id <- ens.id[!is.na(ens.id)]
    }else{
        ens.id <- eg
    }
    ##
    ## obtain the indices of the longest intron for each gene (ens.i)
    ens.i <- tapply( 1:length(intron.stat), ids[,'homo.sapiens'], function(i){
        i[ int.criteria( intron.stat[i] )] })
    ## then get the length of the longest intron
    int.s <- intron.stat[ens.i]
    ## get the human ensembl ids (int stands for intron here)
    int.id <- ids[ens.i, 'homo.sapiens']
    ## very short introns are probably unreasonable; ignore these
    ## removing both the ids and the intron sizes
    int.id <- int.id[int.s >= min.l]
    int.s <- int.s[int.s >= min.l] 
    ## reorder by intron size from small to big
    o <- order(int.s, decreasing=decreasing)
    int.id <- int.id[o]
    int.s <- int.s[o]
    
    ## which of the introns are in eg?
    o.b <- int.id %in% ens.id
    ## the following will give the number of genes in the set (eg)
    ## with an increase in size
    get.stats <- function(b){
        q <- cumsum( b )
        p <- phyper(q, sum(b), sum(!b), 1:length(b))
        p2 <- phyper(q-1, sum(b), sum(!b), 1:length(b), lower.tail=FALSE)
        exp <- 1:length(b) * sum(b) / length(b)
        data.frame('id'=int.id, 'l'=int.s, 'b'=b, 'q'=q, 'p'=p, 'p2'=p2, 'exp'=exp )
    }
    stats <- get.stats(o.b)
    stats.permed <- lapply( 1:perm.n, function(i){
        get.stats( sample(o.b) )
    })
    list('o'=stats, 'o.perm'=stats.permed, 'ens.i'=ens.i, 'ens.id'=ens.id)
}

## takes an object as returned by gene.set.association
plot.association.stats <- function(stats.l, min.sample=50,
                                   cols=c(1, 2, 4, 6), p2=FALSE,
                                   draw.permuted=FALSE,
                                   draw.min.p=FALSE,
                                   legpos='bottomright', leg.inset=c(0.0, 0.1),
                                   ...){
    ## only plot points where a gene is a member;
    ## makes the plot a little bit cleaner
    stats <- stats.l$o
    stats.p <- stats.l$o.perm
    b <- stats$b
    b[1:min.sample] <- FALSE
    ## plot observed vs expected
    par(mar=c(5.1, 4.1, 4.1, 8.1))
    plot( stats$l[b], log2(stats$q[b] / stats$exp[b]),
         xlab='log2 intron size', ylab='', ##  'observed / expected',
         col=cols[1], type='l', axes=FALSE, ...)
    axis(1)
    axis(2, col=cols[1], col.axis=cols[1])
    ## lets do the permuted values
    if(draw.permuted){
        invisible( lapply(stats.p, function(x){
            lines(x$l[ x$b ], log2(x$q[ x$b ] / x$exp[ x$b ]), col=rgb(0.8, 0.8, 0.8))
        }))
    }
    ## then plot p-values
    if(!p2)
        pp <- -log10( stats$p[b] )
    else
        pp <- -log10( stats$p2[b] )
##            
    plot.window( xlim=range( stats$l[b] ), ylim=range(pp) )
    lines( stats$l[b], pp, col=cols[2] )
    axis(4, line=0, col=cols[2], col.axis=cols[2], )
    if(draw.min.p){
        min.p.i <- which.max( pp )
        abline(v=(stats$l[b])[min.p.i], col=cols[2], lty=2)
    }
    ## then plot n (or rather q)
    plot.window( xlim=range( stats$l[b] ), ylim=range(stats$q[b]) )
    lines( stats$l[b], stats$q[b], type='l', col=cols[3] )
    axis(4, line=2.5, col=cols[3], col.axis=cols[3])
    plot.window( xlim=range( stats$l[b] ), ylim=c(0, nrow(stats) ))
    lines( stats$l, 1:nrow(stats), col=cols[4], lty=2, lwd=2 )
    axis(4, line=5, col=cols[4], col.axis=cols[4])
    plot.window(xlim=c(0,1), ylim=c(0,1))
    legend(legpos, legend=c('obs/exp', '-log10 p', 'q', 'k'),
           text.col=cols, box.lty=0, lty=c(1,1,1,2), col=cols, lwd=1, inset=leg.inset,
           bg=NA)
}


bp.sum.eg <- eg.from.enrichment.table(bp.sum, 1e-8 ) ## 95 GO categories..
bp.sum.stats <- lapply( bp.sum.eg, gene.set.association ) ## 95 categories!
bp.sum.stats.e <- lapply( bp.sum.eg, gene.set.association, decreasing=TRUE ) ## 95 categories!

mf.sum.eg <- eg.from.enrichment.table(mf.sum, 1e-8)  ## only 8
mf.sum.stats <- lapply( mf.sum.eg, gene.set.association )
mf.sum.stats.e <- lapply( mf.sum.eg, gene.set.association, decreasing=TRUE )

cc.sum.eg <- eg.from.enrichment.table(cc.sum, 1e-8)  ## 31
cc.sum.stats <- lapply( cc.sum.eg, gene.set.association )
cc.sum.stats.e <- lapply( cc.sum.eg, gene.set.association, decreasing=TRUE )

tel.min.h <- hist( tel.var[,'min'], breaks=30 )
tel.int.l.h <- hist( unlist(tapply( intron.orth.l[,1], intron.orth.id[,'homo.sapiens'], length)) )

##par(mfrow=c(1,1))

cairo_pdf("GO_bp_stats_plots.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
par(mfrow=c(5,3))
par(omi=c(0.5, 0.5, 0.5, 0.5))
for(i in 1:length(bp.sum.stats)){
    plot.association.stats( bp.sum.stats[[i]], draw.min.p=TRUE, draw.permuted=FALSE,
                           main=strwrap(paste(bp.sum[i, c('GOBPID', 'Term')], collapse=" "),
                                        width=45),
                           min.sample=100 )
    b <- intron.orth.id[,'homo.sapiens'] %in%  bp.sum.stats[[i]]$ens.id
    ##inpt <- readline('next: ')
}
dev.off()

cairo_pdf("GO_bp_enr_stats_plots.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
par(mfrow=c(5,3))
par(omi=c(0.5, 0.5, 0.5, 0.5))
for(i in 1:length(bp.sum.stats)){
    plot.association.stats( bp.sum.stats.e[[i]], draw.min.p=TRUE, draw.permuted=FALSE,
                           main=strwrap(paste(bp.sum[i, c('GOBPID', 'Term')], collapse=" "),
                                        width=45),
                           min.sample=100, p2=TRUE, legpos='bottomleft', leg.inset=c(0.2,0) )
    ##inpt <- readline('next: ')
}
dev.off()


cairo_pdf("GO_mf_stats_plots.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
par(mfrow=c(5,3))
par(omi=c(0.5, 0.5, 0.5, 0.5))
for(i in 1:length(mf.sum.stats)){
    plot.association.stats( mf.sum.stats[[i]], draw.min.p=TRUE, draw.permuted=FALSE,
                           main=strwrap(paste(mf.sum[i, c('GOMFID', 'Term')], collapse=" "),
                                        width=45),
                           min.sample=100 )
}
dev.off()

cairo_pdf("GO_mf_enr_stats_plots.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
par(mfrow=c(5,3))
par(omi=c(0.5, 0.5, 0.5, 0.5))
for(i in 1:length(mf.sum.stats)){
    plot.association.stats( mf.sum.stats.e[[i]], draw.min.p=TRUE, draw.permuted=FALSE,
                           main=strwrap(paste(mf.sum[i, c('GOMFID', 'Term')], collapse=" "),
                                        width=45),
                           min.sample=100, p2=TRUE, legpos='bottomleft', leg.inset=c(0.2,0) )
}
dev.off()


cairo_pdf("GO_cc_stats_plots.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
par(mfrow=c(5,3))
par(omi=c(0.5, 0.5, 0.5, 0.5))
for(i in 1:length(cc.sum.stats)){
    plot.association.stats( cc.sum.stats[[i]], draw.min.p=TRUE, draw.permuted=FALSE,
                           main=strwrap(paste(cc.sum[i, c('GOCCID', 'Term')], collapse=" "),
                                   width=45),
                           min.sample=100 )
}
dev.off()

cairo_pdf("GO_cc_enr_stats_plots.pdf", width=a4.w * pdf.m, height=a4.h * pdf.m, onefile=TRUE )
par(mfrow=c(5,3))
par(omi=c(0.5, 0.5, 0.5, 0.5))
for(i in 1:length(cc.sum.stats)){
    plot.association.stats( cc.sum.stats.e[[i]], draw.min.p=TRUE, draw.permuted=FALSE,
                           main=strwrap(paste(cc.sum[i, c('GOCCID', 'Term')], collapse=" "),
                                   width=45),
                           min.sample=100, p2=TRUE, legpos='bottomleft', leg.inset=c(0.2,0) )
}
dev.off()

## And then, how to summarise the plots?
## The easiest is probably to simply plot the p-values from above, allowing
## one to do several lines in one go.

plot.stat.summary <- function(stats, par, n, cols=hsvScale(1:n, val=0.5),
                              stat.labels=NULL, ...){
    values <- lapply( stats[1:n], function(x){ x$o[,c('l',par)] })
    xlim <- range(unlist(lapply(values, function(x){ x[,'l'] })))
    ylim <- range(-log10(unlist(lapply(values, function(x){ x[,par] }))))
    plot(1, type='n', xlim=xlim, ylim=ylim, xlab='log2 intron length', ...)
    grid()
    for(i in 1:length(values)){
        lines( values[[i]][,'l'], -log10(values[[i]][,par]), col=cols[i] )
    }
    if(!is.null(stat.labels))
        legend('topright', stat.labels[1:n], text.col=cols)
}

plot.labels <- function(labels, n, cols=hsvScale(1:n, val=0.5), cex=1,
                        adj=c(1,1)){
    old.mar <- par('mar')
    mar <- old.mar
    mar[c(1,3)] <- c(1.1, 0)
    par('mar'=mar)
    plot.new()
    plot.window(xlim=0:1, ylim=0:1, yaxs='i')
    h <- max(strheight(labels, cex=cex)) * 1.5
    text(1, 1 - 1:n * h, labels[1:n], col=cols, cex=cex,
              adj=adj)
    par(mar=old.mar)
}

lab.cex <- 0.75
cairo_pdf("GO_pval_summary.pdf", width=a4.w * pdf.m * 0.8, height=a4.w * pdf.m * 0.4 )
layout(matrix(1:6, nrow=2), heights=c(1, 0.45))
plot.stat.summary( bp.sum.stats, 'p', 10, ylab='-log10 p')
with(par(), mtext("A", line=1, cex=mt.cex, at=usr[1]))
plot.labels(bp.sum$Term, 10, adj=c(1, 0), cex=lab.cex)
plot.stat.summary( mf.sum.stats, 'p', 8, ylab='-log10 p')
with(par(), mtext("B", line=1, cex=mt.cex, at=usr[1]))
plot.labels(mf.sum$Term, 8, adj=c(1, 0), cex=lab.cex)
plot.stat.summary( cc.sum.stats, 'p', 10, ylab='-log10 p' )
with(par(), mtext("C", line=1, cex=mt.cex, at=usr[1]))
plot.labels(cc.sum$Term, 10, adj=c(1, 0), cex=lab.cex)
dev.off()


lab.cex <- 0.75
cairo_pdf("GO_enr_pval_summary.pdf", width=a4.w * pdf.m * 0.8, height=a4.w * pdf.m * 0.4 )
layout(matrix(1:6, nrow=2), heights=c(1, 0.45))
plot.stat.summary( bp.sum.stats.e, 'p2', 10, ylab='-log10 p')
with(par(), mtext("A", line=1, cex=mt.cex, at=usr[1]))
plot.labels(bp.sum$Term, 10, adj=c(1, 0), cex=lab.cex)
plot.stat.summary( mf.sum.stats.e, 'p2', 8, ylab='-log10 p')
with(par(), mtext("B", line=1, cex=mt.cex, at=usr[1]))
plot.labels(mf.sum$Term, 8, adj=c(1, 0), cex=lab.cex)
plot.stat.summary( cc.sum.stats.e, 'p2', 10, ylab='-log10 p' )
with(par(), mtext("C", line=1, cex=mt.cex, at=usr[1]))
plot.labels(cc.sum$Term, 10, adj=c(1, 0), cex=lab.cex)
dev.off()


### Do the enriched / depleted groups contain genes with an unusually
### large number of introns? And is this what causes them to be
### selected?

## eg is a list of entrez identifiers
eg.set.intron.n.dist <- function(eg=bp.sum.eg, sum.table=bp.sum, sp='homo.sapiens',
                                 cex=1){
    int.n <- tapply(1:nrow(intron.orth.id), intron.orth.id[,sp], length)
    n.max <- max(int.n)
    all.counts <- table(c(1:n.max, int.n)) - 1
    for(i in 1:length(eg)){
        x <- eg[[i]]
        term <- sum.table[i,'Term']
        id <- sum.table[i,2]
        ens.id <- unlist( eg2ens[x] )
        ens.id <- ens.id[!is.na(ens.id)]
        b <- names(int.n) %in% ens.id
        tab.1 <- table(c(1:n.max, int.n[b])) - 1
        tab.2 <- table(c(1:n.max, int.n[!b])) - 1
        plot(1:length(all.counts), as.numeric(all.counts) / sum(all.counts),
             type='l', col='grey', cex=cex, ylim=c(0,
                                                   max(c(all.counts/sum(all.counts),
                                                         tab.1/sum(tab.1),
                                                         tab.2/sum(tab.2)))))
        lines(1:length(all.counts), tab.1/sum(tab.1), type='l', col='red', cex=cex)
        lines(1:length(all.counts), tab.2/sum(tab.2), type='l', col='blue', cex=cex)
    }
}

## I might need this as a supplementary figure
par(mfrow=c(2,2))
eg.set.intron.n.dist( bp.sum.eg[1:4], bp.sum, cex=0.5 )
## there is some relationship, but it is not strong.

## we can maybe more simply ask if there is a relationship between
## max intron length and intron number;

max.intron.l <- tapply( tel.var[,'min'], intron.orth.id[,'danio.rerio'], max )
intron.n <- tapply( 1:nrow(intron.orth.id), intron.orth.id[,'danio.rerio'], length )

par(mfrow=c(1,1))
plot(intron.n, max.intron.l )
abline(h=8, col='red')
## we could then use entropy

require(entropy)

## this needs to be plotted with some sort of heatmap..
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

int.n.l <- hist.2d( intron.n, max.intron.l, numBins=20 )

par(mfrow=c(1,3))
par('mar'=c(5.1, 4.1, 4.1, 2.1))
plot(intron.n, max.intron.l )
image(int.n.l$rb, int.n.l$cb, scale(int.n.l$h))
image(int.n.l$rb, int.n.l$cb, t(scale(t(int.n.l$h))))

image(int.n.l$rb, int.n.l$cb, int.n.l$h)

## There is a weak relationship between maximum intron length and number
## of introns. But this is not likely to have a strong effect on the
## results here. We can at this point simply run a large number of iterations of
## a control sampling like we did before.

