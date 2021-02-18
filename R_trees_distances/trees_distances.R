## read in trees created from the intron orthology.

source("../R/functions.R")
source("~/R/general_functions.R")


## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

half.w <- a4.w * 85 / 210
full.w <- a4.w * 170 / 210

## genome sizes
tmp <- read.table("../R_172_genomes/genome_sizes.txt")
genome.sizes <- tmp[,1]
names(genome.sizes) <-  rownames(tmp)
rm(tmp)

genome.sizes <- sort(genome.sizes)

assemblies <- readRDS( "../R_172_genomes/assemblies.rds" )
## the assemblies are not that useful for us to get an idea of N50s, I
## assembled a set of chromosome / scaffold sequence lengths
seq.region <- readRDS("../R_172_genomes/seq_region.rds" )
## and we can reorder these:
seq.region <- seq.region[ names(genome.sizes) ]

genome.n50 <- sapply( names(genome.sizes), function(nm){
    N50( seq.region[[nm]], genome.sizes[nm], 0.5) })


rm.diagonal <- function(m){
    for(i in 1:nrow(m))
        m[i,i] <- NA
    m
}

## distances by mutual information
int.sp.mi.dr.m <- readRDS("../R_172_genomes/int_sp_mi_dr_m.rds")
## the colnames of the genetic distances are separated by '_' rather than
## by ' '
colnames(int.sp.mi.dr.m) <- sub(' ', '_', colnames(int.sp.mi.dr.m) )
rownames(int.sp.mi.dr.m) <- sub(' ', '_', rownames(int.sp.mi.dr.m) )

int.sp.mi.dr.m <- rm.diagonal( int.sp.mi.dr.m )

## distances by genetic distances (by alignment of exons)
## (Jukes-Cantor, Jukes-Cantor-gap, Kimura two factor)
ex.align.2.jc <-  readRDS("../R_172_genomes/ex_align_2_jc.rds")
ex.align.2.jcg <-  readRDS("../R_172_genomes/ex_align_2_jcg.rds")
ex.align.2.k2 <-  readRDS("../R_172_genomes/ex_align_2_k2.rds")

all(rownames(int.sp.mi.dr.m) == rownames(ex.align.2.jc)) ## TRUE
all(colnames(int.sp.mi.dr.m) == colnames(ex.align.2.jc)) ## TRUE

all(rownames(int.sp.mi.dr.m) == rownames(ex.align.2.k2)) ## TRUE
all(colnames(int.sp.mi.dr.m) == colnames(ex.align.2.k2)) ## TRUE


## We also need to get information about lineages and quality.
## From the basic_stats analysis it looks like we can use the proportion
## of introns that are too small as a proxy for annotation quality.

## we can use stats for all genes.
### for all genes:
all.gene.stats <- read.table("../R_basic_stats/all_genes_exon_intron_stats.csv",
                             sep="\t", stringsAsFactors=FALSE)
colnames(all.gene.stats) <- c('sp', 'db', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

all.intron.s <- lapply( strsplit( all.gene.stats$in.s, ","), as.numeric )
all.exon.s <- lapply( strsplit( all.gene.stats$ex.s, ","), as.numeric )

## a couple of convenience functions
sp.db2sp <- function(x){
    x <- sub( "([^_]+)_([^_]+)_?.*", "\\1 \\2", x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

uc.1 <- function(x){
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

teleost.b <- read.table("../R_172_genomes/teleost_b.txt", sep="\t" )
mammal.b <- read.table("../R_172_genomes/mammal_b.txt", sep="\t" )
sauria.b <- read.table("../R_172_genomes/sauria_b.txt", sep="\t" )

sp.class <- cbind( teleost.b, mammal.b, sauria.b )
colnames(sp.class) <- c('teleost', 'mammal', 'sauria')
rownames(sp.class) <- sub(' ', '_', rownames(teleost.b) )

all(rownames(sp.class) == rownames( ex.align.2.k2 ) )
## [1] TRUE
all(rownames(sp.class) == colnames( ex.align.2.k2 ) )
## [1] TRUE
## this is very convenient.. 

sp.col <- rgb(sp.class[,1] * 0.5 + 0.5, sp.class[,2] * 0.5 + 0.5, sp.class[,3] * 0.5 + 0.5)
names(sp.col) <- rownames(sp.class)

sp.col.2 <- rgb(sp.class[,1] * 0.5 + 0.2, sp.class[,2] * 0.5 + 0.2, sp.class[,3] * 0.5 + 0.2)
names(sp.col.2) <- rownames(sp.class)

class.col <- sort(unique( sp.col ), decreasing=TRUE )
class.col.2 <- sort(unique( sp.col.2 ), decreasing=TRUE )

names(class.col) <- c('teleost', 'mammal', 'sauria', 'others')
names(class.col.2) <- c('teleost', 'mammal', 'sauria', 'others')

## the number of introns with size below min.s,
## 
n.below <- function(x, t.s=c(5,8)){
    y <- c( sapply(t.s, function(t){
        sum(x < t) }), length(x) )
    names(y) <- c(as.character(t.s), 'n')
    y
}

all.intron.s.n <- tapply( all.intron.s, all.gene.stats[,'sp'], function(x){
    n.below( log2(1+unlist(x))) })
all.intron.s.n <- t(sapply(all.intron.s.n, eval))
## reorder
all.intron.s.n <- all.intron.s.n[ rownames( ex.align.2.k2 ), ]

saveRDS(all.intron.s.n, 'all_intron_s_n.rds')

## 75 bp is a better measure than 32..
all.intron.s.n.2 <- tapply( all.intron.s, all.gene.stats[,'sp'], function(x){
    n.below( log2(1+unlist(x)), t.s=c(5, log2(75))) })
all.intron.s.n.2 <- t(sapply(all.intron.s.n.2, eval))
## reorder
all.intron.s.n.2 <- all.intron.s.n.2[ rownames( ex.align.2.k2 ), ]

plot(all.intron.s.n[,1] / all.intron.s.n[,3], col=sp.col, pch=19)
int.s.b <- (all.intron.s.n[,1] / all.intron.s.n[,3]) < 0.025
saveRDS( int.s.b, "int_s_b.rds" )

plot(all.intron.s.n[,1] / all.intron.s.n[,3], col=sp.col, pch=19,
     xlab='species', ylab='intron size < 33 bp')
abline(h=0.025, lty=2)
identify(1:length(sp.n), all.intron.s.n[,1] / all.intron.s.n[,3],
         names(sp.n), cex=0.8)
legend('topleft', names(class.col), pch=19, col=class.col)
plt <- recordPlot()

cairo_pdf("unreasonable_intron_size_props.pdf", width=0.9*pdf.m*a4.w, height=0.45*pdf.m*a4.w)
replayPlot(plt)
dev.off()

rownames( all.intron.s.n) == names( genome.n50 ) ## TRUE

### plot n50 vs unreasonable intron_size_props
plot( log10(1e4 + genome.n50), all.intron.s.n[,1] / all.intron.s.n[,3],
     col=sp.col, pch=19,
     ylab='intron size < 33 bp', xlab='log10(1e4 + N50)' )
abline(h=0.025, lty=2)
##
identify( log10(1e4 + genome.n50), all.intron.s.n[,1] / all.intron.s.n[,3],
         names(sp.n) )


######### To plot mutual information against genetic distance
set.alpha <- function(cols, alpha){
    cc <- col2rgb(cols)
    cc <- rgb( cc[1,]/255, cc[2,]/255, cc[3,]/255, alpha )
    names(cc) <- names(cols)
    cc
}

class.col.a <- set.alpha(class.col, 0.1)

pt.col <- matrix( rgb(1,1,1,0), nrow=nrow(ex.align.2.k2), ncol=ncol(ex.align.2.k2))
pt.col[ sp.class[,'mammal'] & int.s.b, sp.class[,'mammal'] & int.s.b ] <- class.col['mammal']
pt.col[ sp.class[,'teleost'] & int.s.b, sp.class[,'teleost'] & int.s.b ] <- class.col['teleost']
pt.col[ sp.class[,'sauria'] & int.s.b, sp.class[,'sauria'] & int.s.b ] <- class.col['sauria']

## do a similar setup, but here use specifically teleostei, eutheria and aves
## to be consistent with the first figure..
## if we use NA as the default colour then we can use set.alpha when calling
## ordered.plot or rand.plot
pt.col.2 <- matrix( NA, nrow=nrow(ex.align.2.k2), ncol=ncol(ex.align.2.k2))
pt.col.2[ sp.class.3[,'mammalia'] & int.s.b, sp.class.3[,'mammalia'] & int.s.b ] <- class.col.3['mammalia']
pt.col.2[ sp.class.3[,'eutheria'] & int.s.b, sp.class.3[,'eutheria'] & int.s.b ] <- class.col.3['eutheria']
pt.col.2[ sp.class.3[,'teleostei'] & int.s.b, sp.class.3[,'teleostei'] & int.s.b ] <- class.col.3['teleostei']
pt.col.2[ sp.class.3[,'sauria'] & int.s.b, sp.class.3[,'sauria'] & int.s.b ] <- class.col.3['sauria']
pt.col.2[ sp.class.3[,'aves'] & int.s.b, sp.class.3[,'aves'] & int.s.b ] <- class.col.3['aves']

rand.plot <- function(x, y, bg, col, ...){
    i <- sample(1:length(x))
    plot(x[i], y[i], bg=bg[i], col=col[i], ...)
}

ordered.plot <- function(x, y, bg, col, ...){
    i <- tapply( 1:length(col), col, eval )
    i <- i[ order( sapply(i, length), decreasing=TRUE )]
    plot(x, y, type='n', ...)
    for(j in i)
        points(x[j], y[j], bg=bg[j], col=col[j], ...)
}

sp.n <- rownames(sp.class)
names(sp.n) <- uc.1( sp.db2sp(sp.n) )

## The magic for exclusing bad data lies in the pt.col vectore
## which makes invisible points for bad points.
cairo_pdf("k2f_intron_mi.pdf", width=0.5 * a4.w * pdf.m, height=0.3 * a4.w * pdf.m)
par(mar=c(4.1, 4.1, 2.1, 2.1))
rand.plot( ex.align.2.k2, int.sp.mi.dr.m,
          xlab='Kimura two factor', ylab='Intron size mutual information',
          bg=pt.col, col=pt.col, pch=21, xlim=c(0, 0.3))
legend('topright', legend=c('Teleost', 'Mammal', 'Sauria'),
       pch=19, col=class.col[c('teleost', 'mammal', 'sauria')] )
dev.off()

par(mar=c(4.1, 4.1, 2.1, 2.1))
rand.plot( ex.align.2.k2, int.sp.mi.dr.m,
          xlab='Kimura two factor', ylab='Intron size mutual information',
          bg=pt.col, col=pt.col, pch=21, xlim=c(0, 0.3))
legend('topright', legend=c('Teleost', 'Mammal', 'Sauria'),
       pch=19, col=class.col[c('teleost', 'mammal', 'sauria')] )

rand.plot( ex.align.2.jc, int.sp.mi.dr.m,
     xlab='Kimura two factor', ylab='Intron size mutual information',
     bg=pt.col, col=pt.col, pch=21)

rand.plot( ex.align.2.jcg, int.sp.mi.dr.m,
     xlab='Kimura two factor', ylab='Intron size mutual information',
     bg=pt.col, col=pt.col, pch=21)

### Lets try to load up the trees and see if we can make some plots from them.
ex.align.2.k2.nj <- readRDS( "../R_172_genomes/ex_align_2_k2_nj.rds" )
require(ape)
plot(ex.align.2.k2.nj)
usr <- par("usr")
leaf.n <- length( with(ex.align.2.k2.nj, setdiff(edge[,2], edge[,1])) )

## named list containing;
## edge, edge.length, tip.label, Nnode
## the tip nodes are 1->172, and are present only in the second column of edge
## To draw this we can use a function edge.lines
source("/home/lmj/R/R_max_parsimony/R/functions.R")
tree.lines <- edge.lines( ex.align.2.k2.nj )

## it is then possible to draw:
##cairo_pdf("exon_based_nj_tree.pdf", width=a4.w*0.8*pdf.m, height=a4.h*pdf.m)
plot.new()
with(tree.lines, plot.window(xlim=c(0, 1.5 * max(x[,1:2])),
                             ylim=c(0, 1.2 * max(y)), xaxs='i', yaxs='i'))
with(tree.lines, segments(x[,1], y, x[,2], y))
with(tree.lines, segments(v[,1], v[,2], v[,1], v[,3]))
with(tree.lines, text( x[,2] + 0.001, y, labels=nodes, adj=c(0,0.5), cex=0.5))
with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    text( x[b,2] + 0.003, y[b], ex.align.2.k2.nj$tip.label[ nodes[b] ], cex=0.5, pos=4 )
    })
##dev.off()
## That might not be that useful. But we can use the y values for nodes 1->172 as
## the order for various plots.
## From that:
## 266 is the common ancestor of teleosts.
## 268 aves (birds)
## 256 other sauria (crocodile, turtles, tuatara
## 262 episquamata (snakes, lizards, etc)
## 247 sauria
## 239 mammalia
## 87 prototheria (Platypus)
## 237 theria (marsupials and placentals)
## 250 metatheria (marsupials)
## 215 eutheria (placentals)
## classify into:
## teleosts, sauria, aves, episquamata, other sauria,
## theria, metatheria, eutheria

nodes.descend <- function(tree, root){
    b <- rep(FALSE, length(unique(as.integer(tree$edge))))
    descend.tree <- function(root.i){
#        cat("setting: ", tree$edge[root.i,2], " TRUE\n")
        b[tree$edge[root.i,2]] <<- TRUE
        root.i <- which( tree$edge[,1] == tree$edge[root.i,2] )
        for(i in root.i)
            descend.tree(i)
    }
    b[root] <- TRUE
    root.i <- which( tree$edge[,1] == root )
    for(i in 1:length(root.i))
        descend.tree(root.i[i])
    b
}

mammalia.nodes <- nodes.descend( ex.align.2.k2.nj, 239 )
eutheria.nodes <- nodes.descend( ex.align.2.k2.nj, 215 )
teleost.nodes <- nodes.descend( ex.align.2.k2.nj, 266 )
aves.nodes <- nodes.descend( ex.align.2.k2.nj, 268 )
sauria.nodes <- nodes.descend( ex.align.2.k2.nj, 247 )

#nodes.col <- rgb(teleost.nodes | eutheria.nodes, aves.nodes, theria.nodes)
class.col.3 <- rgb( c(0, 0, 0.6, 0.8, 0, 0), c(0, 0, 0, 0, 0.8, 0.6), c(0, 0.8, 0.6, 0, 0, 0.6))
names(class.col.3) <- c('others', 'mammalia', 'eutheria', 'teleostei', 'sauria', 'aves')
## a more reasonable order:
class.col.3 <- class.col.3[ c(4, 5, 6, 2, 3, 1) ]
saveRDS( class.col.3, "class_col_3.rds")

nodes.col <- rep('black', length(mammalia.nodes) )
nodes.col[ mammalia.nodes ] <- rgb(0, 0, 0.8)
nodes.col[ eutheria.nodes ] <- rgb(0.6, 0, 0.6)
nodes.col[ teleost.nodes  ] <- rgb(0.8, 0, 0)
nodes.col[ sauria.nodes ] <- rgb(0, 0.8, 0)
nodes.col[ aves.nodes ] <- rgb(0, 0.6, 0.6)

## I do not remember exaclty why this works, but it is somewher in the details
## of the parsing of the tree.
sp.col.3 <- nodes.col[ 1:length(ex.align.2.k2.nj$tip.label) ]
names(sp.col.3) <- ex.align.2.k2.nj$tip.label
saveRDS( sp.col.3, "sp_col_3.rds" )

## We can in a similar manner make a class tree:
sp.class.3 <- cbind('teleostei'=teleost.nodes, 'mammalia'=mammalia.nodes,
                    'eutheria'=eutheria.nodes, 'aves'=aves.nodes,
                    'sauria'=sauria.nodes)[ 1:length(ex.align.2.k2.nj$tip.label), ]
rownames(sp.class.3) <- ex.align.2.k2.nj$tip.label
saveRDS( sp.class.3, "sp_class_3.rds" )

## it is then possible to draw:
cairo_pdf("exon_based_nj_tree.pdf", width=a4.w*0.8*pdf.m, height=a4.h*0.9*pdf.m)
par('mar'=c(0, 2.1, 0, 2.1))
plot.new()
with(tree.lines, plot.window(xlim=c(0, 1.5 * max(x[,1:2])),
                             ylim=c(0.25, 1.01 * max(y)), xaxs='i', yaxs='i'))
with(tree.lines, segments(x[,1], y, x[,2], y, col=nodes.col[x[,'p'] ] ))
with(tree.lines, segments(v[,1], v[,2], v[,1], v[,3], col=nodes.col[v[,'node']] ))
with(tree.lines, text( x[,2] + strwidth("9", cex=0.6), y, labels=nodes, adj=c(0,0.5), cex=0.6))
with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    offset <- strwidth("9999", cex=0.6) * 1.25
    text( x[b,2] + offset, y[b], uc.1(sp.db2sp(ex.align.2.k2.nj$tip.label[ nodes[b] ])),
         cex=0.6, adj=c(0,0.5) )
    })
legend('bottomleft', legend=names(class.col.3), text.col=class.col.3, box.lty=0,
       inset=c(0,0.05))
dev.off()


# then we can make the y positions
sp.y <- with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    lab <- ex.align.2.k2.nj$tip.label[ nodes[b] ]
    y.pos <- y[b]
    names(y.pos) <- lab
    y.pos
})

saveRDS( sp.y, "sp_y.rds" )

## The inferred intron lengths are in:

int.s.inf <- readRDS("../R_172_genomes/int_s_inf.rds") ## which is a big data set..
## and the changes in intron size and genetic distances
leaf.int.d <- readRDS("../R_172_genomes/leaf_int_d.rds")

## there are a lot of ways in which we can look at this data.
## One way would be to plot the accumulated delta for all lineages;
## One would hope that would separate the groups in some reasonable manner..
## Also nice would be to plot phylogenetic edge length vs delta or abs(delta)
## to see what outliers we have.

## The internal nodes with the biggest reduction in size is 266, which is the
## teleostei common node. This is followed by the node that joins up
## the puffer fishes and some other teleost internal nodes. I believe that this
## probably argues for a general decrease in intron size in teleosts.

### traverse the leaf nodes to get the distance. Use both the tree (for
### phylogenetic differences) and a summarised change in intron size
### We have a total of leaf.n leaf nodes..

## this differs from the version in ../R_172_genomes/ in that we only
## use reasonable intron sizes
leaf.int.d.2 <- lapply( 1:leaf.n, function(i){
    pg.dist <- c()
    int.d <- c()
    node.i <- c()
    ascend.tree <- function(j){
        b <- ex.align.2.k2.nj$edge[,2] == j
        if(sum(b) == 1){
            pg.dist <<- c(ex.align.2.k2.nj$edge.length[b], pg.dist)
            p.j <- ex.align.2.k2.nj$edge[b,1]
            b2 <- (int.s.inf[[j]]$state[,1] >= 6 & int.s.inf[[p.j]]$state[,1] >= 6)
            int.d <<- c(mean(int.s.inf[[j]]$state[b2,2]), int.d)
            node.i <<- c(j, node.i)
            ascend.tree( p.j )
        }
    }
    ascend.tree(i)
    cbind('pg'=cumsum(pg.dist),
          'int.d'=cumsum(int.d),
          'node'=node.i)
})
names(leaf.int.d.2) <- names( leaf.int.d )


draw.delta.lines <- function(m, cols=nodes.col, ...){
    x <- c(0, m[,'pg'])
    y <- c(0, m[,'int.d'])
    nodes <- c(0, m[,'node'])
    ## tel.b <- cumsum( m[,'node'] == 266 )
    ## mammal.b <- cumsum( m[,'node'] == 215 ) ## does not include marsupials
    ## aves.b <- cumsum( m[,'node'] == 268 )
    ## if(is.na(cols))
    ##     cols <- rgb(tel.b, aves.b, mammal.b)
    l <- length(x)
    segments( x[1:(l-1)], y[1:(l-1)], x[2:l], y[2:l], col=cols[nodes[2:l] ], ... )
    nr <- nrow(x)
    c('x'=x[l], 'y'=y[l])
}

xlim <- range( unlist(lapply( leaf.int.d.2, function(x){ x[,'pg'] })) )
xlim[2] <- 1.2 * xlim[2]
xlim[1] <- 0
ylim <- range( unlist(lapply( leaf.int.d.2, function(x){ x[,'int.d'] })) )

## we can put labels at places that are of interest
label.i <- which( names(leaf.int.d) %in%
                  c('danio rerio', 'tetraodon nigroviridis', 'betta splendens',
                    'astyanax mexicanus', 'erpetoichthys calabaricus', 'latimeria chalumnae',
                    'sphenodon punctatus', 'pygocentrus nattereri','eptatretus burgeri',
                    'takifugu rubripes', 'xenopus tropicalis', 'callorhinchus milii',
                    'lepisosteus oculatus', 'gasterosteus aculeatus', 'parambassis ranga',
                    'electrophorus electricus', 'denticeps clupeoides'))

cairo_pdf("Intron_size_evolution.pdf", width=a4.w*0.5*pdf.m, height=a4.h*0.75*pdf.m)
plot.new()
plot.window(xlim=xlim, ylim=ylim, xaxs='i')
abline(h=0, col='grey')
ends <- sapply(leaf.int.d, draw.delta.lines, lwd=2)
axis(2)
with( leaf.int.d, {
    o <- order( ends[2,label.i] )
    x <- (ends[1,label.i])[o]
    y <- sort(ends[2,label.i])
    y.d <- diff(y)
    y.d.i <- which(y.d < strheight( 'A', cex=0.8 ) )
    y[ y.d.i ] <- y[y.d.i] - strheight('A') * 0.25 
    text( x, y, uc.1(colnames(ends)[label.i])[o], pos=4,
     col=nodes.col[label.i[o]], cex=0.8) })
legend('bottomleft', legend=names(class.col.3), text.col=class.col.3, box.lty=0,
       inset=c(0,0.05))
dev.off()

cairo_pdf("Intron_size_evolution_2.pdf", width=a4.w*0.5*pdf.m, height=a4.h*0.75*pdf.m)
plot.new()
plot.window(xlim=xlim, ylim=ylim, xaxs='i')
abline(h=0, col='grey')
ends <- sapply(leaf.int.d.2, draw.delta.lines, lwd=2)
axis(2)
with( leaf.int.d.2, {
    o <- order( ends[2,label.i] )
    x <- (ends[1,label.i])[o]
    y <- sort(ends[2,label.i])
    y.d <- diff(y)
    y.d.i <- which(y.d < strheight( 'A', cex=0.8 ) )
    y[ y.d.i ] <- y[y.d.i] - strheight('A') * 0.25 
    text( x, y, uc.1(colnames(ends)[label.i])[o], pos=4,
     col=nodes.col[label.i[o]], cex=0.8) })
legend('bottomleft', legend=names(class.col.3), text.col=class.col.3, box.lty=0,
       inset=c(0,0.05))
dev.off()


identify( ends[1,], ends[2,], labels=names(sp.n) )

### We can look at the changes at the leaf nodes to see if there is evidence for
### for continuing changes. That would be better done by individual distributions
### for the underlying numbers, but we can first have a look at the mean changes:

leaf.d <- t(sapply( int.s.inf[1:length(leaf.int.d)], function(x){
    b <- x$state[,1] >= 6
    c( mean( x$state[b,2]), median( x$state[b,2] ))
}))
colnames(leaf.d) <- c('mean', 'median')
rownames(leaf.d) <- names(leaf.int.d)

par(mar=c(5.1, 10.1, 4.1, 2.1))
barplot( leaf.d[,'mean'], col=(nodes.col[ 1:length(leaf.d) ]), horiz=T, las=2, cex.names=0.3)
barplot( leaf.d[int.s.b,'mean'], col=(nodes.col[ 1:length(leaf.d) ])[int.s.b], horiz=T, las=2, cex.names=0.3)

## let us do a manual bar plot so that we can add some information more easily
sp.lab <- strsplit(rownames(leaf.d), ' ')
sp.lab <- sapply( sp.lab, function(x){ paste( toupper(substr(x[1], 1, 1)), ". ", x[2], sep="") })

cairo_pdf( 'leaf_mean_intron_changes.pdf', width=a4.w * pdf.m, height=a4.h * pdf.m )
par(omi=c(0.5, 0.5, 0.5, 0.5)) 
par(mar=c(5.1, 4.1, 4.1, 0.1))
layout(matrix(1:2, nrow=1), widths=c(15, 90))
plot.new()
plot.window( xlim=0:1, ylim=c(0,nrow(leaf.d)), xaxs='i')
text(1, 1:nrow(leaf.d)-0.5, sp.lab, cex=0.5, adj=c(1,0.5), col=ifelse(int.s.b, 'black', 'red'))
par(mar=c(5.1, 0.1, 4.1, 2.1))
plot.new()
plot.window( xlim=range(leaf.d[,'mean']), ylim=c(0,nrow(leaf.d)), xaxs='i')
y <- 1:nrow(leaf.d)
with(par(), rect(usr[1], y-1, usr[2], y, col=c(rgb(1, 1, 1), rgb(0.9, 0.9, 0.9)), border=NA ))
with(par(), rect(usr[1], y[!int.s.b]-1, usr[2], y[!int.s.b], col=rgb(1, 0.7, 0.7, 0.3), border=NA))
rect( leaf.d[,'mean'], y-1, 0, y, col=nodes.col[ 1:length(leaf.d) ] )
axis(1)
dev.off()

## median is not interesting as it is almost always 0.
barplot( leaf.d[int.s.b,'median'], col=(nodes.col[ 1:length(leaf.d) ])[int.s.b], horiz=T, las=2, cex.names=0.3)

## more reasonable to go back to the original inferred states. The good thing here is
## that the first 172 nodes are the leaf nodes and have the approprate names.

leaf.d.table <- t(sapply( int.s.inf[1:length(leaf.int.d)], function(x){
    b <- x$state[,1] >= 6
    table( c(-12:15, x$state[b,2]) )
}))
rownames(leaf.d.table) <- names(leaf.int.d)

image( log(leaf.d.table[ int.s.b, ]))

## this gives us the numbers to put in above to get a non-ragged matrix
range(unlist(sapply(leaf.d.table, function(x){ as.numeric(names(x))})))

### lets look at the mutual information and how it distributes across the data set
y <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=TRUE )
x <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=FALSE )


sp.o <- order(sp.y[sp.n])
tmp <- int.sp.mi.dr.m[sp.o, sp.o]
par(mfrow=c(1,1))
cm <- 0.3  ## something to moderate the color
cairo_pdf("species_mutual_information_1.pdf", height=25, width=25)
par(mar=c(1.1,1.1,1.1,1.1))
plot.new()
plot.window(xlim=c(-15,length(sp.n)), ylim=c(1, 15+length(sp.n)))
## rect( x, y, x+1, y+1,
##      col=hsvScale(int.sp.mi.dr.m, val=(cm + int.sp.mi.dr.m)/(cm + max(int.sp.mi.dr.m, na.rm=TRUE))  ),
##      border='grey', lwd=0.25 )
rect( x, y, x+1, y+1,
     col=hsvScale(tmp, val=(cm + tmp)/(cm + max(tmp, na.rm=TRUE))  ),
     border='grey', lwd=0.25 )
## text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", rownames(int.sp.mi.dr.m) ), adj=c(1,0.5), cex=0.8)
## text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", colnames(int.sp.mi.dr.m)), adj=c(0,0.5), srt=90, cex=1)
text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", rownames(tmp) ), adj=c(1,0.5), cex=0.8,
     col=nodes.col[sp.o])
text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", colnames(tmp)), adj=c(0,0.5),
     srt=90, cex=1, col=nodes.col[sp.o])
dev.off()

## lets do this for the subset specified by int.s.b which is defined by the
## proportion of unreasonably small introns
## let us make this into a function

mi2col <- function(v, cm){
    hsvScale(v, val=(cm + v)/(cm + max(v, na.rm=TRUE)))
}

plot.mutual.info <- function(mutual.info=int.sp.mi.dr.m, o.y=sp.y, sp.b=int.s.b, cm=0.3,
                             col=nodes.col[1:nrow(mutual.info)], class.col=class.col.3,
                             border.col=NA, col.f=mi2col, text.labels=TRUE,
                             lab.cex=1, scale.cex=lab.cex, y.margins=c(0,15), x.margins=c(15,0),
                             class.bar.w=2){
    names(col) <- rownames(mutual.info)
    mi <- mutual.info[ sp.b[rownames(mutual.info)], ]
    mi <- mi[ , sp.b[colnames(mi)] ]
    sp.n <- colnames(mi)
    sp.o <- order( o.y[ sp.n ] )
    mi <- mi[ sp.o, sp.o ]
    sp.n <- sp.n[sp.o]
    y <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=TRUE )
    x <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=FALSE )
    par(mar=c(1.1,1.1,1.1,1.1))
    plot.new()
    plot.window(xlim=c(-x.margins[1],length(sp.n)+x.margins[2]),
                ylim=c(-y.margins[1], length(sp.n)+y.margins[2]))
    rect( x, y, x+1, y+1,
     col=col.f(mi, cm),
     border=border.col, lwd=0.1 )
    if(text.labels){
        labels <- uc.1( sp.db2sp( colnames(mi) ))
        text( 0.95, y[1,] + 0.5, labels, adj=c(1,0.5), cex=lab.cex,
             col=col[rownames(mi)])
        text( x[,1]+0.5, length(sp.n) + 1.05, labels, adj=c(0,0.5),
             srt=90, cex=lab.cex, col=col[colnames(mi)])
    }else{
        rect(-0.5-class.bar.w, y[1,], -0.5, y[1,]+1, col=col[ rownames(mi) ], border=NA)
        rect(x[,1], length(sp.n)+1.5, x[,1]+1, length(sp.n)+1.5+class.bar.w,
             col=col[ colnames(mi) ],
             border=NA)
        ## lets put a scale bar at the left top.
        ## we should have this regardless, but it is more difficult to fit with text labels.
        scale.y <- 1:(length(sp.n)/4) ## + length(sp.n)/2
        scale.x <- -6
        mi.range <- range(mutual.info, na.rm=TRUE)
        scale.v <- seq(mi.range[1], mi.range[2], length.out=length(scale.y))
        rect( scale.x, scale.y, scale.x+2, scale.y+1, col=col.f(scale.v, cm), border=NA )
        tick.i <- seq(1, length(scale.y), length.out=5)
        tick.label <- sprintf("%1.1e", scale.v[tick.i])
        text(scale.x, scale.y[tick.i]+0.5, tick.label, adj=c(1,0.5), cex=scale.cex)
        txt.h <- strheight("A", cex=lab.cex)
        txt.y <- length(sp.n)
        A.w <- strheight("A", cex=lab.cex)
        invisible(sapply(1:length(class.col), function(i){
            text(-0.5-class.bar.w-A.w/3, txt.y - 1.5 * txt.h * i, names(class.col)[i],
                 adj=c(1,0), col=class.col[i], cex=lab.cex )
        }))
        ## legend('topleft', names(class.col), text.col=class.col, bty='n',
        ##        inset=c(0,0.1))
    }
    invisible( list('mi'=mi, 'col'=col.f(mi, cm)) )
}

cairo_pdf("species_mutual_information_2.pdf", height=25, width=25)
mi.cols <- plot.mutual.info(text.labels=TRUE)
dev.off()

cairo_pdf("species_mutual_information_3.pdf", height=10, width=10)
mi.cols <- plot.mutual.info(text.labels=FALSE, lab.cex=0.8)
dev.off()

### that looks OK. We will need to put together in a figure; but before that we
### need to ask what level of mutual information is less than expected. For that
### we need to run a whole load of permutations...
### Lets first see how much time this takes..
intron.orth.l <- read.table("../R_172_genomes/dr_intron_orthology_l.txt",
                            sep="\t", stringsAsFactors=FALSE, header=TRUE )
colnames(intron.orth.l) <- sub(".", "_", colnames(intron.orth.l), fixed=TRUE)

intron.orth.id <- read.table("../R_172_genomes/dr_intron_orthology_id.txt",
                            sep="\t", stringsAsFactors=FALSE, header=TRUE )
colnames(intron.orth.id) <- sub(".", "_", colnames(intron.orth.id), fixed=TRUE)


require('entropy')

system.time(
    tmp <- mutual.info(intron.orth.l[,1], intron.orth.l[,2],
                       numBins=20)
)
## felt like that took longer.. but..
## that is probably because the system is very heavily used at the moment
##  user  system elapsed 
## 0.105   0.003   0.139

## this should be parallelisable
require(parallel)

## first let us get the set of real mi values
int.sp.mi <- unlist(sapply( 1:(nrow(int.sp.mi.dr.m)-1), function(i){
    int.sp.mi.dr.m[i, (i+1):ncol(int.sp.mi.dr.m)] }))

mi.perm <- mclapply(1:length(int.sp.mi), function(i){
    sp.i <- sample(1:ncol(intron.orth.l), size=2)
    list(sp=sp.i, mi=mutual.info( intron.orth.l[,sp.i[1]],
                                  sample(intron.orth.l[,sp.i[2]])))
    }, mc.cores=10 )

## that didn't take that much time at all.
## let us harvest the mutual information values
mi.perm.mi <- sapply( mi.perm, function(x){ x$mi$mi })


par(mfrow=c(1,2))
all.mi.h <- hist(c(int.sp.mi, mi.perm.mi))

int.sp.mi.h <- hist(int.sp.mi, breaks=all.mi.h$breaks)
mi.perm.mi.h <- hist(mi.perm.mi, breaks=all.mi.h$breaks)

## anyway, let's plot those two
plot( mi.perm.mi.h$mids, mi.perm.mi.h$density, type='l', col='black', lwd=2 )
lines( int.sp.mi.h$mids, int.sp.mi.h$density, type='l', col='red', lwd=2 )

## more informatively
## I think these two for supplementary as a single figure
cairo_pdf("Mutual_info_distribution.pdf", width=a4.w*pdf.m*0.9, height=a4.w*pdf.m*0.45 )
par(mfrow=c(1,2))
plot( int.sp.mi.h, main='', xlab='Mutual information' )
plot( sort(mi.perm.mi), type='p', ylab='Mutual information', cex=0.5 )
abline(h=min(int.sp.mi), col='red', lty=2)
dev.off()

plot( sort(mi.perm.mi),
     sort(int.sp.mi), type='p', ylab='mi', cex=0.5  )

mi.perm.mi.qnt <- quantile( mi.perm.mi, probs=seq(0,1,0.01) )
## 99% are below 0.006880976
## more reasonable to ask
sum( mi.perm.mi >= min(int.sp.mi) ) / length( mi.perm.mi )
## 0.0002039984  (2/10000)

## I think that the mutual information figure needs only to have
## the image and the mutual information vs kimura two factor

cairo_pdf("Mutual_info_figure.pdf", width=a4.w * 0.9 * pdf.m, height=a4.w * 0.45 * pdf.m )
layout(matrix(c(1,2), nrow=1), widths=c(1,1))
mi.cols <- plot.mutual.info(text.labels=FALSE, lab.cex=0.75, scale.cex=0.6,
                            y.margins=c(8,5))
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=-1))
par(mar=c(4.1, 4.1, 2.1, 2.1))
rand.plot( ex.align.2.k2, int.sp.mi.dr.m,
          xlab='Kimura two factor', ylab='Intron size mutual information',
          bg=pt.col, col=pt.col, pch=21, xlim=c(0, 0.3))
legend('topright', legend=c('Teleost', 'Mammal', 'Sauria'),
       pch=19, col=class.col[c('teleost', 'mammal', 'sauria')] )
with(par(), mtext('B', at=usr[1], cex=mt.cex))
dev.off()

cairo_pdf("Mutual_info_figure_2.pdf", width=full.w * pdf.m, height=full.w * 0.45 * pdf.m )
layout(matrix(c(1,2), nrow=1), widths=c(1,1))
mi.cols <- plot.mutual.info(text.labels=FALSE, lab.cex=0.75, scale.cex=0.6,
                            y.margins=c(8,5), x.margins=c(18, 0))
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=-1))
##
par(mar=c(4.1, 4.1, 2.1, 2.1))
ordered.plot( ex.align.2.k2, int.sp.mi.dr.m,
          xlab='Kimura two factor', ylab='Intron size mutual information',
          bg=set.alpha(pt.col.2, 0.25), col=pt.col.2, pch=21, xlim=c(0, 0.3))
legend('topright', legend=c('Teleost', 'Mammal', 'Eutheria', 'Sauria', 'Aves'),
       pch=19, col=class.col.3[c('teleostei', 'mammalia', 'eutheria', 'sauria', 'aves')] )
with(par(), mtext('B', at=usr[1], cex=mt.cex))
dev.off()


### replot the nreasonable_intron_size_props.pdf using a taxonomy based order
### of species.
par('mfrow'=c(1,1))
x <- 1:nrow(all.intron.s.n)
o <- order( sp.y[ rownames(all.intron.s.n) ] )
plot(x, all.intron.s.n[o,1] / all.intron.s.n[o,3], col=sp.col[o], pch=19,
     xlab='species', ylab='intron size < 33 bp')
abline(h=0.025, lty=2)
pos.1 <- identify(x, all.intron.s.n[o,1] / all.intron.s.n[o,3],
         names(sp.n[o]), cex=0.8, pos=TRUE)
legend('topright', names(class.col), pch=19, col=class.col)
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=-1))

### plot n50 vs unreasonable intron_size_props
plot( log10(1e4 + genome.n50), all.intron.s.n[,1] / all.intron.s.n[,3],
     col=sp.col, pch=19,
     ylab='intron size < 33 bp', xlab='log10(1e4 + N50)' )
abline(h=0.025, lty=2)
##
pos.2 <- identify( log10(1e4 + genome.n50), all.intron.s.n[,1] / all.intron.s.n[,3],
         names(sp.n), pos=TRUE )
legend('topright', names(class.col), pch=19, col=class.col)
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=-1))


cairo_pdf("unreasonable_intron_size_props_2.pdf", width=0.9*pdf.m*a4.w, height=0.9*pdf.m*a4.w)
par('mfrow'=c(2,1))
x <- 1:nrow(all.intron.s.n)
o <- order( sp.y[ rownames(all.intron.s.n) ] )
plot(x, all.intron.s.n[o,1] / all.intron.s.n[o,3], col=sp.col[o], pch=19,
     xlab='species', ylab='intron size < 33 bp')
abline(h=0.025, lty=2)
text(x[pos.1$ind], (all.intron.s.n[o,1] / all.intron.s.n[o,3])[pos.1$ind],
         names(sp.n[o])[pos.1$ind], cex=0.8, pos=pos.1$pos)
legend('topright', names(class.col), pch=19, col=class.col)
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))
##
### plot n50 vs unreasonable intron_size_props
plot( log10(1e4 + genome.n50), all.intron.s.n[,1] / all.intron.s.n[,3],
     col=sp.col, pch=19,
     ylab='intron size < 33 bp', xlab='log10(1e4 + N50)' )
abline(h=0.025, lty=2)
##
text( log10(1e4 + genome.n50)[pos.2$ind], (all.intron.s.n[,1] / all.intron.s.n[,3])[pos.2$ind],
         names(sp.n)[pos.2$ind], pos=pos.2$pos, cex=0.8 )
legend('topright', names(class.col), pch=19, col=class.col)
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
dev.off()

cairo_pdf("unreasonable_intron_size_props_3.pdf", width=0.9*pdf.m*a4.w, height=0.9*pdf.m*a4.w)

par('mfrow'=c(2,1))
x <- 1:nrow(all.intron.s.n.2)
o <- order( sp.y[ rownames(all.intron.s.n.2) ] )
plot(x, all.intron.s.n.2[o,2] / all.intron.s.n.2[o,3], col=sp.col[o], pch=19,
     xlab='species', ylab='intron size < 75 bp')
abline(h=0.05, lty=2)
text(x[pos.1$ind], (all.intron.s.n.2[o,2] / all.intron.s.n.2[o,3])[pos.1$ind],
         names(sp.n[o])[pos.1$ind], cex=0.8, pos=pos.1$pos)
legend('topright', names(class.col), pch=19, col=class.col)
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))

##
### plot n50 vs unreasonable intron_size_props
plot( log10(1e4 + genome.n50), all.intron.s.n.2[,2] / all.intron.s.n.2[,3],
     col=sp.col, pch=19,
     ylab='intron size < 33 bp', xlab='log10(1e4 + N50)' )
abline(h=0.05, lty=2)
##
text( log10(1e4 + genome.n50)[pos.2$ind], (all.intron.s.n.2[,2] / all.intron.s.n.2[,3])[pos.2$ind],
         names(sp.n)[pos.2$ind], pos=pos.2$pos, cex=0.8 )
legend('topright', names(class.col), pch=19, col=class.col)
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))

dev.off()


## ask questions about conservation across large evoluionary distances
## for this we can start off with considering the sizes across different
## clades.

## it seems that I mixed up the column names when calculting the variances in the
## the main analysis. So I need to redo these here.
var.par <- function(x){
    c('n'=sum(!is.na(x)), 'mean'=mean(x, na.rm=TRUE), 'sd'=sd(x, na.rm=TRUE),
      'min'=min(x, na.rm=TRUE), 'max'=max(x, na.rm=TRUE), quantile(x, probs=seq(0,1,0.1), na.rm=TRUE))
}

## we should use int.s.b here
all(names(int.s.b) == rownames(sp.class))  ## well that is true at least..
## 

tel.var <- t(apply( log2(intron.orth.l[ ,rownames(sp.class)[sp.class[,'teleost'] & int.s.b]]),
                   1, var.par ))

mam.var <- t(apply( log2(intron.orth.l[ ,rownames(sp.class)[sp.class[,'mammal'] & int.s.b]]),
                   1, var.par ))

sau.var <- t(apply( log2(intron.orth.l[ ,rownames(sp.class)[sp.class[,'sauria'] & int.s.b]]),
                   1, var.par ))

plot(tel.var[,'min'], mam.var[,'min'])
plot(tel.var[,'min'], sau.var[,'min'])
plot(sau.var[,'min'], mam.var[,'min'])


plot(tel.var[,'10%'], mam.var[,'10%'])
plot(tel.var[,'10%'], sau.var[,'10%'])
plot(sau.var[,'10%'], mam.var[,'10%'])

hist(tel.var[,'10%'])
hist(mam.var[,'10%'])
hist(tel.var[,'50%'])
hist(mam.var[,'50%'])

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

tel.mam.min <- hist.2d( tel.var[,'min'], mam.var[,'min'], numBins=20 )
tel.mam.10 <-  hist.2d( tel.var[,'10%'], mam.var[,'10%'], numBins=20 )
tel.mam.50 <-  hist.2d( tel.var[,'50%'], mam.var[,'50%'], numBins=20 )
tel.mam.10.50 <-  hist.2d( tel.var[,'10%'], mam.var[,'50%'], numBins=20 )

image(tel.mam.10$cb, tel.mam.10$rb, log(tel.mam.10$h))
image(tel.mam.10$cb, tel.mam.10$rb, tel.mam.10$h)

with(tel.mam.10, image(rb, cb, h))
with(tel.mam.10, image(rb, cb, scale((h)), xlab='Teleost', ylab='Mammal' ) )
## scale the rows, so that we can see how teleost size
## predicts mammal size.
with(tel.mam.10, image(rb, cb, t(scale((t(h)))), xlab='Teleost', ylab='Mammal' ) )
with(tel.mam.10, image(rb, cb, log((h)), xlab='Teleost', ylab='Mammal' ) )
plot(tel.var[,'10%'], mam.var[,'10%'], cex=0.5 )

with(tel.mam.50, image(rb, cb, t(scale((t(h)))), xlab='Teleost', ylab='Mammal' ) )
with(tel.mam.50, image(rb, cb, log((h)), xlab='Teleost', ylab='Mammal' ) )
plot(tel.var[,'50%'], mam.var[,'50%'], cex=0.5 )


with(tel.mam.10, image(rb, cb, t(scale(t(h)))))
with(tel.mam.50, image(cb, rb, h))
with(tel.mam.50, image(cb, rb, log(h)))

with(tel.mam.10.50, image(cb, rb, t(scale(t(h)))))

with(tel.mam.10, plot(rb[-1] - diff(rb)/2, rowSums(h) ))
with(tel.mam.10, points(cb[-1] - diff(cb)/2, colSums(h), col='red' ))
## this says that the rows are the teleosts..

sau.mam.10 <-  hist.2d( sau.var[,'10%'], mam.var[,'10%'], numBins=20 )
with(sau.mam.10, image(rb, cb, h)) ## !!!
with(sau.mam.10, image(rb, cb, t(scale(t(h)))) ) ## !!!

tel.sau.10 <-  hist.2d( tel.var[,'10%'], sau.var[,'10%'], numBins=20 )
with(tel.sau.10, image(rb, cb, h)) ## !!!
with(tel.sau.10, image(rb, cb, t(scale(t(h)))) ) ## !!!

v2col <- function(v, cm){
    vv= v - min(v)
    val <- (cm + vv)/(cm + max(vv, na.rm=TRUE))
    hsvScale(v, val=val)
}

v2hcl <- function(v, beg=240, end=360, min.l=0, max.l=100, chroma=75){
    vv <- (v - min(v)) / diff( range(v, na.rm=TRUE))
    hue <- beg + (end - beg) * vv
    luminance <- min.l + (max.l - min.l) * vv
##    hsv( cm, val, val )
    hcl( hue, chroma, luminance )
}

## OK. We will make a plot with these types of values:
## (cm is a colour moderator that moderates the colour value)
## x and y are the breaks of the rectangles, (x[-1] == x2)
## (inner) margins are bottom, left, top, right
rect_image <- function(m, x=0:ncol(m), y=0:nrow(m),
                       col.f=v2col, col.args=NULL, margins=rep(0, 4),
                       border=NA, clear=TRUE, axis=TRUE,
                       xlab=NULL, ylab=NULL, lab.cex=1.5){
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
        cols <- col.f(m)
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


with(tel.mam.10, rect_image(scale(t(h)), cb, rb, col.f=v2hcl,
                            col.args=list(beg=0, end=360, chroma=50, min.l=0)) )
with(tel.mam.10, rect_image(scale(t(h)), cb, rb, col.f=v2sv, cm=0))

## this is supposed to be bad, but it still looks the best of the above alternatives
## but whatever.. we can make better colours at another time
with(tel.mam.10, rect_image(scale(t(h)), cb, rb, col.f=v2col, col.args=list(cm=0.3),
                            xlab='Teleost', ylab='Mammal'))

## the median is even better..
with(tel.mam.50, rect_image(scale(t(h)), cb, rb, col.f=v2col, col.args=list(cm=0.3),
                            xlab='Teleost', ylab='Mammal'))

### I suggest we plot:
## scatter plot (with alpha blending), image with no scaling and scaled image
par(mfrow=c(1,3))
plot( tel.var[,'50%'], mam.var[,'50%'], cex=0.5, xaxs='i', yaxs='i',
     xlab='', ylab='')
mtext('Teleost', side=1, line=2.5)
mtext('Mammal', side=2, line=2.5)
with(tel.mam.50, rect_image( t(log(1 + h)), cb, rb, col.f=v2col, col.args=list(cm=0.3),
                            xlab='Teleost', ylab='Mammal', lab.cex=1))
with(tel.mam.50, rect_image(scale(t(h)), cb, rb, col.f=v2col, col.args=list(cm=0.3),
                            xlab='Teleost', ylab='Mammal', lab.cex=1))
### we can ask:
### as the (quantile) teleost size increases what proportion of introns are above
### median (quantile) mammal length


long.prediction <- function(m1, m2, q1, q2, m2.q=0.5, decreasing=TRUE){
    v1 <- m1[,q1]
    v2 <- m2[,q2]
    b <- !is.na(v1) & !is.na(v2)
    v1 <- v1[b]
    v2 <- v2[b]
    o <- order(v1, decreasing=decreasing)
    v1 <- v1[o]
    v2 <- v2[o]
    v2.m <- sort(v2)[ length(v2) * m2.q ]
    v2.b <- v2 >= v2.m
    v2.cs <- cumsum(v2.b)
    v2.exp <- 1:length(v2) * (sum(v2.b) / length(v2.b))
    v2.fraction <- v2.cs / 1:length(v2)
    if(decreasing)
        v2.p <- phyper( v2.cs-1, sum(v2.b), length(v2.b) - sum(v2.b), 1:length(v2), lower.tail=FALSE )
    else
        v2.p <- phyper( v2.cs, sum(v2.b), length(v2.b) - sum(v2.b), 1:length(v2) )
    list(v1=v1, v2=v2, m=v2.m, cs=v2.cs, exp=v2.exp, fraction=v2.fraction, p=v2.p)
}

long.pred <- lapply( seq(0.5, 0.95, 0.05), function(x){
    long.prediction(tel.var, mam.var, '50%', '50%', m2.q=x) })
names(long.pred) <- sprintf("%.0f%%", seq(0.5, 0.95, 0.05) * 100)

plot.pred.1 <- function(pred=long.pred, use.index=FALSE, use.log=TRUE,
                        plot.exp.r=TRUE, ...){
    exp <- sapply(pred, function(x){ x$exp })
    obs <- sapply(pred, function(x){ x$cs })
    if(use.index)
        v1 <- sapply(pred, function(x){ (length(x$v1):1)/length(x$v1) })
    else
        v1 <- sapply(pred, function(x){ x$v1 })
    if(plot.exp.r)
        ratio <- obs / exp
    else
        ratio <- obs / 1:nrow(obs)
    if(use.log)
        ratio=log2(ratio)
    cols <- hsvScale(1:ncol(ratio), sat=1, val=0.7)
    plot(v1[,1], ratio[,1], xlim=range(v1), ylim=range(ratio),
         type='l', col=cols[1], ...)
    for(i in 2:ncol(ratio))
        lines(v1[,i], ratio[,i], col=cols[i])
    invisible(list('x'=v1, 'y'=ratio, cols=cols))
}

plot.pred.1(use.index=FALSE, xlab='log2 teleost median length', ylab='obs / expected above mammal quantile')
plot.pred.1(use.index=TRUE, xlab='teleost quantile', ylab='obs / expected  above mammal quantile')



par(mfrow=c(1,3))
b <- !is.na(mam.var[,'50%'])
plot(1:sum(b) / sum(b), sort(mam.var[b,'50%']), type='p', xlab='quantile', ylab='log2 intron length' )
b <- !is.na(tel.var[,'50%'])
points(1:sum(b) / sum(b), sort(tel.var[b,'50%']), type='p', col='red' )
#
## let us have some values
i <- sum(tel.var[b,'50%'] <= 10) / sum(b)
with(par(), lines( c(i,i,usr[1], usr[2]), c(usr[3], 10, 10, 10) ))
legend('topleft', c('mammalia', 'teleostia'), pch=1, col=c('black', 'red'))
with(par(), text( usr[1], 10, sprintf("%.1f %% <= 10", 100*(i)), adj=c(-0.1, 1.5), col='red'))
with(par(), text( usr[2], 10, sprintf("%.1f %% > 10", 100*(1-i)), adj=c(1.1, 1.5), col='red' ))
#
##
tmp <- plot.pred.1(use.index=FALSE, xlab='log2 teleost median length', ylab='fraction above mammal quantile', plot.exp.r=FALSE, use.log=FALSE)
i <- sum(tmp$x[,1] >= 10)
i.y <- max( tmp$y[i,] )
with(par(), lines( c(10, 10, usr[2]), c(usr[3], i.y, i.y), col=tmp$cols[1]), )
with(par(), text(usr[2], i.y - 1:ncol(tmp$y) * strheight('A') * 1.5,
                 sprintf('%.1f%%', 100 * tmp$y[i,]), adj=c(1.2,1), col=tmp$cols) )
legend('topleft', legend=names(long.pred), text.col=tmp$col, lty=1, col=tmp$col,
       title='mammalian\nquantiles')
##
plot.pred.1(use.index=TRUE, xlab='teleost median quantile', ylab='log2 obs-expected ratio', plot.exp.r=TRUE, use.log=TRUE)
legend('topleft', legend=names(long.pred), text.col=tmp$col, lty=1, col=tmp$col,
       title='mammalian\nquantiles')

## to look at the minimal intron length:
## (this would be better in basic stats, but since I have the figure here)
estimate.min.length <- function(c.var, param="50%", r=1:100){
    v <- sort(c.var[,param])
    v.d1 <- diff(2^v)
    v.d2 <- diff(v.d1)
    par(mfrow=c(1,1))
    plot(v[r], type='l')
#    plot(2^v[r], v.d1[r], type='l')
#    plot(2^v[r], v.d2[r], type='l')
}

## from playing around with the above it looks like 75 or 76 base pairs
## is a reasonable lower limit.

## That seems ok. Lets try to make a figure together with
## the other plots.

cairo_pdf("quantiles_prediction_figure.pdf", width=full.w * pdf.m, height=full.w * 0.6 * pdf.m )
par(mfrow=c(2,3))
plot( tel.var[,'50%'], mam.var[,'50%'], cex=0.5, xaxs='i', yaxs='i',
     xlab='', ylab='', col=rgb(0,0,0,0.2))
abline(h=log2(75), col='red', lty=2)
abline(v=log2(75), col='red', lty=2)
mtext('Teleost', side=1, line=2.5)
mtext('Mammal', side=2, line=2.5)
with(par(), mtext('A', at=usr[1], cex=mt.cex, line=1))
with(tel.mam.50, rect_image( t(log(1 + h)), rb, cb, col.f=v2col, col.args=list(cm=0.3),
                            xlab='Teleost', ylab='Mammal', lab.cex=1))
with(par(), mtext('B', at=usr[1], cex=mt.cex, line=1))
with(tel.mam.50, rect_image(scale(t(h)), rb, cb, col.f=v2col, col.args=list(cm=0.3),
                            xlab='Teleost', ylab='Mammal', lab.cex=1))
with(par(), mtext('C', at=usr[1], cex=mt.cex, line=1))
### Then repeat the above code (oh, so ugly this makes me feel)
cex.lab <- 1.5
mgp <- c(2.5, 1, 0)
b <- !is.na(mam.var[,'50%'])
plot(1:sum(b) / sum(b), sort(mam.var[b,'50%']), type='p', xlab='quantile', ylab='log2 intron length',
     cex.lab=cex.lab, mgp=mgp)
b <- !is.na(tel.var[,'50%'])
points(1:sum(b) / sum(b), sort(tel.var[b,'50%']), type='p', col='red' )
with(par(), mtext('D', at=usr[1], cex=mt.cex, line=1))
#
## let us have some values
i <- sum(tel.var[b,'50%'] <= 10) / sum(b)
with(par(), lines( c(i,i,usr[1], usr[2]), c(usr[3], 10, 10, 10) ))
legend('topleft', c('mammalia', 'teleostia'), pch=1, col=c('black', 'red'))
with(par(), text( usr[1], 10, sprintf("%.1f %% <= 10", 100*(i)), adj=c(-0.1, 1.5), col='red'))
with(par(), text( usr[2], 10, sprintf("%.1f %% > 10", 100*(1-i)), adj=c(1.1, 1.5), col='red' ))
#
##
tmp <- plot.pred.1(use.index=FALSE, xlab='log2 teleost median length', ylab='fraction above mammal quantile', plot.exp.r=FALSE, use.log=FALSE, cex.lab=cex.lab, mgp=mgp)
i <- sum(tmp$x[,1] >= 10)
i.y <- max( tmp$y[i,] )
with(par(), lines( c(10, 10, usr[2]), c(usr[3], i.y, i.y), col=tmp$cols[1]), )
with(par(), text(usr[2], i.y - 1:ncol(tmp$y) * strheight('A') * 1.5,
                 sprintf('%.1f%%', 100 * tmp$y[i,]), adj=c(1.2,1), col=tmp$cols) )
legend('topleft', legend=names(long.pred), text.col=tmp$col, lty=1, col=tmp$col,
       title='mammalian\nquantiles')
with(par(), mtext('E', at=usr[1], cex=mt.cex, line=1))
##
plot.pred.1(use.index=TRUE, xlab='teleost median quantile', ylab='log2 obs-expected ratio', plot.exp.r=TRUE, use.log=TRUE, cex.lab=cex.lab, mgp=mgp)
legend('topleft', legend=names(long.pred), text.col=tmp$col, lty=1, col=tmp$col,
       title='mammalian\nquantiles')
with(par(), mtext('F', at=usr[1], cex=mt.cex, line=1))
##
## and switch off the device
dev.off()


### Simple plot of the number of orthologous introns for each species.

orth.n <- apply(intron.orth.l, 2, function(x){ sum(!is.na(x)) })

cairo_pdf("orthologue_counts.pdf", width=a4.w * 0.8 * pdf.m, height=a4.h * 0.9 * pdf.m )
par(mar=c(5.1, 8.1, 4.1, 0.1))
layout(matrix(1:2, nrow=1), widths=c(0.85, 0.15))
sp <- names(sp.y)
y <- barplot(orth.n[sp], las=2, names.arg="", ## uc.1(sub("([^_]+)_.+", "\\1",  sp)),
        cex.names=0.5, col=sp.col.3[sp], xaxs='i', yaxs='i', horiz=T )
## mtext(uc.1(sub("([^_]+)_.+", "\\1",  sp)), side=2, at=y, cex=0.5, las=2)
mtext(uc.1(sub("_", " ",  sp)), side=2, at=y, cex=0.5, las=2)
par(mar=c(5.1, 0.1, 4.1, 1.1))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
legend('topleft', legend=names(class.col.3), fill=class.col.3)
dev.off()

quantile(orth.n, probs=seq(0, 1, 0.1))
##      0%     10%     20%     30%     40%     50%     60%     70%     80%     90% 
## 22349.0 47610.2 50025.6 51602.1 52757.8 54007.5 55386.8 56108.3 56720.4 57340.6 
##    100% 
## 63068.0 

orth.sp.rep <- apply(intron.orth.l, 1, function(x){ sum(!is.na(x)) })
quantile(orth.sp.rep, probs=seq(0, 1, 0.1))
## 0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
##  1  111  135  145  151  156  159  162  164  166  172 

## how many genes represented in the orthology?
length(unique(intron.orth.id[,'danio_rerio']))
## 5752 genes
