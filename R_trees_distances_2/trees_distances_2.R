## Redo the Sankoff maximum parsimony, but with much finer grained discretization.
## This is probably rather stupid, and will be massively computationally expensive
## both in terms of CPU and memory, but finding a better way would take longer.

## read in trees created from the intron orthology.

source("../R/functions.R")
source("~/R/general_functions.R")

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


## plotting dimensions
a4.w <- 8.27
a4.h <- 11.69
pdf.m <- 1.6
mt.cex <- 2

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


## distances by genetic distances (by alignment of exons)
## (Jukes-Cantor, Jukes-Cantor-gap, Kimura two factor)
ex.align.2.k2 <-  readRDS("../R_172_genomes/ex_align_2_k2.rds")

ex.align.2.k2.nj <- readRDS( "../R_172_genomes/ex_align_2_k2_nj.rds" )
require(ape)
plot(ex.align.2.k2.nj)

total.nodes <- length(unique(as.numeric(ex.align.2.k2.nj$edge)))
leaf.n <- length( with(ex.align.2.k2.nj, setdiff(edge[,2], edge[,1])) )

## named list containing;
## edge, edge.length, tip.label, Nnode
## the tip nodes are 1->172, and are present only in the second column of edge
## To draw this we can use a function edge.lines
source("/home/lmj/R/R_max_parsimony/R/functions.R")
tree.lines <- edge.lines( ex.align.2.k2.nj )

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

class.col.3 <- readRDS("../R_trees_distances/class_col_3.rds")

nodes.col <- rep('black', length(mammalia.nodes) )
nodes.col[ mammalia.nodes ] <- rgb(0, 0, 0.8)
nodes.col[ eutheria.nodes ] <- rgb(0.6, 0, 0.6)
nodes.col[ teleost.nodes  ] <- rgb(0.8, 0, 0)
nodes.col[ sauria.nodes ] <- rgb(0, 0.8, 0)
nodes.col[ aves.nodes ] <- rgb(0, 0.6, 0.6)

sp.col.3 <- nodes.col[ 1:length(ex.align.2.k2.nj$tip.label) ]
names(sp.col.3) <- ex.align.2.k2.nj$tip.label

## We can in a similar manner make a class tree:
sp.class.3 <- cbind('teleostei'=teleost.nodes, 'mammalia'=mammalia.nodes,
                    'eutheria'=eutheria.nodes, 'aves'=aves.nodes,
                    'sauria'=sauria.nodes)[ 1:length(ex.align.2.k2.nj$tip.label), ]
rownames(sp.class.3) <- ex.align.2.k2.nj$tip.label
saveRDS( sp.class.3, "sp_class_3.rds" )

## it is then possible to draw:
plot.new()
with(tree.lines, plot.window(xlim=c(0, 1.5 * max(x[,1:2])),
                             ylim=c(0, 1.02 * max(y)), xaxs='i', yaxs='i'))
with(tree.lines, segments(x[,1], y, x[,2], y, col=nodes.col[x[,'p'] ] ))
with(tree.lines, segments(v[,1], v[,2], v[,1], v[,3], col=nodes.col[v[,'node']] ))
with(tree.lines, text( x[,2] + strwidth("9", cex=0.4), y, labels=nodes, adj=c(0,0.5), cex=0.4))
with(tree.lines, {
    b <- !(nodes %in% ex.align.2.k2.nj$edge[,1])
    offset <- strwidth("9999", cex=0.4) * 1.25
    text( x[b,2] + offset, y[b], uc.1(sp.db2sp(ex.align.2.k2.nj$tip.label[ nodes[b] ])),
         cex=0.5, adj=c(0,0.5) )
    })
legend('bottomleft', legend=names(class.col.3), text.col=class.col.3, box.lty=0,
       inset=c(0,0.05))

int.s.l <- readRDS("../R_172_genomes/int_s_l.rds")

## This far I have just been repeating myself. Now it gets more interesting
## encode.dist takes a conversion function: We will use
## (encode.dist is part of the sankoff functions)

## The sankoff function I have taken takes the states as character vectors.
## This means that I have a maximum of 256 possible states that I can use.
## Instead I can use the ace function from the ape package.
## That will be horribly expensive as well, but it can be use for continuous characters;
## We will need to convert handle missing states; Unfortunately the ace function is unable
## to handle NAs for discrete states so we will have to replace any such states. The
## only reasonable way in which I can work that out would be to make an unknown state
## from closest neighbour. This means going up the tree for all unknown states and assigning
## from the parent. I would also need to remember how to do that. Which is a bugger.

## I have modified the sankoff encode helper function.
## for this reason I should run it as before first to make sure that it does the
## same thing.


dyn.load("~/R/R_max_parsimony/src/max_parsimony.so")
source("~/R/R_max_parsimony/R/functions.R")

## we then need to make a substitution matrix
al.size <- as.integer( 1 + max(int.s.l, na.rm=TRUE) - min(int.s.l, na.rm=TRUE))
al.offset <- 64L  ## starts at @
sub.matrix <- make.sub.matrix( al.size )
## first row and column indicate missing values, and should have maximum penalties
missing.penalty <- as.integer(max(sub.matrix) + 1)
sub.matrix[1,] <- missing.penalty
sub.matrix[,1] <- missing.penalty

int.s.enc <- encode.dist( int.s.l, offset=al.offset )

## then we can infer the ancestral states.
nj.edge <- matrix( as.integer(ex.align.2.k2.nj$edge), nrow=nrow(ex.align.2.k2.nj$edge))

int.s.inf <- sankoff( nj.edge, c(total.nodes, leaf.n), sub.matrix, c(al.offset, al.size),
                      int.s.enc )
### traverse the leaf nodes to get the distance. Use both the tree (for
### phylogenetic differences) and a summarised change in intron size
### We have a total of leaf.n leaf nodes..

## this differs from the version in ../R_172_genomes/ in that we only
## use reasonable intron sizes
leaf.int.d <- lapply( 1:leaf.n, function(i){
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
names(leaf.int.d) <- colnames( int.s.l )


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

## we can put labels at places that are of interest
label.i <- which( names(leaf.int.d) %in%
                  c('danio rerio', 'tetraodon nigroviridis', 'betta splendens',
                    'astyanax mexicanus', 'erpetoichthys calabaricus', 'latimeria chalumnae',
                    'sphenodon punctatus', 'pygocentrus nattereri','eptatretus burgeri',
                    'takifugu rubripes', 'xenopus tropicalis', 'callorhinchus milii',
                    'lepisosteus oculatus', 'gasterosteus aculeatus', 'parambassis ranga',
                    'electrophorus electricus', 'denticeps clupeoides'))


xlim <- range( unlist(lapply( leaf.int.d, function(x){ x[,'pg'] })) )
xlim[2] <- 1.2 * xlim[2]
xlim[1] <- 0
ylim <- range( unlist(lapply( leaf.int.d, function(x){ x[,'int.d'] })) )

##cairo_pdf("Intron_size_evolution.pdf", width=a4.w*0.5*pdf.m, height=a4.h*0.75*pdf.m)
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
##dev.off()


### That seems to work just as before. Now we want to modify this to use a much larger
### alphabet.

## this will allow us to count differences of 7% in intron length.
## that is still quite a large indel, and so we really would want to have
## a continuous function. But lets try that in any case and see how it goes
al.2.size <- as.integer( 1 + diff( range(int.s.l * 10, na.rm=TRUE )) )
al.2.offset <- 33L  ## we have to have some sort of an offset.
sub.matrix.2 <- make.sub.matrix( al.2.size )
## first row and column indicate missing values, and should have maximum penalties as before
missing.penalty <- as.integer(max(sub.matrix.2) + 1)
sub.matrix.2[1,] <- missing.penalty
sub.matrix.2[,1] <- missing.penalty

int.s.2.enc <- encode.dist( int.s.l * 10, offset=al.2.offset )

int.s.inf.2 <- sankoff( nj.edge, c(total.nodes, leaf.n), sub.matrix.2, c(al.2.offset, al.2.size),
                       int.s.2.enc, check.states=FALSE )

## this differs from the version in ../R_172_genomes/ in that we only
## use reasonable intron sizes
leaf.int.2.d <- lapply( 1:leaf.n, function(i){
    pg.dist <- c()
    int.d <- c()
    node.i <- c()
    ascend.tree <- function(j){
        b <- ex.align.2.k2.nj$edge[,2] == j
        if(sum(b) == 1){
            pg.dist <<- c(ex.align.2.k2.nj$edge.length[b], pg.dist)
            p.j <- ex.align.2.k2.nj$edge[b,1]
            b2 <- (int.s.inf.2[[j]]$state[,1] >= 60 & int.s.inf.2[[p.j]]$state[,1] >= 60)
            int.d <<- c(mean(int.s.inf.2[[j]]$state[b2,2]), int.d)
            node.i <<- c(j, node.i)
            ascend.tree( p.j )
        }
    }
    ascend.tree(i)
    cbind('pg'=cumsum(pg.dist),
          'int.d'=cumsum(int.d),
          'node'=node.i)
})
names(leaf.int.2.d) <- colnames( int.s.l )

## we can put labels at places that are of interest
label.i <- which( names(leaf.int.d) %in%
                  c('danio rerio', 'tetraodon nigroviridis', 'betta splendens',
                    'astyanax mexicanus', 'erpetoichthys calabaricus', 'latimeria chalumnae',
                    'sphenodon punctatus', 'pygocentrus nattereri','eptatretus burgeri',
                    'takifugu rubripes', 'xenopus tropicalis', 'callorhinchus milii',
                    'lepisosteus oculatus', 'gasterosteus aculeatus', 'parambassis ranga',
                    'electrophorus electricus', 'denticeps clupeoides',
                    'myotis lucifugus', 'pteropus vampyrus', 'coturnix japonicus',
                    'apteryx rowi', 'coturnix japonica',
                    'meleagris gallopavo'))


xlim <- range( unlist(lapply( leaf.int.2.d, function(x){ x[,'pg'] })) )
xlim[2] <- 1.2 * xlim[2]
xlim[1] <- 0
ylim <- range( unlist(lapply( leaf.int.2.d, function(x){ x[,'int.d'] })) )

cairo_pdf("Intron_size_evolution.pdf", width=a4.w*0.5*pdf.m, height=a4.h*0.75*pdf.m)
plot.new()
plot.window(xlim=xlim, ylim=ylim, xaxs='i')
abline(h=0, col='grey')
ends <- sapply(leaf.int.2.d, draw.delta.lines, lwd=2)
axis(2)
with( leaf.int.2.d, {
    o <- order( ends[2,label.i] )
    x <- (ends[1,label.i])[o]
    y <- sort(ends[2,label.i])
    y.d <- diff(y)
    y.d.i <- which(y.d < 1.1 * strheight( 'A', cex=0.8 ) )
    count <- 0
    while(length(y.d.i) && count < 20){
        y[ y.d.i ] <- y[y.d.i] - strheight('A') * 0.25
        y[ y.d.i + 1] <- y[y.d.i + 1] + strheight('A') * 0.25
        y.d <- diff(y)
        y.d.i <- which(y.d < 1.1 * strheight( 'A', cex=0.8 ) )
        count <- count + 1
    }
    text( x, y, uc.1(colnames(ends)[label.i])[o], pos=4,
     col=nodes.col[label.i[o]], cex=0.8) })
legend('bottomleft', legend=names(class.col.3), text.col=class.col.3, box.lty=0,
       inset=c(0,0.05))
dev.off()

identify( t(ends), labels=colnames(ends))


### We can look at the changes at the leaf nodes to see if there is evidence for
### for continuing changes. That would be better done by individual distributions
### for the underlying numbers, but we can first have a look at the mean changes:


leaf.d <- t(sapply( int.s.inf.2[1:length(leaf.int.2.d)], function(x){
    b <- x$state[,1] >= 60
    c( mean( x$state[b,2]), median( x$state[b,2] ))
}))
colnames(leaf.d) <- c('mean', 'median')
rownames(leaf.d) <- names(leaf.int.d)

## let us do a manual bar plot so that we can add some information more easily
sp.lab <- strsplit(rownames(leaf.d), ' ')
sp.lab <- sapply( sp.lab, function(x){ paste( toupper(substr(x[1], 1, 1)), ". ", x[2], sep="") })

## we need int.s.b from 172 genomes
int.s.b <- readRDS("../R_trees_distances/int_s_b.rds")

cairo_pdf( 'leaf_mean_intron_changes.pdf', width=a4.w * pdf.m, height=a4.h * pdf.m )
par(omi=c(0.5, 0.5, 0.5, 0.5)) 
par(mar=c(5.1, 4.1, 4.1, 0.1))
layout(matrix(1:2, nrow=1), widths=c(15, 90))
plot.new()
plot.window( xlim=0:1, ylim=c(0,nrow(leaf.d)), xaxs='i', yaxs='i')
text(1, 1:nrow(leaf.d)-0.5, sp.lab, cex=0.5, adj=c(1,0.5), col=ifelse(int.s.b, 'black', 'red'))
par(mar=c(5.1, 0.1, 4.1, 2.1))
plot.new()
plot.window( xlim=range(leaf.d[,'mean']), ylim=c(0,nrow(leaf.d)), xaxs='i', yaxs='i')
y <- 1:nrow(leaf.d)
with(par(), rect(usr[1], y-1, usr[2], y, col=c(rgb(1, 1, 1), rgb(0.9, 0.9, 0.9)), border=NA ))
with(par(), rect(usr[1], y[!int.s.b]-1, usr[2], y[!int.s.b], col=rgb(1, 0.7, 0.7, 0.3), border=NA))
rect( leaf.d[,'mean'], y-1, 0, y, col=nodes.col[ 1:length(leaf.d) ] )
axis(1)
legend('bottomright', legend=names(class.col.3), text.col=class.col.3, box.lty=0,
       inset=c(0,0))
dev.off()



