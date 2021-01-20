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

sp.db2sp.2 <- function(x){
    x <- sub( "([^_]+)_([^_]+)_?.*", "\\1 \\2", x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    xl  <- strsplit(x, " ")
    x  <- sapply( xl, function(y){
        paste( substr(y[1], 1, 1), ". ", paste(y[-1], collapse=" "), sep="" )})
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

half.w <- a4.w * 85 / 210
full.w <- a4.w * 170 / 210
max.h <- a4.w * 225 / 210

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

cairo_pdf("Intron_size_evolution.pdf", width=half.w * pdf.m, height=max.h*0.75*pdf.m)
par(mar=c(1.1, 3.1, 1.1, 1.1))
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
       inset=c(0.01,0.05))
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


#### Have the same sets of introns been minimised independently ?
#### ex.align.2.k2.nj, we have the following nodes:

## overlap between introns that remain long in the two descendants
## long is defined as >= l.min
## loss is defined as becoming <= s.max
long.enrichment  <- function(s.max, l.min, ancestor, d1, d2, i.b=NULL, min.dec=1){
    if(is.null(i.b))
        i.b  <- int.s.inf.2[[d1]]$state[,1] > 60 & int.s.inf.2[[d2]]$state[,1] > 60
    l.b  <- int.s.inf.2[[ancestor]]$state[i.b,1] >= l.min
    l.b1  <- l.b & int.s.inf.2[[d1]]$state[i.b,1] > s.max
    l.b2  <- l.b & int.s.inf.2[[d2]]$state[i.b,1] > s.max
    l.b1.r  <- l.b & int.s.inf.2[[d1]]$state[i.b,1] <= s.max
    l.b2.r  <- l.b & int.s.inf.2[[d2]]$state[i.b,1] <= s.max
    ## which decrease or grow more than min.dec
    d.b1  <- l.b & (int.s.inf.2[[ancestor]]$state[i.b,1] - int.s.inf.2[[d1]]$state[i.b,1]) > min.dec
    d.b2  <- l.b & (int.s.inf.2[[ancestor]]$state[i.b,1] - int.s.inf.2[[d2]]$state[i.b,1]) > min.dec
    g.b1  <- l.b & (int.s.inf.2[[d1]]$state[i.b,1] - int.s.inf.2[[ancestor]]$state[i.b,1]) > min.dec
    g.b2  <- l.b & (int.s.inf.2[[d2]]$state[i.b,1] - int.s.inf.2[[ancestor]]$state[i.b,1]) > min.dec
    ## and we check the enrichment of group 1 in 2
    p  <- phyper( sum(l.b1 & l.b2)-1, sum(l.b1), sum(l.b) - sum(l.b1), sum(l.b2), lower.tail=FALSE )
    p.r  <- phyper( sum(l.b1.r & l.b2.r)-1, sum(l.b1.r), sum(l.b) - sum(l.b1.r), sum(l.b2.r), lower.tail=FALSE )
    p.d  <- phyper( sum(d.b1 & d.b2)-1, sum(d.b1), sum(l.b) - sum(d.b1), sum(d.b2), lower.tail=FALSE )
    p.g  <- phyper( sum(g.b1 & g.b2)-1, sum(g.b1), sum(l.b) - sum(g.b1), sum(g.b2), lower.tail=FALSE )
    c('ll'=sum(l.b1 & l.b2), 'l1'=sum(l.b1), 'l2'=sum(l.b2), 'l1.2'=sum(l.b1 | l.b2), 'n'=sum(l.b), 'p'=p,
      'o.e'=sum(l.b1 & l.b2) / (sum(l.b2) * sum(l.b1) / sum(l.b)),
      'll.r'=sum(l.b1.r & l.b2.r), 'l1.r'=sum(l.b1.r), 'l2.r'=sum(l.b2.r), 'l1.2.r'=sum(l.b1.r | l.b2.r), 'p.r'=p.r,
      'o.e.r'=sum(l.b1.r & l.b2.r) / (sum(l.b2.r) * sum(l.b1.r) / sum(l.b)),
      'p.d'=p.d, 'g.d'=p.g
      )
}

evo.hist  <- function(ancestor, d1, i.b=NULL, ...){
    if(is.null(i.b))
        i.b  <- int.s.inf.2[[d1]]$state[,1] > 60
    hist( int.s.inf.2[[d1]]$state[i.b,1] - int.s.inf.2[[ancestor]]$state[i.b,1], ...)
}

### We can define statistics for all the descendants of node 288
### apart from 1 and 2, what proportions have increased or decreased in size.
### Nodes 1 and 2 are T. rubripes and T. nigroviridis and form a separate branch
### in our tree. Note that it looks like Mola_mola should be on the same
### branch as T. rubripes and T. nigroviridis as it is a Tetraodontiformes.

desc.288  <- which( nodes.descend( ex.align.2.k2.nj, 288 )[ 1:length(ex.align.2.k2.nj$tip.label) ] )

t.fugu.288  <- sapply( desc.288[c(-1,-2)], function(x){ long.enrichment( 80, 80, ancestor=288, d1=2, x ) })
colnames(t.fugu.288)  <- ex.align.2.k2.nj$tip.label[ desc.288[-(1:2)] ]

## It would be nice to give this a group name. We can try to do this by having a look at the NCBI taxonomy:
## Takifugu:      Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Eupercaria; Tetraodontiformes; Tetraodontoidei; Tetradontoidea; Tetraodontidae
## Monopterus     Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Anabantaria; Synbranchiformes; Synbranchoidei; Synbranchidae;
## Mastacembelus  Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Anabantaria; Synbranchiformes; Mastacembeloidei; Mastacembelidae
## Anabas:        Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Anabantaria; Anabantiformes; Anabantoidei; Anabantidae
## Betta:         Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Anabantaria; Anabantiformes; Anabantoidei; Osphronemidae; Macropodinae;
## Parambassis    Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Ovalentaria incertae sedis; Ambassidae
## Oreochromis    Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Cichlomorphae; Cichliformes; Cichlidae; African cichlids; Pseudocrenilabrinae; Oreochromini
## Neolamprologus Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Cichlomorphae; Cichliformes; Cichlidae; African cichlids; Pseudocrenilabrinae; Lamprologini
## Pundamilia     Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Cichlomorphae; Cichliformes; Cichlidae; African cichlids; Pseudocrenilabrinae; Haplochromini
## Haplochromis   Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Cichlomorphae; Cichliformes; Cichlidae; African cichlids; Pseudocrenilabrinae; Haplochromini
## Maylandia      Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Cichlomorphae; Cichliformes; Cichlidae; African cichlids; Pseudocrenilabrinae; Haplochromini
## Astatotilapia  Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Cichlomorphae; Cichliformes; Cichlidae; African cichlids; Pseudocrenilabrinae; Haplochromini
## Stegastes      Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Ovalentaria incertae sedis; Pomacentridae
## Acanthochromis Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Ovalentaria incertae sedis; Pomacentridae
## Amphiprion     Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Ovalentaria incertae sedis; Pomacentridae
## Labrus         Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Eupercaria; Labriformes; Labridae
## Mola           Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Eupercaria; Tetraodontiformes; Moloidei; Molidae  #####****#####
## Scopthalmus    Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Carangaria; Pleuronectiformes; Pleuronectoidei; Scophthalmidae
## Seriola        Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Carangaria; Carangiformes; Carangidae
## Lates          Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Carangaria; Carangaria incertae sedis; Centropomidae
## Gasterosteus   Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Eupercaria; Perciformes; Cottioidei; Gasterosteales; Gasterosteidae
## Larimichthys   Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Eupercaria; Eupercaria incertae sedis; Sciaenidae
## Cottoperca     Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Eupercaria; Perciformes; Notothenioidei; Bovichtidae
## Kryptolebias   Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Atherinomorphae; Cyprinodontiformes; Aplocheiloidei; Rivulidae
## Cyprinodon     Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Atherinomorphae; Cyprinodontiformes; Cyprinodontoidei; Cyprinodontidae; Cyprinodontinae; Cyprinodontini
## Fundulus       Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Atherinomorphae; Cyprinodontiformes; Cyprinodontoidei; Fundulidae
## Poecilia:      Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Atherinomorphae; Cyprinodontiformes; Cyprinodontoidei; Poeciliidae; Poeciliinae
## Gambusia       Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Atherinomorphae; Cyprinodontiformes; Cyprinodontoidei; Poeciliidae; Poeciliinae
## Xiphophorus    Acanthomorphata; Euacanthomorphacea; Percomorphaceae; Ovalentaria; Atherinomorphae; Cyprinodontiformes; Cyprinodontoidei; Poeciliidae; Poeciliinae

## note that some species that are defined as being Ovalentaria are located outside of this
## set, suggesting either errors in the tree (not that unlikely) or issues with the
## standard taxonomy (less likely).

## in order to draw:
tree.err  <- c("mola_mola")
tree.eupercaria  <- c("cottoperca_gobio", "larimichthys_crocea", "gasterosteus_aculeatus",
                      "labrus_bergylta")
tree.f  <- setdiff( colnames(t.fugu.288), c(tree.err, tree.eupercaria))
tree.fb  <- colnames( t.fugu.288 ) %in% tree.f 


cairo_pdf("selective_retention_1.pdf", width=half.w * pdf.m, height=half.w * pdf.m * 1.25 )
par(mfrow=c(3,2))
par(mar=c(5.1, 4.1, 0.4, 1.1))
par(oma=c(0, 0, 2.1, 1.0))
## What proportion are retained as long?
t.fugu.288.col  <- ifelse( tree.fb, 'black', 'gray' )
lab.i  <- match( c("betta_splendens", "gasterosteus_aculeatus", "parambassis_ranga", "anabas_testudineus",
                   "scophthalmus_maximus", "mastacembelus_armatus", "mola_mola"), colnames(t.fugu.288))
lab.lab  <- paste(1:length(lab.i), sp.db2sp.2( colnames(t.fugu.288)[ lab.i ] ) )
plot( genome.sizes[ colnames(t.fugu.288) ], t.fugu.288['l2',] / t.fugu.288['n',], col=t.fugu.288.col,
     xlab='genome size', ylab='proportion retained', ylim=c(0.65, 1), pch=19 )
points( genome.sizes[ colnames(t.fugu.288) ], t.fugu.288['l1',] / t.fugu.288['n',], col=t.fugu.288.col, pch=1 )
legend( 'topleft', legend=c('T. rubripes', 'other'), pch=c(1,19), cex=0.8, box.col='black' )
with(par(), mtext("A", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(t.fugu.288[,lab.i]) ], t.fugu.288['l2',lab.i] / t.fugu.288['n',lab.i], pos=1 )
plot( t.fugu.288['l1.2',] / t.fugu.288['n',], t.fugu.288['o.e',], xlab='proportion retained (union)', ylab='observed/expected',
     col=t.fugu.288.col, pch=19)
legend('topright', legend=c('Eupercaria', 'other'), pch=19, col=c('gray', 'black'), cex=0.8)
text( t.fugu.288['l1.2',lab.i] / t.fugu.288['n',lab.i], t.fugu.288['o.e',lab.i], pos=c(1, 1, 1, 1, 3, 1, 3) )
with(par(), mtext("B", cex=mt.cex, at=usr[1], adj=1, line=0.2))
## and the enrichment plots
plot( genome.sizes[ colnames(t.fugu.288) ], -log10(1e-320 + t.fugu.288['p',]), xlab='genome size', ylab='-log10 p', col=t.fugu.288.col, pch=19 )
text( genome.sizes[ colnames(t.fugu.288[,lab.i]) ], -log10(1e-320 + t.fugu.288['p',lab.i]), pos=1 )
with(par(), mtext("C", cex=mt.cex, at=usr[1], adj=1, line=0.2))
plot( ex.align.2.k2[ 'takifugu_rubripes', colnames(t.fugu.288) ], -log10(1e-320 + t.fugu.288['p',]),
     xlab='Kimura distance', ylab='-log10 p', col=t.fugu.288.col, pch=19 )
with(par(), mtext("D", cex=mt.cex, at=usr[1],, adj=1, line=0.2))
text( ex.align.2.k2[ 'takifugu_rubripes', colnames(t.fugu.288[,lab.i]) ], -log10(1e-320 + t.fugu.288['p',lab.i]), pos=1 )
plot( genome.sizes[ colnames(t.fugu.288) ], t.fugu.288['o.e',], xlab='genome size', ylab='observed / expected', col=t.fugu.288.col, pch=19 )
with(par(), mtext("E", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(t.fugu.288[,lab.i]) ], t.fugu.288['o.e',lab.i], pos=1 )
lab.w  <- max( strwidth( lab.lab, cex=0.8 ) )
with(par(), text(usr[2] - lab.w * 1.1, 1.15, paste(lab.lab, collapse='\n'), adj=c(0, 1), cex=0.8) )
plot( ex.align.2.k2[ 'takifugu_rubripes', colnames(t.fugu.288) ], t.fugu.288['o.e',], xlab='Kimura distance',
     ylab='', col=t.fugu.288.col, pch=19 )
with(par(), mtext("F", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( ex.align.2.k2[ 'takifugu_rubripes', colnames(t.fugu.288[,lab.i]) ], t.fugu.288['o.e',lab.i], pos=1 )
dev.off()


### let us have a look at the MRCA of Danio rerio and others in the same
### manner. This is more complicated, because many introns in D. rerio
### will be long, so the commonality will be less.

desc.266  <- which( nodes.descend( ex.align.2.k2.nj, 266 )[ 1:length(ex.align.2.k2.nj$tip.label) ] )
## that's basically all the teleosts.
## we want to remove those that share a clade with Danio rerio.
desc.267  <- which( nodes.descend( ex.align.2.k2.nj, 267 )[ 1:length(ex.align.2.k2.nj$tip.label) ] )

desc.266.dr  <- setdiff( desc.266, desc.267 )

## we can then ask the question of any node in desc.267 and desc.267.dr
d.rerio.266  <- sapply( desc.266.dr, function(x){ long.enrichment(80, 80, ancestor=266, d1=78, x)})
colnames(d.rerio.266)  <- ex.align.2.k2.nj$tip.label[ desc.266.dr ] 

sp  <- colnames(d.rerio.266)

cairo_pdf("selective_retention_2.pdf", width=half.w * pdf.m, height=half.w * pdf.m * 1.25 )
par(mfrow=c(3,2))
par(mar=c(5.1, 4.1, 0.4, 1.1))
par(oma=c(0, 0, 2.1, 1.0))
## What proportion are retained as long?
d.rerio.266.col  <- 'black'
lab.i  <- match( c("tetraodon_nigroviridis", "takifugu_rubripes", "betta_splendens", "gasterosteus_aculeatus", "parambassis_ranga", "anabas_testudineus",
                   "scophthalmus_maximus", "mastacembelus_armatus"), colnames(d.rerio.266))
lab.lab  <- paste(1:length(lab.i), sp.db2sp.2( colnames(d.rerio.266)[ lab.i ] ) )
plot( genome.sizes[ colnames(d.rerio.266) ], d.rerio.266['l2',] / d.rerio.266['n',], col=d.rerio.266.col,
     xlab='genome size', ylab='proportion retained', ylim=c(0.3, 0.85), pch=19 )
points( genome.sizes[ colnames(d.rerio.266) ], d.rerio.266['l1',] / d.rerio.266['n',], col=d.rerio.266.col, pch=1 )
legend( 'bottomright', legend=c('D. rerio', 'other'), pch=c(1,19), cex=0.8, box.col='black' )
with(par(), mtext("A", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(d.rerio.266[,lab.i]) ], d.rerio.266['l2',lab.i] / d.rerio.266['n',lab.i], pos=1 )
plot( d.rerio.266['l1.2',] / d.rerio.266['n',], d.rerio.266['o.e',], xlab='proportion retained (union)', ylab='observed/expected',
     col=d.rerio.266.col, pch=19)
text( d.rerio.266['l1.2',lab.i] / d.rerio.266['n',lab.i], d.rerio.266['o.e',lab.i], pos=c(1, 3, 1, 1, 3, 3, 3, 1) )
lab.w  <- max( strwidth( lab.lab, cex=0.8 ) )
with(par(), text(usr[1] + lab.w * 0.1, usr[3] + lab.w * 0.1, paste(lab.lab, collapse='\n'), adj=c(0, 0), cex=0.8) )
with(par(), mtext("B", cex=mt.cex, at=usr[1], adj=1, line=0.2))
## and the enrichment plots
plot( genome.sizes[ colnames(d.rerio.266) ], -log10(1e-320 + d.rerio.266['p',]), xlab='genome size', ylab='-log10 p', col=d.rerio.266.col, pch=19 )
text( genome.sizes[ colnames(d.rerio.266[,lab.i]) ], -log10(1e-320 + d.rerio.266['p',lab.i]), pos=1 )
with(par(), mtext("C", cex=mt.cex, at=usr[1], adj=1, line=0.2))
plot( ex.align.2.k2[ 'danio_rerio', colnames(d.rerio.266) ], -log10(1e-320 + d.rerio.266['p',]),
     xlab='Kimura distance', ylab='-log10 p', col=d.rerio.266.col, pch=19 )
with(par(), mtext("D", cex=mt.cex, at=usr[1],, adj=1, line=0.2))
text( ex.align.2.k2[ 'danio_rerio', colnames(d.rerio.266[,lab.i]) ], -log10(1e-320 + d.rerio.266['p',lab.i]), pos=1 )
plot( genome.sizes[ colnames(d.rerio.266) ], d.rerio.266['o.e',], xlab='genome size', ylab='observed / expected', col=d.rerio.266.col, pch=19 )
with(par(), mtext("E", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(d.rerio.266[,lab.i]) ], d.rerio.266['o.e',lab.i], pos=1 )
plot( ex.align.2.k2[ 'danio_rerio', colnames(d.rerio.266) ], d.rerio.266['o.e',], xlab='Kimura distance',
     ylab='', col=d.rerio.266.col, pch=19 )
with(par(), mtext("F", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( ex.align.2.k2[ 'danio_rerio', colnames(d.rerio.266[,lab.i]) ], d.rerio.266['o.e',lab.i], pos=1 )
dev.off()

## Then do the same for a couple of species with smaller genomes
e.elec.266  <- sapply( desc.266.dr, function(x){ long.enrichment(80, 80, ancestor=266, d1=9, x)})
colnames(e.elec.266)  <- ex.align.2.k2.nj$tip.label[ desc.266.dr ] 

sp  <- colnames(e.elec.266)

cairo_pdf("selective_retention_3.pdf", width=half.w * pdf.m, height=half.w * pdf.m * 1.25 )
par(mfrow=c(3,2))
par(mar=c(5.1, 4.1, 0.4, 1.1))
par(oma=c(0, 0, 2.1, 1.0))
## What proportion are retained as long?
e.elec.266.col  <- 'black'
lab.i  <- match( c("tetraodon_nigroviridis", "takifugu_rubripes", "betta_splendens", "gasterosteus_aculeatus", "parambassis_ranga", "anabas_testudineus",
                   "scophthalmus_maximus", "mastacembelus_armatus"), colnames(e.elec.266))
lab.lab  <- paste(1:length(lab.i), sp.db2sp.2( colnames(e.elec.266)[ lab.i ] ) )
plot( genome.sizes[ colnames(e.elec.266) ], e.elec.266['l2',] / e.elec.266['n',], col=e.elec.266.col,
     xlab='genome size', ylab='proportion retained', ylim=c(0.3, 0.75), pch=19 )
points( genome.sizes[ colnames(e.elec.266) ], e.elec.266['l1',] / e.elec.266['n',], col=e.elec.266.col, pch=1 )
legend( 'bottomright', legend=c('E. electrophorus', 'other'), pch=c(1,19), cex=0.8, box.col='black' )
with(par(), mtext("A", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(e.elec.266[,lab.i]) ], e.elec.266['l2',lab.i] / e.elec.266['n',lab.i], pos=1 )
plot( e.elec.266['l1.2',] / e.elec.266['n',], e.elec.266['o.e',], xlab='proportion retained (union)', ylab='observed/expected',
     col=e.elec.266.col, pch=19)
text( e.elec.266['l1.2',lab.i] / e.elec.266['n',lab.i], e.elec.266['o.e',lab.i], pos=c(1, 3, 1, 1, 3, 3, 3, 1) )
lab.w  <- max( strwidth( lab.lab, cex=0.8 ) )
with(par(), text(usr[1] + lab.w * 0.1, usr[3] + lab.w * 0.1, paste(lab.lab, collapse='\n'), adj=c(0, 0), cex=0.8) )
with(par(), mtext("B", cex=mt.cex, at=usr[1], adj=1, line=0.2))
## and the enrichment plots
plot( genome.sizes[ colnames(e.elec.266) ], -log10(1e-320 + e.elec.266['p',]), xlab='genome size', ylab='-log10 p', col=e.elec.266.col, pch=19 )
text( genome.sizes[ colnames(e.elec.266[,lab.i]) ], -log10(1e-320 + e.elec.266['p',lab.i]), pos=1 )
with(par(), mtext("C", cex=mt.cex, at=usr[1], adj=1, line=0.2))
plot( ex.align.2.k2[ 'electrophorus_electricus', colnames(e.elec.266) ], -log10(1e-320 + e.elec.266['p',]),
     xlab='Kimura distance', ylab='-log10 p', col=e.elec.266.col, pch=19 )
with(par(), mtext("D", cex=mt.cex, at=usr[1],, adj=1, line=0.2))
text( ex.align.2.k2[ 'electrophorus_electricus', colnames(e.elec.266[,lab.i]) ], -log10(1e-320 + e.elec.266['p',lab.i]), pos=1 )
plot( genome.sizes[ colnames(e.elec.266) ], e.elec.266['o.e',], xlab='genome size', ylab='observed / expected', col=e.elec.266.col, pch=19 )
with(par(), mtext("E", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(e.elec.266[,lab.i]) ], e.elec.266['o.e',lab.i], pos=1 )
plot( ex.align.2.k2[ 'electrophorus_electricus', colnames(e.elec.266) ], e.elec.266['o.e',], xlab='Kimura distance',
     ylab='', col=e.elec.266.col, pch=19 )
with(par(), mtext("F", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( ex.align.2.k2[ 'electrophorus_electricus', colnames(e.elec.266[,lab.i]) ], e.elec.266['o.e',lab.i], pos=1 )
dev.off()


d.clup.266  <- sapply( desc.266.dr, function(x){ long.enrichment(80, 80, ancestor=266, d1=10, x)})
colnames(d.clup.266)  <- ex.align.2.k2.nj$tip.label[ desc.266.dr ] 

sp  <- colnames(d.clup.266)

cairo_pdf("selective_retention_4.pdf", width=half.w * pdf.m, height=half.w * pdf.m * 1.25 )
par(mfrow=c(3,2))
par(mar=c(5.1, 4.1, 0.4, 1.1))
par(oma=c(0, 0, 2.1, 1.0))
## What proportion are retained as long?
d.clup.266.col  <- 'black'
lab.i  <- match( c("tetraodon_nigroviridis", "takifugu_rubripes", "betta_splendens", "gasterosteus_aculeatus", "parambassis_ranga", "anabas_testudineus",
                   "scophthalmus_maximus", "mastacembelus_armatus"), colnames(d.clup.266))
lab.lab  <- paste(1:length(lab.i), sp.db2sp.2( colnames(d.clup.266)[ lab.i ] ) )
plot( genome.sizes[ colnames(d.clup.266) ], d.clup.266['l2',] / d.clup.266['n',], col=d.clup.266.col,
     xlab='genome size', ylab='proportion retained', ylim=c(0.3, 0.75), pch=19 )
points( genome.sizes[ colnames(d.clup.266) ], d.clup.266['l1',] / d.clup.266['n',], col=d.clup.266.col, pch=1 )
legend( 'bottomright', legend=c('D. clupeoides', 'other'), pch=c(1,19), cex=0.8, box.col='black' )
with(par(), mtext("A", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(d.clup.266[,lab.i]) ], d.clup.266['l2',lab.i] / d.clup.266['n',lab.i], pos=1 )
plot( d.clup.266['l1.2',] / d.clup.266['n',], d.clup.266['o.e',], xlab='proportion retained (union)', ylab='observed/expected',
     col=d.clup.266.col, pch=19)
text( d.clup.266['l1.2',lab.i] / d.clup.266['n',lab.i], d.clup.266['o.e',lab.i], pos=c(1, 3, 1, 1, 3, 3, 3, 1) )
lab.w  <- max( strwidth( lab.lab, cex=0.8 ) )
with(par(), text(usr[1] + lab.w * 0.1, usr[3] + lab.w * 0.1, paste(lab.lab, collapse='\n'), adj=c(0, 0), cex=0.8) )
with(par(), mtext("B", cex=mt.cex, at=usr[1], adj=1, line=0.2))
## and the enrichment plots
plot( genome.sizes[ colnames(d.clup.266) ], -log10(1e-320 + d.clup.266['p',]), xlab='genome size', ylab='-log10 p', col=d.clup.266.col, pch=19 )
text( genome.sizes[ colnames(d.clup.266[,lab.i]) ], -log10(1e-320 + d.clup.266['p',lab.i]), pos=1 )
with(par(), mtext("C", cex=mt.cex, at=usr[1], adj=1, line=0.2))
plot( ex.align.2.k2[ 'denticeps_clupeoides', colnames(d.clup.266) ], -log10(1e-320 + d.clup.266['p',]),
     xlab='Kimura distance', ylab='-log10 p', col=d.clup.266.col, pch=19 )
with(par(), mtext("D", cex=mt.cex, at=usr[1],, adj=1, line=0.2))
text( ex.align.2.k2[ 'denticeps_clupeoides', colnames(d.clup.266[,lab.i]) ], -log10(1e-320 + d.clup.266['p',lab.i]), pos=1 )
plot( genome.sizes[ colnames(d.clup.266) ], d.clup.266['o.e',], xlab='genome size', ylab='observed / expected', col=d.clup.266.col, pch=19 )
with(par(), mtext("E", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( genome.sizes[ colnames(d.clup.266[,lab.i]) ], d.clup.266['o.e',lab.i], pos=1 )
plot( ex.align.2.k2[ 'denticeps_clupeoides', colnames(d.clup.266) ], d.clup.266['o.e',], xlab='Kimura distance',
     ylab='', col=d.clup.266.col, pch=19 )
with(par(), mtext("F", cex=mt.cex, at=usr[1], adj=1, line=0.2))
text( ex.align.2.k2[ 'denticeps_clupeoides', colnames(d.clup.266[,lab.i]) ], d.clup.266['o.e',lab.i], pos=1 )
dev.off()


