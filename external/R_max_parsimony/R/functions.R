## some of these require loading of the ../src/max_parsimony.so
## shared object

## encode numbers as characters.. 
encode.dist <- function(data, offset=64, conversion=as.integer, missing.value=0){
    data <- apply(data, 2, conversion)
    data[ is.na(data) ] <- missing.value
    data <- data + offset
    ##    apply(data, 2, intToUtf8)
    ## intToUtf8 may do strange things with values above
    ## 127. Instead use a simpler function that does
    ## abs(i) % 255
    ## to get char structures.
    if(!is.integer(data))
        data <- matrix( as.integer(data), ncol=ncol(data) )
    apply(data, 2, function(x){ .Call("intToAscii", x) })
}


make.sub.matrix <- function(size){
    m <- matrix(nrow=size, ncol=size)
    for(i in 1:size){
        for(j in 1:size){
            m[i,j] <- as.integer(abs(i-j))
        }
    }
    ## missing value encoding should be handled by the user as it is
    ## difficult to predict how this should be done here.
    m
}

## edges:       A matrix of edge information as created by the ape nj function.
##              This contains two columns giving: 1: the parent node index,
##              and 2: the child index.
## node.info:   a vector containing 1) the total number of nodes and 2) the
##              the number of leaf nodes.
## sub.matrix:  A square substitution matrix.
## al.info:     A vector containing the alphabet offset (minimum value) and the
##              size of the alphabet. The width and height of the substitution matrix
##              should be the size of the alphabet.
## leaf.states: The states of each leaf node encoded as an ASCII string. Note that no
##              letter in the states should exceed offset + al.size - 1
## All arguments should be integral.
sankoff <- function(edges, node.info, sub.matrix, al.info, leaf.states, check.states=TRUE){
    if(!is.character(leaf.states))
        stop("leaf.states should be a character vector\n")
    state.range <- range(unlist(sapply(leaf.states, utf8ToInt)))
    if(check.states && (state.range[1] < al.info[1] || state.range[2] >= sum(al.info)))
        stop("leaf.stats should not have any characters outside of the alphabet")
    .Call("sankoff", edges, as.integer(node.info), sub.matrix, as.integer(al.info), leaf.states)
}

## extracts suitable coordinates for drawing a dendrogram from a table of edges
## as provided by the ape nj() function. This gives the coordinates used by the
## ape plot( tree ) function.

edge.lines <- function(tree){
    y <- vector(mode='numeric', length=nrow(tree$edge))
    nodes <- y
    x <- matrix(0, nrow=nrow(tree$edge), ncol=4)
    colnames(x) <- c('x1', 'x2', 'p', 'c')
    v.lines <- matrix(nrow=0, ncol=4)
    colnames(v.lines) <- c('x', 'y1', 'y2', 'node')
    leaf.b <- tree$edge[,2] < min(tree$edge[,1])
    y[leaf.b] <- 1:sum(leaf.b)
    ## start from the root...
    ## the root has no parent..
    root <- setdiff( tree$edge[,1], tree$edge[,2] )
    visit.tree <- function(root, r.x ){
        ## child.i is wrong. lets call it this.i
        this.i <- which( tree$edge[,2] == root )
        nodes[this.i] <<- root
        root.i <- which( tree$edge[,1] == root )
        if(!length(root.i)){
            return(y[ this.i ])
        }
        x[root.i,1] <<- r.x
        x[root.i,2] <<- x[root.i,1] + tree$edge.length[root.i]
        x[root.i,'p'] <<- tree$edge[root.i,1]
        x[root.i,'c'] <<- tree$edge[root.i,2]
        children <- tree$edge[root.i, 2]
        child.y <- vector(length=length(root.i))
        for(i in 1:length(root.i))
            child.y[i] <- visit.tree( children[i], x[root.i[i], 2] )
        y[ this.i ] <<- mean(child.y)
        v.lines <<- rbind(v.lines, c(x[ root.i[1], 1], min(child.y), max(child.y), root ))
        ##return( y[this.i] )
        return( mean(child.y) )
    }
    y.root <- visit.tree(root, 0)
    nodes <- c(nodes, root)
    x <- rbind(x, c(-max(x)/50, 0, root, root))
    y <- c(y, y.root)
    list('x'=x, 'y'=y, 'nodes'=nodes, 'v'=v.lines)
}

