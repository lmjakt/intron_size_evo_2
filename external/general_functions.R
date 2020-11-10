## a set of general purpose functions

## makes a color for each of level of v
## with low (blue) to high (purple) via, cyan, green, yellow, red.
## this can also be done by reordering the
## output of the rainbow function.
## but not sure how to get the radial shift.
hsvScale <- function(v, sat=1, val=0.75, alpha=1,
                     min.v=min(v, na.rm=TRUE), max.v=max(v, na.rm=TRUE),
                     na.col=rgb(0.6,0.6,0.6)){
    ## handle na values gracefully
    na.v <- is.na(v)
    v[ na.v ] <- mean( v[ !na.v ] ) ## should not change the colour scale
    ## ugly hack
    if(is.na(min.v)){
        min.v <- min(v)
        max.v <- max(v)
    }
    val[ is.na(val) ] <- 0
    sat[ is.na(sat) ] <- 0
    ## run from blue (4/6) -> magenta (5/6)
    v.range <- max.v - min.v
    if(!v.range)
        return(rep(hsv(4/6, sat, val, alpha), length(v)))
    ## flatten the scale outside of the specified range
    v[ v < min.v ] <- min.v
    v[ v > max.v ] <- max.v
    v <- 5 + 5 * (min.v - v) / v.range
    ## now runs from 5 (magenta) -> 0 (red)
    ## convert to 4, 3, 2, 1, 0, 5
    v <- (v - 1) ## and now runs 4, 3, 2, 1, 0, -1
    v[ v < 0 ] <- 6 + v[ v < 0 ] ## -> 4, 3, 2, 1, 0, 5
    cols <- hsv( v/6, sat, val, alpha )
    cols[na.v] <- na.col
    cols
}


## tapply can be used to do something similar.
## the main difference being that using tapply will sort the
## levels through as.numeric(as.factor).
## Note that tapply should be much faster and a better choice for
## anything big.
## returns a vector of numerical factor values that identify
## unique ids for the columns specified by c.id
mulFactor <- function(m, c.id){
    m.u <- unique(m[,c.id, drop=FALSE])
    m.f <- vector(mode='numeric', length=nrow(m))
    for(i in 1:nrow(m.u)){
        b <- rep(TRUE, nrow(m))
        for(j in c.id)
            b <- b & (m[,j] == m.u[i,j])
        m.f[b] <- i
    }
    ## some discrete mathematics (set theory) could probably help this,
    ## but this is written simply
    list('f'=m.f, 'm'=m.u)
}


## calculate f statistics from a table of values
## with grouped values.

## d is a table of numerical values (matrix)
## d.g is a vector of ncol(d) length containing the
## group membership of the columns.
## adjustIntVarFunc should be a function that performs an adjustment on the internal
## variance to compensate for weirdness in this. Normally this, would be used to
## try to get around the problem of too low variance for certain types of data
## where we get close to 0 variance.
fStats <- function(d, d.g, warnOnZeroVar=FALSE, adjustIntVarFunc=eval){
    ## calculate the means of the groups:
    g.ids <- sort(unique(d.g))
    g.sizes <- vector(mode='numeric', length=length(g.ids))
    g.means <- matrix(nrow=nrow(d), ncol=length(g.ids))
    means <- rowSums( d ) / ncol( d ) ## faster than rowMeans, but not as safe
    var.within <- rep(0, nrow(d))
    var.between <- var.within
    N <- ncol( d )        ## the number of samples
    K <- length( g.ids )  ## the number of groups
    for(i in 1:length(g.ids)){
        c.i <- which(d.g == g.ids[i])
        g.sizes[i] <- length(c.i)
        g.means[,i] <- rowSums(d[ ,c.i, drop=FALSE]) / g.sizes[i]
        ## we also need the internal variances, but the way we combine these
        ## means that we can't just cal var if I recall correctly.
        var.within <- var.within + ( rowSums((d[,c.i, drop=FALSE] - g.means[,i])^2) / (N - K) )
        var.between <- var.between + ( g.sizes[i] * ((g.means[,i] - means)^2) / (K - 1) )
    }
    var.within <- adjustIntVarFunc( var.within )
    f <- var.between / var.within
    ## if there is no variance between f should be 0, regardless of variance
    ## within
    f[ var.between == 0 ] <- 0
    ## if there is no variance within we really prefer not to have 
    f[ var.within == 0 ] <- 0
    ## note that var.within and var.between can themselves have NaNs due to the division
    ## we should provide some option to allow the user to handle this.
    if(warnOnZeroVar){
        print("0 variance within for: ")
        print(which(var.within == 0))
    }
    return( list('f'=f, 'df1'=K-1, 'df2'=N-K) )
}

## identical to fStats, but for a maximum of two classes and thus
## gives a directional test.
## Equation taken from https://en.wikipedia.org/wiki/Student%27s_t-test
## Welch's t-test (equal or unequal samples sizes, unequal variances)
tStats <- function(d, d.g, warnOnZeroVar=FALSE, adjustIntVarFunc=eval){
    g.ids <- sort(unique(d.g))
    if(length(g.ids) != 2){
        stop("tStats must specify exactly two groups")
    }
    ## only two groups so we then use 1 and 2 here
    c.b <- (d.g == g.ids[1])
    n1 <- length(which(c.b))
    n2 <- length(which(!c.b))
    g.means.1 <- rowSums(d[,c.b, drop=FALSE]) / n1
    g.means.2 <- rowSums(d[,!c.b, drop=FALSE]) / n2
    var.1 <- rowSums((d[,c.b, drop=FALSE] - g.means.1)^2) / n1
    var.2 <- rowSums((d[,!c.b, drop=FALSE] - g.means.2)^2) / n2
    var.1 <- adjustIntVarFunc(var.1)
    var.2 <- adjustIntVarFunc(var.2)
    tStat <- (g.means.1 - g.means.2) / sqrt( var.1 + var.2 )
    tStat
}

plot.kmeans <- function(km, exp, column.order=NULL, ...){
    if(is.null(column.order))
        column.order = 1:ncol(km$centers)
    for(i in 1:nrow(km$centers)){
        plot(1:ncol(km$centers), km$centers[i,column.order], ylim=range(exp[km$cluster, column.order]),
             type='n', main=paste(i, ":", km$size[i], "genes"), ...)
        b <- km$cluster == i
        b.n <- sum( b )
        e.t <- t(exp[b, ])
        segments( rep( 1:(ncol(exp)-1), b.n ),
                  e.t[ 1:(ncol(exp)-1), ],
                  rep( 2:ncol(exp), b.n ),
                  e.t[ 2:ncol(exp), ],
                  col='grey')
        ## apply(exp[ km$cluster == i, ], 1, function(x){
        ##     lines(1:ncol(km$centers), x[column.order], col='grey') })
        lines(1:ncol(km$centers), km$centers[i,column.order],
              col='red', lwd=2)
        
    }
}


## let's make a table for the klusters..
extract.clusters <- function(km, exp, genes){
        ## e is the expression data,, but the first column is the cluster
        e <- matrix(nrow=length(km$cluster), ncol=(1+ncol(exp)), data=0)
        colnames(e) <- c('cluster', colnames(exp))
        d <- matrix(nrow=length(km$cluster), ncol=2) ## the
        colnames(d) <- c("gene_id", "gene")
        r <- 1
        for(i in 1:nrow(km$centers)){
                    r.b <- km$cluster == i
                            r.r <- r:(r + sum(r.b) - 1)
                            d[ r.r, 1] <- names(km$cluster[r.b])
                            d[ r.r, 2] <- genes[ d[r.r,1], 'gene']

                    e[ r.r, 1] <- i
                            e[ r.r, 2:ncol(e)] = as.matrix(exp[ d[r.r,1] , ])
                            r <- r + sum(r.b)
        }
        data.frame(d, e)
}

rowVars <- function( d ){
    m <- rowMeans( d )
    rowSums( ((m - d)^2) ) / (ncol(d) - 1)
}

## a function to calculate mirrored distance matrices
dyn.load("~/R/cpp/col_dist.so")
colDist <- function( m, scale.columns=FALSE ){
    if(!is.matrix( m ))
        stop('data should be supplied as a matrix of double values')
    if(!is.double(m))
        m <- matrix( as.double(m), nrow=nrow(m) )
    if(scale.columns)
        m <- scale(m)
    
    d <- .Call("col_dist", m )
    if( !is.null(colnames(m)) ){
        rownames(d) <- colnames(m)
        colnames(d) <- colnames(m)
    }
    d
}

## a function that provides a kernel smoothened distribution
## the user must specify a set of values, a set of mid points for
## the distribution, a half-kernel, and the step size of the kernel
## as related to the dimensions of the values and mid sizes.
dyn.load("~/R/cpp/smoothed_hist.so")
smoothHist <- function( values, mids, kernel, kernel.unit ){
    values <- as.double( values )
    mids <- as.double( mids )
    kernel <- as.double( kernel )
    kernel.unit <- as.double( kernel.unit )
    if( !length(values) || !length(mids) || !length(kernel) || !length(kernel.unit))
        stop("empty vectors are not allowed")
    .Call( "smoothed_hist", values, mids, kernel, kernel.unit )
}


## a method to determine clusters of points in linear space by a neareest neighbour method,
## where nearest is defined by a combiantion of proximity and similarity by simply
## allowing links to skip proximal neighbours.
dyn.load("~/R/cpp/linear_nn.so")
linearNN <- function( delta, maxSkip ){
    if(!is.matrix(delta))
        stop("linearNN: delta should be a matrix")
    if(!is.double(delta))
        delta <- matrix( as.double(delta), nrow=nrow(delta))
    maxSkip = as.integer( maxSkip )

    .Call( "linear_nn", delta, maxSkip )
}

## a threshold based method that simply makes clusters on the basis of the longest
## internal distance within a set of adjacent positions.
dyn.load("~/R/cpp/grow_linear_clusters.so")
growLinearClusters <- function( delta, max.delta ){
    if( !is.matrix(delta) )
        stop("delta must be a matrix");
    if( !is.double( delta ) )
        delta <- matrix( as.double(delta), nrow=nrow(delta) )
    max.delta <- as.double( max.delta )

    .Call("grow_linear_clusters", delta, max.delta )
}
