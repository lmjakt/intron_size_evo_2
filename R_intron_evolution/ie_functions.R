## this requires:
dyn.load("~/R/c/convolve_ks.so")
source("~/R/general_functions.R")

## for convenience
## inverts range -> seq
xrange <- function(r){
    r[1]:r[2]
}

## make a simulate function on the basis of the distribution of
## an ancestor and descendant pair
simulate.is.evol <- function(anc, desc, gen.n, min.s=10*log2(76),
                             inf=int.s.inf.2, blur.sd=6, blur.w=10,
                             size.depend=TRUE, norm.sd=NA, use.blurred=TRUE,
                             from.inf=TRUE){
    if(from.inf){
        int.s <- cbind(inf[[anc]]$state[,1], inf[[desc]]$state[,1])
    }else{
        int.s <- cbind(anc, desc)
    }
    int.s <- int.s[ apply(int.s, 1, min) > min.s, ]
    int.d <- int.s[,2] - int.s[,1]
    int.s.r <- range(int.s)
    int.d.r <- range(int.d)
    dummy <- suppressWarnings( cbind(xrange(int.s.r), xrange(int.d.r) ) )
    dummy.tbl <- table(dummy[,1], dummy[,2])
    s.d <- table( c(dummy[,1], int.s[,1]), c(dummy[,2], int.d) )
    s.d <- s.d - dummy.tbl
    s1.d <- table( c(xrange(int.d.r), int.d) ) - 1
    s1.d <- s1.d / sum(s1.d)
    kernel <- dnorm(0:blur.w, mean=0, sd=blur.sd)
    s.d.blurred <- apply(s.d, 2, function(x){ .Call('convolve_ks', as.double(x), kernel) })
    s.d.blurred <- t(apply(s.d.blurred, 1, function(x){ .Call('convolve_ks', as.double(x), kernel) }))
    s.d.blurred[ rowSums(s.d.blurred) == 0, ] <- 1
    s.d.blurred <- s.d.blurred / rowSums(s.d.blurred)
    s.d[ rowSums(s.d) == 0, ] <- 1
    s.d <- s.d / rowSums(s.d)
    x <- xrange(int.s.r)
    y <- xrange(int.d.r)
    ## evolve by sampling from the distributions:
    make.descendants <- function(anc.i, gen, max.i){
        if(gen > gen.n)
            return(max.i)
        a <- nodes[[anc.i]]$size
        a[ a < min(x) ] <- min(x)
        a[ a > max(x) ] <- max(x)
        if(size.depend){
            if(use.blurred){
                d.d <- sapply( a, function(l){
                    sample( y, size=2, replace=TRUE, prob=s.d.blurred[ 1 + l - x[1], ] )})
            }else{
                d.d <- sapply( a, function(l){
                    sample( y, size=2, replace=TRUE, prob=s.d[ 1 + l - x[1], ] )})
            }                
        }else{
            if(is.na(norm.sd)){
                d.d <- matrix( sample(y, size=2*length(a), replace=TRUE, prob=s1.d), nrow=2 )
            }else{
                d.d <- matrix( rnorm(n=2*length(a), sd=norm.sd), nrow=2 )
            }
        }
        d1 <- a + d.d[1,]
        d2 <- a + d.d[2,]
        n1 <- max.i + 1
        n2 <- max.i + 2
        nodes[[n1]] <<- list(parent=anc.i, node=n1, size=d1, gen=gen+1)
        nodes[[n2]] <<- list(parent=anc.i, node=n2, size=d2, gen=gen+1)
        max.i <- make.descendants( n1, gen+1, n2 )
        make.descendants( n2, gen+1, max.i )
    }
    nodes <- vector(mode='list', length=sum( 2^(0:gen.n) ))
    nodes[[1]] <- list(parent=1, node=1, size=int.s[,1], gen=1)
    make.descendants( 1, 1, 1 )
    nodes
}

## If we want to infer back from the leaves of a tree produced by
## simulate.is.evol()
dyn.load("~/R/R_max_parsimony/src/max_parsimony.so")
source("~/R/R_max_parsimony/R/functions.R")

## a wrapper function; reorganises the tree and calls sankoff()
infer.ancestors <- function(desc){
    ## we need an edge table, this should have:
    ## parent, node
    edge <- t(sapply(desc, function(x){ c(x$parent, x$node) }))
    ## remove the first circular root node..
    ## alternatively modify the simulate.is.evol function
    edge <- edge[ edge[,1] != edge[,2], ]
    ## leaf nodes need to be numbered 1..leaf.n
    ## in order that they can correspond to a table of
    ## of states
    leaf.i <- which( !(edge[,2] %in% edge[,1]) )
    ## create a rearranged edge table:-> edge.ra
    edge.ra <- edge + length(leaf.i)
    edge.ra[ leaf.i, 2 ] <- 1:length(leaf.i)
    ##  I also need to make sure that the index of the nodes does not
    ##  exceed the number of nodes..
    total.nodes <- length(unique(as.integer(edge.ra)))
    node.remap <- rep(NA, max(edge.ra))
    node.remap[ sort(unique(as.integer(edge.ra))) ] <- 1:total.nodes
    edge.ra <- cbind( node.remap[ edge.ra[,1] ], node.remap[ edge.ra[,2] ])
    ## then obtain the leaf states and prepare the alphabet.
    states <- sapply( desc[leaf.i + 1], function(x){ x$size })
    al.size <- as.integer( 1 + max(states))
    ##    al.size <- as.integer( 1 + diff( range(states) ))
    ## I should rewrite the max parsimony function so that it can handle
    ## non 0 starts: note that we can run out of sizes here
    ## if the largest intron is larger than 254
    al.offset <- 1L
    sub.matrix <- make.sub.matrix( al.size )
    missing.penalty <- as.integer( max(sub.matrix) + 1 )
    sub.matrix[1,] <- missing.penalty
    sub.matrix[,1] <- missing.penalty
    states.enc <- encode.dist( states, offset=al.offset )
    leaf.n <- sum( !(edge[,2] %in% edge[,1]) )
    total.nodes <- length(unique(as.integer(edge)))
    inf <- sankoff( edge.ra, c(total.nodes, leaf.n), sub.matrix,
                   c(al.offset, al.size), states.enc, check.states=FALSE)
    list(inf=inf, leaves=states, edge=edge, edge.ra=edge.ra)
}

## statistics for common retention of long introns in real and
## simulated data
## it does look like we may have an issue. To deal with this, the simplest
## thing is to divide up into separate sections..

## anc, d1, and d2 are the sizes of the introns in the ancestor and the two descendants
common.retention <- function(anc, d1, d2, thr=80, breaks=seq(80, 200, 10), min.l=10*log2(76)){
    int.l <- cbind(anc, d1, d2)
    int.l <- int.l[ apply(int.l, 1, min) >= min.l, ]
    sapply( 2:length(breaks), function(i){
        b <- int.l[,1] >= breaks[i-1] & int.l[,1] < breaks[i]
        q <- sum( b & int.l[,2] >= thr & int.l[,3] >= thr )
        m <- sum( b & int.l[,2] >= thr )
        n <- sum( b ) - m
        k <- sum(b & int.l[,3] >= thr )
        exp <- k * m / sum(b)
        r <- q/exp
        p <- phyper( q-1, m, n, k, lower.tail=FALSE )
        c(q=q, m=m, n=n, k=k, exp=exp, r=r, p=p)
    })
}

common.retention.tree <- function(anc, d1, d2, inf=int.s.inf.2, thr=80, breaks=seq(80, 200, 10), min.l=10*log2(76),
                             is.simulated=FALSE){
    if(is.simulated){
        int.l <- sapply(inf[c(anc, d1, d2)], function(x){ x$size })
    }else{
        int.l <- sapply( inf[c(anc, d1, d2)], function(x){ x$state[,1] })
    }
    common.retention(int.l[,1], int.l[,2], int.l[,3], thr=thr, breaks=breaks)
}


test.common.retention <- function(anc, d1, d2, inf=int.s.inf.2,
                                  thr=80, breaks=seq(80, 200, 10), min.l=10*log2(76),
                                  sim.n=10){
    int.s <- sapply( c(anc, d1, d2), function(i){ inf[[i]]$state[,1] })
    int.s <- int.s[ apply(int.s, 1, min) > min.l, ]
    obs.ret <- common.retention( int.s[,1], int.s[,2], int.s[,3], thr, breaks, min.l )
    sim.d1 <- lapply(1:sim.n, function(x){
        simulate.is.evol(int.s[,1], int.s[,2], gen.n=1, use.blurred=FALSE, from.inf=FALSE)
    })
    sim.d2 <- lapply(1:sim.n, function(x){
        simulate.is.evol(int.s[,1], int.s[,3], gen.n=1, use.blurred=FALSE, from.inf=FALSE)
    })
    ## then do common retention for all sim.d1 against all sim.d2
    anc.size <- sim.d1[[1]][[1]]$size
    if(any(anc.size != sim.d2[[1]][[1]]$size))
        stop("The ancestral sizes in the simulations do not match")
    ##
    sim.ret <- lapply(sim.d1, function(x){
        lapply(sim.d2, function(y){
            lapply( list( c(2,2), c(2,3), c(3,2), c(3,3)), function(z){
                common.retention(anc.size, x[[z[1]]]$size, y[[z[2]]]$size,
                                 thr=thr, breaks=breaks, min.l=min.l)
            })
            ## at this point I probably have a list of tables.
            ## I probably want to keep all of these for the time being, and extract
            ## the relevant things later on.
        })
    })
    sim.p <- t(do.call(cbind, lapply( sim.ret,
                                     function(x){
                                         do.call(cbind, lapply(x, function(y){
                                             sapply(y, function(z){ z['p',] }) })) })))
    sim.r <- t(do.call(cbind, lapply( sim.ret,
                                     function(x){
                                         do.call(cbind, lapply(x, function(y){
                                             sapply(y, function(z){ z['r',] }) })) })))
    list(obs=obs.ret, sim=sim.ret, sim.p=sim.p, sim.r=sim.r, breaks=breaks)
}

## collect specific rows from the sim part of the data structure
## returned by test.common.retention (tcr) for row rn
collect.tested.rows <- function(tcr, rn){
    t(do.call(cbind, lapply(tcr$sim,
                            function(x){
                                do.call(cbind, lapply(x, function(y){
                                    sapply(y, function(z){ z[rn,] }) })) })))
}

plot.tested.stat <- function(tcr, stat.name='p', sp1, sp2, ...){
    if(!(stat.name %in% c('p', 'r', 'f1', 'f2', 'qf')))
        stop("unsupported stat in plot.tested.stat")
    if(stat.name == 'p'){
        st <- -log10(tcr$obs['p',])
        st.sim <- -log10(tcr$sim.p)
        ylab='-log10(p)'
    }
    if(stat.name == 'r'){
        st <- tcr$obs['r',]
        st.sim <- tcr$sim.r
        ylab='observed / expected'
    }
    ## f for fraction retained in sp 1
    if(stat.name == 'f1'){
        st <- with(tcr, obs['m',] / (obs['m',] + obs['n',]))
        m.sim <- collect.tested.rows( tcr, 'm' )
        n.sim <- collect.tested.rows( tcr, 'n' )
        st.sim <- m.sim / (m.sim + n.sim)
        ylab=paste('retained in', sp1)
    }
    if(stat.name == 'f2'){
        st <- with(tcr, obs['k',] / (obs['m',] + obs['n',]))
        m.sim <- collect.tested.rows( tcr, 'm' )
        n.sim <- collect.tested.rows( tcr, 'n' )
        k.sim <- collect.tested.rows( tcr, 'k' )
        st.sim <- k.sim / (m.sim + n.sim)
        ylab=paste('retained in', sp2)
    }
    if(stat.name == 'qf'){
        st <- with(tcr, obs['q',] / (obs['m',] + obs['n',]))
        m.sim <- collect.tested.rows( tcr, 'm' )
        n.sim <- collect.tested.rows( tcr, 'n' )
        q.sim <- collect.tested.rows( tcr, 'q' )
        st.sim <- q.sim / (m.sim + n.sim)
        ylab=paste('retained in', sp2)
    }
    main <- paste(sp1, "&", sp2)
    xlab <- 'ancestral size'
    x <- tcr$breaks[-1] - diff(tcr$breaks)[1]/2
    plot(x, st, xlim=range(tcr$breaks), ylim=range(c(st, st.sim), na.rm=TRUE),
         xlab=xlab, ylab=ylab, main=main, type='n', ...)
    invisible(apply(st.sim, 1, function(y){ lines(x, y, col='gray', lty=3) }))
    lines(x, st)
}

collect.leaves <- function(node, tree=ex.align.2.k2.nj){
    descend.tree <- function(i){
        j <- which(tree$edge[,1] == tree$edge[i,2])
        if(length(j) == 0){
            leaves <<- c(leaves, tree$edge[i,2])
        }
        for(k in j){
            descend.tree(k)
        }
    }
    leaves <- c()
    i <- which(tree$edge[,2] == node)
    descend.tree(i)
    list(i=leaves, names=tree$tip.label[leaves])
}

collect.descendants <- function(node, tree=ex.align.2.k2.nj, leaf.only=FALSE){
    descend.tree <- function(i){
        j <- which(tree$edge[,1] == tree$edge[i,2])
        descendants <<- c(descendants, tree$edge[i,2])
        for(k in j){
            descend.tree(k)
        }
    }
    descendants <- c()
    i <- which(tree$edge[,2] == node)
    descend.tree(i)
    descendants <- descendants[-1]
    if(leaf.only)
        descendants <- descendants[ descendants <= length(tree$tip.label) ]
    list(i=descendants, names=tree$tip.label[descendants])
}

mk.species.label <- function(labels, short=TRUE){
    labels.sp <- strsplit( labels, "[ _.]+" )
    sapply( labels.sp, function(x){
        gn <- paste(toupper( substring(x[1], 1, 1) ), ".", sep="")
        paste(c(gn, x[-1]), collapse=" ")
    })
}

## lm, rm, bm, tm are the left, right, bottom, and top margins
## expressed as a proportion of the plotted width and height
draw.tree <- function( tree=ex.align.2.k2.nj, t.lines=tree.lines, root=NULL,
                      lm=0.04, rm=0.25, bm=0.04, tm=0.04, label.nodes=TRUE,
                      label.leaves=TRUE, p.new=TRUE,
                      leaf.cex=0.5, node.cex=leaf.cex,
                      lab.f=mk.species.label, node.l.offset=0.001, leaf.l.offset=node.l.offset,
                      leaf.font=3){
    if(!is.null(root)){
        sel.nodes <- collect.descendants(root, tree, leaf.only=FALSE)
        t.lines$y <- with(t.lines, y[ x[,'c'] %in% sel.nodes$i ])
        t.lines$nodes <- with(t.lines, nodes[ x[,'c'] %in% sel.nodes$i ] )
        t.lines$x <- with(t.lines, x[ x[,'c'] %in% sel.nodes$i, ])
        t.lines$v <- with(t.lines, v[ v[,'node'] %in% x[,'p'], ])
    }
    yr <- range( t.lines$y )
    xr <- range( t.lines$x[,1:2] )
    yr.d <- diff(yr)
    xr.d <- diff(xr)
    if(p.new)
        plot.new()
    with(t.lines, plot.window(xlim=xr + c(-lm, rm) * xr.d,
                              ylim=yr + c(-bm, tm) * yr.d, xaxs='i', yaxs='i'))
    lab.r <- with(t.lines, {
        segments( x[,1], y, x[,2], y )
        segments( v[,1], v[,2], v[,1], v[,3])
        node.labels.w <- strwidth( nodes, cex=node.cex )
        if(!label.nodes)
            node.labels.w <- 0
        node.labels.r <- x[,2] + node.l.offset + node.labels.w
        if(label.nodes)
            text( x[,2] + node.l.offset, y, labels=nodes, adj=c(0,0.5), cex=node.cex )
        b <- !(nodes %in% tree$edge[,1])
        leaf.labels <- lab.f(tree$tip.label[nodes[b]])
        leaf.labels.w <- strwidth( leaf.labels, cex=leaf.cex )
        if(label.leaves)
            text( node.labels.r[b] + leaf.l.offset, y[b], lab.f(tree$tip.label[nodes[b]]), cex=leaf.cex, adj=c(0,0.5),
                 font=leaf.font)
        list(node.right=node.labels.r,
             leaf.right=node.labels.r[b] + leaf.l.offset + leaf.labels.w,
             node.w=node.labels.w, leaf.w=leaf.labels.w,
             leaf.labels=leaf.labels, leaf.b=b, sp=tree$tip.label[nodes[b]])
    })
    invisible(list(lines=t.lines, nodes=lab.r, usr=par('usr')))
}

## Visualise the changes between ancestors and descendants:
## instead use my own convolve function
dyn.load("~/R/c/convolve_ks.so")

## for discrete integer values only
hist.2d <- function(x, y, min.x=10*log2(76), min.y=min.x){
    b <- !(is.na(x) | is.na(y))
    if(!is.na(min.x))
        b <- b & x >= min.x 
    if(!is.na(min.y))
        b <- b & y >= min.y 
    x <- x[b]
    y <- y[b]
    x.r <- range(x)
    y.r <- range(y)
    dummy <- suppressWarnings( cbind( x.r[1]:x.r[2], y.r[1]:y.r[2] ))
    dummy.tbl <- table(dummy[,1], dummy[,2])
    h2d <- table( c(dummy[,1], x), c(dummy[,2], y))
    h2d <- h2d - dummy.tbl
    list(x=x.r[1]:x.r[2], y=y.r[1]:y.r[2], h=h2d, xv=x, yv=y)
}

## this uses ancestral values for sd and min.f
## NOT GOOD.
trace.modes.x <- function(x, y, m, lwd=2){
    pk <- apply(m, 1, function(x){ get.peaks( y, x, 3, 0.5) })
    abline(0,1,col=navy, lty=2, lwd=lwd)
    b.i <- which(sapply(pk, function(x){ length(x) > 0 }))
    pts <- lapply(b.i, function(i){ cbind( x[i], pk[[i]]) })
    ridges <- trace.lines(pts, diff(range(y))*0.1)
    ridge.w <- sapply( ridges, function(r.x){ sum( m[ cbind(1+r.x[,1]-x[1], 1+r.x[,2]-y[1]) ])})
    b <- (ridge.w / max(ridge.w)) > 0.05
    for(k in which(b)){
        ##        lines(ridges[[i]], col='blue', lwd=3)
        lines(ridges[[k]][,1]/10, .Call("gs_smooth", ridges[[k]], sd, min.f)/10, col=navy, lwd=lwd)
    }
    invisible(list(pts=ridges, wt=ridge.w, b=b))
}


## for discrete integer values only
## this plots with discrete distributions shown as well
## uses split.screen; cannot be used with par() / layout()
hist.2d.2 <- function(x, y, screen.n=NULL, min.x=10*log2(76), min.y=min.x,
                      dist.w=0.15, dist.h=dist.w, do.plot=TRUE,
                      cols=hcl.colors(256, "YlOrRd", rev=TRUE),
                      mar=c(5.1, 4.1, 4.1, 2.1), main=NA,
                      h.tform=eval, image.after=NA, ax.cex=1,
                      ax.mgp=c(3,1,0), ax.tcl=-ax.mgp[2],
                      xlab=NA, ylab=NA, trace.modes=FALSE, x.div=1, y.div=x.div,
                      cex.main=0.75, cex.xy.lab=0.75,
                      panel.label=NA, panel.cex=2, panel.line=1, panel.padj=0, panel.adj=1, panel.side=2
                      ){
    d.range <- function(x){
        c(min(x), 1.1*max(x))
    }
    h <- hist.2d( x, y, min.x, min.y )
    ## list with x, y, h, xv and yv
    ## h is a table where each row correspondes to one value in x
    ## and each column one value in y.
    h.x <- rowSums(h$h)
    h.y <- colSums(h$h)
    if(!is.null(screen.n)){
        screens <- split.screen( rbind(c(0, dist.w, dist.h, 1),
                                       c(dist.w, 1, 0, dist.h),
                                       c(dist.w, 1, dist.h, 1)),
                                screen.n)
    }else{
        screens <- split.screen( rbind(c(0, dist.w, dist.h, 1),
                                       c(dist.w, 1, 0, dist.h),
                                       c(dist.w, 1, dist.h, 1))
                                )
    }
    screen(screens[3])
    par(mar=c(0,0,mar[3:4]))
    ht <- h.tform(h$h)
    image(h$x/x.div, h$y/y.div, ht, col=cols, axes=FALSE, main=main, cex.main=cex.main)
    eval(quote(image.after))
    ridges <- NULL
    if(trace.modes)
        ridges <- trace.modes.x(h$x, h$y, ht)
    screen(screens[1])
    par(mar=c(0, mar[2:3], 0))
    plot.window(xlim=d.range(h.y), ylim=range(h$y/y.div), yaxs='i', xaxs='i')
    polygon(1.1 * max(h.y) - c(min(h.y), h.y), c(min(h$y), h$y)/y.div, col='grey')
    axis(2, cex.axis=ax.cex, tcl=ax.tcl, mgp=ax.mgp)
    mtext(ylab, side=2, line=ax.mgp[2]+1, cex=cex.xy.lab)
    if(!is.na(panel.label)){
        with(par(), mtext(panel.label, side=panel.side, las=1, line=panel.line,
                          at=usr[4], cex=panel.cex, padj=panel.padj, adj=panel.adj))
    }
    screen(screens[2])
    par(mar=c(mar[1], 0, 0, mar[4]))
    plot.window(xlim=range(h$x/x.div), ylim=d.range(h.x), yaxs='i', xaxs='i')
    polygon(c(h$x[1], h$x)/x.div, 1.1 * max(h.x) - c(min(h.x), h.x), col='grey')
    axis(1, cex.axis=ax.cex, tcl=ax.tcl, mgp=ax.mgp)
    mtext(xlab, side=1, line=ax.mgp[2]+1, cex=cex.xy.lab)
    invisible(list(h=h, ht=ht, screens=screens, ridges=ridges))
}

## x and y, should be breaks that are suitable for images
## z should be a matrix of counts
image.marginal <- function(x, y, z,
                           dist.w=0.15, dist.h=dist.w, do.plot=TRUE,
                           cols=hcl.colors(256, "YlOrRd", rev=TRUE),
                           mar=c(5.1, 4.1, 4.1, 2.1), main=NA,
                           ax.cex=1,
                           ax.mgp=c(3,1,0), ax.tcl=-ax.mgp[2],
                           xlab=NA, ylab=NA, trace.modes=FALSE, x.div=1, y.div=x.div,
                           cex.main=0.75, cex.xy.lab=0.75,
                           panel.label=NA, panel.cex=2, panel.line=1, panel.padj=0, panel.adj=1, panel.side=2,
                           leg.cex=0.5, leg.ticks=5
                           ){
    d.range <- function(x){
        c(min(x), 1.1*max(x))
    }
    h.x <- rowSums(z)
    h.y <- colSums(z)
    grey <- rgb(0.8,0.8,0.8)
    screen.m <- rbind(c(0, dist.w, dist.h, 1),
                      c(dist.w, 1, 0, dist.h),
                      c(dist.w, 1, dist.h, 1),
                      c(0, dist.w, 0, dist.h))
    screens <- split.screen( screen.m )
    screen(screens[3])
    par(mar=c(0,0,mar[3:4]))
    image(x, y, z, col=cols, axes=FALSE, main=main, cex.main=cex.main)
    ## Left histogram
    screen(screens[1])
    par(mar=c(0, mar[2:3], 0))
    plot.window(xlim=d.range(h.y), ylim=range(y), yaxs='i', xaxs='i')
    with(par(), rect(usr[2]-h.y, y[-length(y)], usr[2], y[-1], col=grey, border=NA))
##    polygon(1.1 * max(h.y) - c(min(h.y), h.y), c(min(h$y), h$y)/y.div, col='grey')
    axis(2, cex.axis=ax.cex, tcl=ax.tcl, mgp=ax.mgp)
    mtext(ylab, side=2, line=ax.mgp[2]+1, cex=cex.xy.lab)
    if(!is.na(panel.label)){
        with(par(), mtext(panel.label, side=panel.side, las=1, line=panel.line,
                          at=usr[4], cex=panel.cex, padj=panel.padj, adj=panel.adj))
    }
    screen(screens[2])
    par(mar=c(mar[1], 0, 0, mar[4]))
    plot.window(xlim=range(x), ylim=d.range(h.x), yaxs='i', xaxs='i')
    with(par(), rect(x[-length(x)], usr[4]-h.x, x[-1], usr[4], col=grey, border=NA))
##    polygon(c(h$x[1], h$x)/x.div, 1.1 * max(h.x) - c(min(h.x), h.x), col='grey')
    axis(1, cex.axis=ax.cex, tcl=ax.tcl, mgp=ax.mgp)
    mtext(xlab, side=1, line=ax.mgp[2]+1, cex=cex.xy.lab)
    ## set up a label:
    screen(screens[4])
    par(mar=c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim=c(0,1), c(0,1))
    s.y <- seq(0.15, 0.95, length.out=length(cols)) ## scale y
    s.yd <- diff(s.y)[1]
    s.x <- c(0.9, 1.0)
    rect(s.x[1], s.y, s.x[2], s.y+s.yd, col=cols, border=NA)
    rect(s.x[1], s.y[1], s.x[2], max(s.y+s.yd), col=NA, lwd=0.1)
    ## to determine suitable tick values..
    digits <- floor(2 - log10(max(abs(z))))
    ## constrain the tick ranges
    tick.range <- min(z) + c(0.05, 0.95) * (max(z) - min(z))
    t.v <- unique(round( seq(tick.range[1], tick.range[2], length.out=leg.ticks), digits=digits ))
    t.y <- s.y[1] + s.yd/2 + diff(range(s.y)) * (t.v - min(z)) / diff(range(z))
    segments(0.85, t.y, 0.9, t.y)
    text(0.84, t.y, t.v, adj=c(1,0.5), cex=leg.cex)
    close.screen(screens)
##    invisible(list(h=h, ht=ht, screens=screens, ridges=ridges))
}
        

diff.2d <- function(anc, desc, min.s=10*log2(76),
                    inf=int.s.inf.2, blur.sd=6, blur.w=10,
                    abline.args=NULL, do.plot=TRUE,
                    im.col=hcl.colors(256, "YlOrRd", rev=TRUE)
                    ){
    int.s <- cbind(inf[[anc]]$state[,1], inf[[desc]]$state[,1])
    b <- apply(int.s, 1, min) > min.s
    int.s <- int.s[b,] 
    int.d <- int.s[,2] - int.s[,1]
    ## int.s.r <- range(int.s)
    ## int.d.r <- range(int.d)
    ## dummy <- suppressWarnings( cbind(int.s.r[1]:int.s.r[2], int.d.r[1]:int.d.r[2] ) )
    ## dummy.tbl <- table(dummy[,1], dummy[,2])
    ## s.d <- table( c(dummy[,1], int.s[,1]), c(dummy[,2], int.d) )
    ## s.d <- s.d - dummy.tbl
    s.d <- hist.2d(int.s[,1], int.d, min.x=min.s, min.y=NA)
##    s1.d <- table( c(int.d.r[1]:int.d.r[2], int.d) ) - 1
##    s1.d <- s1.d / sum(s1.d)
    kernel <- dnorm(0:blur.w, mean=0, sd=blur.sd)
    s.d.blurred <- apply(s.d$h, 2, function(x){ .Call('convolve_ks', as.double(x), kernel) })
    s.d.blurred <- t(apply(s.d.blurred, 1, function(x){ .Call('convolve_ks', as.double(x), kernel) }))
    s.d.blurred <- s.d.blurred ## / rowSums(s.d.blurred)
    par(mfrow=c(2,3))
    x <- s.d$x ## int.s.r[1]:int.s.r[2] / 10
    y <- s.d$y ## int.d.r[1]:int.d.r[2] / 10
    if(do.plot){
        plot(int.s[,1]/10, int.d/10, cex=0.5, col=rgb(0,0,0,0.1))
        if(!is.null(abline.args)) do.call( abline, abline.args )
        plot(int.s[,2]/10, int.d/10, cex=0.5, col=rgb(0,0,0,0.1))
        if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, s.d$h, col=im.col )
        if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, s.d.blurred, col=im.col )
        if(!is.null(abline.args)) do.call( abline, abline.args )
        ## image(x=x, y=y, log(s.d.blurred) )
        ## if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, t(scale(t(s.d$h), center=FALSE)), col=im.col)
        ##    image(x=x, y=y, log( s.d / rowSums(s.d) ))
        if(!is.null(abline.args)) do.call( abline, abline.args )
        ## image(x=x, y=y, scale(s.d$h), col=im.col)
        ## if(!is.null(abline.args)) do.call( abline, abline.args )
        image(x=x, y=y, t(scale(t(s.d.blurred))), col=im.col)
        if(!is.null(abline.args)) do.call( abline, abline.args )
    }
    invisible(list(x=x, y=y, sh=s.d$h, bl=s.d.blurred, pt=cbind(as=int.s[,1], ds=int.s[,2], delta=int.d)))
}

## summarises change from x -> y discretized by levels in
## x
## it is assumed that x and y are log2 transformed values.
summarise.delta <- function(x, y, min.x=10*log2(76), min.y=min.x){
    b <- x >= min.x & y >= min.y
    x <- x[b]
    y <- y[b]
    delt <- sapply(sort(unique(x)), function(xv){
        b <- x == xv
        delta <- y[b] - x[b]
        tbl.delta <- table(delta)
        c('x'=xv,
          'sum'=sum(delta),
          'mean'=mean(delta),
          'median'=median(delta),
          'mode'=as.numeric( names(tbl.delta)[ which.max(tbl.delta) ] ),
          'npos'=sum(delta > 0),
          'nneg'=sum(delta < 0),
          'n'=sum(b),
          'x'=xv)
    })
    t(delt)
}

## it is assumed that x and y are log2 transformed values.
## and that the analysis will be for all discrete values
summarise.minimisation <- function(x, y, min.s=10*log2(100),
                                   min.x=10*log2(76), min.y=min.x,
                                   im.col=hcl.colors(256, "YlOrRd", rev=TRUE),
                                   do.plot=FALSE){
    b <- x >= min.x & y >= min.y
    x <- x[b]
    y <- y[b]
    xv <- xrange(range(x))
    yv <- xrange(range(y))
    dummy.v <- suppressWarnings( cbind(xv, yv ) )
    dummy.t <- table(dummy.v[,1], dummy.v[,2])
    tbl <- table( c(dummy.v[,1], x), c(dummy.v[,2], y) ) - dummy.t
    ## the rows of the table represent x, and the columns y
    yb <- yv <= min.s
    min.prop <- apply( tbl, 1, function(z){ sum(z[yb]) / sum(z) })
    min.prop[ is.na(min.prop) ] <-  0
##    par(mfrow=c(1,2))
    if(do.plot){
        image(xrange(xr), xrange(yr), tbl, col=im.col)
        abline(h=min.s)
        plot(xrange(xr), min.prop, type='l')
    }
    stats <- c( 'x.min'=sum(x <= min.s), 'y.min'=sum(y <= min.s),
                'xy.min'= sum(x <= min.s & y <= min.s ),
               'x.l'=sum(x), 'y.l'=sum(y), 'n'=length(x) )
    return(list(x=xv, y=yv, p=min.prop, tbl=tbl, stats=stats))
}

vis.minim <- function(anc, inf=int.s.inf.2, tree=ex.align.2.k2.nj, gs=genome.size,
                      cex.xylab=1, line.xylab=2, cex.leg=1, cex.axis=1, mgp.axis=c(3,1,0), tcl.axis=-mgp.axis[2],
                            mar=c(5.1, 4.1, 2.1, 4.1), max.p=0.7,
                      label.pos=NULL, lab.f=mk.species.label, lab.cex=0.75, sp.b=int.s.b,
                      int.x=0.4, int.y=0.4, int.w=0.6, int.h=0.6, pt.cex=0.75, lwd=1,
                      panel.label=NA, panel.cex=2, panel.line=1, panel.padj=0, panel.adj=1, panel.side=2){
    ## label.pos should contain a matrix with columns: leaf id, label pos
    leaves <- collect.leaves(anc, tree=tree)
    if(is.null(sp.b) || is.na(sp.b))
        sp.b <- rep(TRUE, max(leaves$i))
    i <- leaves$i
    i <- i[ sp.b[ i ] ]
    l.names <- with(leaves, names[ sp.b[i] ] )
    minim.summary <- lapply( i, function(i){
        summarise.minimisation( inf[[anc]]$state[,1], inf[[i]]$state[,1] )
    })
### 
    x <- minim.summary[[1]]$x / 10
    h <- with(inf[[anc]], hist(state[ state[,1] >= 10*log2(76),1] / 10, breaks=60, plot=FALSE ))
    par(mar=mar)
    screens <- split.screen( rbind(c(0,1,0,1), c(int.x, int.x+int.w, int.y, int.y+int.h)))
    screen(screens[1])
    plot.new()
    plot.window(xlim=range(h$breaks), ylim=range(h$count), yaxs='i', xaxs='i')
    mtext("Proportion minimised", side=2, line=line.xylab, cex=cex.xylab)
    mtext("Ancestral size", side=1, line=line.xylab, cex=cex.xylab)
    mtext("Ancestral count", side=4, line=line.xylab, cex=cex.xylab)
    if(!is.na(panel.label)){
        with(par(), mtext(panel.label, side=panel.side, las=1, line=panel.line,
                          at=usr[4], cex=panel.cex, padj=panel.padj, adj=panel.adj))
    }
    with(h, rect( breaks[-length(breaks)], 0, breaks[-1], counts, col=rgb(0.8, 0.8, 0.8), border=NA))
    axis(1, cex.axis=cex.axis, mgp=mgp.axis, tcl=tcl.axis)
    axis(4, cex.axis=cex.axis, mgp=mgp.axis, tcl=tcl.axis)
    plot.window( xlim=range(h$breaks), ylim=c(0,max.p), yaxs='i', xaxs='i' )
    axis(2, cex.axis=cex.axis, mgp=mgp.axis, tcl=tcl.axis)
    abline(v=log2(100), lty=2)
    cols <- hsvScale( log(gs[l.names]) )
    for(j in 1:length(i)){
        with(minim.summary[[j]], {
            y <-  .Call("convolve_ks", p, kernel)
            lines(x/10,
                  y,
                  type='l', col=cols[j], lwd=lwd) ##  nodes.col[j])
            print(paste(c(j, i[j], stats), collapse=" "))
        })
    }
    stats <- sapply(minim.summary, function(x){ x$stats })
    screen(screens[2])
    par(mar=c(0,0,0,0))
    plot( log2(gs[l.names]), stats[2,] / stats[1,], col=cols, pch=19, xlab='log2 genome size',
         ylab='n min. anc / desc', cex.axis=cex.axis, mgp=mgp.axis, tcl=tcl.axis, cex.lab=cex.axis,
         cex=pt.cex)
    abline(h=1)
    if(!is.null(label.pos)){
        lab.i <- match(label.pos[,1], i)
        text( log2(gs[l.names[lab.i]]), stats[2,lab.i] / stats[1,lab.i],
             lab.f(l.names[lab.i]), pos=label.pos[,2], font=3, cex=lab.cex)
    }
    ## plot(stats[1,], stats[2,], col=cols, pch=19,
    ##      xlab='Minimal in ancestor', ylab='Minimal in descendant')
    ## abline(0, 1)
    close.screen(screens)
##    legend('topright', names(class.col), lwd=2, col=class.col, cex=cex.leg)
    invisible(rbind(i, stats))
}


## provide histograms for values in x grouped by values
## in g (there is a better word for this, but I can't remember
## it at the moment).
hist.divided <- function(x, g, x.n, breaks, plot=FALSE){
    ## remove any values that are NA
    ## or infinite
    b <- !(is.na(x) | is.na(g) | is.infinite(x) | is.infinite(g))
    x <- x[b]
    g <- g[b]
    g.r <- range(g)
    g.breaks <- seq( g.r[1], g.r[2], length.out=x.n+1 )
    o <- order(g)
    x <- x[o]
    g <- g[o]
    h <- vector(mode='list', length=x.n-1)
    for(i in 2:length(g.breaks)){
        if(i < length(g.breaks)){
            b <- g >= g.breaks[i-1] & g < g.breaks[i]
        }else{
            b <- g >= g.breaks[i-1]
        }
        if(plot){
            h[[i-1]] <- hist(x[b], breaks=breaks, plot=plot, main=sprintf("%.1f -> %.1f", g.breaks[i-1], g.breaks[i]))
            abline( v=g.breaks[(i-1):i] )
        }else{
            h[[i-1]] <- hist(x[b], breaks=breaks, plot=plot)
        }
    }
    list(h=h, breaks=g.breaks)
}
    
mi2col <- function(v, cm){
    hsvScale(v, val=(cm + v)/(cm + max(v, na.rm=TRUE)))
}

## v is a set of values.
## counts gives the number of instances for each
## level of v
## counts and v should be the same length
un.table <- function(counts, v){
    counts <- as.integer(counts)
    counts[is.na(counts)] <- 0
    unlist(mapply( function(cnt, vv){ rep(vv, cnt) },
           counts, v))
}

## take a two dimensional table of counts
## and two vectors corresponding to the rows
## and column values and calculate variances
## for (A-B)
## rv and cv are the values corresponding
## to the rows and columns.. 
table.to.var <- function(tbl, rv, cv){
    if(nrow(tbl) != length(rv) || ncol(tbl) != length(cv))
        stop("Bad argument dimensions")
    vr <- sapply(1:length(rv), function(i){
        delta <- rv[i] - cv
        sum( (delta^2) * tbl[i,] ) / (sum(tbl[i,])-1)
    })
    list(v=vr, rv=rv, cv=cv, n=rowSums(tbl))
}

table.to.var.2 <- function(tbl, rv, cv, min.n){
    if(nrow(tbl) != length(rv) || ncol(tbl) != length(cv))
        stop("Bad argument dimensions")
    rs <- rowSums(tbl)
    v <- matrix(nrow=length(rs), ncol=3)
    colnames(v) <- c('v', 'rv', 'n')
    b.i <- 1
    row <- 1
    while(b.i <= nrow(tbl)){
        e.i <- b.i
        while( sum(rs[b.i:e.i]) < min.n && e.i < nrow(tbl) )
            e.i <- e.i + 1
        rvv <- t(matrix(rv[b.i:e.i], nrow=1+e.i-b.i, ncol=ncol(tbl)))
        delta <- rvv - cv
        counts <- t(tbl[e.i:b.i, ,drop=FALSE])
        vr <- sum( (delta^2) * counts ) / (sum(counts)-1)
        v[row, ] <- c(vr, mean(rv[b.i:e.i]), sum(counts))
        b.i <- e.i + 1
        row <- row + 1
    }
    v[!is.na(v[,1]), ]
}

table.to.ks <- function(tbl, rv, cv, min.n){
    if(nrow(tbl) != length(rv) || ncol(tbl) != length(cv))
        stop("Bad argument dimensions")
    rs <- rowSums(tbl)
    un.rs <- un.table(rs, rv) ## very long vector...
    ks <- matrix(nrow=length(rs), ncol=4)
    colnames(ks) <- c('d', 'p', 'rv', 'n')
    b.i <- 1
    row <- 1
    while(b.i <= nrow(tbl)){
        e.i <- b.i
        while( sum(rs[b.i:e.i]) < min.n && e.i < nrow(tbl) )
            e.i <- e.i + 1
        counts <- t(tbl[e.i:b.i, ,drop=FALSE])
        v <- un.table( counts, rep(cv, ncol(counts)) )
        k <- suppressWarnings(ks.test( un.rs, v ))
        ks[row, ] <- c(k$statistic, k$p.value, mean(rv[b.i:e.i]), sum(counts))
        b.i <- e.i + 1
        row <- row + 1
    }
    ks[!is.na(ks[,1]), ]
}

## this requires
## "~/R/c/convolve_ks.so"
blur.matrix <- function(m, kernel, na.v=0){
    if(!is.null(na.v)){
        m[is.na(m)] <- na.v
    }
    m <- apply(m, 1, function(x){ .Call("convolve_ks", x, kernel) })
    apply(m, 1, function(x){ .Call("convolve_ks", x, kernel) })
}

get.peaks <- function(x, y, w, min.v){
    d <- diff(y)
    p.i <- which( sapply(1:length(d), function(i){
        br <- (i-w):(i-1)
        ar <- (i+1):(i+w)
        br <- br[ br > 0 ]
        ar <- ar[ ar <= length(d) ]
        all( c(d[br] > 0, d[ar] < 0) )}))
    p.i <- p.i[ (y[p.i] + y[p.i+1]) > 2*min.v ]
    ## p.i <- which(sapply((w+1):(length(d)-w), function(i){
    ##         all(c( d[(i-w):(i-1)] > 0, d[(i+1):(i+w)] < 0))
    ## }))
    if(length(p.i) < 2)
        return(x[p.i])
    ## we expect each peak to be represented by a doublet.
    p.d <- diff(p.i)
    p.db <- p.d == 1
    peak.pos <- c()
    i <- 1
    while(i <= length(p.i)){
        if(i < length(p.i) && p.i[i+1] - p.i[i] == 1){
            peak.pos <- c(peak.pos, (x[p.i[i]] + x[p.i[i+1]]) / 2)
            i <- i + 2
            next
        }
        peak.pos <- c(peak.pos, x[p.i[i]])
        i <- i + 1
    }
    peak.pos
}

## takes as an argument an object returned by
## get.peaks()
## which should be a list with the peak positions for
## each slice of a matrix
## x is the x-value associated with each entry of pts
trace.lines <- function(pts, max.d){
    trace.ridge <- function(i, j, ridge, max.d){
        if(i > length(pts) || pts[[i]][j,'assigned'])
            return(ridge)
        ridge <- rbind(ridge, pts[[i]][j,1:2])
        pts[[i]][j,'assigned'] <<- TRUE
        if(i == length(pts))
            return(ridge)
        d <- abs( pts[[i]][j,'y'] - pts[[i+1]][,'y'] )
        d.i <- which.min(d)
        if(d[d.i] < max.d && !pts[[i+1]][d.i,'assigned'] )
            ridge <- trace.ridge(i+1, d.i, ridge, max.d)
        ridge
    }
    pts <- lapply(1:length(pts), function(i){ cbind('x'=pts[[i]][,1], 'y'=pts[[i]][,2], 'assigned'=0) })
    ridges <- list()
    for(i in 1:length(pts)){
        for(j in 1:nrow(pts[[i]])){
            if(!pts[[i]][j,'assigned']){
                ridge <- matrix(nrow=0, ncol=2)
                ridge <- trace.ridge(i, j, ridge, max.d)
                ridges <- c(ridges, list(ridge))
            }
        }
    }
    ridges
}
    

## width and step are in terms of number of measures
## which are ordered by x
windowed.density <- function(x, y, width, step, max.d=0.02 * diff(range(y))){
    o <- order(x);
    beg <- seq(1, length(x)-width, step)
    end <- beg + width
    wd <- lapply( 1:length(beg), function(i){
        j <- o[beg[i]:end[i]]
        x.q <- quantile(x[j])
        x.m <- mean(x[j])
        d <- density( y[j] )
        p <- get.peaks( d$x, d$y, 3 )
        list(i=j, x.q=x.q, x.m=x.m, d=d, p=p)
    })
    ## merge peaks into lines..
    peaks <- lapply( wd, function(x){ cbind(x=x$x.m, p=x$p, assigned=FALSE) })
    ## a little bit of recursion might work here..
    trace.ridge <- function(i, j, ridge, max.d){
        if(i > length(peaks) || peaks[[i]][j,'assigned'])
            return(ridge)
        ridge <- rbind(ridge, peaks[[i]][j,1:2])
        peaks[[i]][j,'assigned'] <<- TRUE
        if(i == length(peaks))
            return(ridge)
        d <- abs( peaks[[i]][j,'p'] - peaks[[i+1]][,'p'] )
        d.i <- which.min(d)
        if(d[d.i] < max.d && !peaks[[i+1]][d.i,'assigned'] )
            ridge <- trace.ridge(i+1, d.i, ridge, max.d)
        ridge
    }
    ridges <- list()
    for(i in 1:length(peaks)){
        for(j in 1:nrow(peaks[[i]])){
            if(!peaks[[i]][j,'assigned']){
                ridge <- matrix(nrow=0, ncol=2)
                ridge <- trace.ridge(i, j, ridge, max.d)
                ridges <- c(ridges, list(ridge))
            }
        }
    }
    list(wd=wd, ridges=ridges)
}

column.ks <- function(x, y, n, min.x=log2(76), min.y=min.x, rand.n=0){
    b <- x >= min.x & y >= min.y & !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
    x <- x[b]
    y <- y[b]
    x.r <- range(x)
    y.r <- range(y)
    breaks <- seq(x.r[1], x.r[2], length.out=n+1)
    x.lev <- cut(x, breaks=breaks, include.lowest=TRUE)
    ks <- tapply(y, x.lev, function(z){
        rand.ks <- NULL
        if(rand.n > 0)
            rand.ks <- sapply(1:rand.n, function(i){
                k <- ks.test(y, sample(y, length(z)))
                c('D'=k$statistic, 'p'=k$p.value) })
        list(k=ks.test(y, z), n=length(z), rand=rand.ks)})
    Dn <- sapply(ks, function(x){ c(x$n, x$k$statistic, x$k$p.value) })
    Dn.rd = NULL
    Dn.rp = NULL
    if(rand.n > 0){
        Dn.rd <- sapply(ks, function(x){ x$rand[1,] })
        Dn.rp <- sapply(ks, function(x){ x$rand[2,] })
    }
    mids <- breaks[-1] - diff(breaks)/2
    list(ks=ks, D=Dn[2,], n=Dn[1,], p=Dn[3,], breaks=breaks, mids=mids, rD=Dn.rd, rp=Dn.rp)
}

column.mean.dist <- function(x, y, n, min.x=log2(76), min.y=min.x, h.breaks.n=40){
    b <- x >= min.x & y >= min.y & !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
    x <- x[b]
    y <- y[b]
    x.r <- range(x)
    y.r <- range(y)
    breaks <- seq(x.r[1], x.r[2], length.out=n+1)
    x.lev <- cut(x, breaks=breaks, include.lowest=TRUE)
    all.h <- hist(y, breaks=seq(min(y), max(y), length.out=h.breaks.n+1), plot=FALSE)
    all.delta.h <- hist(y - x, breaks=h.breaks.n, plot=FALSE)
    y.p <- all.h$counts / sum(all.h$counts)
    y.l <- cut(y, breaks=all.h$breaks, included.lowest=TRUE)
    H.all <- entropy(all.h$counts, unit='log2')
    c.dist <- tapply(1:length(y), x.lev, function(i){
        n <- length(i)
        d <- sum( (x[i] - y[i])^2 )
        v <- var( (x[i] - y[i]) )
        delta <- sum( y[i] - x[i] )
        delta.h <- hist( y[i] - x[i], breaks=all.delta.h$breaks, plot=FALSE )
        v.x <- var( x[i] )
        v.y <- var( y[i] )
        cr <- cor(x[i], y[i])
        ks.d <- suppressWarnings( ks.test( y, y[i] )$statistic )
        h <- hist(y[i], breaks=all.h$breaks, plot=FALSE)
        kl.d <- KL.plugin( h$counts, all.h$counts, unit='log2' )
        H = entropy( h$counts, unit='log2' )
        y.lh <- sum( log(y.p[ y.l[i] ]) )
        exp.lh <- sum( ifelse(y.p > 0, length(i) * y.p * log(y.p), 0) )
        list(n=n, d=d, v=v, delta=delta, delta.h=delta.h, v.x=v.x, v.y=v.y, ks.d=ks.d, kl.d=kl.d, h=h, H=H,
             lh=y.lh, e.lh=exp.lh, cor=cr)
    })
    mids <- breaks[-1] - diff(breaks)[1]/2
    col.n <- sapply( c.dist, function(x){ x$n })
    col.d <- sapply( c.dist, function(x){ x$d }) / (col.n)
    col.v <- sapply( c.dist, function(x){ x$v })
    col.delta <- sapply( c.dist, function(x){ x$delta })
    col.delta.h <- lapply( c.dist, function(x){ x$delta.h })
    col.xv <- sapply( c.dist, function(x){ x$v.x })
    col.yv <- sapply( c.dist, function(x){ x$v.y })
    col.h <-  lapply( c.dist, function(x){ x$h })
    col.ks.d <- sapply( c.dist, function(x){ x$ks.d })
    col.kl.d <- sapply( c.dist, function(x){ x$kl.d })
    col.lh <- sapply( c.dist, function(x){ x$lh })
    col.e.lh <- sapply( c.dist, function(x){ x$e.lh })
    col.cor <- sapply( c.dist, function(x){ x$cor })
    H = sapply( c.dist, function(x){ x$H })
    list(breaks=breaks, d=col.d, n=col.n, v=col.v, delta=col.delta, delta.h=col.delta.h, all.delta.h=all.delta.h,
         mids=mids, xv=col.xv, yv=col.yv,
         ks.d=col.ks.d, kl.d=col.kl.d, h=col.h, all.h=all.h, H=H, H.all=H.all,
         lh=col.lh, e.lh=col.e.lh, cor=col.cor)
}

column.mean.dist.2 <- function(x, y, n, min.x=log2(76), min.y=min.x){
    b <- x >= min.x & y >= min.y & !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
    x <- x[b]
    y <- y[b]
    x.r <- range(x)
    y.r <- range(y)
    o <- order(x)
    o.beg <- as.integer(seq(1, length(o), length.out=n+1))
    o.end <- o.beg[-1]
    o.beg <- o.beg[-length(o.beg)]
    all.h <- hist(y, breaks=10, plot=FALSE)
    c.dist <- lapply(1:length(o.beg), function(ind){
        i <- o[o.beg[ind]:o.end[ind]]  ## will have overlap of one
        n <- length(i)
        d <- sum( (x[i] - y[i])^2 )
        v <- var( (x[i] - y[i]) )
        v.x <- var( x[i] )
        v.y <- var( y[i] )
        m.x <- mean( x[i] )
        m.y <- mean( y[i] )
        r.x <- range( x[i] )
        r.y <- range( y[i] )
        x.y.cor <- cor(x[i], y[i])
        h <- hist(y[i], breaks=all.h$breaks, plot=FALSE)
        kl.d <- KL.plugin( h$density, all.h$density )
        chi.d <- chi2.plugin( h$density, all.h$density )
        ks.d <- suppressWarnings( ks.test( y[i], y )$statistic )
        list(n=n, d=d, v=v, v.x=v.x, v.y=v.y, m.x=m.x, m.y=m.y, r.x=r.x, r.y=r.y, xy.c=x.y.cor,
             h=h, kl.d=kl.d, chi.d=chi.d, ks.d=ks.d)
    })
    x.m <- sapply( c.dist, function(x){ x$m.x })
    y.m <- sapply( c.dist, function(x){ x$m.y })
    x.r <- sapply( c.dist, function(x){ x$r.x })
    y.r <- sapply( c.dist, function(x){ x$r.y })
    xy.c <- sapply( c.dist, function(x){ x$xy.c })
    kl.d <- sapply( c.dist, function(x){ x$kl.d })
    chi.d <- sapply( c.dist, function(x){ x$chi.d })
    ks.d <-  sapply( c.dist, function(x){ x$ks.d })
    col.n <- sapply( c.dist, function(x){ x$n })
    col.d <- sapply( c.dist, function(x){ x$d }) / (col.n)
    col.v <- sapply( c.dist, function(x){ x$v })
    col.xv <- sapply( c.dist, function(x){ x$v.x })
    col.yv <- sapply( c.dist, function(x){ x$v.y })
    col.h <-  lapply( c.dist, function(x){ x$h })
    list(d=col.d, n=col.n, v=col.v, x.m=x.m, y.m=y.m, x.r=x.r, y.r=y.r, xy.c=xy.c,
         xv=col.xv, yv=col.yv, h=col.h, kl.d=kl.d, chi.d=chi.d, ks.d=ks.d)
}

column.H <- function(x, y, n, min.x=log2(76), min.y=min.x){ ## , rand.n=0){
    b <- x >= min.x & y >= min.y & !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
    x <- x[b]
    y <- y[b]
    x.r <- range(x)
    y.r <- range(y)
    x.breaks <- seq(x.r[1], x.r[2], length.out=n+1)
    x.lev <- cut(x, breaks=x.breaks, include.lowest=TRUE)
    y.breaks <- seq(y.r[1], y.r[2], length.out=n+1)
    y.lev <- cut(y, breaks=y.breaks, include.lowest=TRUE)
    y.h <- entropy( table(y.lev) )
    c.h <- tapply(y, x.lev, function(z){
        rand.mi <- NULL
#        if(rand.n > 0)
#            rand.mi <- sapply(1:rand.n, function(i){
#                k <- ks.test(y, sample(y, length(z)))
#                c('D'=k$statistic, 'p'=k$p.value) })
        z.h <- entropy( table(cut(z, breaks=y.breaks, include.lowest=TRUE)) )
        list(z.h=z.h, n=length(z))})
    zh <- sapply(c.h, function(x){ c(x$n, x$z.h) })
    mids <- x.breaks[-1] - diff(x.breaks)/2
    list(zh=zh[2,], n=zh[1,], breaks=x.breaks, mids=mids, yh=y.h)
}

## how many steps up the tree from node A to node B?
nodes.sep <- function(A, B, tree=ex.align.2.k2.nj){
    ## tree$edge[,1] = parent node
    ## tree$edge[,2] = child node
    ascend.tree <- function(child, target, dist, edge=tree$edge){
        dist <- dist + 1
        p.i <- which( edge[,2] == child )
        if(length(p.i) == 0)
            return(0)
        p <- edge[p.i, 1]
        if(p == child)
            return(0)
        if(p == target)
            return(dist)
        ascend.tree( p, target, dist, edge )
    }
    ascend.tree(A, B, 0)
}

## x is a set of values from which we will pick n
## values with replacement and then calculate the likelihood
## of the sample based on discretizing the values into n.bins
## rep is the number of times that the function will sample n values
## from x
exp.samp.likelihood <- function(x, n, n.bins, rep=1, r=range(x)){
    ## these two lines are essentially copied from the discretize
    ## function of the entropy package.
    breaks <- seq(from=r[1], to=r[2], length.out = n.bins + 1)
    x.l <- cut(x, breaks=breaks, include.lowest=TRUE)
    x.counts <- table(x.l)
    x.p <- x.counts / sum(x.counts)
    ## we will simulate reps by sampling a single time and subdividing
    ## the resulting values into ranges
    ## we could sample from the counts as well by specifying the probabilities
    x.s <- sample(x, size=n*rep, replace=TRUE)
    x.s.l <- cut(x.s, breaks=breaks, include.lowest=TRUE)
    s.beg <- seq(1, n*rep, n)
    s.end <- s.beg + n - 1
    sampled.likelihoods <- sapply(1:length(s.beg), function(i){
        sum( log( x.p[ x.s.l[s.beg[i]:s.end[i]] ] ))
    })
    exp.likelihood <- sum( ifelse(x.p > 0, n*x.p * log(x.p), 0) )
    list(exp=exp.likelihood, obs=sampled.likelihoods)
}

## sizes should be discretized
common.minimisation.by.ancestral.size <- function(anc, d1, d2, max.s=10*log2(100), min.s=10*log2(76),
                                                  max.s2=10*log2(1024), inf=int.s.inf.2, disc.n=1){
    int.s <- cbind( inf[[anc]]$state[,1], inf[[d1]]$state[,1], inf[[d2]]$state[,1] )
    colnames(int.s) <- c('a', 'd1', 'd2')
    row.min <- apply(int.s, 1, min)
    b <- !is.na(row.min) & row.min >= min.s
    int.s <- int.s[b,]
    anc.r <- range(int.s[,'a'])
    anc.breaks <- seq(anc.r[1], anc.r[2], disc.n)
    anc.l <- cut(int.s[,'a'], breaks=anc.breaks, include.lowest=TRUE)
    tbl <- tapply( 1:nrow(int.s), anc.l, function(i){
        b1 <- int.s[i,'d1'] <= max.s
        b2 <- int.s[i,'d2'] <= max.s
        c('as'=mean(int.s[i,'a']), 'p1'=sum(b1)/length(i), 'p2'=sum(b2)/length(i),
          'p1.2'=sum(b1 & b2) / length(i), 'n1'=sum(b1), 'n2'=sum(b2), 'n1.2'=sum(b1 & b2),
          'n.exp'=as.double(sum(b1)) * sum(b2) / length(i), 'n'=length(i),
          p=phyper(sum(b1 & b2)-1, sum(b1), sum(!b1), sum(b2), lower.tail=FALSE))
    })
    df <- as.data.frame(do.call(rbind, tbl))
    observed <- sum( df[,'n1.2'] )
    expected <- sum( df[,'n.exp'] )
    comb.p <- pchisq( -2 * sum(log(df[,'p'])), df=2*nrow(df), lower.tail=FALSE)
    b2 <- df$as <= max.s2
    comb.p2 <- pchisq( -2 * sum(log(df[b2,'p'])), df=2*sum(b2), lower.tail=FALSE)
    list(df=df, comb=c(obs=observed, exp=expected, p=comb.p, p2=comb.p2))
}


## This is the simple version which simply uses a set range of ancestral sizes
## and calculates a single set of statistics
common.minimisation <- function(anc, d1, d2, max.s=10*log2(100), min.s=10*log2(76), inf=int.s.inf.2,
                                                  anc.range=c(max.s, 10*log2(1024))){
    int.s <- cbind( inf[[anc]]$state[,1], inf[[d1]]$state[,1], inf[[d2]]$state[,1] )
    colnames(int.s) <- c('a', 'd1', 'd2')
    row.min <- apply(int.s, 1, min)
    b <- !is.na(row.min) & row.min >= min.s
    b <- b & int.s[,'a'] >= anc.range[1] & int.s[,'a'] <= anc.range[2]
    int.s <- int.s[b,]
##
    b1 <- int.s[,'d1'] <= max.s
    b2 <- int.s[,'d2'] <= max.s
    n <- nrow(int.s)
    c('p1'=sum(b1)/n, 'p2'=sum(b2)/n,
      'p1.2'=sum(b1 & b2) / n, 'n'=n, 'oe'=sum(b1 & b2) / (sum(b1) * sum(b2) / n),
      p=phyper(sum(b1 & b2)-1, sum(b1), sum(!b1), sum(b2), lower.tail=FALSE))
}

## This is the simple version which simply uses a set range of ancestral sizes
## and calculates a single set of statistics
common.long <- function(anc, d1, d2, min.s=10*log2(76), min.s2=10*log2(256), inf=int.s.inf.2,
                                                  anc.range=c(min.s2, 10*log2(1024))){
    int.s <- cbind( inf[[anc]]$state[,1], inf[[d1]]$state[,1], inf[[d2]]$state[,1] )
    colnames(int.s) <- c('a', 'd1', 'd2')
    row.min <- apply(int.s, 1, min)
    b <- !is.na(row.min) & row.min >= min.s
    b <- b & int.s[,'a'] >= anc.range[1] & int.s[,'a'] <= anc.range[2]
    int.s <- int.s[b,]
##
    b1 <- int.s[,'d1'] >= min.s2
    b2 <- int.s[,'d2'] >= min.s2
    n <- nrow(int.s)
    c('p1'=sum(b1)/n, 'p2'=sum(b2)/n,
      'p1.2'=sum(b1 & b2) / n, 'n'=n, 'oe'=sum(b1 & b2) / (sum(b1) * sum(b2) / n),
      p=phyper(sum(b1 & b2)-1, sum(b1), sum(!b1), sum(b2), lower.tail=FALSE))
}


common.long.by.ancestral.size <- function(anc, d1, d2, min.s1=10*log2(76), min.s2=10*log2(1024), inf=int.s.inf.2, disc.n=1){
    int.s <- cbind( inf[[anc]]$state[,1], inf[[d1]]$state[,1], inf[[d2]]$state[,1] )
    colnames(int.s) <- c('a', 'd1', 'd2')
    row.min <- apply(int.s, 1, min)
    b <- !is.na(row.min) & row.min >= min.s1
    int.s <- int.s[b,]
    
    anc.r <- range(int.s[,'a'])
    anc.breaks <- seq(anc.r[1], anc.r[2], disc.n)
    anc.l <- cut(int.s[,'a'], breaks=anc.breaks, include.lowest=TRUE)

    tbl <- tapply( 1:nrow(int.s), anc.l, function(i){
        b1 <- int.s[i,'d1'] >= min.s2
        b2 <- int.s[i,'d2'] >= min.s2
        c('as'=mean(int.s[i[1],'a']), 'p1'=sum(b1)/length(i), 'p2'=sum(b2)/length(i),
          'p1.2'=sum(b1 & b2) / length(i), 'n1'=sum(b1), 'n2'=sum(b2), 'n1.2'=sum(b1 & b2),
          'n.exp'=as.double(sum(b1) * sum(b2) / length(i)), 'n'=length(i),
          p=phyper(sum(b1 & b2)-1, sum(b1), sum(!b1), sum(b2), lower.tail=FALSE) )
    })
    df <- as.data.frame(do.call(rbind, tbl))
    observed <- sum( df[,'n1.2'] )
    expected <- sum( df[,'n.exp'] )
    comb.p <- pchisq( -2 * sum(log(df[,'p'])), df=2*nrow(df), lower.tail=FALSE)
    ## p2 is for those that are at least as long as min.s2 in the ancestral size
    b2 <- df$as >= min.s2
    comb.p2 <- pchisq( -2 * sum(log(df[b2,'p'])), df=2*sum(b2), lower.tail=FALSE)
    list(df=df, comb=c(obs=observed, exp=expected, p=comb.p, p2=comb.p2))
}

common.constant.ancestral.size <- function(anc, d1, d2, min.s1=10*log2(76), max.d=10, inf=int.s.inf.2, disc.n=1){
    int.s <- cbind( inf[[anc]]$state[,1], inf[[d1]]$state[,1], inf[[d2]]$state[,1] )
    colnames(int.s) <- c('a', 'd1', 'd2')
    row.min <- apply(int.s, 1, min)
    b <- !is.na(row.min) & row.min >= min.s1
    int.s <- int.s[b,]
    
    anc.r <- range(int.s[,'a'])
    anc.breaks <- seq(anc.r[1], anc.r[2], disc.n)
    anc.l <- cut(int.s[,'a'], breaks=anc.breaks, include.lowest=TRUE)

    tbl <- tapply( 1:nrow(int.s), anc.l, function(i){
        b1 <- abs(int.s[i,'d1'] - int.s[i,'a']) <= max.d
        b2 <- abs(int.s[i,'d2'] - int.s[i,'a']) <= max.d
        c('as'=mean(int.s[i[1],'a']), 'p1'=sum(b1)/length(i), 'p2'=sum(b2)/length(i),
          'p1.2'=sum(b1 & b2) / length(i), 'n'=length(i),
          p=phyper(sum(b1 & b2)-1, sum(b1), sum(!b1), sum(b2), lower.tail=FALSE) )
    })
    as.data.frame(do.call(rbind, tbl))
}

## return number of introns within range (inclusive) in the root
## the number of the union (i.e. within range in any species)
## some form of quantiles of the above
## note that min.s and max.s should be given as log2(size) units
clade.shared.length <- function(root, min.s, max.s, tree=ex.align.2.k2.nj, int.s=orth$l, inf=int.s.inf.2){
    root.s <- int.s.inf.2[[root]]$state[,1] / 10
    root.b <- root.s >= min.s & root.s <= max.s
    root.n <- sum(root.b, na.rm=TRUE)
    leaves <- collect.leaves(root, tree)
    leaves.s <- log2( int.s[, leaves$names] )
    leaves.b <- apply( leaves.s, 2, function(x){ x >= min.s & x <= max.s } )
    union <- rep(FALSE, nrow(leaves.b))
    for(i in 1:ncol(leaves.b))
        union <- union | leaves.b[,i]
    sp.count <- rowSums(leaves.b, na.rm=TRUE)
    qnt <- quantile( sp.count[union], probs=seq(0,1,0.1), na.rm=TRUE )
    tbl <- table( c(0:ncol(leaves.b), sp.count) ) - 1
    pair.n <- apply(leaves.b, 2, function(x){
        colSums(x & leaves.b, na.rm=TRUE)
    })
    pair.nu <- apply(leaves.b, 2, function(x){
        colSums(x | leaves.b, na.rm=TRUE)
    })
    pair.n.exp <- apply(leaves.b, 2, function(x){
        apply( leaves.b, 2, function(y){
            (sum(x, na.rm=TRUE) * sum(y, na.rm=TRUE)) / sum( !is.na(x) & !is.na(y) )
        })
    })
    list(u=sum(union, na.rm=TRUE),q=qnt, tbl=tbl, sp.n=ncol(leaves.b), leaves=leaves,
         leaf.n=colSums(leaves.b, na.rm=TRUE), pair.n=pair.n, pair.nu=pair.nu, pair.n.exp=pair.n.exp,
         root.n=root.n,
         b=cbind(root=root.b, leaves.b))
}

recycle <- function(v, i){
    v[ 1 + (i-1) %% length(v) ]
}

## the following is copied from the entropy package function. It is different
## in that it returns the breaks used as well as the frequencies
discretize2d.2 <- function(x1, x2, numBins1, numBins2, r1 = range(x1), r2 = range(x2))
{
    b1 = seq(from = r1[1], to = r1[2], length.out = numBins1 + 
        1)
    b2 = seq(from = r2[1], to = r2[2], length.out = numBins2 + 
        1)
    y2d = table(cut(x1, breaks = b1, include.lowest = TRUE), 
        cut(x2, breaks = b2, include.lowest = TRUE))
    list(f=y2d, b1=b1, b2=b2)
}

## the following requires the entropy package.
## equations are taken from the entropy package manual.
KL.matrix <- function(x1, x2, nbins1=20, nbins2=nbins1, unit="log", h.tform=eval){
    require(entropy)
    b <- !(is.na(x1) | is.na(x2))
    f.2 <- discretize2d.2( x1[b], x2[b], numBins1=nbins1, numBins2=nbins2 )
    f.2d <- with(f.2, h.tform(f))
    f.2d <- f.2d/sum(f.2d)
    xf <- rowSums(f.2d)
    yf <- colSums(f.2d)
    freqs.null <- xf %o% yf
    LR <- ifelse( f.2d > 0, log( f.2d / freqs.null ), 0 )
    chi.sq.d <- (f.2d - freqs.null)^2 / freqs.null
    KL <- f.2d * LR
    mi <- mi.plugin( f.2d, unit=unit )
    ## note: sum(KL) should equal mi
    ## included here to confirm my calculations
    list("f2"=f.2d, "null"=freqs.null, "LR"=LR, "KL"=KL, "chi.d"=chi.sq.d, "xf"=xf, "yf"=yf, "mi"=mi, b1=f.2$b1, b2=f.2$b2)
}

## obtain KL.matrix between:
## root and leaves
## leaves against each other..
clade.KL <- function(root, tree=ex.align.2.k2.nj, int.s=orth$l, inf=int.s.inf.2,
                     nbins1=20, nbins2=nbins1, min.s=log2(76), h.tform=eval){
    ## intron sizes as log2 length
    root.s <- int.s.inf.2[[root]]$state[,1] / 10
    leaves <- collect.leaves(root, tree)
    leaves.s <- log2(int.s[,leaves$names])
    l.b <- apply( leaves.s, 2, function(x){ !is.na(x) & x >= min.s })
    root.kl <- lapply(leaves$names, function(sp1){
        b <- l.b[,sp1] & root.s >= min.s
        KL.matrix( root.s[b], leaves.s[b,sp1], h.tform=h.tform, nbins1=nbins1, nbins2=nbins2 )
    })
    names(root.kl) <- leaves$names
    leaves.kl <- lapply(leaves$names, function(sp1){
        kl <- lapply(leaves$names, function(sp2){
            b <- l.b[,sp1] & l.b[,sp2]
            KL.matrix(leaves.s[b,sp1], leaves.s[b,sp2], h.tform=h.tform, nbins1=nbins1, nbins2=nbins2)
        })
        names(kl) <- leaves$names
        kl
    })
    names(leaves.kl) <- leaves$names
    list('root'=root.kl, 'leaves'=leaves.kl)
}
