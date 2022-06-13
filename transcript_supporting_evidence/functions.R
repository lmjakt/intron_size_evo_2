space.text.y <- function(y, words, margin=1.5, y.range=par("usr")[3:4], ...){
                                        #    words <- words[ !is.na(y) ]
                                        #    y <- y[ !is.na(y) ]
    y.r <- rank(y)
    y <- sort(y)
    h <- max(strheight(words, ...), na.rm=TRUE) * margin
    while(min(diff(y)) < h){ ## potentially infinite loop
        ## do the stupid way, iterating through each time
        b.i <- 0
        e.i <- 0
        for(i in 2:length(y)){
            if(!b.i && y[i] - y[i-1] < h)
                b.i <- i-1
            if(b.i && (is.na(y[i]) || y[i] - y[i-1] > h)){
                e.i <- i-1
                break
            }
        }
        if(!e.i) e.i <- i
        ## space out from the region starting at the middle
        ## the middle position is biassed towards the top
        c.i <- b.i + as.integer( (e.i - b.i) / 2 )
        ## first go downwards..
        i <- c.i
        kludge <- 1.001 ## or we can end up in infinite loop
        while(i > b.i){
            y[i-1] <- y[i-1] - (h - (y[i]-y[i-1]))*kludge
            i <- i - 1
        }
        if(y[b.i] < (y.range[1] + h/2)){
            e.j <- b.i + min( which( diff(y[b.i:length(y)]) / h) - kludge) > 0.01
            e.j <- ifelse(length(e.j) > 1, e.j, b.i)
            y[b.i:e.j] <- y[b.i:e.j] + kludge * ((y.range[1] + h/2) - y[b.i])
        }
        i <- c.i
        while(i < e.i){
            y[i+1] <- y[i+1] + (h - (y[i+1]-y[i]))*kludge
            i <- i + 1
        }
        if(y[e.i] > (y.range[2] - h/2)){
            b.j <- max( which( ((diff(y[1:e.i]) / h) - kludge) > 0.01 ))
            b.j <- ifelse( length(b.j) > 0, b.j, e.i )
            y[b.j:e.i] <- y[b.j:e.i] + kludge * ((y.range[2] - h/2) - y[e.i])
        }
    }
    return(y[y.r])
}
