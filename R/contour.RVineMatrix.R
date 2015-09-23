## -----------------------------------------------------------------------------
## contour generic for RVineMatrix objects
contour.RVineMatrix <- function(x, tree = "ALL", xylim = NULL, cex.nums = 1, ...) {
    
    ## check input
    d <- nrow(x$Matrix)
    if (all(tree == "ALL"))
        tree <- seq.int(d-1)
    n.tree <- length(tree)
    if (!is.null(list(...)$type)) 
        stop("Only contour plots allowed. Don't use the type argument!")
    
    ## set up for plotting windows (restore settings on exit)
    usr <- par(mfrow = c(n.tree, d - min(tree)), mar = rep(0, 4))  # dimensions of contour matrix
    on.exit(par(usr))
    
    ## default style --------------------------------------------------
    # headings: create blue color scale
    TUMblue <- rgb(0, 103/255, 198/255)
    tint.seq <- seq(0.5, 0.5, length.out = d - 1)
    clrs <- rev(sapply(tint.seq, function(x) tint(TUMblue, x, 0.7)))
    
    # contours: set limits for plots
    if (!is.null(list(...)$margins)) {
        margins <- list(...)$margins
        if (!(margins %in% c("norm", "unif")))
            stop("margins not supported")
    } else {
        margins <- "norm"
    }
    if (is.null(xylim))
        xylim <- switch(margins,
                        "norm" = c(-3, 3),
                        "unif" = c(1e-1, 1 - 1e-1))
    xlim <- ylim <- xylim
    
    # contours: adjust limits for headings
    offs <- 0.25 
    mult <- 1.5
    ylim[2] <- ylim[2] + offs*diff(ylim)
    
    
    ## run through trees -----------------------------------------------
    # initialize check variables 
    cnt <- 0
    k <- d
    e <- numeric(0)
    class(e) <- "try-error"
    
    while ("try-error" %in% class(e)) {
        e <- try({
            maxnums <- get_num(1, tree = max(tree), RVM = x)
            for (i in tree) {
                for (j in 1:(d - min(tree))) {
                    if (d - i >= j) {
                        # set up list of contour arguments
                        args <- list(x = BiCop(family=x$family[d-i+1,j],
                                               par=x$par[d-i+1,j],
                                               par2=x$par2[d-i+1,j]),
                                     drawlabels = FALSE,
                                     xlab = "",
                                     ylab = "",
                                     xlim = xlim, 
                                     ylim = ylim,
                                     xaxt = "n",
                                     yaxt = "n",
                                     add  = TRUE)
                        
                        # create empty plot
                        plot.new()
                        plot.window(xlim = xlim, ylim = ylim, 
                                    xaxs = "i",  yaxs = "i")
                        
                        # call plot.BiCop with ... arguments 
                        do.call(plot, modifyList(args, list(...)))
                        
                        # draw area for headings
                        abline(h = ylim[2] - diff(ylim)/mult*offs)
                        ci <- min(length(clrs) + 1 - i, 10)
                        polygon(x = c(xlim[1] - diff(xlim),
                                      xlim[1] - diff(xlim), 
                                      xlim[2] + diff(xlim),
                                      xlim[2] + diff(xlim)),
                                y = c(ylim[2] + diff(ylim)/mult*offs,
                                      ylim[2] - diff(ylim)/mult*offs, 
                                      ylim[2] - diff(ylim)/mult*offs, 
                                      ylim[2] + diff(ylim)/mult*offs),
                                col = clrs[ci])
                        
                        # add separating lines
                        abline(v = xlim)
                        abline(h = ylim)
                        
                        # add pair-copula ID
                        cx1 <- 0.95 * diff(xlim) / strwidth(maxnums)
                        cx1 <- cx1
                        ty <- ylim[2] - diff(ylim)/mult*offs
                        cx2 <- 0.95 * (ylim[2] - ty) / strheight(maxnums)
                        cx2 <- cx2
                        cx <- min(cx1, cx2)
                        text(x = sum(xlim)/2,
                             y = ty + 0.225 / cex.nums * (ylim[2] - ty), 
                             cex    = cex.nums * cx,
                             labels = get_num(j, tree = i, RVM = x), 
                             pos    = 3,
                             offset = 0)
                    } else {
                        plot.new()
                    }
                }
            }
        }
        , silent = TRUE)
        
        ## adjust to figure margins if necessary
        if (length(tree) < 1)
            stop("Error in plot.new() : figure margins too large")
        if ("try-error" %in% class(e)) {
            cnt <- cnt + 1
            tree <- tree[-which(tree == max(tree))]
            par(mfrow = c(n.tree - cnt, d - min(tree)))
        }        
    }
    
    ## message for the user if not all trees could be plotted -----------
    if (length(tree) != n.tree) {
        nmbr.msg <- as.character(tree[1])
        if (length(tree) > 2) {
            for (i in tree[-c(1, length(tree))]) {
                nmbr.msg <- paste(nmbr.msg, i, sep=", ")
            }
        }
        if (length(tree) > 1) {
            s.msg <- "s "
            nmbr.msg <- paste(nmbr.msg,
                              "and",
                              tree[length(tree)],
                              "were plotted. ")
        } else {
            s.msg <- " "
            nmbr.msg <- paste(nmbr.msg, "was plotted. ", sep=" ")
        }    
        msg.space <- "There is not enough space."
        msg.tree <- paste("Only Tree",
                          s.msg,
                          nmbr.msg,
                          "Use the 'tree' argument or enlarge figure margins",
                          " to see the others.",
                          sep = "")
        message(paste(msg.space, msg.tree))
    }
}

tint <- function(x, fac, alpha = 1) {
    x <- c(col2rgb(x))
    x <- (x + (255 - x) * fac) / 255
    rgb(x[1], x[2], x[3], alpha)
}
