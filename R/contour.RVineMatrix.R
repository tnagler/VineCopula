#' @method contour RVineMatrix
#' @rdname plot.RVineMatrix
contour.RVineMatrix <- function(x, tree = "ALL", xylim = NULL, cex.nums = 1,
                                data = NULL, ...) {

    ## check input
    d <- nrow(x$Matrix)
    if (all(tree == "ALL"))
        tree <- seq.int(d-1)
    n.tree <- length(tree)
    if (!is.null(list(...)$type))
        stop("Only contour plots allowed. Don't use the type argument!")

    ## set up for plotting windows (restore settings on exit)
    usr <- par(mfrow = c(n.tree, d - min(tree)), mar = rep(0, 4))
    on.exit(par(usr))

    ## calculate pseudo-observations (if necessary)
    psobs <- if (!is.null(data)) vine_psobs(data, x) else NULL

    ## default style --------------------------------------------------
    # headings: create blue color scale
    TUMblue <- rgb(0, 103/255, 198/255)
    tint.seq <- seq(0.5, 0.5, length.out = d - 1)
    clrs <- rev(sapply(tint.seq, function(x) tint(TUMblue, x, 0.7)))

    # contours: set limits for plots
    if (!is.null(list(...)$margins)) {
        margins <- list(...)$margins
        if (!(margins %in% c("norm", "unif", "exp", "flexp")))
            contour(BiCop(0), margins = c(0, 10))
    } else {
        margins <- "norm"
    }
    if (is.null(xylim))
        xylim <- switch(margins,
                        "norm"  = c(-3, 3),
                        "unif"  = c(0, 1 - 1e-2),
                        "exp"   = c(0, 10),
                        "flexp" = c(-10, 0))
    xlim <- ylim <- xylim

    # contours: adjust limits for headings
    offs <- 0.25
    mult <- 1.35
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
                        if (is.null(psobs)) {
                            pcfit <- BiCop(family = x$family[d - i + 1, j],
                                           par    = x$par[d - i + 1, j],
                                           par2   = x$par2[d - i + 1, j],
                                           check.pars = FALSE)
                        } else {
                            pcfit <- kdecopula::kdecop(psobs[[i]][[j]])
                        }

                        # set up list of contour arguments
                        args <- list(x = pcfit,
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

                        # call contour.BiCop with ... arguments
                        do.call(contour, modifyList(args, list(...)))

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
                        cx1 <- 0.75 * diff(xlim) / strwidth(maxnums)
                        ty <- ylim[2] - diff(ylim)/mult*offs
                        cx2 <- 0.75 * (ylim[2] - ty) / strheight(maxnums)
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

vine_psobs <- function (uev, object) {
    uev <- as.matrix(uev)
    if (ncol(uev) == 1)
        uev <- matrix(uev, 1, nrow(uev))
    if (any(uev > 1) || any(uev < 0))
        stop("Data has be in the interval [0,1].")
    n <- ncol(uev)
    N <- nrow(uev)
    if (ncol(uev) != ncol(object$Matrix))
        stop("Dimensions of 'data' and 'object' do not match.")
    if (!is(object,"RVineMatrix"))
        stop("'object' has to be an RVineMatrix object")

    o <- diag(object$Matrix)
    oldobject <- object
    if (any(o != length(o):1)) {
        object <- normalizeRVineMatrix(object)
        uev <- matrix(uev[, o[length(o):1]], N, n)
    }

    ## initialize objects
    CondDistr <- neededCondDistr(object$Matrix)
    val <- array(1, dim = c(n, n, N))
    out <- lapply(1:(n - 1), list)
    V <- list()
    V$direct <- array(NA, dim = c(n, n, N))
    V$indirect <- array(NA, dim = c(n, n, N))
    V$direct[n, , ] <- t(uev[, n:1])

    for (i in (n - 1):1) {
        for (k in n:(i + 1)) {
            ## extract data for current tree
            m <- object$MaxMat[k, i]
            zr1 <- V$direct[k, i, ]
            if (m == object$Matrix[k, i]) {
                zr2 <- V$direct[k, (n - m + 1), ]
            } else {
                zr2 <- V$indirect[k, (n - m + 1), ]
            }

            ## store data
            out[[n - k + 1]][[i]] <- cbind(zr2, zr1)

            ## calculate pseudo-observations for next tree
            if (CondDistr$direct[k - 1, i])
                V$direct[k - 1, i, ] <- BiCopHfunc1(zr2, zr1,
                                                    object$family[k, i],
                                                    object$par[k, i],
                                                    object$par2[k, i],
                                                    check.pars = FALSE)
            if (CondDistr$indirect[k - 1, i])
                V$indirect[k - 1, i, ] <- BiCopHfunc2(zr2, zr1,
                                                      object$family[k, i],
                                                      object$par[k, i],
                                                      object$par2[k, i],
                                                      check.pars = FALSE)
        }
    }

    ## return list of pseudo-observations
    out
}

