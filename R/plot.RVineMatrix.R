#' Plotting `RVineMatrix` objects.
#'
#' There are two plotting generics for `RVineMatrix` objects.
#' `plot.RVineMatrix` plots one or all trees of a given R-vine copula
#' model. Edges can be labeled with information about the corresponding
#' pair-copula. `contour.RVineMatrix` produces a matrix of contour plots
#' (using [VineCopula::plot.BiCop()]).
#'
#' If you want the contour boxes to be perfect squares, the plot height should
#' be `1.25/length(tree)*(d - min(tree))` times the plot width.
#'
#' @method plot RVineMatrix
#'
#' @aliases plot.RVineMatrix contour.RVineMatrix
#'
#' @param x `RVineMatrix` object.
#' @param tree `"ALL"` or integer vector; specifies which trees are
#' plotted.
#' @param type integer; specifies how to make use of variable names: \cr
#' `0` = variable names are ignored, \cr `1` = variable names are
#' used to annotate vertices, \cr `2` = uses numbers in plot and adds a
#' legend for variable names.
#' @param edge.labels character; either a vector of edge labels or one of the
#' following: \cr `"family"` = pair-copula family abbreviation (see
#' [VineCopula::BiCopName()]), \cr `"par"` =
#' pair-copula parameters, \cr `"tau"` = pair-copula Kendall's tau (by
#' conversion of parameters) \cr `"family-par"` = pair-copula family and
#' parameters \cr `"family-tau"` = pair-copula family and Kendall's tau.
#' @param legend.pos the `x` argument for
#' [graphics::legend()].
#' @param interactive logical; if TRUE, the user is asked to adjust the
#' positioning of vertices with his mouse.
#' @param xylim numeric vector of length 2; sets `xlim` and `ylim`
#' for the contours
#' @param cex.nums numeric; expansion factor for font of the numbers.
#' @param data a data matrix for creating kernel density contours of each pair.
#' @param \dots Arguments passed to
#' [network::plot.network()] or
#' [VineCopula::plot.BiCop()] respectively.
#'
#' @author Thomas Nagler, Nicole Barthel
#'
#' @seealso [VineCopula::RVineMatrix()],
#' [network::plot.network()],
#' [VineCopula::plot.BiCop()],
#' [VineCopula::BiCopName()],
#' [graphics::legend()]
#'
#' @keywords plot
#'
#' @examples
#'
#' ## build vine model
#' strucmat <- matrix(c(3,   1, 2, 0, 2, 1, 0, 0, 1), 3, 3)
#' fammat   <- matrix(c(0,   1, 6, 0, 0, 3, 0, 0, 0), 3, 3)
#' parmat   <- matrix(c(0, 0.3, 3, 0, 0, 1, 0, 0, 0), 3, 3)
#' par2mat  <- matrix(c(0,   0, 0, 0, 0, 0, 0, 0, 0), 3, 3)
#' RVM  <- RVineMatrix(strucmat, fammat, parmat, par2mat)
#'
#' # plot trees
#' \dontrun{plot(RVM)}
#'
#' # show contour plots
#' contour(RVM)
#'
plot.RVineMatrix <- function(x, tree = "ALL", type = 0, edge.labels = NULL,
                             legend.pos = "bottomleft", interactive = FALSE, ...) {
    if (!requireNamespace("network", quietly = TRUE))
        stop("The 'network' package must be installed.")

    M <- x$Matrix
    d <- nrow(M)

    ## sanity checks
    if (!inherits(x, "RVineMatrix"))
        stop("'x' has to be an RVineMatrix object.")
    if (tree != "ALL" && tree > d - 1)
        stop("Selected tree does not exist.")
    if (any(tree == "ALL") )
        tree <- 1:(d - 1)
    if (!all(type %in% c(0, 1, 2)))
        stop("type not implemented")
    stopifnot(is.logical(interactive))

    ## set names if empty
    if (is.null(x$names))
        x$names <- paste("V", 1:d, sep = "")

    #### set up plotting options ----------------------------
    # reduce default margins of plot range
    usr <- par(mar = c(1.1,0.1,3.1,0.1))
    on.exit(par(usr))

    # set plot.network options
    TUMblue <- rgb(0, 103/255, 198/255)
    TUMlightblue <- tint(TUMblue, 0.5)
    dflt <- list(interactive = interactive,
                 displaylabels = TRUE,
                 pad = 1.5e-1,
                 edge.lwd = 0.25,
                 edge.col = gray(0.3),
                 boxed.labels = TRUE,
                 label.pad = 1.5,
                 label.bg = TUMlightblue,
                 label.pos = 7,
                 label.col = "gray97",
                 label.cex = 1.3,
                 vertex.cex = 0,
                 object.scale = 0.05)
    # Same color for edges, edge labels and label borders
    dflt <- append(dflt, list(label.border   = dflt$edge.col,
                              edge.label.col = dflt$edge.col,
                              edge.label.cex = dflt$label.cex - 0.2))

    ## overwrite defaults with ... argument
    temp.args <- modifyList(dflt, list(...))

    #### loop through the trees -----------------------------
    for (i in tree) {

        main <- list(main = paste("Tree ", i, sep = ""),
                     col.main = ifelse("col.main" %in% names(temp.args),
                                       temp.args$col.main,
                                       temp.args$edge.col))
        final.args <- append(temp.args, main)

        ## create network object
        g <- makeNetwork(x, i, !(type %in% c(0, 2)))
        final.args$x = g$nw

        ## set edge labels
        if (!is.null(edge.labels))
            final.args$edge.label <- set_edge_labels(tree = i,
                                                     RVM = x,
                                                     edge.labels = edge.labels,
                                                     type = type)

        do.call(plot, final.args)

        ## add legend
        if (type == 2) {
            legend(legend.pos,
                   legend = paste(1:d, x$name, sep = " \U002194 "),
                   bty = "n",
                   xjust = 0,
                   text.col = final.args$edge.col,
                   cex = final.args$label.cex)
        }

        ## wait for key stroke
        if (i != max(tree)) {
            par(ask = TRUE)
        } else {
            par(ask = FALSE)
        }
    }
}


## creates a network object for a tree in a given RVineMatrix ------------------
makeNetwork <- function(RVM, tree, use.names = FALSE) {
    M <- RVM$Matrix
    d <- ncol(M)

    I <- matrix(0, d - tree + 1, d - tree + 1)

    ## extract node and edge labels as numbers
    if (tree > 1) {
        node.lab <- sapply(1:(d - tree + 1),
                           get_num,
                           tree = tree - 1,
                           RVM = RVM)
    } else {
        node.lab <- paste(diag(M))
    }
    edge.lab <- sapply(seq.int(d - tree),
                       get_num,
                       tree = tree,
                       RVM = RVM)

    ## convert to numeric matrices V and E
    V <- t(sapply(strsplit(node.lab,  " *[,;] *"), as.numeric))
    V <- matrix(V, ncol = tree)
    E <- t(sapply(strsplit(edge.lab,  " *[,;] *"), as.numeric))

    ## build incident matrix by matching V and E
    for (i in 1:nrow(E)) {
        ind.i <- which(apply(V, 1, function(x) all(x %in% E[i, ])))
        I[ind.i[1], ind.i[2]] <- I[ind.i[1], ind.i[2]] <- 1
    }

    ## convert to variable names (if asked for)
    if (use.names) {
        if (tree > 1) {
            node.lab <- sapply(1:(d - tree + 1),
                               get_name,
                               tree = tree - 1,
                               RVM = RVM)
        } else {
            node.lab <- RVM$names[diag(M)]
        }
    }

    ## create network
    colnames(I) <- rownames(I) <- node.lab
    nw <- network::network(I, directed = FALSE)

    ## return network and labels
    list(nw = nw, vlabs = node.lab)
}


## finds appropriate edge labels for the plot ----------------------------------
set_edge_labels <- function(tree, RVM, edge.labels, type) {
    d <- nrow(RVM$Matrix)
    if (edge.labels[1] == "family") {
        elabel <- sapply(1:(d - tree + 1),
                         get_family,
                         tree = tree,
                         RVM = RVM)
        elabel <- BiCopName(as.numeric(elabel))
    } else if (edge.labels[1] == "par") {
        elabel <- sapply(1:(d - tree + 1),
                          get_par,
                          tree = tree,
                          RVM = RVM)
    } else if (edge.labels[1] == "tau") {
        elabel <- sapply(1:(d - tree + 1),
                         get_tau,
                         tree = tree,
                         RVM = RVM)
    } else if (edge.labels[1] == "family-par") {
        elabel1 <- sapply(1:(d - tree + 1),
                          get_family,
                          tree = tree,
                          RVM = RVM)
        elabel1 <- BiCopName(as.numeric(elabel1))
        elabel2 <- sapply(1:(d - tree + 1),
                          get_par,
                          tree = tree,
                          RVM = RVM)
        elabel <- paste0(elabel1, "(", elabel2, ")")
        elabel <- sapply(elabel,
                         function(x){
                             tmp <- gsub("((", "(", x, fixed = TRUE)
                             gsub("))", ")", tmp, fixed = TRUE)
                         })
    } else if (edge.labels[1] == "family-tau") {
        elabel1 <- sapply(1:(d - tree + 1),
                          get_family,
                          tree = tree,
                          RVM = RVM)
        elabel1 <- BiCopName(as.numeric(elabel1))
        elabel2 <- sapply(1:(d - tree + 1),
                          get_tau,
                          tree = tree,
                          RVM = RVM)
        elabel <- paste0(elabel1, "(", elabel2, ")")
    } else if (length(edge.labels) > 1) {
        # user may provide own labels
        if (length(edge.labels) == d - tree) {
            elabel <- as.character(edge.labels)
        } else {
            stop("length of edge.labels does not equal the number of edges in the tree")
        }
    } else if (edge.labels[1] == "pair"){
        if (type %in% c(0, 2)) {
            elabel <- sapply(1:(d - tree + 1),
                             get_num,
                             tree = tree,
                             RVM = RVM)
        } else {
            elabel <- sapply(1:(d - tree + 1),
                             get_name,
                             tree = tree,
                             RVM = RVM)
        }
    } else {
        stop("edge.labels not implemented")
    }

    elabel
}


## get info for a pair-copula from RVineMatrix object --------------------------
get_num <-  function(j, tree, RVM) {
    M <- RVM$Matrix
    d <- nrow(M)
    # get numbers from structure matrix
    nums <- as.character(M[c(j, (d - tree + 1):d), j])
    # conditioned set
    bef <- paste(nums[2],
                 nums[1],
                 sep = ",",
                 collapse = "")
    # conditioning set
    aft <- if (length(nums) > 2) {
        gsub(" ",
             ",",
             do.call(paste, as.list(as.character(nums[3:length(nums)]))))
    }  else ""
    # paste together
    sep <- if (length(nums) > 2) " ; " else ""
    paste(bef, aft, sep = sep, collapse = "")
}

get_name <-  function(j, tree, RVM) {
    M <- RVM$Matrix
    d <- nrow(M)
    # variable names
    nams <- RVM$names[M[c(j, (d - tree + 1):d), j]]
    # conditioned set
    bef <- paste(nams[2],
                 nams[1],
                 sep = ",",
                 collapse = "")
    # conditioning set
    aft <- if (length(nams) > 2) {
        gsub(" ",  ",", do.call(paste, as.list(nams[3:length(nams)])))
    }  else ""
    # paste together
    sep <- if (length(nams) > 2) " ; " else ""
    paste(bef, aft, sep = sep, collapse = "")
}

get_family <- function(j, tree, RVM) {
    d <- nrow(RVM$family)
    paste(RVM$family[d - tree + 1, j])
}

get_par <- function(j, tree, RVM) {
    d <- nrow(RVM$family)
    # get parameters
    par  <- round(RVM$par[d - tree + 1, j], digits = 2)
    par2 <- round(RVM$par2[d - tree + 1, j], digits = 2)
    # add brackets if par2 != 0
    apply(cbind(par, par2), 1, join_par)
}

join_par <- function(x) {
    if (x[2] != 0)
        return(paste0("(", x[1], ",", x[2], ")"))
    x[1]
}

get_tau <- function(j, tree, RVM) {
    round(RVM$tau[nrow(RVM$tau) - tree + 1, j], digits = 2)
}
