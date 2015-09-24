#' Visualisation of R-Vine Tree Structure
#' 
#' Function is deprecated since \code{VineCopula 2.0}. Use
#' \code{\link[VineCopula:plot.RVineMatrix]{plot.RVineMatrix}} instead.
#' 
#' 
#' @param x \code{RVineMatrix} object.
#' @param tree \code{"ALL"} or integer vector; specifies which trees are
#' plotted.
#' @param type integer; specifies how to make use of variable names: \cr
#' \code{0} = variable names are ignored, \cr \code{1} = variable names are
#' used to annotate vertices, \cr \code{2} = uses numbers in plot and adds a
#' legend for variable names.
#' @param edge.labels character; either a vector of edge labels or one of the
#' following: \cr \code{"family"} = pair-copula family abbreviation (see
#' \code{\link[VineCopula:BiCopName]{BiCopName}}), \cr \code{"par"} =
#' pair-copula parameters, \cr \code{"tau"} = pair-copula Kendall's tau (by
#' conversion of parameters) \cr \code{"family-par"} = pair-copula family and
#' parameters \cr \code{"family-tau"} = pair-copula family and Kendall's tau.
#' @param legend.pos the \code{x} argument for
#' \code{\link[graphics:legend]{legend}}.
#' @param interactive logical; if TRUE, the user is asked to adjust the
#' positioning of vertices with his mouse.
#' @param \dots Arguments passed to
#' \code{\link[network:plot.network]{plot.network}}.
#' @author Thomas Nagler
#' @seealso \code{\link[VineCopula:plot.RVineMatrix]{plot.RVineMatrix}}
#' @export RVineTreePlot
RVineTreePlot <- function(x, tree = "ALL", type = 0, edge.labels = NULL, legend.pos = "bottomleft", interactive = FALSE, ...) {
    if (!inherits(x, "RVineMatrix")) {
        stop("'x' has to be an RVineMatrix object.")
    }
    warning("RVineTreePlot is deprecated and behaves differently compared to versions < 2.0. Use plot.RVineMatrix instead.")
    plot(x, tree, type, edge.labels, legend.pos, interactive, ...)
}


# RVineTreePlot <- function(data = NULL, RVM, method = "mle", max.df = 30,
#                           max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)), 
#                           tree = "ALL", edge.labels = c("family"), P = NULL, legend = FALSE) {
#     
#     if (is(RVM)[1] != "RVineMatrix") 
#         stop("'RVM' has to be an RVineMatrix object.")
#     
#     if (edge.labels[1] != FALSE & !all(edge.labels %in% c("family", "par", "par2", "theotau", "emptau", "pair"))) 
#         stop("Edge label not implemented.")
#     if (is.null(data) & any(edge.labels == "emptau")) 
#         stop("Empirical Kendall's tau values cannot be obtained if no data is provided.")
#     
#     if (is.null(data) == FALSE && (any(data > 1) || any(data < 0))) 
#         stop("Data has be in the interval [0,1].")
#     d <- dim(RVM)
#     
#     if (is.null(RVM$names)) 
#         RVM$names <- paste("V", 1:d, sep = "")
#     
#     empTauMat <- matrix(NA, d, d)
#     
#     if (!is.null(data)) {
#         
#         seqpar <- RVineSeqEstTau(data,
#                                  RVM, 
#                                  method = method,
#                                  se = FALSE,
#                                  max.df = max.df,
#                                  max.BB = max.BB, 
#                                  progress = FALSE)
#         
#         RVM$par <- seqpar$par
#         RVM$par2 <- seqpar$par2
#         
#         if (any(edge.labels == "emptau")) 
#             empTauMat <- round(seqpar$tau, 2)
#         
#     }
#     
#     theoTauMat <- round(RVinePar2Tau(RVM), 2)
#     if (any(edge.labels == "par")) 
#         parMat <- round(RVM$par, 2)
#     if (any(edge.labels == "par2")) 
#         parMat2 <- round(RVM$par2, 2)
#     
#     if (is.null(P)) {
#         P <- list()
#         for (i in 1:(d - 1)) P[[i]] <- 0
#     }
#     
#     if (tree != "ALL" && tree > d - 1) 
#         stop("Selected tree does not exist.")
#     
#     if (tree == "ALL") 
#         tree <- 1:(d - 1)
#     
#     M <- RVM$Matrix
#     
#     # cp. Alg. 3.1 in Dissmann
#     
#     edges <- list()
#     for (j in 1:(d - 1)) edges[[j]] <- array(NA, dim = c(d - j, 2, j))
#     
#     weight <- list()
#     for (j in 1:(d - 1)) weight[[j]] <- rep(NA, d - j)
#     
# #     # label the nodes
# #     for (j in 1:(d - 1))
# #         for (k in 1:d) edges[[j]][edges[[j]] == k] <- RVM$names[k]
#     
#     
#     if (edge.labels[1] != FALSE) {
#         numlabels <- length(edge.labels)
#         elabels <- list()
#         for (j in 1:(d - 1)) elabels[[j]] <- matrix(NA, d - j, numlabels)
#     }
#     
#     # initial edge
#     edges[[1]][1, , ] <- sort(c(M[d - 1, d - 1], M[d, d - 1]))
#     weight[[1]][1] <- ifelse(is.null(data), theoTauMat[d, d - 1], empTauMat[d, d - 1])
#     if (edge.labels[1] != FALSE) {
#         for (jj in 1:numlabels) {
#             if (edge.labels[jj] == "family") 
#                 elabels[[1]][1, jj] <- BiCopName(RVM$family[d, d - 1],
#                                                  short = TRUE)
#             if (edge.labels[jj] == "par") 
#                 elabels[[1]][1, jj] <- parMat[d, d - 1]
#             if (edge.labels[jj] == "par2") 
#                 elabels[[1]][1, jj] <- parMat2[d, d - 1]
#             if (edge.labels[jj] == "theotau") 
#                 elabels[[1]][1, jj] <- theoTauMat[d, d - 1]
#             if (edge.labels[jj] == "emptau") 
#                 elabels[[1]][1, jj] <- empTauMat[d, d - 1]
#             if (edge.labels[jj] == "pair") 
#               if (legend == TRUE) {
#                 elabels[[1]][1, jj] <- paste(RVM$Matrix[d - 1, d - 1],
#                                              RVM$Matrix[d, d - 1],
#                                              sep = ",")
#               } else {
#                 elabels[[1]][1, jj] <- paste(RVM$names[RVM$Matrix[d - 1, d - 1]],
#                                              RVM$names[RVM$Matrix[d, d - 1]],
#                                              sep = ",")
#               }
#         }
#     }
#     
#     for (i in (d - 2):1) {
#         
#         # new edge in first tree
#         ee <- sort(c(M[i, i], M[d, i]))
#         edges[[1]][d - i, , ] <- ee
#         weight[[1]][d - i] <- ifelse(is.null(data), theoTauMat[d, i], empTauMat[d, i])
#         if (edge.labels[1] != FALSE) {
#             for (jj in 1:numlabels) {
#                 if (edge.labels[jj] == "family") 
#                     elabels[[1]][d - i, jj] <- BiCopName(RVM$family[d, i], 
#                                                          short = TRUE)
#                 if (edge.labels[jj] == "par") 
#                     elabels[[1]][d - i, jj] <- parMat[d, i]
#                 if (edge.labels[jj] == "par2") 
#                     elabels[[1]][d - i, jj] <- parMat2[d, i]
#                 if (edge.labels[jj] == "theotau") 
#                     elabels[[1]][d - i, jj] <- theoTauMat[d, i]
#                 if (edge.labels[jj] == "emptau") 
#                     elabels[[1]][d - i, jj] <- empTauMat[d, i]
#                 if (edge.labels[jj] == "pair") 
#                   if (legend == TRUE) {
#                     elabels[[1]][d - i, jj] <- paste(RVM$Matrix[i, i],
#                                                      RVM$Matrix[d, i],
#                                                      sep = ",")
#                   } else {
#                     elabels[[1]][d - i, jj] <- paste(RVM$names[RVM$Matrix[i, i]],
#                                                      RVM$names[RVM$Matrix[d, i]],
#                                                      sep = ",")
#                   }
#             }
#         }
#         # edges in further trees
#         for (k in 1:(d - i - 1)) {
#             edges[[k + 1]][d - i - k, 1, ] <- ee
#             
#             # identify conditioned and conditioning sets
#             if (length(M[(d - k):d, i]) >= 3) {
#                 if (setequal(M[(d - k):d, i], ee_old)) {
#                     edges[[k + 1]][d - i - k, 2, ] <- ee_old
#                 } else {
#                     for (j in 1:(d - i - k)) {
#                         if (setequal(M[(d - k):d, i], edges[[k + 1]][j, 1, ])) 
#                             edges[[k + 1]][d - i - k, 2, ] <- edges[[k + 1]][j, 1, ]
#                         if (setequal(M[(d - k):d, i], edges[[k + 1]][j, 2, ])) 
#                             edges[[k + 1]][d - i - k, 2, ] <- edges[[k + 1]][j, 2, ]
#                     }
#                 }
#             } else {
#                 edges[[k + 1]][d - i - k, 2, ] <- sort(M[(d - k):d, i])
#             }
#             
#             # create edge lables
#             weight[[k + 1]][d - i - k] <- ifelse(is.null(data), theoTauMat[d - k, i], empTauMat[d - k, i])
#             if (edge.labels[1] != FALSE) {
#                 for (jj in 1:numlabels) {
#                     if (edge.labels[jj] == "family") 
#                         elabels[[k + 1]][d - i - k, jj] <- BiCopName(RVM$family[d - k, i], short = TRUE)
#                     if (edge.labels[jj] == "par") 
#                         elabels[[k + 1]][d - i - k, jj] <- parMat[d - k, i]
#                     if (edge.labels[jj] == "par2") 
#                         elabels[[k + 1]][d - i - k, jj] <- parMat2[d - k, i]
#                     if (edge.labels[jj] == "theotau") 
#                         elabels[[k + 1]][d - i - k, jj] <- theoTauMat[d - k, i]
#                     if (edge.labels[jj] == "emptau") 
#                         elabels[[k + 1]][d - i - k, jj] <- empTauMat[d - k, i]
#                     if (edge.labels[jj] == "pair") {
#                       if (legend == TRUE) {
#                         handle1 <- paste(RVM$Matrix[i, i], 
#                                          RVM$Matrix[d - k, i],
#                                          sep = ",")
#                         handle2 <- paste(RVM$Matrix[(d - k + 1):d, i],
#                                          collapse = ",")
#                         handle3 <- paste(handle1, 
#                                          handle2, 
#                                          sep = ";")
#                       } else {
#                         handle1 <- paste(RVM$names[RVM$Matrix[i, i]], 
#                                          RVM$names[RVM$Matrix[d - k, i]],
#                                          sep = ",")
#                         handle2 <- paste(RVM$names[RVM$Matrix[(d - k + 1):d, i]],
#                                          collapse = ",")
#                         handle3 <- paste(handle1, 
#                                          handle2, 
#                                          sep = ";")
#                       }
#                       elabels[[k + 1]][d - i - k, jj] <- handle3  #paste(handle1,handle2,sep=';')
#                     }
#                 }
#             }
#             
#             # identify conditioned and conditioning sets
#             ee <- c(sort(c(setdiff(ee, M[(d - k):d, i]), 
#                            setdiff(M[(d - k):d, i], ee))),
#                     sort(intersect(ee, M[(d - k):d, i])))
#         }
#         
#         ee_old <- ee
#         
#     }
#     
#     # label the nodes
#     if (legend == FALSE) {
#       for (j in 1:(d - 1)) for (k in 1:d) edges[[j]][edges[[j]] == k] <- RVM$names[k]
#     }
#     
#     # convert to edge lists
#     edgelist <- list()
#     for (j in 1:(d - 1)) edgelist[[j]] <- matrix(NA, d - j, 2)
#     
#     edgelist[[1]] <- matrix(as.character(edges[[1]][, , 1]), d - 1, 2)
#     
#     for (j in 1:(d - 2)) edgelist[[2]][j, ] <- c(paste(edges[[2]][j, 1, ], collapse = ","), 
#                                                  paste(edges[[2]][j, 2, ], collapse = ","))
#     
#     # separate conditioned and conditioning sets
#     if (d > 3) {
#         for (i in 3:(d - 1)) {
#             for (j in 1:(d - i)) {
#                 edgelist[[i]][j, 1] <- paste(paste(edges[[i]][j, 1, 1:2], collapse = ","),
#                                              paste(edges[[i]][j, 1, 3:i], collapse = ","), sep = ";")
#                 edgelist[[i]][j, 2] <- paste(paste(edges[[i]][j, 2, 1:2], collapse = ","),
#                                              paste(edges[[i]][j, 2, 3:i], collapse = ","), sep = ";")
#             }
#         }
#     }
#     
#     # combine edge lables
#     if (edge.labels[1] != FALSE) {
#         elabels2 <- list()
#         for (j in 1:(d - 1)) {
#             elabels2[[j]] <- rep(NA, d - j)
#             for (i in 1:(d - j)) elabels2[[j]][i] <- paste(elabels[[j]][i, ], collapse = ",")
#         }
#     }
#     
#     # create graphs
#     gg <- list()
#     for (i in 1:(d - 1)) {
#         gg[[i]] <- graph_from_edgelist(edgelist[[i]], directed = FALSE)
#         E(gg[[i]])$weight <- weight[[i]]
#         if (edge.labels[1] != FALSE) 
#             E(gg[[i]])$name <- elabels2[[i]]
#     }
#     
#     # loop through the trees
#     for (i in tree) {
#         
#         g <- gg[[i]]
#         
#         if (edge.labels[1] != FALSE) {
#             elabel <- E(g)$name
#         } else {
#             elabel <- NULL
#         }
#         
#         ## specify layout for plotting
#         if (all(P[[i]] == 0)) {
#             P[[i]] <- layout_in_circle(g)
#             P[[i]] <- layout_with_graphopt(g, start = P[[i]], niter = 50, spring.length = 1)
#         }
#         
#         ## initialize plotting
#         main <- paste("Tree ", i, sep = "")
#         if (legend == TRUE) {
#           vwidth <- max(strwidth(V(g)$name, units = "figure")) * 800
#           vheight <- 20
#         } else {
#           vwidth <- max(strwidth(V(g)$name, units = "figure")) * 1000
#           vheight <- 20
#         }
#         
#         ## plot tree
#         plot(g, layout = P[[i]],
#              vertex.label = V(g)$name,
#              vertex.shape = "rectangle",
#              vertex.size = vwidth,
#              vertex.size2 = vheight,
#              edge.label.family = "sans",
#              edge.label = elabel,
#              edge.width = (10*abs(E(g)$weight) + 0.5),
#              edge.arrow.size = 0,
#              main = main)
#         if (legend == TRUE) {
#           legend("bottomleft", legend = paste(1:d, RVM$name, sep = " = "),
#                  bty = "n", xjust = 0)
#         }
#         
#         if (i != max(tree)) {
#             par(ask = TRUE)
#         } else {
#             par(ask = FALSE)
#         }
#         
#     }
#     
#     return(P)
#     
# }
# 
# RVineSeqEstTau <- function(data, RVM, method = "mle", se = FALSE, max.df = 30, max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)), 
#                            progress = FALSE) {
#     data <- as.matrix(data)
#     n <- dim(RVM)
#     N <- nrow(data)
#     if (dim(data)[2] != dim(RVM)) 
#         stop("Dimensions of 'data' and 'RVM' do not match.")
#     if (N < 2) 
#         stop("Number of observations has to be at least 2.")
#     if (!("RVineMatrix" %in% is(RVM))) 
#         stop("'RVM' has to be an RVineMatrix object.")
#     
#     if (method != "mle" && method != "itau") 
#         stop("Estimation method has to be either 'mle' or 'itau'.")
#     
#     if (max.df <= 2) 
#         stop("The upper bound for the degrees of freedom parameter has to be larger than 2.")
#     if (!is.list(max.BB)) 
#         stop("'max.BB' has to be a list.")
#     if (max.BB$BB1[1] < 0.001) 
#         stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
#     if (max.BB$BB1[2] < 1.001) 
#         stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
#     if (max.BB$BB6[1] < 1.001) 
#         stop("The upper bound for the first parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
#     if (max.BB$BB6[2] < 1.001) 
#         stop("The upper bound for the second parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
#     if (max.BB$BB7[1] < 1.001) 
#         stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
#     if (max.BB$BB7[2] < 0.001) 
#         stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")
#     if (max.BB$BB8[1] < 1.001) 
#         stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
#     if (max.BB$BB8[2] < 0.001 || max.BB$BB8[2] > 1) 
#         stop("The upper bound for the second parameter of the BB1 copula should be in the interval [0,1].")
#     
#     o <- diag(RVM$Matrix)
#     
#     if (any(o != length(o):1)) {
#         oldRVM <- RVM
#         RVM <- normalizeRVineMatrix(RVM)
#         data <- data[, o[length(o):1]]
#     }
#     
#     Params <- RVM$par
#     Params2 <- RVM$par2
#     
#     if (se == TRUE) {
#         seMat1 <- matrix(0, nrow = n, ncol = n)
#         seMat2 <- matrix(0, nrow = n, ncol = n)
#     }
#     
#     empTauMat <- matrix(0, nrow = n, ncol = n)
#     
#     V <- list()
#     V$direct <- array(NA, dim = c(n, n, N))
#     V$indirect <- array(NA, dim = c(n, n, N))
#     
#     V$direct[n, , ] <- t(data[, n:1])
#     
#     for (i in (n - 1):1) {
#         
#         for (k in n:(i + 1)) {
#             
#             m <- RVM$MaxMat[k, i]
#             zr1 <- V$direct[k, i, ]
#             
#             if (m == RVM$Matrix[k, i]) {
#                 zr2 <- V$direct[k, (n - m + 1), ]
#             } else {
#                 zr2 <- V$indirect[k, (n - m + 1), ]
#             }
#             
#             
#             if (RVM$family[k, i] == 2 | RVM$family[k, i] == 7 | RVM$family[k, i] == 8 | RVM$family[k, i] == 9) {
#                 if (progress == TRUE) {
#                     if (k == n) {
#                         message(oldRVM$Matrix[i, i], ",", oldRVM$Matrix[k, i]) 
#                     } else { 
#                         message(oldRVM$Matrix[i, i], ",", oldRVM$Matrix[k, i],  ";",
#                                 paste(oldRVM$Matrix[(k + 1):n, i], collapse = ","))
#                     }
#                 }
#                 par.out <- BiCopEst(zr2, 
#                                     zr1,
#                                     RVM$family[k, i], 
#                                     method,
#                                     se,
#                                     max.df, 
#                                     max.BB)
#                 # par1 <- out.par$par
#                 Params[k, i] <- par.out$par
#                 Params2[k, i] <- par.out$par2
#                 # empTauMat[k,i] = cor(zr2,zr1,method='kendall')
#                 empTauMat[k, i] <- fasttau(zr2, zr1)
#                 if (se == TRUE) {
#                     # se1 <- par.out$se
#                     seMat1[k, i] <- par.out$se
#                     seMat2[k, i] <- par.out$se2
#                 }
#             } else {
#                 if (progress == TRUE) {
#                     if (k == n) {
#                         message(oldRVM$Matrix[i, i], ",", oldRVM$Matrix[k, i]) 
#                     } else { 
#                         message(oldRVM$Matrix[i, i], ",", oldRVM$Matrix[k, i], ";",
#                                 paste(oldRVM$Matrix[(k + 1):n, i], collapse = ","))
#                     }
#                 }
#                 par.out <- BiCopEst(zr2, zr1, RVM$family[k, i], method, se, max.df, max.BB)
#                 Params[k, i] <- par.out$par
#                 empTauMat[k, i] <- cor(zr2, zr1, method = "kendall")
#                 empTauMat[k, i] <- fasttau(zr2, zr1)
#                 if (se == TRUE) {
#                     seMat1[k, i] <- par.out$se
#                 }
#             }
#             
#             
#             if (RVM$CondDistr$direct[k - 1, i]) {
#                 V$direct[k - 1, i, ] <- .C("Hfunc1",
#                                            as.integer(RVM$family[k, i]), 
#                                            as.integer(length(zr1)),
#                                            as.double(zr1),
#                                            as.double(zr2), 
#                                            as.double(Params[k, i]),
#                                            as.double(Params2[k, i]),
#                                            as.double(rep(0, length(zr1))),
#                                            PACKAGE = "VineCopula")[[7]]
#             }
#             if (RVM$CondDistr$indirect[k - 1, i]) {
#                 V$indirect[k - 1, i, ] <- .C("Hfunc2",
#                                              as.integer(RVM$family[k, i]),
#                                              as.integer(length(zr2)),
#                                              as.double(zr2), 
#                                              as.double(zr1),
#                                              as.double(Params[k, i]), 
#                                              as.double(Params2[k, i]),
#                                              as.double(rep(0, length(zr1))),
#                                              PACKAGE = "VineCopula")[[7]]
#             }
#             
#         }
#     }
#     
#     if (se == FALSE) {
#         return(list(par = Params, 
#                     par2 = Params2, 
#                     tau = empTauMat)) 
#     } else {
#         return(list(par = Params, 
#                     par2 = Params2, 
#                     tau = empTauMat, 
#                     se = seMat1, 
#                     se2 = seMat2))
#     }
# }
