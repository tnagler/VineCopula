#' Pairs Plot of Copula Data
#'
#' This function provides pair plots for copula data. Using default setting it
#' plots bivariate contour plots on the lower panel, scatter plots and
#' correlations on the upper panel and histograms on the diagonal panel.
#'
#'
#' @param x \code{copuladata} object.
#' @param labels variable names/labels.
#' @param \dots other graphical parameters (see \code{\link[graphics]{par}}) or
#' options passed to \code{\link{BiCopKDE}}.
#' @param lower.panel panel function to be used on the lower diagonal panels
#' (if not supplied, a default function is used)
#' @param upper.panel panel function to be used on the upper diagonal panels
#' (if not supplied, a default function is used)
#' @param diag.panel panel function to be used on the diagonal panels (if not
#' supplied, a default function is used)
#' @param label.pos y position of labels in the diagonal panel; default:
#' \code{label.pos = 0.85}.
#' @param cex.labels magnification to be used for the labels of the diagonal
#' panel; default: \code{cex.labels = 1}.
#' @param gap distance between subplots, in margin lines; default: \code{gap =
#' 0}.
#' @param method a character string indicating which correlation coefficients
#' are computed. One of \code{"pearson"}, \code{"kendall"} (default), or
#' \code{"spearman"}
#' @param ccols colour to be used for the contour plots; default: \code{ccols =
#' terrain.colors(30)}.
#' @param margins character; margins for the contour plots. Options are:\cr
#' \code{"unif"} for the original copula density,
#' \code{"norm"} for the transformed density with standard normal margins,
#' \code{"exp"} with standard exponential margins, and  \code{"flexp"} with
#' flipped exponential margins.
#' @note If the default panel functions are used \cr \itemize{ \item \code{col}
#' changes only the colour of the points in the scatter plot
#' (\code{upper.panel}) \cr \item \code{cex} changes only the magnification of
#' the points in the scatter plot (\code{upper.panel}) }
#' @author Tobias Erhardt
#' @seealso \code{\link[graphics]{pairs}}, \code{\link{as.copuladata}},
#' \code{\link{BiCopKDE}}
#' @examples
#'
#' data(daxreturns)
#' \dontshow{daxreturns <- daxreturns[1:50, ]}
#' data <- as.copuladata(daxreturns)
#' sel <- c(4,5,14,15)
#'
#' ## pairs plot with default settings
#' pairs(data[sel])
#'
#' ## pairs plot with custom settings
#' nlevels <- 20
#' pairs(data[sel], cex = 2, pch = 1, col = "black",
#'       diag.panel = NULL, label.pos = 0.5,
#'       cex.labels = 2.5, gap = 1,
#'       method = "pearson", ccols = heat.colors(nlevels),
#'       margins = "flexp")
#'
#' ## pairs plot with own panel functions
#' up <- function(x, y) {
#'   # upper panel: empirical contour plot
#'   op <- par(usr = c(-3, 3, -3, 3), new = TRUE)
#'   BiCopKDE(x, y,
#'            levels = c(0.01, 0.05, 0.1, 0.15, 0.2),
#'            margins = "exp",
#'            axes = FALSE)
#'   on.exit(par(op))
#' }
#'
#' lp <- function(x, y) {
#'   # lower panel: scatter plot (copula data) and correlation
#'   op <- par(usr = c(0, 1, 0, 1), new = TRUE)
#'   points(x, y, pch = 1, col = "black")
#'   r <- cor(x, y, method = "spearman") # Spearman's rho
#'   txt <- format(x = r, digits = 3, nsmall = 3)[1]
#'   text(x = 0.5, y = 0.5, labels = txt, cex = 1 + abs(r) * 2, col = "blue")
#'   on.exit(par(op))
#' }
#'
#' dp <- function(x) {
#'   # diagonal panel: histograms (copula data)
#'   op <- par(usr = c(0, 1, 0, 1.5), new = TRUE)
#'   hist(x, freq = FALSE, add = TRUE, col = "brown", border = "black", main = "")
#'   abline(h = 1, col = "black", lty = 2)
#'   on.exit(par(op))
#' }
#'
#' nlevels <- 20
#' pairs(data[sel],
#'       lower.panel = lp, upper.panel = up, diag.panel = dp, gap = 0.5)
#'
#' @export pairs.copuladata
pairs.copuladata <- function(x,
                             labels = names(x),
                             ...,
                             lower.panel = lp.copuladata,
                             upper.panel = up.copuladata,
                             diag.panel = dp.copuladata,
                             label.pos = 0.85,
                             cex.labels = 1,
                             gap = 0,
                             method = "kendall",
                             ccols = terrain.colors(11),
                             margins = "norm") {
  ## pairs plot for 'copuladata'

  # panel functions
  if (!is.null(lower.panel) && missing(lower.panel)) {
    lp <- function(x, y, cc = ccols, mar = margins, ...) {
      lower.panel(x, y, cc = cc, mar = mar, ...)
    }
  } else {
    lp <- lower.panel
  }
  if (!is.null(upper.panel) && missing(upper.panel)) {
    up <- function(x, y, mthd = method, ...) {
      upper.panel(x, y, mthd = mthd, ...)
    }
  } else {
    up <- upper.panel
  }
  if (!is.null(diag.panel) && missing(diag.panel)) {
    dp <- function(x, ...) {
      diag.panel(x, ...)
    }
  } else {
    dp <- diag.panel
  }

  # pairs plot (with panel functions as defined below or as provided by user)
  op <- par(xaxt = "n", yaxt = "n")
  #do.call(pairs, modifyList(default, list(...)))
  pairs(x = as.matrix(x),
        labels = labels,
        ...,
        lower.panel = lp,
        upper.panel = up,
        diag.panel = dp,
        label.pos = label.pos,
        cex.labels = cex.labels,
        gap = gap)
  on.exit(par(op))
}


## lower panel: empirical contour plot
lp.copuladata <- function(x, y, cc, mar, marpar, ...) {
  # set number of levels and if contour labels are drawn
  if (length(cc) == 1) {
    nlvls <- 11
    dl <- TRUE
  } else {
    nlvls <- length(cc)
    dl <- FALSE
  }
  # set levels according to margins
  if (mar %in% c("norm", "exp", "flexp")) {
    lvls <- seq(0.0, 0.25, length.out = nlvls)
  } else if (mar == "unif") {
    lvls <- seq(0.3, 4, length.out = nlvls)
  }
  # set default parameters
  pars <- list(u1 = x,
               u2 = y,
               size = 100,
               levels = lvls,
               margins = mar,
               axes = FALSE,
               drawlabels = dl)
  # get non-default parameters
  pars <- modifyList(pars, list(...))
  pars <- modifyList(modifyList(pars, list(col = NULL)), list(col = cc))
  op <- par(usr = c(-3, 3, -3, 3), new = TRUE)
  # call BiCopMetaContour
  do.call(BiCopKDE, pars)
  on.exit(par(op))
}


## upper panel: scatter plot (copula data) and correlation
up.copuladata <- function(x, y, mthd, ...) {
  # set default parameters
  pars <- list(x = x,
               y = y,
               pch = ".",
               cex = 1,
               col = "grey"
  )
  # get non-default parameters
  pars <- modifyList(pars, list(...))
  op <- par(usr = c(0, 1, 0, 1), new = TRUE)
  # call points (to produce scatter plot)
  do.call(points, pars)
  r <- cor(x = x, y = y, method = mthd)
  txt <- format(x = r, digits = 2, nsmall = 2)[1]
  # call text
  pars.txt <- list(x = 0.5,
                   y = 0.5,
                   labels = txt)
  pars.txt <- modifyList(pars.txt, list(...))
  do.call(text, modifyList(pars.txt, list(cex = 1 + abs(r) * 3, col = "red")))
  on.exit(par(op))
}


## diagonal panel: histograms (copula data)
dp.copuladata <- function(x, ...) {
  # set default parameters
  pars <- list(x = x,
               freq = FALSE,
               add = TRUE,
               border = "white",
               main = "")
  # get non-default parameters
  pars <- modifyList(modifyList(pars, list(...)), list(col = "grey"))
  op <- par(usr = c(0, 1, 0, 1.6), new = TRUE)
  # call hist
  do.call(hist, pars)
  box()
  if (pars$freq == FALSE)
    abline(h = 1, col = "black", lty = 3)
  on.exit(par(op))
}
