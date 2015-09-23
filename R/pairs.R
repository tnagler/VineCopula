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
                             ccols = terrain.colors(30),
                             margins = "norm",
                             margins.par = 0) {
  ## pairs plot for 'copuladata'
  
  # panel functions
  if (!is.null(lower.panel) && missing(lower.panel)) {
    lp <- function(x, y, cc = ccols, mar = margins, marpar = margins.par, ...) {
      lower.panel(x, y, cc = cc, mar = mar, marpar = marpar, ...)
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
  if (mar %in% c("norm", "t", "exp")) {
    lvls <- seq(0.0, 0.2, length.out = nlvls)
  } else if (mar == "gamma") {
    lvls <- seq(0.0, 0.1, length.out = nlvls)
  } else if (mar == "unif") {
    lvls <- seq(0.1, 1.5, length.out = nlvls)
  }
  # set default parameters
  pars <- list(u1 = x,
               u2 = y,
               bw = 2, 
               size = 100, 
               levels = lvls,
               margins = mar, 
               margins.par = marpar, 
               xylim = NA,
               axes = FALSE,
               drawlabels = dl)
  # get non-default parameters
  pars <- modifyList(pars, list(...))
  pars <- modifyList(modifyList(pars, list(col = NULL)), list(col = cc))
  op <- par(usr = c(-3, 3, -3, 3), new = TRUE)
  # call BiCopMetaContour
  do.call(BiCopMetaContour, pars)
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
