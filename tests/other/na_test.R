library("VineCopula")


## BiCopCDF -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopCDF(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopCDF(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopCDF(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopCDF(s[, 1], s[, 2], BiCop(5, 4))
BiCopCDF(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopCDF(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopCDF(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopCDF(s[, 1], s[, 2], 2, 1:10/11, 4))
suppressWarnings(BiCopCDF(0.5, 0.5, 3, 3))
suppressWarnings(BiCopCDF(NA, u2 = 0.5, family = 3, par = 3))

## BiCopChiPlot --------------------------------------
s <- BiCopSim(500, 3 ,3)
op <- par(mfrow = c(1, 3))
BiCopChiPlot(s[,1], s[,2], xlim = c(-1,1), ylim = c(-1,1),
             main="General chi-plot")
BiCopChiPlot(s[,1], s[,2], mode = "lower", xlim = c(-1,1),
             ylim = c(-1,1), main = "Lower chi-plot")
BiCopChiPlot(s[,1], s[,2], mode = "upper", xlim = c(-1,1),
             ylim = c(-1,1), main = "Upper chi-plot")
s[1, 2] <- NA
BiCopChiPlot(s[,1], s[,2], xlim = c(-1,1), ylim = c(-1,1),
             main="General chi-plot")
BiCopChiPlot(s[,1], s[,2], mode = "lower", xlim = c(-1,1),
             ylim = c(-1,1), main = "Lower chi-plot")
BiCopChiPlot(s[,1], s[,2], mode = "upper", xlim = c(-1,1),
             ylim = c(-1,1), main = "Upper chi-plot")
par(op)

## BiCopKPlot  --------------------------------------
s <- BiCopSim(500, 3 ,3)
BiCopKPlot(s[,1], s[,2])
s[1, 2] <- NA
BiCopKPlot(s[,1], s[,2])

## BiCopCompare ---------------------------------------
# s <- BiCopSim(100, 3, 3)
# fit <- BiCopCompare(s[, 1], s[, 2])
# s[2, 1] <- s[10, 2] <-- NA
# fit <- BiCopCompare(s[, 1], s[, 2])

## BiCopDeriv -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopDeriv(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopDeriv(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopDeriv(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopDeriv(u2 = 0.5, family = 3, par = 3), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopDeriv(s[, 1], s[, 2], BiCop(5, 4))
BiCopDeriv(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopDeriv(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- s[10, 2] <-- NA
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 3, 1:10/11))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 3, 1:10/11, deriv = "u1"))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 3, 1:10/11, deriv = "u2"))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 2, 1:10/11, 4))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par2"))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "u1"))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "u2"))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 3, 1:10/11, log = TRUE))
suppressWarnings(BiCopDeriv(s[, 1], s[, 2], 2, 1:10/11, 4, log = TRUE))

## BiCopDeriv2 -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopDeriv2(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopDeriv2(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopDeriv2(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopDeriv2(u2 = 0.5, family = 3, par = 3), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopDeriv2(s[, 1], s[, 2], BiCop(5, 4))
BiCopDeriv2(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopDeriv2(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- s[10, 2] <-- NA
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 3, 1:10/11))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 3, 1:10/11, deriv = "u1"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 3, 1:10/11, deriv = "u2"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 3, 1:10/11, deriv = "u1"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 3, 1:10/11, deriv = "par1u1"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 3, 1:10/11, deriv = "par1u2"))

suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par2"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "u1"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "u2"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par1u1"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par1u2"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par1par2"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par2u1"))
suppressWarnings(BiCopDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par2u2"))

## BiCopEst -------------------------------------------
s <- BiCopSim(500, 1, 0.7)
s[2, 1] <- s[10, 2] <-- NA

BiCopEst(s[, 1], s[, 2], family = 3, method = "itau")
BiCopEst(s[, 1], s[, 2], family = 2, method = "mle" ,
         max.BB = list(BB1 = c(1, 6), BB6 = c(6, 6), BB7 = c(5, 6), BB8 = c(6, 1)))

## BiCopEstList -------------------------------------------
s <- BiCopSim(500, 1, 0.7)
s[, 1] <- 1 - s[, 1]
e <- try(BiCopEstList(s[, 1], s[, 2], familyset = 3:4, method = "itau", rotations = FALSE))
stopifnot(inherits(e, "try-error"))
BiCopEstList(s[, 1], s[, 2], familyset = 3:4, method = "mle")
s[2, 1] <- s[10, 2] <-- NA
BiCopEstList(s[, 1], s[, 2], familyset = 1:4, method = "itau")
BiCopEstList(s[, 1], s[, 2], familyset = 1:4, method = "mle")

## BiCopGofTest --------------------------
set.seed(123)
simdata <- BiCopSim(300, 3, 2)
u1 <- simdata[,1]
u2 <- simdata[,2]
u1[1] <- u2[3] <- NA
BiCopGofTest(u1, u2, family = 3)
BiCopGofTest(u1, u2, family = 5)
gof <- BiCopGofTest(u1, u2, family = 3, method = "kendall", B=50)
gof <- BiCopGofTest(u1, u2, family = 5, method = "kendall", B=50)

## BiCopHfunc -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHfunc(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfunc(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfunc(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHfunc(s[, 1], s[, 2], BiCop(5, 4))
BiCopHfunc(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHfunc(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopHfunc(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHfunc(s[, 1], s[, 2], 2, 1:10/11, 4))

## BiCopHfunc1 -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHfunc1(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfunc1(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfunc1(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHfunc1(s[, 1], s[, 2], BiCop(5, 4))
BiCopHfunc1(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHfunc1(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopHfunc1(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHfunc1(s[, 1], s[, 2], 2, 1:10/11, 4))

## BiCopHfunc2 -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHfunc2(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfunc2(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfunc2(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHfunc2(s[, 1], s[, 2], BiCop(5, 4))
BiCopHfunc2(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHfunc2(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopHfunc2(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHfunc2(s[, 1], s[, 2], 2, 1:10/11, 4))

## BiCopHfuncDeriv -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHfuncDeriv(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfuncDeriv(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfuncDeriv(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfuncDeriv(u2 = 0.5, family = 3, par = 3), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHfuncDeriv(s[, 1], s[, 2], BiCop(5, 4))
BiCopHfuncDeriv(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHfuncDeriv(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- s[10, 2] <-- NA
suppressWarnings(BiCopHfuncDeriv(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHfuncDeriv(s[, 1], s[, 2], 3, 1:10/11))
suppressWarnings(BiCopHfuncDeriv(s[, 1], s[, 2], 3, 1:10/11, deriv = "u2"))
suppressWarnings(BiCopHfuncDeriv(s[, 1], s[, 2], 2, 1:10/11, 4))
suppressWarnings(BiCopHfuncDeriv(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par2"))
suppressWarnings(BiCopHfuncDeriv(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "u2"))

## BiCopHfuncDeriv2 -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHfuncDeriv2(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfuncDeriv2(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfuncDeriv2(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHfuncDeriv2(u2 = 0.5, family = 3, par = 3), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHfuncDeriv2(s[, 1], s[, 2], BiCop(5, 4))
BiCopHfuncDeriv2(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHfuncDeriv2(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- s[10, 2] <-- NA
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 3, 1:10/11))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 3, 1:10/11, deriv = "u2"))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 3, 1:10/11, deriv = "par1u2"))

suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par2"))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "u2"))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par1u2"))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par1par2"))
suppressWarnings(BiCopHfuncDeriv2(s[, 1], s[, 2], 2, 1:10/11, 4, deriv = "par2u2"))

## BiCopHinv -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHinv(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHinv(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHinv(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHinv(s[, 1], s[, 2], BiCop(5, 4))
BiCopHinv(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHinv(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopHinv(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHinv(s[, 1], s[, 2], 2, 1:10/11, 4))

## BiCopHinv1 -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHinv1(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHinv1(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHinv1(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHinv1(s[, 1], s[, 2], BiCop(5, 4))
BiCopHinv1(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHinv1(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopHinv1(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHinv1(s[, 1], s[, 2], 2, 1:10/11, 4))

## BiCopHinv2 -------------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopHinv2(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHinv2(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopHinv2(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopHinv2(s[, 1], s[, 2], BiCop(5, 4))
BiCopHinv2(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopHinv2(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopHinv2(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopHinv2(s[, 1], s[, 2], 2, 1:10/11, 4))

## BiCopIndTest -------------------------
dat <- BiCopSim(500, BiCop(3, 3))
dat[1, 1] <- dat[3, 2] <- NA
BiCopIndTest(dat[, 1], dat[, 2])

## BiCopKDE --------------------------------------
cop <- BiCop(3, tau = 0.3)
u <- BiCopSim(1000, cop)
u[1, 1] <- u[3, 2] <- NA
contour(cop)  # true contours
BiCopKDE(u[, 1], u[, 2])
BiCopKDE(u[, 1], u[, 2], kde.pars = list(mult = 0.5))  # undersmooth
BiCopKDE(u[, 1], u[, 2], kde.pars = list(mult = 2))  # oversmooth
BiCopKDE(u[, 1], u[, 2], type = "surface", zlim = c(0, 4))
plot(cop, zlim = c(0, 4))  # true density

## BiCopLambda ------------------------------------
cop <- BiCop(3, tau = 0.5)
dat <- BiCopSim(1000, cop)
dat[1, 2] <- dat[3, 2] <- NA
op <- par(mfrow = c(1, 3))
BiCopLambda(dat[, 1], dat[, 2])  # empirical lambda-function
BiCopLambda(cop)	# theoretical lambda-function
BiCopLambda(dat[, 1], dat[, 2], cop)	# both
par(op)

## BiCopMetaContour -------------------------
cop <- BiCop(family = 1, tau = 0.5)
BiCopMetaContour(obj = cop, main = "Clayton - normal margins")
dat <- BiCopSim(1000, cop)
dat[1, 1] <- dat[3, 2] <- NA
BiCopMetaContour(dat[, 1], dat[, 2], bw = 2, family = "emp",
                 main = "empirical - normal margins")

## BiCopPar2Beta ----------------------------
BiCopPar2Beta(family = 3, par = 2)
BiCopPar2Beta(BiCop(3, 3))
BiCopPar2Beta(family = c(3,4,6), par = 2:4)

## BiCopPar2TailDep ----------------------------
BiCopPar2TailDep(1, 0.7)
BiCopPar2TailDep(BiCop(3, 3))
BiCop(1, 0.7)$taildep  # alternative
BiCopPar2TailDep(2, c(0.6, 0.7, 0.8), 4)
BiCopPar2TailDep(c(3, 4, 6), 2)

## BiCopPar2Tau ----------------------------
BiCopPar2Tau(1, 0.7)
BiCopPar2Tau(BiCop(3, 3))
BiCop(1, 0.7)$tau  # alternative
BiCopPar2Tau(2, c(0.6, 0.7, 0.8), 4)
BiCopPar2Tau(c(3, 4, 6), 2)

## BiCopPDF -------------------------------
s <- BiCopSim(10, 3, 3)
e <- try(BiCopPDF(s[, 1], s[, 2], 1:5, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopPDF(s[, 1], s[, 2], 1:10, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
e <- try(BiCopPDF(s[, 1], s[, 2], 11, 0.5), silent = TRUE)
stopifnot(inherits(e, "try-error"))
BiCopPDF(s[, 1], s[, 2], BiCop(5, 4))
BiCopPDF(s[, 1], s[, 2], obj = BiCop(5, 4))
BiCopPDF(s[, 1], s[, 2], 2, 0.5, 4)
s[2, 1] <- NA
suppressWarnings(BiCopPDF(s[, 1], s[, 2], 3, 3))
suppressWarnings(BiCopPDF(s[, 1], s[, 2], 2, 1:10/11, 4))
suppressWarnings(BiCopPDF(0.5, 0.5, 3, 3))
suppressWarnings(BiCopPDF(NA, u2 = 0.5, family = 3, par = 3))

## BiCopSelect -------------------------------------------
s <- BiCopSim(500, 1, 0.7)
s[, 1] <- 1 - s[, 1]
e <- try(BiCopSelect(s[, 1], s[, 2], familyset = 3:4, rotations = FALSE))
stopifnot(inherits(e, "try-error"))
BiCopSelect(s[, 1], s[, 2], familyset = 3:4)
s[2, 1] <- s[10, 2] <-- NA
BiCopSelect(s[, 1], s[, 2], familyset = 1:4, selectioncrit = "BIC")
BiCopSelect(s[, 1], s[, 2], familyset = 1:4, selectioncrit = "logLik")

## BiCopSim ----------------------------------
set.seed(123)
cop <- BiCop(family = 2, par = -0.7, par2 = 4)
simdata <- BiCopSim(100, cop)
simdata <- BiCopSim(3, 1:3, 1:3/4, 4)

## BiCopTau2Par ------------------------------
tau0 <- 0.5
rho <- BiCopTau2Par(family = 1, tau = tau0)
rho <- BiCopTau2Par(family = 2, tau = c(0.4, 0.5, 0.6))
vtau <- seq(from = 0.1, to = 0.8, length.out = 100)
thetaC <- BiCopTau2Par(family = 3, tau = vtau)
thetaG <- BiCopTau2Par(family = 4, tau = vtau)
thetaF <- BiCopTau2Par(family = 5, tau = vtau)
thetaJ <- BiCopTau2Par(family = 6, tau = vtau)
plot(thetaC ~ vtau, type = "l", ylim = range(thetaF))
lines(thetaG ~ vtau, col = 2)
lines(thetaF ~ vtau, col = 3)
lines(thetaJ ~ vtau, col = 4)
theta <- BiCopTau2Par(family = c(3,4,6), tau = c(0.4, 0.5, 0.6))
BiCopPar2Tau(family = c(3,4,6), par = theta)

## BiCopVuongClarke ------------------------------------
dat <- BiCopSim(500, 2, 0.7, 5)
BiCopVuongClarke(dat[,1], dat[,2], familyset = 1:6)



#-----------------------------------------------------------------#
Matrix <- c(5, 2, 3, 1, 4,
            0, 2, 3, 4, 1,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 1)
Matrix <- matrix(Matrix, 5, 5)
family <- c(0, 1, 3, 4, 4,
            0, 0, 3, 4, 1,
            0, 0, 0, 4, 1,
            0, 0, 0, 0, 3,
            0, 0, 0, 0, 0)
family <- matrix(family, 5, 5)
par <- c(0, 0.2, 0.9, 1.5, 3.9,
         0, 0, 1.1, 1.6, 0.9,
         0, 0, 0, 1.9, 0.5,
         0, 0, 0, 0, 4.8,
         0, 0, 0, 0, 0)
par <- matrix(par, 5, 5)
par2 <- matrix(0, 5, 5)


## RVineAIC / RVineBIC ----------------------------
RVM <- RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2)
simdata <- RVineSim(10, RVM)
simdata[2, 2] <- NA
RVineAIC(simdata, RVM)
RVineBIC(simdata, RVM)
RVM$family[5, 1] <- -3
e <- try(RVineAIC(simdata, RVM))
stopifnot(inherits(e, "try-error"))
e <- try(RVineBIC(simdata, RVM))
stopifnot(inherits(e, "try-error"))
simdata[1, 1] <- 2
e <- try(RVineAIC(simdata, RVM))
stopifnot(inherits(e, "try-error"))
e <- try(RVineBIC(simdata, RVM))
stopifnot(inherits(e, "try-error"))

## RVineClarkeTest ----------------------------
data(daxreturns)
RVM <- RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2)
CVM <- RVineStructureSelect(daxreturns[,1:5], c(1:6), type = "CVine")
daxreturns[2, 2] <- NA
clarke <- RVineClarkeTest(daxreturns[,1:5], RVM, CVM)
clarke$statistic
clarke$statistic.Schwarz
clarke$p.value
clarke$p.value.Schwarz

## RVineLogLik --------------------------
RVM <- RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2)
simdata <- RVineSim(10, RVM)
simdata[2, 2] <- NA
RVineLogLik(simdata, RVM, separate = TRUE)
RVineLogLik(simdata, RVM, separate = FALSE)

## RVineCopSelect ----------------------------
RVM <- RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2)
simdata <- RVineSim(100, RVM)
RVM1 <- RVineCopSelect(simdata, familyset = c(1, 3, 4, 5 ,6), Matrix)
simdata[1:99, 2] <- NA
RVM1 <- RVineCopSelect(simdata, familyset = c(1, 3, 4, 5 ,6), Matrix)

## RVineGofTest ----------------------------------
data(daxreturns)
daxreturns[1, 2] <- NA
RVineGofTest(daxreturns[1:100,1:5], RVM, B = 0)
RVineGofTest(daxreturns[1:100,1:5], RVM, method = "ECP2",
             statistic = "CvM", B = 200)

## RVineGrad / RVineHessian / RVineStdError ----------------------------
simdata <- RVineSim(300, RVM)
RVineGrad(simdata, RVM)
RVineHessian(simdata[1,], RVM)
simdata[1, 1] <- NA
RVineGrad(simdata, RVM)
out2 <- RVineHessian(simdata, RVM)
RVineStdError(out2$hessian, RVM)

## RVineMatrix -----------------------------------------
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3", "V4", "V5"))

## RVineMLE -------------------------------------------
simdata <- RVineSim(100, RVM)
RVineMLE(simdata, RVM, grad = FALSE, trace = 0)
simdata[1:2, 2] <- NA
RVineMLE(simdata, RVM, grad = FALSE, trace = 0)

## RVinePIT -----------------------------------------
simdata <- RVineSim(10, RVM)
RVinePIT(simdata, RVM)
simdata[1:2, 2] <- NA
RVinePIT(simdata, RVM)

## RVineSeqEst --------------------------------------
simdata <- RVineSim(50, RVM)
RVineSeqEst(simdata, RVM, method = "itau", se = FALSE)
simdata[1:49, 2] <- NA
RVineSeqEst(simdata, RVM, method = "itau", se = TRUE)




