context("C++ class BiCop")

call_BiCops <- function(fun_name, family, par, par2) {
    .C(fun_name,
       family = as.integer(family),
       par = as.double(par),
       par2 = as.double(par2),
       npars = as.integer(0),
       PACKAGE = "VineCopula")
}

test_that("default constructor works", {
    obj <- call_BiCops("test_BiCop", -1, 234, 2341)
    expect_equal(obj$family, 0)
    expect_equal(obj$par, 0)
    expect_equal(obj$par2, 0)
})

test_that("number of parameters are calculated correctly", {
    objs <- list(
        BiCop(0, 0),
        BiCop(1, 0.5),
        BiCop(2, -0.5, 4),
        BiCop(33, -0.5),
        BiCop(24, -2),
        BiCop(5, 0.5),
        BiCop(16, 2),
        BiCop(7, 0.5, 2),
        BiCop(18, 0.5, 2),
        BiCop(28, -2, -2),
        BiCop(39, -2, -2),
        BiCop(10, 2, 0.5),
        BiCop(114, 2, 0.4),
        BiCop(224, -2, 0.4)
    )
    lapply(
        objs,
        function(x) {
            expect_equal(
                call_BiCops("test_BiCop", x$family, x$par, x$par2)$npars,
                x$npars
            )
        }
    )
})

test_that("parametrized constructor works", {
    obj <- call_BiCops("test_BiCop", 1, 0.3, 7.0)
    expect_equal(obj$family, 1)
    expect_equal(obj$par, 0.3)
    expect_equal(obj$par2, 7.0)
})

test_that("setter works", {
    obj <- call_BiCops("test_BiCop_set", 1, 0.3, 7.0)
    expect_equal(obj$family, 1)
    expect_equal(obj$par, 0.3)
    expect_equal(obj$par2, 7.0)
})


library(VineCopula)
## check that asymmetries are handled correctly by using a family with
## asymmetric tails and negative dependence
set.seed(5)
family <- 23
n <- 10
u <- BiCopSim(n, 23, -1)
u1 <- u[, 1]
u2 <- u[, 2]
fit <- BiCopEst(u[, 1], u[, 2], 23)
par <- fit$par
par2 <- fit$par2

call_BiCop_funcs <- function(fun_name, u1, u2, n, family, par, par2) {
    .C(fun_name,
       u1 = as.double(u1),
       u2 = as.double(u2),
       out = as.double(rep(0, n)),
       n = as.integer(n),
       family = as.integer(family),
       par = as.double(par),
       par2 = as.double(par2),
       PACKAGE = "VineCopula")
}

test_that("hFunc1 works", {
    obj <- call_BiCop_funcs("test_BiCop_hFunc1", u1, u2, n, family, par, par2)
    expect_equal(obj$out, BiCopHfunc1(u1, u2, family, par))
})

test_that("hFunc2 works", {
    obj <- call_BiCop_funcs("test_BiCop_hFunc2", u1, u2, n, family, par, par2)
    expect_equal(obj$out, BiCopHfunc2(u1, u2, family, par))
})

test_that("PDF works", {
    obj <- call_BiCop_funcs("test_BiCop_PDF", u1, u2, n, family, par, par2)
    expect_equal(obj$out, BiCopPDF(u1, u2, family, par))
})

test_that("logLik works", {
    obj <- call_BiCop_funcs("test_BiCop_logLik", u1, u2, n, family, par, par2)
    expect_equal(obj$out[1], sum(log(BiCopPDF(u1, u2, family, par))))
})

test_that("AIC works", {
    obj <- call_BiCop_funcs("test_BiCop_AIC", u1, u2, n, family, par, par2)
    expect_equal(obj$out[1], fit$AIC)
})

test_that("BIC works", {
    obj <- call_BiCop_funcs("test_BiCop_BIC", u1, u2, n, family, par, par2)
    expect_equal(obj$out[1], fit$BIC)
})

