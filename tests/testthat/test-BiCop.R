context("C++ class BiCop")

call_BiCops <- function(fun_name, family, par, par2) {
    .C(fun_name,
       family = as.integer(family),
       par = as.double(par),
       par2 = as.double(par2),
       PACKAGE = "VineCopula")
}

test_that("default constructor works", {
    obj <- call_BiCops("test_BiCop", -1, 234, 2341)
    expect_equal(obj$family, 0)
    expect_equal(obj$par, 0)
    expect_equal(obj$par2, 0)
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
n <- 10
u1 <- runif(n)
u2 <- runif(n)

## check that asymmetries are handled correctly by using a family with
## asymmetric tails and negative dependence
family <- 23
par <- -1
par2 <- 0
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
