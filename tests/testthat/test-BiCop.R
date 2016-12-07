context("C++ class BiCop")

call_bicops <- function(fun_name, family, par, par2) {
    .C(fun_name,
       family = as.integer(family),
       par = as.double(par),
       par2 = as.double(par2),
       PACKAGE = "VineCopula")
}

test_that("default constructor works", {
    obj <- call_bicops("test_bicop", -1, 234, 2341)
    expect_equal(obj$family, 0)
    expect_equal(obj$par, 0)
    expect_equal(obj$par2, 0)
})

test_that("parametrized constructor works", {
    obj <- call_bicops("test_bicop", 1, 0.3, 7.0)
    expect_equal(obj$family, 1)
    expect_equal(obj$par, 0.3)
    expect_equal(obj$par2, 7.0)
})

test_that("setter works", {
    obj <- call_bicops("test_bicop_set", 1, 0.3, 7.0)
    expect_equal(obj$family, 1)
    expect_equal(obj$par, 0.3)
    expect_equal(obj$par2, 7.0)
})

call_bicop_funcs <- function(fun_name, u1, u2, n, family, par, par2) {
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

test_that("hfunc1 works", {
    n <- 10
    u1 <- runif(n)
    u2 <- runif(n)
    obj <- call_bicop_funcs("test_bicop_hfunc1", u1, u2, n, 3, 1, 0)
    expect_equal(obj$out, BiCopHfunc1(u1, u2, 3, 1))
    expect_length(obj$u1, n)
    expect_length(obj$u2, n)
    expect_length(obj$out, n)
})


library(VineCopula)

test_that("hfunc2 works", {
    n <- 10
    u1 <- runif(n)
    u2 <- runif(n)
    obj <- call_bicop_funcs("test_bicop_hfunc2", u1, u2, n, 3, 1, 0)
    expect_equal(obj$out, BiCopHfunc2(u1, u2, 3, 1))
    expect_length(obj$u1, n)
    expect_length(obj$u2, n)
    expect_length(obj$out, n)
})
