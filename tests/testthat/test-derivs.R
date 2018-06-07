context("rotations")

test_that("90 deg rotation is handled correctly in evaluators", {
    ## model setup
    u1 <- 0.4
    u2 <- 0.7
    cop <- BiCop(family = 23, par = -3)
    cop_unrot <- BiCop(family = 3, par = 3)
    theta <- c(u1, u2, -3)

    ## pdf
    expect_equal(
        BiCopPDF(u1, u2, cop),
        BiCopPDF(1 - u1, u2, cop_unrot),
        label = "pdf"
    )

    ## h-functions
    expect_equal(
        BiCopHfunc1(u1, u2, cop),
        BiCopHfunc1(1 - u1, u2, cop_unrot),
        label = "hfunc1"
    )
    expect_equal(
        BiCopHfunc2(u1, u2, cop),
        1 - BiCopHfunc2(1 - u1, u2, cop_unrot),
        label = "hfunc1"
    )

    ## 1st pdf derivative
    pdf_fun <- function(theta) {
        cop$par <- theta[3]
        BiCopPDF(theta[1], theta[2], cop)
    }
    grad <- numDeriv::grad(pdf_fun, theta)
    derivs_1st <- list("u1", "u2", "par")
    for (i in seq_along(derivs_1st) ) {
        expect_equal(
            BiCopDeriv(theta[1], theta[2], cop, deriv = derivs_1st[[i]]),
            grad[i],
            label = paste("1st pdf derivative w.r.t.", derivs_1st[[i]]),
            tol = 1e-2
        )
    }

    ## 1st hfunc derivative
    for (deriv in c("u2", "par")) {
        expect_equal(
            BiCopHfuncDeriv(u1, u2, cop, deriv = deriv),
            numDeriv::grad(BiCopHfunc2, u1, u2 = u2, obj = cop),
            label = paste("1st hfunc derivative w.r.t.", deriv)
        )
    }


    derivs_2nd <- list(
        c("u1", "u1"),
        c("u2", "u2"),
        c("par1", "u1"),
        c("par1", "u2")
    )

    ## 2nd pdf derivates
    hess <- numDeriv::hessian(pdf_fun, c(u1, u2, -3))
    colnames(hess) <- rownames(hess) <- c("u1", "u2", "par1")

    for (deriv in derivs_2nd) {
        expect_equal(
            BiCopDeriv2(u1, u2, cop, deriv = paste(unique(deriv), collapse = "")),
            hess[deriv[1], deriv[2]],
            label = paste("2nd pdf derivative w.r.t.", deriv, collapse = ", ")
        )
    }

    ## 2nd hfunc derivatives
    hfunc_fun <- function(theta) {
        cop$par <- theta[3]
        BiCopHfunc2(theta[1], theta[2], cop)
    }
    hess <- numDeriv::hessian(hfunc_fun, c(u1, u2, -3))
    colnames(hess) <- rownames(hess) <- c("u1", "u2", "par1")

    for (deriv in derivs_2nd[c(2, 4)]) {
        expect_equal(
            BiCopHfuncDeriv2(u1, u2, cop, deriv = paste(unique(deriv), collapse = "")),
            hess[deriv[1], deriv[2]],
            label = paste("2nd pdf derivative w.r.t.", deriv, collapse = ", ")
        )
    }

})
