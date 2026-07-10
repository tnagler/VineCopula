test_that("vine derivatives handle 90 and 270 degree rotations", {
    RVM <- D2RVine(
        1:4,
        family = c(23, 24, 26, 33, 34, 36),
        par = c(-2.5, -1.5, -2, -2.5, -1.5, -2)
    )
    data <- matrix(c(0.2, 0.4, 0.65, 0.8), nrow = 1)
    par_index <- lower.tri(RVM$par)
    parameters <- RVM$par[par_index]

    loglik <- function(parameters) {
        par <- RVM$par
        par[par_index] <- parameters
        RVineLogLik(data, RVM, par = par, separate = FALSE)$loglik
    }

    expected_gradient <- numDeriv::grad(loglik, parameters)
    expected_hessian <- numDeriv::hessian(loglik, parameters)
    hessian <- RVineHessian(data, RVM)$hessian
    reverse_order <- nrow(hessian):1

    expect_equal(
        RVineGrad(data, RVM)$gradient,
        expected_gradient,
        tolerance = 1e-5
    )
    expect_equal(
        hessian[reverse_order, reverse_order],
        expected_hessian,
        tolerance = 1e-5
    )
})
