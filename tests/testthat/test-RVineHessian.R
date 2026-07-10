test_that("RVineHessian uses all observations", {
    Matrix <- matrix(c(2, 1, 0, 1), 2, 2)
    family <- par <- par2 <- matrix(0, 2, 2)
    family[2, 1] <- 1
    par[2, 1] <- 0.5
    RVM <- RVineMatrix(Matrix, family, par, par2)
    data <- matrix(c(0.1, 0.3, 0.7, 0.2, 0.8, 0.4), ncol = 2)

    per_observation <- lapply(seq_len(nrow(data)), function(i) {
        RVineHessian(data[i, ], RVM)$hessian
    })

    expect_equal(
        RVineHessian(data, RVM)$hessian,
        Reduce(`+`, per_observation)
    )
})
