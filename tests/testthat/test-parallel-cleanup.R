context("Parallel cleanup")

## The functions that register a socket cluster call rm(list = ls()) before
## returning.  Because ls() does not list dot-prefixed names, a cluster stored
## in a plain variable is removed before the on.exit handler runs; the handler
## then fails on the lookup, try(silent = TRUE) hides it, and the workers stay
## alive for the rest of the session.  Repeated calls used to accumulate
## `cores` worker processes each.

n_workers <- function() {
    out <- suppressWarnings(
        system("pgrep -u $USER -f workRSOC[K]", intern = TRUE,
               ignore.stderr = TRUE)
    )
    length(out)
}

test_that("fitting with cores > 1 leaves no worker processes behind", {
    skip_on_cran()
    skip_on_os("windows")
    if (Sys.which("pgrep") == "") {
        skip("pgrep not available")
    }

    set.seed(5)
    u <- matrix(runif(400 * 4), 400, 4)
    fams <- c(0, 1, 3, 4, 5)

    before <- n_workers()
    RVM <- RVineStructureSelect(u, familyset = fams, method = "itau",
                                indeptest = FALSE, progress = FALSE,
                                cores = 2)
    Sys.sleep(1)
    expect_equal(n_workers(), before)

    RVineCopSelect(u, familyset = fams, Matrix = RVM$Matrix,
                   indeptest = FALSE, cores = 2)
    Sys.sleep(1)
    expect_equal(n_workers(), before)

    RVineSeqEst(u, RVM, method = "itau", cores = 2)
    Sys.sleep(1)
    expect_equal(n_workers(), before)
})

test_that("no stale default cluster is left registered", {
    skip_on_cran()

    set.seed(6)
    u <- matrix(runif(300 * 3), 300, 3)
    RVineStructureSelect(u, familyset = c(0, 1, 5), method = "itau",
                         indeptest = FALSE, progress = FALSE, cores = 2)

    cl <- tryCatch(parallel::getDefaultCluster(), error = function(e) NULL)
    expect_null(cl)
})

test_that("the tree criterion does not capture the caller's frame", {
    ## set_treecrit() returns a closure carrying its own frame as environment.
    ## Only the AIC/BIC branches read `famset`, so under the default "tau" the
    ## promise went unforced and kept a live reference to the frame of
    ## RVineStructureSelect() -- data matrix included. The closure is shipped to
    ## every worker on every parallel dispatch, so it must stay small.
    tc <- VineCopula:::set_treecrit("tau", famset = 1:40)
    expect_lt(length(serialize(tc, NULL)), 100000)
})
