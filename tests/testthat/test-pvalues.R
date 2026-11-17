test_that("zero and complete exceedance counts stay strictly inside the endpoints", {
  set.seed(123)
  x <- rlnorm(40, 1, 0.5)
  # These seeded, one-replicate runs give zero or complete exceedances for
  # every statistic, exercising both real bootstrap paths without mocks.
  cases <- data.frame(method = c("mult", "mult", "pb", "pb"),
                      seed = c(16, 2, 2, 3),
                      count = c(0, 1, 0, 1),
                      expected = c(0.25, 0.75, 0.25, 0.75))
  for (i in seq_len(nrow(cases))) {
    set.seed(cases$seed[i])
    fit <- gof.test(x, "lnorm", 1, list(meanlog = 1, sdlog = 0.5),
                    simulation = cases$method[i], ng = 17)
    counts <- rowSums(fit$stat.sim >= fit$statistic)
    expect_equal(unname(counts), rep(cases$count[i], 4))
    expect_identical(unname(fit$p.value), rep(cases$expected[i], 4))
    expect_true(all(fit$p.value > 0 & fit$p.value < 1))
  }
})
