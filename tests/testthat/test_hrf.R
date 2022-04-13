library(autohrf)

# create_boynton_hrf
test_that("create_boynton_hrf", {
  hrf <- create_boynton_hrf()
  expect_equal(mean(hrf), 0.1328353)
})

# create_spm_hrf
test_that("create_spm_hrf", {
  hrf <- create_spm_hrf()
  expect_equal(mean(hrf), 0.1467475)
})

# convolve_hrf
test_that("convolve_hrf", {
  hrf <- create_boynton_hrf()
  x <- seq(1:10)
  c_hrf <- convolve_hrf(x, hrf)
  expect_equal(mean(c_hrf), 0.4297225)
})

