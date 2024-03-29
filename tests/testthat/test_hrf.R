# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

library(autohrf)

# set tolerance
tol <- 0.05

# create_boynton_hrf
test_that("create_boynton_hrf", {
  hrf <- create_boynton_hrf(tr = 2.5)
  expect_equal(mean(hrf), 0.133, tolerance = tol)
})

# create_spm_hrf
test_that("create_spm_hrf", {
  hrf <- create_spm_hrf(tr = 2.5)
  expect_equal(mean(hrf), 0.147, tolerance = tol)
})

# convolve_hrf
test_that("convolve_hrf", {
  hrf <- create_boynton_hrf(tr = 2.5)
  y <- seq(1:10)
  c_hrf <- convolve_hrf(y, hrf)
  expect_equal(mean(c_hrf), 0.43, tolerance = tol)
})
