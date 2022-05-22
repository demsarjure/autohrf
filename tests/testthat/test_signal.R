# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

library(autohrf)

# set tolerance
tol <- 0.05

# downsample
test_that("downsample", {
  y <- seq(1:1000)
  y_ds <- downsample(y)
  expect_equal(mean(y_ds), 501, tolerance = tol)
})
