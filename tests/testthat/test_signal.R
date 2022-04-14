# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

library(autohrf)

# set tolerance
tol <- 0.01

# downsample
test_that("downsample", {
  x <- seq(1:1000)
  x_ds <- downsample(x)
  expect_equal(mean(x_ds), 501, tolerance = tol)
})
