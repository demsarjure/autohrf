# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

library(autohrf)

# set tolerance
tol <- 0.01

# set seed
set.seed(27)

# prepare model specs
# 3 events: encoding, delay, response
model3 <- data.frame(
  event        = c("encoding", "delay", "response"),
  start_time   = c(0,          2.65,     12.5),
  end_time     = c(3,          12.5,     16)
)

# 4 events: fixation, target, delay, response
model4 <- data.frame(
  event        = c("fixation", "target", "delay", "response"),
  start_time   = c(0,          2.5,      2.65,    12.5),
  end_time     = c(2.5,        3,        12.5,    15.5)
)

model_specs <- list(model3, model4)

# run autohrf
df <- swm
autofit <- autohrf(df, model_specs, tr = 2.5, population = 2, iter = 2)

# convolve_events
test_that("convolve_events", {
  # create the model
  m <- data.frame(event = c("encoding", "delay", "response"),
                  start_time = c(0, 2.5, 12.5), duration = c(2.5, 10, 5))

  # convolve
  r <- convolve_events(m, tr = 2.5)

  # test
  expect_equal(mean(r$m), 0.111, tolerance = tol)
  expect_equal(mean(r$x), 0.124, tolerance = tol)
  expect_equal(mean(r$ts), 0.291, tolerance = tol)
})

# plot_events
test_that("plot_events", {
  plot <- plot_events(autofit[[1]])
  vdiffr::expect_doppelganger("plot_events output", plot)
})
