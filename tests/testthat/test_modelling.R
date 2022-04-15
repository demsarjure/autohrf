# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

library(autohrf)

# set tolerance
tol <- 0.01

# set seed
set.seed(27)

# load data
df <- swm

# run_model --------------------------------------------------------------------
model <- data.frame(event      = c("encoding", "delay", "response"),
                    start_time = c(0,           2.65,    12.5),
                    duration   = c(2.65,        9.85,    3))

# convolve
ce <- convolve_events(model)

# set some weights
roi_weights <-
  data.frame(roi = c("L_LIPv_ROI", "L_SCEF_ROI", "R_p32pr_ROI"),
             weight = c(2, 2, 4))

# run_model
test_that("run_model", {
  # run_model
  res <- run_model(df, ce, model, roi_weights = roi_weights)

  # test
  expect_equal(res$r2$mean, 0.897, tolerance = tol)
  expect_equal(res$r2$median, 0.92, tolerance = tol)
  expect_equal(res$r2$min, 0.765, tolerance = tol)
  expect_equal(res$r2$weighted, 0.9, tolerance = tol)
})

# run_model
test_that("plot_model", {
  plot <- plot_model(model, ce)
  vdiffr::expect_doppelganger("plot_model output", plot)
})


# autohrf ----------------------------------------------------------------------
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
autofit <- autohrf(df,
                  model_specs,
                  roi_weights = roi_weights,
                  population = 2,
                  iter = 2)

# autohrf
test_that("autohrf", {
  expect_equal(mean(autofit[[1]]$fitness), 0.802, tolerance = tol)
  expect_equal(mean(autofit[[2]]$fitness), 0.849, tolerance = tol)
})

# plot_best_models
test_that("plot_best_models", {
  plot <- plot_best_models(autofit)
  vdiffr::expect_doppelganger("plot_best_models output", plot)
})

# print_best_models
test_that("print_best_models", {
  expect_output(print_best_models(autofit))
})
