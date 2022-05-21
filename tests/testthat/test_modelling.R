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

# evaluate_model ---------------------------------------------------------------
# 3 events
model <- data.frame(event      = c("encoding", "delay", "response"),
                    start_time = c(0,           2.65,    12.5),
                    duration   = c(2.65,        9.85,    3))

# single event
model_single <- data.frame(event      = c("encoding"),
                           start_time = c(0),
                           duration   = c(3))

# set some weights
roi_weights <-
  data.frame(roi = c("L_LIPv_ROI", "L_SCEF_ROI", "R_p32pr_ROI"),
             weight = c(2, 2, 4))

em <- evaluate_model(df, model, tr = 2.5, roi_weights = roi_weights)

# evaluate_model
test_that("evaluate_model", {
  expect_output(evaluate_model(df, model, tr = 2.5, roi_weights = roi_weights))
})

# evaluate_model single event
test_that("evaluate_model_single", {
  expect_output(evaluate_model(df, model_single, tr = 2.5))
})


# plot_model
test_that("plot_model", {
  plot <- plot_model(em)
  vdiffr::expect_doppelganger("plot_model output", plot)
})

# plot_model by roi
test_that("plot_model_by_roi", {
  plot <- plot_model(em, by_roi = TRUE)
  vdiffr::expect_doppelganger("plot_model_by_roi output", plot)
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

model_constraints <- list(model3, model4)

# run autohrf
autofit <- autohrf(df,
                   model_constraints,
                   tr = 2.5,
                   roi_weights = roi_weights,
                   population = 2,
                   iter = 2)

# autohrf
test_that("autohrf", {
  expect_equal(mean(autofit[[1]]$fitness), 0.909, tolerance = tol)
  expect_equal(mean(autofit[[2]]$fitness), 0.905, tolerance = tol)
})

# plot_best_models
test_that("plot_best_models", {
  plot <- plot_best_models(autofit)
  vdiffr::expect_doppelganger("plot_best_models output", plot)
})

# get_best_models
test_that("get_best_models", {
  expect_output(get_best_models(autofit))
})