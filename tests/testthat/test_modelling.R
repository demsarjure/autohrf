# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

library(autohrf)

# set tolerance
tol <- 0.5

# set seed
set.seed(27)

# load data
df <- flanker

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
  data.frame(roi = c("1", "2", "3"),
             weight = c(2, 2, 4))

em <- evaluate_model(df, model, tr = 2.5, roi_weights = roi_weights)

# evaluate_model
test_that("evaluate_model", {
  expect_output(evaluate_model(df, model, tr = 2.5, roi_weights = roi_weights))
})

# evaluate_model by roi
test_that("evaluate_model_by_roi", {
  expect_equal(mean(em$by_roi$r2), 0.28, tolerance = tol)
  expect_equal(mean(em$by_roi$r2w), 0.43, tolerance = tol)
  expect_equal(mean(em$by_roi$bic), 236.62, tolerance = tol)
  expect_equal(mean(em$by_roi$bicw), 446.47, tolerance = tol)
})

# evaluate_model single event
test_that("evaluate_model_single", {
  expect_output(evaluate_model(df, model_single, tr = 2.5))
})

# plot_model
test_that("plot_model", {
  plot <- plot_model(em)
  expect_s3_class(plot, "ggplot")
})

# plot_model by roi
test_that("plot_model_by_roi", {
  plot <- plot_model(em, by_roi = TRUE)
  expect_s3_class(plot, "ggplot")
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
                   iter = 2,
                   cores = 1)

# autohrf
test_that("autohrf", {
  expect_equal(mean(autofit[[1]]$fitness), 0.26, tolerance = tol)
  expect_equal(mean(autofit[[2]]$fitness), 0.27, tolerance = tol)
})

# autohrf_parallel
test_that("autohrf_parallel", {
  autofit_p <- autohrf(df,
                       model_constraints,
                       tr = 2.5,
                       population = 2,
                       iter = 2,
                       cores = 2)
  expect_equal(mean(autofit_p[[1]]$fitness), 0.30, tolerance = tol)
  expect_equal(mean(autofit_p[[2]]$fitness), 0.32, tolerance = tol)
})

# plot_best_models
test_that("plot_best_models", {
  plot <- plot_best_models(autofit)
  expect_s3_class(plot, "ggplot")
})

# get_best_models_return_fitness
test_that("get_best_models_return_fitness", {
  expect_output(get_best_models(autofit, return_fitness = TRUE))
})
