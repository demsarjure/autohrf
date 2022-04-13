library(autohrf)

# set tolerance
tol <- 0.01

# set seed
set.seed(27)

# prepare model specs
# 3 events: encoding, delay, response
model3 <- data.frame(
  event        = c("encoding", "delay", "response"),
  start_time   = c(0,          2.65,     12.5     ),
  end_time     = c(3,          12.5,     16       ),
  min_duration = c(1,          5,        1        )
)

# 4 events: fixation, target, delay, response
model4 <- data.frame(
  event        = c("fixation", "target", "delay", "response"),
  start_time   = c(0,          2.5,      2.65,    12.5),
  end_time     = c(2.5,        3,        12.5,    15.5),
  min_duration = c(1,          0.1,      5,       1)
)

model_specs <- list(model3, model4)

# run autohrf
df <- swm
autofit <- autohrf(df, model_specs, population=2, iter=2)

# convolve_events
test_that("convolve_events", {
  # create the model
  m <- data.frame(event = c("encoding", "delay", "response"),
                  start_time = c(0, 2.5, 12.5), duration = c(2.5, 10, 5))

  # convolve
  r <- convolve_events(m)

  # test
  expect_equal(mean(r$m), 0.111, tolerance=tol)
  expect_equal(mean(r$x), 0.135, tolerance=tol)
  expect_equal(mean(r$ts), 0.333, tolerance=tol)
})

# plot_events
test_that("plot_events", {
  plot <- plot_events(autofit[[1]])
  vdiffr::expect_doppelganger("plot_events output", plot)
})
