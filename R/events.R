# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

# A helper function for convolving events of a model with a generated HRF
# signal.
convolve_events <- function(model,
                            tr = 2.5,
                            f = 100,
                            hrf = "boynton",
                            t = 32,
                            delta = 2.25,
                            tau = 1.25,
                            alpha = 2,
                            p = c(6, 16, 1, 1, 6, 0, 32)) {

  len <- max(model$start_time) + max(model$duration) + 30
  len <- ceiling(len / tr) * tr
  e_tr <- tr / f
  len_ts <- ceiling(len / e_tr)
  n_events <- dim(model)[1]

  ts <- c(1:len_ts) * 0
  m <- matrix(0, len_ts, n_events)

  for (i in 1:n_events) {
    start <- round(model$start_time[i] / e_tr) + 1
    end <- round((model$start_time[i] + model$duration[i]) / e_tr)
    for (j in start:end) {
      ts[j] <- ts[j] + 1
      m[j, i] <- 1
    }
  }

  # hrf
  if (hrf == "spm") {
    hrf_s <- create_spm_hrf(tr = e_tr, t = t, p = p)
  } else {
    hrf_s <- create_boynton_hrf(
      tr = e_tr,
      t = t,
      delta = delta,
      tau = tau,
      alpha = alpha)
  }
  x <- downsample(convolve_hrf(m, hrf_s), f)

  # time series
  ts <- downsample(convolve_hrf(ts, hrf_s), f)

  # model
  m <- downsample(m, f)

  # return results as a list
  return(list(m = m, x = x, ts = ts))
}


# A helper function for plotting events of a fitted model.
plot_events <- function(af) {
  # init local variables for CRAN check
  start_time <- NULL
  duration <- NULL
  x <- NULL
  y <- NULL

  # prepare the data
  model <- af$models[[1]]
  ce <- af$best
  len_ts <- dim(ce$x)[1]
  time_ts <- c(0:(len_ts - 1)) * af$tr
  d <- data.frame(y = ce$ts, x = time_ts, ts = "ts")

  for (n in 1:dim(ce$x)[2]) {
    d <- rbind(d, data.frame(y = ce$x[, n], x = time_ts, ts = model$event[n]))
  }
  model$ts <- model$event

  # construct the plot
  ggplot() +
    geom_line(data = d[d$ts == "ts", ], aes(x = x, y = y),
              color = "black", size = 1, alpha = 0.75) +
    geom_line(data = d[d$ts != "ts", ],
              aes(x = x, y = y, color = ts, group = ts)) +
    geom_rect(data = model, aes(xmin = start_time,
                                xmax = start_time + duration,
                                ymax = 0,
                                ymin = -0.1,
                                color = ts,
                                fill = ts)) +
    scale_color_brewer(type = "qual", palette = "Set1", name = "Event") +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "Event") +
    ylab("") +
    xlab("Time")
}