# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

#' @title convolve_events
#' @description A helper function for convolving events of a model with a
#' generated HRF signal.
#' @export
#'
#' @param model A data frame containing information about the model to use
#' and its events (event, start_time and duration).
#' @param tr MRI's repetition time.
#' @param max_duration Maximum duration of the signal.
#' @param hrf Method to use for HRF generation, can be "boynton" or "spm".
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param p_boynton Parameters for the Boynton's HRF.
#' @param p_spm Parameters for the SPM HRF.
#' @param f Upsampling factor.
#'
#' @return Returns a list with the convolved signal and time series.
convolve_events <- function(model,
                            tr,
                            max_duration,
                            hrf = "spm",
                            t = 32,
                            p_boynton = c(2.25, 1.25, 2),
                            p_spm = c(6, 16, 1, 1, 6, 0),
                            f = 100) {

  len <- max(model$start_time) + max(model$duration) + max_duration
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
    hrf_s <- autohrf::create_spm_hrf(tr = e_tr, t = t, p = p_spm)
  } else {
    hrf_s <- autohrf::create_boynton_hrf(
      tr = e_tr,
      t = t,
      p = p_boynton)
  }

  # signal
  y <- autohrf::downsample(autohrf::convolve_hrf(m, hrf_s), f)

  # time series
  ts <- autohrf::downsample(autohrf::convolve_hrf(ts, hrf_s), f)

  # return results as a list
  return(list(y = y, ts = ts))
}


#' @title plot_events
#' @description A helper function for plotting events of a fitted model.
#' @export
#'
#' @param af The output from the autohrf function.
#' @param i Model index.
#'
#' @return Returns a plot of the events.
plot_events <- function(af, i = NULL) {
  # init local variables for CRAN check
  start_time <- NULL
  duration <- NULL
  x <- NULL
  y <- NULL

  # prepare the data
  model <- af$models[[1]]
  ce <- af$ce
  len_ts <- dim(ce$y)[1]
  time_ts <- c(0:(len_ts - 1)) * af$tr
  d <- data.frame(y = ce$ts, x = time_ts, ts = "ts")

  for (n in 1:dim(ce$y)[2]) {
    d <- rbind(d, data.frame(y = ce$y[, n], x = time_ts, ts = model$event[n]))
  }

  # convert events to factors
  d$ts <- factor(d$ts, levels = c(model$event, "ts"))
  model$event <- factor(model$event, levels = model$event)
  model$ts <- model$event

  # construct the plot
  p <- ggplot() +
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

  if (!is.null(i))
    p <- p + ggtitle(paste0("Model ", i))

  p
}
