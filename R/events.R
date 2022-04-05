# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

#' @title convolve_events
#' @description TODO.
#' @export
#'
#' @param model TODO.
#' @param tr MRI's repetition time.
#' @param method TODO.
#' @param f TODO.
#' @param hrf TODO.
#' @param t TODO.
#' @param delta TODO.
#' @param tau TODO.
#' @param alpha TODO.
#' @param p TODO.
#'
#' @return TODO.
#'
#' @examples
#' # TODO
#' x <- "TODO"
#'
convolve_events <- function(model,
                            tr=2.5,
                            method="middle",
                            f=100,
                            hrf="boynton",
                            t=32,
                            delta=2.25, tau=1.25, alpha=2,
                            p=c(6, 16, 1, 1, 6, 0, 32)) {

  len <- max(model$time)+max(model$duration)+30
  len <- ceiling(len/tr)*tr
  e_tr <- tr/f
  len_ts <- ceiling(len/e_tr)
  n_events <- dim(model)[1]

  ts <- c(1:len_ts)*0
  m <- matrix(0, len_ts, n_events)

  for (i in 1:n_events) {
    start <- round(model$time[i]/e_tr) + 1
    end <- round((model$time[i] + model$duration[i])/e_tr)
    for (j in start:end) {
      ts[j] <- ts[j] + model$value[i]
      m[j,i] <- model$value[i]
    }
  }

  # hrf
  if (hrf == "spm") {
    hrf_s <- create_spm_hrf(tr=e_tr, t=t, p=p)
  } else {
    hrf_s <- create_boynton_hrf(tr=e_tr, t=t, delta=delta, tau=tau, alpha=alpha)
  }
  x <- downsample(convolve_hrf(m, hrf_s), f, method)

  # time series
  ts <- downsample(convolve_hrf(ts, hrf_s), f, method)

  # model
  m <- downsample(m, f, method)

  # return results as a list
  return(list(m=m, x=x, ts=ts))
}


#' @title plot_events
#' @description TODO.
#' @import ggplot2
#' @export
#'
#' @param autofit TODO.
#'
#' @return TODO.
#'
#' @examples
#' # TODO
#' x <- "TODO"
#'
plot_events <- function(autofit) {
  # init local variables for CRAN check
  duration <- NULL
  event <- NULL
  value <- NULL
  x <- NULL
  y <- NULL

  # prepare the data
  model <- autofit$models[[1]]
  ce <- autofit$best
  len_ts <- dim(ce$x)[1]
  time_ts <- c(0:(len_ts-1)) * autofit$tr
  d <- data.frame(y=ce$ts, x=time_ts, ts="ts")

  for (n in 1:dim(ce$x)[2]) {
    d <- rbind(d, data.frame(y=ce$x[,n], x=time_ts, ts=model$event[n]))
  }
  model$ts <- model$event

  # construct the plot
  ggplot() +
    geom_line(data=d[d$ts=="ts",], aes(x=x,y=y), color="black", size=1, alpha=0.5) +
    geom_line(data=d[d$ts != "ts",], aes(x=x, y=y, color=ts, group=ts)) +
    geom_rect(data=model, aes(xmin=time, xmax=time+duration, ymax=0, ymin=-0.1*value, color=ts, fill=ts)) +
    scale_color_brewer(type="qual", palette="Set1", name="Event") +
    scale_fill_brewer(type="qual", palette="Set1", name="Event") +
    ylab("") +
    xlab("Time")
}

