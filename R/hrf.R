# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

#' @title create_boynton_hrf
#' @description A helper function for creating a Boynton HRF.
#' @export
#'
#' @param tr MRI's repetition time.
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param p Parameters for the Boynton's HRF.
#'
#' @return Returns a Boynton HRF function.
create_boynton_hrf <- function(tr,
                               t = 32,
                               p = c(2.25, 1.25, 2)) {
  # Boynton: p[1] = delta, p[2] = tau, p[3] = alpha
  y <- c(0:(t / tr)) * tr
  r <- (y - p[1]) / p[2]
  h <- ((r ^ p[3]) * exp(-r))
  h[y < p[1]] <- 0
  h <- h / ((p[3] ^ p[3]) * exp(-p[3]))
  return(h)
}

#' @title create_boynton_hrf
#' @description A helper function for creating a SPM HRF.
#' @export
#'
#' @param tr MRI's repetition time.
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param p Parameters for the SPM HRF.
#'
#' @return Returns a SPM HRF function.
create_spm_hrf <- function(tr,
                           t = 32,
                           p = c(6, 16, 1, 1, 6, 0)) {
  dt <- tr / t
  u <- (c(0:(t / tr)) * t) - p[6] / dt
  h <- dgamma(u * (dt / p[3]), p[1] / p[3]) -
         dgamma(u * (dt / p[4]), p[2] / p[4]) / p[5]
  h <- h / max(h)
  return(h)
}

#' @title convolve_hrf
#' @description A helper function for convolving HRF with a signal.
#' @export
#'
#' @param y The signal.
#' @param hrf_s The HRF.
#' @return Returns the convolution between HRF and the signal.
convolve_hrf <- function(y, hrf_s) {
  hrf_s <- as.matrix(hrf_s)
  pad <- length(hrf_s) + 20
  if (is.matrix(y)) {
    m  <- matrix(0, dim(y)[1] + 2 * pad, dim(y)[2])
    m[(pad + 1):(pad + dim(y)[1]), ] <- y
    hrf_f <- stats::filter(m, filter = hrf_s, method = "convolution", sides = 1)
    hrf_s <- as.matrix(hrf_f[(pad + 1):(pad + dim(y)[1]), ])
    maxv <- apply(hrf_s, 2, FUN = function(x) max(abs(x)))
    hrf_s <- hrf_s /
      matrix(maxv, nrow = dim(hrf_s)[1], ncol = dim(hrf_s)[2], byrow = TRUE)
  } else {
    hrf_s <- stats::filter(c(c(1:pad) * 0,
                    y,
                    c(1:pad) * 0),
                    filter = hrf_s,
                    method = "convolution",
                    sides = 1)
    hrf_s <- hrf_s[(pad + 1):(pad + length(y))]
    hrf_s <- hrf_s / max(abs(hrf_s))
  }
  return(hrf_s)
}
