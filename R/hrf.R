# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

#' @title create_boynton_hrf
#' @description Generate an HRF using the Boynton method.
#' @export
#'
#' @param tr MRI's repetition time.
#' @param t TODO
#' @param delta the delta parameter of Boynton's HRF.
#' @param tau the tau parameter of Boynton's HRF.
#' @param alpha the alpha parameter of Boynton's HRF.
#'
#' @return TODO.
#'
#' @examples
#' # create the HRF with a smaller tr
#' x <- create_boynton_hrf(tr=1)
#'
create_boynton_hrf <- function(tr=2.5, t=32, delta=2.25, tau=1.25, alpha=2) {
  t <- c(0:(t/tr)) * tr
  r <- (t - delta) / tau
  h <- ((r ^ alpha) * exp(-r))
  h[t < delta] <- 0
  h <- h / ((alpha ^ alpha) * exp(-alpha))
  return (h)
}


#' @title create_spm_hrf
#' @description Generate an HRF using the SPM method.
#' @import stats
#' @export
#'
#' @param tr MRI's repetition time.
#' @param t TODO
#' @param p TODO
#'
#' @return TODO.
#'
#' @examples
#' # create the HRF with a smaller tr
#' x <- create_spm_hrf(tr=1)
#'
create_spm_hrf <- function(tr=2.5, t=32, p=c(6, 16, 1, 1, 6, 0, 32)) {
  dt <- tr/t
  u <- c(0:(p[7]/dt)) - p[6]/dt;
  hrf <- dgamma(u*(dt/p[3]), p[1]/p[3]) - dgamma(u*(dt/p[4]), p[2]/p[4])/p[5]
  hrf <- hrf[c(0:(p[7]/tr)) * t + 1];
  hrf <- hrf / max(hrf)
  return(hrf)
}


#' @title convolve_hrf
#' @description Convolve HRF with a signal.
#' @export
#'
#' @param x the input signal.
#' @param hrf HRF to convolve the signal with.
#'
#' @return TODO.
#'
#' @examples
#' # create the HRF
#' hrf <- create_boynton_hrf(tr=1)
#'
#' # our signal
#' x <- sin(1:1000)
#'
#' ts <- convolve_hrf(x, hrf)
#'
convolve_hrf <- function(x, hrf) {
  pad <- length(hrf)+20
  if (is.matrix(x)) {
    m  <- matrix(0, dim(x)[1]+2*pad, dim(x)[2])
    m[(pad+1):(pad+dim(x)[1]),] <- x
    hrf <- filter(m, filter=hrf, method="convolution", sides=1)[(pad+1):(pad+dim(x)[1]),]
    maxv <- apply(hrf, 2, FUN=function(x) max(abs(x)))
    hrf <- hrf / matrix(maxv, nrow=dim(hrf)[1], ncol=dim(hrf)[2], byrow=TRUE)
  } else {
    hrf <- filter(c(c(1:pad)*0, x, c(1:pad)*0), filter=hrf, method="convolution", sides=1)
    hrf <- hrf[(pad+1):(pad+length(x))]
    hrf <- hrf / max(abs(hrf))
  }
  return(hrf)
}
