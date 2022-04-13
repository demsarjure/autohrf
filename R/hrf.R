# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

# A helper function for creating a Boynton HRF.
create_boynton_hrf <- function(tr=2.5, t=32, delta=2.25, tau=1.25, alpha=2) {
  t <- c(0:(t/tr)) * tr
  r <- (t - delta) / tau
  h <- ((r ^ alpha) * exp(-r))
  h[t < delta] <- 0
  h <- h / ((alpha ^ alpha) * exp(-alpha))
  return (h)
}

# A helper function for creating SPM HRF.
create_spm_hrf <- function(tr=2.5, t=32, p=c(6, 16, 1, 1, 6, 0, 32)) {
  dt <- tr/t
  u <- c(0:(p[7]/dt)) - p[6]/dt;
  hrf <- dgamma(u*(dt/p[3]), p[1]/p[3]) - dgamma(u*(dt/p[4]), p[2]/p[4])/p[5]
  hrf <- hrf[c(0:(p[7]/tr)) * t + 1];
  hrf <- hrf / max(hrf)
  return(hrf)
}


# A helper function for convolving HRF with a signal.
convolve_hrf <- function(x, hrf_s) {
  pad <- length(hrf_s)+20
  if (is.matrix(x)) {
    m  <- matrix(0, dim(x)[1]+2*pad, dim(x)[2])
    m[(pad+1):(pad+dim(x)[1]),] <- x
    hrf_s <- filter(m, filter=hrf_s, method="convolution", sides=1)[(pad+1):(pad+dim(x)[1]),]
    maxv <- apply(hrf_s, 2, FUN=function(x) max(abs(x)))
    hrf_s <- hrf_s / matrix(maxv, nrow=dim(hrf_s)[1], ncol=dim(hrf_s)[2], byrow=TRUE)
  } else {
    hrf_s <- filter(c(c(1:pad)*0, x, c(1:pad)*0), filter=hrf_s, method="convolution", sides=1)
    hrf_s <- hrf_s[(pad+1):(pad+length(x))]
    hrf_s <- hrf_s / max(abs(hrf_s))
  }
  return(hrf_s)
}
