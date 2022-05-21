# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

# A helper function for downsampling a given signal.
downsample <- function(y, f = 100) {
  # prepare the input data
  y <- as.matrix(y)
  l <- dim(y)[1]
  n <- ceiling(l / f)
  m <- matrix(0, n, dim(y)[2])

  # iterate over the series
  for (i in 1:n) {
    # downsample
    start <- (i - 1) * f + 1
    end <- i * f
    if (end > l) end <- l
    middle <- round(start + f / 2)
    if (middle > l) middle <- l
    m[i, ] <- y[middle, ]
  }

  # return the downsampled results
  return(m)
}
