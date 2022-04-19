# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

# A helper function for downsampling a given signal.
downsample <- function(x, f=100, method="middle") {
  # prepare the input data
  x <- as.matrix(x)
  l <- dim(x)[1]
  n <- ceiling(l / f)
  m <- matrix(0, n, dim(x)[2])

  # iterate over the series
  for (i in 1:n) {
    start <- (i - 1) * f + 1
    end <- i * f
    if (end > l) end <- l

    # middle or mean
    if (method == "middle") {
      middle <- round(start + f / 2)
      if (middle > l) middle <- l
      m[i, ] <- x[middle, ]
    }
    else if (method == "mean") {
      m[i, ] <- apply(as.matrix(x[start:end, ]), 2, FUN = mean)
    }
  }

  # return the downsampled results
  return(m)
}
