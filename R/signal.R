#' @title downsample
#' @description Downsample a time series by filtering its time stamps/indexes.
#' @export
#' @param x The indexes/time stamps to downsample.
#' @param f Downsampling frequency.
#' @param method Can be "middle" or "mean". Middle will return integer results, mean will return floats.
#' @return A vector representing a downsampled input signal.
#'
#' @examples
#' # time stamps
#' x <- seq(1:1000)
#'
#' # downsample
#' x_d <- downsample(x)
#'
downsample <- function(x, f=100, method="middle") {
  # prepare the input data
  x <- as.matrix(x)
  l <- dim(x)[1]
  n <- ceiling(l / f)
  m <- matrix(0, n, dim(x)[2])

  # iterate over the series
  for (i in 1:n) {
    start = (i-1)*f + 1
    end   = i*f
    if (end > l) end <- l

    # middle or mean
    if (method == "middle") {
      middle <- round(start + f/2)
      if (middle > l) middle <- l
      m[i,] = x[middle, ]
    }
    else if (method == "mean") {
      m[i,] = apply(as.matrix(x[start:end,]), 2, FUN=mean)
    }
  }

  # return the downsampled results
  return(m)
}
