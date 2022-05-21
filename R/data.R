# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

#' Datasets for autohrf examples
#' Example datasets for use in \pkg{autohrf} examples and vignettes.
#' The datasets were extracted from the internal
#' Mind and Brain Lab's (MBLab, \url{http://www.mblab.si} repository.
#' MBLab is a research lab at the
#' Faculty of Arts, Department of Psychology,
#' University of Ljubljana, Slovenia.
#'
#' @name autohrf-datasets
#' @aliases swm
#'
#' @format
#' \describe{
#' \item{\code{swm}}{
#' fMRI dataset for a spatial working memory experiment.
#'
#' Source: Internal MBLab repository.
#'
#' 504 obs. of 3 variables
#' \itemize{
#' \item \code{roi} region of interest.
#' \item \code{time} time stamp.
#' \item \code{y} BOLD value.
#' }
#' }
#' }
#'
#' @examples
#' # load swm data
#' data_swm <- swm
#'
NULL
