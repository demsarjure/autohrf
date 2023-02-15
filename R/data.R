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
#' @aliases swm swm_autofit swm_autofit1 swm_autofit2 flanker flanker_autofit
#'
#' @format
#' \describe{
#'
#' \item{\code{swm}}{
#' fMRI dataset for a spatial working memory experiment.
#'
#' Source: Internal MBLab repository.
#'
#' 11520 obs. of 3 variables
#' \itemize{
#' \item \code{roi} region of interest.
#' \item \code{time} time stamp.
#' \item \code{y} BOLD value.
#' }
#'
#' }
#'
#' \item{\code{swm_autofit}}{
#' Stored results from a pre-completed autohrf run.
#'
#' Source: Internal MBLab repository.
#' }
#'
#' \item{\code{swm_autofit1}}{
#' Stored results from a pre-completed autohrf run.
#'
#' Source: Internal MBLab repository.
#' }
#'
#' \item{\code{swm_autofit2}}{
#' Stored results from a pre-completed autohrf run.
#'
#' Source: Internal MBLab repository.
#' }
#'
#' \item{\code{flanker}}{
#' fMRI dataset for a flanker experiment.
#'
#' Source: Internal MBLab repository.
#'
#' 192 obs. of 3 variables
#' \itemize{
#' \item \code{roi} region of interest.
#' \item \code{time} time stamp.
#' \item \code{y} BOLD value.
#' }
#'
#' }
#'
#' \item{\code{flanker_autofit}}{
#' Stored results from a pre-completed autohrf run.
#'
#' Source: Internal MBLab repository.
#' }
#'
#' }
#'
#' @examples
#' # load swm data
#' data_swm <- swm
#'
#' # load the previously completed autofits
#' autofit <- swm_autofit
#' autofit1 <- swm_autofit1
#' autofit2 <- swm_autofit2
#'
#' # load flanker data
#' data_flanker <- flanker
#'
#' # load the previously completed autofits
#' autofit3 <- flanker_autofit
NULL
