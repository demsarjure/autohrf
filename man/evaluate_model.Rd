% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelling.R
\name{evaluate_model}
\alias{evaluate_model}
\title{evaluate_model}
\usage{
evaluate_model(
  d,
  model,
  tr,
  roi_weights = NULL,
  hrf = "spm",
  t = 32,
  p_boynton = c(2.25, 1.25, 2),
  p_spm = c(6, 16, 1, 1, 6, 0),
  f = 100,
  verbose = TRUE
)
}
\arguments{
\item{d}{A dataframe with the signal data: roi, t and y. ROI is the name of
the region, t is the timestamp and y the value of the signal.}

\item{model}{A data frame containing information about the model to use
and its events (event, start_time and duration).}

\item{tr}{MRI's repetition time.}

\item{roi_weights}{A data frame with ROI weights: roi, weight. ROI is the
name of the region, weight a number that defines the importance of that roi,
the default weight for a ROI is 1. If set to 2 for a particular ROI that ROI
will be twice as important.}

\item{hrf}{Method to use for HRF generation, can be "boynton" or "spm".}

\item{t}{The t parameter for Boynton or SPM HRF generation.}

\item{p_boynton}{Parameters for the Boynton's HRF.}

\item{p_spm}{Parameters for the SPM HRF.}

\item{f}{Upsampling factor.}

\item{verbose}{Whether to print a report of the evaluation results.}
}
\value{
Returns a list that contains the model, fits of events for
each ROI, convolved events, TR and evaluation scores for each ROI.
}
\description{
A function for evaluating the model against the data.
}
\examples{
# create the model
m <- data.frame(event = c("encoding", "delay", "response"),
start_time = c(0, 2.5, 12.5), duration = c(2.5, 10, 5))

# evaluate
df <- flanker
res <- evaluate_model(df, m, tr = 2.5)

}
