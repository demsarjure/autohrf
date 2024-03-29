# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

#' @title evaluate_model
#' @description A function for evaluating the model against the data.
#' @import ggplot2
#' @export
#'
#' @param d A dataframe with the signal data: roi, t and y. ROI is the name of
#' the region, t is the timestamp and y the value of the signal.
#' @param model A data frame containing information about the model to use
#' and its events (event, start_time and duration).
#' @param tr MRI's repetition time.
#' @param roi_weights A data frame with ROI weights: roi, weight. ROI is the
#' name of the region, weight a number that defines the importance of that roi,
#' the default weight for a ROI is 1. If set to 2 for a particular ROI that ROI
#' will be twice as important.
#' @param hrf Method to use for HRF generation, can be "boynton" or "spm".
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param p_boynton Parameters for the Boynton's HRF.
#' @param p_spm Parameters for the SPM HRF.
#' @param f Upsampling factor.
#' @param verbose Whether to print a report of the evaluation results.
#'
#' @return Returns a list that contains the model, fits of events for
#' each ROI, convolved events, TR and evaluation scores for each ROI.
#'
#' @examples
#' # create the model
#' m <- data.frame(event = c("encoding", "delay", "response"),
#' start_time = c(0, 2.5, 12.5), duration = c(2.5, 10, 5))
#'
#' # evaluate
#' df <- flanker
#' res <- evaluate_model(df, m, tr = 2.5)
#'
evaluate_model <- function(d,
                           model,
                           tr,
                           roi_weights = NULL,
                           hrf = "spm",
                           t = 32,
                           p_boynton = c(2.25, 1.25, 2),
                           p_spm = c(6, 16, 1, 1, 6, 0),
                           f = 100,
                           verbose = TRUE) {

  ce <- autohrf::convolve_events(model = model,
                                 tr = tr,
                                 max_duration = max(d$t),
                                 hrf = hrf,
                                 t = t,
                                 p_boynton = p_boynton,
                                 p_spm = p_spm,
                                 f = f)

  rm <- run_model(d, ce, model, roi_weights)

  by_roi <- data.frame(roi = rm$c$r, r2 = rm$c$r2, r2w = rm$c$r2w,
                       bic = rm$c$bic, bicw = rm$c$bicw)

  em <- list(model = model, rm = rm, ce = ce,
             tr = tr, coefficients = rm$c, by_roi = by_roi)

  # report
  if (verbose) {
    cat("\nMean R2: ", rm$r2$mean)
    cat("\nMedian R2: ", rm$r2$median)
    cat("\nMin R2: ", rm$r2$min)
    cat("\nWeighted R2: ", rm$r2$weighted)
    cat("\n\nMean BIC: ", rm$bic$mean)
    cat("\nMedian BIC: ", rm$bic$median)
    cat("\nMin BIC: ", rm$bic$min)
    cat("\nWeighted BIC: ", rm$bic$weighted, "\n")
  }

  return(em)
}


#' @title plot_model
#' @description Plots a manually constructed model.
#' @importFrom magrittr %>%
#' @importFrom dplyr case_when mutate
#' @import ggplot2
#' @import RColorBrewer
#' @export
#'
#' @param model_evaluation The output from the evaluate_model function.
#' @param by_roi Whether to plot the fit for each ROI independently.
#' @param ncol Number of columns in the facet wrap.
#' @param nrow Number of rows in the facet wrap.
#' @param scales Whether to free certain axes of the facet wrap.
#' @param rois A subset of ROIs to visualize.
#'
#' @return A ggplot visualization of the model.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(event      = c("encoding", "delay", "response"),
#'                      start_time = c(0,           2.65,    12.5),
#'                      duration   = c(2.65,        9.85,    3))
#'
#'
plot_model <- function(model_evaluation,
                       by_roi = FALSE,
                       ncol = NULL,
                       nrow = NULL,
                       scales = "free_y",
                       rois = NULL) {
  # init local variables for CRAN check
  event <- NULL
  roi <- NULL
  y <- NULL
  signal <- NULL

  # set events
  events  <- model_evaluation$rm$events

  # prepare af
  af <- list(models = list(model_evaluation$model),
             fitness = 0,
             ce = model_evaluation$ce,
             tr = model_evaluation$tr)

  # plot
  if (by_roi) {
    # fit
    fit <- model_evaluation$rm$fit

    # filter
    if (!is.null(rois)) {
      fit <- fit %>% dplyr::filter(roi %in% rois)
      fit$roi <- factor(fit$roi, levels = rois)
    }

    # visualization
    event_data <- fit[fit$event %in% events, ]
    signal_data <- fit[!(fit$event %in% events), ]
    signal_data <- signal_data %>%
        mutate(signal = case_when(event == "y" ~ "bold", TRUE ~ "model"))
    legend_order <- c("bold", "model", events)
    my_palette <- c("#000000", "#FF7E79", brewer.pal(length(events), "Set1"))
    p <- ggplot() +
      geom_line(data = signal_data,
                aes(x = t, y = y, color = signal, group = signal), size = 1) +
      geom_line(data = event_data,
                aes(x = t, y = y, color = event, group = event)) +
      ylab("y") +
      xlab("time") +
      scale_fill_manual(values = my_palette, breaks = legend_order) +
      scale_color_manual(values = my_palette, breaks = legend_order) +
      theme(legend.title = element_blank()) +
      facet_wrap(~ roi, scales = scales)

      if (!is.null(ncol) && !is.null(nrow)) {
        cat("\nWARNING: Both ncol and nrow are provided, using only ncol!\n")
        nrow <- NULL
      }

      if (is.null(ncol) && !is.null(nrow)) {
        p <- p + facet_wrap(~ roi, nrow = nrow, scales = scales)
      } else if (!is.null(ncol) && is.null(nrow)) {
        p <- p + facet_wrap(~ roi, ncol = ncol, scales = scales)
      } else {
        p <- p + facet_wrap(~ roi, scales = scales)
      }
      # plot
      p
  } else {
    autohrf::plot_events(af)
  }
}


#' @title run_model
#' @description A helper function for evaluating a model.
#' @export
#'
#' @param d A dataframe with the signal data: roi, t and y. ROI is the name of
#' the region, t is the timestamp and y the value of the signal.
#' @param ce Result of the convolve_events function.
#' @param model A data frame containing information about the model to use
#' and its events (event, start_time and duration).
#' @param roi_weights A data frame with ROI weights: roi, weight. ROI is the
#' name of the region, weight a number that defines the importance of that roi,
#' the default weight for a ROI is 1. If set to 2 for a particular ROI that ROI
#' will be twice as important.
#'
#' @return Returns the model's evaluation.
run_model <- function(d,
                      ce,
                      model,
                      roi_weights = NULL) {

  # init local variables for CRAN check
  event <- NULL
  y <- NULL

  # set up variables
  coeffs <- c()
  fit <- c()
  rois <- unique(d$roi)
  events <- as.character(model$event)
  n_events <- length(events)

  # create weights if not defined
  default_roi_weights <- data.frame(roi = unique(d$roi), weight = 1)
  if (!is.null(roi_weights)) {
    roi_weights <- merge(default_roi_weights, roi_weights,
                         by = "roi", all.x = TRUE)
    roi_weights <-
      data.frame(roi = roi_weights$roi,
                 weight = ifelse(is.na(roi_weights$weight.y),
                                 roi_weights$weight.x, roi_weights$weight.y))
  } else {
    roi_weights <- default_roi_weights
  }

  # expand the dataframe
  d[, c("(Intercept)", events, "y_m", "r")] <- 0
  d[, "(Intercept)"] <- 1
  l <- dim(d[d$roi == rois[1], ])[1]

  # normalize
  ce[1:l, ] <- ce[1:l, ] /
    matrix(apply(as.matrix(ce[1:l, ]), 2, FUN = function(x) max(abs(x))),
            nrow = l,
            ncol = n_events,
            byrow = TRUE)

  # run through rois
  for (roi in rois) {
    # compute the linear model
    d[d$roi == roi, events] <- ce[1:l, ]
    m <- lm(formula(d[, c("y", events)]), d[d$roi == roi, ])

    # save component timeseries
    d[d$roi == roi, events] <- ce[1:l, ] *
                               matrix(m$coefficients[events],
                                      l, length(model$event), byrow = TRUE)
    d[d$roi == roi, "(Intercept)"] <- m$coefficients["(Intercept)"][[1]]
    d[d$roi == roi, "y_m"] <-
      apply(as.matrix(d[d$roi == roi, c("(Intercept)", events)]), 1, FUN = sum)
    d[d$roi == roi, "r"] <- m$residuals
    d[d$roi == roi, events] <-
      ce[1:l, ] *
      matrix(m$coefficients[events], l, length(model$event), byrow = TRUE) +
      m$coefficients["(Intercept)"][[1]]

    # calculate r2
    r2 <- 1 - var(m$residuals) / var(d[d$roi == roi, "y"])

    # calculate r2w
    r2w <- r2 * roi_weights[roi_weights$roi == roi, "weight"]

    # calculate bic
    bic <- BIC(m)

    # calculate bicw
    bicw <- bic * roi_weights[roi_weights$roi == roi, "weight"]

    coeffs <-
      rbind(coeffs,
            data.frame(c(r = roi,
                       as.list(m$coefficients[events]),
                       r2 = r2,
                       r2w = r2w,
                       bic = bic,
                       bicw = bicw)))
  }

  for (v in c("y", events, "y_m")) {
    t <- d[, c("roi", "t", v)]
    t$event <- v
    names(t) <- c("roi", "t", "y", "event")
    fit <- rbind(fit, t)
  }

  fit$event <- factor(fit$event, levels = c("y", events, "y_m"), ordered = TRUE)

  # weighted r2
  r2w <- sum(coeffs$r2w) / sum(roi_weights$weight)

  # weighted r2
  bicw <- sum(coeffs$bicw) / sum(roi_weights$weight)

  # store evaluation
  r2 <- list(mean = mean(coeffs$r2),
               median = median(coeffs$r2),
               min = min(coeffs$r2),
               weighted = r2w)

  bic <- list(mean = mean(coeffs$bic),
               median = median(coeffs$bic),
               min = min(coeffs$bic),
               weighted = bicw)

  return(list(fit = fit, d = d, c = coeffs, r2 = r2, bic = bic, events = events))
}


#' @title autohrf
#' @description A function that automatically finds the parameters of model's
#' that best match the underlying data.
#' @importFrom lubridate day hour minute second seconds_to_period
#' @import doParallel
#' @import foreach
#' @import gtools
#' @import stats
#' @export
#'
#' @param d A dataframe with the signal data: roi, t and y. ROI is the name of
#' the region, t is the timestamp and y the value of the signal.
#' @param model_constraints A list of model specifications to use for fitting.
#' Each specification is represented as a data frame containing information
#' about it (event, start_time, end_time, min_duration and max_duration).
#' @param tr MRI's repetition time.
#' @param roi_weights A data frame with ROI weights: roi, weight. ROI is the
#' name of the region, weight a number that defines the importance of that roi,
#' the default weight for a ROI is 1. If set to 2 for a particular ROI that ROI
#' will be twice as important.
#' @param allow_overlap Whether to allow overlap between events.
#' @param population The size of the population in the genetic algorithm.
#' @param iter Number of iterations in the genetic algorithm.
#' @param mutation_rate The mutation rate in the genetic algorithm.
#' @param mutation_factor The mutation factor in the genetic algorithm.
#' @param elitism The degree of elitism (promote a percentage of the best
#' solutions) in the genetic algorithm.
#' @param hrf Method to use for HRF generation.
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param p_boynton Parameters for the Boynton's HRF.
#' @param p_spm Parameters for the SPM HRF.
#' @param f Upsampling factor.
#' @param cores Number of cores to use for parallel processing. Set to the
#' number of provided model constraints by default.
#' @param autohrf Results of a previous autohrf run to continue.
#' @param verbose Whether to print progress of the fitting process.
#'
#' @return A list containing model fits for each of the provided model
#' specifications.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(
#'   event        = c("encoding", "delay", "response"),
#'   start_time   = c(0,          2.65,     12.5),
#'   end_time     = c(3,          12.5,     16)
#' )
#'
#' model4 <- data.frame(
#'   event        = c("fixation", "target", "delay", "response"),
#'   start_time   = c(0,          2.5,      2.65,    12.5),
#'   end_time     = c(2.5,        3,        12.5,    15.5)
#' )
#'
#' model_constraints <- list(model3, model4)
#'
#' # run autohrf
#' df <- flanker
#' autofit <- autohrf(df, model_constraints, tr = 2.5,
#'                    population = 2, iter = 2, cores = 1)
#'
autohrf <- function(d,
                    model_constraints,
                    tr,
                    roi_weights = NULL,
                    allow_overlap = FALSE,
                    population = 100,
                    iter = 100,
                    mutation_rate = 0.1,
                    mutation_factor = 0.05,
                    elitism = 0.1,
                    hrf = "spm",
                    t = 32,
                    p_boynton = c(2.25, 1.25, 2),
                    p_spm = c(6, 16, 1, 1, 6, 0),
                    f = 100,
                    cores = NULL,
                    autohrf = NULL,
                    verbose = TRUE) {

  # parameters
  n_models <- length(model_constraints)
  elitism <- ceiling(population * elitism)

  # set cores to number of model constraints
  if (is.null(cores)) {
    cores <- n_models

    # decrease cores if too large
    if (cores > parallel::detectCores() - 1) {
      cores <- parallel::detectCores() - 1
    }
  }

  # setup parallelism
  cl <- NULL
  if (cores != 1) {
    cl <- parallel::makeCluster(cores, outfile = "")
    doParallel::registerDoParallel(cl)
  }

  i <- 1
  results <- foreach(i = 1:n_models) %dopar% {
    autohrf::fit_to_constraints(i,
                                d,
                                model_constraints,
                                tr,
                                roi_weights,
                                allow_overlap,
                                population,
                                iter,
                                mutation_rate,
                                mutation_factor,
                                elitism,
                                hrf,
                                t,
                                p_boynton,
                                p_spm,
                                f,
                                autohrf,
                                verbose)
  }

  # close cluster
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
  }

  return(results)
}


#' @title fit_to_constraints
#' @description A helper function for fitting a model to constraints.
#' @export
#'
#' @param model_id ID of the model.
#' @param d A dataframe with the signal data: roi, t and y. ROI is the name of
#' the region, t is the timestamp and y the value of the signal.
#' @param model_constraints A list of model specifications to use for fitting.
#' Each specification is represented as a data frame containing information
#' about it (event, start_time, end_time, min_duration and max_duration).
#' @param tr MRI's repetition time.
#' @param roi_weights A data frame with ROI weights: roi, weight. ROI is the
#' name of the region, weight a number that defines the importance of that roi,
#' the default weight for a ROI is 1. If set to 2 for a particular ROI that ROI
#' will be twice as important.
#' @param allow_overlap Whether to allow overlap between events.
#' @param population The size of the population in the genetic algorithm.
#' @param iter Number of iterations in the genetic algorithm.
#' @param mutation_rate The mutation rate in the genetic algorithm.
#' @param mutation_factor The mutation factor in the genetic algorithm.
#' @param elitism The degree of elitism (promote a percentage of the best
#' solutions) in the genetic algorithm.
#' @param hrf Method to use for HRF generation.
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param p_boynton Parameters for the Boynton's HRF.
#' @param p_spm Parameters for the SPM HRF.
#' @param f Upsampling factor.
#' @param autohrf Results of a previous autohrf run to continue.
#' @param verbose Whether to print progress of the fitting process.
#'
#' @return Returns the best model given provided constraints.
fit_to_constraints <- function(model_id,
                               d,
                               model_constraints,
                               tr,
                               roi_weights,
                               allow_overlap,
                               population,
                               iter,
                               mutation_rate,
                               mutation_factor,
                               elitism,
                               hrf,
                               t,
                               p_boynton,
                               p_spm,
                               f,
                               autohrf = NULL,
                               verbose = TRUE) {

  # iterate over all models
  execution_time <- Sys.time()

  # get model
  current_model <- model_constraints[[model_id]]

  # get number of events
  n_events <- nrow(current_model)

  # some basic checks
  if (!"event" %in% colnames(current_model)) {
    cat("\nERROR: Missing event column in one of provided model constraints!")
  }
  if (!"start_time" %in% colnames(current_model)) {
    cat("\nERROR: Missing start_time column in one of provided model
         constraints!")
  }
  if (!"end_time" %in% colnames(current_model)) {
    cat("\nERROR: Missing end_time column in one of provided model
        constraints!")
  }
  for (j in 1:n_events) {
    # get event
    event <- current_model[j, ]

    # check if end_time larger than start_time
    if (event$end_time < event$start_time) {
        cat("\nERROR: end_time is smaller than start_time in one of provided
            model constraints!")
    }
  }

  # set min duration to default if not set
  if (!"min_duration" %in% colnames(current_model)) {
    current_model$min_duration <- rep(0.01, nrow(current_model))
  }

  # set max duration to default if not set
  if (!"max_duration" %in% colnames(current_model)) {
    current_model$max_duration <-
      current_model$end_time - current_model$start_time
  }

  # create first generation if not continuing
  if (is.null(autohrf)) {
    times <- autohrf::create_first_generation(current_model,
                                              n_events,
                                              population,
                                              allow_overlap)
    start_time <- times[[1]]
    end_time <- times[[2]]
    start_iter <- 1
  } else {
    ah <- autohrf[[model_id]]
    start_time <- ah$start_time
    end_time <- ah$end_time
    iter <- iter + ah$iter
    start_iter <- ah$iter + 1

    # check if populations sizes are matcing
    if (population != length(start_time)) {
      cat(paste0("ERROR: population in the previous run, [", length(start_time),
                "] does not match the population in the new run [",
                population, "]!"))
    }
  }

  # iterate over generations
  max_fitness <- vector()
  for (i in start_iter:iter) {
    # calculate eta
    seconds <- difftime(Sys.time(), execution_time, units = "secs")
    eta <- round(seconds * (iter - i + 1), 1)
    eta <- seconds_to_period(eta)
    eta <- sprintf("%02d-%02d:%02d:%02d",
                    day(eta), hour(eta), minute(eta), round(second(eta)))

    # time of execution
    execution_time <- Sys.time()

    # print
    if (verbose) {
      cat("Model:", model_id, "\t|\t",
          "Iteration:", i, "/", iter, "\t|\t",
          "ETA:", eta, "\n")
    }

    # evaluate each model
    fitness <- vector()
    for (j in 1:population) {
      # round
      start_time[[j]] <- round(start_time[[j]], 2)
      end_time[[j]] <- round(end_time[[j]], 2)

      # create model from data
      model <- data.frame(event = current_model$event,
                          start_time = start_time[[j]],
                          duration = end_time[[j]] - start_time[[j]])

      em <- autohrf::evaluate_model(d = d,
                                    model = model,
                                    roi_weights = roi_weights,
                                    tr = tr,
                                    f = f,
                                    hrf = hrf,
                                    t = t,
                                    p_boynton = p_boynton,
                                    p_spm = p_spm,
                                    verbose = FALSE)

      fitness <- append(fitness, em$rm$r2$weighted)
    }

    # sort
    start_time <- start_time[mixedorder(fitness, decreasing = TRUE)]
    end_time <- end_time[mixedorder(fitness, decreasing = TRUE)]
    fitness <- sort(fitness, decreasing = TRUE)
    max_fitness <- append(max_fitness, max(fitness))

    # create next gen (skip in last generation)
    if (i != iter) {
      times <- autohrf::create_new_generation(elitism,
                                              population,
                                              start_time,
                                              end_time,
                                              fitness,
                                              n_events,
                                              mutation_factor,
                                              mutation_rate,
                                              current_model,
                                              allow_overlap)

      start_time <- times[[1]]
      end_time <- times[[2]]
    }
  }

  # construct new models
  new_models <- list()
  for (j in 1:population) {
    # round
    start_time[[j]] <- round(start_time[[j]], 2)
    end_time[[j]] <- round(end_time[[j]], 2)

    # create model from data
    model <- data.frame(event = current_model$event,
                        start_time = start_time[[j]],
                        duration = end_time[[j]] - start_time[[j]])

    new_models[[j]] <- model
  }

  # store convolved events
  ce <- autohrf::convolve_events(model = new_models[[1]],
                                 tr = tr,
                                 max_duration = max(d$t),
                                 hrf = hrf,
                                 t = t,
                                 p_boynton = p_boynton,
                                 p_spm = p_spm,
                                 f = f)

  result <- list(models = new_models,
                 fitness = max_fitness,
                 ce = ce,
                 tr = tr,
                 start_time = start_time,
                 end_time = end_time,
                 iter = iter)

  return(result)
}


#' @title create_first_generation
#' @description A helper function for creating the first generation.
#' @export
#'
#' @param current_model The constraints of the current model.
#' @param n_events Number of events in the model.
#' @param population The size of the population in the genetic algorithm.
#' @param allow_overlap Whether to allow overlap between events.
#'
#' @return Returns the first generation of models.
create_first_generation <- function(current_model,
                                    n_events,
                                    population,
                                    allow_overlap) {

  start_time <- list()
  end_time <- list()
  for (i in 1:population) {
    # variables for start and end times
    starts <- NULL
    ends <- NULL

    for (j in 1:n_events) {
      # get event
      event <- current_model[j, ]

      # create random and add to variable
      duration <- runif(1, event$min_duration, event$max_duration)
      start <- runif(1, event$start_time, event$end_time - duration)
      end <- start + duration

      # round to 2 digits
      start <- round(start, 2)
      end <- round(end, 2)

      starts <- append(starts, start)
      ends <- append(ends, end)
    }

    # sort in case of overlaps
    if (!allow_overlap) {
      starts <- sort(starts)
      ends <- sort(ends)
    }

    # store
    start_time[[i]] <- starts
    end_time[[i]] <- ends
  }

  times <- list()
  times[[1]] <- start_time
  times[[2]] <- end_time
  return(times)
}


#' @title create_new_generation
#' @description A helper function for creating a new generation of possible
#' solutions.
#' @importFrom lubridate day hour minute second seconds_to_period
#' @import doParallel
#' @import foreach
#' @import gtools
#' @import stats
#' @export
#'
#' @param elitism The degree of elitism (promote a percentage of the best
#' solutions) in the genetic algorithm.
#' @param population The size of the population in the genetic algorithm.
#' @param start_time A list with model's event start times.
#' @param end_time A list with model's event end times.
#' @param fitness A fitness score of all candidate models.
#' @param n_events Number of events in the model.
#' @param mutation_factor The mutation factor in the genetic algorithm.
#' @param mutation_rate The mutation rate in the genetic algorithm.
#' @param current_model The constraints of the current model.
#' @param allow_overlap Whether to allow overlap between events.
#'
#' @return A new generation of candidate models.
create_new_generation <- function(elitism,
                                  population,
                                  start_time,
                                  end_time,
                                  fitness,
                                  n_events,
                                  mutation_factor,
                                  mutation_rate,
                                  current_model,
                                  allow_overlap) {

  new_start <- list()
  new_end  <- list()

  # copy the best ones (elitism)
  for (e in 1:elitism) {
    new_start[[e]] <- start_time[[e]]
    new_end[[e]] <- end_time[[e]]
  }

  # create new ones with genetic algorithms
  for (j in (elitism + 1):population) {
    # get parents
    parents <- autohrf::get_parents(fitness)
    p1 <- parents[[1]]
    p2 <- parents[[2]]

    # create child
    child <- autohrf::create_child(start_time,
                                   end_time,
                                   n_events,
                                   mutation_rate,
                                   mutation_factor,
                                   current_model,
                                   p1,
                                   p2,
                                   allow_overlap)

    # store
    new_start[[j]] <- child[[1]]
    new_end[[j]] <- child[[2]]
  }

  # replace old generation with new
  times <- list()
  times[[1]] <- new_start
  times[[2]] <- new_end

  return(times)
}


#' @title get_parents
#' @description A helper function for getting parents for the child model.
#' @export
#'
#' @param fitness A fitness score of all candidate models.
#'
#' @return Parents for the child model.
get_parents <- function(fitness) {
  sum_fitness <- sum(fitness)

  # get parents
  p1 <- 1
  p2 <- 1

  for (k in 1:2) {
    # random number for parent lottery
    rand <- runif(1, 0, sum_fitness)
    sum <- 0
    p_temp <- 1
    # iterate until sum larger than random
    while (TRUE) {
      sum <- sum + fitness[p_temp]
      if (sum >= rand) {
        break
      }
      p_temp <- p_temp + 1
    }
    # set parents
    if (k == 1) {
      p1 <- p_temp
    } else {
      p2 <- p_temp
    }
  }

  # store parents
  parents <- list()
  parents[[1]] <- p1
  parents[[2]] <- p2

  return(parents)
}


#' @title create_child
#' @description A helper function for creating a child from parents.
#' @export
#'
#' @param start_time A list with model's event start times.
#' @param end_time A list with model's event end times.
#' @param n_events Number of events in the model.
#' @param mutation_rate The mutation rate in the genetic algorithm.
#' @param mutation_factor The mutation factor in the genetic algorithm.
#' @param current_model The constraints of the current model.
#' @param p1 The first selected parent.
#' @param p2 The second selected parent.
#' @param allow_overlap Whether to allow overlap between events.
#'
#' @return A child model created from two parents.
create_child <- function(start_time,
                         end_time,
                         n_events,
                         mutation_rate,
                         mutation_factor,
                         current_model,
                         p1,
                         p2,
                         allow_overlap) {
  # crossover
  start1 <- start_time[[p1]]
  end1 <- end_time[[p1]]
  start2 <- start_time[[p2]]
  end2 <- end_time[[p2]]

  # take half of p1
  half1 <- ceiling(n_events / 2)
  half2 <- n_events - half1

  # only 1 event take start from one, end from the other
  if (half2 == 0) {
    start <- start1
    end <- end2
  } else {
    # take first half from p1 and second from p2
    start <- append(start1[1:half1], start2[half1 + 1:half2])
    end <- append(end1[1:half1], end2[half1 + 1:half2])
  }

  # mutations modify values
  for (k in 1:n_events) {
    # add some random value (depends on mutation_factor)
    duration <- end[k] - start[k]
    mutation <- duration * mutation_factor

    # start
    rand <- runif(1)
    if (rand < mutation_rate) {
      change <- max(0.01, runif(1, -mutation, mutation))
      start[k] <- start[k] + change
    }

    # end
    rand <- runif(1)
    if (rand < mutation_rate) {
      change <- max(0.01, runif(1, -mutation, mutation))
      end[k] <- end[k] + change
    }

    if (end[k] < start[k]) {
      temp <- start[k]
      start[k] <- end[k]
      end[k] <- temp
    }
  }

  # clamp model events to boundaries
  for (k in 1:n_events) {
    # get model
    event <- current_model[k, ]

    # extract info
    st <- event$start_time
    et <- event$end_time
    min_d <- event$min_duration
    max_d <- event$max_duration

    # start time needs to be larger than start_time
    # and lower than end_time - min_duration
    start[k] <- max(st, min(et - min_d, start[k]))
    start[k] <- round(start[k], 2)

    # end time needs to be lower than end_time
    # and larger than start_time + min_duration
    # and lower than start_time + max_duration
    end[k] <- min(et, max(st + min_d, min(st + max_d, end[k])))
    end[k] <- round(end[k], 2)
  }

  child <- list()
  child[[1]] <- start
  child[[2]] <- end
  return(child)
}


#' @title plot_fitness
#' @description Plots how fitness changed through iterations of autohrf.
#' Use this to investigate whether your solution converged.
#' @import ggplot2
#' @export
#'
#' @param autofit Output of the autohrf function.
#'
#' @return A ggplot visualization of fitness through time.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(
#'   event        = c("encoding", "delay", "response"),
#'   start_time   = c(0,          2.65,     12.5),
#'   end_time     = c(3,          12.5,     16)
#' )
#'
#' model4 <- data.frame(
#'   event        = c("fixation", "target", "delay", "response"),
#'   start_time   = c(0,          2.5,      2.65,    12.5),
#'   end_time     = c(2.5,        3,        12.5,    15.5)
#' )
#'
#' model_constraints <- list(model3, model4)
#'
#' # run autohrf
#' df <- flanker
#' autofit <- autohrf(df, model_constraints, tr = 2.5,
#'                    population = 2, iter = 2, cores = 1)
#'
#' # plot fitness
#' plot_fitness(autofit)
#'
plot_fitness <- function(autofit) {
  # init local variables for CRAN check
  index <- NULL
  fitness <- NULL
  model <- NULL

  # empty variables for
  fitness <- NULL

  # iterate over all fits and prepare the data frame
  for (i in seq_len(length(autofit))) {
    fit <- data.frame(fitness = autofit[[i]]$fitness,
                      index = seq(length(autofit[[i]]$fitness)),
                      model = as.factor(i))

    fitness <- rbind(fitness, fit)
  }

  # plot the results
  ggplot(data = fitness, aes(x = index, y = fitness, color = model)) +
    geom_line(size = 1) +
    ylab("Fitness") +
    xlab("Iteration") +
    labs(color = "Model") +
    scale_color_brewer(type = "qual", palette = "Set1")
}

#' @title plot_best_models
#' @description Plots the best fitted model for each of the specs in autohrf.
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @param autofit Output of the autohrf function.
#' @param ncol Number of columns in the plot.
#' @param nrow Number of rows in the plot.
#'
#' @return Plots the grid containing a visualization of the best models for each
#' of the provided constraints.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(
#'   event        = c("encoding", "delay", "response"),
#'   start_time   = c(0,          2.65,     12.5),
#'   end_time     = c(3,          12.5,     16)
#' )
#'
#' model4 <- data.frame(
#'   event        = c("fixation", "target", "delay", "response"),
#'   start_time   = c(0,          2.5,      2.65,    12.5),
#'   end_time     = c(2.5,        3,        12.5,    15.5)
#' )
#'
#' model_constraints <- list(model3, model4)
#'
#' # run autohrf
#' df <- flanker
#' autofit <- autohrf(df, model_constraints, tr = 2.5,
#'                    population = 2, iter = 2, cores = 1)
#'
#' # plot best models
#' plot_best_models(autofit)
#'
plot_best_models <- function(autofit, ncol = NULL, nrow = NULL) {
  # plot list storage
  graphs <- list()
  i <- 1

  # iterate over models
  for (af in autofit) {
    graphs[[i]] <- autohrf::plot_events(af, i)
    i <- i + 1
  }

  if (is.null(ncol) && is.null(nrow)) {
    nrow <- i - 1
    ncol <- 1
  } else if (!is.null(ncol) && !is.null(nrow)) {
    cat("\nWARNING: Both ncol and nrow are provided, using only ncol!\n")
    nrow <- NULL
  }

  # plot grid
  cowplot::plot_grid(plotlist = graphs, nrow = nrow, ncol = ncol, scale = 0.95)
}


#' @title get_best_models
#' @description Returns and prints the best fitted model for each of the specs
#' used in autohrf.
#' @import utils
#' @export
#'
#' @param autofit Output of the autohrf function.
#' @param return_fitness Whether to return models or fitness.
#' @param verbose Whether to print information or only return the result.
#'
#' @return Returns a list containing the best models for each of the provided
#' constraints.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(
#'   event        = c("encoding", "delay", "response"),
#'   start_time   = c(0,          2.65,     12.5),
#'   end_time     = c(3,          12.5,     16)
#' )
#'
#' model4 <- data.frame(
#'   event        = c("fixation", "target", "delay", "response"),
#'   start_time   = c(0,          2.5,      2.65,    12.5),
#'   end_time     = c(2.5,        3,        12.5,    15.5)
#' )
#'
#' model_constraints <- list(model3, model4)
#'
#' # run autohrf
#' df <- flanker
#' autofit <- autohrf(df, model_constraints, tr = 2.5,
#'                    population = 2, iter = 2, cores = 1)
#'
#' # print best models
#' get_best_models(autofit)
#'
get_best_models <- function(autofit, return_fitness = FALSE, verbose = TRUE) {
  # best models storage
  models <- list()
  fitness <- vector()

  # iterate over models
  i <- 1
  for (af in autofit) {
    models[[i]] <- af$models[[1]]
    fitness <- c(fitness, tail(af$fitness, 1))

    if (verbose) {
      cat("\n----------------------------------------\n")
      cat("\nModel", i, "\n\n")
      cat("Fitness: ", tail(af$fitness, 1), "\n\n")
      print(af$models[[1]])
      cat("\n----------------------------------------\n")
    }

    i <- i + 1
  }

  if (return_fitness) {
    return(fitness)
  }

  return(models)
}
