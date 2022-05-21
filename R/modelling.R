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
#' @param report Whether to print a report of the evaluation results.
#'
#' @return Returns a list that contains the model, fits of events for
#' each ROI, convolved events and the TR.
#'
#' @examples
#' # create the model
#' m <- data.frame(event = c("encoding", "delay", "response"),
#' start_time = c(0, 2.5, 12.5), duration = c(2.5, 10, 5))
#'
#' # evaluate
#' df <- swm
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
                           report = TRUE) {

  ce <- convolve_events(model, tr, f, hrf, t, p_boynton, p_spm)
  rm <- run_model(d, ce, model, roi_weights)

  em <- list(model = model, rm = rm, ce = ce, tr = tr, coefficients = rm$c)

  # report
  if (report) {
    cat("\nMean R2: ", rm$r2$mean)
    cat("\nMedian R2: ", rm$r2$median)
    cat("\nMin R2: ", rm$r2$min)
    cat("\nWeighted R2: ", rm$r2$weighted, "\n")
  }

  return(em)
}


#' @title plot_model
#' @description Plots a manually constructed model.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @import ggplot2
#' @export
#'
#' @param model_evaluation The output from the evaluate_model function.
#' @param by_roi Whether to plot the fit for each ROI independently.
#' @param ncol Number of columns in the facet wrap.
#' @param nrow Number of rows in the facet wrap.
#' @param scales Whether to free certain axes of the facet wrap.
#' @param rois A subset of ROIs to visualize.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(event      = c("encoding", "delay", "response"),
#'                      start_time = c(0,           2.65,    12.5),
#'                      duration   = c(2.65,        9.85,    3))
#'
#'
plot_model <- function(model_evaluation,
                       by_roi = TRUE,
                       ncol = NULL,
                       nrow = NULL,
                       scales = "free_y",
                       rois = NULL) {
  # init local variables for CRAN check
  event <- NULL
  roi <- NULL
  y <- NULL

  # set events
  events  <- model_evaluation$rm$events

  # prepare af
  af <- list(models = list(model_evaluation$model),
             fitness = 0,
             best = model_evaluation$ce,
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
    p <- ggplot() +
      geom_line(data = fit[fit$event == "y_m", ],
                aes(x = t, y = y), color = "black", size = 1, alpha = 0.5) +
      geom_line(data = fit[fit$event %in% events, ],
                aes(x = t, y = y, color = event, group = event)) +
      geom_line(data = fit[fit$event == "y", ],
                aes(x = t, y = y), color = "red", size = 1, alpha = 0.3) +
      ylab("") +
      xlab("time") +
      scale_fill_discrete(name = "event") +
      scale_color_discrete(name = "event") +
      facet_wrap(~ roi, scales = scales)

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
    plot_events(af)
  }
}


# a helper function for evaluating a model
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
  ce$y[1:l, ] <- ce$y[1:l, ] /
    matrix(apply(ce$y[1:l, ], 2, FUN = function(x) max(abs(x))),
            nrow = l,
            ncol = n_events,
            byrow = TRUE)

  # run through rois
  for (roi in rois) {
    # compute the linear model
    d[d$roi == roi, events] <- ce$y[1:l, ]
    m <- lm(formula(d[, c("y", events)]), d[d$roi == roi, ])

    # save component timeseries
    d[d$roi == roi, events] <- ce$y[1:l, ] *
                               matrix(m$coefficients[events],
                                      l, length(model$event), byrow = TRUE)
    d[d$roi == roi, "(Intercept)"] <- m$coefficients["(Intercept)"][[1]]
    d[d$roi == roi, "y_m"] <-
      apply(as.matrix(d[d$roi == roi, c("(Intercept)", events)]), 1, FUN = sum)
    d[d$roi == roi, "r"] <- m$residuals
    d[d$roi == roi, events] <-
      ce$y[1:l, ] *
      matrix(m$coefficients[events], l, length(model$event), byrow = TRUE) +
      m$coefficients["(Intercept)"][[1]]

    # calculate r2
    r2 <- 1 - var(m$residuals) / var(d[d$roi == roi, "y_m"])

    # calculate r2w
    r2w <- r2 * roi_weights[roi_weights$roi == roi, "weight"]

    coeffs <-
      rbind(coeffs,
            data.frame(c(r = roi,
                       as.list(m$coefficients[events]),
                       r2 = r2,
                       r2w = r2w)))
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

  # store r2
  r2 <- list(mean = mean(coeffs$r2),
             median = median(coeffs$r2),
             min = min(coeffs$r2),
             weighted = r2w)

  return(list(fit = fit, d = d, c = coeffs, r2 = r2, events = events))
}


#' @title autohrf
#' @description A function that automatically finds the parameters of model's
#' that best match the underlying data.
#' @importFrom lubridate day hour minute second seconds_to_period
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
#' @param iter Number of iterations in the genetic algorithm.
#' @param population The size of the population in the genetic algorithm.
#' @param mutation_rate The mutation rate in the genetic algorithm.
#' @param mutation_factor The mutation factor in the genetic algorithm.
#' @param elitism The degree of elitism (promote a percentage of the best
#' solutions) in the genetic algorithm.
#' @param hrf Method to use for HRF generation.
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param p_boynton Parameters for the Boynton's HRF.
#' @param p_spm Parameters for the SPM HRF.
#' @param f Upsampling factor.
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
#' df <- swm
#' autofit <- autohrf(df, model_constraints, tr = 2.5, population = 2, iter = 2)
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
                    f = 100) {

  # parameters
  pop <- population
  m_rate <- mutation_rate
  m_factor <- mutation_factor * pop
  elitism <- ceiling(pop * elitism)

  # results
  results <- list()

  # iterate over all models
  n_models <- length(model_constraints)
  total_iterations <- n_models * iter
  execution_time <- Sys.time()
  for (m in 1:n_models) {
    # get model
    current_model <- model_constraints[[m]]
    n_events <- nrow(current_model)

    # set min duration to default if not set
    if (!"min_duration" %in% colnames(current_model)) {
      current_model$min_duration <- rep(0.1, nrow(current_model))
    }

    # set max duration to default if not set
    if (!"max_duration" %in% colnames(current_model)) {
      current_model$max_duration <-
        current_model$end_time - current_model$start_time
    }

    # create first generation
    times <- create_first_generation(current_model,
                                     n_events,
                                     pop,
                                     allow_overlap)
    start_time <- times[[1]]
    end_time <- times[[2]]

    # iterate over generations
    max_fitness <- vector()
    for (i in 1:iter) {
      # calculate eta
      current_iteration <- i + ((m - 1) * iter)
      seconds <- difftime(Sys.time(), execution_time, units = "secs")
      eta <- round(seconds * (total_iterations - current_iteration + 1), 1)
      eta <- seconds_to_period(eta)
      eta <- sprintf("%02d-%02d:%02d:%02d",
                      day(eta), hour(eta), minute(eta), round(second(eta)))

      # time of execution
      execution_time <- Sys.time()

      # print
      cat("Progress:\t", current_iteration, "/",
          total_iterations, "\teta:", eta, "\n")

      # evaluate each model
      fitness <- vector()
      for (j in 1:pop) {
        # create model from data
        model <- data.frame(event = current_model$event,
                            start_time = start_time[[j]],
                            duration = end_time[[j]] - start_time[[j]])

        em <- evaluate_model(d = d,
                             model = model,
                             roi_weights = roi_weights,
                             tr = tr,
                             f = f,
                             hrf = hrf,
                             t = t,
                             p_boynton = p_boynton,
                             p_spm = p_spm,
                             report = FALSE)

        fitness <- append(fitness, em$rm$r2$weighted)
      }

      # sort
      start_time <- start_time[mixedorder(fitness, decreasing = TRUE)]
      end_time <- end_time[mixedorder(fitness, decreasing = TRUE)]
      fitness <- sort(fitness, decreasing = TRUE)
      max_fitness <- append(max_fitness, max(fitness))

      # create next gen (skip in last generation)
      if (i != iter) {
        times <- create_new_generation(elitism,
                                       pop,
                                       start_time,
                                       end_time,
                                       fitness,
                                       n_events,
                                       m_factor,
                                       m_rate,
                                       current_model,
                                       allow_overlap)

        start_time <- times[[1]]
        end_time <- times[[2]]
      }
    }

    # construct new models
    new_models <- list()
    for (j in 1:pop) {
      # create model from data
      model <- data.frame(event = current_model$event,
                          start_time = start_time[[j]],
                          duration = end_time[[j]] - start_time[[j]])

      new_models[[j]] <- model
    }

    #evaluate the best model
    ce <- convolve_events(model = new_models[[1]],
                          tr = tr,
                          f = f,
                          hrf = hrf,
                          t = t,
                          p_boynton = c(2.25, 1.25, 2),
                          p_spm = c(6, 16, 1, 1, 6, 0))

    results[[m]] <- list(models = new_models,
                         fitness = max_fitness,
                         best = ce,
                         tr = tr)
  }

  return(results)
}


# a helper function for creating the first generation
create_first_generation <- function(current_model,
                                    n_events,
                                    pop,
                                    allow_overlap) {

  start_time <- list()
  end_time <- list()
  for (i in 1:pop) {
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


# a helper function for creating a new generation of possible solutions
create_new_generation <- function(elitism,
                                  pop,
                                  start_time,
                                  end_time,
                                  fitness,
                                  n_events,
                                  m_factor,
                                  m_rate,
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
  for (j in (elitism + 1):pop) {
    # get parents
    parents <- get_parents(fitness)
    p1 <- parents[[1]]
    p2 <- parents[[2]]

    # create child
    child <- create_child(start_time,
                          end_time,
                          n_events,
                          m_rate,
                          m_factor,
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

# a helper function for getting parents
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

# a helper function for creating a child from parents
create_child <- function(start_time,
                         end_time,
                         n_events,
                         m_rate,
                         m_factor,
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
  half1 <- round(n_events / 2)
  half2 <- n_events - half1

  # take first half from p1 and second from p2
  start <- append(start1[1:half1], start2[half1 + 1:half2])
  end <- append(end1[1:half1], end2[half1 + 1:half2])

  # mutations modify values
  for (k in 1:n_events) {
  # add some random value (depends on m_factor)
    duration <- end[k] - start[k]
    mutation <- duration * m_factor

    # start
    rand <- runif(1)
    if (rand < m_rate) {
      start[k] <- start[k] + runif(1, -mutation, mutation)
    }

    # end
    rand <- runif(1)
    if (rand < m_rate) {
      end[k] <- end[k] + runif(1, -mutation, mutation)
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

    # end time needs to be lower than end_time
    # and larger than start_time + min_duration
    # and lower than start_time + max_duration
    end[k] <- min(et, max(st + min_d, min(st + max_d, end[k])))
  }

  child <- list()
  child[[1]] <- start
  child[[2]] <- end
  return(child)
}


#' @title plot_fitness
#' @description Plots how fitness changed through iterations of autohrf.
#' Use this to invesitage whether your solution converged.
#' @import ggplot2
#' @export
#'
#' @param autofit Output of the autohrf function.
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
#' df <- swm
#' autofit <- autohrf(df, model_constraints, tr = 2.5, population = 2, iter = 2)
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
#' df <- swm
#' autofit <- autohrf(df, model_constraints, tr = 2.5, population = 2, iter = 2)
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
    graphs[[i]] <- plot_events(af, i)
    i <- i + 1
  }

  if (is.null(ncol) && is.null(nrow)) {
    nrow <- i - 1
    ncol <- 1
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
#' df <- swm
#' autofit <- autohrf(df, model_constraints, tr = 2.5, population = 2, iter = 2)
#'
#' # print best models
#' get_best_models(autofit)
#'
get_best_models <- function(autofit) {
  # best models storage
  models <- list()

  # iterate over models
  i <- 1
  cat("\n----------------------------------------\n")
  for (af in autofit) {
    cat("\nModel", i, "\n\n")
    models[[i]] <- af$models[[1]]
    cat("Fitness: ", tail(af$fitness, 1), "\n\n")
    print(af$models[[1]])
    cat("\n----------------------------------------\n")

    i <- i + 1
  }

  return(models)
}
