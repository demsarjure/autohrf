# SPDX-FileCopyrightText: 2022 Jure Demšar, Nina Purg, Grega Repovš
#
# SPDX-License-Identifier: GPL-3.0-or-later

#' @title run_model
#' @description A function for evaluating the model against the data.
#' @import ggplot2
#' @export
#'
#' @param d A dataframe with the signal data: roi, t and x. ROI is the name of the region, t timestamps and x values of the signal.
#' @param r The output from the convolve_events function.
#' @param model A data frame containing information about the model to use and its events (event, start_time and duration).
#' @param normalize Whether to normalize the signal.
#' @param report Wheter to plot the report of once done.
#'
#' @return Returns a list that contains the HRF function, fits of events for each ROI, estimates of fit quality for each ROI and a summary of model's fits.
#'
#' @examples
#' # create the model
#' m <- data.frame(event = c("encoding", "delay", "response"),
#' start_time = c(0, 2.5, 12.5), duration = c(2.5, 10, 5))
#'
#' # convolve
#' r <- convolve_events(m)
#'
#' # evaluate
#' df <- swm
#' res <- run_model(df, r, m)
#'
run_model <- function(d, r, model, normalize=TRUE, report=TRUE) {
  # init local variables for CRAN check
  event <- NULL
  y <- NULL

  # set up variables
  coeffs <- c()
  fit <- c()
  rois <- unique(d$roi)
  events <- as.character(model$event)
  n_events <- length(events)

  # expand the dataframe
  d[, c("(Intercept)", events, "y", "r")] <- 0
  d[, "(Intercept)"] <- 1
  l = dim(d[d$roi==rois[1],])[1]

  # run through rois
  for (roi in rois) {
    # normalize if needed
    if (normalize) {
      r$x[1:l,] <- r$x[1:l,] /
        matrix(apply(r$x[1:l,], 2, FUN=function(x) max(abs(x))), nrow=l, ncol=n_events, byrow=TRUE)
    }

    # compute the linear model
    d[d$roi == roi, events] <- r$x[1:l,]
    m <- lm(formula(d[, c("x", events)]), d[d$roi == roi,] )

    # save component timeseries
    d[d$roi == roi, events] <- r$x[1:l,] * matrix(m$coefficients[events], l, length(model$event), byrow=TRUE)
    d[d$roi == roi, "(Intercept)"] <- m$coefficients["(Intercept)"][[1]]
    d[d$roi == roi, "y"] <- apply(as.matrix(d[d$roi == roi, c("(Intercept)", events)]), 1, FUN=sum)
    d[d$roi == roi, "r"] <- m$residuals
    d[d$roi == roi, events] <- r$x[1:l,] * matrix(m$coefficients[events], l, length(model$event), byrow=TRUE) + m$coefficients["(Intercept)"][[1]]

    r2 <- 1 - var(m$residuals)/var(d[d$roi == roi, "x"])

    coeffs <- rbind(coeffs, data.frame(c(r=roi, as.list(m$coefficients[events]), r2=r2)))
  }

  for (v in c("x", events, "y")) {
    t <- d[, c("roi", "t", v)]
    t$event <- v
    names(t) <- c("roi", "t", "y", "event")
    fit <- rbind(fit, t)
  }

  fit$event <- factor(fit$event, levels=c("x", events, "y"), ordered=TRUE)

  r = list(mean=mean(coeffs$r2), median=median(coeffs$r2), min=min(coeffs$r2))

  # print the report
  if (report) {
    print(r)

    p <- ggplot() +
      geom_line(data=fit[fit$event == "x",],      aes(x=t, y=y), color="black", size=1, alpha=0.5) +
      geom_line(data=fit[fit$event %in% events,], aes(x=t, y=y, color=event, group=event)) +
      geom_line(data=fit[fit$event == "y",],      aes(x=t, y=y), color="red", size=1, alpha=0.3) +
      ylab("") + xlab("time") + scale_fill_discrete(name="event") + scale_color_discrete(name="event") +
      facet_wrap(~ roi, scales="free_y")
    print(p)

    print(ggplot(fit, aes(t, y, color=event)) + geom_line() + facet_wrap(~ roi, scales="free_y"))
  }

  return(list(fit=fit, d=d, c=coeffs, r=r))
}


#' @title autohrf
#' @description A function that automatically finds the parameters of model's that best match the underlying data.
#' @import gtools lubridate stats
#' @export
#'
#' @param d A dataframe with the signal data: roi, t and x. ROI is the name of the region, t timestamps and x values of the signal.
#' @param model_specs A list of model specifications to use for fitting. Each specification is represented as a data frame containing information about it (event, start_time, end_time, min_duration and max_duration).
#' @param population The size of the population for use in the genetic algorithm.
#' @param iter Number of iterations in the genetic algorithm.
#' @param mutation_rate The mutation rate in the genetic algorithm.
#' @param mutation_factor The mutation factor in the genetic algorithm.
#' @param elitism Whether to use elitism (promote a number of the best solutions) in the genetic algorithm.
#' @param tr MRI's repetition time.
#' @param f Downsampling frequency.
#' @param method Can be "middle" or "mean". Middle will return integer results, mean will return floats.
#' @param hrf Method to use for HRF generation.
#' @param t The t parameter for Boynton or SPM HRF generation.
#' @param delta The delta parameter of Boynton's HRF.
#' @param tau The tau parameter of Boynton's HRF.
#' @param alpha The alpha parameter of Boynton's HRF.
#' @param p The p parameter of SPM's HRF.
#'
#' @return A list containing model fits for each of the provided model specifications.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(event = c("encoding", "delay", "response"),
#'                      start_time = c(0, 2.65, 12.5),
#'                      end_time = c(3, 12.5, 16),
#'                      min_duration = c(1, 5, 1))
#'
#' model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
#'                      start_time = c(0, 2.5, 2.65, 12.5),
#'                      end_time = c(2.5, 3, 12.5, 15.5),
#'                      min_duration = c(1, 0.1, 5, 1))
#'
#' model_specs <- list(model3, model4)
#'
#' # run autohrf
#' df <- swm
#' autofit <- autohrf(df, model_specs, population=2, iter=2)
#'
autohrf <- function(d,
                    model_specs,
                    population = 100,
                    iter = 100,
                    mutation_rate = 0.1,
                    mutation_factor = 0.05,
                    elitism = 0.1,
                    tr=2.5,
                    f=100,
                    method = "middle",
                    hrf = "boynton",
                    t=32,
                    delta=2.25, tau=1.25, alpha=2,
                    p=c(6, 16, 1, 1, 6, 0, 32)) {

  # parameters
  pop <- population
  m_rate <- mutation_rate
  m_factor <- mutation_factor * pop
  elitism <- ceiling(pop * elitism)

  # results
  results <- list()

  # iterate over all models
  n_models <- length(model_specs)
  total_iterations <- n_models * iter
  execution_time <- NULL
  for (m in 1:n_models) {
    # get model
    current_model <- model_specs[[m]]

    # generate starting
    start_time <- list()
    end_time <- list()
    n_events <- nrow(current_model)
    for (i in 1:pop) {
      # variables for start and end times
      starts <- NULL
      ends <- NULL
      for (j in 1:n_events) {
        # get event
        event <- current_model[j, ]

        # create random and add to variable
        start <- runif(1, event$start_time, event$end_time - event$min_duration)
        end <- runif(1, start + event$min_duration, event$end_time)
        starts <- append(starts, start)
        ends <- append(ends, end)
      }

      # sort in case of overlaps
      starts <- sort(starts)
      ends <- sort(ends)

      # store
      start_time[[i]] <- starts
      end_time[[i]] <- ends
    }

    # iterate over generations
    max_fitness <- vector()
    for (i in 1:iter) {
      # calculate eta
      current_iteration <- i + ((m-1) * iter)
      if (!is.null(execution_time)) {
        seconds <- difftime(Sys.time(), execution_time, units="secs")
        eta <- round(seconds * (total_iterations - current_iteration + 1), 1)
        eta <- seconds_to_period(eta)
        eta <- sprintf('%02d-%02d:%02d:%02d',
                       day(eta), hour(eta), minute(eta), round(second(eta)))
      } else {
        eta <- "x"
      }

      # time of execution
      execution_time <- Sys.time()

      # print
      cat("Progress:\t", current_iteration, "/", total_iterations, "\teta:", eta, "\n")

      # evaluate each model
      fitness <- vector()
      for (j in 1:pop) {
        # create model from data
        model <- data.frame(event = current_model$event,
                            start_time = start_time[[j]],
                            duration = end_time[[j]] - start_time[[j]])

        r <- convolve_events(model=model,
                             tr=tr,
                             f=f,
                             method=method,
                             hrf=hrf,
                             t=t,
                             delta=delta, tau=tau, alpha=alpha,
                             p=p)

        rm <- run_model(d=d, r=r, model=model, report=FALSE)
        r2 <- rm$r$mean
        fitness <- append(fitness, r2)
      }

      # sort
      start_time <- start_time[mixedorder(fitness, decreasing=TRUE)]
      end_time <- end_time[mixedorder(fitness, decreasing=TRUE)]
      fitness <- sort(fitness, decreasing=TRUE)
      max_fitness <- append(max_fitness, max(fitness))
      sum_fitness <- sum(fitness)

      # create next gen (skip in last generation)
      if (i != iter) {
        new_start <- list()
        new_end  <- list()

        # copy the best ones (elitism)
        for (e in 1:elitism) {
          new_start[[e]] <- start_time[[e]]
          new_end[[e]] <- end_time[[e]]
        }

        # create new ones with genetic algorithms
        for (j in (elitism+1):pop) {
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

          # crossover
          start1 <- start_time[[p1]]
          end1 <- end_time[[p1]]
          start2 <- start_time[[p2]]
          end2 <- end_time[[p2]]

          # take half of p1
          half1 <- round(n_events/2)
          half2 <- n_events - half1

          # take first half from p1 and second from p2
          start <- append(start1[1:half1], start2[half1+1:half2])
          end <- append(end1[1:half1], end2[half1+1:half2])

          # mutations modify values
          for (k in 1:n_events) {
            # time
            rand <- runif(1)
            if (rand < m_rate) {
              # add some random value (depends on m_factor)
              start[k] <- start[k] + runif(1, -m_factor, m_factor)
              end[k] <- end[k] + runif(1, -m_factor, m_factor)

              if (end[k] < start[k]) {
                temp <- start[k]
                start[k] <- end[k]
                end[k] <- temp
              }
            }
          }

          # clamp events to boundaries
          for (k in 1:n_events) {
            # get event
            event <- current_model[k, ]

            end_limit <- event$end_time
            # get boundaries
            if (k == 1) {
              start_limit <- event$start_time
            } else {
              start_limit <- max(event$start_time, end[k-1])
            }

            # clamp to boundaries
            start[k] <- max(start_limit, min(end_limit-event$min_duration, start[k]))
            end[k] <- max(start[k]+event$min_duration, min(end_limit, end[k]))

            # prevent too short events
            if ((end[k] - start[k]) < event$min_duration) {
              end[k] = start[k] + event$min_duration
            }
          }

          # store
          new_start[[j]] <- start
          new_end[[j]] <- end
        }

        # replace old generation with new and repeat
        start_time <- new_start
        end_time <- new_end
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
    r <- convolve_events(model=new_models[[1]],
                         tr=tr,
                         f=f,
                         method=method,
                         hrf=hrf,
                         t=t,
                         delta=delta, tau=tau, alpha=alpha,
                         p=p)

    results[[m]] <- list(models=new_models,
                         fitness=max_fitness,
                         best=r,
                         tr=tr)
  }

  return(results)
}


#' @title plot_fitness
#' @description Plots how fitness changed through iterations of autohrf. Use this to invesitage whether your solution converged.
#' @import ggplot2
#' @export
#'
#' @param autofit Output of the autohrf function.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(event = c("encoding", "delay", "response"),
#'                      start_time = c(0, 2.65, 12.5),
#'                      end_time = c(3, 12.5, 16),
#'                      min_duration = c(1, 5, 1))
#'
#' model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
#'                      start_time = c(0, 2.5, 2.65, 12.5),
#'                      end_time = c(2.5, 3, 12.5, 15.5),
#'                      min_duration = c(1, 0.1, 5, 1))
#'
#' model_specs <- list(model3, model4)
#'
#' # run autohrf
#' df <- swm
#' autofit <- autohrf(df, model_specs, population=2, iter=2)
#'
#' # plot fitness
#' plot_fitness(autofit)
#'
plot_fitness <- function(autofit) {
  # init local variables for CRAN check
  Index <- NULL
  Fitness <- NULL
  Model <- NULL

  # empty variables for
  fitness <- NULL

  # iterate over all fits and prepare the data frame
  for (i in 1:length(autofit)) {
    fit <- data.frame(Fitness=autofit[[i]]$fitness,
                      Index=seq(length(autofit[[i]]$fitness)),
                      Model=as.factor(i))

    fitness <- rbind(fitness, fit)
  }

  # plot the results
  ggplot(data=fitness, aes(x=Index, y=Fitness, color=Model)) +
    geom_line(size=1) +
    ylab("Fitness") +
    xlab("Iteration") +
    scale_color_brewer(type="qual", palette="Set1")
}


#' @title plot_best_models
#' @description Plots the best fitted model for each of the specs used in autohrf.
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @param autofit Output of the autohrf function.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(event = c("encoding", "delay", "response"),
#'                      start_time = c(0, 2.65, 12.5),
#'                      end_time = c(3, 12.5, 16),
#'                      min_duration = c(1, 5, 1))
#'
#' model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
#'                      start_time = c(0, 2.5, 2.65, 12.5),
#'                      end_time = c(2.5, 3, 12.5, 15.5),
#'                      min_duration = c(1, 0.1, 5, 1))
#'
#' model_specs <- list(model3, model4)
#'
#' # run autohrf
#' df <- swm
#' autofit <- autohrf(df, model_specs, population=2, iter=2)
#'
#' # plot best models
#' plot_best_models(autofit)
#'
plot_best_models <- function(autofit) {
  # plot list storage
  graphs <- list()
  i <- 1

  # iterate over models
  for (af in autofit) {
    graphs[[i]] <- plot_events(af)
    i <- i + 1
  }

  # plot grid
  cowplot::plot_grid(plotlist=graphs, nrow=i-1, ncol=1, scale=0.95)
}


#' @title print_best_models
#' @description Prints the best fitted model for each of the specs used in autohrf.
#' @import utils
#' @export
#'
#' @param autofit Output of the autohrf function.
#'
#' @examples
#' # prepare model specs
#' model3 <- data.frame(event = c("encoding", "delay", "response"),
#'                      start_time = c(0, 2.65, 12.5),
#'                      end_time = c(3, 12.5, 16),
#'                      min_duration = c(1, 5, 1))
#'
#' model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
#'                      start_time = c(0, 2.5, 2.65, 12.5),
#'                      end_time = c(2.5, 3, 12.5, 15.5),
#'                      min_duration = c(1, 0.1, 5, 1))
#'
#' model_specs <- list(model3, model4)
#'
#' # run autohrf
#' df <- swm
#' autofit <- autohrf(df, model_specs, population=2, iter=2)
#'
#' # print best models
#' print_best_models(autofit)
#'
print_best_models <- function(autofit) {
  # iterate over models
  i <- 1
  cat("\n----------------------------------------\n")
  for (af in autofit) {
    cat("\nModel", i, "\n\n")
    i <- i + 1
    cat("Fitness: ", tail(af$fitness, 1), "\n\n")
    print(af$models[[1]])
    cat("\n----------------------------------------\n")
  }
}
