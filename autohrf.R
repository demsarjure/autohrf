# libraries
library(scales)
library(ez)
library(ggplot2)
library(gtools)


# Support functions ------------------------------------------------------------
downsample <- function(x, f=100, method="middle") {
  x <- as.matrix(x)
  l <- dim(x)[1]
  n <- ceiling(l / f)
  m <- matrix(0, n, dim(x)[2])
  for (i in 1:n) {
    start = (i-1)*f + 1
    end   = i*f
    if (end > l) end <- l
    if (method == "middle") {
      middle <- round(start + f/2)
      if (middle > l) middle <- l
      m[i,] = x[middle, ]
    }
    else if (method == "mean") {
      m[i,] = apply(as.matrix(x[start:end,]), 2, FUN=mean)
    }
  }
  return(m)
}


hrf_boynton <- function(TR=2.5, t=32, delta=2.25, tau=1.25, alpha=2) {
  t <- c(0:(t/TR)) * TR
  r <- (t - delta) / tau
  h <- ((r ^ alpha) * exp(-r))
  h[t < delta] <- 0
  h <- h / ((alpha ^ alpha) * exp(-alpha))
  return (h)
}


hrf_spm <- function(TR=2.5, t=16, p=c()) {
  if (length(p)==0) {
    p <- c(6, 16, 1, 1, 6, 0, 32)
  }
  dt    <- TR/t
  u     <- c(0:(p[7]/dt)) - p[6]/dt;
  hrf   <- dgamma(u*(dt/p[3]), p[1]/p[3]) - dgamma(u*(dt/p[4]), p[2]/p[4])/p[5]
  hrf   <- hrf[c(0:(p[7]/TR)) * t + 1];
  hrf   <- hrf / max(hrf)
  return (hrf)
}


conv_hrf <- function(x, TR=2.5, p=c(), hrf="spm", value="peak") {
  if (hrf=="spm") {
    ts <- hrf_spm(TR, p=p)
  }
  else if (hrf=="boynton") {
    ts <- hrf_boynton(TR)
  }
  pad <- length(ts)+20
  if (is.matrix(x)) {
    m  <- matrix(0, dim(x)[1]+2*pad, dim(x)[2])
    m[(pad+1):(pad+dim(x)[1]),] <- x
    ts <- filter(m, filter=ts, method="convolution", sides=1)[(pad+1):(pad+dim(x)[1]),]
    if (value == "peak") {
      maxv <- apply(ts, 2, FUN=function(x) max(abs(x)))
      ts <- ts / matrix(maxv, nrow=dim(ts)[1], ncol=dim(ts)[2], byrow=TRUE)
    }
  } else {
    ts <- filter(c(c(1:pad)*0, x, c(1:pad)*0), filter=ts, method="convolution", sides=1)[(pad+1):(pad+length(x))]
    if (value == "peak") {
      ts <- ts / max(abs(ts))
    }
  }
  return (ts)
}


conv_events <- function(events, TR=2.5, length=max(events$time)+max(events$duration)+30, p=c(), hrf="boynton", method="middle", f=100, value="scale") {
  eTR = TR/f
  length  <- ceiling(length/TR)*TR
  tslen   <- ceiling(length/eTR)
  nevents <- dim(events)[1]
  
  ts <- c(1:tslen)*0
  m  <- matrix(0, tslen, nevents)
  
  for (r in 1:nevents) {
    start = round(events$time[r]/eTR) + 1
    end   = round((events$time[r] + events$duration[r])/eTR)
    for (t in start:end) {
      ts[t] = ts[t] + events$value[r]
      m[t,r] = events$value[r]
    }
  }
  return (list(m=downsample(m, f, method), X=downsample(conv_hrf(m, eTR, p=p, hrf=hrf, value=value), f, method), ts=downsample(conv_hrf(ts, eTR, p=p, hrf=hrf, value=value), f, method), TR=TR, events=events))
}


plot_events <- function(r) {
  plot(r$ts, type="l", col="red")
  for (n in 1:dim(r$m)[2]) {
    lines(r$m[,n]*0.1, col="gray")
    lines(r$X[,n], col="green")
  }
  lines(r$ts, type="l", col="red")
}


ggplot_events <- function(r) {
  tslen  <- dim(r$X)[1]
  tstime <- c(0:(tslen-1)) * r$TR
  d <- data.frame(y=r$ts, x=tstime, ts="ts")
  if ("event" %in% names(r$events)) {
    for (n in 1:dim(r$X)[2]) {
      d <- rbind(d, data.frame(y=r$X[,n], x=tstime, ts=r$events$event[n]))
    }
    r$events$ts <- r$events$event
  } else {
    for (n in 1:dim(r$X)[2]) {
      d <- rbind(d, data.frame(y=r$X[,n], x=tstime, ts=paste("e", n, sep="")))
    }
    r$events$ts <- paste("e", c(1:dim(r$events)[1]), sep="")
  }
  p <- ggplot() +
    geom_line(data=d[d$ts=="ts",], aes(x=x,y=y), color="black", size=1, alpha=0.5) +
    geom_line(data=d[d$ts != "ts",], aes(x=x, y=y, color=ts, group=ts)) +
    geom_rect(data=r$events, aes(x=1, y=1, xmin=time, xmax=time+duration, ymax=0, ymin=-0.1*value, color=ts, fill=ts)) +
    ylab("") + xlab("time")
  if (!"event" %in% names(r$events)) p <- p + theme(legend.position = "none")  else p <- p + scale_fill_discrete(name="event") + scale_color_discrete(name="event")
  
  print(p)
}


# Modelling --------------------------------------------------------------------
runmodel <- function(d, r, normalize=TRUE, report=TRUE) {
  
  # --- Set up variables
  
  coeffs  <- c()
  fit     <- c()
  rois    <- unique(d$roi)
  events  <- as.character(r$events$event)
  nevents <- length(events)
  
  # --- Expand dataframe
  
  d[, c("(Intercept)", events, "y", "r")] <- 0
  d[, "(Intercept)"] <- 1
  
  l = dim(d[d$roi==rois[1],])[1]
  
  # --- Run through ROI
  
  for (roi in rois) {
    
    # --- Normalize if needed
    
    if (normalize) {
      r$X[1:l,] <- r$X[1:l,] / matrix(apply(r$X[1:l,], 2, FUN=function(x) max(abs(x))), nrow=l, ncol=nevents, byrow=TRUE)
    }
    
    # --- Compute linear model
    
    d[d$roi == roi, events] <- r$X[1:l,]
    m <- lm(formula(d[, c("x", events)]), d[d$roi == roi,] )
    
    # --- Save component timeseries
    
    d[d$roi == roi, events]        <- r$X[1:l,] * matrix(m$coefficients[events], l, length(r$events$event), byrow=TRUE)
    d[d$roi == roi, "(Intercept)"] <- m$coefficients["(Intercept)"][[1]]
    d[d$roi == roi, "y"]           <- apply(as.matrix(d[d$roi == roi, c("(Intercept)", events)]), 1, FUN=sum)
    d[d$roi == roi, "r"]           <- m$residuals
    d[d$roi == roi, events]        <- r$X[1:l,] * matrix(m$coefficients[events], l, length(r$events$event), byrow=TRUE) + m$coefficients["(Intercept)"][[1]]
    
    R2 <- 1 - var(m$residuals)/var(d[d$roi == roi, "x"])
    
    coeffs <- rbind(coeffs, data.frame(c(r=roi, as.list(m$coefficients[events]), R2=R2)))
  }
  for (v in c("x", events, "y")) {
    t <- d[, c("roi", "t", v)]
    t$event <- v
    names(t) <- c("roi", "t", "y", "event")
    fit <- rbind(fit, t)
  }
  
  fit$event <- factor(fit$event, levels=c("x", events, "y"), ordered=TRUE)

  r = list(mean=mean(coeffs$R2), median=median(coeffs$R2), min=min(coeffs$R2))
  
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

# Automated parameter search ---------------------------------------------------
autohrf <- function(d,
                    models,
                    population = 100,
                    iter = 100,
                    mutation_rate = 0.1,
                    mutation_factor = 0.05,
                    elitism = 0.1,
                    method = "middle",
                    value = "peak",
                    hrf = "spm") {
  
  # parameters
  pop <- population
  m_rate <- mutation_rate
  m_factor <- mutation_factor * population
  elitism <- round(population * elitism)
  
  # results
  results <- list()
  
  # iterate over all models
  n_models <- length(models)
  total_iterations <- n_models * iter
  execution_time <- NULL
  for (m in 1:n_models) {
    # get model
    current_model <- models[[m]]
    
    # generate starting
    start_time <- list()
    end_time <- list()
    n_events <- nrow(current_model)
    for (p in 1:pop) {
      # varibales to start start and end times
      starts <- NULL
      ends <- NULL
      for (i in 1:n_events) {
        # get event
        event <- current_model[i, ]
        
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
      start_time[[p]] <- starts
      end_time[[p]] <- ends
    }
    
    # iterate over generations
    max_fitness <- vector()
    for (i in 1:iter) {
      # calculate eta 
      current_iteration <- i + ((m-1) * iter)
      if (!is.null(execution_time)) {
        minutes <- difftime(Sys.time(), execution_time, units="mins")
        eta <- round(minutes * (total_iterations - current_iteration + 1), 1)
      } else {
        eta <- "x"
      }
      
      # time of execution
      execution_time <- Sys.time()
      
      # print
      cat("Progress:\t", current_iteration, "/", total_iterations, "\teta:", eta, "min\n")
      
      # evaluate each model
      fitness <- vector()
      for (p in 1:pop) {
        # create model from data
        model <- data.frame(event = current_model$event,
                            time = start_time[[p]],
                            duration = end_time[[p]] - start_time[[p]],
                            value = rep(1, n_events))
        
        r <- conv_events(model, 2.5, method=method, value=value, hrf=hrf)
        r2 <- runmodel(d, r, report=FALSE)$r$mean
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
        for (p in 1:elitism) {
          new_start[[p]] <- start_time[[p]]
          new_end[[p]] <- end_time[[p]]
        }
        
        # create new ones with genetic algorithms
        for (p in (elitism+1):pop) {
          # get parents
          p1 <- 1
          p2 <- 1
          for (j in 1:2) {
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
            if (j == 1) {
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
          for (j in 1:n_events) {
            # time
            rand <- runif(1)
            if (rand < m_rate) {
              # add some random value (depends on m_factor)
              start[j] <- start[j] + runif(1, -m_factor, m_factor)
              end[j] <- end[j] + runif(1, -m_factor, m_factor)
              
              if (end[j] < start[j]) {
                temp <- start[j]
                start[j] <- end[j]
                end[j] <- temp
              }
            }
          }
          
          # clamp events to boundaries
          for (j in 1:n_events) {
            # get event
            event <- current_model[j, ]
            
            end_limit <- event$end_time
            # get boundaries
            if (j == 1) {
              start_limit <- event$start_time
            } else {
              start_limit <- max(event$start_time, end[j-1])
            }
            
            # clamp to boundaries
            start[j] <- max(start_limit, min(end_limit-event$min_duration, start[j]))
            end[j] <- max(start[j]+event$min_duration, min(end_limit, end[j]))
            
            # prevent too short events
            if ((end[j] - start[j]) < event$min_duration) {
              end[j] = start[j] + event$min_duration
            }
          }

          # store
          new_start[[p]] <- start
          new_end[[p]] <- end
        }
        
        # replace old generation with new and repeat
        start_time <- new_start
        end_time <- new_end
      }
    }
    
    # construct models
    new_models <- list()
    for (p in 1:pop) {
      # create model from data
      model <- data.frame(event = current_model$event,
                          time = start_time[[p]],
                          duration = end_time[[p]] - start_time[[p]],
                          value = rep(1, n_events))
      
      new_models[[p]] <- model
    }
    
    results[[m]] <- list(models=new_models, fitness=max_fitness)
  }

  return(results)
}
