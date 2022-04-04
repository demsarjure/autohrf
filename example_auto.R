# Unassumed modelling

# source functions
source("autohrf.R")

# load data
d <- read.delim("./data/sWMsm_all_U_HCP_raw_long.txt")

# select specific ROIs
d <- d[d$roicode %in% c(1, 2, 10, 17, 20, 43, 45, 48, 50, 52, 60, 80, 84, 96, 116, 117, 138, 169,
                    181, 182, 190, 197, 200, 223, 225, 228, 230, 232, 240, 260, 264, 276, 296, 297, 318, 349),] # 78, 11

# exclude error/incorrect trials
d <- d[!d$event %in% c("e", "i"),]

# join data
dm <- aggregate(d$mean, list(roi=d$roi, t=(d$frame-1)*2.5), function(x) mean(x, na.rm=TRUE))

# prepare data frames
model3 <- data.frame(event = c("encoding", "delay", "response"),
                     start_time = c(0, 2.65, 12.5),
                     end_time = c(3, 12.5, 16),
                     min_duration = c(1, 5, 1))

model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
                     start_time = c(0, 2.5, 2.65, 12.5),
                     end_time = c(2.5, 3, 12.5, 15.5),
                     min_duration = c(1, 0.1, 5, 1))

model5a <- data.frame(event = c("fixation", "target", "delayA", "delayB", "response"),
                     start_time = c(0, 2.5, 2.65, 5, 12.5),
                     end_time = c(2.5, 3, 10, 12.5, 15.5),
                     min_duration = c(1, 0.1, 2.5, 2.5, 1))

model5b <- data.frame(event = c("fixation", "target", "delay", "responseA", "responseB"),
                      start_time = c(0, 2.5, 2.65, 12.5, 13),
                      end_time = c(2.5, 3, 12.5, 15, 15.5),
                      min_duration = c(1, 0.1, 5, 1, 1))

model6 <- data.frame(event = c("fixation", "target", "delayA", "delayB", "responseA", "responseB"),
                     start_time = c(0, 2.5, 2.65, 5, 12.5, 13),
                     end_time = c(2.5, 3, 10, 12.5, 15, 15.5),
                     min_duration = c(1, 0.1, 2.5, 2.5, 1, 1))

# merge models in a list
#models <- list(model3)
models <- list(model4)
#models <- list(model3, model4, model5b)
#models <- list(model3, model4, model5a, model5b, model6)

# auto fit
autofit <- autohrf(dm, models)

# fast fit - do not use this for final research
autofit <- autohrf(dm, models, population=25, iter=50)

# fast fit - boynton
autofit <- autohrf(dm, models, population=25, iter=50, hrf = "boynton")

# plot fitness (did the solution converge?)
fitness <- NULL
for (i in 1:length(autofit)) {
  fit <- data.frame(fitness=autofit[[i]]$fitness,
                    index=seq(length(autofit[[i]]$fitness)),
                    model=as.factor(i))
  
  fitness <- rbind(fitness, fit)
}
ggplot(data=fitness, aes(x=index, y=fitness, color=model)) +
  geom_line()

 # test best models
model1 <- autofit[[1]]$models[[1]]
model1
r1 <- conv_events(model1, 2.5, method="middle", value="peak")
ggplot_events(r1)
dm1 <- runmodel(dm, r1)

model2 <- autofit[[2]]$models[[1]]
model2
r2 <- conv_events(model2, 2.5, method="middle", value="peak")
ggplot_events(r2)
dm2 <- runmodel(dm, r2)

model3 <- autofit[[3]]$models[[1]]
model3
r3 <- conv_events(model3, 2.5, method="middle", value="peak")
ggplot_events(r3)
dm3 <- runmodel(dm, r3)

model4 <- autofit[[4]]$models[[1]]
model4
r4 <- conv_events(model4, 2.5, method="middle", value="peak")
ggplot_events(r4)
dm4 <- runmodel(dm, r4)
