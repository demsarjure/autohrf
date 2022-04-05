# select specific ROIs
d <- swm[d$roicode %in% c(1, 2, 10, 17, 20, 43, 45, 48, 50, 52, 60, 80, 84, 96, 116, 117, 138, 169,
                        181, 182, 190, 197, 200, 223, 225, 228, 230, 232, 240, 260, 264, 276, 296, 297, 318, 349),]

d <- d[!d$event %in% c("e", "i"),]

# join data
tr <- 2.5
dm <- aggregate(d$mean, list(roi=d$roi, t=(d$frame-1)*tr), function(x) mean(x, na.rm=TRUE))

# prepare data frames
model3 <- data.frame(event = c("encoding", "delay", "response"),
                     start_time = c(0, 2.65, 12.5),
                     end_time = c(3, 12.5, 16),
                     min_duration = c(1, 5, 1))

model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
                     start_time = c(0, 2.5, 2.65, 12.5),
                     end_time = c(2.5, 3, 12.5, 15.5),
                     min_duration = c(1, 0.1, 5, 1))

models <- list(model3, model4)

autofit <- autohrf(dm, models, population=20, iter=10, hrf = "boynton")

plot_fitness(autofit)

plot_best_models(autofit)

print_best_models(autofit)
