# Load d
df <- swm

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

autofit <- autohrf(df, models, population=2, iter=2)

plot_fitness(autofit)

plot_best_models(autofit)

print_best_models(autofit)
