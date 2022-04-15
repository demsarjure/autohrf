# this is a development test script

# Load d
df <- swm

# prepare data frames
model3 <- data.frame(event = c("encoding", "delay", "response"),
                     start_time = c(0, 2.65, 12.5),
                     end_time = c(3, 12.5, 16))

model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
                     start_time = c(0, 2.5, 2.65, 12.5),
                     end_time = c(2.5, 3, 12.5, 15.5))

models <- list(model3, model4)

# auto
autofit <- autohrf(df, models, population = 2, iter = 2)

plot_fitness(autofit)

plot_best_models(autofit)

print_best_models(autofit)

# manual
d <- df
model_specs <- models
allow_overlap <- FALSE
population <- 100
iter <- 100
mutation_rate <- 0.1
mutation_factor <- 0.05
elitism <- 0.1
tr <- 2.5
f <- 100
method <- "middle"
hrf <- "boynton"
t <- 32
delta <- 2.25
tau <- 1.25
alpha <- 2
p <- c(6, 16, 1, 1, 6, 0, 32)
roi_weights <- data.frame(roi = c("L_LIPv_ROI", "L_SCEF_ROI", "R_p32pr_ROI"),
                          weight = c(10, 20, 30))
