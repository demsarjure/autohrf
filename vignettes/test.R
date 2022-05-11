# this is a development test script
library(ggplot2)

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
hrf <- "boynton"
t <- 32
delta <- 2.25
tau <- 1.25
alpha <- 2
p <- c(6, 16, 1, 1, 6, 0, 32)
roi_weights <- data.frame(roi = c("L_LIPv_ROI", "L_SCEF_ROI", "R_p32pr_ROI"),
                          weight = c(5, 5, 5))

model <- data.frame(event      = c("encoding", "delay", "response"),
                    start_time = c(0,           2.65,    12.5),
                    duration   = c(2.65,        9.85,    3))


# installation
install.packages("devtools")
library(devtools)
install_github("demsarjure/autohrf")
library(autohrf)

# load data
d <- swm

# manual evaluation
model <- data.frame(event      = c("encoding", "delay", "response"),
                    start_time = c(0,           2.65,    12.5),
                    duration   = c(2.65,        9.85,    3))
em <- evaluate_model(d, model)
plot_model(em)

ggsave(paste0("plot_model.pdf"),
       width = 3840,
       height = 1920,
       dpi = 500,
       units = "px")

plot_model(em, by_roi = TRUE)

ggsave(paste0("plot_model_by_roi.pdf"),
       width = 3840,
       height = 3840,
       dpi = 400,
       units = "px")

# autohrf
model3 <- data.frame(event = c("enconvolveding", "delay", "response"),
                     start_time = c(0, 2.65, 12.5),
                     end_time = c(3, 12.5, 16))
model4 <- data.frame(event = c("fixation", "target", "delay", "response"),
                     start_time = c(0, 2.5, 2.65, 12.5),
                     end_time = c(2.5, 3, 12.5, 15.5))
models <- list(model3, model4)
autofit <- autohrf(df, models, population = 2, iter = 2)
plot_fitness(autofit)
plot_best_models(autofit)
print_best_models(autofit)
