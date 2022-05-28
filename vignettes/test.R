# this is a development test script
library(ggplot2)

# libs
library(devtools)
install_github("demsarjure/autohrf")
library(autohrf)

# Load data
df <- swm

# prepare data frames
model1 <- data.frame(event = c("encoding", "delay", "response"),
                     start_time = c(0, 2.5, 12.5),
                     end_time = c(2.5, 12.5, 15))

model2 <- data.frame(event = c("encoding", "delay", "response"),
                     start_time = c(0, 3, 10),
                     end_time = c(3, 10, 15))

model3 <- data.frame(event = c("encoding", "delay", "response"),
                     start_time = c(0, 5, 10),
                     end_time = c(5, 10, 15))

model4 <- data.frame(event = c("encoding", "test", "delay", "response"),
                     start_time = c(0, 2.5, 5, 10),
                     end_time = c(2.5, 5, 10, 15))

# autohrf
models <- list(model1, model2, model3, model4)
autofit1 <- autohrf(df, models, tr = 2.5, population = 10, iter = 5)
p1 <- plot_fitness(autofit1)
p1
ggsave(paste0("plot_fitness1.pdf"),
       width = 3840,
       height = 1920,
       dpi = 500,
       units = "px")

autofit <- autohrf(df, models, tr = 2.5, population = 100, iter = 100)
p <- plot_fitness(autofit)
p
ggsave(paste0("plot_fitness.pdf"),
       width = 3840,
       height = 1920,
       dpi = 500,
       units = "px")

plot_best_models(autofit1)

# extract models
m <- get_best_models(autofit1)
m1 <- m[[1]]
m2 <- m[[2]]

# evaluate m1 manually
em <- evaluate_model(df, m1, tr = 2.5, hrf = "spm")
plot_model(em)
plot_model(em, by_roi = TRUE)

# manual
d <- df
roi_weights <- NULL
allow_overlap <- FALSE
population <- 10
iter <- 10
mutation_rate <- 0.1
mutation_factor <- 0.05
elitism <- 0.1
tr <- 2.5
hrf <- "spm"
t <- 32
p_boynton <- c(2.25, 1.25, 2)
p_spm <- c(6, 16, 1, 1, 6, 0)
f <- 100

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

# single event
model <- data.frame(event      = "encoding",
                    start_time = 0,
                    duration   = 10)
em <- evaluate_model(swm, model, tr = 2.5)
plot_model(em)
plot_model(em, by_roi = TRUE,
           rois = c("L_LO1_ROI", "R_PFt_ROI"))

# manual evaluation
model <- data.frame(event      = c("encoding", "delay", "response"),
                    start_time = c(0,           2.65,    12.5),
                    duration   = c(2.65,        9.85,    3))
em <- evaluate_model(swm, model, tr = 2.5)
plot_model(em, by_roi = TRUE,
           rois = c("R_2_ROI", "L_LIPv_ROI", "L_SCEF_ROI", "R_p32pr_ROI"))

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

# manual example 2
# load sample data
d <- swm

# test model
model <- data.frame(event = c("encoding", "delay", "response"),
                    start_time = c(0, 0.15, 10),
                    duration = c(0.15, 9.85, 3))

# evaluate the model
em <- evaluate_model(d, model, hrf = "spm")

# plot fits
dataset <- "1A"
dataset <- "1B"
dataset <- "1C"

d <- read.delim(sprintf("data/cluster_ts_%s.tsv", dataset))

d$y <- d$mean
d$t <- (d$frame - 1) * 2.5

model1 <- data.frame(event        = c("onset", "block", "rest"),
                     start_time   = c(0,        0,       60),
                     end_time     = c(2,        60,      65),
                     min_duration = c(0.1,      55,      0.1))

model2 <- data.frame(event        = c("onset", "block"),
                     start_time   = c(0,       0),
                     end_time     = c(2,       60),
                     min_duration = c(0.1,     55))

model3 <- data.frame(event        = c("onset", "block", "rest"),
                     start_time   = c(0,        0,       55),
                     end_time     = c(5,        60,      65),
                     min_duration = c(0.1,      50,      0.1))

model_constraints <- list(model1, model2, model3)
library(autohrf)
autofit <- autohrf(d, model_constraints, tr = 2.5, population = 5, iter = 20)

plot_best_models(autofit)

model2 <- data.frame(event        = "block",
                     start_time   = 0,
                     end_time     = 60,
                     min_duration = 10)

model_constraints <- list(model2)
autofit <- autohrf(d, model_constraints, tr = 2.5, population = 10, iter = 10)
plot_best_models(autofit)

model2 <- data.frame(event        = "block",
                     start_time   = 0,
                     end_time     = 60,
                     min_duration = 55)
model_constraints <- list(model2)
autofit <- autohrf(d, model_constraints, tr = 2.5, population = 10, iter = 10)
plot_best_models(autofit)
