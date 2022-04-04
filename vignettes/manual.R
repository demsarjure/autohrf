# Unassumed modelling

# source functions
source("autohrf.R")

# load data
d <- read.delim("./data/sWMsm_all_U_HCP_raw_long.txt")

# select specific ROIs
d <- d[d$roicode %in% c(1, 2, 10, 17, 20, 43, 45, 48, 50, 52, 60, 80, 84, 96, 116, 117, 138, 169,
                    181, 182, 190, 197, 200, 223, 225, 228, 230, 232, 240, 260, 264, 276, 296, 297, 318, 349),] # 78, 11


# view data
ggplot(d, aes(frame, mean)) + stat_summary(geom="line", fun.y="mean", size=0.8) + stat_summary(geom="errorbar", fun.data="mean_se", width=0.4)
ggplot(d, aes(frame, mean)) + stat_summary(geom="line", fun.y="mean", size=0.8) + stat_summary(geom="errorbar", fun.data="mean_se", width=0.4) + facet_wrap(~ roi)
ggplot(d, aes(frame, mean)) + stat_summary(geom="line", fun.y="mean", size=0.8) + stat_summary(geom="errorbar", fun.data="mean_se", width=0.4) + facet_wrap(~ roicode)

ggplot(d, aes(frame, mean, color=event)) + stat_summary(geom="line", fun.y="mean", size=0.8) + stat_summary(geom="errorbar", fun.data="mean_se", width=0.4)
ggplot(d, aes(frame, mean, color=event)) + stat_summary(geom="line", fun.y="mean", size=0.8) + stat_summary(geom="errorbar", fun.data="mean_se", width=0.4) + facet_wrap(~ roi)

# exclude error/incorrect trials
d <- d[!d$event %in% c("e", "i"),]

# join data
dm <- aggregate(d$mean, list(roi=d$roi, t=(d$frame-1)*2.5), function(x) mean(x, na.rm=TRUE))
de <- aggregate(d$mean, list(roi=d$roi, event=d$event, t=(d$frame-1)*2.5), function(x) mean(x, na.rm=TRUE))
ds <- aggregate(d$mean, list(roi=d$roi, subject=d$subject, t=(d$frame-1)*2.5), function(x) mean(x, na.rm=TRUE))
des <- aggregate(d$mean, list(roi=d$roi, event=d$event, subject=d$subject, t=(d$frame-1)*2.5), function(x) mean(x, na.rm=TRUE))

ggplot(dm, aes(t, x)) + geom_line() + facet_wrap(~ roi)
ggplot(de, aes(t, x, color=event)) + geom_line() + facet_wrap(~ roi)
ggplot(ds, aes(t, x)) + stat_summary(geom="line", fun.y="mean", size=0.8) + stat_summary(geom="errorbar", fun.data="mean_se", width=0.4) + facet_wrap(~ roi)
ggplot(des, aes(t, x, color=event)) + stat_summary(geom="line", fun.y="mean", size=0.8) + stat_summary(geom="errorbar", fun.data="mean_se", width=0.4) + facet_wrap(~ roi)


# Testing assumed models
# events: encoding, delay, response
modelA <- data.frame(event    = c("encoding", "delay", "response"),
                     time     = c(0,           2.65,    12.5     ),
                     duration = c(2.65,        9.85,    3        ),
                     value    = c(1,           1,       1        ))
rA <- conv_events(modelA, 2.5, method="middle", value="peak", hrf="boynton")
ggplot_events(rA)
dmA <- runmodel(dm, rA) # mean R: 0.897232, median R: 0.9204496


# events: fixation, target, delay, response
modelB <- data.frame(event    = c("fixation", "target", "delay", "response"),
                     time     = c(0,           2.5,      2.65,    12.5     ),
                     duration = c(2.5,         0.1,      9.85,    3        ),
                     value    = c(1,           1,        1,       1        ))
rB <- conv_events(modelB, 2.5, method="middle", value="peak")
ggplot_events(rB)
dmB <- runmodel(dm, rB) # mean R: 0.9214711, median R: 0.940774


# events: encoding, early delay, late delay, response
modelC <- data.frame(event    = c("encoding", "delay1", "delay2", "response"),
                     time     = c(0,           2.65,     7.5,      12.5     ),
                     duration = c(2.65,        4.85,     5,        3        ),
                     value    = c(1,           1,        1,        1        ))
rC <- conv_events(modelC, 2.5, method="middle", value="peak")
ggplot_events(rC)
dmC <- runmodel(dm, rC) # mean R: 0.9312625, median R: 0.9517234


# events: encoding, delay, probe, response
modelD <- data.frame(event    = c("encoding", "delay", "probe", "response"),
                     time     = c(0,           2.65,    12.5,    13       ),
                     duration = c(2.65,        9.85,    0.5,     2.5      ),
                     value    = c(1,           1,       1,       1        ))
rD <- conv_events(modelD, 2.5, method="middle", value="peak")
ggplot_events(rD)
dmD <- runmodel(dm, rD) # mean R: 0.9427942, median R: 0.9614632


# events: fixation, target, delay, probe, response

modelE <- data.frame(event    = c("fixation", "target", "delay", "probe", "response"),
                     time     = c(0,           2.5,      2.65,    12.5,    13        ),
                     duration = c(2.5,         0.1,      9.85,    0.5,     2.5       ),
                     value    = c(1,           1,        1,       1,       1         ))
rE <- conv_events(modelE, 2.5, method="middle", value="peak")
ggplot_events(rE)
dmE <- runmodel(dm, rE) # mean R: 0.9517975, median R: 0.9708709


# events: fixation, target, early and late delay, response
modelF <- data.frame(event    = c("fixation", "target", "delay1", "delay2", "response"),
                     time     = c(0,           2.5,      2.65,     7.5,      12.5     ),
                     duration = c(2.5,         0.1,      4.85,     5,        3        ),
                     value    = c(1,           1,        1,        1,        1        ))
rF <- conv_events(modelF, 2.5, method="middle", value="peak")
ggplot_events(rF)
dmF <- runmodel(dm, rF) # mean R: 0.9348051, median R: 0.9523274


# events: fixation, target, early and late delay, probe, response
modelG <- data.frame(event    = c("fixation", "target", "delay1", "delay2", "probe", "response"),
                     time     = c(0,           2.5,      2.65,     7.5,      12.5,    13       ),
                     duration = c(2.5,         0.1,      4.85,     5,        0.5,     2.5      ),
                     value    = c(1,           1,        1,        1,        1,       1        ))
rG <- conv_events(modelG, 2.5, method="middle", value="peak")
ggplot_events(rG)
dmG <- runmodel(dm, rG) # mean R: 0.9667715, median R: 0.9805873


# events: encoding, early and late delay, probe, response
modelH <- data.frame(event    = c("encoding", "delay1", "delay2", "probe", "response"),
                     time     = c(0,           2.65,     7.5,      12.5,    13        ),
                     duration = c(2.5,         4.85,     5,        0.5,     2.5       ),
                     value    = c(1,           1,        1,        1,       1         ))
rH <- conv_events(modelH, 2.5, method="middle", value="peak")
ggplot_events(rH)
dmH <- runmodel(dm, rH) # mean R: 0.9475901, median R: 0.9682923


# events: fixation, early and late delay, probe, response
modelI <- data.frame(event    = c("fixation", "delay1", "delay2", "probe", "response"),
                     time     = c(0,           2.5,      7.5,      12.5,     15        ),
                     duration = c(2.5,         5,        5,        1,        1         ),
                     value    = c(1,           1,        1,        1,        1         ))
rI <- conv_events(modelI, 2.5, method="middle", value="peak")
ggplot_events(rI)
dmI <- runmodel(dm, rI) # mean R: 0.9713018, median R: 0.9766944
