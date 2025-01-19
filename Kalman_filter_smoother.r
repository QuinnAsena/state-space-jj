source("./tony/TVARSS_14Nov24.R")

library(readr)
d <- read_csv("./data/Sunfish_dataset.csv")

# remove rows with deltas = NA
d <- d[!is.na(d$Hemlock_d13c),]

# create dummy variable for a step change in 5000
# d$shift <- d$Age < 5000

# create time going forward and reorder dataset
d$Time <- max(d$Age) - d$Age
d <- d[order(d$Time),]

# scale variables
d$CHAR <- scale(d$CHAR)
d$LL <- scale(d$LL)

# Adjust plot margins to create space on the right for the legend
par(mar = c(5, 4, 4, 12))  # Increases right margin

# Plot data
plot(scale(Hemlock_d13c) ~ Age, data = d, typ = "l", ylim = c(-3,3), ylab = "Scaled Values", xlab = "Time")
lines(scale(Beech_d13c) ~ Age, data = d, col = "darkgreen")
lines(LL ~ Age, data = d, col = "blue")
lines(CHAR ~ Age, data = d, col = "brown")

legend("topright", inset = c(-0.55, 0), legend = c("Hemlock_d13c", "Beech_d13c", "LL", "CHAR"),
       col = c("black", "darkgreen", "blue", "brown"), lty = 1, cex = 0.8, xpd = TRUE, bty = "n")

# Set up parameters -------------------------------------------------------

# number of lags in the autocorrelation function
p <- 3
# standard error of the measurement error. Setting su.fixed = 1 assumes that the measurement error has SD = 1. Setting su.fixed = NA will get the fitting to estimate su.
su.fixed <- 1

# It is possible to have the autoregression coefficients change through time by setting sb0.fixed = NA and/or sb.fixed = matrix(NA,1,p). Setting them to zero means the auoregression coeffients don't change.
sb0.fixed <- 0
sb.fixed <- matrix(0,1,p)


# Hemlock -----------------------------------------------------------------

ME <- rep(1, nrow(d))

mod0_hem <- TVARSS(X = d$Hemlock_d13c, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed
)

LL <- d$LL
temp <- d$Temperature

# Create a list of U datasets
U_list <- list(
  LL = scale(as.matrix(LL)),
  temp = scale(as.matrix(temp)),
  ll_temp = scale(cbind(LL = LL, temp = temp))
)

# Create and empty list that will be filled with models
mod_list_hem <- vector(mode = "list", length = length(U_list))
mod_fit_list_hem <- vector(mode = "list", length = length(U_list))

# Run a loop over the U datasets saving them into the empty list
for (i in seq_along(U_list)) {
  print(i)

  ssmod <- TVARSS(X = d$Hemlock_d13c, p = p, ME = ME, U = U_list[[i]], Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod0_hem$b0,
                  b.start = mod0_hem$b
  )
  mod_list_hem[[i]] <- ssmod

  mod_fit_list_hem[[i]] <- TVARSS_KalmanSmoother(ssmod)
}

# There is now a list object, each element contains a model
# To indext a model use double square brackets
# Select model from list
mod_list_hem[[1]]
mod_list_hem[[2]]
mod_list_hem[[3]]

# There is also a list of fitted models for the plots
mod_fit_list_hem[[1]]
mod_fit_list_hem[[2]]
mod_fit_list_hem[[3]]


plot(d$Hemlock_d13c ~ d$Age, xlab = "Cal Years BP", ylab = "Data (black), Filter (blue), and Smoother (red)")
lines(d$Hemlock_d13c ~ d$Age, lty = 2, col = "black")  # Original data as dashed line
lines(mod_fit_list_hem[[1]]$Y.filtered[,1] ~ d$Age, type = "l", col = "blue")  # Kalman filter
lines(mod_fit_list_hem[[1]]$Y.smoothed[,1] ~ d$Age, col = "red")  # Kalman smoother

# Likelihood ratio test loop
# This loop is not saving any outputs, it just prints the likelihood ratio
# of each model against the mod0
for (i in seq_along(mod_list_hem)) {
  print(c(
    deviance = 2 * (mod_list_hem[[i]]$logLik - mod0_hem$logLik),
    P = pchisq(
      2 * (mod_list_hem[[i]]$logLik - mod0_hem$logLik),
      df = ncol(U_list[[i]]),
      lower.tail = F
    )
  ))
}


### lake-level pre/post 5ka
# create dummy variable for a step change in 5000
d$shift <- d$Age < 5000
U_shift <- cbind(d$LL, d$LL * d$shift); colnames(U_shift) <- c("LL", "LL.after.5000")
U0 <- as.matrix(d$LL); colnames(U0) <- c("LL")

mod0_shift_hem <- TVARSS(X = d$Hemlock_d13c, p = p, ME = ME, U = U0, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed
)
mod0_shift_hem


# mod0_list <- vector(mode = "list", length = 5)
# mod0_list[[1]] <- mod0_shift_hem
#
# for (i in 1:(length(mod0_list)-1)) {
#   mod0_list[[i]]$sb0.fixed
#
#   mod0_list[[i+1]] <- TVARSS(X = d$Hemlock_d13c, p = p, ME = ME, U = U0, Tsamplefract = .9, show.fig = F, annealing = F,
#                            sb0.fixed = mod0_list[[i]]$sb0.fixed,
#                            sb.fixed = mod0_list[[i]]$sb.fixed,
#                            su.fixed = su.fixed
#   )
# }


mod1_shift_hem <- TVARSS(X = d$Hemlock_d13c, p = p, ME = ME, U = U_shift, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               #	su.fixed = su.fixed,
               b0.start = mod0_shift_hem$b0,
               b.start = mod0_shift_hem$b
)
mod1_shift_hem
# statistical test of mod1 vs. mod0. This is a standard likelihood ratio test.
c(deviance = 2*(mod1_shift_hem$logLik - mod0_shift_hem$logLik), P = pchisq(2*(mod1_shift_hem$logLik - mod0_shift_hem$logLik), df = ncol(U_shift) - ncol(U0), lower.tail = F))






# Beech -------------------------------------------------------------------


mod0_beech <- TVARSS(X = d$Beech_d13c, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
                   sb0.fixed = sb0.fixed,
                   sb.fixed = sb.fixed,
                   su.fixed = su.fixed
)


# Create and empty list that will be filled with models
mod_list_beech <- vector(mode = "list", length = length(U_list))
mod_fit_list_beech <- vector(mode = "list", length = length(U_list))

# Run a loop over the U datasets saving them into the empty list
for (i in seq_along(U_list)) {
  print(i)

  ssmod <- TVARSS(X = d$Beech_d13c, p = p, ME = ME, U = U_list[[i]], Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod0_beech$b0,
                  b.start = mod0_beech$b
  )
  mod_list_beech[[i]] <- ssmod

  mod_fit_list_beech[[i]] <- TVARSS_KalmanSmoother(ssmod)
}

# Select model from list
mod_list_beech[[1]]
mod_list_beech[[2]]
mod_list_beech[[3]]

# Select fitted models
mod_fit_list_beech[[1]]


plot(d$Beech_d13c ~ d$Age, xlab = "Cal Years BP", ylab = "Data (black), Filter (blue), and Smoother (red)")
lines(d$Beech_d13c ~ d$Age, lty = 2, col = "black")  # Original data as dashed line
lines(mod_fit_list_beech[[2]]$Y.filtered[,1] ~ d$Age, type = "l", col = "blue")  # Kalman filter
lines(mod_fit_list_beech[[2]]$Y.smoothed[,1] ~ d$Age, col = "red")  # Kalman smoother

# Likelihood ratio test
for (i in seq_along(mod_list_beech)) {
  print(c(
    deviance = 2 * (mod_list_beech[[i]]$logLik - mod0_beech$logLik),
    P = pchisq(
      2 * (mod_list_beech[[i]]$logLik - mod0_beech$logLik),
      df = ncol(U_list[[i]]),
      lower.tail = F
    )
  ))
}


### lake-level pre/post 5ka
# create dummy variable for a step change in 5000
d$shift <- d$Age < 5000
U_shift <- cbind(d$LL, d$LL * d$shift); colnames(U_shift) <- c("LL", "LL.after.5000")
U0 <- as.matrix(d$LL); colnames(U0) <- c("LL")

mod0_shift_beech <- TVARSS(X = d$Beech_d13c, p = p, ME = ME, U = U0, Tsamplefract = .9, show.fig = F, annealing = F,
                         sb0.fixed = sb0.fixed,
                         sb.fixed = sb.fixed,
                         #	su.fixed = su.fixed
)
mod0_shift_beech


mod1_shift_beech <- TVARSS(X = d$Beech_d13c, p = p, ME = ME, U = U_shift, Tsamplefract = .9, show.fig = F, annealing = F,
                         sb0.fixed = sb0.fixed,
                         sb.fixed = sb.fixed,
                         #	su.fixed = su.fixed,
                         b0.start = mod0_shift_beech$b0,
                         b.start = mod0_shift_beech$b
)
mod1_shift_beech
# statistical test of mod1 vs. mod0. This is a standard likelihood ratio test.
c(deviance = 2*(mod1_shift_beech$logLik - mod0_shift_beech$logLik), P = pchisq(2*(mod1_shift_beech$logLik - mod0_shift_beech$logLik), df = ncol(U_shift) - ncol(U0), lower.tail = F))



# Plotting ----------------------------------------------------------------

library(patchwork)
library(ggtext)
library(RColorBrewer)
library(tidyverse)

dat <- read_csv("./data/Sunfish_dataset.csv")
dat <- dat[!is.na(dat$Hemlock_d13c),]
dat$Time <- max(dat$Age) - dat$Age
dat <- dat[order(dat$Time),]


# kalman_filt <- lapply(seq_along(mod_fit_list_hem), \(f) {
#   data.frame(Y.filtered = mod_fit_list_hem[[f]]$Y.smoothed[,1],
#              age = dat$Age,
#              obs = dat$Hemlock_d13c)
# })
# kalman_filt
#
# names(kalman_filt) <- paste0("mod_", 1:length(kalman_filt))
# kalman_filt_long <- bind_rows(kalman_filt, .id = "mod") |>
#   mutate(cat = "filt") |>
#   pivot_longer(-c(age, mod, cat))
#
#
# ggplot(kalman_filt_long, aes(x = age, y = value, shape = name, colour = name)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap( ~ mod) +
#   theme_minimal()


hem_mod_1 <- bind_cols(
  age = dat$Age,
  ll = dat$LL,
  hem_d13 = dat$Hemlock_d13c,
  Y_smoothed = mod_fit_list_hem[[1]]$Y.smoothed[, 1]) |>
  pivot_longer(-age) |>
  mutate(facet = ifelse(name == "ll", "facet_ll", "facet_mod"))


ggplot(hem_mod_1, aes(x = age, y = value, color = name)) +
  geom_point() +
  geom_line() +
  facet_grid(facet ~ ., scales = "free_y") +
  # facet_wrap(~facet, scales = "free_y", nrow = 2) +
  labs(x = "Age", y = "Value", color = "Variable") +
  theme_minimal()


mod_1 <- bind_cols(
  age = dat$Age,
  temp = dat$Temperature,
  ll = -dat$LL,
  hem_d13 = dat$Hemlock_d13c,
  beech_d13 = dat$Beech_d13c,
  Y_smoothed_beech = mod_fit_list_beech[[1]]$Y.smoothed[, 1],
  Y_smoothed_hem = mod_fit_list_hem[[1]]$Y.smoothed[, 1]) |>
  pivot_longer(-age) |>
  mutate(facet = case_when(
    name == "ll" ~ "facet_ll",
    name == "temp" ~ "facet_temp",
    # name %in% c("ll", "temp") ~ "facet_cov",
    name %in% c("hem_d13", "Y_smoothed_hem") ~ "facet_hem",
    name %in% c("beech_d13", "Y_smoothed_beech") ~ "facet_beech"
  ),
  facet_1 = case_when(name == "hem_d13" ~ "observed",
                      name == "beech_d13" ~ "observed",
                      name == "ll" ~ "facet_ll",
                      name == "temp" ~ "facet_temp",
                      .default = "smoothed"))

ggplot(mod_1, aes(x = age, y = value, color = name)) +
  geom_point() +
  geom_line() +
  facet_grid(facet ~ ., scales = "free_y") +
  theme_minimal()



plot_hem <- ggplot(mod_1 %>% filter(facet == "facet_hem"), aes(x = age, y = value, color = name)) +
  geom_line() +
  geom_point() +
  labs(x = NULL,
       y = expression(paste(delta^{13}, "which element was it?")),
       colour = NULL) +
  theme_minimal() +
  ylim(c(-31, -21)) +
  # scale_colour_brewer(palette="Paired", labels = c("Observed", "Kalman smoother")) +
  scale_colour_manual(values = brewer.pal(11, "BrBG")[c(3, 10)],
                      labels = c("Observed", "Kalman smoother")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_markdown(),
        legend.position = c(0.1, 0.8))


plot_beech <- ggplot(mod_1 %>% filter(facet == "facet_beech"), aes(x = age, y = value, color = name)) +
  geom_line() +
  geom_point() +
  ylim(c(-31, -21)) +
  labs(x = NULL,
       y = expression(paste(delta^{13}, "which element was it?")),
       colour = NULL) +
  theme_minimal() +
  # scale_colour_brewer(palette="Paired", labels = c("Observed", "Kalman smoother")) +
  scale_colour_manual(values = brewer.pal(11, "BrBG")[c(3, 10)], labels = c("Observed", "Kalman smoother")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_markdown(),
        legend.position = c(0.1, 0.8))


plot_ll <- ggplot(mod_1 %>% filter(facet == "facet_ll"), aes(x = age, y = value)) +
  geom_line() +
  labs(x = NULL, y = "cm relative to modern") +
  scale_y_reverse() +
  theme_minimal() +
  theme(axis.text.x = element_blank())


plot_temp <- ggplot(mod_1 %>% filter(facet == "facet_temp"), aes(x = age, y = value)) +
  geom_line() +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  labs(x = "Age (kabp)", y = "degrees C") +
  theme_minimal()

patchplot <- (plot_hem / plot_beech / plot_ll / plot_temp) +
  plot_layout(heights = c(2, 2, 1, 1))

patchplot


# Best to specify dimensions in the save pathway
# Vector formats like eps or svg are best for plots
# You can open up an SVG in inkscape or photoshop and have control of all the elements
ggsave(filename = "a_file_pathway.svg", plot = patchplot,
       width = 9, hieght = 9,
       units = "in",
       device = "svg")



ggplot(mod_1 %>% filter(facet == "facet_hem" | facet == "facet_beech"), aes(x = age, y = value, color = name, linetype = facet)) +
  geom_line() +
  geom_point() +
  labs(x = NULL,
       y = expression(paste(delta^{13}, "which element was it?")),
       colour = NULL) +
  theme_minimal() +
  ylim(c(-31, -21)) +
  scale_colour_brewer(palette="Paired", labels = c("Observed", "Kalman smoother")) +
  # scale_colour_manual(values = brewer.pal(11, "BrBG")[c(3, 10)],
  #                     labels = c("Observed", "Kalman smoother")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_markdown(),
        legend.position = c(0.1, 0.8))



ggplot() +
  geom_line(data = mod_1 %>% filter(facet == "facet_hem" | facet == "facet_beech", facet_1 == "smoothed"),
            aes(x = age, y = value, color = name)) +
  geom_point(data = mod_1 %>% filter(facet == "facet_hem" | facet == "facet_beech", facet_1 == "observed"),
             aes(x = age, y = value, color = name)) +
  labs(x = NULL,
       y = expression(paste(delta^{13}, "which element was it?")),
       colour = NULL) +
  theme_minimal() +
  ylim(c(-31, -21)) +
  # scale_colour_brewer(palette="Paired", labels = c("Observed", "Kalman smoother")) +
  scale_colour_manual(values = brewer.pal(11, "BrBG")[c(3, 10, 3, 10)],
                      labels = c("Observed", "Kalman smoother")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_markdown(),
        legend.position = c(0.1, 0.8))


ggplot(mod_1 %>% filter(facet_1 == "smoothed"), aes(x = age, y = value, color = name)) +
  geom_line() +
  geom_point() +
  labs(x = NULL,
       y = expression(paste(delta^{13}, "which element was it?")),
       colour = NULL) +
  theme_minimal() +
  ylim(c(-31, -21)) +
  # scale_colour_brewer(palette="Paired", labels = c("Observed", "Kalman smoother")) +
  scale_colour_manual(values = brewer.pal(11, "BrBG")[c(3, 10)],
                      labels = c("Kalman smoother hemlock?", "Kalman smoother beech?")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_markdown(),
        legend.position = c(0.1, 0.8))
