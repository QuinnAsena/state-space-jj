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
ME <- rep(1, nrow(d))
su.fixed <- 1

# It is possible to have the autoregression coefficients change through time by setting sb0.fixed = NA and/or sb.fixed = matrix(NA,1,p). Setting them to zero means the auoregression coeffients don't change.
sb0.fixed <- 0
sb.fixed <- matrix(0,1,p)

LL <- scale(d$LL)
temp <- scale(d$Temperature)
Age <- d$Age

# Create a list of U variables that can be used in the analyses
U.list <- data.frame(
  LL = LL,
  temp = temp,
  LL_temp = LL * temp,
  Age5000 = Age < 5000,
  Age5500 = Age < 5500,
  Age4800 = Age < 4800,
  LL_Age5000 = LL * Age5000,
  LL_Age5500 = LL * Age5500,
  LL_Age4800 = LL * Age4800,
  temp_Age5000 = temp * Age5000,
  temp_Age5500 = temp * Age5500,
  temp_Age4800 = temp * Age4800,
  LL_temp_Age5000 = LL * temp * Age5000,
  LL_temp_Age5500 = LL * temp * Age5500,
  LL_temp_Age4800 = LL * temp * Age4800
)

# -----------------------------------------------------------------
# Analyses: I set this up so you have to pick and choose the
# X and U variables. I like this better than doing all combinations,
# because you can go through each model separately to see 
# what is going on.
# -----------------------------------------------------------------

# Hemlock
X <- d$Hemlock_d13c

# Beech
X <- d$Beech_d13c

###############################################
# tests of LL and temp

# pick whichever one you want to test
U <- U.list[,c("LL")]
U <- U.list[,c("temp")]
U <- U.list[,c("LL", "temp")]

mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed)
mod0
mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod0$b0,
                  b.start = mod0$b)
mod
c(dev = 2*(mod$logLik - mod0$logLik), P = pchisq(2*(mod$logLik - mod0$logLik), df = ncol(as.matrix(U)), lower.tail = F))


###############################################
# tests of interactions

U1 <- U.list[,c("LL","temp")]
U <- cbind(U1, inter = U.list[,"LL_temp"])

mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed)
               
mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod0$b0,
                  b.start = mod0$b)
                  
mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = T, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod1$b0,
                  b.start = mod1$b)
mod

# the first test is for all variables in U together
dev <- 2*(mod$logLik - mod0$logLik)
df <- ncol(U)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))

# the second test is for the variables in U that aren't in U1
dev <- 2*(mod$logLik - mod1$logLik)
df <- ncol(U) - ncol(U1)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))


###############################################
# tests of changes in time

# pick whichever pair you want to test; the code is set up to test the significance of the
# variable(s) added to U1 to produce U.
U1 <- U.list[,c("LL", "Age5000")]
U <- cbind(U1, inter = U.list[,"LL_Age5000"])

U1 <- U.list[,c("temp", "Age5000")]
U <- cbind(U1, inter = U.list[,"temp_Age5000"])

U1 <- U.list[,c("LL", "temp", "Age5000")]
U <- cbind(U1, inter = U.list[,c("LL_Age5000", "temp_Age5000")])

mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed)
               
mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod0$b0,
                  b.start = mod0$b)
                  
mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod1$b0,
                  b.start = mod1$b)
mod

# the first test is for all variables in U together
dev <- 2*(mod$logLik - mod0$logLik)
df <- ncol(U)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))

# the second test is for the variables in U that are not in U1
dev <- 2*(mod$logLik - mod1$logLik)
df <- ncol(U) - ncol(U1)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))

###############################################
# tests of interactions with a change point

U0 <- U.list[,c("LL","temp", "Age5000", "LL_temp")]
U1 <- U.list[,c("LL","temp", "Age5000", "LL_temp","LL_Age5000", "temp_Age5000")]
U <- cbind(U1, inter_change = U.list[,"LL_temp_Age5000"])

mod0 <- TVARSS(X = X, p = p, ME = ME, U = U0, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed)
               
mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod0$b0,
                  b.start = mod0$b)
                  
mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = T, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  b0.start = mod1$b0,
                  b.start = mod1$b)
mod

# the first test is for all variables in U1 that are not in U0
dev <- 2*(mod1$logLik - mod0$logLik)
df <- ncol(U1) - ncol(U0)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))

# the second test is for all variables in U that are not in U0
dev <- 2*(mod$logLik - mod0$logLik)
df <- ncol(U) - ncol(U0)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))

# the third test is for the variables in U that are not in U1
dev <- 2*(mod$logLik - mod1$logLik)
df <- ncol(U) - ncol(U1)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))




###############################################
# plot the Kalman-filtered and smoothed data
fit <- TVARSS_KalmanSmoother(mod)

# blue = Kalman filter, red = Kalman smoother
plot(X ~ time, data = fit, xlab = "Time", ylab="data (black), filter (blue), and smoother (red)")
lines(X ~ time, data = fit, lty=2)
lines(X.filtered ~ time, data = fit, typ = "l", col = "blue")
lines(X.smoothed ~ time, data = fit, col = "red")


