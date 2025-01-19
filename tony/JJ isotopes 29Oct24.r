source("./tony/TVARSS_29Oct24.R")

d <- read.csv("./data/Sunfish_dataset.csv")
summary(d)

# remove rows with deltas = NA
d <- d[!is.na(d$Hemlock_d13c),]

# create dummy variable for a step change in 5000
d$shift <- d$Age < 5000

# create time going forward and reorder dataset
d$Time <- max(d$Age) - d$Age
d <- d[order(d$Time),]

# scale variables
d$CHAR <- scale(d$CHAR)
d$LL <- scale(d$LL)

# plot data
plot(scale(Hemlock_d13c) ~ Time, data = d, typ = "l", ylim = c(-3,3))
lines(scale(Beech_d13c) ~ Time, data = d, col = "darkgreen")
lines(LL ~ Time, data = d, col = "blue")
lines(CHAR ~ Time, data = d, col = "brown")


########################################
# set up parameters

# number of lags in the autocorrelation function
p <- 3

# standard error of the measurement error. Setting su.fixed = 1 assumes that the measurement error has SD = 1. Setting su.fixed = NA will get the fitting to estimate su.
su.fixed <- 1

# It is possible to have the autoregression coefficients change through time by setting sb0.fixed = NA and/or sb.fixed = matrix(NA,1,p). Setting them to zero means the auoregression coeffients don't change.
sb0.fixed <- 0
sb.fixed <- matrix(0,1,p)

# response variable X
X <- d$Hemlock_d13c
# X <- d$Beech_d13c

# explanayory variable U. U can contain multiple explanatory variables or can be omitted.
U <- cbind(d$CHAR, d$LL); colnames(U) <- c("CHAR","LL")
# U <- as.matrix(d$CHAR)
# U <- as.matrix(d$LL)

########################################

# fit the models with a Kalman filter
ME <- rep(1, nrow(d))
mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
	sb0.fixed = sb0.fixed,
	sb.fixed = sb.fixed,
	su.fixed = su.fixed
)
mod1 <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
	sb0.fixed = sb0.fixed,
	sb.fixed = sb.fixed,
	su.fixed = su.fixed,
	b0.start = mod0$b0,
	b.start = mod0$b
)
mod1
# statistical test of mod1 vs. mod0. This is a standard likelihood ratio test.
c(deviance = 2*(mod1$logLik - mod0$logLik), P = pchisq(2*(mod1$logLik - mod0$logLik), df = ncol(U), lower.tail = F))

# measure of the strength of autocorrelation (1 is very highly autocorrelated)
mod1$eigen.fitted[1]

# fit the trajectory with a Kalman smoother
fit <- TVARSS_KalmanSmoother(mod1)

# blue = Kalman filter, red = Kalman smoother
plot(X ~ d$Time, xlab = "Time", ylab="data (black), filter (blue), and smoother (red)")
lines(X ~ d$Time, lty=2)
lines(fit$y[1,] ~ d$Time, typ = "l", col = "blue")
lines(fit$ysmooth[1,] ~ d$Time, col = "red")

########################################

# check to see if there is a change in the relationship with LL before/after 5000BP
# compare these
U <- cbind(d$LL, d$LL * d$shift); colnames(U) <- c("LL", "LL.after.5000")
U0 <- as.matrix(d$LL); colnames(U0) <- c("LL")

mod0 <- TVARSS(X = X, p = p, ME = ME, U = U0, Tsamplefract = .9, show.fig = F, annealing = F,
	sb0.fixed = sb0.fixed,
	sb.fixed = sb.fixed,
	su.fixed = su.fixed
)
mod0
mod1 <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
	sb0.fixed = sb0.fixed,
	sb.fixed = sb.fixed,
	su.fixed = su.fixed,
	b0.start = mod0$b0,
	b.start = mod0$b
)
mod1
# statistical test of mod1 vs. mod0. This is a standard likelihood ratio test.
c(deviance = 2*(mod1$logLik - mod0$logLik), P = pchisq(2*(mod1$logLik - mod0$logLik), df = ncol(U) - ncol(U0), lower.tail = F))
