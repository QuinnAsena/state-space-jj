source("TVARSS_28Jan25.R")

d <- read.csv("./data/Sunfish_dataset.csv")

# create time going forward and reorder dataset
d$Time <- round(max(d$Age) - d$Age)
d <- d[order(d$Time),]

# remove rows with deltas = NA
d <- d[!is.na(d$Hemlock_d13c),]

# create dummy variable for a step change in 5000
# d$shift <- d$Age < 5000

# scale variables
d$CHAR <- scale(d$CHAR)
d$LL <- scale(d$LL)

# bin data into 100-year bins starting in 2376
d$bin <- 1 + (d$Age - 2376) %/% 100
d$bin <- max(d$bin) - d$bin +1

b <- data.frame(bin = 1:55)
d <- merge(d, b, all = T)

# interpolate mising data for independent variables
var.list <- c("Age", "Time","LL","Temperature","CHAR")
for(i.var in var.list) d[,i.var] <- approx(x = d$bin, y = d[,i.var], xout = d$bin)$y

# check
par(mfrow = c(length(var.list),1))
for(i.var in var.list) plot(d$bin, d[,i.var], col = 1 + is.na(d$Hemlock_d13c), main = i.var)


# Adjust plot margins to create space on the right for the legend
par(mar = c(5, 4, 4, 12))  # Increases right margin

# Plot data
par(mfrow = c(1,1))
plot(scale(Hemlock_d13c) ~ Age, data = d, typ = "l", ylim = c(-3,3), ylab = "Scaled Values", xlab = "Age")
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


LL <- scale(d$LL)
temp <- scale(d$Temperature)
Age <- d$Age

# Create a list of U variables that can be used in the analyses
Age5000 <- as.numeric(Age < 5000)
Age5500 <- as.numeric(Age < 5500)
Age4800 <- as.numeric(Age < 4800)

U.list <- data.frame(
  LL = LL,
  temp = temp,
  LL_temp = LL * temp,
  Age5000 = Age5000,
  Age5500 = Age5500,
  Age4800 = Age4800,
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

# test
source("TVARSS_11Feb25.r")
p <- 2

# It is possible to have the autoregression coefficients change through time by setting sb0.fixed = NA and/or sb.fixed = matrix(NA,1,p). Setting them to zero means the auoregression coeffients don't change.
sb0.fixed <- 0
sb.fixed <- matrix(0,1,p)

X <- d[,"Hemlock_d13c"]

pick.list <- list(c("LL"), c("temp"), c("LL", "temp"))
U <- U.list[,pick.list[[1]], drop = F]

mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed)
mod0

mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  c.fixed = rep(0, ncol(U)),               
                  d.fixed = rep(NA, ncol(U)),
                  d.start = rep(.01, ncol(U)),
                  b0.start = mod0$b0,
                  b.start = mod0$b)
mod
dev <- 2*(mod$logLik - mod0$logLik)
df <- ncol(U)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))

P_indep_vars_TVARSS(mod)


run <- function(p){
	
	w <- list()
	count <- 0
	for(i.var in c("Hemlock_d13c", "Beech_d13c")){
		count <- count + 1
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print(i.var)
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	
		X <- d[,i.var]
		
		###############################################
		print("tests of LL and temp")
		
		pick.list <- list(c("LL"), c("temp"), c("LL", "temp"))
		for(i in 1:3){
			
			U <- U.list[,pick.list[[i]], drop = F]
		
			mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
			               sb0.fixed = sb0.fixed,
			               sb.fixed = sb.fixed,
			               su.fixed = su.fixed)
			mod0
#browser()
			mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
			                  sb0.fixed = sb0.fixed,
			                  sb.fixed = sb.fixed,
			                  su.fixed = su.fixed,
			                  c.fixed = rep(0, ncol(U)),               
			                  d.fixed = rep(NA, ncol(U)),
			                  b0.start = mod0$b0,
			                  b.start = mod0$b)
			show(pick.list[[i]])
			show(mod)
			
			dev <- 2*(mod$logLik - mod0$logLik)
			df <- ncol(U)
			show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		}
		
		###############################################
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("tests of interactions")
		
		U1 <- U.list[,c("LL","temp"), drop = F]
		U <- cbind(U1, inter = U.list[,"LL_temp"])
		
		mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
		               sb0.fixed = sb0.fixed,
		               sb.fixed = sb.fixed,
		               su.fixed = su.fixed)
		               
		mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  d.fixed = rep(NA, ncol(U1)),
		                  c.fixed = rep(0, ncol(U1)),               
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)
		                  
		mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U)),               
		                  d.fixed = rep(NA, ncol(U)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		
		show(mod)
		
		# the first test is for all variables in U together
		dev <- 2*(mod$logLik - mod0$logLik)
		df <- ncol(U)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		# the second test is for the variables in U that aren't in U1
		dev <- 2*(mod$logLik - mod1$logLik)
		df <- ncol(U) - ncol(U1)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		
		###############################################
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("tests of changes in time")
		
		# pick whichever pair you want to test; the code is set up to test the significance of the
		# variable(s) added to U1 to produce U.
		for(i in c("LL","temp")){
			if(i == "LL"){
				U1 <- U.list[,c("LL", "Age5000")]
				U <- cbind(U1, inter = U.list[,"LL_Age5000"])
			}else{
				U1 <- U.list[,c("temp", "Age5000")]
				U <- cbind(U1, inter = U.list[,"temp_Age5000"])
			}
			mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
			               sb0.fixed = sb0.fixed,
			               sb.fixed = sb.fixed,
			               su.fixed = su.fixed)
			               
			mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
			                  sb0.fixed = sb0.fixed,
			                  sb.fixed = sb.fixed,
			                  su.fixed = su.fixed,
			                  c.fixed = rep(0, ncol(U1)),               
			                  d.fixed = rep(NA, ncol(U1)),
			                  b0.start = mod0$b0,
			                  b.start = mod0$b)		
			                  
			mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
			                  sb0.fixed = sb0.fixed,
			                  sb.fixed = sb.fixed,
			                  su.fixed = su.fixed,
			                  c.fixed = rep(0, ncol(U)),               
			                  d.fixed = rep(NA, ncol(U)),
			                  b0.start = mod0$b0,
			                  b.start = mod0$b)		
			show(i)
			show(mod)
			
			# the first test is for all variables in U together
			dev <- 2*(mod$logLik - mod0$logLik)
			df <- ncol(U)
			show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
			
			# the second test is for the variables in U that are not in U1
			dev <- 2*(mod$logLik - mod1$logLik)
			df <- ncol(U) - ncol(U1)
			show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		}
		
		###############################################
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("tests of interactions with a change point")
		
		U1 <- U.list[,c("LL","temp", "Age5000", "LL_temp")]
		U2 <- U.list[,c("LL","temp", "Age5000", "LL_temp", "LL_Age5000", "temp_Age5000")]
		U <- cbind(U2, inter_change = U.list[,"LL_temp_Age5000"])
		
		mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
		               sb0.fixed = sb0.fixed,
		               sb.fixed = sb.fixed,
		               su.fixed = su.fixed)
		               
		mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U1)),               
		                  d.fixed = rep(NA, ncol(U1)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		

		mod2 <- TVARSS(X = X, p = p, ME = ME, U = U2, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U2)),               
		                  d.fixed = rep(NA, ncol(U2)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		
		                  
		mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U)),               
		                  d.fixed = rep(NA, ncol(U)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		
		show(mod)
		
		# test for all variables in the full mod
		dev <- 2*(mod$logLik - mod0$logLik)
		df <- ncol(U)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		# test for all variables in U2 that are not in U1
		dev <- 2*(mod2$logLik - mod1$logLik)
		df <- ncol(U2) - ncol(U1)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		# test for the variables in U that are not in U2
		dev <- 2*(mod$logLik - mod2$logLik)
		df <- ncol(U) - ncol(U2)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		show(c(mod$logLik, mod2$logLik, mod1$logLik, mod0$logLik))
		
		###############################################
		# plot the Kalman-filtered and smoothed data
		fit <- TVARSS_KalmanSmoother(mod)
		
		# blue = Kalman filter, red = Kalman smoother
		pdf(paste0(i.var," figure.pdf"), height = 4, width = 6)
			par(mfrow = c(1,1), mai = c(.8,.8,.1,.1))
			plot(X ~ time, data = fit, xlab = "Time", ylab="data (blk), filter (blu), smooth (red), pred (gr)")
			lines(X ~ time, data = fit, lwd=2)
			lines(X.filtered ~ time, data = fit, col = "blue", lwd=2)
			lines(X.smoothed ~ time, data = fit, col = "red", lwd=2)
			lines(X.mean ~ time, data = fit, col = "green", lwd=2)
			points(X ~ time, data = fit, lty=2)
			points(X.filtered ~ time, data = fit, typ = "l", col = "blue")
			points(X.smoothed ~ time, data = fit, col = "red")
			points(X.mean ~ time, data = fit, col = "green")
		dev.off()
		
		w[[count]] <- data.frame(var = i.var, bin = fit$time, X = fit$X, X.filtered = fit$X.filtered, X.smoothed = fit$X.smoothed)
	}
	return(w)
}
# run analyses
p <- 3
sb.fixed <- matrix(0,1,p)
output <- run(p)
write.csv(output, "isotope Kalmen smoother output p = 3.csv", row.names = F)

p <- 2
sb.fixed <- matrix(0,1,p)
output <- run(p)
write.csv(output, "isotope Kalmen smoother output p = 2.csv", row.names = F)
