## TVARSS.R 

##   Time-Varying Autoregressive State Space model

## Copyright 2016 Anthony R. Ives

## This file is part of the R-package `earlywarnings'.
## See the file ../COPYING for licensing issues.
## TVARSS.R 

##   Time-Varying Autoregressive State Space model

## Copyright 2016 Anthony R. Ives

## This file is part of the R-package `earlywarnings'.
## See the file ../COPYING for licensing issues.

require("GenSA")

TVARSS <- function(X, p = 1, ME = NULL, U = NULL, b0.start = NA, b.start = array(NA, dim = p), se.start = NA, su.start = 0.01, sb0.start = .05, sb.start = array(0.05, dim = p), c.start = if(is.null(U)) NULL else array(NA, dim=dim(U)[2]), b0.fixed = NA, b.fixed = array(NA, dim = p), se.fixed = NA, su.fixed = NA, sb0.fixed = NA, sb.fixed = array(NA, dim = p), c.fixed = if(is.null(U)) NULL else array(NA, dim=dim(U)[2]), Tsamplefract = .25, annealing = F, show.count = 10^10, show.fig = T, maxit.BFGS = 10^4, maxit.SANN = 10^2, optim.BFGS.control = NULL, optim.SANN.control = NULL) {
	
	require(GenSA)
	
	####################################################
	# Begin TVARSS.ml
	TVARSS_ml <- function(par, X, ME, U, p, par.fixed, show.count = 10^10, show.fig = T) {

		Tmax <- dim(X)[1]

		par.full <- par.fixed
		par.full[is.na(par.fixed)] <- par

		b0 <- par.full[1]
		b <- par.full[2:(p+1)]
		se <- par.full[p+2]
		su <- par.full[p+3]
		sb <- par.full[(p+4):(p+4+p)]
		if(!is.null(U)){
			nu <- dim(U)[2]
			cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
		}		

		B0 <- b0		
		B <- as.matrix(rbind(array(0, dim=p), diag(1, p - 1, p)))
		B[1, ] <- b
		
		Se <- diag(0, p)
		Se[1, 1] <- se^2
		
		Su <- su^2	
		Sb <- diag(sb^2)
		
		S <- as.matrix(rbind(cbind(Se, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		Z <- matrix(0, 1, 2*p+1)
		Z[1,1] <- 1
		
		# Initial unconditional values
		x <- X[p:1]
		
		PP <- solve(diag(p*p)-kronecker(B,B)) %*% matrix(Se, nrow = p*p)		
		PP <- matrix(PP, p, p)
		PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		if(counter %% show.count == 0){
			# information used for plotting
			eiglist <- matrix(NA, Tmax, 1)
			eigB <- eigen(B)$values
			eiglist[p:1] <- max(abs(eigB))
			
			xlist <- matrix(NA, nrow=Tmax, ncol=1)
			xlist[p:1] <- x
			
			b0list <- matrix(NA, nrow=Tmax, ncol=1)
			b0list[p:1] <- b0
			
			blist <- matrix(NA, nrow=Tmax, ncol=p)
			blist[p:1] <- sum(b)
			
			PPlist <- matrix(NA, nrow=Tmax, ncol=dim(PP)[1]^2)
			sdX <- sd(X)
		}
		
		logFt <- 0
		vFv <- 0
		for(t in (p+1):Tmax) {
				
			# PREDICTION EQUATIONS
			
			B12 <- as.matrix(rbind(1 - sum(B[1,]), matrix(0, p - 1, 1)))
			B13 <- as.matrix(rbind(t(x) - B0, matrix(0, p-1, p)))
			
			BB <- as.matrix(rbind(cbind(B, B12, B13), cbind(matrix(0, 1, p), 1, matrix(0, 1, p)), cbind(matrix(0, p, p+1), diag(p))))
			PP <- BB %*% PP %*% t(BB) + S
			if(is.null(U)){
				x <- B0 + B %*% (x-B0)
			}else{				
				if(ncol(U) == 1) {
					x <- B0 + B %*% (x-B0) + as.numeric(U[t] * cc)
				}else{
					x <- B0 + B %*% (x-B0) + as.numeric(U[t, ] %*% cc)
				}
			}		
						
			# UPDATING EQUATIONS
			if(!any(is.na(X[t]))){
				FF <- Z %*% PP %*% t(Z) + Su * ME[t]
				invF <- 1/FF
				
				y <- matrix(c(x, B0, B[1,]), ncol=1)		
				v <- X[t] - Z %*% y
				
				y <- y + PP %*% t(Z) %*% invF %*% v
				PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
				
				x <- y[1:p]
				B0 <- y[p+1]
				B[1,] <- y[(p+2):length(y)]
				
				# TERMS OF LIKELIHOOD FUNCTION
				if(FF > 0) {
					logFt <- logFt + log(FF)
				}else{
					logFt <- 10^10
				}
				
				vFv <- vFv + t(v) %*% invF %*% v
				if(counter %% show.count == 0){
					eigB <- eigen(B)$values
					eiglist[t] <- max(abs(eigB))
					xlist[t] <- x[1]
					PPlist[t,] <- matrix(PP, nrow=1)
					b0list[t] <- B0
					blist[t,] <- B[1,]
				}
			}
		}
		
		LL <- logFt + vFv
		if(is.complex(LL)) LL <- 10^10

		if((counter %% show.count == 0) & (show.fig == T)){
			npar <- length(par)
			logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2
			
			par(mfrow=c(2,1), mar=c(4, 4, 2, .5), cex.axis=1, cex.lab=1, cex.main=1)
			plot(1:Tmax, X, typ="p", xlab="Time", ylab="X (o), X.fitted (b), and b0.fitted (g)", main=paste('logLik=', .001*round(1000*logLik)))
			lines(1:Tmax, xlist, col = "blue")
			if(any(is.na(X))) points(1:Tmax, xlist, col="blue", pch = 20, cex = .5)
			
			if(is.null(U)){
				lines(1:Tmax, b0list, col = "green")
				if(any(is.na(X))) points(1:Tmax, b0list, col="green", pch = 20, cex = .5)
			}else{				
				if(ncol(U) == 1) {
					lines(1:Tmax, b0list + U * as.numeric(cc), col = "green")
					if(any(is.na(X))) points(1:Tmax, b0list + U * as.numeric(cc), col = "green", pch = 20, cex = .5)
				}else{
					lines(1:Tmax, b0list + U %*% cc, col = "green")
					if(any(is.na(X))) points(1:Tmax, b0list + U %*% cc, col = "green", pch = 20, cex = .5)
				}
			}		
			
			
			main.text <- NULL
			for(i in 0:p) main.text <- paste(main.text, "   sb[", i, "]=", abs(.001*round(1000*sb[i+1])), sep="")
			
			plot(1:Tmax, eiglist, typ="l", col="red", xlab="Time", ylab=expression(lambda), ylim=c(0,1.3), main=eval(expression(main.text)))
			if(any(is.na(X))) points(1:Tmax, eiglist, col="red", pch = 20, cex = .5)
			lines(c(0, Tmax), c(1, 1), col="black")
		}		
		
		if((counter %% show.count == 0)){
			X.fitted <<- xlist
			b0.fitted <<- b0list
			b.fitted <<- blist
			PP.fitted <<- PPlist
			eigen.fitted <<- eiglist
		}		
		counter <<- counter + 1
		return(LL)
	}			
	# End TVARSS_ml
	####################################################

	if (var(X, na.rm = TRUE) == 0) {
		stop("The response (dependent variable) has no variation.")
	}
	
	X <- as.matrix(X)
	if(!is.null(U)) {
		U <- as.matrix(U)
		nu <- dim(U)[2]
	}
	
	if (is.null(ME)) {
		ME <- matrix(1, dim(X))
	}else{
		ME <- as.matrix(ME)
		if (length(ME) != length(X)) {
			stop("The measurement error matrix ME should have the same length as X.")
		}
	}

	Tmax <- dim(X)[1]	
	Tsample <- floor(Tsamplefract * Tmax)
	
	if(length(sb.start) != p) stop("The length of sb.start must be p.")

	ar.init <- arima(X, xreg = U, order = c(p, 0, 0))
	b0.init <- ar.init$coef[p+1]
	b.init <- ar.init$coef[1:p]
	s2 <- ar.init$sigma2
	se.init <- (s2/2)^.5
	su.init <- (s2/2)^.5
	sb0.init <- sb0.start
	sb.init <- sb.start
	if(!is.null(U)) {
		c.init <- ar.init$coef[(p+2):(p+1+nu)]
	}else{
		c.init <- NULL
	}
	par.init <- c(b0.init, b.init, se.init, su.init, sb0.init, sb.init, c.init) 
	par.start <- c(b0.start, b.start, se.start, su.start, sb0.start, sb.start, c.start)
	par.fixed <- c(b0.fixed, b.fixed, se.fixed, su.fixed, sb0.fixed, sb.fixed, c.fixed)

	# set up variables for fitting
	par.full <- par.init
	par.full[!is.na(par.start)] <- par.start[!is.na(par.start)]
	
	par.full[!is.na(par.fixed)] <- par.fixed[!is.na(par.fixed)]	
	par <- par.full[is.na(par.fixed)]

	if(is.null(optim.BFGS.control)) optim.BFGS.control = list(maxit = maxit.BFGS)
	if(is.null(optim.SANN.control)) optim.SANN.control = list(temp = 0.1, tmax = 10, maxit = maxit.SANN)
		
	counter <- 1
	fitted.values <- NULL

	if(annealing == T){
		if(!is.null(U)) {
			par.upper <- c(100, array(2, dim=p), 10, 10, 10, array(10, dim=p), array(100, dim=nu))
			par.lower <- c(-100, array(-2, dim=p), 0, 0, 0, array(0, dim=p), array(-100, dim=nu))
		}else{
			par.upper <- c(100, array(2, dim=p), 10, 10, 10, array(10, dim=p))
			par.lower <- c(-100, array(-2, dim=p), 0, 0, 0, array(0, dim=p))
		}
			
		par.upper <- par.upper[is.na(par.fixed)]
		par.lower <- par.lower[is.na(par.fixed)]
		
		optSANN <- GenSA(fn = TVARSS_ml, par = par, lower = par.lower, upper = par.upper, X = X, U = U, p = p, par.fixed = par.fixed, show.count = show.count, show.fig = show.fig, control=list(smooth = F, maxit = maxit.SANN))	
		par <- optSANN$par
	}

	opt <- optim(fn = TVARSS_ml, par = par, X = X, U = U, ME = ME, p = p, par.fixed = par.fixed, show.count = show.count, show.fig = show.fig, method = "Nelder-Mead", control = optim.BFGS.control)
	
	if(opt$convergence != 0) cat("/nNelder-Mead failed to converge")
	
	# retrieve final fitted values
	TVARSS_ml(opt$par, X = X, U = U, ME = ME, p, par.fixed, show.count = 1, show.fig = show.fig)

	par.full <- par.fixed
	par.full[is.na(par.fixed)] <- opt$par

	b0 <- par.full[1]
	b <- par.full[2:(p+1)]
	se <- abs(par.full[p+2])
	su <- abs(par.full[p+3])
	sb0 <- abs(par.full[p+4])
	sb <- abs(par.full[(p+5):(p+4+p)])
	if(!is.null(U)){
		cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
	}else{
		cc <- NULL
	}		

	LL <- opt$value
	npar <- length(par)
	logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2
	AIC <- -2*logLik + 2*npar;

	results <- list(X = X, p = p, ME = ME, U = U, se = se, su = su, sb0 = sb0, sb = sb, cc = cc, b0 = b0, b = b, logLik = logLik, AIC = AIC, npar = npar, X.fitted = X.fitted, b0.fitted = b0.fitted, b.fitted = b.fitted, PP.fitted = PP.fitted, eigen.fitted = eigen.fitted, b0.start = b0.start, b.start = b.start, se.start = se.start, su.start = su.start, sb.start = sb.start, c.start = c.start, b0.fixed = b0.fixed, b.fixed = b.fixed, se.fixed = se.fixed, su.fixed = su.fixed, sb0.fixed = sb0.fixed, sb.fixed = sb.fixed, c.fixed = c.fixed, opt.par = opt$par, par.full = par.full, par.fixed = par.fixed, Tsamplefract = Tsamplefract, annealing = annealing, optim.BFGS.control = optim.BFGS.control, optim.SANN.control = optim.SANN.control)
	class(results) <- "TVARSS"
	
	rm(list=c("X.fitted", "b0.fitted", "b.fitted", "PP.fitted", "eigen.fitted"), envir = .GlobalEnv)

	return(results)
}


######################################################
######################################################
# summary.TVARSS
######################################################
######################################################

summary.TVARSS <- function(x, ...) {

	cat("\nCall: TVARSS with p = ", x$p, "\n")

	cat("\nlogLik =", x$logLik)
	cat(",  AIC = ", x$AIC, " [df = ",x$npar, "]\n", sep="")

	cat("\nAutoregression coefficients\n")
	cat("\tb[0] = ", x$b0, "\n", sep="")
	for(i in 1:length(x$b))
		cat("\tb[",i,"] = ", x$b[i], "\n", sep="")

	cat("\nVariation in autoregression coefficients (as standard deviations)\n")
	cat("\tsb[0] = ", x$sb0, "\n", sep="")
	for(i in 1:length(x$sb))
		cat("\tsb[",i,"] = ", x$sb[i], "\n", sep="")

	cat("\nProcess and Measurement errors (as standard deviations)")
	cat("\n\tse =", x$se)
	if(!is.null(x$su)) cat("\n\tsu =", x$su)
	
	if(!is.null(x$cc)){
		cat("\n\nCoefficients for U\n")
		for(i in 1:length(x$cc))
			cat("\tc[",i-1,"] = ", x$cc[i], "\n", sep="")	
	}
	
	cat("\n")
}
print.TVARSS <- summary.TVARSS

######################################################
######################################################
# plot.TVARSS
######################################################
######################################################

plot.TVARSS <- function(x, ...) {

	####################################################
	# Begin TVARSS.fig
	TVARSS.fig <- function(par, X, ME, U, p, par.fixed) {

		Tmax <- dim(X)[1]

		par.full <- par.fixed
		par.full[is.na(par.fixed)] <- par

		b0 <- par.full[1]
		b <- par.full[2:(p+1)]
		se <- par.full[p+2]
		su <- par.full[p+3]
		sb <- par.full[(p+4):(p+4+p)]
		if(!is.null(U)){
			nu <- dim(U)[2]
			cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
		}		

		B0 <- b0		
		B <- as.matrix(rbind(array(0, dim=p), diag(1, p - 1, p)))
		B[1, ] <- b
		
		Se <- diag(0, p)
		Se[1, 1] <- se^2
		
		Su <- su^2	
		Sb <- diag(sb^2)
		
		S <- as.matrix(rbind(cbind(Se, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		Z <- matrix(0, 1, 2*p+1)
		Z[1,1] <- 1
		
		# Initial unconditional values
		x <- X[p:1]
		
		PP <- solve(diag(p*p)-kronecker(B,B)) %*% matrix(Se, nrow = p*p)		
		PP <- matrix(PP, p, p)
		PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		# information used for plotting
		eiglist <- matrix(NA, Tmax, 1)
		eigB <- eigen(B)$values
		eiglist[1:p] <- max(abs(eigB))
		
		xlist <- matrix(NA, nrow=Tmax, ncol=1)
		xlist[p:1] <- x
		
		b0list <- matrix(NA, nrow=Tmax, ncol=1)
		b0list[p:1] <- b0
		
		blist <- matrix(NA, nrow=Tmax, ncol=p)
		blist[p:1] <- sum(b)
		
		PPlist <- matrix(NA, nrow=Tmax, ncol=dim(PP)[1]^2)
		sdX <- sd(X)
		
		logFt <- 0
		vFv <- 0
		for(t in (p+1):Tmax) {
				
			# PREDICTION EQUATIONS
			
			B12 <- as.matrix(rbind(1 - sum(B[1,]), matrix(0, p - 1, 1)))
			B13 <- as.matrix(rbind(t(x) - B0, matrix(0, p-1, p)))
			
			BB <- as.matrix(rbind(cbind(B, B12, B13), cbind(matrix(0, 1, p), 1, matrix(0, 1, p)), cbind(matrix(0, p, p+1), diag(p))))
			PP <- BB %*% PP %*% t(BB) + S
			if(is.null(U)){
				x <- B0 + B %*% (x-B0)
			}else{				
				if(ncol(U) == 1) {
					x <- B0 + B %*% (x-B0) + as.numeric(U[t] * cc)
				}else{
					x <- B0 + B %*% (x-B0) + as.numeric(U[t, ] %*% cc)
				}
			}		
			
			# UPDATING EQUATIONS
			if(!any(is.na(X[t]))){
				FF <- Z %*% PP %*% t(Z) + Su * ME[t]
				invF <- 1/FF
				
				y <- matrix(c(x, B0, B[1,]), ncol=1)		
				v <- X[t] - Z %*% y
				
				y <- y + PP %*% t(Z) %*% invF %*% v
				PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
				
				x <- y[1:p]
				B0 <- y[p+1]
				B[1,] <- y[(p+2):length(y)]
				
				# TERMS OF LIKELIHOOD FUNCTION
				if(FF > 0) {
					logFt <- logFt + log(FF)
				}else{
					logFt <- 10^10
				}
				
				vFv <- vFv + t(v) %*% invF %*% v
	
				eigB <- eigen(B)$values
				eiglist[t] <- max(abs(eigB))
				xlist[t] <- x[1]
				PPlist[t,] <- matrix(PP, nrow=1)
				b0list[t] <- B0
				blist[t,] <- B[1,]
			}
		}
		
		LL <- logFt + vFv
		if(is.complex(LL)) LL <- 10^10

		npar <- length(par)
		logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2
		
		par(mfrow=c(2,1), mar=c(4, 4, 2, .5), cex.axis=1, cex.lab=1, cex.main=1)
		plot(1:Tmax, X, typ="p", xlab="Time", ylab="X (o), X.fitted (b), and b0.fitted (g)", main=paste('logLik=', .001*round(1000*logLik)))
		lines(1:Tmax, xlist, col = "blue")
		if(any(is.na(X))) points(1:Tmax, xlist, col="blue", pch = 20, cex = .5)
		
		if(is.null(U)){
			lines(1:Tmax, b0list, col = "green")
			if(any(is.na(X))) points(1:Tmax, b0list, col="green", pch = 20, cex = .5)
		}else{				
			if(ncol(U) == 1) {
				lines(1:Tmax, b0list + U * as.numeric(cc), col = "green")
				if(any(is.na(X))) points(1:Tmax, b0list + U * as.numeric(cc), col = "green", pch = 20, cex = .5)
			}else{
				lines(1:Tmax, b0list + U %*% cc, col = "green")
				if(any(is.na(X))) points(1:Tmax, b0list + U %*% cc, col = "green", pch = 20, cex = .5)
			}
		}		
		
		
		main.text <- NULL
		for(i in 0:p) main.text <- paste(main.text, "   sb[", i, "]=", abs(.001*round(1000*sb[i+1])), sep="")
		
		plot(1:Tmax, eiglist, typ="l", col="red", xlab="Time", ylab=expression(lambda), ylim=c(0,1.3), main=eval(expression(main.text)))
		if(any(is.na(X))) points(1:Tmax, eiglist, col="red", pch = 20, cex = .5)
		lines(c(0, Tmax), c(1, 1), col="black")
		
	}		
	# End TVARSS.fig
	####################################################

	# plot for par = opt.par
	TVARSS.fig(par=x$opt.par, X=x$X, U=x$U, p=x$p, par.fixed=x$par.fixed)
}

TVARSS_KalmanSmoother <- function(mod){
	
	require(MASS)
	
	####################################################
	# Begin TVARSSsmoother_ml
	TVARSSsmoother_ml <- function(par, X, ME, U, p, par.fixed) {

		Tmax <- dim(X)[1]

		par.full <- par.fixed
		par.full[is.na(par.fixed)] <- par

		b0 <- par.full[1]
		b <- par.full[2:(p+1)]
		se <- par.full[p+2]
		su <- par.full[p+3]
		sb <- par.full[(p+4):(p+4+p)]
		if(!is.null(U)){
			nu <- dim(U)[2]
			cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
		}		

		B0 <- b0		
		B <- as.matrix(rbind(array(0, dim=p), diag(1, p - 1, p)))
		B[1, ] <- b
		
		Se <- diag(0, p)
		Se[1, 1] <- se^2
		
		Su <- su^2	
		Sb <- diag(sb^2)
		
		S <- as.matrix(rbind(cbind(Se, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		Z <- matrix(0, 1, 2*p+1)
		Z[1,1] <- 1
		
		# Initial unconditional values
		x <- X[p:1]
		
		PP <- solve(diag(p*p)-kronecker(B,B)) %*% matrix(Se, nrow = p*p)		
		PP <- matrix(PP, p, p)
		PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		B12 <- as.matrix(rbind(1 - sum(B[1,]), matrix(0, p - 1, 1)))
		B13 <- as.matrix(rbind(t(x) - B0, matrix(0, p-1, p)))
		
		BB <- as.matrix(rbind(cbind(B, B12, B13), cbind(matrix(0, 1, p), 1, matrix(0, 1, p)), cbind(matrix(0, p, p+1), diag(p))))

		y.t.t <- array(0, dim = c(2*p+1, Tmax))
		y.tp1.t <- array(0, dim = c(2*p+1, Tmax))
		BB.tp1.t <- array(0, dim = c(2*p+1, 2*p+1, Tmax))
		PP.t.t <- array(0, dim = c(2*p+1, 2*p+1, Tmax))
		PP.tp1.t <- array(0, dim = c(2*p+1, 2*p+1, Tmax))
		
		for(tt in 1:p){
			y.tp1.t[1:tt,tt] <- X[tt:1]
			y.t.t[1:tt,tt] <- X[tt:1]
		}
		BB.tp1.t[,,1:p] <- BB
		PP.t.t[,,1:p] <- PP
		PP.tp1.t[,,1:p] <- PP
		
		logFt <- 0
		vFv <- 0
		for(t in (p+1):Tmax) {
				
			# PREDICTION EQUATIONS
			
			B12 <- as.matrix(rbind(1 - sum(B[1,]), matrix(0, p - 1, 1)))
			B13 <- as.matrix(rbind(t(x) - B0, matrix(0, p-1, p)))
			
			BB <- as.matrix(rbind(cbind(B, B12, B13), cbind(matrix(0, 1, p), 1, matrix(0, 1, p)), cbind(matrix(0, p, p+1), diag(p))))
			PP <- BB %*% PP %*% t(BB) + S
			if(is.null(U)){
				x <- B0 + B %*% (x-B0)
			}else{				
				if(ncol(U) == 1) {
					x <- B0 + B %*% (x-B0) + as.numeric(U[t] * cc)
				}else{
					x <- B0 + B %*% (x-B0) + as.numeric(U[t, ] %*% cc)
				}
			}
			y <- matrix(c(x, B0, B[1,]), ncol=1)		
						
			y.tp1.t[,t-1] <- y
			BB.tp1.t[,,t-1] <- BB
			PP.tp1.t[,,t-1] <- PP
			
			# UPDATING EQUATIONS
			if(!any(is.na(X[t]))){
				FF <- Z %*% PP %*% t(Z) + Su * ME[t]
				invF <- 1/FF
				
				v <- X[t] - Z %*% y
				
				y <- y + PP %*% t(Z) %*% invF %*% v
				PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
				
				y.t.t[,t] <- y
				PP.t.t[,,t] <- PP

				x <- y[1:p]
				B0 <- y[p+1]
				B[1,] <- y[(p+2):length(y)]
				
				# TERMS OF LIKELIHOOD FUNCTION
				if(FF > 0) {
					logFt <- logFt + log(FF)
				}else{
					logFt <- 10^10
				}
				
				vFv <- vFv + t(v) %*% invF %*% v
			}			
		}
		
		LL <- logFt + vFv
		logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2
			
		return(list(y.t.t = y.t.t, y.tp1.t = y.tp1.t, BB.tp1.t = BB.tp1.t, PP.t.t = PP.t.t, PP.tp1.t = PP.tp1.t))
	}			
	# End TVARSSsmoother_ml
	####################################################
	
	
	Tmax <- dim(mod$X)[1]

	est <- TVARSSsmoother_ml(par = mod$opt.par, mod$X, mod$ME, mod$U, mod$p, mod$par.fixed)
	
	y.tp1.t <- est$y.tp1.t
	y.t.t <- est$y.t.t
	PP.t.t <- est$PP.t.t
	PP.tp1.t <- est$PP.tp1.t
	BB.tp1.t <- est$BB.tp1.t
		
	ys <- array(dim = dim(y.t.t))
	ys[,Tmax] <- y.t.t[,Tmax]
	#ys[,Tmax] <- y.t.tm1[,Tmax]
	for(t in (Tmax-1):1){
		J <- PP.t.t[,,t] %*% t(BB.tp1.t[,,t]) %*% ginv(PP.tp1.t[,,t])
		ys[,t] <- y.t.t[,t] + J %*% (ys[,t+1] - y.tp1.t[,t])
	}
	
	return(list(y = y.t.t, ysmoothed = ys, y.tp1.t = y.tp1.t))
		
}

