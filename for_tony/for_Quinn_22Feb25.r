devtools::install_github("https://github.com/QuinnAsena/multinomialTS")
library(multinomialTS)

# Y contains all target species, X contains all drivers unscaled (as a check on lake level values)
# top of matrix is the oldest time/bin
Y <- readRDS("./for_tony/Y_dat")
X <- readRDS("./for_tony/X_dat")
X <- scale(X)

Tsample <- which(rowSums(Y) != 0)


# To aggregate everything but tsuga in Y:
Y[ ,1] <- rowSums(Y[ , c("other", "P.strobu", "Fagus", "Quercus", "Betula")])
Y <- Y[ , c("other", "Tsuga")]

# look only at LL
XX <- X[,"mean_ll", drop = F]

# include all independent variables
XX <- X

# TEMPERATURE IS POLLEN INFERRED, check final decision

p <- ncol(XX) + 1 # Number of independent variables plus intercept
n <- ncol(Y)

V.fixed = diag(n) # Covariance matrix of environmental variation in process eq
# V.fixed = matrix(NA, n, n)
# V.fixed[1] = 1

B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

glmm_mod <- multinomialTS::mnGLMM(Y = Y[which(rowSums(Y) != 0),],
                                  X = XX[which(rowSums(Y) != 0), ,drop = F],
                                  B.start = B.start, B.fixed = B.fixed,
                                  V.fixed = V.fixed)
summary(glmm_mod)

B0.start <- glmm_mod$B[1, , drop = F]
B.start <- glmm_mod$B[2:p, , drop = F]

sigma.start <- glmm_mod$sigma

V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
V.fixed[1] = 1

V.start <- glmm_mod$V
# V.start <- diag(diag(V.start))

B.fixed <- matrix(NA, ncol(XX), n)
B.fixed[,1] <- 0
B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

# Set-up C without interactions
C.start.diag = .5 * diag(n)
C.fixed.diag <- C.start.diag
C.fixed.diag[C.fixed.diag != 0] <- NA

# Model with no interactions
start_time <- Sys.time()
mnTS_mod <- mnTS(Y = Y,
                 X = XX, 
                 Tsample = 1:nrow(Y),
                 B0.start = B0.start, B0.fixed = B0.fixed,
                 B.start = B.start, B.fixed = B.fixed,
                 C.start = C.start.diag, C.fixed = C.fixed.diag,
                 V.start = V.start, V.fixed = V.fixed,
                 dispersion.fixed = 1, maxit.optim = 1e+6)
# maxit.optim is the max number of iterations the optimiser will complete before stopping.
# increase maxit.optim if the model needs a lot of time to fit.
end_time <- Sys.time()
end_time - start_time
summary(mnTS_mod)


#########################
# plotting

df <- data.frame(Y = Y, X = XX, tsuga_prop = (Y[ ,2] / rowSums(Y)), y.fitted = mnTS_mod$y, time = 1:nrow(Y))
plot(scale(tsuga_prop) ~ time, data = df, type = "l")
points(scale(tsuga_prop) ~ time, data = df)
lines(scale(XX[,"mean_ll"]) ~ time, data = df, col = "blue")

arima(tsuga_prop, xreg = scale(XX[,"mean_ll"]), order = c(1,0,0))

# test
source("TVARSS_11Feb25.r")

# It is possible to have the autoregression coefficients change through time by setting sb0.fixed = NA and/or sb.fixed = matrix(NA,1,p). Setting them to zero means the auoregression coeffients don't change.
sb0.fixed <- 0
sb.fixed <- matrix(0,1,1)

U <- scale(XX[,"mean_ll"])
mod.c <- TVARSS(X = tsuga_prop, p = 1, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  c.fixed = rep(NA, ncol(U)),               
                  d.fixed = rep(0, ncol(U)),
                  c.start = rep(.01, ncol(U)))
mod.d <- TVARSS(X = tsuga_prop, p = 1, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  c.fixed = rep(0, ncol(U)),               
                  d.fixed = rep(NA, ncol(U)),
                  d.start = rep(.01, ncol(U)))
mod.c
mod.d
summary(mnTS_mod)

## setup pairwise interactions -------------------------------------------

c_idx <- t(combn(2:n, 2))
mod_list <- vector(mode = "list", length = nrow(c_idx))

start_time <- Sys.time()

for (idx in 1:nrow(c_idx)) {

  C.start <- .5 * diag(n)
  C.start[c_idx[idx, 1], c_idx[idx, 2]] = C.start[c_idx[idx, 2], c_idx[idx, 1]] = .001
  C.fixed <- C.start
  C.fixed[C.fixed != 0] <- NA

  print(C.start)
  print(C.fixed)

  mnTS_mod_int <- tryCatch(
    {
      multinomialTS::mnTS(Y = Y, X = X, Tsample = Tsample,
                          B0.start = mnTS_mod$B0, B0.fixed = B0.fixed,
                          B.start = mnTS_mod$B, B.fixed = B.fixed,
                          C.start = C.start, C.fixed = C.fixed,
                          V.start = mnTS_mod$V, V.fixed = V.fixed,
                          dispersion.fixed = 1, maxit.optim = 1e+7)
    },
    error=function(cond) {
      message("Here's the original error message:")
      message(conditionMessage(cond))
      return(cond)
    }
  )
  mod_list[[idx]] <- mnTS_mod_int
}

end_time <- Sys.time()
end_time - start_time
