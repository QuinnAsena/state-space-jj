devtools::install_github("https://github.com/QuinnAsena/multinomialTS")

# Y contains all target species, X contains all drivers unscaled (as a check on lake level values)
# top of matrix is the oldest time/bin
Y <- readRDS("./Y_dat")
X <- readRDS("./X_dat")
X <- scale(X)


# To aggregate everything but tsuga in Y:
# Y[ ,1] <- rowSums(Y[ , c("other", "P.strobu", "Fagus", "Quercus", "Betula")])
# Y <- Y[ , c("other", "Tsuga")]

# TEMPERATURE IS POLLEN INFERRED, check final decision

p <- ncol(X) + 1 # Number of independent variables plus intercept
n <- ncol(Y)

V.fixed = diag(n) # Covariance matrix of environmental variation in process eq
# V.fixed = matrix(NA, n, n)
# V.fixed[1] = 1

B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

glmm_mod <- multinomialTS::mnGLMM(Y = Y[which(rowSums(Y) != 0),],
                                  X = X[which(rowSums(Y) != 0), ,drop = F],
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

B.fixed <- matrix(NA, ncol(X), n)
B.fixed[,1] <- 0
B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

# Set-up C without interactions
C.start.diag = .5 * diag(n)
C.fixed.diag <- C.start.diag
C.fixed.diag[C.fixed.diag != 0] <- NA


# Model with no interactions
start_time <- Sys.time()
mnTS_mod <- mnTS(Y = Y,
                 X = X, Tsample = Tsample,
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

tsoga_prop <- mnTS_mod$Y[ ,2] / rowSums(mnTS_mod$Y)
mnTS_mod$X
lizt <- scale(cbind(mnTS_mod$Y[ ,2] / rowSums(mnTS_mod$Y),
mnTS_mod$X))
matplot(lizt, type = 'p', pch = 19)

arima(tsoga_prop, xreg = mnTS_mod$X, order = c(1,0,0))

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
