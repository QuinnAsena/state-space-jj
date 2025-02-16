# Data wrangling ---------------------------------------------------------

if (!require("pacman")) install.packages("pacman", repos="http://cran.r-project.org")
pacman::p_load(tidyverse, multinomialTS, forecast, mgcv)

sunfish_pollen <- read_csv("./data/Sunfish_585_pollen_counts_bchron062722B_cleaned.csv")
sunfish_combined <- read_csv("./data/Sunfish_dataset.csv")
sunfish_rioja <- read_csv("./data/Sunfish_pollen_rioja_2.csv")
sunfish_ll <- read_csv("./data/Sunfish_LL_Recon_Feb22.csv")
ns_temp <- read_csv("./data/Shuman24_NS_Temperatures.csv")
isotopes <- read_csv("./data/isotope_Kalmen_smoother_output.csv")

dim(sunfish_pollen)
dim(sunfish_rioja)

# CHECK THAT THIS IS THE CORRECT AGE MODEL
sunfish_pollen$age <- sunfish_rioja$Age

colnames(sunfish_pollen)
sunfish_pollen_only <- sunfish_pollen[-c(1:8)]
ncol(sunfish_pollen_only)
colSums(sunfish_pollen_only)
rowSums(sunfish_pollen_only)

# Remove spp with less than 5 count in whole core
filter_5_names <- names(which(colSums(sunfish_pollen_only) <= 5))

# Remove some aquatic/semi-aquatic spp as suggested
aquatic_spp <- c(
  "Brasenia", "Nuphar",
  "Nym", "Nym.cell", 
  "Potamoge","Myrica", 
  "C.stolon", "Botrich",
  "Botrioc")

filter_5_names[which(aquatic_spp %in% filter_5_names)]
spp_remove <- c(aquatic_spp, setdiff(filter_5_names, aquatic_spp))
colSums(sunfish_pollen_only[ ,spp_remove])

# clip age
sunfish_pollen_clean <- sunfish_pollen |>
  select(-c(all_of(spp_remove), DepthCM.Vania, Depth.predict,
            Age.bchron062722B, `_2.5`, `_10`, `_90`, `_975`)) |>
  rename("depth" = "Depth.Correct") |> 
  filter(age >= 2374 & age <= 8092) # Clip to age-span of isotopes

head(sunfish_pollen_clean)
tail(sunfish_pollen_clean)

colSums(sunfish_pollen_clean[ ,c("Pinus", "P.strobu", "P.diplo")])

sunfish_pollen_assigned <- sunfish_pollen_clean %>%
  mutate(
    total_pollen = P.strobu + P.diplo,
    P.strobu = if_else(total_pollen == 0, P.strobu, round(P.strobu + (Pinus * (P.strobu / total_pollen)))),
    P.diplo = if_else(total_pollen == 0, P.diplo, round(P.diplo + (Pinus * (P.diplo / total_pollen))))
  ) %>%
  select(-c(total_pollen, Pinus))

sum(colSums(sunfish_pollen_clean[ ,c("Pinus", "P.strobu", "P.diplo")]))
sum(colSums(sunfish_pollen_assigned[ ,c("P.strobu", "P.diplo")]))
colnames(sunfish_pollen_assigned)


sunfish_spp_long <- pivot_longer(sunfish_pollen_assigned, cols = -c(age, depth),
                                 names_to = "variablename")
#  mutate(variablename = replace(variablename, stringr::str_detect(variablename, "P.diplo*"), "Pinus")) |> 
#  group_by(variablename, age) |>
#  summarise(value = sum(value), .groups = 'keep')

# P.diplo now being aggregated with other
target_spp <- c("P.strobu", "Tsuga", "Fagus", "Quercus", "Betula")
colSums(sunfish_pollen_assigned[ ,target_spp])
colSums(sunfish_pollen_assigned[, !colnames(sunfish_pollen_assigned) %in% c(target_spp, "age", "depth"), drop = FALSE])
# High Isoetes counts. Remove?
# What is the peak in other from?

# Pull out target species
sunfish_target_spp <- sunfish_spp_long |>
  filter(variablename %in% target_spp)

# Group all other spp
sunfish_other <- sunfish_spp_long |>
  filter(!variablename %in% target_spp) |>
  mutate(variablename = "other") |>
  group_by(variablename, age, depth) |>
  summarise(value = sum(value), .groups='keep')

# Plot spp
bind_rows(sunfish_other, sunfish_target_spp) |>
  arrange(age) |>
  ggplot(aes(x = age, y = value)) +
    geom_point() +
    geom_line() +
    scale_x_reverse() +
    facet_wrap(~ variablename, scales = "free")

# Stack and pivot to wide
# arrange from oldest to yongest
sunfish_spp_wide <- bind_rows(sunfish_other, sunfish_target_spp) |>
  pivot_wider(id_cols = age, names_from = variablename, values_from = value) |>
  arrange(desc(age))

# Create bins
bin_width <- 100
sunfish_bins <-
  cut(
    sunfish_spp_wide$age,
    breaks = seq(
      from = min(sunfish_spp_wide$age),
      to = max(sunfish_spp_wide$age + bin_width),
      by = bin_width
    ), include.lowest = T, labels = F)

# Check if data fall within a bin and need to be summed
diff(sunfish_bins) # bin 22

nrow(sunfish_spp_wide)

sunfish_spp_binned <- bind_cols(bins = sunfish_bins, sunfish_spp_wide) |>
  group_by(bins) |> # Group the data by the bins so that we calculate per time bin
  summarise(across(c(other, all_of(target_spp)), sum, na.rm = TRUE),
            age = mean(age, na.rm = TRUE)) |> 
              arrange(desc(age))
            
# lost one row to binning, not bad
nrow(sunfish_spp_binned)


## Wrangling covariates --------------------------------------------------

isotopes_clean <- isotopes |>
  (\(x) split(x, x$var))() |>
  map(\(df) df |>
        select(-var) |>
        rename_with(~ paste(unique(df$var), ., sep = "_"), starts_with("X"))
  ) |>
  reduce(full_join, by = "bin") |> 
  select(bin, Beech_d13c_X.smoothed, Hemlock_d13c_X.smoothed) |> 
  rename(
    hem_d13 = Hemlock_d13c_X.smoothed,
    beech_d13 = Beech_d13c_X.smoothed)

dim(isotopes_clean)
print(isotopes_clean, n = 58)
diff(isotopes_clean$bin)

dim(sunfish_spp_binned)
print(sunfish_spp_binned, n = 41)
diff(sunfish_spp_binned$bins)

ll_bins <-
  cut(
    sunfish_ll$Age,
    breaks = seq(
      from = min(sunfish_spp_wide$age),
      to = max(sunfish_spp_wide$age + bin_width),
      by = bin_width
    ), include.lowest = T, labels = F)

ll_mean <- cbind(bins = ll_bins, sunfish_ll) |> 
  drop_na(bins) |> 
  select(bins, Lake_Level_Below_Modern) |> 
  group_by(bins) |> 
  summarise(mean_ll = mean(Lake_Level_Below_Modern)) |> 
  arrange(desc(bins))

dim(ll_mean)
head(ll_mean)
tail(ll_mean)
ll_mean$bins
diff(ll_mean$bins)

ns_temp_bins <-
  cut(
    ns_temp$Age,
    breaks = seq(
      from = min(sunfish_spp_wide$age),
      to = max(sunfish_spp_wide$age + bin_width),
      by = bin_width
    ), include.lowest = T, labels = F)

ns_temp_mean <- cbind(bins = ns_temp_bins, ns_temp) |> 
  drop_na(bins) |> 
  select(bins, Smean) |> 
  group_by(bins) |> 
  summarise(mean_ns_temp = mean(Smean)) |> 
  arrange(desc(bins))

dim(ns_temp_mean)
head(ns_temp_mean)
tail(ns_temp_mean)
ns_temp_mean$bins
diff(ns_temp_mean$bins)


sunfish_all <- sunfish_spp_binned |> 
  full_join(isotopes_clean, by = c("bins" = "bin")) |> 
  full_join(ll_mean) |> 
  full_join(ns_temp_mean) |> 
  arrange(desc(bins))

print(sunfish_all, n = 58)

sunfish_all |>
  select(bins, beech_d13, hem_d13, mean_ll, mean_ns_temp) |> 
  pivot_longer(-bins) |>
  ggplot(aes(x = bins, y = value)) +
  geom_point() +
  geom_line() +
  scale_x_reverse() +
  facet_wrap(~ name, scales = "free")

sunfish_all_interp <- sunfish_all |> 
  mutate(across(c(beech_d13, hem_d13), forecast::na.interp))

sunfish_all_interp |>
  select(bins, beech_d13, hem_d13, mean_ll, mean_ns_temp) |> 
  pivot_longer(-bins) |>
  ggplot(aes(x = bins, y = value)) +
  geom_point() +
  geom_line() +
  scale_x_reverse() +
  facet_wrap(~ name, scales = "free")

# X_ll <- data.frame(y = blanding_ll$LakeElev.cm.below.mod, x = blanding_ll$age)
# X_ll_gam <- mgcv::gam(y ~ s(x, bs = "bs", k = nrow(X_ll)), method = "REML", data =  X_ll[nrow(X_ll):1, , drop = F])
# pred <- predict(X_ll_gam, newdata = data.frame(x = sunfish_pollen$sunfish_years[63:69]))


# Set-up models ----------------------------------------------------------

Y <- sunfish_all_interp |>
  select(other, all_of(target_spp)) |>
  as.matrix()

Tsample <- which(rowSums(Y) != 0)

# set up X with temp
X <- sunfish_all_interp |>
  select(beech_d13, hem_d13, mean_ll, mean_ns_temp) |>
  as.matrix() |>
  scale()

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

saveRDS(mod_list, "./results/pairwise.rds")

lapply(mod_list, \(x){x$C})
sapply(mod_list, \(x){x$opt.convergence})

mod_list <- readRDS("./results/pairwise.rds")

mod_list[[5]]$C

# tsuga:fagus
# tsuga:betula
# fagus:quercus

# fagus:pinus
# tsuga:quercus
# tsuga:pinus
# betula:pinus
# betula:quercus
