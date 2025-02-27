
if (!require("pacman")) install.packages("pacman", repos="http://cran.r-project.org")
pacman::p_load(tidyverse, multinomialTS, forecast, mgcv, future, furrr, ggtext, patchwork)

sunfish_pollen <- read_csv("./data/Sunfish_585_pollen_counts_bchron062722B_cleaned.csv")
sunfish_combined <- read_csv("./data/Sunfish_dataset.csv")
sunfish_ll <- read_csv("./data/Sunfish_LL_Recon_Feb22.csv")

sunfish_pollen_only <- sunfish_pollen[-c(1:8)]
# Remove spp with less than 5 count in whole core
filter_5_names <- names(which(colSums(sunfish_pollen_only) <= 20))

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
             `_2.5`, `_10`, `_90`, `_975`)) |>
  rename("depth" = "Depth.Correct",
          "age" = "Age.bchron062722B")

dim(sunfish_pollen_clean) 
dim(sunfish_combined) 

# Some minor inconsistencies in age-depth model likely due to repeatedly running Bchron without a seed
# Ages replace with lates version (consistent with the one Tony used for isotopes) 


# Distribute undifined pine counts between P.strobus and P.diplo based on their relative proportions
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

sunfish_grouped_long <- bind_rows(sunfish_other, sunfish_target_spp) |> 
  mutate(variablename = fct(variablename, levels = c("other", target_spp)))

# Plot spp
sunfish_grouped_long |>
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
  summarise(across(c(other, all_of(target_spp)), \(x) sum(x ,na.rm = TRUE)),
            age = mean(age, na.rm = TRUE)) |> 
              arrange(desc(age))
 

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
  select(bins, LakeElev.cm.below.mod) |> 
  group_by(bins) |> 
  summarise(mean_ll = mean(LakeElev.cm.below.mod)) |> 
  arrange(desc(bins))


sunfish_all <- sunfish_spp_binned |> 
  full_join(ll_mean) |> 
  arrange(desc(bins))

sunfish_all_interp <- sunfish_all |> 
  mutate(across(c(mean_ll), forecast::na.interp))



Y <- sunfish_all_interp |>
  select(other, all_of(target_spp)) |>
  as.matrix()

# Y[ ,1] <- rowSums(Y[ , c("other", "P.strobu", "Fagus", "Quercus", "Betula")])
# Y <- Y[ , c("other", "Tsuga")]

Tsample <- which(rowSums(Y) != 0)

# set up X with temp
X <- sunfish_all_interp |>
  select(mean_ll) |>
  as.matrix() |>
  scale()

matplot(X, type = 'l')

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
mnTS_mod <- mnTS(Y = Y[which(rowSums(Y) != 0),],
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
