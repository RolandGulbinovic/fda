# install.packages("fda")
# install.packages("refund")
# install.packages("tidyverse") 
#install.packages("refund")
library(fda)
library(dplyr)
library(tibble)
library(readxl)
library(refund)
library(funFEM)
df <- read_excel("galutinis_nenormalizuotas - Copy.xlsx")
consumption_data <- df[, 2:(ncol(df) - 2)]
n_clients <- ncol(consumption_data)
n_timepoints <- nrow(consumption_data)

# Ensure datetime column is POSIXct
datetime_col <- as.POSIXct(df$data_valanda)

start_time <- min(datetime_col)
hours_since_start <- as.numeric(difftime(datetime_col, start_time, units = "hours"))

# Normalize to [0, 1]
time_grid <- hours_since_start / max(hours_since_start)

# === 3. Create basis ===
n_basis <- 30
basis <- create.fourier.basis(rangeval = c(0, 1), nbasis = n_basis)

# === 4. Smooth each client time series ===
consumption_matrix <- as.matrix(consumption_data)
consumption_matrix <- as.matrix(sapply(consumption_data, as.numeric))
lambda_value <- 1e-4  # or even 0
fd_list <- lapply(1:ncol(consumption_matrix), function(i) {
  smooth.basis(argvals = time_grid, y = consumption_matrix[, i],
               fdParobj = fdPar(basis, Lfdobj = int2Lfd(2), lambda = lambda_value))$fd
})

coef_matrix <- do.call(cbind, lapply(fd_list, function(fdobj) fdobj$coefs))

common_basis <- fd_list[[1]]$basis

# Create the combined fd object
fd_smooth <- fd(coef_matrix, common_basis)


###########
# PCA
###########
fpca_result <- pca.fd(fd_smooth, nharm = 3)

# View variance explained
fpca_result$varprop  # proportion of variance for each component

# Scree plot
plot(fpca_result$values, type = "b", main = "Scree Plot", xlab = "Component", ylab = "Eigenvalue")

# Get variance proportions (already in decimal form, e.g., 0.43)
var_props <- round(100 * fpca_result$varprop[1:3], 2)

# Plot with variance in titles
par(mfrow = c(1, 3))  # 3 plots side by side

plot(fpca_result$harmonics[1],
     main = paste0("PC 1 (", var_props[1], "% variance)"))

plot(fpca_result$harmonics[2],
     main = paste0("PC 2 (", var_props[2], "% variance)"))

plot(fpca_result$harmonics[3],
     main = paste0("PC 3 (", var_props[3], "% variance)"))



##################
######### Regression - Attempt 1
##################
Y_mat <- eval.fd(time_grid, fd_smooth) %>% t()  # each row is a client


# Get avg price and temp.
# avg_price <- mean(df$`price EUR/MWH`, na.rm = TRUE)
# avg_temp  <- mean(df$`temperature C`, na.rm = TRUE)
# 
# 
# price_scalar <- rep(avg_price, n_clients)
# temp_scalar  <- rep(avg_temp, n_clients)


###### Regression where the X is a scalar
# calcualte the total consumption for each client
total_consumption <- apply(eval.fd(time_grid, fd_smooth), 2, sum)

xfdlist <- list(
  const = rep(1, n_clients),
  total = total_consumption
)

betalist <- list(
  fdPar(basis, int2Lfd(2), lambda = 1e-4),
  fdPar(basis, int2Lfd(2), lambda = 1e-4)
)
fos_model <- fRegress(fd_smooth, xfdlist, betalist)

par(mfrow = c(1, 2))

plot(fos_model$betaestlist[[1]], 
     main = "Intercept β₀(t)",
     xlab = "Normalized time", 
     ylab = "Energy consumption")

# Total consumption effect β₁(t)
plot(fos_model$betaestlist[[2]], 
     main = "Total Consumption Effect β₁(t)",
     xlab = "Normalized time", 
     ylab = "Effect size")


#########
# Regression where we use Principal Components as X's
scores <- fpca_result$scores

fpca1 <- scores[, 1]
fpca2 <- scores[, 2]
fpca3 <- scores[, 3]

xfdlist <- list(
  const = rep(1, n_clients),
  fpca1 = fpca1,
  fpca2 = fpca2,
  fpca3 = fpca3
)

betalist <- list(
  fdPar(basis, int2Lfd(2), lambda = 1e-4),
  fdPar(basis, int2Lfd(2), lambda = 1e-4),
  fdPar(basis, int2Lfd(2), lambda = 1e-4),
  fdPar(basis, int2Lfd(2), lambda = 1e-4)
)

fpca_model <- fRegress(fd_smooth, xfdlist, betalist)

par(mfrow = c(2, 2))

plot(fpca_model$betaestlist[[1]], main = "Intercept β₀(t)")
plot(fpca_model$betaestlist[[2]], main = "FPCA1 Effect β₁(t)")
plot(fpca_model$betaestlist[[3]], main = "FPCA2 Effect β₂(t)")
plot(fpca_model$betaestlist[[4]], main = "FPCA3 Effect β₃(t)")


par(mfrow = c(1, 1))


####################
####### Clustering
####################
# K is 2 because there are only two possible clusters for our data

res_v = funFEM(fd_smooth, K = 2, model = "AkjBk", init = "kmeans", lambda = 0, disp = TRUE)

# plotting
fdmeans_v = fd_smooth
fdmeans_v$coefs = t(res_v$prms$my)
plot(fdmeans_v, lwd = 2, lty = 1, col = 1:2, 
     xaxt = "n",  # suppress default x-axis
     xlab = "Date", ylab = "Mean Consumption",
     main = "Cluster Mean Curves (funFEM)")
datetime_col <- as.POSIXct(df$data_valanda)

tick_idx <- seq(1, length(datetime_col), length.out = 12)
tick_labels <- format(datetime_col[tick_idx], "%b-%d")
tick_positions <- time_grid[tick_idx]  # since x-axis is normalized [0,1]

axis(1, at = tick_positions, labels = tick_labels, las = 2)

legend("topleft", legend = paste("Cluster", 1:2), col = 1:2, lty = 1, lwd = 2)

# we can see that there are 2 clusters, one for high consumption, one for low

table(res_v$cls)
# most of the clients are in Low Cluster.


# Evaluate consumption curves
consumption_mat <- eval.fd(time_grid, fd_smooth)
total_consumption <- colSums(consumption_mat)

# Aggregate by cluster
aggregate(total_consumption, by = list(cluster = res_v$cls), FUN = summary)

##################
#### ANOVA
##################

library(fda.usc)


# Pairwise FANOVA - we check if the clustered groups really vary by a lot

# Took this function from Teams FDA channel
fANOVA.pointwise <- function(data, groups, t.seq, alpha = 0.05) {
  if (length(groups) != ncol(data)) stop("Length of groups must match number of columns in data")
  if (!is.factor(groups)) groups <- factor(groups)
  library(dplyr)
  n_time <- nrow(data)
  n_groups <- length(levels(groups))
  group_levels <- levels(groups)
  pvals <- numeric(n_time)
  mean.p <- matrix(NA, ncol = n_groups, nrow = n_time)
  # Determine number of pairwise comparisons
  combs <- combn(group_levels, 2)
  perm <- ncol(combs)
  Tukey.posthoc <- matrix(NA, ncol = perm, nrow = n_time)
  colnames(Tukey.posthoc) <- apply(combs, 2, paste, collapse = " - ")
  # Main loop: pointwise ANOVA and Tukey HSD
  for (i in 1:n_time) {
    dt <- data.frame(values = data[i, ], groups = groups)
    av <- aov(values ~ groups, data = dt)
    pvals[i] <- summary(av)[[1]]["Pr(>F)"][1,1]
    mean.p[i, ] <- dt %>%
      group_by(groups) %>%
      summarise(mean_val = mean(values), .groups = "drop") %>%
      pull(mean_val)
    tukey_res <- TukeyHSD(av)$groups[, "p adj"]
    Tukey.posthoc[i, ] <- tukey_res
  }
  overall_mean <- rowMeans(data)
  ## --- Plotting --- ##
  opar <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  # Plot ANOVA p-values
  plot(t.seq, pvals, type = "l", lwd = 2, col = "darkred",
       main = "Pointwise ANOVA p-values", xlab = "Time", ylab = "p-value", ylim = c(0, 1))
  abline(h = alpha, col = "blue", lty = 2, lwd = 2)
  # Plot group means
  col_set <- rainbow(n_groups)
  ylim_range <- range(c(mean.p, overall_mean)) * 1.05
  plot(t.seq, overall_mean, type = "l", lwd = 2, col = "black",
       main = "Group Mean Functions", xlab = "Time", ylab = "Mean", ylim = ylim_range)
  for (j in 1:n_groups) {
    lines(t.seq, mean.p[, j], col = col_set[j], lty = j + 1, lwd = 1.5)
  }
  legend("topright", legend = c("Overall", group_levels),
         col = c("black", col_set), lty = c(1, 2:(n_groups + 1)),
         lwd = c(2, rep(1.5, n_groups)))
  par(opar)
  ## --- Post-hoc Tukey plots --- ##
  opar2 <- par(mfrow = c(1, 1), ask = TRUE)
  for (i in 1:perm) {
    plot(t.seq, Tukey.posthoc[, i], type = "l", col = "purple", lwd = 2,
         main = paste("Tukey HSD p-values for", colnames(Tukey.posthoc)[i]),
         xlab = "Time", ylab = "p-value", ylim = c(0, 1))
    abline(h = alpha, col = "blue", lty = 2, lwd = 2)
  }
  par(opar2)
  ## --- Return --- ##
  return(list(
    p.values = pvals,
    TukeyHSD = Tukey.posthoc,
    group.means = mean.p,
    overall.mean = overall_mean,
    comparisons = colnames(Tukey.posthoc)
  )
  )
}

consumption_matrix <- eval.fd(time_grid, fd_smooth)  # [n_timepoints x n_clients]

groups <- as.factor(res_v$cls)

fanova_results <- fANOVA.pointwise(
  data = consumption_matrix,
  groups = groups,
  t.seq = time_grid,
  alpha = 0.05
)

# this gives plots and shows us that at ALL timepoints, the low and high clusters vary



############
### Regression - Attempt 2
###########

# Evaluate functional objects
cons_matrix <- eval.fd(time_grid, fd_smooth)
price_vec <- eval.fd(time_grid, price_fd)
temp_vec  <- eval.fd(time_grid, temp_fd)

# Average consumption across clients (scalar response per timepoint)
avg_consumption <- rowMeans(cons_matrix)

# Fit a time-level scalar regression
reg_df <- data.frame(cons = avg_consumption, price = price_vec, temp = temp_vec)

model <- lm(cons ~ price + temp, data = reg_df)
summary(model)


# This does get us a model and I think it still counts as FDA, because we are 
# using functional objects that are smoothened. We are just aligning it on time and not client-level
# (because our regressor vary over time and not across clients)

# from the result we can see that both price and temp are significant.
# Our model also explains 75% of variance.
