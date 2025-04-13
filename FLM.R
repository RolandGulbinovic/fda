# install.packages("fda")
# install.packages("refund")
# install.packages("tidyverse") 

library(fda)
library(dplyr)
library(tibble)
library(readxl)

df <- read_excel("galutinis_nenormalizuotas - Copy.xlsx")
head(consumption_matrix[, 2], 10)
consumption_data <- df[, 2:(ncol(df) - 2)]
n_clients <- ncol(consumption_data)
n_timepoints <- nrow(consumption_data)

# === 2. Create normalized time grid ===
# Assuming you have a datetime column or rownames
# Ensure datetime column is POSIXct
datetime_col <- as.POSIXct(df$data_valanda)

# Calculate hours since the start
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
