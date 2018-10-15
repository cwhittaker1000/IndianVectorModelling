# Plotting Parameters
blue <- rgb(0, 0, 1, alpha=0.25)
red <- rgb(1, 0, 0, alpha=0.35)

# Loading Raw Rainfall Data and Processing It
raw_rainfall <- read.csv("C:/Users/Charlie Whittaker/Documents/Fluv Tester/precipitation-(chirps).csv")
raw_rainfall2 <- read.csv("C:/Users/Charlie Whittaker/Downloads/precipitation-(chirps) (1) (bloop).csv")
raw_rainfall_props <- c(rep(0, 120), raw_rainfall$X23.3534N..85.3943W)
raw_rainfall_2_props <- raw_rainfall2$X23.3525N..85.3935W

plot(raw_rainfall_2_props, type = "l")
lines(raw_rainfall_props, type = "l", col = "blue")


colnames(raw_rainfall) = c("Date", "Rainfall")
raw_rainfall <- raw_rainfall[, 2]
rainfall <- c()
for(i in 1:length(raw_rainfall)) {
  temp_rainfall <- rep(raw_rainfall[i], 5)
  rainfall <- c(rainfall, temp_rainfall)
}
plot(rainfall, type = "l", col = blue, lwd = 2.5, ylab = "Rainfall (mm)")
plot(0, 0, xlim = c(0, 3500), ylim = c(0, 70))
polygon(x = c(seq(1, 3505, 1), rev(seq(1, 3505, 1))),
        y = c(rainfall, rep(0, 3505)), col = blue, border = NA)

# Monthly Rainfall Data and Observed Mosquito Densities for Plotting
monthly_rainfall_data <- c(70.9633, 151.3844, 231.8042, 240.04, 212.767, 84.9041, 10, 10, 17.297, 10, 18.3768, 51.4621,
                           98.7449, 265.5596, 305.7442, 194.4492, 29.9026, 66.5369, 10, 10, 19.0457, 10, 10)
fluv <- c(4, 12, 4, 0, 0, 240, 307, 464, 304, 82, 64, 8, 8, 4, 2, 0, 12, 155, 295, 356, 243, 94, 51)
cul <- c(148, 310, 213, 206, 241, 7, 5, 5, 5, 19, 36, 69, 126, 214, 238, 218, 186, 19, 4, 5, 5, 13, 33)
ann <- c(170, 152, 238, 209, 242, 117, 104, 114, 90, 121, 103, 102, 53, 144, 200, 173, 170, 124, 116, 120,
               98, 123, 72)
timepoints <- seq(150, 3450, 150)
library(IndianVectorModelling)
devtools::document()
# Running and Plotting the Culicifacies Model
static_parameters <- c("dt" = 0.2, "dd_pow" = 1)
model_parameters_cul <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96, "muE0" = 0.035, "muL0" = 0.035,
                          "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 50, "beta" = 21.19,
                          "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 200, "z" = 200, "K_static" = 0,
                          "K_max" = 100, "hill_1" = 10, "hill_2" = -1, "E" = 0, "L" = 0, "P" = 0, "M" = 0)

culicifacies <- mosquito_population_model_cul(0, 3500, model_parameters_cul, static_parameters,
                                              rainfall, "linear", "mean", "raw")
plot(culicifacies$K_total, type = "l", lwd = 4, col = "black", xlab = "", ylab = "K")
par(new = TRUE)
plot(rainfall, type = "l", col = blue, lwd = 2.5, ylab = "", axes = F)

plot(culicifacies$M_Output, type = "l", lwd = 2, col = "red", axes = F)

barplot(monthly_rainfall_data, axes = F, col = blue, xlab = "", ylab = "", ylim = c(0, 350), border= 0)
axis(side=4, at = pretty(range(monthly_rainfall_data)))
plot(rainfall, type = "l", col = blue, lwd = 2, axes = F, xlim = c(0, 3500), xlab = "", ylab = "")
par(new = TRUE)
plot(culicifacies$M_Output, type = "l", ylab = "Mosquito Catch", col = "#64B2D1", lwd = 4)
par(new = TRUE)
plot(timepoints, cul, xlim = c(0, 3500), ylab = "", xlab = "", pch = 20, col = "#64B2D1", cex = 2, axes = F)

# Bit of a fudge to quickly plot relatively nicely K and Rainfall
timepoints <- seq(0, 3505/(30/0.2), (1/30 * 0.2))
timepoints <- timepoints[1:3500]
plot(0, 0, xlim = c(0, 3505), ylim = c(0, 6000), pch = 20, cex = 0.01, axes = F, xlab = "", ylab = "")
polygon(x = c(seq(1, 3505, 1), rev(seq(1, 3505, 1))),
        y = c(rainfall*100, rep(0, 3505)), col = "red", border = NA)
par(new = TRUE)
plot(timepoints, culicifacies$K_total, type = "l", lwd = 3, col = "black", xlab = "Time (Months)", ylab = "K")

# Running and Plotting the Fluviatilis Model
static_parameters <- c("dt" = 0.2, "dd_pow" = 1)
model_parameters_fluv <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 0.65, "muE0" = 0.035, "muL0" = 0.035,
                           "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 40, "beta" = 21.19,
                           "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 200, "z" = 200, "K_static" = 200,
                           "K_max" = 0, "hill_1" = 70, "hill_2" = -1, "E" = 0, "L" = 0, "P" = 0, "M" = 0,
                           "Washout_threshold" = 4, "Washout_Scaling_Factor" = 0.004,
                           "washout_hill_one" = 400, "washout_hill_two" = 3)

fluviatilis <- mosquito_population_model_fluv(0, 3500, model_parameters_fluv, static_parameters,
                                        rainfall, "linear", "mean", "exponential")
plot(fluviatilis$K_total, type = "l", lwd = 2, col = "black")
par(new = TRUE)
plot(rainfall, type = "l", col = blue, lwd = 2.5, ylab = "", axes = F)

plot(fluviatilis$M_Output, type = "l", lwd = 2, col = "red", axes = F)

barplot(monthly_rainfall_data, axes = F, col = blue, xlab = "", ylab = "", ylim = c(0, 350), border= 0)
axis(side=4, at = pretty(range(monthly_rainfall_data)))
plot(rainfall, type = "l", col = blue, lwd = 2, axes = F, xlim = c(0, 3500), xlab = "", ylab = "")
par(new = TRUE)
plot(fluviatilis$M_Output, type = "l", ylab = "Mosquito Catch", col = "dark grey", lwd = 4)
par(new = TRUE)
plot(timepoints, fluv, xlim = c(0, 3500), ylab = "", xlab = "", pch = 20, col = "dark grey", cex = 2, axes = F)

# Bit of a fudge to quickly plot relatively nicely K and Rainfall
timepoints <- seq(0, 3505/(30/0.2), (1/30 * 0.2))
timepoints <- timepoints[1:3500]
plot(0, 0, xlim = c(0, 3505), ylim = c(0, 6000), pch = 20, cex = 0.01, axes = F, xlab = "", ylab = "")
polygon(x = c(seq(1, 3505, 1), rev(seq(1, 3505, 1))),
        y = c(rainfall*100, rep(0, 3505)), col = "red", border = NA)
par(new = TRUE)
plot(timepoints, fluviatilis$K_total, type = "l", lwd = 3, col = "black", xlab = "Time (Months)", ylab = "K")

# Running and Plotting the Annularis Model
static_parameters <- c("dt" = 0.2, "dd_pow" = 1)
model_parameters_ann <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 0.65, "muE0" = 0.035, "muL0" = 0.035,
                           "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 50, "beta" = 21.19,
                           "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 210, "z" = 200, "K_static" = 13,
                           "K_max" = 0, "hill_1" = 70, "hill_2" = -1, "E" = 10000, "L" = 500, "P" = 500, "M" = 1000)

annularis <- mosquito_population_model_ann(0, 3500, model_parameters_ann, static_parameters,
                                              rainfall, "linear", "mean", "raw")
plot(annularis$K_total, type = "l", lwd = 2, col = "black", ylim = c(0, max(annularis$K_total)))
par(new = TRUE)
plot(rainfall, type = "l", col = blue, lwd = 2.5, ylab = "", axes = F)

plot(annularis$M_Output, type = "l", lwd = 2, col = "red", axes = F)

barplot(monthly_rainfall_data, axes = F, col = blue, xlab = "", ylab = "", ylim = c(0, 350), border= 0)
axis(side=4, at = pretty(range(monthly_rainfall_data)))
plot(rainfall, type = "l", col = blue, lwd = 2, axes = F, xlim = c(0, 3500), xlab = "", ylab = "")
par(new = TRUE)
plot(annularis$M_Output, type = "l", ylab = "Mosquito Catch", col = "dark green", lwd = 4, ylim = c(0, max(annularis$M_Output)))
par(new = TRUE)
plot(timepoints, ann, xlim = c(0, 3500), ylim = c(0, max(ann)), ylab = "", xlab = "", pch = 20, col = "dark green", cex = 2, axes = F)

plot(annularis$M_Output, type = "l", ylab = "Mosquito Catch", col = "dark grey", lwd = 2, ylim = c(0, max(annularis$M_Output)), axes = F)
axis(side=4)
par(new = TRUE)
plot(timepoints, ann, xlim = c(0, 3500), ylim = c(0, max(ann)), ylab = "", xlab = "", pch = 20, col = "dark green", cex = 2)

# Bit of a fudge to quickly plot relatively nicely K and Rainfall
timepoints <- seq(0, 3505/(30/0.2), (1/30 * 0.2))
timepoints <- timepoints[1:3500]
plot(0, 0, xlim = c(0, 3505), ylim = c(0, 6000), pch = 20, cex = 0.01, axes = F, xlab = "", ylab = "")
polygon(x = c(seq(1, 3505, 1), rev(seq(1, 3505, 1))),
        y = c(rainfall*100, rep(0, 3505)), col = "red", border = NA)
par(new = TRUE)
plot(timepoints, annularis$K_total, type = "l", lwd = 3, col = "black", xlab = "Time (Months)", ylab = "K",
     ylim = c(0, max(annularis$K_total)))




# Testing the C++ Hill Function Out
hill_output <- matrix(nrow = 500, ncol = 2)
for (i in 1:nrow(hill_output)) {
  hill_output[i, 1] <- i
  hill_output[i, 2] <- Hill_Function(i, 10, 200, 3)
}
plot(hill_output, type = "l", xlim = c(0, 500), ylim = c(0, max(hill_output[, 2])), col = "black")



