# Loading the Indian Vector Modelling Library
library(IndianVectorModelling)
library(tictoc)

# Plotting Parameters
blue <- rgb(0, 0, 1, alpha=0.25)
red <- rgb(1, 0, 0, alpha=0.35)

# Loading Raw Rainfall Data, Checking It's Alright
raw_rainfall_old <- read.csv("C:/Users/Charlie Whittaker/Documents/Fluv Tester/precipitation-(chirps).csv")
raw_rainfall_full_length <- read.csv("C:/Users/Charlie Whittaker/Downloads/precipitation-(chirps) (1) (longer).csv")
raw_rainfall_old_props <- c(rep(0, 242), raw_rainfall_old$X23.3534N..85.3943W)
raw_rainfall_full_length_props <- raw_rainfall_full_length$X23.3525N..85.3935W
plot(raw_rainfall_full_length_props, type = "l")
lines(raw_rainfall_old_props, type = "l", col = "blue")

# Processing Old Raw Rainfall Data and Plotting It
colnames(raw_rainfall_old) = c("Date", "Rainfall")
raw_rainfall_old <- raw_rainfall_old[, 2]
rainfall_old <- c()
for(i in 1:length(raw_rainfall_old)) {
  temp_rainfall <- rep(raw_rainfall_old[i], 5)
  rainfall_old <- c(rainfall_old, temp_rainfall)
}
plot(rainfall_old, type = "l", col = blue, lwd = 2.5, ylab = "Rainfall (mm)")

# Processing Proper Raw Rainfall Data and Plotting It
colnames(raw_rainfall_full_length) = c("Date", "Rainfall")
raw_rainfall <- raw_rainfall_full_length[, 2]
rainfall <- c()
for(i in 1:length(raw_rainfall)) {
  temp_rainfall <- rep(raw_rainfall[i], 5)
  rainfall <- c(rainfall, temp_rainfall)
}
plot(rainfall, type = "l", col = blue, lwd = 2.5, ylab = "Rainfall (mm)")

# Observed Mosquito Densities for Plotting
fluv <- c(4, 12, 4, 0, 0, 240, 307, 464, 304, 82, 64, 8, 8, 4, 2, 0, 12, 155, 295, 356, 243, 94, 51)
cul <- c(148, 310, 213, 206, 241, 7, 5, 5, 5, 19, 36, 69, 126, 214, 238, 218, 186, 19, 4, 5, 5, 13, 33)
ann <- c(170, 152, 238, 209, 242, 117, 104, 114, 90, 121, 103, 102, 53, 144, 200, 173, 170, 124, 116, 120,
               98, 123, 72)

# Running the General Mosquito Model
timepoints <- seq(150, 3450, 150)
static_parameters <- c("dt" = 0.2, "dd_pow" = 1)
model_parameters <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96, "muE0" = 0.035, "muL0" = 0.035, "muP" = 0.211,
                      "muM" = 0.096, "lambda" = 13, "tau_rain" = 30, "beta" = 21.19, "overdisp" = 1.4788, "pop_frac" = 0.020,
                      "scaling_factor_rainfall" = 0,
                      "z" = 200,
                      "K_static" = 200, "K_Max_Hill_Rainfall" = 0, "Hill_Rainfall_1" = 0, "Hill_Rainfall_2" = 0,
                      "E" = 0, "L" = 0, "P" = 0, "M" = 0,
                      "Washout_Threshold" = 4, "Washout_Exp_Decline" = 0.004,
                      "Washout_Hill_1" = 0, "Washout_Hill_2" = 0, "tau_static" = 40,
                      "offset" = 1200,
                      "scaling_factor_static" = 200)

tic()
general_model_output <- general_mosquito_population_model(0, 3500, model_parameters, static_parameters, rainfall,
                                                          "linear", # Density Dependence Functional Form
                                                          "mean", # K_Rain Functional Form
                                                          "exponential", # K_Static Decline Functional Form
                                                          "raw") # K_Rain Using Raw Average or Hill Function
toc()

plot(general_model_output$M_Output, type = "l", lwd = 2)
plot(general_model_output$K_Rain, type = "l", lwd = 2, col = "red")
plot(general_model_output$K_Static, type = "l", lwd = 2, col = "red")
plot(general_model_output$K_Total, type = "l", lwd = 2, col = "red")
par(new = T)
plot(rainfall_old, type = "l", col = blue, lwd = 2.5, axes = F, ylab = "", xlab = "")

plot(general_model_output$`prior K`, type = "l")
par(new = T)
plot(rainfall[1:1200], type = "l", col = blue, lwd = 2.5, axes = F, ylab = "", xlab = "")
kloop <- general_model_output$`prior K`

plot(general_model_output$rainfallaverage_Kstatic, type = "l")
lines(fluviatilis$rainfallaverage_Kstatic, type = "l", col = "red")


# Running and Plotting the Fluviatilis Model
model_parameters_fluv <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96, "muE0" = 0.035, "muL0" = 0.035,
                           "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 40, "beta" = 21.19,
                           "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 200, "z" = 200, "K_static" = 200,
                           "K_max" = 0, "hill_1" = 70, "hill_2" = -1, "E" = 0, "L" = 0, "P" = 0, "M" = 0,
                           "Washout_threshold" = 4, "Washout_Scaling_Factor" = 0.004,
                           "washout_hill_one" = 400, "washout_hill_two" = 3)
fluviatilis <- mosquito_population_model_fluv(0, 3500, model_parameters_fluv, static_parameters, rainfall_old, "linear", "mean", "exponential")

plot(general_model_output$M_Output, type = "l", lwd = 2, col = "black")
lines(fluviatilis$M_Output, type = "l", lwd = 2, col = "red")

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "black")
lines(fluviatilis$K_total, type = "l", lwd = 2, col = "red", xlab = "", ylab = "K")
par(new = T)
plot(rainfall_old, type = "l", col = blue, lwd = 2.5, ylab = "", axes = F)

plot(fluviatilis$M_Output, type = "l", lwd = 2, col = "red", axes = F)


# Running and Plotting the Culicifacies Model
model_parameters_cul <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96, "muE0" = 0.035, "muL0" = 0.035,
                          "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 30, "beta" = 21.19,
                          "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 200, "z" = 200, "K_static" = 0,
                          "K_max" = 0, "hill_1" = 0, "hill_2" = 0, "E" = 0, "L" = 0, "P" = 0, "M" = 0)
culicifacies <- mosquito_population_model_cul(0, 3500, model_parameters_cul, static_parameters, rainfall_old, "linear", "mean", "raw")

plot(general_model_output$M_Output, type = "l", lwd = 2, col = "black")
lines(culicifacies$M_Output, type = "l", lwd = 2, col = "red")

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "black")
lines(culicifacies$K_total, type = "l", lwd = 2, col = "red", xlab = "", ylab = "K")



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






################ Miscellaneous Checking Stuff #################


### Testing the Culicifacies Modality for a Mean Functional Form of K_Rain ###
# Calculating Where the Rainfall Average Should Be Calculated From
start_point <- 600 - (model_parameters["tau_rain"]/static_parameters["dt"]) + 1
end_point <- model_parameters["offset"]
# Calculating By Hand the
rainfall_snippet <- model_parameters["scaling_factor_rainfall"] * mean(rainfall[start_point : end_point])
rainfall_snippet1 <- model_parameters["scaling_factor_rainfall"] * mean(rainfall[(start_point + 1): (end_point + 1)])
rainfall_snippet2 <- model_parameters["scaling_factor_rainfall"] * mean(rainfall[(start_point + 2): (end_point + 2)])
rainfall_snippet3 <- model_parameters["scaling_factor_rainfall"] * mean(rainfall[(start_point + 3): (end_point + 3)])
# Comparing Out of Model Calculations to Those Done In the Model
output$K_Rain[1:4]
rainfall_snippet; rainfall_snippet1; rainfall_snippet2; rainfall_snippet3;

