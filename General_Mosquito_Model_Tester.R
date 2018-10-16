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
lines(raw_rainfall_old_props, type = "l", col = "red")

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

model_parameters <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96,
                      "muE0" = 0.035, "muL0" = 0.035, "muP" = 0.211, "muM" = 0.096,
                      "lambda" = 13, "beta" = 21.19,
                      "overdisp" = 1.4788, "pop_frac" = 0.020, "z" = 200,
                      "tau_rain" = 50, "scaling_factor_rainfall" = 0,
                      "K_Max_Hill_Rainfall" = 0, "Hill_Rainfall_1" = 0, "Hill_Rainfall_2" = 0,
                      "tau_static" = 50, "scaling_factor_static" = 2000, "K_Max_Static" = 50,
                      "Washout_Threshold" = 4, "Washout_Exp_Decline" = 0.004, "Washout_Hill_1" = 400, "Washout_Hill_2" = 3,
                      "E" = 0, "L" = 0, "P" = 0, "M" = 0, "offset" = 0)

# Note that when using start and endtimes != 0 that you need to use actual initial conditions
# from the previous timepoint. Starting with 0s again will give you diff values and won't be comparable.
start_time <- 1500
end_time <- 2500

tic()
general_model_output <- general_mosquito_population_model(start_time, end_time, model_parameters, static_parameters, rainfall,
                                                          "linear", # Density Dependence Functional Form
                                                          "mean", # K_Rain Functional Form
                                                          "raw", # K_Rain Using Raw Average or Hill Function
                                                          "exponential") # K_Static Decline Functional Form
toc()

plot(general_model_output$M_Output, type = "l", lwd = 2, col = "blue")
lines(seq(1501, 2500), general_model_output$M_Output, type = "l", lwd = 2, col = "orange")

general_model_output$E_Output[500]
general_model_output$L_Output[500]
general_model_output$P_Output[500]
general_model_output$M_Output[500]

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "blue")
lines(general_model_output$K_Total, type = "l", lwd = 2, col = "red")



plot(general_model_output$`prior K`, type = "l", col = "red")

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "red")

length(seq(start_time, end_time-1, 1))
lines(seq(start_time, end_time-1, 1), general_model_output$M_Output, type = "l", lwd = 2, col = "blue")
plot(seq(start_time, end_time-1, 1), general_model_output$M_Output, type = "l", lwd = 2, col = "black")

plot(general_model_output$K_Rain, type = "l", lwd = 2, col = "red")
plot(general_model_output$K_Static, type = "l", lwd = 2, col = "red")

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "red", ylim = c(0, max(general_model_output$K_Total)))
par(new = T)
plot(rainfall_old, type = "l", col = blue, lwd = 2.5, axes = F, ylab = "", xlab = "")

plot(general_model_output$rainfallaverage_Kstatic, type = "l")

plot(general_model_output$`prior K`, type = "l")
par(new = T)
plot(rainfall[1:1200], type = "l", col = blue, lwd = 2.5, axes = F, ylab = "", xlab = "")

# Running and Plotting the Annularis Model
static_parameters <- c("dt" = 0.2, "dd_pow" = 1)
model_parameters_ann <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96, "muE0" = 0.035, "muL0" = 0.035,
                           "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 50, "beta" = 21.19,
                           "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 200, "z" = 200, "K_static" = 10,
                           "K_max" = 0, "hill_1" = 70, "hill_2" = -1, "E" = 0, "L" = 0, "P" = 0, "M" = 0)

annularis <- mosquito_population_model_ann(0, 3500, model_parameters_ann, static_parameters,
                                              rainfall_old, "linear", "mean", "raw")

plot(general_model_output$M_Output, type = "l", lwd = 2, col = "black")
lines(annularis$M_Output, type = "l", lwd = 2, col = "red")

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "black", ylim = c(0, 6000))
lines(annularis$K_total, type = "l", lwd = 2, col = "red")
par(new = T)
plot(rainfall_old, type = "l", col = blue, lwd = 2.5, ylab = "", axes = F)

# Running and Plotting the Culicifacies Model
model_parameters_cul <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96, "muE0" = 0.035, "muL0" = 0.035,
                          "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 30, "beta" = 21.19,
                          "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 200, "z" = 200, "K_static" = 0,
                          "K_max" = 2000, "hill_1" = 5, "hill_2" = -5, "E" = 0, "L" = 0, "P" = 0, "M" = 0)
culicifacies <- mosquito_population_model_cul(0, 3500, model_parameters_cul, static_parameters, rainfall_old, "linear", "exponential", "hill")

plot(general_model_output$M_Output, type = "l", lwd = 2, col = "black")
lines(culicifacies$M_Output, type = "l", lwd = 2, col = "red")

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "black")
lines(culicifacies$K_total, type = "l", lwd = 2, col = "red", xlab = "", ylab = "K")
par(new = T)
plot(rainfall_old, type = "l", col = blue, lwd = 2.5, ylab = "", axes = F)

# Running and Plotting the Fluviatilis Model
model_parameters_fluv <- c("dE" = 6.96, "dL" = 6.0185, "dP" = 1.96, "muE0" = 0.035, "muL0" = 0.035,
                           "muP" = 0.211, "muM" = 0.096, "lambda" = 13, "tau" = 50, "beta" = 21.19,
                           "overdisp" = 1.4788, "pop_frac" = 0.020, "scaling_factor" = 200, "z" = 200, "K_static" = 200,
                           "K_max" = 0, "hill_1" = 70, "hill_2" = -1, "E" = 0, "L" = 0, "P" = 0, "M" = 0,
                           "Washout_threshold" = 4, "Washout_Scaling_Factor" = 0.004,
                           "washout_hill_one" = 400, "washout_hill_two" = 3)
fluviatilis <- mosquito_population_model_fluv(0, 3500, model_parameters_fluv, static_parameters, rainfall_old,
                                              "linear", "exponential", "raw")

plot(general_model_output$M_Output, type = "l", lwd = 2, col = "black")
lines(fluviatilis$M_Output, type = "l", lwd = 2, col = "red")

plot(general_model_output$K_Total, type = "l", lwd = 2, col = "black")
lines(fluviatilis$K_total, type = "l", lwd = 2, col = "red", xlab = "", ylab = "K")
par(new = T)
plot(rainfall_old, type = "l", col = blue, lwd = 2.5, ylab = "", axes = F)


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

