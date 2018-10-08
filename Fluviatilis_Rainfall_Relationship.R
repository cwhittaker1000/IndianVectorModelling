raw_rainfall <- read.csv("C:/Users/Charlie Whittaker/Documents/Fluv Tester/precipitation-(chirps).csv")
colnames(raw_rainfall) = c("Date", "Rainfall")
raw_rainfall <- raw_rainfall[, 2]
fluv <- c(4, 12, 4, 0, 0, 240, 307, 464, 304, 82, 64, 8, 8, 4, 2, 0, 12, 155, 295, 356, 243, 94, 51)
timepoints <- seq(150, 3450, 150)

rainfall <- c()
for(i in 1:length(raw_rainfall)) {
  temp_rainfall <- rep(raw_rainfall[i], 5)
  rainfall <- c(rainfall, temp_rainfall)
}

plot(rainfall, type = "l", col = "red", lwd = 2)

# Average Rainfall and Washout Threshold - Exponential Decline Post-Washout Cessation
# Play around with marker resetting everytime we have washout vs the first time we lose washout. Probably the former.
# That's what's included for now.
# NEED TO SORT THE FIRST FEW DAYS WHERE WE DON'T HAVE PREVIOUS RAINFALL!!!
tau <- 200
K <- c()
Washout_Threshold <- 5
Kmax <- 20
marker <- 0
b <- 500 # sets the midpoint of the decline
a <- 10

for (i in 1:length(rainfall)) {

  if (i <= tau) {
    rainfall_selection <- rainfall[1:i]
    rain_average <- (1/i) * sum(rainfall_selection)
    if (rain_average > Washout_Threshold) {
      K[i] <- 0
      marker <- 0
    }
    else {
      K[i] <- Kmax
    }
  }

  else {

    # Selecting the relevant previous day's rainfall and calculating the average
    rainfall_selection <- rainfall[(i-tau) : i]
    rain_average <- (1/tau) * sum(rainfall_selection)

    # If washout occurring, set K to 0 and keep marker at 0
    if (rain_average > Washout_Threshold) {
      K[i] <- 0
      marker <- 0
    }

    # When washout stops: if it's the first timepoint since washout,
    # marker will still be 0, K set to K max, and washout stop time will be stored as marker_stop.
    else {
      if (marker == 0) {
        K[i] <- Kmax
        marker_stop <- i
        print(marker_stop)
        marker <- 1
      }

      # Exponential decline from Kmax as time progresses
      else {
        K[i] <- Kmax / (1 + (((i - marker_stop) / b) ^ a))
      }
    }
  }
}

plot(seq(0, 3500, 1), rep(5, 3501), type = "l", col = "red", ylim = c(0, 20), lwd = 3)
par(new = T)
plot(average, type = "l", col = "black", lwd = 1, ylim = c(0, 20), xlab = "", ylab = "")
par(new = T)
plot(K, type = "l", ylim = c(0, 20), xlab = "", ylab = "")
par(new = T)
plot(timepoints, fluv, xlim = c(0, 3500), axes = F, ylab = "", xlab = "", pch = 20, col = "dark grey", cex = 3)







### Miscellaneous Code ###
# Average Rainfall
average <- c()
for (i in 1:length(rainfall)) {
  if (i <= tau) {
    rainfall_selection <- rainfall[1:i]
    average[i] <- (1/i) * sum(rainfall_selection)
  }
  else {
    rainfall_selection <- rainfall[(i-tau) : i]
    average[i] <- (1/tau) * sum(rainfall_selection)
  }
}
plot(timepoints, fluv, xlim = c(0, 3500), axes = F, ylab = "", xlab = "")
par(new = T)
plot(average, type = "l")

# Average Rainfall and Washout Threshold - No Decline Post-Washout Cessation
for (i in 1:length(rainfall)) {
  if (i <= tau) {
    rainfall_selection <- rainfall[1:i]
    rain_average <- (1/i) * sum(rainfall_selection)
    if (rain_average > Washout_Threshold) {
      K[i] <- 0
    }
    else {
      K[i] <- Kmax
    }
  }
  else {
    rainfall_selection <- rainfall[(i-tau) : i]
    rain_average <- (1/tau) * sum(rainfall_selection)
    if (rain_average > Washout_Threshold) {
      K[i] <- 0
    }
    else {
      K[i] <- Kmax
    }
  }
}

# Average Rainfall and Washout Threshold - Exponential Decline Post-Washout Cessation
# Play around with marker resetting everytime we have washout vs the first time we lose washout. Probably the former.
# That's what's included for now.
# NEED TO SORT THE FIRST FEW DAYS WHERE WE DON'T HAVE PREVIOUS RAINFALL!!!
tau <- 200
K <- c()
Washout_Threshold <- 5
Kmax <- 20
sf <- 0.002
marker <- 0

for (i in 1:length(rainfall)) {

  if (i <= tau) {
    rainfall_selection <- rainfall[1:i]
    rain_average <- (1/i) * sum(rainfall_selection)
    if (rain_average > Washout_Threshold) {
      K[i] <- 0
      marker <- 0
    }
    else {
      K[i] <- Kmax
    }
  }

  else {

    # Selecting the relevant previous day's rainfall and calculating the average
    rainfall_selection <- rainfall[(i-tau) : i]
    rain_average <- (1/tau) * sum(rainfall_selection)

    # If washout occurring, set K to 0 and keep marker at 0
    if (rain_average > Washout_Threshold) {
      K[i] <- 0
      marker <- 0
    }

    # When washout stops: if it's the first timepoint since washout,
    # marker will still be 0, K set to K max, and washout stop time will be stored as marker_stop.
    else {
      if (marker == 0) {
        K[i] <- Kmax
        marker_stop <- i
        print(marker_stop)
        marker <- 1
      }

      # Exponential decline from Kmax as time progresses
      else {
        K[i] <- Kmax * exp(-sf * (i - marker_stop))
      }
    }
  }
}

plot(seq(0, 3500, 1), rep(5, 3501), type = "l", col = "red", ylim = c(0, 20), lwd = 3)
par(new = T)
plot(average, type = "l", col = "black", lwd = 1, ylim = c(0, 20), xlab = "", ylab = "")
par(new = T)
plot(K, type = "l", ylim = c(0, 20), xlab = "", ylab = "")
par(new = T)
plot(timepoints, fluv, xlim = c(0, 3500), axes = F, ylab = "", xlab = "", pch = 20, col = "dark grey", cex = 3)





