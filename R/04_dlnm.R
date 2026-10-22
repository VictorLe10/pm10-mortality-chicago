# ---------------------------------------------------------------
# 04_dlnm.R
# Distributed lag models:
#   (a) PM10: linear exposure-response, lag structure over 0-10 days
#   (b) Temperature: non-linear exposure-lag-response over 0-21 days
# ---------------------------------------------------------------

library(dlnm)
library(splines)

data <- chicagoNMMAPS
n_years <- length(unique(data$year))

# ---- Cross-basis definitions ----------------------------------
cb_pm <- crossbasis(data$pm10, lag = 10,
                    argvar = list(fun = "lin"),
                    arglag = list(fun = "ns", knots = logknots(10, 3)))

cb_temp <- crossbasis(data$temp, lag = 21,
                      argvar = list(fun = "ns",
                                    knots = quantile(data$temp, c(0.10, 0.75, 0.90), na.rm = TRUE)),
                      arglag = list(fun = "ns", knots = logknots(21, 3)))

model_dlnm <- glm(death ~ cb_pm + cb_temp + ns(time, df = 7 * n_years) + dow,
                  family = quasipoisson, data = data)

# ---- (a) PM10 lag-response ------------------------------------
pred_pm <- crosspred(cb_pm, model_dlnm, at = 10, cumul = TRUE)

pm_overall <- 100 * (c(pred_pm$allRRfit, pred_pm$allRRlow, pred_pm$allRRhigh) - 1)
cat(sprintf("Cumulative effect over lags 0-10 per 10 ug/m3 PM10: %.2f%% (95%% CI %.2f to %.2f)\n",
            pm_overall[1], pm_overall[2], pm_overall[3]))

png("figures/05_pm10_lag_response.png", width = 1200, height = 750, res = 150)
plot(pred_pm, var = 10, type = "p", pch = 19, cex = 1.1, ci = "bars", col = "#12505e",
     ylab = "RR per 10 ug/m3", xlab = "Lag (days)",
     main = "PM10: lag-response association")
dev.off()

# ---- (b) Temperature exposure-response ------------------------
# Centre on the minimum mortality temperature
# (searched between the 1st and 99th percentile, where estimates are stable)
p01 <- round(quantile(data$temp, 0.01, na.rm = TRUE) * 2) / 2
p99 <- round(quantile(data$temp, 0.99, na.rm = TRUE) * 2) / 2
pred_tmp0 <- crosspred(cb_temp, model_dlnm, at = seq(p01, p99, by = 0.5), cen = 10)
mmt <- pred_tmp0$predvar[which.min(pred_tmp0$allRRfit)]
pred_tmp <- crosspred(cb_temp, model_dlnm, cen = mmt, by = 0.5)
cat(sprintf("Minimum mortality temperature: %.1f C\n", mmt))

# Cumulative RR at the extremes, over the full lag period (0-21)
# and over the first days only (0-3), to separate immediate from delayed effects
extreme_full  <- crosspred(cb_temp, model_dlnm, cen = mmt, at = c(p01, p99))
extreme_short <- crosspred(cb_temp, model_dlnm, cen = mmt, at = c(p01, p99), lag = c(0, 3))
temp_res <- data.frame(
  percentile  = rep(c("1st (cold)", "99th (heat)"), 2),
  lag_period  = rep(c("0-21", "0-3"), each = 2),
  temperature = rep(c(p01, p99), 2),
  RR          = c(extreme_full$allRRfit,  extreme_short$allRRfit),
  lower       = c(extreme_full$allRRlow,  extreme_short$allRRlow),
  upper       = c(extreme_full$allRRhigh, extreme_short$allRRhigh),
  row.names   = NULL
)
temp_res[, 4:6] <- round(temp_res[, 4:6], 3)
print(temp_res)
write.csv(temp_res, "results/temperature_extremes.csv", row.names = FALSE)

png("figures/06_temperature_overall.png", width = 1200, height = 750, res = 150)
plot(pred_tmp, "overall", col = "#c0392b", lwd = 2,
     xlab = "Temperature (C)", ylab = "Cumulative RR (lag 0-21)",
     main = "Temperature-mortality association")
abline(v = mmt, lty = 2, col = "grey40")
dev.off()

png("figures/07_temperature_3d.png", width = 1100, height = 900, res = 150)
plot(pred_tmp, xlab = "Temperature (C)", ylab = "Lag (days)", zlab = "RR",
     theta = 210, phi = 30, lphi = 30, col = "#d6eaf0",
     main = "Temperature exposure-lag-response surface")
dev.off()

# Short (lag 0) vs delayed (lag 10) effects show heat acting fast and cold acting slowly
# (centred at 20 C, a mild reference temperature, so both tails are visible)
pred_slices <- crosspred(cb_temp, model_dlnm, cen = 20, by = 0.5)
png("figures/08_temperature_by_lag.png", width = 1300, height = 600, res = 150)
par(mfrow = c(1, 2))
plot(pred_slices, "slices", lag = 0, col = "#c0392b", lwd = 2,
     xlab = "Temperature (C)", ylab = "RR (ref. 20 C)", main = "Lag 0: same-day effect")
plot(pred_slices, "slices", lag = 10, col = "#12505e", lwd = 2,
     xlab = "Temperature (C)", ylab = "RR (ref. 20 C)", main = "Lag 10: delayed effect")
dev.off()
