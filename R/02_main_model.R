# ---------------------------------------------------------------
# 02_main_model.R
# Short-term effect of PM10 on daily mortality
# Quasi-Poisson time-series regression with natural splines
# ---------------------------------------------------------------

library(dlnm)
library(splines)

dir.create("results", showWarnings = FALSE)

data <- chicagoNMMAPS
n_years <- length(unique(data$year))

# ---- Exposure: average of lag 0 and lag 1 (lag01) -------------
lag_vec <- function(x, k) c(rep(NA, k), head(x, -k))

data$pm10_lag01 <- rowMeans(cbind(data$pm10, lag_vec(data$pm10, 1)))
data$temp_lag13 <- rowMeans(cbind(lag_vec(data$temp, 1),
                                  lag_vec(data$temp, 2),
                                  lag_vec(data$temp, 3)))

# ---- Model -----------------------------------------------------
# Confounder control:
#  - long-term trend and seasonality: ns(time) with 7 df per year
#  - same-day temperature and lagged (1-3 days) temperature: ns, 6 and 3 df
#  - day of the week
model <- glm(
  death ~ pm10_lag01 +
    ns(time, df = 7 * n_years) +
    ns(temp, df = 6) +
    ns(temp_lag13, df = 3) +
    dow,
  family = quasipoisson,
  data   = data,
  na.action = na.exclude
)

# ---- Effect per 10 ug/m3 increase ------------------------------
rr_per10 <- function(fit, term, increment = 10) {
  b  <- coef(fit)[term]
  se <- sqrt(vcov(fit)[term, term])
  pc <- 100 * (exp(increment * (b + c(0, -1.96, 1.96) * se)) - 1)
  data.frame(estimate = pc[1], lower = pc[2], upper = pc[3])
}

main_res <- cbind(
  analysis = "Main model (lag01, 7 df/year)",
  rr_per10(model, "pm10_lag01"),
  n_days = sum(!is.na(residuals(model))),
  dispersion = summary(model)$dispersion
)

print(round(main_res[, -1], 3))
cat(sprintf(
  "\nA 10 ug/m3 increase in PM10 (lag01) is associated with a %.2f%% (95%% CI %.2f to %.2f) change in daily mortality.\n",
  main_res$estimate, main_res$lower, main_res$upper))

write.csv(main_res, "results/main_model.csv", row.names = FALSE)

# ---- Residual diagnostics -------------------------------------
png("figures/03_residual_diagnostics.png", width = 1350, height = 500, res = 150)
par(mfrow = c(1, 2))
res <- residuals(model, type = "deviance")
plot(data$date, res, pch = ".", xlab = NULL, ylab = "Deviance residuals",
     main = "Residuals over time")
abline(h = 0, col = "red")
acf(na.omit(res), main = "Autocorrelation of residuals", lag.max = 30)
dev.off()

saveRDS(model, "results/main_model.rds")
