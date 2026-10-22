# ---------------------------------------------------------------
# 03_sensitivity.R
# Sensitivity analyses for the PM10-mortality association:
#   (a) single-day lags 0-5
#   (b) degrees of freedom per year for the time spline
#   (c) cause-specific mortality (cardiovascular, respiratory)
# ---------------------------------------------------------------

library(dlnm)
library(splines)
library(ggplot2)

data <- chicagoNMMAPS
n_years <- length(unique(data$year))

lag_vec <- function(x, k) if (k == 0) x else c(rep(NA, k), head(x, -k))

data$pm10_lag01 <- rowMeans(cbind(data$pm10, lag_vec(data$pm10, 1)))
data$temp_lag13 <- rowMeans(cbind(lag_vec(data$temp, 1),
                                  lag_vec(data$temp, 2),
                                  lag_vec(data$temp, 3)))

fit_pm <- function(exposure, outcome = "death", df_year = 7) {
  data$exposure <- exposure
  f <- as.formula(paste0(
    outcome, " ~ exposure + ns(time, df = ", df_year * n_years, ")",
    " + ns(temp, df = 6) + ns(temp_lag13, df = 3) + dow"))
  m  <- glm(f, family = quasipoisson, data = data, na.action = na.exclude)
  b  <- coef(m)["exposure"]
  se <- sqrt(vcov(m)["exposure", "exposure"])
  pc <- 100 * (exp(10 * (b + c(0, -1.96, 1.96) * se)) - 1)
  data.frame(estimate = pc[1], lower = pc[2], upper = pc[3])
}

# (a) Single-day lags
lag_res <- do.call(rbind, lapply(0:5, function(k)
  cbind(analysis = "Single-day lag", setting = paste0("lag ", k),
        fit_pm(lag_vec(data$pm10, k)))))

# (b) Degrees of freedom per year for seasonality control
df_res <- do.call(rbind, lapply(c(4, 6, 7, 8, 10, 12), function(d)
  cbind(analysis = "Time spline df/year", setting = paste0(d, " df/year"),
        fit_pm(data$pm10_lag01, df_year = d))))

# (c) Cause-specific mortality
cause_res <- do.call(rbind, lapply(c("death", "cvd", "resp"), function(o)
  cbind(analysis = "Cause of death",
        setting = c(death = "All-cause", cvd = "Cardiovascular", resp = "Respiratory")[[o]],
        fit_pm(data$pm10_lag01, outcome = o))))

sens <- rbind(lag_res, df_res, cause_res)
rownames(sens) <- NULL
sens[, 3:5] <- round(sens[, 3:5], 3)
print(sens)
write.csv(sens, "results/sensitivity_analyses.csv", row.names = FALSE)

# ---- Forest plot ----------------------------------------------
sens$setting <- factor(sens$setting, levels = rev(unique(sens$setting)))

p <- ggplot(sens, aes(estimate, setting)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, colour = "#12505e") +
  geom_point(size = 2, colour = "#12505e") +
  facet_grid(analysis ~ ., scales = "free_y", space = "free_y") +
  labs(x = "% change in daily mortality per 10 ug/m3 PM10 (95% CI)", y = NULL,
       title = "Sensitivity analyses") +
  theme_minimal(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))

ggsave("figures/04_sensitivity_forest.png", p, width = 9, height = 6, dpi = 150)
