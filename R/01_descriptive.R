# ---------------------------------------------------------------
# 01_descriptive.R
# Descriptive analysis of daily mortality, PM10 and temperature
# Data: chicagoNMMAPS (dlnm package), Chicago 1987-2000
# ---------------------------------------------------------------

library(dlnm)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

data <- chicagoNMMAPS

# ---- Summary table --------------------------------------------
vars <- c(death = "All-cause deaths (daily)",
          cvd   = "Cardiovascular deaths (daily)",
          resp  = "Respiratory deaths (daily)",
          pm10  = "PM10 (ug/m3)",
          temp  = "Temperature (C)",
          rhum  = "Relative humidity (%)")

summary_tab <- do.call(rbind, lapply(names(vars), function(v) {
  x <- data[[v]]
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  data.frame(
    variable  = vars[[v]],
    mean      = round(mean(x, na.rm = TRUE), 1),
    sd        = round(sd(x, na.rm = TRUE), 1),
    p25       = round(q[1], 1),
    median    = round(q[2], 1),
    p75       = round(q[3], 1),
    missing   = sum(is.na(x)),
    row.names = NULL
  )
}))

print(summary_tab)
write.csv(summary_tab, "results/summary_table.csv", row.names = FALSE)

# ---- Time-series plots ----------------------------------------
long <- rbind(
  data.frame(date = data$date, value = data$death, series = "Daily deaths"),
  data.frame(date = data$date, value = data$pm10,  series = "PM10 (ug/m3)"),
  data.frame(date = data$date, value = data$temp,  series = "Temperature (C)")
)

p_ts <- ggplot(long, aes(date, value)) +
  geom_line(linewidth = 0.2, colour = "grey40", na.rm = TRUE) +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Chicago 1987-2000: daily mortality, PM10 and temperature") +
  theme_minimal(base_size = 11)

ggsave("figures/01_time_series.png", p_ts, width = 9, height = 7, dpi = 150)

# ---- Seasonality of deaths by month ---------------------------
p_month <- ggplot(data, aes(factor(month), death)) +
  geom_boxplot(outlier.size = 0.4, fill = "#cfe3e8") +
  labs(x = "Month", y = "Daily deaths",
       title = "Seasonal pattern of daily mortality") +
  theme_minimal(base_size = 11)

ggsave("figures/02_deaths_by_month.png", p_month, width = 8, height = 4.5, dpi = 150)

# Note: the extreme values in July 1995 correspond to the Chicago heat wave.
cat("\nHighest daily death counts:\n")
print(head(data[order(-data$death), c("date", "death", "temp")], 5))
