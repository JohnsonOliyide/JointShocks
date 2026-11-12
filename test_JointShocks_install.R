## JointShocks package - verification script
## Each section below is the exact @examples code from the reference
## manual, in the same order the manual lists them (its table of
## contents). Only the package loading line has been changed
## (library() instead of devtools).

library(JointShocks)

# ---- bootstrap_gIRF ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
boot_model <- bootstrap_VAR(VAR_model, B = 500)
gIRF_CI <- bootstrap_gIRF(boot_model, response = "BAC", shock = "BCS", H = 8, CI = 0.9)
summary(gIRF_CI)

# Plot the original gIRF with its bootstrapped confidence interval
horizon <- 0:gIRF_CI$H
plot(horizon, gIRF_CI$original, type = "n",
     ylim = range(c(gIRF_CI$lower, gIRF_CI$upper)),
     xlab = "Horizon", ylab = "Response",
     main = "gIRF: Response of BAC due to a shock from BCS")
polygon(c(horizon, rev(horizon)), c(gIRF_CI$lower, rev(gIRF_CI$upper)),
        col = adjustcolor("purple", alpha.f = 0.25), border = NA)
lines(horizon, gIRF_CI$original, col = "purple", lwd = 2.5)
abline(h = 0, col = "gray40", lty = 2)
grid()

# ---- bootstrap_gIRF_sum ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
boot_model <- bootstrap_VAR(VAR_model, B = 500)
gIRF_sum_CI <- bootstrap_gIRF_sum(boot_model, response = "BAC",
                                   shocks = c("BCS", "CS", "DB"), H = 8, CI = 0.9)
summary(gIRF_sum_CI)

# Plot the original sum of gIRFs with its bootstrapped confidence interval
horizon <- 0:gIRF_sum_CI$H
plot(horizon, gIRF_sum_CI$original, type = "n",
     ylim = range(c(gIRF_sum_CI$lower, gIRF_sum_CI$upper)),
     xlab = "Horizon", ylab = "Response",
     main = "gIRF sum: Response of BAC to individual shocks from BCS, CS, and DB")
polygon(c(horizon, rev(horizon)), c(gIRF_sum_CI$lower, rev(gIRF_sum_CI$upper)),
        col = adjustcolor("purple", alpha.f = 0.25), border = NA)
lines(horizon, gIRF_sum_CI$original, col = "purple", lwd = 2.5)
abline(h = 0, col = "gray40", lty = 2)
grid()

# ---- bootstrap_jIRF ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
boot_model <- bootstrap_VAR(VAR_model, B = 500)
jIRF_CI <- bootstrap_jIRF(boot_model, response = "BAC", shocks = c("BCS", "CS", "DB"),
                           H = 8, CI = 0.9)
summary(jIRF_CI)

# Plot the original jIRF with its bootstrapped confidence interval
horizon <- 0:jIRF_CI$H
plot(horizon, jIRF_CI$original, type = "n",
     ylim = range(c(jIRF_CI$lower, jIRF_CI$upper)),
     xlab = "Horizon", ylab = "Response",
     main = "jIRF: Response of BAC due to joint shocks from BCS, CS, and DB")
polygon(c(horizon, rev(horizon)), c(jIRF_CI$lower, rev(jIRF_CI$upper)),
        col = adjustcolor("black", alpha.f = 0.2), border = NA)
lines(horizon, jIRF_CI$original, col = "black", lwd = 2.5, lty = 2)
abline(h = 0, col = "gray40", lty = 2)
grid()

# ---- bootstrap_VAR ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
boot_model <- bootstrap_VAR(VAR_model, B = 500)

# ---- genspill ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
genspill_results <- genspill(VAR_model, H = 30)
summary(genspill_results)

# ---- genspill_roll ----
data(Volatility_data_pre)
date <- attr(Volatility_data_pre, "date")
rolling_model <- VAR_roll(data = Volatility_data_pre, p = 1, win = 150, const = TRUE, dates = date)
gen_spill_result_roll <- genspill_roll(rolling_model, H = 30)

# Plot the generalized spillover index (gSOI) over rolling windows
plot(gen_spill_result_roll$end_of_window_date, gen_spill_result_roll$gSOI,
     type = "l", col = "blue", lwd = 1,
     xlab = "End of window date", ylab = "gSOI",
     main = "gSOI over rolling windows")

# ---- jointspill ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
jointspill_results <- jointspill(VAR_model, H = 30)
summary(jointspill_results)

# ---- jointspill_roll ----
data(Volatility_data_pre)
date <- attr(Volatility_data_pre, "date")
rolling_model <- VAR_roll(data = Volatility_data_pre, p = 1, win = 150, const = TRUE, dates = date)
joint_spill_results_roll <- jointspill_roll(rolling_model, H = 30)

# Plot the joint spillover index (jSOI) over rolling windows
plot(joint_spill_results_roll$end_of_window_date, joint_spill_results_roll$jSOI,
     type = "l", col = "blue", lwd = 1,
     xlab = "End of window date", ylab = "jSOI",
     main = "jSOI over rolling windows")

# ---- read_data ----
data(Volatility_data_pre)
path <- tempfile(fileext = ".csv")
write.csv(data.frame(date = format(attr(Volatility_data_pre, "date"), "%m/%d/%Y"),
                      Volatility_data_pre, check.names = FALSE),
          path, row.names = FALSE)
Data_date <- read_data(path, date_format = "%m/%d/%Y")
Data <- Data_date$Data
date <- Data_date$date

# ---- run_full_gFEVD ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
gFEVD <- run_full_gFEVD(VAR_model, H = 30)

# ---- run_full_gIRF ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
gIRF <- run_full_gIRF(VAR_model, H = 8)

# Plot the response of BAC to a shock in BCS
horizon <- 0:8
plot(horizon, gIRF["BAC", "BCS", ], type = "l", lwd = 3, col = "#1f78b4",
     xlab = "Horizon", ylab = "Response",
     main = "Response of BAC to Shock in BCS")
abline(h = 0, col = "gray40", lty = 2)
grid()

# ---- run_full_jFEVD ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
jfevd <- run_full_jFEVD(VAR_model, from_set = c("BCS", "CS", "DB"), H = 30)

# ---- run_full_jIRF ----
data(Volatility_data_pre)
VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
jIRF <- run_full_jIRF(VAR_model, from_set = c("BCS", "CS", "DB"), H = 8)

# Plot the jIRF response of a single variable (e.g. BAC) to the joint shocks
horizon <- 0:8
plot(horizon, jIRF["BAC", ], type = "l", lwd = 3, col = "#33a02c",
     xlab = "Horizon", ylab = "Response",
     main = "Joint IRF: Response of BAC to Joint Shock")
abline(h = 0, col = "gray40", lty = 2)
grid()

# ---- VAR_roll ----
data(Volatility_data_pre)
date <- attr(Volatility_data_pre, "date")
rolling_model <- VAR_roll(data = Volatility_data_pre, p = 1, win = 150, const = TRUE, dates = date)
