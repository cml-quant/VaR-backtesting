# Load required libraries
library(rugarch)
library(lubridate)
library(reshape2)
library(ggplot2)


data <- read.csv("project2_data_42416.csv", sep=';')
data$date <- as.Date(data$date)

# we. create data frames for Returns and Prices
tickers <- unique(data$TICKER)
dates <- unique(data$date)
Returns <- data.frame(date = dates)
Prices <- data.frame(date = dates)

for (ticker in tickers) {
  stock_data <- data[data$TICKER == ticker, ]
  stock_data <- stock_data[order(stock_data$date), ]
  Returns[, ticker] <- stock_data$RET
  Prices[, ticker] <- stock_data$PRC
}

# we remove the date column and use rownames
rownames(Returns) <- as.character(Returns$date)
rownames(Prices) <- as.character(Prices$date)
Returns$date <- NULL
Prices$date <- NULL

save(Returns, file = "Returns.RData")
save(Prices, file = "Prices.RData")


# Parameters
p <- 0.01  # probability level
value <- 1  # portfolio value
T <- nrow(Returns)  # total observations
WE <- 100  # estimation window
WT <- T - WE  # testing window
assets <- c("JPM", "BCS", "DB")  # Modify if assets differ
y <- Returns[, assets]

# Historical Simulation VaR
calc_HS <- function(data, window, p, value) {
  T <- length(data)
  var <- rep(NA, T)
  for (t in (window + 1):T) {
    rolling_window <- data[(t - window):(t - 1)]
    var[t] <- -quantile(rolling_window, probs = p, na.rm = TRUE) * value
  }
  return(var)
}

# GARCH VaR (Normal and Student-t)
calc_GARCH <- function(data, window, p, value, distribution = "norm") {
  T <- length(data)
  var <- rep(NA, T)
  spec <- ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    distribution.model = distribution
  )
  
  for (t in (window + 1):T) {
    rolling_window <- data[(t - window):(t - 1)]
    tryCatch({
      fit <- ugarchfit(spec = spec, data = rolling_window, solver = "hybrid")
      if (fit@fit$convergence == 0) {
        sigma <- sigma(fit)
        if (distribution == "norm") {
          var[t] <- -sigma * qnorm(p)
        } else if (distribution == "std") {
          shape <- coef(fit)["shape"]
          var[t] <- -sigma * qdist("std", p, shape = shape)
        }
      }
    }, error = function(e) {
      warning(paste("GARCH fitting failed at time", t, ":", e$message))
    })
  }
  return(var)
}

# we thenalculate VaR for each asset
VaR <- list()
for (asset in assets) {
  VaR[[asset]] <- list(
    HS = calc_HS(y[, asset], WE, p, value),
    GARCH_norm = calc_GARCH(y[, asset], WE, p, value, "norm"),
    GARCH_t = calc_GARCH(y[, asset], WE, p, value, "std")
  )
}

# Use parameter extraction for JPM (GARCH-t model)
library(xts)

jpm_data <- xts(y[, "JPM"], order.by=as.Date(rownames(y)))

# GARCH spec for Student-t
spec_t <- ugarchspec(
  variance.model = list(model="sGARCH", garchOrder=c(1,1)),
  mean.model = list(armaOrder=c(0,0), include.mean=FALSE),
  distribution.model="std"
)

# we create storage for parameters and volatility
params <- matrix(NA, nrow=length(y[,1])-WE, ncol=5)
colnames(params) <- c("omega", "alpha", "beta", "shape", "volatility")
dates <- as.Date(rownames(y)[(WE+1):length(y[,1])])

# we ecxtract new parameters and volatility by using a rolling window
for(t in 1:(length(y[,1])-WE)) {
  window_data <- y[t:(t+WE-1), "JPM"]
  
  tryCatch({
    fit <- ugarchfit(spec=spec_t, data=window_data, solver="hybrid")
    if(fit@fit$convergence == 0) {
      params[t,1:4] <- coef(fit)
      params[t,5] <- tail(sigma(fit), 1)
    }
  }, error=function(e) {
    warning(paste("Fitting failed at time", t))
  })
}

# Create data frame with dates for parameters
params_df <- data.frame(
  date = dates,
  omega = params[,"omega"],
  alpha = params[,"alpha"],
  beta = params[,"beta"],
  volatility = params[,"volatility"]
)

# Plot of theparameters
params_melt <- melt(params_df, id.vars="date")
param_plot <- ggplot(params_melt, aes(x=date, y=value, color=variable)) +
  geom_line() +
  facet_wrap(~variable, scales="free_y", ncol=2) +
  theme_minimal() +
  labs(title="t-GARCH Parameters and Volatility Over Time - JPM",
       x="Date",
       y="Parameter Value") +
  theme(legend.position="none")
print(param_plot)

# We create our functions for violation test and Kupiec
calc_violations <- function(returns, var, p) {
  violations <- returns < -var
  violations[is.na(violations)] <- FALSE
  return(violations)
}

kupiec_test <- function(violations, p, T) {
  x <- sum(violations, na.rm = TRUE)
  if (x == 0) return(list(statistic = NA, p_value = NA))
  
  pi_hat <- x/T
  LR <- -2 * log(((1-p)^(T-x)) * (p^x)) + 
    2 * log(((1-pi_hat)^(T-x)) * (pi_hat^x))
  p_value <- 1 - pchisq(LR, df = 1)
  
  return(list(statistic = LR, p_value = p_value))
}

# we apply backtesting for each asset
backtest_results <- list()
for (asset in assets) {
  returns <- y[, asset]
  
  # Calculate actual losses for backtest
  losses <- -returns * value
  
  # Get the different VaR estimates
  var_hs <- VaR[[asset]]$HS
  var_garch_norm <- VaR[[asset]]$GARCH_norm
  var_garch_t <- VaR[[asset]]$GARCH_t
  
  # wecalculate violations
  violations_hs <- calc_violations(returns, var_hs/value, p)
  violations_garch_norm <- calc_violations(returns, var_garch_norm/value, p)
  violations_garch_t <- calc_violations(returns, var_garch_t/value, p)
  
  # we then do the Kupiec test
  test_period <- (WE+1):T
  kupiec_hs <- kupiec_test(violations_hs[test_period], p, length(test_period))
  kupiec_garch_norm <- kupiec_test(violations_garch_norm[test_period], p, length(test_period))
  kupiec_garch_t <- kupiec_test(violations_garch_t[test_period], p, length(test_period))
  
  # we store thetest results
  backtest_results[[asset]] <- list(
    violations = list(
      HS = sum(violations_hs, na.rm = TRUE),
      GARCH_norm = sum(violations_garch_norm, na.rm = TRUE),
      GARCH_t = sum(violations_garch_t, na.rm = TRUE)
    ),
    violation_rate = list(
      HS = mean(violations_hs, na.rm = TRUE),
      GARCH_norm = mean(violations_garch_norm, na.rm = TRUE),
      GARCH_t = mean(violations_garch_t, na.rm = TRUE)
    ),
    kupiec_test = list(
      HS = kupiec_hs,
      GARCH_norm = kupiec_garch_norm,
      GARCH_t = kupiec_garch_t
    )
  )
  
  #  VaR comparison plots 
  dates <- as.Date(rownames(y))
  plot_data <- data.frame(
    Date = dates[test_period],
    Returns = returns[test_period],
    VaR_HS = var_hs[test_period]/value,
    VaR_GARCH_norm = var_garch_norm[test_period]/value,
    VaR_GARCH_t = var_garch_t[test_period]/value
  )
  
  var_plot <- ggplot(plot_data, aes(x = Date)) +
    geom_line(aes(y = Returns), color = "black", alpha = 0.5) +
    geom_line(aes(y = -VaR_HS, color = "HS")) +
    geom_line(aes(y = -VaR_GARCH_norm, color = "GARCH-norm")) +
    geom_line(aes(y = -VaR_GARCH_t, color = "GARCH-t"), linetype = "dashed") +
    scale_color_manual(values = c("HS" = "blue", "GARCH-norm" = "red", "GARCH-t" = "green")) +
    labs(title = paste("VaR Backtesting -", asset),
         x = "Date",
         y = "Returns",
         color = "VaR Model") +
    theme_minimal()
  
  print(var_plot)
}

# Print backtesting results
cat("\nBACKTESTING RESULTS\n")
cat("===================\n")

for (asset in assets) {
  cat("\nResults for", asset, "\n")
  cat("-------------------\n")
  
  cat("Number of Violations:\n")
  print(backtest_results[[asset]]$violations)
  
  cat("\nViolation Rates:\n")
  print(backtest_results[[asset]]$violation_rate)
  
  cat("\nKupiec Test p-values:\n")
  p_values <- sapply(backtest_results[[asset]]$kupiec_test, function(x) x$p_value)
  print(p_values)
  
  cat("\nInterpretation:\n")
  for (model in c("HS", "GARCH_norm", "GARCH_t")) {
    cat(sprintf("%s model: ", model))
    p_val <- backtest_results[[asset]]$kupiec_test[[model]]$p_value
    if (is.na(p_val)) {
      cat("Insufficient violations for testing\n")
    } else if (p_val < 0.05) {
      cat("Reject H0 - model maybe misspecified (p < 0.05)\n")
    } else {
      cat("Fail to reject H0 - model appears adequate (p >= 0.05)\n")
    }
  }
  cat("\n")
}
