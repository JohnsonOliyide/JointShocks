# JointShocks

**JointShocks** is an R package for estimating joint and generalized impulse response functions, forecast error variance decompositions, and spillover metrics for vector autoregressive (VAR) models, over both full-sample and rolling-window applications.

It provides tools to compute:

- The **joint spillover metrics** of Lastrapes and Wiesen (2021, 2025)
- The **joint impulse response function** and **joint forecast error variance decomposition** of Wiesen and Beaumont (2024)
- The **generalized spillover index**, **generalized impulse response function**, and **generalized forecast error variance decomposition** of Diebold and Yilmaz (2012, 2014) and Koop, Pesaran, and Potter (1996) / Pesaran and Shin (1998)

## Step 1: Install the devtools package

To install an R package from GitHub, start by installing the `devtools` package. The best way to do this is from CRAN, by typing:

```r
install.packages("devtools")
```

## Step 2: Install JointShocks from GitHub

```r
library(devtools)
install_github("JohnsonOliyide/JointShocks")
```

## Step 3: Go through the tutorial

JointShocks Tutorial (link to come)

The tutorial walks through estimating a VAR model, computing joint and generalized impulse response functions, joint and generalized forecast error variance decompositions, and spillover metrics over both the full sample and rolling windows, with full worked examples using the bundled bank volatility dataset.

## Authors

Johnson Oliyide and Thomas F. P. Wiesen

## References

Diebold, F. X., and Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. *International Journal of Forecasting*, 28(1), 57-66.

Diebold, F. X., and Yilmaz, K. (2014). On the network topology of variance decompositions: Measuring the connectedness of financial firms. *Journal of Econometrics*, 182(1), 119-134.

Koop, G., Pesaran, M. H., and Potter, S. M. (1996). Impulse response analysis in nonlinear multivariate models. *Journal of Econometrics*, 74(1), 119-147.

Lastrapes, W. D., and Wiesen, T. F. P. (2021). The joint spillover index. *Economic Modelling*, 94, 681-691.

Lastrapes, W. D., and Wiesen, T. F. P. (2025). Regional bank failures and volatility transmission. *Journal of Financial Stability*, 78, 101404.

Pesaran, H. H., and Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. *Economics Letters*, 58(1), 17-29.

Wiesen, T. F. P., and Beaumont, P. M. (2024). A joint impulse response function for vector autoregressive models. *Empirical Economics*, 66(4), 1553-1585.
