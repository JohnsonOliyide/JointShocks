\name{Volatility_data_fullsample}
\alias{Volatility_data_fullsample}
\docType{data}
\title{Full-Sample Bank Volatility Data (Wiesen and Beaumont, 2024)}
\description{
Daily realized volatility of eleven major banks over the full sample period,
1/5/2017 to 2/3/2021 (spanning both the pre-pandemic and pandemic periods), used
throughout the package examples and vignette to replicate the results of Wiesen and
Beaumont (2024).
}
\usage{data(Volatility_data_fullsample)}
\format{
A matrix with 1027 rows (dates) and 11 columns, one per bank: Bank of America (BAC),
Bank of New York Mellon (BK), Barclays (BCS), Citigroup (C), Credit Suisse (CS),
Deutsche Bank (DB), Goldman Sachs (GS), JP Morgan Chase (JPM), Morgan Stanley (MS),
State Street (STT), and Wells Fargo (WFC). Row names are the trading dates.
}
\source{
Wiesen, T. F. P., and Beaumont, P. M. (2024). A joint impulse response function for
vector autoregressive models. Empirical Economics, 66(4), 1553-1585.
}
\examples{
data(Volatility_data_fullsample)
VAR_model <- vars::VAR(Volatility_data_fullsample, p = 1, type = "const")
}
\keyword{datasets}
