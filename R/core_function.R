## JointShocks R Package
## List of functions (core_function file)

#' Read Data
#'
#' This utility loads a CSV or Excel file and preserves variable names exactly as in the file.
#' @param path File path to the dataset (.csv or .xlsx)
#' @param date_col Logical; whether to include a date column (default = TRUE)
#' @param date_format Optional format string for parsing dates (e.g. "\%m/\%d/\%Y")
#' @return A list with named `Data` (matrix of variables) and `date` (optional date vector)
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @importFrom utils read.csv tail
#' @importFrom stats coef median sd residuals
#' @importFrom vars VAR Phi
#' @examples
#' data(Volatility_data_pre)
#' path <- tempfile(fileext = ".csv")
#' write.csv(data.frame(date = format(attr(Volatility_data_pre, "date"), "%m/%d/%Y"),
#'                       Volatility_data_pre, check.names = FALSE),
#'           path, row.names = FALSE)
#' Data_date <- read_data(path, date_format = "%m/%d/%Y")
#' Data <- Data_date$Data
#' date <- Data_date$date
#' @export
read_data <- function(path, date_col = TRUE, date_format = "%m/%d/%Y") {
  if (!file.exists(path)) stop("File not found")

  if (grepl("\\.xlsx?$", path)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Please install the 'readxl' package to use Excel support.")
    }
    raw_data <- readxl::read_excel(path, .name_repair = "minimal")
  } else {
    raw_data <- read.csv(path, check.names = FALSE)
  }

  if (date_col) {
    dates <- as.Date(unlist(raw_data[[1]]), format = date_format)
    if (all(is.na(dates))) {
      stop("Failed to parse the date column using date_format = '", date_format,
           "'. Check that date_format matches the format of the dates in your file ",
           "(e.g. \"%m/%d/%Y\"), or set date_col = FALSE if your file has no date column.")
    }
    data <- as.matrix(raw_data[, -1])
    colnames(data) <- colnames(raw_data)[-1]
  } else {
    dates <- NULL
    data <- as.matrix(raw_data)
    colnames(data) <- colnames(raw_data)
  }

  return(list(Data = data, date = dates))
}

# Extract parameters from a VAR model
#
# This function extracts the required parameters from a fitted VAR model, including
# the covariance matrix of shocks, the VMA coefficient matrices, lagged coefficients,
# and constant terms.
#
# @param VAR_model A VAR object estimated using the vars package
# @param H number of A matrices (beyond the immediate h=0 one) to save.  If H=5, then there will be a total of 6 A matrices (for h=0,1,2,3,4,5).
# @return A list with sigma_u, A (VMA matrices), coefficients, and constants
extract_VAR_parameters <- function(VAR_model, H) {
  if (!(VAR_model$type %in% c("const", "none"))) {
    stop("extract_VAR_parameters() currently only supports VAR models estimated with type = 'const' or type = 'none'. ",
         "Re-estimate your VAR using vars::VAR(..., type = 'const') or vars::VAR(..., type = 'none').")
  }
  if (H < 1) {
    stop("H must be a strictly positive integer (H >= 1) in extract_VAR_parameters().")
  }
  summary_model <- summary(VAR_model)
  sigma_u <- summary_model$covres
  A <- Phi(VAR_model, nstep = H)    #Note: H must be greater than or equal to 1 for this to work (setting H=0 yields and error message for the extraction of the VMA coefficient matrices).

  coef_matrix <- coef(VAR_model)
  K <- ncol(sigma_u)
  p_lag <- VAR_model$p

  var_coef_matrices <- array(NA, dim = c(K, K, p_lag))
  for (i in 1:K) {
    for (lag in 1:p_lag) {
      start_index <- (lag - 1) * K + 1
      end_index <- lag * K
      var_coef_matrices[i, , lag] <- coef_matrix[[i]][start_index:end_index, 1]
    }
  }

  if (VAR_model$type == "const") {
    constant <- sapply(coef_matrix, function(mat) tail(mat[,1], 1))
  } else {
    constant <- rep(0, K)   #type == "none": the VAR has no constant/intercept term, so it is exactly zero.
  }

  return(list(sigma_u = sigma_u, A = A, coef_matrices = var_coef_matrices, constant = constant))
}



# Calculate Generalized Impulse Response Function (gIRF)
#
# @param A VMA coefficient matrices
# @param sigma_u Covariance matrix of residuals
# @param H time steps ahead (produces gIRF at h=0,1,2,...,H).  Thus, if H=5, then we see the gIRF at h=0,1,2,3,4,5.
# @return 3D array of gIRFs
calculate_gIRF <- function(A, sigma_u, H) {
  if (H < 0) {
    stop("H must be a non-negative integer in calculate_gIRF().")
  }
  if (dim(A)[1] != nrow(sigma_u) || dim(A)[2] != ncol(sigma_u)) {
    stop("The dimensions of A and sigma_u do not match in calculate_gIRF(). ",
         "Check that A and sigma_u both correspond to the same K-variable VAR model.")
  }
  if (dim(A)[3] < H + 1) {
    stop("A does not contain enough VMA coefficient matrices for the requested horizon H. ",
         "A must have at least H+1 matrices along its third dimension.")
  }
  G <- sqrt(diag(diag(sigma_u)))
  gIRF <- array(NA, dim = c(dim(A)[1], dim(A)[2], H+1))
  for (h in 1:(H+1)) {
    gIRF[,,h] <- A[,,h] %*% sigma_u %*% solve(G)
  }
  return(gIRF)
}


# Calculate Joint Impulse Response Function (jIRF)
#
# @param A VMA coefficient matrices
# @param sigma_u Covariance matrix of residuals
# @param e_sel Selection matrix (columns for shocked variables)
# @param shock_std_vec magnitude of the originating shocks
# @param H time steps ahead (produces jIRF at h=0,1,2,...,H).  Thus, if H=5, then we see the jIRF at h=0,1,2,3,4,5.
# @return Matrix of jIRFs
calculate_jIRF <- function(A, sigma_u, e_sel, shock_std_vec, H) {
  if (H < 0) {
    stop("H must be a non-negative integer in calculate_jIRF().")
  }
  if (dim(A)[1] != nrow(sigma_u) || dim(A)[2] != ncol(sigma_u)) {
    stop("The dimensions of A and sigma_u do not match in calculate_jIRF(). ",
         "Check that A and sigma_u both correspond to the same K-variable VAR model.")
  }
  if (dim(A)[3] < H + 1) {
    stop("A does not contain enough VMA coefficient matrices for the requested horizon H. ",
         "A must have at least H+1 matrices along its third dimension.")
  }
  if (ncol(e_sel) != length(shock_std_vec)) {
    stop("The number of columns in e_sel must match the length of shock_std_vec in calculate_jIRF(). ",
         "Each shocked variable needs exactly one column in e_sel and one entry in shock_std_vec.")
  }
  K <- dim(A)[1]
  jIRF <- matrix(NA, nrow = K, ncol = H + 1)
  for (h in 1:(H + 1)) {
    jIRF[, h] <- A[, , h] %*% sigma_u %*% e_sel %*% solve(t(e_sel) %*% sigma_u %*% e_sel) %*% shock_std_vec
  }
  return(jIRF)
}



# Calculate Generalized Forecast Error Variance Decomposition
#
# @param A VMA coefficient matrices
# @param sigma_u Covariance matrix
# @param XiCov H-step Forecast error covariance matrix
# @param H Forecast horizon
# @return List with FEVD
calculate_gFEVD <- function(A, sigma_u, XiCov, H) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  if (dim(A)[1] != nrow(sigma_u) || dim(A)[2] != ncol(sigma_u) ||
      dim(A)[1] != nrow(XiCov) || dim(A)[2] != ncol(XiCov)) {
    stop("The dimensions of A, sigma_u, and XiCov do not match in calculate_gFEVD(). ",
         "Check that A, sigma_u, and XiCov all correspond to the same K-variable VAR model.")
  }
  if (dim(A)[3] < H) {
    stop("A does not contain enough VMA coefficient matrices for the requested horizon H. ",
         "A must have at least H matrices along its third dimension.")
  }
  K <- ncol(sigma_u)
  A_sigma_sq <- array(0, dim = c(K, K, H))
  for (h in 1:H) {
    A_sigma_sq[,,h] <- (A[,,h] %*% sigma_u)^2    #Element by element square.
  }
  A_sigma_sum <- apply(A_sigma_sq, c(1, 2), sum)

  FEVD_gen_num <- matrix(0, K, K)
  FEVD_gen <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      FEVD_gen_num[i, j] <- A_sigma_sum[i, j] / sigma_u[j, j]
      FEVD_gen[i, j] <- FEVD_gen_num[i, j] / XiCov[i, i]
    }
  }
    return(FEVD_gen)
}



# Calculate Joint Forecast Error Variance Decomposition (jFEVD)
#
# @param A VMA coefficient matrices
# @param sigma_u Covariance matrix
# @param e_sel Selection matrix for shocks
# @param XiCov Forecast error covariance matrix
# @param H Forecast horizon
# @return Named vector of jFEVD
calculate_jFEVD <- function(A, sigma_u, e_sel, XiCov, H) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  if (dim(A)[1] != nrow(sigma_u) || dim(A)[2] != ncol(sigma_u) ||
      dim(A)[1] != nrow(XiCov) || dim(A)[2] != ncol(XiCov) ||
      dim(A)[1] != nrow(e_sel)) {
    stop("The dimensions of A, sigma_u, e_sel, and XiCov do not match in calculate_jFEVD(). ",
         "Check that they all correspond to the same K-variable VAR model.")
  }
  if (dim(A)[3] < H) {
    stop("A does not contain enough VMA coefficient matrices for the requested horizon H. ",
         "A must have at least H matrices along its third dimension.")
  }
  K <- dim(A)[1]
  idK <- diag(K)
  jFEVD_num <- numeric(K)  #Vector of all zeros the length of K.
  for (i in 1:K) {
    num_h <- numeric(H)   #Vector of all zeros the length of H. (numerator at each h)
    for (h in 1:H) {
      num_h[h] <- idK[i, ] %*% A[,, h] %*% sigma_u %*% e_sel %*% solve(t(e_sel) %*% sigma_u %*% e_sel) %*% t(e_sel) %*% sigma_u %*% t(A[,, h]) %*% idK[, i]
    }
    jFEVD_num[i] <- sum(num_h)
  }
  jFEVD <- jFEVD_num / diag(XiCov)   #Element by element division.
  names(jFEVD) <- paste0("Var", 1:K)
  return(jFEVD)
}



#' Run Joint IRF Calculation
#'
#' This is a wrapper function. It is "full" in the sense that it outputs how the whole
#' Y vector responds due to the from_set shocks. It computes the joint IRF directly
#' from a VAR model and a specified set of shock variables.
#'
#' @param VAR_model Fitted VAR model
#' @param from_set Character vector of variable names or list of integers for the variables where the shocks originate.
#' @param H Forecast horizon (can be equal to zero)
#' @return Matrix of joint IRFs with labeled row names
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Wiesen, T. F. P., and Beaumont, P. M. (2024). A joint impulse response function for
#' vector autoregressive models. Empirical Economics, 66(4), 1553-1585.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' jIRF <- run_full_jIRF(VAR_model, from_set = c("BCS", "CS", "DB"), H = 8)
#'
#' # Plot the jIRF response of a single variable (e.g. BAC) to the joint shocks
#' horizon <- 0:8
#' plot(horizon, jIRF["BAC", ], type = "l", lwd = 3, col = "#33a02c",
#'      xlab = "Horizon", ylab = "Response",
#'      main = "Joint IRF: Response of BAC to Joint Shock")
#' abline(h = 0, col = "gray40", lty = 2)
#' grid()
#' @export
run_full_jIRF <- function(VAR_model, from_set, H = 30) {
  setup <- extract_VAR_parameters(VAR_model, H+1)   #We want to allow people to set H=0 if they only want the immediate (h=0) IRF.  But if we enter H=0 into extract_VAR_parameters, then we get an error message.  Thus, let us use the extract_VAR_parameters command to get one more (H+1) VMA coefficient matrix than we need and not use it when we compute the IRF.
  var_names <- colnames(VAR_model$y)
  if (is.character(from_set)==1){from_indices = which(var_names %in% from_set)}   #If from_set is list of characters, then convert character names in from_set to indices (e.g., 2,4,6).
  if (is.character(from_set)==1  & sum(!from_set %in% var_names)>0){stop("One or more of your variable names are not recognized.  Check the spelling of your from_set variable names.  Or use integers to define the from_set (e.g., from_set =c(2,4) for the second and fourth variables).")}
  #The above line checks if the variable names in from_set are characters that do not exist in var_names and it yields an error message in that case.
  if (is.numeric(from_set)==1){from_indices=from_set}   #If from_set is list of integers, then simply set from_set equal to from_indices.
  if (is.numeric(from_set)==1 & (any(from_indices < 1) | any(from_indices > length(var_names)))) {stop("One or more of your from_set integers are out of range. from_set must contain integers between 1 and ", length(var_names), " (the number of variables in your VAR model).")}
  if (length(from_indices) < 2) {stop("from_set must contain at least two variables for the joint IRF (jIRF), since it measures the response to joint, simultaneous shocks from multiple variables. For a single shock variable, use run_full_gIRF() instead.")}
  e_sel <- diag(ncol(setup$sigma_u))[, from_indices]
  shock_std_vec <- sqrt(diag(setup$sigma_u[from_indices, from_indices]))
  jIRF <- calculate_jIRF(setup$A, setup$sigma_u, e_sel, shock_std_vec, H)   #Remember that H can be equal to zero here.
  rownames(jIRF) <- var_names
  return(jIRF)
}




#' Run full Generalized IRF Calculation
#'
#' This is a wrapper function. It is "full" in the sense that it outputs how the whole
#' Y vector responds due to shocks from all K variables (one at a time).
#'
#' @param VAR_model Fitted VAR object from vars::VAR()
#' @param H Forecast horizon
#' @return A labeled 3D array of generalized IRFs
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Koop, G., Pesaran, M. H., and Potter, S. M. (1996). Impulse response analysis in
#' nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.
#'
#' Pesaran, H. H., and Shin, Y. (1998). Generalized impulse response analysis in linear
#' multivariate models. Economics Letters, 58(1), 17-29.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' gIRF <- run_full_gIRF(VAR_model, H = 8)
#'
#' # Plot the response of BAC to a shock in BCS
#' horizon <- 0:8
#' plot(horizon, gIRF["BAC", "BCS", ], type = "l", lwd = 3, col = "#1f78b4",
#'      xlab = "Horizon", ylab = "Response",
#'      main = "Response of BAC to Shock in BCS")
#' abline(h = 0, col = "gray40", lty = 2)
#' grid()
#' @export
run_full_gIRF <- function(VAR_model, H = 30) {
  if (H < 0) {stop("H must be a non-negative integer in run_full_gIRF().")}
  setup <- extract_VAR_parameters(VAR_model, H+1)   #We want to allow people to set H=0 if they only want the immediate (h=0) IRF.  But if we enter H=0 into extract_VAR_parameters, then we get an error message.  Thus, let us use the extract_VAR_parameters command to get one more (H+1) VMA coefficient matrix than we need and not use it when we compute the IRF.
  gIRF <- calculate_gIRF(setup$A, setup$sigma_u, H)   #Remember that H can be equal to zero here.
  var_names <- colnames(VAR_model$y)
  if (!is.null(var_names)) {
    dimnames(gIRF) <- list(Response = var_names, Shock = var_names, Horizon = 0:H)
  }
  return(gIRF)
}



#' Run Full Generalized FEVD Calculation
#'
#' This is a wrapper for the gFEVD computation; it uses the calculate_gFEVD subfunction.
#' It is "full" in the sense that it computes the full K by K matrix (the i, j element
#' of the K by K matrix is the proportion of the forecast error variance of variable i
#' that can be explained by the shocks of variable j).
#' @param VAR_model Fitted VAR model
#' @param H Forecast horizon
#' @return A K x K matrix of generalized FEVD values (proportions), with dimnames
#' set to the VAR's variable names
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Koop, G., Pesaran, M. H., and Potter, S. M. (1996). Impulse response analysis in
#' nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.
#'
#' Pesaran, H. H., and Shin, Y. (1998). Generalized impulse response analysis in linear
#' multivariate models. Economics Letters, 58(1), 17-29.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' gFEVD <- run_full_gFEVD(VAR_model, H = 30)
#' @export
run_full_gFEVD <- function(VAR_model, H = 30) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  setup <- extract_VAR_parameters(VAR_model, H)
  K <- dim(setup$A)[1]
  var_names <- colnames(VAR_model$y)
  XiCov_h <- array(NA, dim = c(K, K, H))
  for (h in 1:H) {
    XiCov_h[,,h] <- setup$A[,,h] %*% setup$sigma_u %*% t(setup$A[,,h])
  }
  XiCov <- apply(XiCov_h, c(1,2), sum)
  result <- calculate_gFEVD(setup$A, setup$sigma_u, XiCov, H)
  dimnames(result) <- list(var_names, var_names)
  return(result)
}




#' Run Full Joint FEVD Calculation
#'
#' This function computes the joint forecast error variance decomposition (jFEVD)
#' directly from a VAR model and a specified set of jointly-conditioning shock
#' variables. For each variable in the system, it reports the proportion of that
#' variable's H-step forecast error variance that can be jointly explained by
#' simultaneous shocks from the variables in from_set.
#'
#' @param VAR_model Fitted VAR model
#' @param from_set Character vector of variable names or list of integers where the joint conditioning set
#' @param H Forecast horizon
#' @return Named vector of joint FEVD for each variable
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Wiesen, T. F. P., and Beaumont, P. M. (2024). A joint impulse response function for
#' vector autoregressive models. Empirical Economics, 66(4), 1553-1585.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' jfevd <- run_full_jFEVD(VAR_model, from_set = c("BCS", "CS", "DB"), H = 30)
#' @export
run_full_jFEVD <- function(VAR_model, from_set, H = 30) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  setup <- extract_VAR_parameters(VAR_model, H)
  var_names <- colnames(VAR_model$y)
  if (is.character(from_set)==1){from_indices = which(var_names %in% from_set)}   #If from_set is list of characters, then convert character names in from_set to indices (e.g., 3,5,6).
  if (is.character(from_set)==1  & sum(!from_set %in% var_names)>0){stop("One or more of your variable names are not recognized.  Check the spelling of your from_set variable names.  Or use integers to define the from_set (e.g., from_set =c(2,4) for the second and fourth variables).")}
  #The above line checks if the variable names in from_set are characters that do not exist in var_names and it yields an error message in that case.
  if (is.numeric(from_set)==1){from_indices=from_set}   #If from_set is list of integers, then simply set from_set equal to from_indices.
  if (is.numeric(from_set)==1 & (any(from_indices < 1) | any(from_indices > length(var_names)))) {stop("One or more of your from_set integers are out of range. from_set must contain integers between 1 and ", length(var_names), " (the number of variables in your VAR model).")}
  if (length(from_indices) < 2) {stop("from_set must contain at least two variables for the joint FEVD (jFEVD), since it measures the forecast error variance jointly explained by simultaneous shocks from multiple variables. For a single shock variable, use run_full_gFEVD() instead.")}
  e_sel <- diag(ncol(setup$sigma_u))[, from_indices]
  K <- ncol(setup$sigma_u)
  XiCov_h <- array(NA, dim = c(K, K, H))
  for (h in 1:H) {
    XiCov_h[,,h] <- setup$A[,,h] %*% setup$sigma_u %*% t(setup$A[,,h])
  }
  XiCov <- apply(XiCov_h, c(1,2), sum)
  jFEVD <- calculate_jFEVD(setup$A, setup$sigma_u, e_sel, XiCov, H)
  names(jFEVD) <- var_names
  return(jFEVD)
}



#' Bootstrap a VAR model
#'
#' Bootstraps the VAR model B times using a residual resampling procedure, saving the
#' constant vector, the P VAR coefficient matrices, and the shock covariance matrix for
#' each of the B iterations. The first iteration is always the original (non-resampled)
#' fitted parameters, so the original estimates can be recovered alongside the bootstrap
#' distribution.
#'
#' @param VAR_model Fitted VAR model from vars::VAR(), estimated with type = "const"
#' @param B Number of bootstrap iterations (default = 10000)
#' @return A list of class "bootstrap_VAR", where K is the number of variables in
#' VAR_model and p is its lag order, with: constant (a K x B matrix of constant terms,
#' one column per bootstrap iteration), coef_matrices (a K x K x p x B array of VAR
#' coefficient matrices, one K x K x p set per iteration), Sigma_u (a K x K x B array of
#' shock covariance matrices, one per iteration), p (the lag order), var_names (the K
#' variable names), and B (the number of bootstrap iterations)
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Wiesen, T. F. P., and Beaumont, P. M. (2024). A joint impulse response function for
#' vector autoregressive models. Empirical Economics, 66(4), 1553-1585.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' boot_model <- bootstrap_VAR(VAR_model, B = 500)
#' @export
bootstrap_VAR <- function(VAR_model, B = 10000) {
  if (B < 1) {stop("B must be a strictly positive integer in bootstrap_VAR().")}
  setup <- extract_VAR_parameters(VAR_model, 1)   #H=1 is the minimum allowed by extract_VAR_parameters(); we only need coef_matrices, constant, and sigma_u here, not the VMA matrices.
  Y <- t(VAR_model$y)
  Tt <- ncol(Y)
  K <- nrow(Y)
  p_lag <- VAR_model$p
  var_names <- colnames(VAR_model$y)
  resid <- t(residuals(VAR_model))

  constant <- matrix(NA, nrow = K, ncol = B)
  coef_matrices <- array(NA, dim = c(K, K, p_lag, B))
  Sigma_u <- array(NA, dim = c(K, K, B))

  #Force the first iteration to be the original (non-resampled) data and parameters.
  constant[, 1] <- setup$constant
  coef_matrices[, , , 1] <- setup$coef_matrices
  Sigma_u[, , 1] <- setup$sigma_u

  for (b in 2:B) {
    if (b %% 200 == 0) cat("Bootstrap iteration:", b, "\n")   #This provides a progress indicator on the bootstrap iterations.

    idx <- sample(1:(Tt - p_lag), replace = TRUE)
    boot_resid <- resid[, idx]
    boot_resid <- boot_resid - rowMeans(boot_resid)

    Ynew <- matrix(0, nrow = K, ncol = Tt)
    Ynew[, 1:p_lag] <- Y[, 1:p_lag]
    for (t in (p_lag + 1):Tt) {
      yhat <- setup$constant
      for (lag in 1:p_lag) {
        yhat <- yhat + setup$coef_matrices[, , lag] %*% Ynew[, t - lag]
      }
      Ynew[, t] <- yhat + boot_resid[, t - p_lag]
    }

    model_boot <- VAR(as.data.frame(t(Ynew)), p = p_lag, type = "const")
    setup_boot <- extract_VAR_parameters(model_boot, 1)

    constant[, b] <- setup_boot$constant
    coef_matrices[, , , b] <- setup_boot$coef_matrices
    Sigma_u[, , b] <- setup_boot$sigma_u
  }

  structure(list(
    constant = constant,
    coef_matrices = coef_matrices,
    Sigma_u = Sigma_u,
    p = p_lag,
    var_names = var_names,
    B = B
  ), class = "bootstrap_VAR")
}



# Convert VAR coefficient matrices to VMA (Wold) coefficient matrices
#
# Internal helper that converts a set of K x K x p VAR coefficient matrices into
# K x K x (H+1) VMA coefficient matrices (A[,,1] is the h=0 identity matrix,
# A[,,2] is the h=1 matrix, and so on). This is the same recursive conversion
# used inside genspill_roll() and jointspill_roll(), extracted here so that
# bootstrap_gIRF(), bootstrap_jIRF(), and bootstrap_gIRF_sum() can reuse it on
# bootstrapped VAR coefficient matrices (which, unlike a fitted VAR_model object,
# have no Phi() method available).
#
# @param coef_matrices K x K x p array of VAR coefficient matrices
# @param H Forecast horizon (produces VMA matrices for h=0,1,...,H)
# @return K x K x (H+1) array of VMA coefficient matrices
coef_matrices_to_VMA <- function(coef_matrices, H) {
  K <- dim(coef_matrices)[1]
  p <- dim(coef_matrices)[3]
  A <- array(NA, dim = c(K, K, H + 1))
  A[,,1] <- diag(K)
  coef_long <- array(0, dim = c(K, K, H + 1))
  for (lag in 1:p) {
    coef_long[,,lag] <- coef_matrices[,,lag]
  }
  for (h in 1:H) {
    A_times_coef <- array(0, dim = c(K, K, h))
    for (subh in 1:h) {
      A_times_coef[,,subh] <- A[,,h + 1 - subh] %*% coef_long[,,subh]
    }
    A[,,h + 1] <- apply(A_times_coef, c(1, 2), sum)
  }
  return(A)
}
# Callers must wrap a per-iteration/per-window slice in array(x, dim = dim(...)[1:3])
# before passing it here: indexing a K x K x p x N array at a single position silently
# drops the p dimension when p=1, which would otherwise make dim(coef_matrices)[3] NA.



# Resolve a response/shock variable name or index against a bootstrap_VAR object
#
# Internal helper shared by bootstrap_gIRF(), bootstrap_jIRF(), and bootstrap_gIRF_sum()
# to convert a response/shock argument (character variable name(s) or integer index/indices)
# into integer indices, with clear error messages for unrecognized names or out-of-range
# integers.
#
# @param x Character name(s) or integer index/indices
# @param var_names Character vector of variable names
# @param arg_name Name of the argument being validated, used in error messages
# @return Integer index/indices
resolve_var_indices <- function(x, var_names, arg_name) {
  if (is.character(x)) {
    if (sum(!x %in% var_names) > 0) {
      stop("One or more of your ", arg_name, " variable names are not recognized. ",
           "Check the spelling, or use integers instead (e.g., ", arg_name, " = 2).")
    }
    indices <- which(var_names %in% x)
  } else if (is.numeric(x)) {
    if (any(x < 1) | any(x > length(var_names))) {
      stop("One or more of your ", arg_name, " integers are out of range. ", arg_name,
           " must contain integers between 1 and ", length(var_names), ".")
    }
    indices <- x
  } else {
    stop(arg_name, " must be a character variable name or an integer index.")
  }
  return(indices)
}



#' Bootstrapped Generalized Impulse Response Function (gIRF) with confidence intervals
#'
#' Uses the bootstrapped VAR parameters from bootstrap_VAR() to compute B generalized
#' impulse response functions for a single response variable and a single shock variable,
#' and summarizes them with a confidence interval at each horizon.
#'
#' @param boot_model A "bootstrap_VAR" object returned by bootstrap_VAR()
#' @param response Name or integer index of the response variable
#' @param shock Name or integer index of the variable where the shock originates
#' @param H Forecast horizon (produces the gIRF at h=0,1,...,H)
#' @param CI Confidence interval level (default = 0.9)
#' @details
#' Each component of the returned list can be extracted directly with `$`, e.g.
#' `gIRF_CI$original`, and is a vector over the horizons h=0,...,H. original is the
#' gIRF computed from the original (non-bootstrapped) fitted VAR parameters, the
#' single best point estimate. mean, median, and sd are the mean, median, and standard
#' deviation of the gIRF across all bootstrap iterations. lower and upper are the
#' bounds of the CI-percent bootstrap confidence interval.
#' @return A list of class "bootstrap_gIRF" holding the following elements:
#' \item{lower}{CI lower bound, by horizon}
#' \item{mean}{Bootstrap mean, by horizon}
#' \item{median}{Bootstrap median, by horizon}
#' \item{original}{Original (non-bootstrapped) gIRF, by horizon}
#' \item{upper}{CI upper bound, by horizon}
#' \item{sd}{Bootstrap standard deviation, by horizon}
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Koop, G., Pesaran, M. H., and Potter, S. M. (1996). Impulse response analysis in
#' nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.
#'
#' Pesaran, H. H., and Shin, Y. (1998). Generalized impulse response analysis in linear
#' multivariate models. Economics Letters, 58(1), 17-29.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' boot_model <- bootstrap_VAR(VAR_model, B = 500)
#' gIRF_CI <- bootstrap_gIRF(boot_model, response = "BAC", shock = "BCS", H = 8, CI = 0.9)
#' summary(gIRF_CI)
#'
#' # Plot the original gIRF with its bootstrapped confidence interval
#' horizon <- 0:gIRF_CI$H
#' plot(horizon, gIRF_CI$original, type = "n",
#'      ylim = range(c(gIRF_CI$lower, gIRF_CI$upper)),
#'      xlab = "Horizon", ylab = "Response",
#'      main = "gIRF: Response of BAC due to a shock from BCS")
#' polygon(c(horizon, rev(horizon)), c(gIRF_CI$lower, rev(gIRF_CI$upper)),
#'         col = adjustcolor("purple", alpha.f = 0.25), border = NA)
#' lines(horizon, gIRF_CI$original, col = "purple", lwd = 2.5)
#' abline(h = 0, col = "gray40", lty = 2)
#' grid()
#' @export
bootstrap_gIRF <- function(boot_model, response, shock, H, CI = 0.9) {
  if (!inherits(boot_model, "bootstrap_VAR")) {stop("boot_model must be a 'bootstrap_VAR' object returned by bootstrap_VAR().")}
  if (H < 0) {stop("H must be a non-negative integer in bootstrap_gIRF().")}
  var_names <- boot_model$var_names
  response_idx <- resolve_var_indices(response, var_names, "response")
  shock_idx <- resolve_var_indices(shock, var_names, "shock")
  if (length(response_idx) != 1 | length(shock_idx) != 1) {stop("response and shock must each be a single variable in bootstrap_gIRF().")}

  B <- boot_model$B
  gIRF_boot <- matrix(NA, nrow = H + 1, ncol = B)
  for (b in 1:B) {
    coef_matrices_b <- array(boot_model$coef_matrices[, , , b], dim = dim(boot_model$coef_matrices)[1:3])
    A_b <- coef_matrices_to_VMA(coef_matrices_b, H)
    gIRF_b <- calculate_gIRF(A_b, boot_model$Sigma_u[, , b], H)
    gIRF_boot[, b] <- gIRF_b[response_idx, shock_idx, ]
  }

  alpha <- 1 - CI
  lower_idx <- round(B * (alpha / 2))
  upper_idx <- round(B * (1 - alpha / 2))
  gIRF_sorted <- apply(gIRF_boot, 1, sort)

  structure(list(
    lower = gIRF_sorted[lower_idx, ],
    mean = rowMeans(gIRF_boot),
    median = apply(gIRF_boot, 1, median),
    original = gIRF_boot[, 1],   #Recall that bootstrap_VAR() forces iteration 1 to be the original (non-resampled) parameters.
    upper = gIRF_sorted[upper_idx, ],
    sd = apply(gIRF_boot, 1, sd),
    response = var_names[response_idx],
    shock = var_names[shock_idx],
    var_names = var_names,
    H = H
  ), class = "bootstrap_gIRF")
}



#' Bootstrapped Joint Impulse Response Function (jIRF) with confidence intervals
#'
#' Uses the bootstrapped VAR parameters from bootstrap_VAR() to compute B joint impulse
#' response functions for a single response variable due to joint, simultaneous shocks
#' from multiple shock variables, and summarizes them with a confidence interval at each
#' horizon.
#'
#' @param boot_model A "bootstrap_VAR" object returned by bootstrap_VAR()
#' @param response Name or integer index of the response variable
#' @param shocks Names or integer indices of the variables where the joint shocks originate
#' (must contain at least two variables)
#' @param H Forecast horizon (produces the jIRF at h=0,1,...,H)
#' @param CI Confidence interval level (default = 0.9)
#' @details
#' Each component of the returned list can be extracted directly with `$`, e.g.
#' `jIRF_CI$original`, and is a vector over the horizons h=0,...,H. original is the
#' jIRF computed from the original (non-bootstrapped) fitted VAR parameters, the
#' single best point estimate. mean, median, and sd are the mean, median, and standard
#' deviation of the jIRF across all bootstrap iterations. lower and upper are the
#' bounds of the CI-percent bootstrap confidence interval.
#' @return A list of class "bootstrap_jIRF" holding the following elements:
#' \item{lower}{CI lower bound, by horizon}
#' \item{mean}{Bootstrap mean, by horizon}
#' \item{median}{Bootstrap median, by horizon}
#' \item{original}{Original (non-bootstrapped) jIRF, by horizon}
#' \item{upper}{CI upper bound, by horizon}
#' \item{sd}{Bootstrap standard deviation, by horizon}
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Wiesen, T. F. P., and Beaumont, P. M. (2024). A joint impulse response function for
#' vector autoregressive models. Empirical Economics, 66(4), 1553-1585.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' boot_model <- bootstrap_VAR(VAR_model, B = 500)
#' jIRF_CI <- bootstrap_jIRF(boot_model, response = "BAC", shocks = c("BCS", "CS", "DB"),
#'                            H = 8, CI = 0.9)
#' summary(jIRF_CI)
#'
#' # Plot the original jIRF with its bootstrapped confidence interval
#' horizon <- 0:jIRF_CI$H
#' plot(horizon, jIRF_CI$original, type = "n",
#'      ylim = range(c(jIRF_CI$lower, jIRF_CI$upper)),
#'      xlab = "Horizon", ylab = "Response",
#'      main = "jIRF: Response of BAC due to joint shocks from BCS, CS, and DB")
#' polygon(c(horizon, rev(horizon)), c(jIRF_CI$lower, rev(jIRF_CI$upper)),
#'         col = adjustcolor("black", alpha.f = 0.2), border = NA)
#' lines(horizon, jIRF_CI$original, col = "black", lwd = 2.5, lty = 2)
#' abline(h = 0, col = "gray40", lty = 2)
#' grid()
#' @export
bootstrap_jIRF <- function(boot_model, response, shocks, H, CI = 0.9) {
  if (!inherits(boot_model, "bootstrap_VAR")) {stop("boot_model must be a 'bootstrap_VAR' object returned by bootstrap_VAR().")}
  if (H < 0) {stop("H must be a non-negative integer in bootstrap_jIRF().")}
  var_names <- boot_model$var_names
  response_idx <- resolve_var_indices(response, var_names, "response")
  shock_idx <- resolve_var_indices(shocks, var_names, "shocks")
  if (length(response_idx) != 1) {stop("response must be a single variable in bootstrap_jIRF().")}
  if (length(shock_idx) < 2) {stop("shocks must contain at least two variables for the joint IRF (jIRF), since it measures the response to joint, simultaneous shocks from multiple variables. For a single shock variable, use bootstrap_gIRF() instead.")}

  K <- length(var_names)
  e_sel <- diag(K)[, shock_idx]
  B <- boot_model$B
  jIRF_boot <- matrix(NA, nrow = H + 1, ncol = B)
  for (b in 1:B) {
    coef_matrices_b <- array(boot_model$coef_matrices[, , , b], dim = dim(boot_model$coef_matrices)[1:3])
    A_b <- coef_matrices_to_VMA(coef_matrices_b, H)
    sigma_u_b <- boot_model$Sigma_u[, , b]
    shock_std_b <- sqrt(diag(sigma_u_b[shock_idx, shock_idx]))
    jIRF_b <- calculate_jIRF(A_b, sigma_u_b, e_sel, shock_std_b, H)
    jIRF_boot[, b] <- jIRF_b[response_idx, ]
  }

  alpha <- 1 - CI
  lower_idx <- round(B * (alpha / 2))
  upper_idx <- round(B * (1 - alpha / 2))
  jIRF_sorted <- apply(jIRF_boot, 1, sort)

  structure(list(
    lower = jIRF_sorted[lower_idx, ],
    mean = rowMeans(jIRF_boot),
    median = apply(jIRF_boot, 1, median),
    original = jIRF_boot[, 1],   #Recall that bootstrap_VAR() forces iteration 1 to be the original (non-resampled) parameters.
    upper = jIRF_sorted[upper_idx, ],
    sd = apply(jIRF_boot, 1, sd),
    response = var_names[response_idx],
    shocks = var_names[shock_idx],
    var_names = var_names,
    H = H
  ), class = "bootstrap_jIRF")
}



#' Bootstrapped sum of Generalized Impulse Response Functions (gIRF sum) with confidence intervals
#'
#' Uses the bootstrapped VAR parameters from bootstrap_VAR() to compute, at each of the B
#' bootstrap iterations, the sum of the individual generalized impulse response functions
#' for a single response variable due to separate (not joint) shocks from multiple shock
#' variables, and summarizes them with a confidence interval at each horizon.
#'
#' @param boot_model A "bootstrap_VAR" object returned by bootstrap_VAR()
#' @param response Name or integer index of the response variable
#' @param shocks Names or integer indices of the variables where the individual shocks
#' originate (must contain at least two variables)
#' @param H Forecast horizon (produces the gIRF sum at h=0,1,...,H)
#' @param CI Confidence interval level (default = 0.9)
#' @details
#' Each component of the returned list can be extracted directly with `$`, e.g.
#' `gIRF_sum_CI$original`, and is a vector over the horizons h=0,...,H. original is
#' the gIRF sum computed from the original (non-bootstrapped) fitted VAR parameters,
#' the single best point estimate. mean, median, and sd are the mean, median, and
#' standard deviation of the gIRF sum across all bootstrap iterations. lower and upper
#' are the bounds of the CI-percent bootstrap confidence interval.
#' @return A list of class "bootstrap_gIRF_sum" holding the following elements:
#' \item{lower}{CI lower bound, by horizon}
#' \item{mean}{Bootstrap mean, by horizon}
#' \item{median}{Bootstrap median, by horizon}
#' \item{original}{Original (non-bootstrapped) gIRF sum, by horizon}
#' \item{upper}{CI upper bound, by horizon}
#' \item{sd}{Bootstrap standard deviation, by horizon}
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Koop, G., Pesaran, M. H., and Potter, S. M. (1996). Impulse response analysis in
#' nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.
#'
#' Pesaran, H. H., and Shin, Y. (1998). Generalized impulse response analysis in linear
#' multivariate models. Economics Letters, 58(1), 17-29.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' boot_model <- bootstrap_VAR(VAR_model, B = 500)
#' gIRF_sum_CI <- bootstrap_gIRF_sum(boot_model, response = "BAC",
#'                                    shocks = c("BCS", "CS", "DB"), H = 8, CI = 0.9)
#' summary(gIRF_sum_CI)
#'
#' # Plot the original sum of gIRFs with its bootstrapped confidence interval
#' horizon <- 0:gIRF_sum_CI$H
#' plot(horizon, gIRF_sum_CI$original, type = "n",
#'      ylim = range(c(gIRF_sum_CI$lower, gIRF_sum_CI$upper)),
#'      xlab = "Horizon", ylab = "Response",
#'      main = "gIRF sum: Response of BAC to individual shocks from BCS, CS, and DB")
#' polygon(c(horizon, rev(horizon)), c(gIRF_sum_CI$lower, rev(gIRF_sum_CI$upper)),
#'         col = adjustcolor("purple", alpha.f = 0.25), border = NA)
#' lines(horizon, gIRF_sum_CI$original, col = "purple", lwd = 2.5)
#' abline(h = 0, col = "gray40", lty = 2)
#' grid()
#' @export
bootstrap_gIRF_sum <- function(boot_model, response, shocks, H, CI = 0.9) {
  if (!inherits(boot_model, "bootstrap_VAR")) {stop("boot_model must be a 'bootstrap_VAR' object returned by bootstrap_VAR().")}
  if (H < 0) {stop("H must be a non-negative integer in bootstrap_gIRF_sum().")}
  var_names <- boot_model$var_names
  response_idx <- resolve_var_indices(response, var_names, "response")
  shock_idx <- resolve_var_indices(shocks, var_names, "shocks")
  if (length(response_idx) != 1) {stop("response must be a single variable in bootstrap_gIRF_sum().")}
  if (length(shock_idx) < 2) {stop("shocks must contain at least two variables in bootstrap_gIRF_sum(), since it sums the individual gIRFs due to shocks from multiple variables. For a single shock variable, use bootstrap_gIRF() instead.")}

  B <- boot_model$B
  gIRF_sum_boot <- matrix(NA, nrow = H + 1, ncol = B)
  for (b in 1:B) {
    coef_matrices_b <- array(boot_model$coef_matrices[, , , b], dim = dim(boot_model$coef_matrices)[1:3])
    A_b <- coef_matrices_to_VMA(coef_matrices_b, H)
    gIRF_b <- calculate_gIRF(A_b, boot_model$Sigma_u[, , b], H)
    gIRF_sum_boot[, b] <- apply(gIRF_b[response_idx, shock_idx, , drop = FALSE], 3, sum)
  }

  alpha <- 1 - CI
  lower_idx <- round(B * (alpha / 2))
  upper_idx <- round(B * (1 - alpha / 2))
  gIRF_sum_sorted <- apply(gIRF_sum_boot, 1, sort)

  structure(list(
    lower = gIRF_sum_sorted[lower_idx, ],
    mean = rowMeans(gIRF_sum_boot),
    median = apply(gIRF_sum_boot, 1, median),
    original = gIRF_sum_boot[, 1],   #Recall that bootstrap_VAR() forces iteration 1 to be the original (non-resampled) parameters.
    upper = gIRF_sum_sorted[upper_idx, ],
    sd = apply(gIRF_sum_boot, 1, sd),
    response = var_names[response_idx],
    shocks = var_names[shock_idx],
    var_names = var_names,
    H = H
  ), class = "bootstrap_gIRF_sum")
}



# Internal helper to print a bootstrap IRF confidence interval table
#
# Shared formatting helper used by summary.bootstrap_gIRF(), summary.bootstrap_jIRF(),
# and summary.bootstrap_gIRF_sum(). Prints a title line followed by a table with one row
# per horizon (h=0,1,...,H) and one column for each of: CI lower bound, mean, median,
# original, CI upper bound, and standard deviation.
#
# @param object One of "bootstrap_gIRF", "bootstrap_jIRF", or "bootstrap_gIRF_sum"
# @param title Character string describing the IRF, printed above the table
print_bootstrap_IRF_summary <- function(object, title) {
  H <- object$H
  full_table <- matrix(NA, nrow = H + 1, ncol = 6)
  full_table[, 1] <- format(round(object$lower, 4), nsmall = 4)
  full_table[, 2] <- format(round(object$mean, 4), nsmall = 4)
  full_table[, 3] <- format(round(object$median, 4), nsmall = 4)
  full_table[, 4] <- format(round(object$original, 4), nsmall = 4)
  full_table[, 5] <- format(round(object$upper, 4), nsmall = 4)
  full_table[, 6] <- format(round(object$sd, 4), nsmall = 4)

  row_labels <- paste0("h=", 0:H)
  col_labels <- c("CI lower bound", "mean", "median", "original", "CI upper bound", "stand. dev.")
  dimnames(full_table) <- list(row_labels, col_labels)

  cat(title, "\n\n")
  df <- as.data.frame(full_table, check.names = FALSE)
  print(df, right = TRUE, row.names = TRUE)
  invisible(full_table)
}



#' Summary method for bootstrap_gIRF object
#'
#' Prints a table with one row per horizon and one column for each of the CI lower
#' bound, mean, median, original (non-bootstrapped), CI upper bound, and standard
#' deviation of the bootstrapped generalized impulse response function.
#'
#' @param object A "bootstrap_gIRF" object returned by bootstrap_gIRF()
#' @param ... Not used
#' @rdname bootstrap_gIRF
#' @export
summary.bootstrap_gIRF <- function(object, ...) {
  if (!inherits(object, "bootstrap_gIRF")) {stop("Object is not of class 'bootstrap_gIRF'")}
  title <- paste0("gIRF: The generalized response of ", object$response,
                   " due to a shock from ", object$shock, ".")
  print_bootstrap_IRF_summary(object, title)
}



#' Summary method for bootstrap_jIRF object
#'
#' Prints a table with one row per horizon and one column for each of the CI lower
#' bound, mean, median, original (non-bootstrapped), CI upper bound, and standard
#' deviation of the bootstrapped joint impulse response function.
#'
#' @param object A "bootstrap_jIRF" object returned by bootstrap_jIRF()
#' @param ... Not used
#' @return A formatted table of the bootstrapped jIRF confidence interval
#' @rdname bootstrap_jIRF
#' @export
summary.bootstrap_jIRF <- function(object, ...) {
  if (!inherits(object, "bootstrap_jIRF")) {stop("Object is not of class 'bootstrap_jIRF'")}
  shocks_text <- if (length(object$shocks) > 1) {
    paste0(paste(object$shocks[-length(object$shocks)], collapse = ", "), ", & ", object$shocks[length(object$shocks)])
  } else {
    object$shocks
  }
  title <- paste0("jIRF: The response of ", object$response,
                   " due to joint simultaneous shocks from ", shocks_text, ".")
  print_bootstrap_IRF_summary(object, title)
}



#' Summary method for bootstrap_gIRF_sum object
#'
#' Prints a table with one row per horizon and one column for each of the CI lower
#' bound, mean, median, original (non-bootstrapped), CI upper bound, and standard
#' deviation of the bootstrapped sum of generalized impulse response functions.
#'
#' @param object A "bootstrap_gIRF_sum" object returned by bootstrap_gIRF_sum()
#' @param ... Not used
#' @return A formatted table of the bootstrapped gIRF sum confidence interval
#' @rdname bootstrap_gIRF_sum
#' @export
summary.bootstrap_gIRF_sum <- function(object, ...) {
  if (!inherits(object, "bootstrap_gIRF_sum")) {stop("Object is not of class 'bootstrap_gIRF_sum'")}
  shocks_text <- if (length(object$shocks) > 1) {
    paste0(paste(object$shocks[-length(object$shocks)], collapse = ", "), ", & ", object$shocks[length(object$shocks)])
  } else {
    object$shocks
  }
  title <- paste0("gIRF sum: The sum of the generalized responses of ", object$response,
                   " due to individual shocks from ", shocks_text, ".")
  print_bootstrap_IRF_summary(object, title)
}




#' Generalized Spillover Table and Metrics (Diebold & Yilmaz, 2012/2014)
#'
#' Computes the generalized spillover table (gSOT) and related total/net spillover
#' metrics of Diebold and Yilmaz (2012, 2014) from a fitted VAR model, over the full
#' sample. The gSOT is a row-normalized version of the generalized forecast error
#' variance decomposition (gFEVD), from which the total spillover to/from each
#' variable, the net spillover for each variable, and the overall generalized
#' spillover index (gSOI) are derived.
#'
#' @param VAR_model Estimated VAR object from vars::VAR()
#' @param H Forecast horizon (H >= 1)
#' @details
#' Each component of the returned list can be extracted directly with `$`, e.g.
#' `genspill_results$gSOI`. gFEVD is the unscaled K x K generalized forecast error
#' variance decomposition matrix, as a proportion, where element (i, j) is the
#' proportion of variable i's forecast error variance explained by shocks from
#' variable j alone. gSOT is gFEVD row-normalized to sum to 100 across each row, as
#' percentages - Diebold and Yilmaz's row-sum-normalized bilateral spillover measure.
#' gen_from_others_to_i gives, for each variable, the generalized total spillover it
#' receives from all other variables. gen_to_others_from_i gives the generalized total
#' spillover each variable transmits to all others. gen_net is gen_to_others_from_i
#' minus gen_from_others_to_i, where a positive value means the variable transmits
#' more shocks than it receives. gSOI, the generalized spillover index, summarizes
#' overall generalized connectedness across the whole system in one number.
#' @return A list of class "genspill" holding the following elements:
#' \item{gFEVD}{Unscaled K x K generalized FEVD matrix (proportions)}
#' \item{gSOT}{Row-normalized generalized spillover table (percentages)}
#' \item{gen_from_others_to_i}{Generalized spillover received by each variable}
#' \item{gen_to_others_from_i}{Generalized spillover transmitted by each variable}
#' \item{gen_net}{Net generalized spillover for each variable}
#' \item{gSOI}{Generalized spillover index}
#' \item{var_names}{Variable names, in order}
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Diebold, F. X., and Yilmaz, K. (2012). Better to give than to receive: Predictive
#' directional measurement of volatility spillovers. International Journal of
#' Forecasting, 28(1), 57-66.
#'
#' Diebold, F. X., and Yilmaz, K. (2014). On the network topology of variance
#' decompositions: Measuring the connectedness of financial firms. Journal of
#' Econometrics, 182(1), 119-134.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' genspill_results <- genspill(VAR_model, H = 30)
#' summary(genspill_results)
#' @export
genspill <- function(VAR_model, H=30) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  setup <- extract_VAR_parameters(VAR_model, H)
  sigma_u <- setup$sigma_u
  A <- setup$A
  K <- ncol(sigma_u)

  var_names <- colnames(VAR_model$y)
  if (is.null(var_names)) var_names <- paste0("Var", 1:K)

  XiCov_h <- array(NA, dim = c(K, K, H))
  for (h in 1:H) {
    XiCov_h[,,h] <- A[,,h] %*% sigma_u %*% t(A[,,h])
  }
  XiCov <- apply(XiCov_h, c(1,2), sum)

  FEVD_gen <- calculate_gFEVD(A, sigma_u, XiCov, H)

  row_sum <- rowSums(FEVD_gen)
  gSOT <- 100 * (FEVD_gen / row_sum)
  gSOI <- 100 - sum(diag(gSOT)) / K

  gen_from_others_to_i <- 100 - diag(gSOT)
  gen_to_others_from_i <- colSums(gSOT) - diag(gSOT)
  gen_net <- gen_to_others_from_i - gen_from_others_to_i

  rownames(gSOT) <- colnames(gSOT) <- var_names
  names(gen_from_others_to_i) <- var_names
  names(gen_to_others_from_i) <- var_names
  names(gen_net) <- var_names
  rownames(FEVD_gen) <- colnames(FEVD_gen) <- var_names

  structure(list(
    gSOT = gSOT,
    gen_from_others_to_i = gen_from_others_to_i,
    gen_to_others_from_i = gen_to_others_from_i,
    gen_net = gen_net,
    gSOI = gSOI,
    gFEVD = FEVD_gen,
    var_names = var_names
  ), class = "genspill")
}




#' Summary method for genspill object
#'
#' Prints the gSOT matrix together with the generalized total spillover to/from each
#' variable, the generalized net spillover for each variable, and the overall
#' generalized spillover index (gSOI), in a single formatted table.
#'
#' @param object A genspill object returned by `genspill()`
#' @param ... Not used
#' @return A formatted generalized spillover table
#' @rdname genspill
#' @export
summary.genspill <- function(object, ...) {
  if (!inherits(object, "genspill")) {
    stop("Object is not of class 'genspill'")
  }

  K <- nrow(object$gSOT)
  var_names <- object$var_names

  full_table <- matrix(NA, nrow = K + 2, ncol = K + 1)
  full_table[1:K, 1:K] <- format(round(object$gSOT, 2), nsmall = 2)
  full_table[1:K, K + 1] <- format(round(object$gen_from_others_to_i, 2), nsmall = 2)
  full_table[K + 1, 1:K] <- format(round(object$gen_to_others_from_i, 2), nsmall = 2)
  full_table[K + 2, 1:K] <- format(round(object$gen_net, 2), nsmall = 2)
  full_table[K + 2, K + 1] <- format(round(object$gSOI, 2), nsmall = 2)
  full_table[K + 1, K+1] <- "gSOI"

  row_labels <- c(var_names, "To others from i", "NET")
  col_labels <- c(var_names, "From others to i")
  dimnames(full_table) <- list(row_labels, col_labels)

  df <- as.data.frame(full_table, check.names = FALSE)
  print(df, right = TRUE, row.names = TRUE)
  invisible(full_table)
}



#' Joint Spillover Table and Metrics (Lastrapes & Wiesen, 2021, 2025)
#'
#' Computes the joint spillover metrics of Lastrapes and Wiesen (2021, 2025) from a
#' fitted VAR model, over the full sample. Unlike the generalized spillover metrics,
#' which condition on one shock variable at a time, these metrics measure the
#' proportion of each variable's forecast error variance that can be jointly explained
#' by simultaneous shocks from all other variables at once, and report the joint
#' spillover index (jSOI) summarizing this across the whole system.
#'
#' @param VAR_model Estimated VAR object from vars::VAR()
#' @param H Forecast horizon (default=30)
#' @details
#' Each component of the returned list can be extracted directly with `$`, e.g.
#' `jointspill_results$jSOI`. gFEVD is the unscaled K x K generalized forecast error
#' variance decomposition matrix, as percentages, where element (i, j) is the
#' percentage of variable i's forecast error variance explained by shocks from
#' variable j alone. gFEVD_tilde is that same matrix rescaled by lambda so it is on a
#' comparable scale to the joint spillover metrics. jnt_from_others_to_i gives, for
#' each variable, the percentage of its forecast error variance jointly explained by
#' simultaneous shocks from all other variables at once. jnt_to_others_from_i gives
#' the joint total spillover each variable transmits to all others. jnt_net is
#' jnt_to_others_from_i minus jnt_from_others_to_i, where a positive value means the
#' variable transmits more shocks than it receives. jSOI, the joint spillover index,
#' is the average of jnt_from_others_to_i across all variables.
#' @return A list of class "jointspill" holding the following elements:
#' \item{gFEVD}{Unscaled K x K generalized FEVD matrix (percentages)}
#' \item{gFEVD_tilde}{gFEVD rescaled by lambda}
#' \item{lambda}{Scaling factor applied to gFEVD}
#' \item{jnt_from_others_to_i}{Joint spillover received by each variable}
#' \item{jnt_to_others_from_i}{Joint spillover transmitted by each variable}
#' \item{jnt_net}{Net joint spillover for each variable}
#' \item{jSOI}{Joint spillover index}
#' \item{var_names}{Variable names, in order}
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Lastrapes, W. D., and Wiesen, T. F. P. (2021). The joint spillover index. Economic
#' Modelling, 94, 681-691.
#'
#' Lastrapes, W. D., and Wiesen, T. F. P. (2025). Regional bank failures and volatility
#' transmission. Journal of Financial Stability, 78, 101404.
#' @examples
#' data(Volatility_data_pre)
#' VAR_model <- vars::VAR(Volatility_data_pre, p = 1, type = "const")
#' jointspill_results <- jointspill(VAR_model, H = 30)
#' summary(jointspill_results)
#' @export
jointspill <- function(VAR_model, H=30) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  setup <- extract_VAR_parameters(VAR_model, H)
  sigma_u <- setup$sigma_u
  A <- setup$A
  K <- ncol(sigma_u)
  e <- diag(1, K)

  XiCov_h <- array(NA, dim = c(K, K, H))
  for (h in 1:H) {
    XiCov_h[,,h] <- A[,,h] %*% sigma_u %*% t(A[,,h])
  }
  XiCov <- apply(XiCov_h, c(1,2), sum)

  FEVD_gen <- calculate_gFEVD(A, sigma_u, XiCov, H)

  M <- vector("list", K)
  M[[1]] <- e[, -1]
  for (i in 2:(K-1)) {
    M[[i]] <- cbind(e[, 1:(i-1)], e[, (i+1):K])
  }
  M[[K]] <- e[, -K]

  Sjnt_DOT_to_i_numerator <- numeric(K)
  for (i in 1:K) {
    for (h in 1:H) {
      Mi <- M[[i]]
      term <- solve(t(Mi) %*% sigma_u %*% Mi)
      Sjnt_DOT_to_i_numerator[i] <- Sjnt_DOT_to_i_numerator[i] +t(e[, i]) %*% A[,,h] %*% sigma_u %*% Mi %*% term %*% t(Mi) %*% sigma_u %*% t(A[,,h]) %*% e[, i]
    }
  }
  jnt_from_others_to_i <- 100 * (Sjnt_DOT_to_i_numerator / diag(XiCov))
  jSOI <- mean(jnt_from_others_to_i)

  FEVD_gen_percent <- 100 * FEVD_gen
  lambda <- jSOI / ((1 / K) * (sum(FEVD_gen_percent) - sum(diag(FEVD_gen_percent))))
  gFEVD_tilde <- lambda * FEVD_gen_percent

  jnt_to_others_from_i <- numeric(K)
  for (i in 1:K) {
    jnt_to_others_from_i[i] <- sum(gFEVD_tilde[, i]) - gFEVD_tilde[i, i]
  }
  jnt_net <- jnt_to_others_from_i - jnt_from_others_to_i

  var_names <- colnames(VAR_model$y)
  if (is.null(var_names)) {
    var_names <- paste0("Var", 1:K)
  }

  structure(list(
    gFEVD = FEVD_gen_percent,
    jnt_from_others_to_i = jnt_from_others_to_i,
    jnt_to_others_from_i = jnt_to_others_from_i,
    jnt_net = jnt_net,
    jSOI = jSOI,
    lambda = lambda,
    gFEVD_tilde = gFEVD_tilde,
    var_names = var_names
  ), class = "jointspill")
}



#' Summary method for jointspill object
#'
#' Prints the gFEVD (or gFEVD_tilde) matrix together with the joint total spillover
#' to/from each variable, the joint net spillover for each variable, and the overall
#' joint spillover index (jSOI), in a single formatted table.
#'
#' @param object A jointspill object returned by `jointspill()`
#' @param report_gFEVD Logical. If TRUE (default), the table reports the unscaled
#' gFEVD matrix (as percentages) - the proportion of each variable's forecast error
#' variance explained by shocks from one other variable at a time. If FALSE, the table
#' reports the gFEVD_tilde matrix instead - the same information rescaled by lambda so
#' it is on a comparable scale to the joint spillover metrics (jnt_from_others_to_i,
#' jnt_to_others_from_i, jnt_net, jSOI) shown in the rest of the table.
#' @param ... Not used
#' @return A formatted joint spillover table
#' @rdname jointspill
#' @export
summary.jointspill <- function(object, report_gFEVD=TRUE,...) {
  if (!inherits(object, "jointspill")) {
    stop("Object is not of class 'jointspill'")
  }
  if (!is.logical(report_gFEVD)) {stop("Input argument `report_gFEVD' should be TRUE or FALSE.  Set report_gFEVD=TRUE to report the unscaled gFEVD matrix (reported as percentages).  Set report_gFEVD=FALSE to report the scaled gFEVD_tilde matrix.")}

  K = length(object$jnt_from_others_to_i)
  full_table <- matrix(NA, nrow = K + 2, ncol = K + 1)
  if (report_gFEVD) {full_table[1:K, 1:K] = format(round(object$gFEVD, 2), nsmall = 2)}
  if (!report_gFEVD) {full_table[1:K, 1:K] = format(round(object$gFEVD_tilde, 2), nsmall = 2)}
  full_table[1:K, K + 1] = format(round(object$jnt_from_others_to_i, 2), nsmall = 2)
  full_table[K + 1, 1:K] = format(round(object$jnt_to_others_from_i, 2), nsmall = 2)
  full_table[K + 2, 1:K] = format(round(object$jnt_net, 2), nsmall = 2)
  full_table[K + 2, K + 1] = format(round(object$jSOI, 2), nsmall = 2)
  full_table[K+1, K+1]= "jSOI"

  var_names = object$var_names
  row_labels = c(var_names, "To others from i", "NET")
  col_labels = c(var_names, "From others to i")
  dimnames(full_table) = list(row_labels, col_labels)

  df = as.data.frame(full_table, check.names = FALSE)
  print(df, right = TRUE, row.names = TRUE)
  invisible(full_table)
}




#' Rolling Window VAR Estimation
#'
#' Estimates a separate VAR model over each rolling window of a fixed length win,
#' sliding the window forward one observation at a time across the full sample, and
#' saves the fitted VAR coefficients and shock covariance matrix for each window. The
#' output of this function is used as the input to genspill_roll() and
#' jointspill_roll() to compute spillover metrics over rolling windows.
#'
#' @param data A K-variable dataset (data frame or matrix)
#' @param p Number of lags in each VAR
#' @param win Rolling window size
#' @param const Logical; whether to include a constant term (default = TRUE)
#' @param dates Date vector of the original VAR data.
#' @return A list, where K is the number of variables in data, p is the lag order, and N
#' is the number of rolling windows, with: VARcoeff (a K x K x p x N array of VAR
#' coefficient matrices, one K x K x p set per window), constant (a K x N matrix of
#' constant terms, one column per window; only included if const = TRUE), Sigma (a K x K
#' x N array of shock covariance matrices, one per window), and end_of_window_date (if
#' dates was supplied, the date at the end of each of the N windows)
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @examples
#' data(Volatility_data_pre)
#' date <- attr(Volatility_data_pre, "date")
#' rolling_model <- VAR_roll(data = Volatility_data_pre, p = 1, win = 150, const = TRUE, dates = date)
#' @export
VAR_roll <- function(data, p, win, const = TRUE, dates = NULL) {
  Tt <- nrow(data)
  K <- ncol(data)
  wnum <- Tt + 1 - win
  if (p < 1) {stop("p (the number of lags) must be a strictly positive integer in VAR_roll().")}
  if (win <= p) {stop("The rolling window length (win) must be greater than the number of lags (p). Please select a larger win or a smaller p.")}
  if (win>Tt) {stop("The length of the rolling windows (win) cannot exceed the sample size of the time series data (T).  Please select a smaller rolling window length (win).")  }
  if (!is.null(dates)) { if (Tt != length(dates)) {stop("The length of the time series data (T) does not match the length of the dates.  Check data or dates input.  Or set dates equal to NULL.")}  }   #Checks if the length of the data and the length of the dates match.
  VARcoeff <- array(NA, dim = c(K, K, p, wnum))
  constant <- if (const) matrix(NA, nrow = K, ncol = wnum) else NULL
  Sigma <- array(NA, dim = c(K, K, wnum))

  for (tau in 1:wnum) {
    window_data <- data[tau:(tau + win - 1), ]
    VAR_model <- VAR(window_data, p = p, type = if (const) "const" else "none")
    coef_list <- coef(VAR_model)

    for (i in 1:K) {
      for (lag in 1:p) {
        lag_start <- (lag - 1) * K + 1
        lag_end <- lag_start + K - 1
        VARcoeff[i, , lag, tau] <- coef_list[[i]][lag_start:lag_end, 1]
      }
      if (const) {
        constant[i, tau] <- coef_list[[i]][nrow(coef_list[[i]]), 1]
      }
    }
    Sigma[,,tau] <- summary(VAR_model)$covres
    if ((tau %% 50) == 0) {cat("rolling window progress (%):", round(100*(tau/wnum),2), "\n")}
  }
  cat("rolling window progress (%):", 100, "\n")
  var_names <- colnames(data)
  if (!is.null(var_names)) {
    dimnames(VARcoeff) <- list(var_names, var_names, NULL, NULL)
    dimnames(Sigma) <- list(var_names, var_names, NULL)
    if (const) {rownames(constant) <- var_names}
  }
  result <- list(VARcoeff = VARcoeff,
                 Sigma = Sigma,
                 end_of_window_date = if (!is.null(dates)) dates[(length(dates) - wnum + 1):length(dates)] else NULL
                 )
  if (const) {result$constant <- constant}
  return(result)
}




#' Generalized Spillover Metrics Over Rolling Windows (Diebold & Yilmaz, 2012/2014)
#'
#' Computes the same generalized spillover metrics as genspill() (the gSOT, total/net
#' spillovers, and gSOI), but separately for each rolling window of a VAR_roll() object,
#' so that the evolution of generalized spillovers over time can be tracked and plotted.
#'
#' @param rolling_VAR Estimated rolling window VAR object from VAR_roll
#' @param H Forecast horizon (H >= 1)
#' @return A list of generalized spillover metrics over rolling windows, where N is the
#' number of rolling windows, holding the following elements:
#' \item{gSOT}{Row-normalized generalized spillover table (percentages), one K x K
#' matrix per window}
#' \item{gFEVD}{Unscaled K x K generalized FEVD matrix (proportions), one per window}
#' \item{gen_from_others_to_i}{Generalized spillover received by each variable, one K x
#' N matrix}
#' \item{gen_to_others_from_i}{Generalized spillover transmitted by each variable, one K
#' x N matrix}
#' \item{gen_net}{Net generalized spillover for each variable, one K x N matrix}
#' \item{gSOI}{Generalized spillover index, one value per window}
#' \item{end_of_window_date}{The end-of-window date for each rolling window, carried
#' through from rolling_VAR}
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Diebold, F. X., and Yilmaz, K. (2012). Better to give than to receive: Predictive
#' directional measurement of volatility spillovers. International Journal of
#' Forecasting, 28(1), 57-66.
#'
#' Diebold, F. X., and Yilmaz, K. (2014). On the network topology of variance
#' decompositions: Measuring the connectedness of financial firms. Journal of
#' Econometrics, 182(1), 119-134.
#' @examples
#' data(Volatility_data_pre)
#' date <- attr(Volatility_data_pre, "date")
#' rolling_model <- VAR_roll(data = Volatility_data_pre, p = 1, win = 150, const = TRUE, dates = date)
#' gen_spill_result_roll <- genspill_roll(rolling_model, H = 30)
#'
#' # Plot the generalized spillover index (gSOI) over rolling windows
#' plot(gen_spill_result_roll$end_of_window_date, gen_spill_result_roll$gSOI,
#'      type = "l", col = "blue", lwd = 1,
#'      xlab = "End of window date", ylab = "gSOI",
#'      main = "gSOI over rolling windows")
#' @export
genspill_roll = function(rolling_VAR, H) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  VARcoeff_tau = rolling_VAR$VARcoeff
  Sigma_tau = rolling_VAR$Sigma
  K = dim(Sigma_tau)[1]
  p = dim(VARcoeff_tau)[3]
  wnum = dim(Sigma_tau)[3]
  var_names = rownames(Sigma_tau)

  gFEVD_tau = array(NA, dim = c(K, K, wnum))
  gSOT_tau = array(NA, dim = c(K, K, wnum))
  gen_from_others_to_i_tau = matrix(NA, nrow = K, ncol = wnum)
  gen_to_others_from_i_tau = matrix(NA, nrow = K, ncol = wnum)
  gen_net_tau = matrix(NA, nrow = K, ncol = wnum)
  gSOI_tau = numeric(wnum)

  for (tau in 1:wnum) {
    VARcoeff_tau_slice = array(VARcoeff_tau[,,,tau], dim = dim(VARcoeff_tau)[1:3])
    A = coef_matrices_to_VMA(VARcoeff_tau_slice, H)   #Note that A is recalculated for each tau

    #compute Xi
    XiCov_h = array(NA, dim = c(K, K, H))
    for (h in 1:H) {
      XiCov_h[,,h] = A[,,h] %*% Sigma_tau[,,tau] %*% t(A[,,h])
    }
    XiCov = apply(XiCov_h, c(1, 2), sum)

    #Compute gFEVD
    FEVD_gen = calculate_gFEVD(A, Sigma_tau[,,tau], XiCov, H)
    gFEVD_tau[,,tau]=FEVD_gen  #Save the gFEVD for iteration tau.

    #Compute the gSOT and gSOI
    row_sum = rowSums(FEVD_gen)
    gSOT = 100 * (FEVD_gen / row_sum)
    gSOI = 100 - sum(diag(gSOT)) / K
    gSOT_tau[,,tau] = gSOT #Save the gSOT for iteration tau.
    gSOI_tau[tau] = gSOI   #Save the gSOI for iteration tau.

    #Compute the generalized total to other, from other, and net spillovers.
    gen_from_others_to_i = 100-diag(gSOT)
    gen_to_others_from_i = colSums(gSOT)-diag(gSOT)
    gen_net = gen_to_others_from_i-gen_from_others_to_i
    gen_from_others_to_i_tau[,tau]=gen_from_others_to_i  #Save the gen total from others to i spillover for iteration tau.
    gen_to_others_from_i_tau[,tau]=gen_to_others_from_i  #Save the gen total to others from i spillover for iteration tau.
    gen_net_tau[,tau]=gen_net                            #Save the gen net total spillover for iteration tau.
  }
  if (!is.null(var_names)) {
    dimnames(gSOT_tau) <- list(var_names, var_names, NULL)
    dimnames(gFEVD_tau) <- list(var_names, var_names, NULL)
    rownames(gen_from_others_to_i_tau) <- var_names
    rownames(gen_to_others_from_i_tau) <- var_names
    rownames(gen_net_tau) <- var_names
  }
  return(list(
    gSOT = gSOT_tau,
    gFEVD = gFEVD_tau,
    gen_from_others_to_i = gen_from_others_to_i_tau,
    gen_to_others_from_i = gen_to_others_from_i_tau,
    gen_net = gen_net_tau,
    gSOI = gSOI_tau,
    end_of_window_date = rolling_VAR$end_of_window_date
  ))
}




#' Joint Spillover Metrics Over Rolling Windows (Lastrapes and Wiesen, 2021, 2025)
#'
#' Computes the same joint spillover metrics as jointspill() (joint total/net
#' spillovers, lambda, gFEVD_tilde, and jSOI), but separately for each rolling window of
#' a VAR_roll() object, so that the evolution of joint spillovers over time can be
#' tracked and plotted.
#'
#' @param rolling_VAR Estimated rolling window VAR object from VAR_roll
#' @param H Forecast horizon (H >= 1)
#' @return A list of joint spillover metrics over rolling windows, where N is the number
#' of rolling windows, holding the following elements:
#' \item{gFEVD}{Unscaled K x K generalized FEVD matrix (percentages), one per window}
#' \item{jnt_from_others_to_i}{Joint spillover received by each variable, one K x N
#' matrix}
#' \item{jnt_to_others_from_i}{Joint spillover transmitted by each variable, one K x N
#' matrix}
#' \item{jnt_net}{Net joint spillover for each variable, one K x N matrix}
#' \item{jSOI}{Joint spillover index, one value per window}
#' \item{lambda}{Scaling factor applied to gFEVD, one value per window}
#' \item{gFEVD_tilde}{gFEVD rescaled by lambda, one K x K matrix per window}
#' \item{end_of_window_date}{The end-of-window date for each rolling window, carried
#' through from rolling_VAR}
#' @author Johnson Oliyide, Thomas F. P. Wiesen
#' @references
#' Lastrapes, W. D., and Wiesen, T. F. P. (2021). The joint spillover index. Economic
#' Modelling, 94, 681-691.
#'
#' Lastrapes, W. D., and Wiesen, T. F. P. (2025). Regional bank failures and volatility
#' transmission. Journal of Financial Stability, 78, 101404.
#' @examples
#' data(Volatility_data_pre)
#' date <- attr(Volatility_data_pre, "date")
#' rolling_model <- VAR_roll(data = Volatility_data_pre, p = 1, win = 150, const = TRUE, dates = date)
#' joint_spill_results_roll <- jointspill_roll(rolling_model, H = 30)
#'
#' # Plot the joint spillover index (jSOI) over rolling windows
#' plot(joint_spill_results_roll$end_of_window_date, joint_spill_results_roll$jSOI,
#'      type = "l", col = "blue", lwd = 1,
#'      xlab = "End of window date", ylab = "jSOI",
#'      main = "jSOI over rolling windows")
#' @export
jointspill_roll <- function(rolling_VAR, H) {
  if (H<1) {stop("The forecast horizon (H) must be a strictly positive integer.")}
  VARcoeff_tau = rolling_VAR$VARcoeff
  Sigma_tau = rolling_VAR$Sigma
  K = dim(Sigma_tau)[1]
  p = dim(VARcoeff_tau)[3]
  wnum = dim(Sigma_tau)[3]
  var_names = rownames(Sigma_tau)

  gFEVD_percent_tau = array(NA, dim = c(K, K, wnum))
  jSOI_tau = numeric(wnum)
  lambda_tau = numeric(wnum)
  gFEVD_tilde_tau = array(NA, dim = c(K, K, wnum))
  jnt_from_others_to_i_tau = matrix(NA, nrow = K, ncol = wnum)
  jnt_to_others_from_i_tau = matrix(NA, nrow = K, ncol = wnum)
  jnt_net_tau = matrix(NA, nrow = K, ncol = wnum)

  e = diag(1, K)
  M = vector("list", length = K)
  M[[1]] = e[,-1]
  for (i in 2:(K-1)) {
    M[[i]] = cbind(e[,1:(i-1)], e[,(i+1):K])
  }
  M[[K]] = e[,-K]

  for (tau in 1:wnum) {
    VARcoeff_tau_slice = array(VARcoeff_tau[,,,tau], dim = dim(VARcoeff_tau)[1:3])
    A = coef_matrices_to_VMA(VARcoeff_tau_slice, H)   #Note that A is recalculated for each tau

    #Compute Xi
    XiCov_h = array(NA, dim = c(K, K, H))
    for (h in 1:H) {
      XiCov_h[,,h] = A[,,h] %*% Sigma_tau[,,tau] %*% t(A[,,h])
    }
    XiCov = apply(XiCov_h, c(1, 2), sum)

    #Compute gFEVD
    FEVD_gen = calculate_gFEVD(A, Sigma_tau[,,tau], XiCov, H)
    gFEVD_percent=100*FEVD_gen              #The gFEVD (as a percent)
    gFEVD_percent_tau[,,tau]=gFEVD_percent  #Save the gFEVD (as a percent) for iteration tau.

    # Joint Total Spillover FROM others TO i
    Sjnt_DOT_to_i_numerator = numeric(K)
    for (i in 1:K) {
      for (h in 1:H) {
        Mi = M[[i]]
        term = solve(t(Mi) %*% Sigma_tau[,,tau] %*% Mi)
        Sjnt_DOT_to_i_numerator[i] = Sjnt_DOT_to_i_numerator[i] +t(e[, i]) %*% A[,,h] %*% Sigma_tau[,,tau] %*% Mi %*% term %*% t(Mi) %*% Sigma_tau[,,tau] %*% t(A[,,h]) %*% e[, i]
      }
    }
    jnt_from_others_to_i = 100 * (Sjnt_DOT_to_i_numerator / diag(XiCov))
    jnt_from_others_to_i_tau[,tau]=jnt_from_others_to_i   #Save the the joint total spillover from others to i for window tau.

    jSOI = mean(jnt_from_others_to_i)
    jSOI_tau[tau] = jSOI     #save the jSOI for window tau.

    lambda = jSOI / ((1 / K) * (sum(gFEVD_percent) - sum(diag(gFEVD_percent))))
    gFEVD_tilde = lambda * gFEVD_percent
    lambda_tau[tau]=lambda      #Save lambda for window tau.
    gFEVD_tilde_tau[,,tau]=gFEVD_tilde      #Save gFEVD_tilde for window tau.

    for (i in 1:K) {
      jnt_to_others_from_i_tau[i, tau] = sum(gFEVD_tilde_tau[, i, tau]) - gFEVD_tilde_tau[i, i, tau]
      jnt_net_tau[i, tau] = jnt_to_others_from_i_tau[i, tau] - jnt_from_others_to_i_tau[i, tau]
    }
  }

  if (!is.null(var_names)) {
    dimnames(gFEVD_percent_tau) <- list(var_names, var_names, NULL)
    dimnames(gFEVD_tilde_tau) <- list(var_names, var_names, NULL)
    rownames(jnt_from_others_to_i_tau) <- var_names
    rownames(jnt_to_others_from_i_tau) <- var_names
    rownames(jnt_net_tau) <- var_names
  }
  return(list(
    gFEVD = gFEVD_percent_tau,
    jnt_from_others_to_i = jnt_from_others_to_i_tau,
    jnt_to_others_from_i = jnt_to_others_from_i_tau,
    jnt_net = jnt_net_tau,
    jSOI = jSOI_tau,
    lambda = lambda_tau,
    gFEVD_tilde = gFEVD_tilde_tau,
    end_of_window_date = rolling_VAR$end_of_window_date
  ))
}

