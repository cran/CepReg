#' @import MASS
#' @import rrpack
#' @import psych
#' @import Renvlp
#' @importFrom stats mvfft quantile rnorm

##################################################
#' @title Compute the Periodogram of Multivariate Time Series
#'
#' @description This function computes the periodogram for each time series in the input matrix.
#'
#' @param Y A numeric matrix of dimension \code{T x N}, where each column is a univariate time series.
#'
#' @returns A numeric matrix of dimension \code{N x L}, where each row is the periodogram of a time series.
#' @export
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(20), ncol = 4)
#' perd <- perd_get(Y)
#'
perd_get <- function(Y) {
  T <- nrow(Y)
  N <- ncol(Y)
  L <- floor(T / 2) - 1
  Y_centered <- apply(Y, 2, function(x)
    x - mean(x))
  fft_result <- mvfft(Y_centered) / sqrt(T)
  periodogram <- abs(fft_result)^2
  perd_matrix <- t(periodogram[2:(L + 1), ])
  return(perd_matrix)
}




##################################################
#' Generate a Fourier Cosine Basis Matrix for Log-Spectral Modeling
#'
#' Constructs a matrix of Fourier cosine basis functions evaluated at a given frequency grid.
#' Used in cepstral smoothing of log-spectra.
#'
#' @param k0 Number of cepstral basis function.
#' @param frq A vector of frequencies in \code{[0,1]}.
#'
#' @returns A \code{k0 x length(frq)} matrix of basis function.
#' @export
#'
#' @examples
#' set.seed(123)
#' frq<-seq(0,1, length.out=5)
#' psi<-psi_get(k0=3, frq)

psi_get <- function(k0, frq) {
  psi = matrix(1, k0, length(frq))
  if (k0 >= 2) {
    for (j in 2:k0) {
      psi[j, ] = sqrt(2) * cos(2 * pi * (j - 1) * frq)
    }
  }
  return(psi)
}


##################################################
#' Estimate Cepstral Coefficients from Periodogram
#'
#' Estimates replicate-specific cepstral coefficients and smoothed log-spectra
#' using a Fourier cosine basis and Whittle-type approximation.
#'
#' @param perd An matrix of periodogram.
#' @param k0 Number of cepstral coefficients.
#' @param frq A vector of frequencies in \code{[0,1]}.
#'
#' @returns A list with:
#' \describe{
#'   \item{\code{f}}{An \code{N × k0} matrix of estimated cepstral coefficients.}
#'   \item{\code{ff}}{An \code{N × K} matrix of smoothed log-spectra.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
#' len <- nrow(Y)
#' L <- floor(len/2)-1
#' frq <- (1:L)/len
#' perd <- perd_get(Y)
#' result <- cep_get(perd = perd, k0 = 3, frq = frq)

cep_get = function(perd,k0,frq){
  # get basis function
  psi = psi_get(k0,frq)
  # Estimate the terms from the penalized sum-of-squares
  Wmatin = matrix(0,k0,k0)
  for(j in 1:length(frq)){
    Wmatin = Wmatin + psi[,j]%*%t(psi[,j])
  }
  Wmat = solve(Wmatin)
  tmp = spec_regress(perd, psi, Wmat, k0)
  f = tmp$f
  ff = tmp$ff
  return(list(f=f,ff=ff))
}


##################################################
#' Fisher Scoring Algorithm For Estimating Cepstral Coefficients
#'
#' Estimates replicate-specific cepstral coefficients and corresponding smoothed log-spectra
#' using a Whittle likelihood approximation.
#'
#' @param perd An N x K matrix of periodogram.
#' @param psi A matrix of cepstral basis functions of dimension \code{k0 × K}.
#' @param Wmat The inverse Gram matrix of the basis functions.
#' @param k0 Number of cepstral basis function
#'
#' @returns A list with:
#' \describe{
#'   \item{\code{f}}{An \code{N × k0} matrix of estimated cepstral coefficients.}
#'   \item{\code{ff}}{An \code{N × K} matrix of smoothed log-spectra.}
#'   }
#' @export
#'
#' @examples
#' set.seed(123)
#' N <- 5
#' len <- 20
#' L <- floor(len/2) - 1
#' frq <- (1:L) / len
#'
#' Y <- matrix(rnorm(len * N), nrow = len, ncol = N)
#'
#' perd <- perd_get(Y)
#'
#' k0 <- 3
#' psi <- psi_get(k0, frq)
#'
#' Wmatin <- matrix(0, k0, k0)
#' for (j in 1:ncol(psi)) {
#'   Wmatin <- Wmatin + psi[, j] %*% t(psi[, j])
#' }
#' Wmat <- solve(Wmatin)
#'
#' out <- spec_regress(perd, psi, Wmat, k0)


spec_regress = function(perd, psi, Wmat, k0) {
  curr_tol = 1
  max_it = 100
  fish_run = 1
  dimen = dim(perd)
  N = dimen[1]
  K = dimen[2]
  # initialize
  fold = cbind(log(rowMeans(perd)), matrix(0, N, k0 - 1))
  Fold = fold %*% psi
  # run iterations
  while (curr_tol > 0.001 && fish_run < max_it) {
    fwork = psi %*% (matrix(1, K, N) - t(perd) * exp(-t(Fold)))
    f = fold - t(fwork) %*% Wmat
    ff = f %*% psi
    curr_tol = mean((f - fold)^2) / mean(fold^2)
    fish_run = fish_run + 1
    fold = f
    Fold = ff
  }
  return(list(f = f, ff = ff))
}



##################################################
#The Reduce-Rank Regression
#' Reduced-Rank Regression on Cepstral Coefficients
#'
#' Fits a reduced-rank regression (RRR) between covariates and cepstral coefficients
#' using a specified maximum rank, and reconstructs log-spectra.
#'
#'
#' @param X A numeric matrix of predictors (N x P).
#' @param f A numeric matrix of cepstral coefficients.
#' @param frq A vector of frequencies in \code{[0,1]}.
#' @param nbase Number of Fourier basis functions.
#' @param nrank Fixed Rank for the reduced-rank regression.
#'
#' @returns A list containing:
#' \describe{
#'   \item{\code{alph}}{Estimated intercept vector.}
#'   \item{\code{bet}}{Estimated coefficient matrix.}
#'   \item{\code{spechat}}{Estimated log-spectra.}
#'   \item{\code{res}}{Matrix of residuals.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' frq <- seq(0, 1, length.out = 16)[2:8]
#' n <- 5
#' p <- 2
#' nbase <- 2
#'
#' X <- matrix(rnorm(n * p), n, p)
#'
#' psi <- psi_get(nbase, frq)
#'
#' true_beta <- matrix(rnorm(p * nbase), p, nbase)
#' alph <- rnorm(nbase)
#' f <- X %*% true_beta + matrix(alph, n, nbase, byrow = TRUE) +
#'      matrix(rnorm(n * nbase), n, nbase)
#'
#' rrr <- rrr_get(X, f, frq, nbase = nbase, nrank = 1)

rrr_get <- function(X, f, frq, nbase, nrank) {
  # get basis function
  psi = psi_get(nbase, frq)
  rr_results = rrr(f, scale(X, scale = FALSE), maxrank = nrank)
  bet = rr_results$coef
  alph = colMeans(f) - t(bet) %*% colMeans(X)
  res = f - X %*% (bet) - matrix(alph, nrow(X), nbase, byrow = TRUE)
  fit = t(alph %*% matrix(1, 1, dim(X)[1])) + X %*% bet
  spechat =  fit %*% psi
  return(list(
    alph = alph,
    bet = bet,
    spechat = spechat,
    res = res
  ))
}


##################################################
#The Envelope Estimator

#' Envelope Estimator for Log-Spectral Regression
#'
#' Fits an envelope regression model to predict cepstral coefficients from covariates.
#'
#' @param X A numeric matrix of predictors (N × p).
#' @param f A numeric matrix of cepstral coefficients (N × nbase).
#' @param frq Numeric vector of frequencies in \code{[0,1]}.
#' @param nbase Number of Fourier basis functions.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{alph}}{Intercept vector.}
#'   \item{\code{bet}}{Envelope regression coefficient matrix.}
#'   \item{\code{spechat}}{Estimated smoothed log-spectra.}
#'   \item{\code{res}}{Residuals from envelope model.}
#' }
#' @export
#'
#' @examples
#' library(Renvlp)
#'
#' set.seed(123)
#' frq <- seq(0, 1, length.out = 16)[2:8]
#' n <- 20
#' p <- 3
#' nbase <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' f <- matrix(rnorm(n * nbase), n, nbase)
#'
#' u_max <- min(ncol(X), ncol(f))
#' cv_errors <- numeric(u_max)
#' for (j in 1:u_max) {
#'   cv_errors[j] <- cv.xenv(X, f, j, m = 5, nperm = 50)
#' }
#' optimal_u <- which.min(cv_errors)
#'
#' env_result <- env_get(X, f, frq, nbase = nbase)


env_get <- function(X, f, frq, nbase) {
  # get basis function
  psi = psi_get(nbase, frq)
  ee = c()
  for (j in 1:min(ncol(f), ncol(X))) {
    ee[j] = cv.xenv(X, f, j, m = 5, nperm = 50)
  }
  u = which(ee == min(ee))
  env_results = xenv(X, f, u)
  bet = env_results$beta
  alph = colMeans(f) - t(bet) %*% colMeans(X)
  res = f - X %*% (bet) - matrix(alph, nrow(X), nbase, byrow = TRUE)
  fit = t(alph %*% matrix(1, 1, dim(X)[1])) + X %*% bet
  spechat =  fit %*% psi
  return(list(
    alph = alph,
    bet = bet,
    spechat = spechat,
    res = res
  ))
}



##################################################
#The Ordinary Least Square Estimator
#' Ordinary Least Squares Estimator for Log-Spectral Regression
#'
#' Performs OLS regression to estimate the association between covariates
#' and cepstral coefficients.
#'
#' @param X A numeric matrix of predictors (N x P).
#' @param f A numeric matrix of cepstral coefficients.
#' @param frq A vector of frequencies in \code{[0,1]}.
#' @param nbase Number of Fourier basis functions.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{alph}}{Intercept vector.}
#'   \item{\code{bet}}{OLS coefficient matrix.}
#'   \item{\code{spechat}}{Estimated smoothed log-spectra.}
#'   \item{\code{res}}{Matrix of residuals.}
#' }
#'
#' @export
#'
#' @examples
#' frq <- seq(0, 1, length.out = 16)[2:8]
#' n <- 10
#' p <- 3
#' nbase <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' f <- matrix(rnorm(n * nbase), n, nbase)
#'
#' ols_result <- ols_get(X, f, frq, nbase)
#'


ols_get = function(X, f, frq, nbase) {
  # get basis function
  psi = psi_get(nbase, frq)
  ols_results = xenv(X, f, dim(X)[2])
  bet = ols_results$beta
  alph = colMeans(f) - t(bet) %*% colMeans(X)
  res = f - X %*% (bet) - matrix(alph, nrow(X), nbase, byrow = TRUE)
  fit = t(alph %*% matrix(1, 1, dim(X)[1])) + X %*% bet
  spechat =  fit %*% psi
  return(list(
    alph = alph,
    bet = bet,
    spechat = spechat,
    res = res
  ))
}



##################################################
#' Compute Functional Effects of Intercept and Covariates
#'
#' Projects cepstral coefficient intercept and covariate effects onto the frequency domain
#' using the cepstral basis functions.
#'
#' @param alpha A numeric vector of cepstral intercept coefficients.
#' @param beta A numeric matrix of regression coefficients.
#' @param frq Numeric vector of frequency points in \code{[0,1]}.
#' @param nbase Number of Fourier basis functions.
#' @param ind An integer vector indicating the indices of covariates to be included in the model.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{alpha_effect}}{Functional intercept across frequency.}
#'   \item{\code{beta_effect}}{Matrix of functional covariate effects.}
#' }
#' @export
#'
#' @examples
#' frq <- seq(0, 1, length.out = 16)[2:8]
#' alpha <- rnorm(3)
#' beta <- matrix(rnorm(2 * 3), 2, 3)
#' result <- effect_get(alpha, beta, frq, nbase = 3, ind = c(1, 2))
#'
effect_get = function(alpha, beta, frq, nbase, ind) {
  # get basis function
  psi = psi_get(nbase, frq)
  alpha_effect = as.vector(t(alpha) %*% psi)

  beta_effect = matrix(0, length(frq), length(ind))
  for (j in 1:length(ind)) {
    beta_effect[, j] = t(beta[j, ]) %*% psi
  }
  return(list(alpha_effect = alpha_effect, beta_effect = beta_effect))
}



##################################################
#Bootstrap Confidence Interval for Effect Function
#' Bootstrap Confidence Intervals for Functional Effect Curves
#'
#' Computes bias-corrected percentile bootstrap confidence intervals for the intercept
#' and covariate effect functions in cepstral-based regression models.
#'
#' @param logspect Matrix of estimated log-spectra.
#' @param res Matrix of residuals from cepstral regression.
#' @param alpha_effect Vector of estimated intercept effect function.
#' @param beta_effect Matrix of estimated covariate effect functions.
#' @param X Covariate matrix.
#' @param nbase Number of cepstral basis functions.
#' @param frq1 Frequency grid used for cepstral modeling.
#' @param frq2 Frequency grid used for reconstructing spectra.
#' @param nrank Rank for reduced-rank regression.
#' @param ind A vector of indices indicating which covariates to compute effect confidence intervals for.
#' @param level Confidence level.
#' @param nboot Number of bootstrap iterations.
#' @param method Regression method: "rrr", "ols", or "env".
#' @param verb Logical; if TRUE, prints bootstrap iteration number.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{alpha_ci}}{Matrix with lower and upper CI for intercept effect.}
#'   \item{\code{beta_ci}}{Array or matrix of lower and upper CI for covariate effects.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' N <- 10
#' len <- 12
#' nbase <- 2
#' nrank <- 1
#' nboot <- 10
#' level <- 0.95
#' p <- 2
#' ind <- 1:p
#'
#' Y <- matrix(rnorm(N * nbase), nrow = N, ncol = nbase)
#'
#' X <- matrix(rnorm(N * p), nrow = N, ncol = p)
#'
#'
#' frq <- seq(1, nbase) / len
#'
#' rrr_out <- rrr_get(Y, X, frq, nbase, nrank)
#'
#' eff <- effect_get(rrr_out$alph, rrr_out$bet, frq, nbase, ind)
#' alpha_eff <- eff$alpha_effect
#' beta_eff <- eff$beta_effect
#'
#' logspect <- matrix(rnorm(N * length(frq)), nrow = N, ncol = length(frq))
#'
#' boot_ci <- boot_effect(
#'   logspect = logspect,
#'   res = rrr_out$res,
#'   alpha_effect = alpha_eff,
#'   beta_effect = beta_eff,
#'   X = X,
#'   nbase = nbase,
#'   frq1 = frq,
#'   frq2 = frq,
#'   nrank = nrank,
#'   ind = ind,
#'   level = level,
#'   nboot = nboot,
#'   method = "rrr",
#'   verb = TRUE
#' )
#'
#' plot(frq, beta_eff[, 1], type = "l", col = "blue", lwd = 2,
#'      ylab = "Effect", xlab = "Frequency",
#'      main = paste("Effect Function and", level*100, "% CI for Covariate", ind[1]))
#' lines(frq, boot_ci[[ind[1]]][, 1], col = "red", lty = 2)
#' lines(frq, boot_ci[[ind[1]]][, 2], col = "red", lty = 2)
#' legend("topright", legend = c("Effect", "Bootstrap CI"), col = c("blue", "red"),
#'        lty = c(1, 2), lwd = c(2, 1))

boot_effect <- function(logspect,
                        res,
                        alpha_effect,
                        beta_effect,
                        X,
                        nbase,
                        frq1,
                        frq2,
                        nrank,
                        ind,
                        level,
                        nboot,
                        method = "rrr",
                        verb = FALSE) {
  len = ncol(logspect)
  N = nrow(X)
  psi = psi_get(nbase, frq2)
  Zboot = matrix(0, len, N)
  alpha_effect_boot = matrix(0, length(frq1), nboot)
  beta_effect_boot = array(0, dim = c(length(frq1), length(ind), nboot))
  for (boot in 1:nboot) {
    if (verb == TRUE) {
      print(boot)
    }
    indx = sample(1:N, N, replace = TRUE)
    res_boot <- res[indx, ]
    res_fit =  res_boot %*% psi
    specthat = exp(logspect + res_fit)
    Zboot = data_generater(N, len, t(sqrt(specthat)))
    perd = perd_get(Zboot)
    f = cep_get(perd, nbase, frq1)$f
    if (method == "rrr") {
      output = rrr_get(X, f, frq1, nbase, nrank)
    } else if (method == "ols") {
      output = ols_get(X, f, frq1, nbase)
    } else if (method == "env") {
      output = env_get(X, f, frq1, nbase)
    }
    alph = output$alph
    bet = output$bet

    effect_output = effect_get(alph, bet, frq1, nbase, ind)
    alpha_effect_boot[, boot] =  effect_output$alpha_effect
    beta_effect_boot[, , boot] = effect_output$beta_effect
  }


  if (length(ind) == 1) {
    bias_alpha = rowMeans(alpha_effect_boot) - alpha_effect

    bias_beta =  rowMeans(beta_effect_boot) - beta_effect

    bias_corrected_alpha = alpha_effect_boot - bias_alpha

    bias_corrected_beta = beta_effect_boot - array(bias_beta, dim = dim(beta_effect_boot))

    alpha_ci = cbind(
      apply(bias_corrected_alpha, 1, quantile, probs = (1 - level) / 2),
      apply(bias_corrected_alpha, 1, quantile, probs = 1 -
              (1 - level) / 2)
    )
    beta_ci = cbind(
      apply(bias_corrected_beta, 1, quantile, probs = (1 - level) / 2),
      apply(bias_corrected_beta, 1, quantile, probs = 1 -
              (1 - level) / 2)
    )
  } else{
    bias_alpha = rowMeans(alpha_effect_boot) - alpha_effect

    bias_beta = apply(beta_effect_boot, c(1, 2), mean) - beta_effect

    bias_corrected_alpha = alpha_effect_boot - bias_alpha

    bias_corrected_beta = beta_effect_boot - array(bias_beta, dim = dim(beta_effect_boot))

    alpha_ci = cbind(
      apply(bias_corrected_alpha, 1, quantile, probs = (1 - level) / 2),
      apply(bias_corrected_alpha, 1, quantile, probs = 1 -
              (1 - level) / 2)
    )
    beta_ci = cbind(
      apply(bias_corrected_beta, c(1, 2), quantile, probs = (1 - level) / 2),
      apply(bias_corrected_beta, c(1, 2), quantile, probs =
              1 - (1 - level) / 2)
    )
  }


  return(list(alpha_ci = alpha_ci, beta_ci = beta_ci))
}



##################################################
#' Generate Exponential Correlation Covariance Matrix
#'
#' Creates an n × n covariance matrix with entries \eqn{\rho^{|i - j|}}.
#'
#' @param n Dimension of the covariance matrix.
#' @param rho Correlation decay parameter.
#'
#' @returns An n × n positive definite covariance matrix.
#' @export
#'
#' @examples
#' S <- generate_sig(5, 0.5)
#'
#'
generate_sig = function(n, rho) {
  if (rho == 0) {
    sig = diag(n)
  } else{
    a = seq(1:n)
    b = matrix(a, n, n)
    c = t(b)
    q = abs(b - c)
    tmp = q * log(rho)
    sig = exp(tmp)
  }
  return(sig)
}


##################################################
#Generate data from spectral representation
#' Generate Time Series
#'
#' Simulates real-valued time series using the Cramér spectral representation and inverse FFT.
#'
#' @param N Number of time series to generate.
#' @param nobs Number of time points.
#' @param spec Spetral density matrix.
#'
#' @returns Matrix of size \code{nobs × N} of generated time series
#' @export
#'
#' @examples
#' set.seed(123)
#' N    <- 3
#' nobs <- 20
#' freqs <- (1:nobs) / nobs
#'
#' spec <- matrix(NA, nrow = nobs, ncol = N)
#' for (i in 1:N) {
#'   spec[, i] <- exp(2 * cos(2 * pi * freqs) + rnorm(1, sd = 0.1))
#' }
#'
#' data_generater(N = N, nobs = nobs, spec = spec)
#'
data_generater <- function(N, nobs, spec) {
  x <- matrix(0, nrow = nobs, ncol = N)
  for (i in 1:N) {
    Rez <- rnorm((nobs / 2) - 1, 0, 1 / sqrt(nobs))
    Imz <- rnorm((nobs / 2) - 1, 0, 1 / sqrt(nobs))
    z <- Rez + (1i) * Imz
    z[nobs / 2] = rnorm(1, 0, 1 / sqrt(nobs))
    z[nobs] = rnorm(1, 0, 1 / sqrt(nobs))
    z[(nobs / 2 + 1):(nobs - 1)] <- Conj(z[(nobs / 2 - 1):1])
    for (j in 1:nobs) {
      tmp <- exp(1i * 2 * pi * 1 * j / nobs) * spec[1, i] %*% z[1]
      for (k in 2:nobs) {
        tmp <- tmp + exp(1i * 2 * pi * k * j / nobs) * spec[k, i] %*% z[k]
      }
      x[j, i] <- Re(tmp)
    }
  }
  return(x)
}


##################################################
#' Cepstral Regression
#'
#' Performs cepstral regression to model frequency domain relationships between a
#' functional response and scalar covariates. Supports ordinary least squares (OLS),
#' reduced-rank regression (RRR), and envelope regression (ENV) methods. Automatically
#' selects the number of cepstral basis functions via AIC.
#'
#' @param y Numeric matrix of dimension (time points) × (samples).
#' @param x Numeric matrix of scalar covariates with dimensions (samples) × (covariates).
#' @param method One of "ols", "rrr", or "env" specifying the regression method.
#' @param number_of_K Maximum number of cepstral basis functions to consider for AIC selection.
#' @param if_bootstrap Logical; whether to compute bootstrap confidence intervals (default FALSE).
#' @param level Confidence level for bootstrap intervals. Required if `if_bootstrap = TRUE`.
#' @param nboot Integer; the number of bootstrap samples. Required if `if_bootstrap = TRUE`.
#' @param ind Integer vector; indices of covariates for which the effect functions are to be estimated and plotted.
#'            Required if `if_bootstrap = TRUE`.
#' @param nrank Integer; the rank used for reduced-rank regression. Required when `method = "rrr"` or when bootstrapping with `"rrr"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{eff}{A list of estimated effect functions (e.g., \code{alpha_effect}, \code{beta_effect}).}
#'   \item{boot}{A list of bootstrap results including confidence intervals; `NULL` if `if_bootstrap = FALSE`.}
#'   \item{fit}{A list containing regression coefficients, residuals, smoothed spectral estimates, and other model outputs.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' niter <- 5
#' len <- 10
#' N <- 3
#' p <- 2
#' L <- floor(len/2)-1
#' frq <- (1:L)/len
#' mu <- rep(0, p)
#' rho <- 0
#' Sigma <- generate_sig(p, rho)
#'
#' X <- MASS::mvrnorm(N, mu, Sigma)
#' X[,1] <- runif(N, 0, 1)
#'
#' spec <- matrix(0,len,N)
#' for(j in 1:N){
#'   eta1 <- rnorm(1,0,0.5)
#'   eta2 <- rnorm(1,0,0.5)
#'   eta3 <- rnorm(1,0,0.5)
#'   spec[,j] <- exp(
#'     2*cos(2*pi*(1:len)/len) +
#'     X[j,1]*(2*cos(4*pi*(1:len)/len)) +
#'     eta1 + eta2*cos(2*pi*(1:len)/len) +
#'     eta3*(cos(4*pi*(1:len)/len))
#'     )
#'  }
#'
#' Z <- data_generater(N,len,sqrt(spec))
#'
#' res_ols <- CepReg(Z, X, method = "ols", number_of_K = 2,
#'          if_bootstrap = TRUE, level = 0.95,
#'          nboot = 2, ind = 1)
#'
#' eff_ols <- res_ols$eff
#' boot_ols <- res_ols$boot
#'
#' plot(frq, eff_ols$alpha_effect, type = 'l', col = "black", xlab = "Frequency", ylab = "",
#'      ylim = range(c(boot_ols$alpha_ci,
#'      eff_ols$alpha_effect, 2*cos(2*pi*frq)+0.577)))
#' title(ylab = expression(alpha(omega)), line = 2, cex.lab = 1.2)
#' lines(frq, boot_ols$alpha_ci[, 1], col = "black")
#' lines(frq, boot_ols$alpha_ci[, 2], col = "black")
#' lines(frq, 2 * cos(2 * pi * frq) + 0.577, col = "red", lty = 1, lwd = 2)
#' legend("topright",legend = c("True", "Estimated", "CI"),
#'      col = c("red", "black", "black"), ncol = 1)
#'
CepReg <- function(y,
                   x,
                   method = c("ols", "rrr", "env"),
                   number_of_K,
                   if_bootstrap = FALSE,
                   level = NULL,
                   nboot = NULL,
                   ind = NULL,
                   nrank = NULL) {
  method <- match.arg(method)
  n <- nrow(x)
  p <- ncol(x)

  if (if_bootstrap) {
    if (is.null(level) || is.null(nboot) || is.null(ind)) {
      stop("When bootstrap=TRUE, you must provide 'level', 'nboot', and 'ind'.")
    }
    if (method == "rrr" && is.null(nrank)) {
      stop("When bootstrap=TRUE and method='rrr', 'nrank' must be provided.")
    }
  }

  # Prepare frequency and basis
  len <- nrow(y)
  L <- floor(len / 2) - 1
  frq <- (1:L) / len
  frq2 <- (1:len) / len

  # Compute periodogram
  perd <- perd_get(y)

  # Select parameter nbase by AIC
  Ks <- number_of_K
  AICs <- numeric(Ks)
  for (j in 1:Ks) {
    tmp <- cep_get(perd, j, frq)
    AICs[j] <- sum(tmp$ff + perd * exp(-tmp$ff)) + 2 * n * j
  }
  nbase <- which.min(AICs)

  # Estimate cepstral coefficients
  cep_obj <- cep_get(perd, nbase, frq)
  f <- cep_obj$f

  # Fit regression based on method
  if (method == "ols") {
    fit <- ols_get(as.matrix(x), f, frq2, nbase)
  } else if (method == "rrr") {
    if (is.null(nrank))
      stop("Please specify nrank for RRR method")
    fit <- rrr_get(as.matrix(x), f, frq2, nbase, nrank)
  } else if (method == "env") {
    fit <- env_get(as.matrix(x), f, frq2, nbase)
  }

  # Compute effect functions
  if (is.null(ind))
    ind <- 1:p

  eff <- effect_get(fit$alph, fit$bet, frq, nbase, ind)

  # Bootstrap CI if requested
  if (if_bootstrap) {
    boot_res <- boot_effect(
      fit$spechat,
      fit$res - mean(fit$res),
      eff$alpha_effect,
      eff$beta_effect,
      as.matrix(x),
      nbase,
      frq,
      frq2,
      nrank,
      ind,
      level,
      nboot,
      method = method,
      verb = FALSE
    )
  } else {
    boot_res <- NULL
  }

  # Return results
  list(eff = eff, boot = boot_res, fit = fit)
}

