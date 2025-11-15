##############################################
# ------------ HELPER FUNCTIONS ------------ #
##############################################
### TODO: move helper functions into their own file

#' Extracts row-wise dot products
#'
#' @param mat numeric matrix; represents covariates/features for each individual
#' @param coef numeric matrix; represents coefficients corresponding to each of the covariates in `mat`
#' @param n integer; sample size
#'
#' @return Returns a numeric vector of linear predictors, one per individual. Returns a vector of zeroes if there are no predictors.
#'
#' @keywords internal
.get_linear_predictor <- function(mat, coef, n) {
  if (is.null(mat) || length(coef) == 0L) {
    return(rep(0, n))
  }
  else{
    drop(mat %*% coef)
  }
}

#' Logistic function as used in confounded treatment assignment
#'
#' @param x numeric vector
#'
#' @return Element-wise logistic function.
#'
#' @keywords internal
.logistic_function <- function(x){
  return(1 / (1 + exp(- x )))
}

#' Gives the inverse generator function for the Gumbel-Hougaard Copula with parameter `theta`.
#'
#' @param t numeric; dummy variable of the inverse generator \eqn{1/\phi}
#' @param theta numeric \eqn{>=} 1; Gumbel-Hougaard association parameter
#'
#' @return
#' Gives the inverse generator.
#'
#' @details
#' Recreated from the `gumbel` package. See package documentation for details.
#'
#' @keywords internal
.invphigumbel <- function(t, theta = 1) {
  exp(-t^(1/theta))
}

#' Draw bivariate uniforms from a Gumbel--Hougaard copula
#' @param n integer; number of observations
#' @param theta numeric \eqn{>=} 1; Gumbel-Hougaard association parameter
#'
#' @return n x 2 matrix, each row (U1, U2), both in (0,1), with GH copula parameter = theta.
#'
#' @details
#' This mirrors the construction of the `gumbel.jeong` function used in the `cWR` package. For even more details, see the manual for R package `gumbel`.
#'
#' @keywords internal
.draw_gumbel_uniforms <- function(n, theta = 1) {

  # We are only interested in generating a bivariate vector for now
  # Will look to expand to more than two outcomes in the future
  dim <- 2L

  # Generate positive-stable latent term
  beta <- 1 / theta
  unifpirand <- runif(n, 0, pi)
  exprand <- matrix(rexp(dim * n), nrow = n, ncol = dim)
  exprand2   <- rexp(n)
  stablerand <- sin((1 - beta) * unifpirand)^((1 - beta) / beta) *
    (sin(beta * unifpirand)) /
    (sin(unifpirand))^(1 / beta)
  stablerand <- stablerand / (exprand2^(theta - 1))

  # Inverse Gumbel-Hougaard transformation
  # This returns an n x 2 matrix of uniforms with GH(theta) dependence.
  u_mat <- .invphigumbel(exprand / stablerand, theta)

  return(u_mat)
}


############################################
# ------------ MAIN FUNCTIONS ------------ #
############################################

#' @description
#' Generate a single dataset of i.i.d. subjects under semi-competing risks, with:
#' \itemize{
#' \item Measured covariates \code{X} w/ dimension \code{p}.
#' \item Unmeasured covariates \code{U} w/ dimension \code{q}.
#' \item Treatment assignment \code{A}. Randomized or confounded.
#' \item Observed semi-competing outcomes: \code{Y_D}. \code{delta_D}, \code{Y_H}, \code{delta_H}.
#' \item Derived time-to-first-event outcomes \code{Y_tfe}, \code{delta_D}.
#' }
#' Under:
#' \itemize{
#' \item Latent time to fatal and non-fatal event: \code{T_D}, \code{T_H}.
#' \item Gumbel–Hougaard copula dependence with parameter \code{theta_copula}.
#' \item Random censoring time \code{C}.
#' \item Administrative censoring time \code{phi_admin}.
#' }
#'
#' @param N integer; sample size for this data set.
#' @param p integer; number of measured covariates.
#' @param q Integer; number of unmeasured covariates.
#' @param mu_X numeric vector length p; means for measured covariates.
#' @param sd_X numeric vector length p; SDs for measured covariates.
#' @param mu_U numeric vector length q; means for unmeasured covariates.
#' @param sd_U numeric vector length q; SDs for unmeasured covariates.
#' @param treat_assign character; one of "randomized" or "confounded".
#' @param alpha_0 numeric; intercept in treatment model if confounded.
#' @param alpha_X numeric vector length p; coefficients for X in treatment model.
#' @param alpha_U numeric vector length q; coefficients for U in treatment model.
#' @param lambda_D numeric; baseline hazard rate for fatal outcome.
#' @param kappa_D numeric; Weibull shape parameter for marginal time-to-fatal outcome. Becomes Exponentially distributed when equal to 1.
#' @param beta_A_D numeric; treatment log-HR effect on fatal outcome.
#' @param beta_X_D numeric vector length p; measured covariate effects on time-to-fatal outcome.
#' @param beta_U_D numeric vector length q; unmeasured covariate effects on time-to-fatal outcome.
#' @param lambda_H numeric; baseline hazard rate for time-to-non-fatal outcome.
#' @param kappa_H numeric; Weibull shape for time-to-non-fatal outcome. Becomes Exponentially distributed when equal to 1.
#' @param beta_A_H numeric; Treatment log-HR effect on time-to-non-fatal outcome.
#' @param beta_X_H numeric vector length p; Measured covariate effects on time-to-non-fatal outcome.
#' @param beta_U_H numeric vector length q; Unmeasured covariate effects on time-to-non-fatal outcome.
#' @param theta_copula numeric \eqn{>=} 1; Gumbel–Hougaard dependence parameter. No association between latent competing event times when equal to 1.
#' @param lambda_C numeric; baseline censoring rate.
#' @param beta_A_C numeric; treatment effect on censoring rate.
#' @param phi_admin numeric; administrative censoring time.
#' @param incl_latent boolean; default FALSE; additionally includes latent event times in the output dataset.
#' @param incl_tfe boolean: default FALSE; additionally includes calculated time-to-first-event variables in the output dataset.
#'
#' @return
#' A \code{data.frame} object with columns:
#' \itemize{
#' \item \code{id}: unique subject identifier
#' \item \code{A}: treatment assignment
#' \item \code{Y_D}: observed time to fatal event
#' \item \code{delta_D}: oberved fatal event indicator
#' \item \code{Y_H}: observed time to non-fatal event
#' \item \code{delta_H}: observed non-fatal event indicator
#' \item \code{X_1}, ..., \code{X_p}: measured covariates
#' \item \code{U_1}, ..., \code{U_q}: unmeasured covariates
#' \item \code{Y_tfe}: observed time to first event (if \code{incl_tfe} is set to \code{TRUE})
#' \item \code{delta_tfe}: observed time to first event indicator (if \code{incl_tfe} is set to \code{TRUE})
#' \item\code{T_D}: latent time to fatal event (if \code{incl_latent} is set to \code{TRUE})
#' \item\code{T_H}: latent time to non-fatal event (if \code{incl_latent} is set to \code{TRUE})
#' }
#'
#' @details
#' Weibull PH form:
#'   S_D(t|A,X,U) = exp( - (Lambda_D * t)^kappa_D ),
#'   Lambda_D = lambda_D * exp(beta_A_D*A + beta_X_D^T X + beta_U_D^T U).
#'
#' Event time simulation via inverse transform:
#'   T_D = ( -log(V_D) / Lambda_D )^(1/kappa_D ),
#'   given V_D ~ U(0,1).
#'
#' Joint dependence across (T_H, T_D) is induced by sampling
#' a bivariate (V_H, V_D) from a Gumbel–Hougaard copula with parameter
#' theta_copula, then applying the inverse Weibull transformations
#' separately for fatal and non-fatal events.
#'
#' @export
simulate_dataset <- function(
    N,
    p, q,
    mu_X, sd_X,
    mu_U, sd_U,
    treat_assign = c("randomized", "confounded"),
    alpha_0 = 0,
    alpha_X = NULL,
    alpha_U = NULL,
    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
    theta_copula,
    lambda_C, beta_A_C,
    phi_admin,
    incl_tfe = FALSE,
    incl_latent = FALSE
) {

  # Validate treatment assignment argument
  treat_assign <- match.arg(treat_assign)

  #---------------------------
  # 1. GENERATE COVARIATES
  #---------------------------

  # Generate measured covariate matrix X (N x p)
  if (p > 0L) {
    X_mat <- sapply(1:p, function(j) {
      rnorm(N, mean = mu_X[j], sd = sd_X[j])
    })
    colnames(X_mat) <- paste0("X_", 1:p)
  } else {
    X_mat <- NULL
  }

  # Generate unmeasured covariate matrix U (N x q)
  if (q > 0L) {
    U_mat <- sapply(1:q, function(j) {
      rnorm(N, mean = mu_U[j], sd = sd_U[j])
    })
    colnames(U_mat) <- paste0("U_", 1:q)
  } else {
    U_mat <- NULL
  }

  #---------------------------
  # 2. TREATMENT ASSIGNMENT
  #---------------------------

  # Generate treatment vector A (N x 1)
  if (treat_assign == "randomized") {
    # Simulate randomized "coin-flip" assignment
    A <- rbinom(N, size = 1, prob = 0.5)
  } else if (treat_assign == "confounded") {
    # Simulate "confounded" assignment via logit model with X and U
    lp_treat <- alpha_0 +
      .get_linear_predictor(X_mat, alpha_X, N) +
      .get_linear_predictor(U_mat, alpha_U, N)
    A <- rbinom(N, size = 1, prob = .logistic_function(lp_treat))
  }

  #---------------------------
  # 3. COPULA
  #---------------------------

  # Sample (V_H, V_D) ~ GH(theta_copula), where V_H,V_D are marginal Uniform(0,1)
  V_mat <- .draw_gumbel_uniforms(N, theta_copula)
  V_H <- V_mat[, 1]
  V_D <- V_mat[, 2]

  #---------------------------
  # 4. WEIBULL PARAMETERS
  #---------------------------

  # Fatal Weibull scale parameter
  Lambda_D <- lambda_D *
    exp(beta_A_D * A +
          .get_linear_predictor(X_mat, beta_X_D, N) +
          .get_linear_predictor(U_mat, beta_U_D, N))
  # Non-fatal Weibull scale parameter
  Lambda_H <- lambda_H *
    exp(beta_A_H * A +
          .get_linear_predictor(X_mat, beta_X_H, N) +
          .get_linear_predictor(U_mat, beta_U_H, N))

  #---------------------------
  # 5. DERIVE EVENT TIMES
  #---------------------------

  # Latent fatal event time
  T_D <- ( -log(V_D) / Lambda_D )^(1 / kappa_D)
  # Latent non-fatal event time
  T_H <- ( -log(V_H) / Lambda_H )^(1 / kappa_H)

  #---------------------------
  # 6. GENERATE CENSORING
  #---------------------------

  # If lambda_C is non-positive...
  if (isTRUE(lambda_C <= 0L)) {
    # Interpret as "no random censoring" (C is infinite)...
    C <- rep.int(Inf, N)
  } else {
    # Otherwise, calculate censoring rate and censoring variable...
    Lambda_C <- lambda_C * exp(-beta_A_C * A)
    # ... keeping Lambda_C positive in case non-positive rates are accidentally generated
    C <- rexp(N, rate = pmax(Lambda_C, .Machine$double.eps))
  }

  #---------------------------
  # 7. OBSERVED EVENT DATA
  #---------------------------

  # Observed fatal event data
  Y_D <- pmin(T_D, C, phi_admin)
  delta_D <- as.integer(T_D <= pmin(C, phi_admin))
  # Observed non-fatal event data:
  Y_H <- pmin(T_H, T_D, C, phi_admin)
  delta_H <- as.integer(
    (T_H <= T_D) & (T_H <= C) & (T_H <= phi_admin)
  )

  #---------------------------
  # 8. ASSEMBLE DATA
  #---------------------------

  # Create output data-frame
  output_df <- data.frame(
    id = seq_len(N),
    A = A,
    Y_D = Y_D,
    delta_D = delta_D,
    Y_H = Y_H,
    delta_H = delta_H
  )

  # Attach measured covariates
  if (p > 0L) {
    output_df <- cbind(output_df, as.data.frame(X_mat, optional = TRUE))
  }
  # Attach unmeasured covariates
  if (q > 0L) {
    output_df <- cbind(output_df, as.data.frame(U_mat, optional = TRUE))
  }
  # Attach time-to-first-event variables if incl_tfe is set to TRUE
  if (incl_tfe == TRUE){
    # Calculate time-to-first-event endpoint
    Y_tfe <- pmin(T_H, T_D, C, phi_admin)
    delta_tfe <- as.integer(pmin(T_H, T_D) <= pmin(C, phi_admin))
    # Attach it to the output
    output_df <- cbind(output_df, Y_tfe, delta_tfe)
  }

  # Attach latent event times if incl_latent is set to TRUE
  if (incl_latent == TRUE){
    output_df <- cbind(output_df, T_D, T_H)
  }

  rownames(output_df) <- NULL

  return(output_df)
}
