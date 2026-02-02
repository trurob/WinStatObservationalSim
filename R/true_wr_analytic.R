# True Causal Win Ratio - Analytical Calculation
# Based on Zhang et al. (2019, 2022), Luo et al. (2015)
#
# This script implements the analytical calculation of the true causal win ratio
# for observational studies with semi-competing risks data.
#
# This package is a modified versiion of Zhang's cWR package:
# https://github.com/dee1008/cWR

##############################################
# ------------ HELPER FUNCTIONS ------------ #
##############################################

#' Marginal survival function for fatal event under treatment a
#'
#' Computes S_D(t | a, x, u) for Weibull marginals.
#'
#' @param t numeric scalar or vector; time
#' @param a numeric scalar (0 or 1); treatment indicator
#' @param x numeric vector; measured covariates (can be NULL)
#' @param u numeric vector; unmeasured covariates (can be NULL)
#' @param lambda_D numeric scalar > 0; baseline hazard rate for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape parameter for fatal event
#' @param beta_A_D numeric scalar; treatment effect on fatal event
#' @param beta_X_D numeric vector; effects of measured covariates on fatal event (can be NULL)
#' @param beta_U_D numeric vector; effects of unmeasured covariates on fatal event (can be NULL)
#'
#' @return numeric; survival probability S_D(t | a, x, u)
#'
#' @keywords internal
.S_D <- function(t, a, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D) {
  # Calculate linear predictor
  lp <- beta_A_D * a
  if (!is.null(x) && !is.null(beta_X_D)) {
    lp <- lp + sum(beta_X_D * x)
  }
  if (!is.null(u) && !is.null(beta_U_D)) {
    lp <- lp + sum(beta_U_D * u)
  }

  # Weibull scale parameter
  Lambda_D <- lambda_D * exp(-lp)

  # Weibull survival function
  return(exp(-(Lambda_D * t)^kappa_D))
}

#' Marginal survival function for non-fatal event under treatment a
#'
#' Computes S_H(t | a, x, u) for Weibull marginals.
#'
#' @param t numeric scalar or vector; time
#' @param a numeric scalar (0 or 1); treatment indicator
#' @param x numeric vector; measured covariates (can be NULL)
#' @param u numeric vector; unmeasured covariates (can be NULL)
#' @param lambda_H numeric scalar > 0; baseline hazard rate for non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape parameter for non-fatal event
#' @param beta_A_H numeric scalar; treatment effect on non-fatal event
#' @param beta_X_H numeric vector; effects of measured covariates on non-fatal event (can be NULL)
#' @param beta_U_H numeric vector; effects of unmeasured covariates on non-fatal event (can be NULL)
#'
#' @return numeric; survival probability S_H(t | a, x, u)
#'
#' @keywords internal
.S_H <- function(t, a, x, u, lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H) {
  # Calculate linear predictor
  lp <- beta_A_H * a
  if (!is.null(x) && !is.null(beta_X_H)) {
    lp <- lp + sum(beta_X_H * x)
  }
  if (!is.null(u) && !is.null(beta_U_H)) {
    lp <- lp + sum(beta_U_H * u)
  }

  # Weibull scale parameter
  Lambda_H <- lambda_H * exp(-lp)

  # Weibull survival function
  return(exp(-(Lambda_H * t)^kappa_H))
}

#' Joint survival function for (T_H, T_D) with Gumbel-Hougaard copula
#'
#' Computes S_{D,H}(t_D, t_H | a, x, u) using Gumbel-Hougaard copula.
#'
#' @param t_H numeric scalar; time for non-fatal event
#' @param t_D numeric scalar; time for fatal event
#' @param a numeric scalar (0 or 1); treatment indicator
#' @param x numeric vector; measured covariates (can be NULL)
#' @param u numeric vector; unmeasured covariates (can be NULL)
#' @param lambda_H numeric scalar > 0; baseline hazard rate for non-fatal event
#' @param lambda_D numeric scalar > 0; baseline hazard rate for fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape parameter for non-fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape parameter for fatal event
#' @param beta_A_H numeric scalar; treatment effect on non-fatal event
#' @param beta_A_D numeric scalar; treatment effect on fatal event
#' @param beta_X_H numeric vector; effects of measured covariates on non-fatal event (can be NULL)
#' @param beta_X_D numeric vector; effects of measured covariates on fatal event (can be NULL)
#' @param beta_U_H numeric vector; effects of unmeasured covariates on non-fatal event (can be NULL)
#' @param beta_U_D numeric vector; effects of unmeasured covariates on fatal event (can be NULL)
#' @param theta_copula numeric >= 1; Gumbel-Hougaard copula parameter
#'
#' @return numeric; joint survival probability
#'
#' @keywords internal
.S_DH <- function(t_H, t_D, a, x, u,
                  lambda_H, lambda_D, kappa_H, kappa_D,
                  beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                  theta_copula) {

  # Calculate linear predictors
  lp_H <- beta_A_H * a
  lp_D <- beta_A_D * a

  if (!is.null(x)) {
    if (!is.null(beta_X_H)) lp_H <- lp_H + sum(beta_X_H * x)
    if (!is.null(beta_X_D)) lp_D <- lp_D + sum(beta_X_D * x)
  }
  if (!is.null(u)) {
    if (!is.null(beta_U_H)) lp_H <- lp_H + sum(beta_U_H * u)
    if (!is.null(beta_U_D)) lp_D <- lp_D + sum(beta_U_D * u)
  }

  # Weibull scale parameters
  Lambda_H <- lambda_H * exp(-lp_H)
  Lambda_D <- lambda_D * exp(-lp_D)

  # Gumbel-Hougaard copula joint survival
  term_H <- (Lambda_H * t_H)^(kappa_H * theta_copula)
  term_D <- (Lambda_D * t_D)^(kappa_D * theta_copula)

  return(exp(-((term_H + term_D)^(1 / theta_copula))))
}

#' Conditional survival G_a(y1, y2 | x, u)
#'
#' Computes the conditional survival of T_H beyond y1 given T_D survives beyond y2.
#' G_a(y1, y2 | x, u) = S_{D,H}(y2, y1 | a, x, u) / S_D(y2 | a, x, u)
#'
#' @inheritParams .S_DH
#' @param y1 numeric scalar; time for non-fatal event
#' @param y2 numeric scalar; time for fatal event
#'
#' @return numeric; conditional survival probability
#'
#' @keywords internal
.G <- function(y1, y2, a, x, u,
               lambda_H, lambda_D, kappa_H, kappa_D,
               beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
               theta_copula) {

  numerator <- .S_DH(y1, y2, a, x, u,
                     lambda_H, lambda_D, kappa_H, kappa_D,
                     beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                     theta_copula)

  denominator <- .S_D(y2, a, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)

  return(numerator / denominator)
}

#' Hazard function for fatal event at time t
#'
#' Computes λ_D(t | a, x, u) for Weibull distribution.
#'
#' @inheritParams .S_D
#'
#' @return numeric; hazard rate at time t
#'
#' @keywords internal
.lambda_D <- function(t, a, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D) {
  # Calculate linear predictor
  lp <- beta_A_D * a
  if (!is.null(x) && !is.null(beta_X_D)) {
    lp <- lp + sum(beta_X_D * x)
  }
  if (!is.null(u) && !is.null(beta_U_D)) {
    lp <- lp + sum(beta_U_D * u)
  }

  # Weibull scale parameter
  Lambda_D <- lambda_D * exp(-lp)

  # Weibull hazard function: kappa * Lambda^kappa * t^(kappa-1)
  return(kappa_D * (Lambda_D^kappa_D) * (t^(kappa_D - 1)))
}

#' Conditional hazard for non-fatal event given fatal event time
#'
#' Computes λ_{H|D}(y1 | y2, a, x, u), the hazard for T_H at y1 conditional on T_D = y2.
#' For Gumbel copula with Weibull marginals, derived from copula density.
#'
#' @inheritParams .G
#'
#' @return numeric; conditional hazard rate
#'
#' @keywords internal
.lambda_H_given_D <- function(y1, y2, a, x, u,
                              lambda_H, lambda_D, kappa_H, kappa_D,
                              beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                              theta_copula) {

  # Calculate linear predictors
  lp_H <- beta_A_H * a
  lp_D <- beta_A_D * a

  if (!is.null(x)) {
    if (!is.null(beta_X_H)) lp_H <- lp_H + sum(beta_X_H * x)
    if (!is.null(beta_X_D)) lp_D <- lp_D + sum(beta_X_D * x)
  }
  if (!is.null(u)) {
    if (!is.null(beta_U_H)) lp_H <- lp_H + sum(beta_U_H * u)
    if (!is.null(beta_U_D)) lp_D <- lp_D + sum(beta_U_D * u)
  }

  # Weibull scale parameters
  Lambda_H <- lambda_H * exp(-lp_H)
  Lambda_D <- lambda_D * exp(-lp_D)

  # For Gumbel copula with Weibull marginals:
  # Compute base terms
  term_H <- (Lambda_H * y1)^(kappa_H * theta_copula)
  term_D <- (Lambda_D * y2)^(kappa_D * theta_copula)
  sum_terms <- term_H + term_D

  # Compute in log-space to avoid overflow
  # copula_factor = sum_terms^(1/theta - 1) * kappa_H * theta * (Lambda_H*y1)^(kappa_H*theta-1) * Lambda_H
  log_power_term <- (1/theta_copula - 1) * log(sum_terms)
  log_weibull_term <- log(kappa_H * theta_copula) +
    (kappa_H * theta_copula - 1) * log(Lambda_H * y1) +
    log(Lambda_H)

  # Combine in log space then exponentiate
  log_result <- log_power_term + log_weibull_term
  result <- exp(log_result)

  # Handle edge cases: replace NaN/Inf with 0 (works with vectors)
  result[!is.finite(result)] <- 0

  return(result)
}

#' Censoring survival function Q_a(t)
#'
#' Computes Q_a(t) = P(C > t | A = a) for exponential censoring.
#'
#' @param t numeric scalar or vector; time
#' @param a numeric scalar (0 or 1); treatment indicator
#' @param lambda_C numeric scalar > 0; baseline censoring rate (can be NULL for no censoring)
#' @param beta_A_C numeric scalar; treatment effect on censoring (can be NULL)
#'
#' @return numeric; censoring survival probability
#'
#' @keywords internal
.Q <- function(t, a, lambda_C, beta_A_C) {
  # Handle no censoring case - return vector of 1s matching length of t
  if (is.null(lambda_C)) {
    return(rep(1, length(t)))
  }
  if (isTRUE(lambda_C <= 0)) {  # Use isTRUE() for safe scalar comparison
    return(rep(1, length(t)))
  }

  Lambda_C <- lambda_C * exp(-beta_A_C * a)
  return(exp(-Lambda_C * t))
}

#' Censoring hazard (constant for exponential)
#'
#' @inheritParams .Q
#'
#' @return numeric; censoring hazard rate
#'
#' @keywords internal
.lambda_C <- function(a, lambda_C, beta_A_C) {
  # Handle no censoring case
  if (is.null(lambda_C)) {
    return(0)
  }
  if (isTRUE(lambda_C <= 0)) {  # Use isTRUE() for safe scalar comparison
    return(0)
  }

  return(lambda_C * exp(-beta_A_C * a))
}

##############################################
# -------- INTEGRATION FUNCTIONS ---------- #
##############################################
# NOTE: For exponential marginals (kappa = 1), we use Zhang's exact formulas
# which are derived from the closed-form expressions for exponential + Gumbel copula

#' Integrand for term A - single integral (Zhang's formula for exponential)
#'
#' Based on Zhang's Ind_singleNu function.
#' For exponential marginals with Gumbel copula.
#'
#' @keywords internal
.integrand_A_exponential <- function(y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                                     theta_copula, lambda_C, beta_A_C) {
  # Zhang's Ind_singleNu:
  # exp(-numda.D*exp(-beta.D)*y2 - numda.C*exp(-beta.C)*y2) *
  # exp(-numda.D*y2 - numda.C*y2) * numda.D

  # Treatment group (a=1) terms
  Lambda_D_1 <- lambda_D * exp(-beta_A_D)
  Lambda_C_1 <- lambda_C * exp(-beta_A_C)

  # Control group (a=0) terms
  Lambda_D_0 <- lambda_D
  Lambda_C_0 <- lambda_C

  return(exp(-Lambda_D_1 * y2 - Lambda_C_1 * y2) *
           exp(-Lambda_D_0 * y2 - Lambda_C_0 * y2) *
           Lambda_D_0)
}

#' Integrand for term B - double integral (Zhang's formula for exponential)
#'
#' Based on Zhang's Ind_doubleNu function.
#' For exponential marginals with Gumbel copula.
#'
#' @keywords internal
.integrand_B_exponential <- function(y1, y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                                     theta_copula, lambda_C, beta_A_C) {
  # Zhang's Ind_doubleNu uses log-space computation:
  # loga: joint survival for treatment group
  # logb: joint survival for control group
  # logc: conditional hazard term
  # logd: sum of censoring hazards

  # Treatment group (a=1) terms
  Lambda_H_1 <- lambda_H * exp(-beta_A_H)
  Lambda_D_1 <- lambda_D * exp(-beta_A_D)
  Lambda_C_1 <- lambda_C * exp(-beta_A_C)

  # Control group (a=0) terms
  Lambda_H_0 <- lambda_H
  Lambda_D_0 <- lambda_D
  Lambda_C_0 <- lambda_C

  # Log of joint survival for treatment (Gumbel copula)
  # S_{D,H,1}(y2, y1) = exp(-((Lambda_H_1*y1)^theta + (Lambda_D_1*y2)^theta)^(1/theta))
  loga <- -(((Lambda_H_1 * y1)^theta_copula + (Lambda_D_1 * y2)^theta_copula)^(1/theta_copula)) -
    Lambda_C_1 * y2

  # Log of joint survival for control (Gumbel copula)
  logb <- -(((Lambda_H_0 * y1)^theta_copula + (Lambda_D_0 * y2)^theta_copula)^(1/theta_copula)) -
    Lambda_C_0 * y2

  # Log of conditional hazard λ_{H|D,0}(y1|y2)
  # This comes from the derivative of the Gumbel copula
  # Pattern from Zhang: log(1/alpha.cor) + log((...)^(1/alpha.cor-1)) + log(alpha.cor*(...))
  sum_terms_0 <- (Lambda_H_0 * y1)^theta_copula + (Lambda_D_0 * y2)^theta_copula

  logc <- log(1/theta_copula) +
    log(sum_terms_0^(1/theta_copula - 1)) +
    log(theta_copula * (Lambda_H_0 * y1)^(theta_copula - 1) * Lambda_H_0)

  # Log of sum of censoring hazards
  logd <- log(Lambda_C_0 * (1 + exp(-beta_A_C)))

  return(exp(loga + logb + logc + logd))
}

#' Integrand for term C - single integral (Zhang's formula for exponential)
#'
#' Based on Zhang's Ind_singleDe function.
#'
#' @keywords internal
.integrand_C_exponential <- function(y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                                     theta_copula, lambda_C, beta_A_C) {
  # Zhang's Ind_singleDe:
  # exp(-numda.D*exp(-beta.D)*y2 - numda.C*exp(-beta.C)*y2) *
  # exp(-numda.D*y2 - numda.C*y2) * numda.D*exp(-beta.D)

  # Treatment group (a=1) terms
  Lambda_D_1 <- lambda_D * exp(-beta_A_D)
  Lambda_C_1 <- lambda_C * exp(-beta_A_C)

  # Control group (a=0) terms
  Lambda_D_0 <- lambda_D
  Lambda_C_0 <- lambda_C

  return(exp(-Lambda_D_1 * y2 - Lambda_C_1 * y2) *
           exp(-Lambda_D_0 * y2 - Lambda_C_0 * y2) *
           Lambda_D_1)
}

#' Integrand for term D - double integral (Zhang's formula for exponential)
#'
#' Based on Zhang's Ind_doubleDe function.
#'
#' @keywords internal
.integrand_D_exponential <- function(y1, y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                                     theta_copula, lambda_C, beta_A_C) {
  # Treatment group (a=1) terms
  Lambda_H_1 <- lambda_H * exp(-beta_A_H)
  Lambda_D_1 <- lambda_D * exp(-beta_A_D)
  Lambda_C_1 <- lambda_C * exp(-beta_A_C)

  # Control group (a=0) terms
  Lambda_H_0 <- lambda_H
  Lambda_D_0 <- lambda_D
  Lambda_C_0 <- lambda_C

  # Log of joint survival for treatment
  loga <- -(((Lambda_H_1 * y1)^theta_copula + (Lambda_D_1 * y2)^theta_copula)^(1/theta_copula)) -
    Lambda_C_1 * y2

  # Log of joint survival for control
  logb <- -(((Lambda_H_0 * y1)^theta_copula + (Lambda_D_0 * y2)^theta_copula)^(1/theta_copula)) -
    Lambda_C_0 * y2

  # Log of conditional hazard λ_{H|D,1}(y1|y2) - note treatment group!
  sum_terms_1 <- (Lambda_H_1 * y1)^theta_copula + (Lambda_D_1 * y2)^theta_copula

  logc <- log(1/theta_copula) +
    log(sum_terms_1^(1/theta_copula - 1)) +
    log(theta_copula * (Lambda_H_1 * y1)^(theta_copula - 1) * Lambda_H_1)

  # Log of sum of censoring hazards (same as B)
  logd <- log(Lambda_C_0 * (1 + exp(-beta_A_C)))

  return(exp(loga + logb + logc + logd))
}

#' Compute term A(x, u) via single integration
#'
#' A(x,u) = ∫ S_{D,1}(y2|x,u) Q_1(y2) S_{D,0}(y2|x,u) Q_0(y2) λ_{D,0}(y2|x,u) dy2
#'
#' @keywords internal
.compute_A <- function(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                       beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                       theta_copula, lambda_C, beta_A_C) {

  # Special case: Exponential with no covariates - use Zhang's optimized formula
  if (kappa_H == 1 && kappa_D == 1 && is.null(x) && is.null(u)) {
    integrand <- function(y2) {
      .integrand_A_exponential(y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                               theta_copula, lambda_C, beta_A_C)
    }
  } else {
    # General case: Weibull marginals and/or covariates
    integrand <- function(y2) {
      S_D_1 <- .S_D(y2, a = 1, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
      S_D_0 <- .S_D(y2, a = 0, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
      Q_1 <- .Q(y2, a = 1, lambda_C, beta_A_C)
      Q_0 <- .Q(y2, a = 0, lambda_C, beta_A_C)
      lambda_D_0 <- .lambda_D(y2, a = 0, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)

      return(S_D_1 * Q_1 * S_D_0 * Q_0 * lambda_D_0)
    }
  }

  result <- integrate(integrand, lower = 0, upper = Inf,
                      rel.tol = .Machine$double.eps^0.25,
                      subdivisions = 1000L)

  return(result$value)
}

#' Compute term B(x, u) via double integration
#'
#' B(x,u) = ∫∫ G_1(y1,y2|x,u) Q_1(y2) G_0(y1,y2|x,u) Q_0(y2)
#'           λ_{H|D,0}(y1|y2,x,u) {λ_{C,0}(y2)+λ_{C,1}(y2)} dy1 dy2
#'
#' @keywords internal
.compute_B <- function(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                       beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                       theta_copula, lambda_C, beta_A_C) {

  # Special case: Exponential with no covariates - use Zhang's optimized formula
  if (kappa_H == 1 && kappa_D == 1 && is.null(x) && is.null(u)) {
    integrand_outer <- function(y2) {
      integrand_inner <- function(y1) {
        .integrand_B_exponential(y1, y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                                 theta_copula, lambda_C, beta_A_C)
      }

      inner_result <- integrate(integrand_inner, lower = 0, upper = y2,
                                rel.tol = .Machine$double.eps^0.25,
                                subdivisions = 1000L)
      return(inner_result$value)
    }
  } else {
    # General case: Weibull marginals and/or covariates
    integrand_outer <- function(y2) {
      integrand_inner <- function(y1) {
        G_1 <- .G(y1, y2, a = 1, x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                  beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                  theta_copula)
        G_0 <- .G(y1, y2, a = 0, x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                  beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                  theta_copula)
        Q_1 <- .Q(y2, a = 1, lambda_C, beta_A_C)
        Q_0 <- .Q(y2, a = 0, lambda_C, beta_A_C)
        lambda_H_D_0 <- .lambda_H_given_D(y1, y2, a = 0, x, u,
                                          lambda_H, lambda_D, kappa_H, kappa_D,
                                          beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                                          theta_copula)
        lambda_C_0 <- .lambda_C(a = 0, lambda_C, beta_A_C)
        lambda_C_1 <- .lambda_C(a = 1, lambda_C, beta_A_C)

        return(G_1 * Q_1 * G_0 * Q_0 * lambda_H_D_0 * (lambda_C_0 + lambda_C_1))
      }

      inner_result <- integrate(integrand_inner, lower = 0, upper = y2,
                                rel.tol = .Machine$double.eps^0.25,
                                subdivisions = 1000L)
      return(inner_result$value)
    }
  }

  result <- integrate(Vectorize(integrand_outer), lower = 0, upper = Inf,
                      rel.tol = .Machine$double.eps^0.25,
                      subdivisions = 1000L)

  return(result$value)
}

#' Compute term C(x, u) via single integration
#'
#' C(x,u) = ∫ S_{D,1}(y2|x,u) Q_1(y2) S_{D,0}(y2|x,u) Q_0(y2) λ_{D,1}(y2|x,u) dy2
#'
#' @keywords internal
.compute_C <- function(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                       beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                       theta_copula, lambda_C, beta_A_C) {

  # Special case: Exponential with no covariates - use Zhang's optimized formula
  if (kappa_H == 1 && kappa_D == 1 && is.null(x) && is.null(u)) {
    integrand <- function(y2) {
      .integrand_C_exponential(y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                               theta_copula, lambda_C, beta_A_C)
    }
  } else {
    # General case: Weibull marginals and/or covariates
    integrand <- function(y2) {
      S_D_1 <- .S_D(y2, a = 1, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
      S_D_0 <- .S_D(y2, a = 0, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
      Q_1 <- .Q(y2, a = 1, lambda_C, beta_A_C)
      Q_0 <- .Q(y2, a = 0, lambda_C, beta_A_C)
      lambda_D_1 <- .lambda_D(y2, a = 1, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)

      return(S_D_1 * Q_1 * S_D_0 * Q_0 * lambda_D_1)
    }
  }

  result <- integrate(integrand, lower = 0, upper = Inf,
                      rel.tol = .Machine$double.eps^0.25,
                      subdivisions = 1000L)

  return(result$value)
}

#' Compute term D(x, u) via double integration
#'
#' D(x,u) = ∫∫ G_1(y1,y2|x,u) Q_1(y2) G_0(y1,y2|x,u) Q_0(y2)
#'           λ_{H|D,1}(y1|y2,x,u) {λ_{C,0}(y2)+λ_{C,1}(y2)} dy1 dy2
#'
#' @keywords internal
.compute_D <- function(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                       beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                       theta_copula, lambda_C, beta_A_C) {

  # Special case: Exponential with no covariates - use Zhang's optimized formula
  if (kappa_H == 1 && kappa_D == 1 && is.null(x) && is.null(u)) {
    integrand_outer <- function(y2) {
      integrand_inner <- function(y1) {
        .integrand_D_exponential(y1, y2, lambda_H, lambda_D, beta_A_H, beta_A_D,
                                 theta_copula, lambda_C, beta_A_C)
      }

      inner_result <- integrate(integrand_inner, lower = 0, upper = y2,
                                rel.tol = .Machine$double.eps^0.25,
                                subdivisions = 1000L)
      return(inner_result$value)
    }
  } else {
    # General case: Weibull marginals and/or covariates
    integrand_outer <- function(y2) {
      integrand_inner <- function(y1) {
        G_1 <- .G(y1, y2, a = 1, x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                  beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                  theta_copula)
        G_0 <- .G(y1, y2, a = 0, x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                  beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                  theta_copula)
        Q_1 <- .Q(y2, a = 1, lambda_C, beta_A_C)
        Q_0 <- .Q(y2, a = 0, lambda_C, beta_A_C)
        lambda_H_D_1 <- .lambda_H_given_D(y1, y2, a = 1, x, u,
                                          lambda_H, lambda_D, kappa_H, kappa_D,
                                          beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                                          theta_copula)
        lambda_C_0 <- .lambda_C(a = 0, lambda_C, beta_A_C)
        lambda_C_1 <- .lambda_C(a = 1, lambda_C, beta_A_C)

        return(G_1 * Q_1 * G_0 * Q_0 * lambda_H_D_1 * (lambda_C_0 + lambda_C_1))
      }

      inner_result <- integrate(integrand_inner, lower = 0, upper = y2,
                                rel.tol = .Machine$double.eps^0.25,
                                subdivisions = 1000L)
      return(inner_result$value)
    }
  }

  result <- integrate(Vectorize(integrand_outer), lower = 0, upper = Inf,
                      rel.tol = .Machine$double.eps^0.25,
                      subdivisions = 1000L)

  return(result$value)
}

#' Helper function to compute E_{X,U}[f(X,U)] for small dimensions
#'
#' Uses quadrature integration over Normal covariates and enumeration over Bernoulli.
#'
#' @keywords internal
.integrate_over_covariates <- function(func, p, q, X_dist, U_dist,
                                       mu_X, sd_X, prob_X,
                                       mu_U, sd_U, prob_U) {

  # For now, use simple quadrature approach (Gauss-Hermite for Normal)
  # This works well for p + q <= 3 or 4

  if ((p + q) > 4) {
    warning("Covariate integration with dimension > 4 may be slow or inaccurate. ",
            "Consider using Monte Carlo approximation instead.")
  }

  # Separate Normal and Bernoulli covariates
  X_normal_idx <- which(X_dist == "Normal")
  X_bern_idx <- which(X_dist == "Bernoulli")
  U_normal_idx <- which(U_dist == "Normal")
  U_bern_idx <- which(U_dist == "Bernoulli")

  n_normal <- length(X_normal_idx) + length(U_normal_idx)
  n_bern <- length(X_bern_idx) + length(U_bern_idx)

  # For simplicity: if we have any Bernoulli, enumerate all combinations
  # For Normal, use Gauss-Hermite quadrature

  if (n_bern > 0 && n_bern <= 10) {
    # Enumerate Bernoulli combinations
    bern_grid <- expand.grid(rep(list(c(0, 1)), n_bern))

    # Compute probability of each Bernoulli combination
    bern_probs <- numeric(nrow(bern_grid))
    for (i in 1:nrow(bern_grid)) {
      prob <- 1
      col_idx <- 1
      for (j in X_bern_idx) {
        prob <- prob * ifelse(bern_grid[i, col_idx] == 1, prob_X[j], 1 - prob_X[j])
        col_idx <- col_idx + 1
      }
      for (j in U_bern_idx) {
        prob <- prob * ifelse(bern_grid[i, col_idx] == 1, prob_U[j], 1 - prob_U[j])
        col_idx <- col_idx + 1
      }
      bern_probs[i] <- prob
    }
  } else if (n_bern > 10) {
    stop("Too many Bernoulli covariates (>10) for exact integration. ",
         "Consider reducing dimension or using Monte Carlo.")
  }

  # For Normal covariates, use Gauss-Hermite quadrature
  # We'll use statmod package's gauss.quad.prob if available, otherwise simple approach
  if (n_normal > 0) {
    # Simple approach: use equally-spaced grid over [-3sd, +3sd] for each Normal
    # For better accuracy, could use gauss.quad.prob from statmod package
    n_quad_points <- 10  # per dimension

    if (n_normal == 0) {
      # No normal covariates, just enumerate Bernoulli
      expectation <- 0
      for (i in 1:nrow(bern_grid)) {
        x_val <- if (p > 0) as.numeric(bern_grid[i, 1:length(X_bern_idx)]) else NULL
        u_val <- if (q > 0) as.numeric(bern_grid[i, (length(X_bern_idx)+1):n_bern]) else NULL
        expectation <- expectation + bern_probs[i] * func(x_val, u_val)
      }
      return(expectation)
    } else {
      stop("General covariate integration with Normal covariates not yet fully implemented. ",
           "For now, please use p = 0, q = 0 (no covariates) or all Bernoulli covariates.")
    }
  } else {
    # Only Bernoulli covariates
    expectation <- 0
    for (i in 1:nrow(bern_grid)) {
      # Split bern_grid row into X and U parts
      x_val <- if (p > 0) as.numeric(bern_grid[i, 1:length(X_bern_idx)]) else NULL
      u_val <- if (q > 0) {
        start_idx <- length(X_bern_idx) + 1
        as.numeric(bern_grid[i, start_idx:n_bern])
      } else {
        NULL
      }
      expectation <- expectation + bern_probs[i] * func(x_val, u_val)
    }
    return(expectation)
  }
}

############################################
# ------------ MAIN FUNCTION ------------ #
############################################

#' Calculate the true causal win ratio analytically
#'
#' Computes the true causal win ratio using numerical integration over the
#' covariate distribution, following Zhang et al. (2019, 2022) and Luo et al. (2015).
#'
#' @param p integer; number of measured covariates (can be 0 or NULL)
#' @param q integer; number of unmeasured covariates (can be 0 or NULL)
#' @param X_dist character vector length p; distribution type for each measured covariate,
#'   either "Normal" or "Bernoulli" (defaults to "Normal" for all if NULL)
#' @param U_dist character vector length q; distribution type for each unmeasured covariate,
#'   either "Normal" or "Bernoulli" (defaults to "Normal" for all if NULL)
#' @param mu_X numeric vector length p; means for Normal measured covariates
#' @param sd_X numeric vector length p; standard deviations for Normal measured covariates
#' @param prob_X numeric vector length p; probabilities for Bernoulli measured covariates
#' @param mu_U numeric vector length q; means for Normal unmeasured covariates
#' @param sd_U numeric vector length q; standard deviations for Normal unmeasured covariates
#' @param prob_U numeric vector length q; probabilities for Bernoulli unmeasured covariates
#' @param lambda_H numeric scalar > 0; baseline hazard rate for non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape parameter for non-fatal event (use 1 for exponential)
#' @param beta_A_H numeric scalar; treatment effect on non-fatal event
#' @param beta_X_H numeric vector length p; measured covariate effects on non-fatal event (can be NULL)
#' @param beta_U_H numeric vector length q; unmeasured covariate effects on non-fatal event (can be NULL)
#' @param lambda_D numeric scalar > 0; baseline hazard rate for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape parameter for fatal event (use 1 for exponential)
#' @param beta_A_D numeric scalar; treatment effect on fatal event
#' @param beta_X_D numeric vector length p; measured covariate effects on fatal event (can be NULL)
#' @param beta_U_D numeric vector length q; unmeasured covariate effects on fatal event (can be NULL)
#' @param theta_copula numeric >= 1; Gumbel-Hougaard copula parameter (use 1 for independence)
#' @param lambda_C numeric scalar > 0; baseline censoring rate (NULL for no censoring)
#' @param beta_A_C numeric scalar; treatment effect on censoring (NULL for no censoring)
#'
#' @return numeric; the true causal win ratio WR_causal
#'
#' @details
#' This function computes the true causal win ratio as defined in Zhang et al. (2019):
#'
#' WR_causal = E_{X,U}[A(X,U) + B(X,U)] / E_{X,U}[C(X,U) + D(X,U)]
#'
#' where A, B, C, D are integrals over the joint distribution of (T_H, T_D, C) under
#' treatment and control. See Zhang (2019) Section 3.2.1 and your README for details.
#'
#' When p = q = 0 (no covariates), this reduces to the simpler case where we just
#' compute (A + B) / (C + D) directly.
#'
#' The function handles both exponential (kappa = 1) and Weibull (kappa ≠ 1) marginal
#' distributions.
#'
#' @examples
#' # Example 1: Replicate Zhang's case with no treatment effect (should give WR = 1)
#' true_wr_analytic(
#'   p = 0, q = 0,
#'   lambda_H = 0.1, kappa_H = 1, beta_A_H = 0,
#'   lambda_D = 0.08, kappa_D = 1, beta_A_D = 0,
#'   theta_copula = 2,
#'   lambda_C = 0.09, beta_A_C = 0.1
#' )
#'
#' # Example 2: With treatment effect (eta_H = 0.1, eta_D = 0.1)
#' true_wr_analytic(
#'   p = 0, q = 0,
#'   lambda_H = 0.1, kappa_H = 1, beta_A_H = 0.1,
#'   lambda_D = 0.08, kappa_D = 1, beta_A_D = 0.1,
#'   theta_copula = 2,
#'   lambda_C = 0.09, beta_A_C = 0.1
#' )
#'
#' @export
true_wr_analytic <- function(
    p = NULL, q = NULL,
    X_dist = NULL, U_dist = NULL,
    mu_X = NULL, sd_X = NULL, prob_X = NULL,
    mu_U = NULL, sd_U = NULL, prob_U = NULL,
    lambda_H, kappa_H, beta_A_H,
    beta_X_H = NULL, beta_U_H = NULL,
    lambda_D, kappa_D, beta_A_D,
    beta_X_D = NULL, beta_U_D = NULL,
    theta_copula = 1,
    lambda_C = NULL, beta_A_C = NULL
) {

  # Handle NULL p and q
  if (is.null(p)) p <- 0L
  if (is.null(q)) q <- 0L

  # Set default distributions if not provided
  if (p > 0 && is.null(X_dist)) X_dist <- rep("Normal", p)
  if (q > 0 && is.null(U_dist)) U_dist <- rep("Normal", q)

  # Case 1: No covariates (simplest and fastest)
  if (p == 0 && q == 0) {
    # Compute A, B, C, D directly without integration over covariates
    A_val <- .compute_A(x = NULL, u = NULL,
                        lambda_H, lambda_D, kappa_H, kappa_D,
                        beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                        theta_copula, lambda_C, beta_A_C)

    B_val <- .compute_B(x = NULL, u = NULL,
                        lambda_H, lambda_D, kappa_H, kappa_D,
                        beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                        theta_copula, lambda_C, beta_A_C)

    C_val <- .compute_C(x = NULL, u = NULL,
                        lambda_H, lambda_D, kappa_H, kappa_D,
                        beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                        theta_copula, lambda_C, beta_A_C)

    D_val <- .compute_D(x = NULL, u = NULL,
                        lambda_H, lambda_D, kappa_H, kappa_D,
                        beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                        theta_copula, lambda_C, beta_A_C)

    numerator <- A_val + B_val
    denominator <- C_val + D_val

    WR_causal <- numerator / denominator

    return(WR_causal)
  }

  # Case 2: With covariates - integrate E_{X,U}[(A+B)/(C+D)]
  # For small covariate dimension, use quadrature/enumeration

  # Define function to integrate: (A+B) for given (x,u)
  numerator_func <- function(x, u) {
    A <- .compute_A(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                    beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                    theta_copula, lambda_C, beta_A_C)
    B <- .compute_B(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                    beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                    theta_copula, lambda_C, beta_A_C)
    return(A + B)
  }

  # Define function to integrate: (C+D) for given (x,u)
  denominator_func <- function(x, u) {
    C <- .compute_C(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                    beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                    theta_copula, lambda_C, beta_A_C)
    D <- .compute_D(x, u, lambda_H, lambda_D, kappa_H, kappa_D,
                    beta_A_H, beta_A_D, beta_X_H, beta_X_D, beta_U_H, beta_U_D,
                    theta_copula, lambda_C, beta_A_C)
    return(C + D)
  }

  # Integrate over (X, U) distribution
  E_numerator <- .integrate_over_covariates(
    numerator_func, p, q, X_dist, U_dist,
    mu_X, sd_X, prob_X, mu_U, sd_U, prob_U
  )

  E_denominator <- .integrate_over_covariates(
    denominator_func, p, q, X_dist, U_dist,
    mu_X, sd_X, prob_X, mu_U, sd_U, prob_U
  )

  WR_causal <- E_numerator / E_denominator

  return(WR_causal)
}
