# The functions in this package are modified versions of code used in Di Zhang's
# `cWR` package, found here: https://github.com/dee1008/cWR


##############################################
# ------------ HELPER FUNCTIONS ------------ #
##############################################
# TODO: Move helper functions to their own file

#' Linear predictor for fatal event
#'
#' Compute the linear predictor for the fatal event hazard given treatment, measured covariates x, and unmeasured covariates u.
#'
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param beta_A_D numeric scalar; treatment coefficient; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u; can also be \code{NULL}
#'
#' @return numeric scalar; linear predictor for the fatal event
#'
#' @details One of x or u (and the corresponding beta parameter) can be \code{NULL}, but NOT both.
#'
#' @keywords internal
.linpred_D <- function(a, x, u, beta_A_D, beta_X_D, beta_U_D) {
  out <- 0
  if (!is.null(beta_A_D)){
    out <- out + beta_A_D * a
  }
  if (!is.null(beta_X_D) && !is.null(x) && length(beta_X_D) == length(x)) {
    out <- out + sum(beta_X_D * x)
  }
  if (!is.null(beta_U_D) && !is.null(u) && length(beta_U_D) == length(u)) {
    out <- out + sum(beta_U_D * u)
  }
  return(out)
}

#' Linear predictor for non-fatal event
#'
#' Compute the linear predictor for the non-fatal event hazard given treatment, measured covariates x, and unmeasured covariates u.
#'
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param beta_A_H numeric scalar; treatment coefficient; can also be \code{NULL}
#' @param beta_X_H numeric vector; coefficients for x; can also be \code{NULL}
#' @param beta_U_H numeric vector; coefficients for u; can also be \code{NULL}
#'
#' @return numeric scalar; linear predictor for the non-fatal event
#'
#' @details One of x or u (and the corresponding beta parameter) can be \code{NULL}, but NOT both.
#'
#' @keywords internal
.linpred_H <- function(a, x, u, beta_A_H, beta_X_H, beta_U_H) {
  out <- 0
  if (!is.null(beta_A_H)){
    out <- out + beta_A_H * a
  }
  if (!is.null(beta_X_H) && !is.null(x) && length(beta_X_H) == length(x)) {
    out <- out + sum(beta_X_H * x)
  }
  if (!is.null(beta_U_H) && !is.null(u) && length(beta_U_H) == length(u)) {
    out <- out + sum(beta_U_H * u)
  }
  return(out)
}

#' Weibull scale parameter for fatal event
#'
#' Compute the scale parameter for the Weibull fatal event distribution.
#'
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard rate
#' @param beta_A_D numeric scalar; treatment coefficient; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u; can also be \code{NULL}
#'
#' @return numeric scalar; Weibull scale parameter for fatal event
#'
#' @keywords internal
.Lambda_D <- function(a, x, u, lambda_D, beta_A_D, beta_X_D, beta_U_D) {
  lp <- .linpred_D(a, x, u, beta_A_D, beta_X_D, beta_U_D)
  return(lambda_D * exp(-lp))
}

#' Weibull scale parameter for non-fatal event
#'
#' Compute the scale parameter for the Weibull non-fatal event distribution.
#'
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard rate
#' @param beta_A_H numeric scalar; treatment coefficient; can also be \code{NULL}
#' @param beta_X_H numeric vector; coefficients for x; can also be \code{NULL}
#' @param beta_U_H numeric vector; coefficients for u; can also be \code{NULL}
#'
#' @return numeric scalar; Weibull scale parameter for non-fatal event
#'
#' @keywords internal
.Lambda_H <- function(a, x, u, lambda_H, beta_A_H, beta_X_H, beta_U_H) {
  lp <- .linpred_H(a, x, u, beta_A_H, beta_X_H, beta_U_H)
  return(lambda_H * exp(-lp))
}

#' Censoring hazard for random censoring
#'
#' Compute the hazard rate for random censoring. If \code{lambda_C} is NULL,
#' there is no random censoring and the hazard is zero.
#'
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param lambda_C numeric scalar > 0; baseline censoring rate; can also be \code{NULL}
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#'
#' @return numeric scalar; censoring hazard
#'
#' @keywords internal
.Lambda_C <- function(a, lambda_C, beta_A_C) {
  if (is.null(lambda_C)){
    return(0)
  }
  else{
    return(lambda_C * exp(-beta_A_C * a))
  }
}

#' Censoring survival function Q(t | A=a)
#'
#' Compute the censoring survival probability Q(t | A=a). If
#' \code{lambda_C} is NULL, there is no random censoring and Q(t|a) = 1.
#'
#' @param t numeric scalar; time
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param lambda_C numeric scalar > 0; baseline censoring rate; can also be \code{NULL}
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#'
#' @return numeric scalar; censoring survival probability
#'
#' @keywords internal
.Q <- function(t, a, lambda_C, beta_A_C) {
  if (is.null(lambda_C)){
    return(1)
  }
  else{
    return(exp(-.Lambda_C(a, lambda_C, beta_A_C) * t))
  }

}

#' Marginal Weibull survival for fatal event
#'
#' Compute the marginal survival function S_D(t | A,X,U) for the fatal event.
#'
#' @param t numeric scalar; time
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard rate
#' @param kappa_D numeric scalar > 0; Weibull shape parameter
#' @param beta_A_D numeric scalar; treatment coefficient; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u; can also be \code{NULL}
#'
#' @return numeric scalar; survival probability
#'
#' @keywords internal
.S_D <- function(t, a, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D) {
  lam <- .Lambda_D(a, x, u, lambda_D, beta_A_D, beta_X_D, beta_U_D)
  return(exp(-(lam * t)^kappa_D))
}

#' Marginal Weibull survival for non-fatal event
#'
#' Compute the marginal survival function S_H(t | A,X,U) for the non-fatal event.
#'
#' @param t numeric scalar; time
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard rate
#' @param kappa_H numeric scalar > 0; Weibull shape parameter
#' @param beta_A_H numeric scalar; treatment coefficient; can also be \code{NULL}
#' @param beta_X_H numeric vector; coefficients for x; can also be \code{NULL}
#' @param beta_U_H numeric vector; coefficients for u; can also be \code{NULL}
#'
#' @return numeric scalar; survival probability
#'
#' @keywords internal
.S_H <- function(t, a, x, u, lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H) {
  lam <- .Lambda_H(a, x, u, lambda_H, beta_A_H, beta_X_H, beta_U_H)
  return(exp(-(lam * t)^kappa_H))
}

#' Joint survival S_{D,H}(t_D, t_H | A,X,U) via Gumbel–Hougaard copula
#'
#' Compute the joint survival for (T_D, T_H) under a Gumbel–Hougaard copula.
#' When \code{theta_copula = 1}, this reduces to independence
#' S_D(t_D) * S_H(t_H).
#'
#' @param y2 numeric scalar; time for fatal event T_D
#' @param y1 numeric scalar; time for non-fatal event T_H
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard for non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape for non-fatal event
#' @param beta_A_H numeric scalar; treatment coefficient for non-fatal event; can also be \code{NULL}
#' @param beta_X_H numeric vector; coefficients for x on non-fatal event; can also be \code{NULL}
#' @param beta_U_H numeric vector; coefficients for u on non-fatal event; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard parameter
#'
#' @return numeric scalar; joint survival probability
#'
#' @keywords internal
.S_DH <- function(y2, y1, a, x, u,
                  lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                  lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                  theta_copula) {
  lamD <- .Lambda_D(a, x, u, lambda_D, beta_A_D, beta_X_D, beta_U_D)
  lamH <- .Lambda_H(a, x, u, lambda_H, beta_A_H, beta_X_H, beta_U_H)

  # If we have independent events
  if (theta_copula == 1) {
    return( exp(-(lamD * y2)^kappa_D) * exp(-(lamH * y1)^kappa_H) )
  }
  else{
    term <- (lamD * y2)^(kappa_D * theta_copula) + (lamH * y1)^(kappa_H * theta_copula)
    return(exp(-term^(1 / theta_copula)))
  }
}

#' Conditional survival of T_H > y1 given T_D > y2
#'
#' Compute G(y1, y2 | A,X,U) = P(T_H > y1 | T_D > y2, A,X,U).
#'
#' @param y1 numeric scalar; time for non-fatal event
#' @param y2 numeric scalar; time for fatal event
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard for non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape for non-fatal event
#' @param beta_A_H numeric scalar; treatment coefficient for non-fatal event; can also be \code{NULL}
#' @param beta_X_H numeric vector; coefficients for x on non-fatal event; can also be \code{NULL}
#' @param beta_U_H numeric vector; coefficients for u on non-fatal event; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard parameter
#'
#' @return numeric scalar; conditional survival probability
#'
#' @keywords internal
.G_cond <- function(y1, y2, a, x, u,
                    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                    theta_copula) {
  num <- .S_DH(y2, y1, a, x, u,
               lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
               lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
               theta_copula)
  den <- .S_D(y2, a, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
  return(num/den)
}

#' Density of T_D
#'
#' Compute the density f_D(y2 | A,X,U) for the fatal event time T_D at y2.
#'
#' @param y2 numeric scalar; time for fatal event
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#'
#' @return numeric scalar; density value
#'
#' @keywords internal
.f_D <- function(y2, a, x, u,
                 lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D) {
  lam <- .Lambda_D(a, x, u, lambda_D, beta_A_D, beta_X_D, beta_U_D)
  t_exp <- (lam * y2)^kappa_D
  return(kappa_D * lam^kappa_D * y2^(kappa_D - 1) * exp(-t_exp))
}

#' Hazard of T_D
#'
#' Compute the hazard lambda_D(y2 | A,X,U) for the fatal event time T_D.
#'
#' @param y2 numeric scalar; time for fatal event
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#'
#' @return numeric scalar; hazard value
#'
#' @keywords internal
.lambda_D_fun <- function(y2, a, x, u,
                          lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D) {
  return(
    .f_D(y2, a, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D) /
      .S_D(y2, a, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
  )
}

#' Conditional hazard lambda_{H|D}
#'
#' Compute the conditional hazard lambda_{H|D}(y1 | y2, A,X,U) via finite differences of the conditional survival G(y1, y2 | A,X,U).
#' This is basically a numerical method of calculating the numerator in the conditional hazard.
#'
#' @param y1 numeric scalar; time for non-fatal event
#' @param y2 numeric scalar; time for fatal event
#' @param a numeric scalar, takes on value 1 or 0; treatment indicator
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard for non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape for non-fatal event
#' @param beta_A_H numeric scalar; treatment coefficient for non-fatal event; can also be \code{NULL}
#' @param beta_X_H numeric vector; coefficients for x on non-fatal event; can also be \code{NULL}
#' @param beta_U_H numeric vector; coefficients for u on non-fatal event; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard parameter
#' @param fd_step numeric scalar > 0; finite difference step size
#'
#' @return numeric scalar; conditional hazard value
#'
#' @keywords internal
.lambda_H_given_D <- function(y1, y2, a, x, u,
                              lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                              lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                              theta_copula,
                              fd_step) {

  if (y1 <= 0){
    return(0)
  }

  h   <- fd_step
  y1p <- y1 + h
  G1  <- .G_cond(y1,  y2, a, x, u,
                 lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                 lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                 theta_copula)
  G1p <- .G_cond(y1p, y2, a, x, u,
                 lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                 lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                 theta_copula)
  logG  <- log(G1)
  logGp <- log(G1p)
  d_logG <- (logGp - logG) / h
  return(-d_logG)
}

#' Nested double integral helper
#'
#' Compute a nested double integral over 0 < y1 < y2 < y2_max using a supplied integrand function of (y1, y2).
#'
#' @param myfun function(y1, y2) returning numeric
#' @param y2_max numeric scalar > 0; upper limit for y2
#' @param rel.tol numeric scalar; relative tolerance for \code{integrate()}
#' @param abs.tol numeric scalar; absolute tolerance for \code{integrate()}
#'
#' @return an object of class "integrate" with components including \code{value}
#'
#' @keywords internal
.double_integral <- function(myfun, y2_max, rel.tol, abs.tol) {
  integrate(
    Vectorize(function(y2) {
      inner <- integrate(
        Vectorize(function(y1) myfun(y1, y2)),
        lower = 0, upper = y2,
        rel.tol = rel.tol, abs.tol = abs.tol
      )
      inner$value
    }),
    lower = 0, upper = y2_max,
    rel.tol = rel.tol, abs.tol = abs.tol
  )
}

#' Individual-setting single integral numerator A(x,u)
#'
#' Compute the integrand for the single integral component A(x,u) in the numerator of the true win ratio.
#'
#' @param y2 numeric scalar; variable of integration, time for fatal event
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard for non-fatal event (unused)
#' @param kappa_H numeric scalar > 0; Weibull shape for non-fatal event (unused)
#' @param beta_A_H numeric scalar; treatment coefficient for non-fatal event; can also be \code{NULL} (unused)
#' @param beta_X_H numeric vector; coefficients for x on non-fatal event; can also be \code{NULL} (unused)
#' @param beta_U_H numeric vector; coefficients for u on non-fatal event; can also be \code{NULL} (unused)
#' @param lambda_C numeric scalar > 0; baseline censoring rate; can also be \code{NULL}
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard parameter (unused)
#' @param fd_step numeric scalar > 0; finite difference step size
#'
#' @return numeric scalar; value of the integrand at y2
#'
#' @keywords internal
.Ind_singleNu <- function(
    y2, x, u,
    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
    lambda_C, beta_A_C,
    theta_copula,
    fd_step
) {
  a1 <- 1
  a0 <- 0

  S_D1 <- .S_D(y2, a1, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
  S_D0 <- .S_D(y2, a0, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
  Q1   <- .Q(y2, a1, lambda_C, beta_A_C)
  Q0   <- .Q(y2, a0, lambda_C, beta_A_C)
  lamD0 <- .lambda_D_fun(y2, a0, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)

  loga <- log(S_D1) + log(Q1)
  logb <- log(S_D0) + log(Q0)
  logc <- log(lamD0)

  return(exp(loga + logb + logc))
}

#' Individual-setting single integral denominator C(x,u)
#'
#' Compute the integrand for the single integral component C(x,u) in the denominator of the true win ratio.
#'
#' @param y2 numeric scalar; variable of integration, time for fatal event
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard for non-fatal event (unused)
#' @param kappa_H numeric scalar > 0; Weibull shape for non-fatal event (unused)
#' @param beta_A_H numeric scalar; treatment coefficient for non-fatal event; can also be \code{NULL} (unused)
#' @param beta_X_H numeric vector; coefficients for x on non-fatal event; can also be \code{NULL} (unused)
#' @param beta_U_H numeric vector; coefficients for u on non-fatal event; can also be \code{NULL} (unused)
#' @param lambda_C numeric scalar > 0; baseline censoring rate; can also be \code{NULL}
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard parameter (unused)
#' @param fd_step numeric scalar > 0; finite difference step size (unused)
#'
#' @return numeric scalar; value of the integrand at y2
#'
#' @keywords internal
.Ind_singleDe <- function(
    y2, x, u,
    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
    lambda_C, beta_A_C,
    theta_copula,
    fd_step
) {
  a1 <- 1
  a0 <- 0

  S_D1 <- .S_D(y2, a1, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
  S_D0 <- .S_D(y2, a0, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)
  Q1   <- .Q(y2, a1, lambda_C, beta_A_C)
  Q0   <- .Q(y2, a0, lambda_C, beta_A_C)

  lamD1 <- .lambda_D_fun(y2, a1, x, u, lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D)

  loga <- log(S_D1) + log(Q1)
  logb <- log(S_D0) + log(Q0)
  logc <- log(lamD1)

  exp(loga + logb + logc)
}

#' Individual-setting double integral numerator B(x,u)
#'
#' Compute the integrand for the double integral component B(x,u) in the numerator of the true win ratio.
#'
#' @param y1 numeric scalar; variable of integration, time for non-fatal event
#' @param y2 numeric scalar; variable of integration, time for fatal event
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard for non-fatal event (unused)
#' @param kappa_H numeric scalar > 0; Weibull shape for non-fatal event (unused)
#' @param beta_A_H numeric scalar; treatment coefficient for non-fatal event; can also be \code{NULL} (unused)
#' @param beta_X_H numeric vector; coefficients for x on non-fatal event; can also be \code{NULL} (unused)
#' @param beta_U_H numeric vector; coefficients for u on non-fatal event; can also be \code{NULL} (unused)
#' @param lambda_C numeric scalar > 0; baseline censoring rate; can also be \code{NULL}
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard parameter (unused)
#' @param fd_step numeric scalar > 0; finite difference step size (unused)
#'
#' @return numeric scalar; value of the integrand at (y1, y2)
#'
#' @keywords internal
.Ind_doubleNu <- function(
    y1, y2, x, u,
    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
    lambda_C, beta_A_C,
    theta_copula,
    fd_step
) {
  a1 <- 1
  a0 <- 0

  G1 <- .G_cond(y1, y2, a1, x, u,
                lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                theta_copula)
  G0 <- .G_cond(y1, y2, a0, x, u,
                lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                theta_copula)

  Q1 <- .Q(y2, a1, lambda_C, beta_A_C)
  Q0 <- .Q(y2, a0, lambda_C, beta_A_C)

  lamH0 <- .lambda_H_given_D(
    y1, y2, a0, x, u,
    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
    theta_copula,
    fd_step
  )

  lamC0 <- .Lambda_C(a0, lambda_C, beta_A_C)
  lamC1 <- .Lambda_C(a1, lambda_C, beta_A_C)
  sumC  <- lamC0 + lamC1

  # If we have no random censoring, then this component vanishes
  if (sumC <= 0) {
    return(0)
  }
  else{
    loga <- log(G1) + log(Q1)
    logb <- log(G0) + log(Q0)
    logc <- log(lamH0)
    logd <- log(sumC)

    return(exp(loga + logb + logc + logd))
  }
}

#' Individual-setting double integral denominator D(x,u)
#'
#' Compute the integrand for the double integral component D(x,u) in the denominator of the true win ratio.
#'
#' @param y1 numeric scalar; variable of integration, time for non-fatal event
#' @param y2 numeric scalar; variable of integration, time for fatal event
#' @param x numeric vector; measured covariates; can also be \code{NULL}
#' @param u numeric vector; unmeasured covariates; can also be \code{NULL}
#' @param lambda_D numeric scalar > 0; baseline hazard for fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape for fatal event
#' @param beta_A_D numeric scalar; treatment coefficient for fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector; coefficients for x on fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector; coefficients for u on fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard for non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape for non-fatal event
#' @param beta_A_H numeric scalar; treatment coefficient for non-fatal event; can also be \code{NULL}
#' @param beta_X_H numeric vector; coefficients for x on non-fatal event; can also be \code{NULL}
#' @param beta_U_H numeric vector; coefficients for u on non-fatal event; can also be \code{NULL}
#' @param lambda_C numeric scalar > 0; baseline censoring rate; can also be \code{NULL}
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard parameter
#' @param fd_step numeric scalar > 0; finite difference step size
#'
#' @return numeric scalar; value of the integrand at (y1, y2)
#'
#' @keywords internal
.Ind_doubleDe <- function(
    y1, y2, x, u,
    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
    lambda_C, beta_A_C,
    theta_copula,
    fd_step
) {
  a1 <- 1
  a0 <- 0

  G1 <- .G_cond(y1, y2, a1, x, u,
                lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                theta_copula)
  G0 <- .G_cond(y1, y2, a0, x, u,
                lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
                lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
                theta_copula)

  Q1 <- .Q(y2, a1, lambda_C, beta_A_C)
  Q0 <- .Q(y2, a0, lambda_C, beta_A_C)

  lamH1 <- .lambda_H_given_D(
    y1, y2, a1, x, u,
    lambda_D, kappa_D, beta_A_D, beta_X_D, beta_U_D,
    lambda_H, kappa_H, beta_A_H, beta_X_H, beta_U_H,
    theta_copula,
    fd_step
  )

  lamC0 <- .Lambda_C(a0, lambda_C, beta_A_C)
  lamC1 <- .Lambda_C(a1, lambda_C, beta_A_C)
  sumC  <- lamC0 + lamC1

  # If we have no random censoring, then this component vanishes
  if (sumC <= 0) {
    return(0)
  }
  else{
    loga <- log(G1) + log(Q1)
    logb <- log(G0) + log(Q0)
    logc <- log(lamH1)
    logd <- log(sumC)

    return(exp(loga + logb + logc + logd))
  }
}

#' True causal win ratio for a single covariate pattern
#'
#' Compute the analytic "true causal win ratio" for fixed measured and unmeasured covariates x, u.
#' This corresponds to
#' WR(x, u) = {A(x,u) + B(x,u)} / {C(x,u) + D(x,u)}
#' where A, B, C, D are the components defined via integrals over the joint distribution of (T_H, T_D, C) under treatment and control.
#'
#' @param x numeric vector length p; measured covariates (or NULL, as long as \code{u} is not null)
#' @param u numeric vector length q; unmeasured covariates (or NULL, as long as \code{x} is not null)
#' @param lambda_D numeric scalar > 0; baseline hazard rate for the fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape parameter for the fatal event; defaults to 1 for exponential
#' @param beta_A_D numeric scalar; treatment coefficient for the fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector length p; coefficients for x on the fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector length q; coefficients for u on the fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard rate for the non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape parameter for the non-fatal event; defaults to 1 for exponential
#' @param beta_A_H numeric scalar; treatment coefficient for the non-fatal event; can also be \code{NULL}
#' @param beta_X_H numeric vector length p; coefficients for x on the non-fatal event; can also be \code{NULL}
#' @param beta_U_H numeric vector length q; coefficients for u on the non-fatal event; can also be \code{NULL}
#' @param lambda_C numeric scalar > 0; baseline hazard rate for random censoring; can also be \code{NULL}, denoting no random censoring
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard copula parameter; defaults to 1 for no association
#' @param phi_admin numeric scalar > 0; administrative censoring time; can also be \code{NULL} for no administrative censoring
#' @param rel.tol numeric scalar; relative tolerance for \code{integrate()}, default 1e-5
#' @param abs.tol numeric scalar; absolute tolerance for \code{integrate()}, default 1e-8
#' @param y2_max numeric scalar > 0; upper limit for integration in y2; defaults to Inf, or truncated at \code{phi_admin} if non-NULL
#' @param fd_step numeric scalar > 0; finite difference step for derivative in y1, default 1e-4
#'
#' @return list with elements \code{WR}, \code{A}, \code{B}, \code{C}, \code{D}
#'
#' @keywords internal
.true_wr_single_covariate <- function(
    x = NULL,
    u = NULL,
    lambda_D,
    kappa_D = 1,
    beta_A_D = NULL,
    beta_X_D = NULL,
    beta_U_D = NULL,
    lambda_H,
    kappa_H = 1,
    beta_A_H = NULL,
    beta_X_H = NULL,
    beta_U_H = NULL,
    lambda_C = NULL,
    beta_A_C = NULL,
    theta_copula = 1,
    phi_admin = NULL,
    rel.tol = 1e-5,
    abs.tol = 1e-8,
    y2_max = Inf,
    fd_step = 1e-4
) {

  # Make sure there exists at least one set of covariates that aren't NULL
  if (is.null(x) && is.null(u)) {
    stop("At least one of 'x' or 'u' must be non-NULL.")
  }
  if (!is.null(x) && !is.numeric(x)) stop("'x' must be numeric or NULL.")
  if (!is.null(u) && !is.numeric(u)) stop("'u' must be numeric or NULL.")

  # Enforce length compatibility with the coefficients, silent if NULL
  if (!is.null(beta_X_D) && !is.null(x) && length(beta_X_D) != length(x)) {
    stop("length(beta_X_D) must match length(x).")
  }
  if (!is.null(beta_X_H) && !is.null(x) && length(beta_X_H) != length(x)) {
    stop("length(beta_X_H) must match length(x).")
  }
  if (!is.null(beta_U_D) && !is.null(u) && length(beta_U_D) != length(u)) {
    stop("length(beta_U_D) must match length(u).")
  }
  if (!is.null(beta_U_H) && !is.null(u) && length(beta_U_H) != length(u)) {
    stop("length(beta_U_H) must match length(u).")
  }

  # Adjust y2_max using phi_admin if provided
  if (!is.null(phi_admin)) {
    if (!is.numeric(phi_admin) || phi_admin <= 0) {
      stop("'phi_admin' must be a numeric scalar > 0 or NULL.")
    }
    if (!is.infinite(y2_max)) {
      y2_max <- min(y2_max, phi_admin)
    } else {
      y2_max <- phi_admin
    }
  }

  # Calculate the four integrals for the true win ratio
  A_int <- integrate(
    Vectorize(function(y2) {
      .Ind_singleNu(
        y2 = y2, x = x, u = u,
        lambda_D = lambda_D, kappa_D = kappa_D,
        beta_A_D = beta_A_D, beta_X_D = beta_X_D, beta_U_D = beta_U_D,
        lambda_H = lambda_H, kappa_H = kappa_H,
        beta_A_H = beta_A_H, beta_X_H = beta_X_H, beta_U_H = beta_U_H,
        lambda_C = lambda_C, beta_A_C = beta_A_C,
        theta_copula = theta_copula,
        fd_step = fd_step
      )
    }),
    lower = 0, upper = y2_max,
    rel.tol = rel.tol, abs.tol = abs.tol
  )

  C_int <- integrate(
    Vectorize(function(y2) {
      .Ind_singleDe(
        y2 = y2, x = x, u = u,
        lambda_D = lambda_D, kappa_D = kappa_D,
        beta_A_D = beta_A_D, beta_X_D = beta_X_D, beta_U_D = beta_U_D,
        lambda_H = lambda_H, kappa_H = kappa_H,
        beta_A_H = beta_A_H, beta_X_H = beta_X_H, beta_U_H = beta_U_H,
        lambda_C = lambda_C, beta_A_C = beta_A_C,
        theta_copula = theta_copula,
        fd_step = fd_step
      )
    }),
    lower = 0, upper = y2_max,
    rel.tol = rel.tol, abs.tol = abs.tol
  )

  B_int <- .double_integral(
    myfun = function(y1, y2) {
      .Ind_doubleNu(
        y1 = y1, y2 = y2, x = x, u = u,
        lambda_D = lambda_D, kappa_D = kappa_D,
        beta_A_D = beta_A_D, beta_X_D = beta_X_D, beta_U_D = beta_U_D,
        lambda_H = lambda_H, kappa_H = kappa_H,
        beta_A_H = beta_A_H, beta_X_H = beta_X_H, beta_U_H = beta_U_H,
        lambda_C = lambda_C, beta_A_C = beta_A_C,
        theta_copula = theta_copula,
        fd_step = fd_step
      )
    },
    y2_max = y2_max,
    rel.tol = rel.tol,
    abs.tol = abs.tol
  )

  D_int <- .double_integral(
    myfun = function(y1, y2) {
      .Ind_doubleDe(
        y1 = y1, y2 = y2, x = x, u = u,
        lambda_D = lambda_D, kappa_D = kappa_D,
        beta_A_D = beta_A_D, beta_X_D = beta_X_D, beta_U_D = beta_U_D,
        lambda_H = lambda_H, kappa_H = kappa_H,
        beta_A_H = beta_A_H, beta_X_H = beta_X_H, beta_U_H = beta_U_H,
        lambda_C = lambda_C, beta_A_C = beta_A_C,
        theta_copula = theta_copula,
        fd_step = fd_step
      )
    },
    y2_max = y2_max,
    rel.tol = rel.tol,
    abs.tol = abs.tol
  )

  # Grab the results
  A <- A_int$value
  B <- B_int$value
  C <- C_int$value
  D <- D_int$value

  # Put them together to get the true win ratio
  WR <- (A + B) / (C + D)

  return(
    list(
      WR = WR,
      A  = A,
      B  = B,
      C  = C,
      D  = D
    )
  )

}

#' True causal win ratio averaged over covariate patterns
#'
#' Compute the true causal win ratio by averaging A(x,u)+B(x,u) and C(x,u)+D(x,u) over covariate patterns (x,u).
#' The averaging is done over the rows of X_mat and U_mat, treating the empirical distribution of (X,U) as the target distribution.
#'
#' @param X_mat numeric matrix N x p; measured covariates; can be NULL if there are no measured covariates
#' @param U_mat numeric matrix N x q; unmeasured covariates; can be NULL if there are no unmeasured covariates
#' @param lambda_D numeric scalar > 0; baseline hazard rate for the fatal event
#' @param kappa_D numeric scalar > 0; Weibull shape parameter for the fatal event; defaults to 1 for exponential
#' @param beta_A_D numeric scalar; treatment coefficient for the fatal event; can also be \code{NULL}
#' @param beta_X_D numeric vector length p; coefficients for x on the fatal event; can also be \code{NULL}
#' @param beta_U_D numeric vector length q; coefficients for u on the fatal event; can also be \code{NULL}
#' @param lambda_H numeric scalar > 0; baseline hazard rate for the non-fatal event
#' @param kappa_H numeric scalar > 0; Weibull shape parameter for the non-fatal event; defaults to 1 for exponential
#' @param beta_A_H numeric scalar; treatment coefficient for the non-fatal event; can also be \code{NULL}
#' @param beta_X_H numeric vector length p; coefficients for x on the non-fatal event; can also be \code{NULL}
#' @param beta_U_H numeric vector length q; coefficients for u on the non-fatal event; can also be \code{NULL}
#' @param lambda_C numeric scalar > 0; baseline hazard rate for random censoring; can also be \code{NULL} for no random censoring
#' @param beta_A_C numeric scalar; treatment coefficient for censoring; can also be \code{NULL}
#' @param theta_copula numeric scalar >= 1; Gumbel–Hougaard copula parameter (default 1 for independence)
#' @param phi_admin numeric scalar > 0; administrative censoring time; can also be \code{NULL} for no administrative censoring
#' @param rel.tol numeric scalar; relative tolerance for \code{integrate()}, default 1e-5
#' @param abs.tol numeric scalar; absolute tolerance for \code{integrate()}, default 1e-8
#' @param y2_max numeric scalar > 0; upper limit for integration in y2; defaults to Inf, or truncated at \code{phi_admin} if non-NULL
#' @param fd_step numeric scalar > 0; finite difference step for derivative in y1, default 1e-4
#'
#' @return list with elements:
#'   \describe{
#'     \item{WR}{estimated true causal win ratio}
#'     \item{num}{Monte Carlo estimate of E[A(X,U)+B(X,U)]}
#'     \item{den}{Monte Carlo estimate of E[C(X,U)+D(X,U)]}
#'     \item{details}{list with vectors A, B, C, D for each row}
#'   }
#'
#' @export
true_wr_analytic <- function(
    X_mat = NULL,
    U_mat = NULL,
    lambda_D,
    kappa_D = 1,
    beta_A_D = NULL,
    beta_X_D = NULL,
    beta_U_D = NULL,
    lambda_H,
    kappa_H = 1,
    beta_A_H = NULL,
    beta_X_H = NULL,
    beta_U_H = NULL,
    lambda_C = NULL,
    beta_A_C = NULL,
    theta_copula = 1,
    phi_admin = NULL,
    rel.tol = 1e-5,
    abs.tol = 1e-8,
    y2_max = Inf,
    fd_step = 1e-4
) {

  # Require at least one of X_mat or U_mat to be non-NULL
  if (is.null(X_mat) && is.null(U_mat)) {
    stop("At least one of X_mat or U_mat must be non-NULL.")
  }
  if (!is.null(X_mat)) {
    X_mat <- as.matrix(X_mat)
  }
  if (!is.null(U_mat)) {
    U_mat <- as.matrix(U_mat)
  }

  # Enforce consistency in dimensions and number of samples N
  if (!is.null(X_mat) && !is.null(U_mat)) {
    if (nrow(X_mat) != nrow(U_mat)) {
      stop("X_mat and U_mat must have the same number of rows if both are non-NULL.")
    }
    N <- nrow(X_mat)
  } else if (!is.null(X_mat)) {
    N <- nrow(X_mat)
  } else {
    N <- nrow(U_mat)
  }

  # Result vectors storing the integrals for individual covaraites
  # We will average these at the end
  res_A <- numeric(N)
  res_B <- numeric(N)
  res_C <- numeric(N)
  res_D <- numeric(N)

  # For each individual
  for (i in seq_len(N)) {
    # Grab their covariates
    x_i <- ifelse(!is.null(X_mat), X_mat[i, ], NULL)
    u_i <- ifelse(!is.null(U_mat), U_mat[i, ], NULL)

    # Calculate the true causal WR integral given this individual's covariates
    out_i <- .true_wr_single_covariate(
      x = x_i, u = u_i,
      lambda_D = lambda_D,
      kappa_D  = kappa_D,
      beta_A_D = beta_A_D,
      beta_X_D = beta_X_D,
      beta_U_D = beta_U_D,
      lambda_H = lambda_H,
      kappa_H  = kappa_H,
      beta_A_H = beta_A_H,
      beta_X_H = beta_X_H,
      beta_U_H = beta_U_H,
      lambda_C = lambda_C,
      beta_A_C = beta_A_C,
      theta_copula = theta_copula,
      phi_admin = phi_admin,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      y2_max  = y2_max,
      fd_step = fd_step
    )

    # Store the results
    res_A[i] <- out_i$A
    res_B[i] <- out_i$B
    res_C[i] <- out_i$C
    res_D[i] <- out_i$D
  }

  # Average out the numerator and denominator, completing the expectation
  num <- mean(res_A + res_B)
  den <- mean(res_C + res_D)
  # Calculate the true win ratio
  WR  <- num / den
  return(
    list(
      WR  = WR,
      num = num,
      den = den,
      details = list(A = res_A, B = res_B, C = res_C, D = res_D)
    )
  )
}
