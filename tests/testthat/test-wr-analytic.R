test_that("Weibull helpers reduce to exponential when kappa = 1", {

  a      <- 1
  x      <- NULL
  u      <- NULL
  lambda <- 0.2
  beta_A <- 0.3
  t      <- 2

  # Scale parameter: Lambda = lambda * exp(-beta_A * a)
  lam_scale <- .Lambda_D(
    a      = a,
    x      = x,
    u      = u,
    lambda_D = lambda,
    beta_A_D = beta_A,
    beta_X_D = NULL,
    beta_U_D = NULL
  )
  expect_equal(lam_scale, lambda * exp(-beta_A * a))

  # Survival function for death: S_D(t) = exp(-(Lambda * t)^kappa) with kappa = 1
  S_val <- .S_D(
    t       = t,
    a       = a,
    x       = x,
    u       = u,
    lambda_D = lambda,
    kappa_D  = 1,
    beta_A_D = beta_A,
    beta_X_D = NULL,
    beta_U_D = NULL
  )
  expect_equal(S_val, exp(-lam_scale * t))

  # Density for time to death: f_D(t) = Lambda * exp(-Lambda * t) when kappa = 1
  f_val <- .f_D(
    y2      = t,
    a       = a,
    x       = x,
    u       = u,
    lambda_D = lambda,
    kappa_D  = 1,
    beta_A_D = beta_A,
    beta_X_D = NULL,
    beta_U_D = NULL
  )
  expect_equal(f_val, lam_scale * exp(-lam_scale * t), tolerance = 1e-10)

  # Hazard for death: lambda_D_fun = f / S
  h_val <- .lambda_D_fun(
    y2      = t,
    a       = a,
    x       = x,
    u       = u,
    lambda_D = lambda,
    kappa_D  = 1,
    beta_A_D = beta_A,
    beta_X_D = NULL,
    beta_U_D = NULL
  )
  expect_equal(h_val, lam_scale, tolerance = 1e-10)

})

test_that("Joint survival reduces to product under independence via theta_copula = 1", {

  a <- 0
  x <- NULL
  u <- NULL

  lambda_D <- 0.1
  lambda_H <- 0.2
  kappa_D  <- 1.5
  kappa_H  <- 0.7

  y2 <- 1.1
  y1 <- 0.8

  # Survival function for death
  SD <- .S_D(
    t        = y2,
    a        = a,
    x        = x,
    u        = u,
    lambda_D = lambda_D,
    kappa_D  = kappa_D,
    beta_A_D = 0,
    beta_X_D = NULL,
    beta_U_D = NULL
  )

  # Survival function for hospitalization
  SH <- .S_H(
    t        = y1,
    a        = a,
    x        = x,
    u        = u,
    lambda_H = lambda_H,
    kappa_H  = kappa_H,
    beta_A_H = 0,
    beta_X_H = NULL,
    beta_U_H = NULL
  )

  # Joint survival function of death and hospitalization
  SDH <- .S_DH(
    y2        = y2,
    y1        = y1,
    a         = a,
    x         = x,
    u         = u,
    lambda_D  = lambda_D,
    kappa_D   = kappa_D,
    beta_A_D  = 0,
    beta_X_D  = NULL,
    beta_U_D  = NULL,
    lambda_H  = lambda_H,
    kappa_H   = kappa_H,
    beta_A_H  = 0,
    beta_X_H  = NULL,
    beta_U_H  = NULL,
    theta_copula = 1
  )
  expect_equal(SDH, SD * SH, tolerance = 1e-10)

})

test_that("No random censoring makes B(x,u) and D(x,u) vanish", {

  x <- 0
  u <- NULL

  lambda_D <- 0.08
  lambda_H <- 0.10

  out <- .true_wr_single_covariate(
    x = x,
    u = u,
    lambda_D = lambda_D,
    kappa_D  = 1,
    beta_A_D = 0.2,
    beta_X_D = 0,
    beta_U_D = NULL,
    lambda_H = lambda_H,
    kappa_H  = 1,
    beta_A_H = 0.3,
    beta_X_H = 0,
    beta_U_H = NULL,
    # No random censoring
    lambda_C = NULL,
    beta_A_C = NULL,
    theta_copula = 1,
    phi_admin = 5,
    y2_max   = Inf,
    rel.tol  = 1e-5,
    abs.tol  = 1e-8,
    fd_step  = 1e-4
  )
  expect_equal(out$B, 0, tolerance = 1e-8)
  expect_equal(out$D, 0, tolerance = 1e-8)

})

test_that("true_wr_analytic averages over identical covariate patterns correctly", {

  set.seed(1)

  # Single measured covariate value, no unmeasured covariates
  x_val <- 0.5
  X_mat <- matrix(x_val, nrow = 50, ncol = 1)
  U_mat <- NULL

  lambda_D <- 0.08
  lambda_H <- 0.10

  # Single-covariate WR
  single <- .true_wr_single_covariate(
    x = x_val,
    u = NULL,
    lambda_D = lambda_D,
    kappa_D  = 1,
    beta_A_D = 0.2,
    beta_X_D = 0.1,
    beta_U_D = NULL,
    lambda_H = lambda_H,
    kappa_H  = 1,
    beta_A_H = 0.3,
    beta_X_H = 0.05,
    beta_U_H = NULL,
    lambda_C = 0.09,
    beta_A_C = 0.1,
    theta_copula = 1,
    phi_admin = 5,
    y2_max   = Inf,
    rel.tol  = 1e-5,
    abs.tol  = 1e-8,
    fd_step  = 1e-4
  )

  # Averaged WR over identical rows
  avg <- true_wr_analytic(
    X_mat = X_mat,
    U_mat = U_mat,
    lambda_D = lambda_D,
    kappa_D  = 1,
    beta_A_D = 0.2,
    beta_X_D = 0.1,
    beta_U_D = NULL,
    lambda_H = lambda_H,
    kappa_H  = 1,
    beta_A_H = 0.3,
    beta_X_H = 0.05,
    beta_U_H = NULL,
    lambda_C = 0.09,
    beta_A_C = 0.1,
    theta_copula = 1,
    phi_admin = 5,
    y2_max     = Inf,
    rel.tol    = 1e-5,
    abs.tol    = 1e-8,
    fd_step    = 1e-4
  )

  expect_equal(avg$WR, single$WR, tolerance = 1e-4)
})


