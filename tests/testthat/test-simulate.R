test_that("simulate_dataset returns expected output object; case: no covariates, no random censoring, no association between events", {
  set.seed(1)
  N <- 1000
  dat <- simulate_dataset(
    N = N,
    p = 0, q = 0,
    mu_X = numeric(0), sd_X = numeric(0),
    mu_U = numeric(0), sd_U = numeric(0),
    treat_assign = "randomized",
    lambda_H = 0.1, kappa_H = 1, beta_A_H = 0,
    beta_X_H = numeric(0), beta_U_H = numeric(0),
    lambda_D = 0.2, kappa_D = 1, beta_A_D = 0,
    beta_X_D = numeric(0), beta_U_D = numeric(0),
    theta_copula = 1,
    lambda_C = 0, beta_A_C = 0,
    phi_admin = 10,
    incl_tfe = TRUE,
    incl_latent = TRUE
  )
  # returns expected class
  expect_s3_class(dat, "data.frame")
  # returns expected number of rows
  expect_equal(nrow(dat), N)
  # returns expected variables
  expect_true(all(c("A","Y_D","delta_D","Y_H","delta_H","Y_tfe","delta_tfe","T_D","T_H") %in% names(dat)))
  # Semi-competing risk paradigm in tact
  expect_true(all(dat$Y_H <= pmin(dat$Y_D, dat$Y_tfe)))
  # Latent event times happen before true event times
  expect_true(all(dat$T_D >= pmin(dat$Y_D, dat$Y_tfe)))
  expect_true(all(dat$T_H >= pmin(dat$Y_D, dat$Y_H, dat$Y_tfe)))
})
