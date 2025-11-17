### Test helper functions
test_that(".logistic_function maps to (0,1)", {
  x <- c(-10, -1, 0, 1, 10)
  p <- .logistic_function(x)
  expect_true(all(p > 0 & p < 1))
  expect_equal(p[3], 0.5)
})
test_that(".get_linear_predictor works with a NULL matrix", {
  n <- 10L
  lp_null <- .get_linear_predictor(NULL, coef = numeric(0), n = n)
  expect_length(lp_null, n)
  expect_true(all(lp_null == 0))
})
test_that(".get_linear_predictor returns expected number of covariates", {
  n <- 10L
  set.seed(1)
  X <- matrix(rnorm(n * 2), ncol = 2)
  beta <- c(0.5, -0.25)
  lp <- .get_linear_predictor(X, beta, n)
  expect_length(lp, n)
  expect_equal(lp, drop(X %*% beta))
})
### Test `.generate_covariate_matrix`
test_that(".generate_covariate_matrix handles p being NULL or 0", {
  M1 <- .generate_covariate_matrix(N = 100, p = NULL)
  M2 <- .generate_covariate_matrix(N = 100, p = 0)
  expect_null(M1)
  expect_null(M2)
})
test_that(".generate_covariate_matrix generates Normal covariates by default", {
  set.seed(1)
  N <- 1000
  p <- 2
  mu <- c(0, 5)
  sd <- c(1, 2)

  X <- .generate_covariate_matrix(
    N = N, p = p,
    mu = mu, sd = sd,
    dist = NULL, prob = NULL,
    prefix = "X"
  )
  # We get expected dimensions and column names
  expect_equal(dim(X), c(N, p))
  expect_equal(colnames(X), c("X_1", "X_2"))
  # Sample means of the respective covariates should be close to the population mean
  expect_equal(mean(X[, 1]), mu[1], tolerance = 0.1)
  expect_equal(mean(X[, 2]), mu[2], tolerance = 0.1)
})
test_that(".generate_covariate_matrix generates expected number of Bernoulli covariates", {
  set.seed(1)
  N <- 5000
  p <- 3
  dist <- c("Bernoulli", "Bernoulli", "Bernoulli")
  prob <- c(0.2, 0.5, 0.8)
  X <- .generate_covariate_matrix(
    N = N, p = p,
    mu = NULL, sd = NULL,
    dist = dist, prob = prob,
    prefix = "X"
  )
  expect_equal(dim(X), c(N, p))
  # Generated sample proportions should be close to prob
  prop <- colMeans(X)
  expect_equal(unname(prop), prob, tolerance = 0.02)
})
test_that(".generate_covariate_matrix handles a vector composed of both Normal and Bernoulli random variables", {
  set.seed(1)
  N <- 5000
  p <- 3
  dist <- c("Normal", "Bernoulli", "Normal")
  mu <- c(10, NA, 0)
  sd <- c(2,  NA, 1)
  prob <- c(NA, 0.3, NA)

  X <- .generate_covariate_matrix(
    N = N, p = p,
    mu = mu, sd = sd,
    dist = dist, prob = prob,
    prefix = "X"
  )
  expect_equal(dim(X), c(N, p))
  # Check Normal columns
  expect_equal(mean(X[, 1]), 10, tolerance = 0.15)
  expect_equal(mean(X[, 3]), 0, tolerance = 0.15)
  # Check Bernoulli column
  expect_equal(mean(X[, 2]), 0.3, tolerance = 0.02)
})
### TODO: Investigate the Bernoulli prob error
#test_that(".generate_covariate_matrix throws informative errors", {
#  N <- 10
#  p <- 2
#  # dist length mismatch
#  expect_error(
#    .generate_covariate_matrix(N, p, dist = "Normal"),
#    "Length of 'dist'"
#  )
#  # Bernoulli with prob NULL
#  expect_error(
#    .generate_covariate_matrix(
#      N, p,
#      dist = c("Bernoulli", "Normal"),
#      mu = c(0, 0), sd = c(1, 1),
#      prob = NULL
#    ),
#    "prob is NULL"
#  )
#  # Normal with mu and sd NULL
#  expect_error(
#    .generate_covariate_matrix(
#      N, p,
#      dist = c("Normal", "Normal"),
#      mu = NULL, sd = NULL,
#      prob = NULL
#    ),
#    "mu or 'sd' is NULL"
#  )
#  # Bernoulli with invalid probability
#  expect_error(
#    .generate_covariate_matrix(
#      N, p,
#      dist = c("Bernoulli", "Normal"),
#      mu = c(0, 0), sd = c(1, 1),
#      prob = c(-0.1, NA)
#    ),
#    "requires a probability in \\(0, 1\\)"
#  )
#})
### Test `.simulate_dataset`
test_that("simulate_dataset works with no covariates and no censoring", {
  set.seed(1)
  N <- 2000
  dat <- simulate_dataset(
    N = N,
    lambda_H = 0.1, kappa_H = 1.2, beta_A_H = -0.3,
    lambda_D = 0.08, kappa_D = 1.1, beta_A_D = -0.25
  )
  # Expected object type and dimensions
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), N)
  # Columns present
  expect_true(all(c("id", "A", "Y_D", "delta_D", "Y_H", "delta_H") %in% names(dat)))
  # No covariates
  expect_false(any(grepl("^X_", names(dat))))
  expect_false(any(grepl("^U_", names(dat))))
  # All times positive
  expect_true(all(dat$Y_D > 0))
  expect_true(all(dat$Y_H > 0))
  # Some events should occur (not all censored)
  expect_true(any(dat$delta_D == 1))
  expect_true(any(dat$delta_H == 1))
})
test_that("simulate_dataset respects semi-competing risk structure", {
  set.seed(1)
  N <- 2000
  dat <- simulate_dataset(
    N = N,
    lambda_H = 0.1, kappa_H = 1.2, beta_A_H = -0.3,
    lambda_D = 0.08, kappa_D = 1.1, beta_A_D = -0.25,
    phi_admin = 10,
    incl_latent = TRUE
  )
  # latent times T_D and T_H present in final columns
  expect_true(all(c("T_D", "T_H") %in% names(dat)))

  # Check semi-competing risks: Y_H = min(T_H, T_D, C, phi)
  # We can't see C here, but we can check certain cases.
  # Case 1: non-fatal observed (delta_H == 1) and fatal observed later (delta_D == 1)
  idx <- which(dat$delta_H == 1 & dat$delta_D == 1)
  # For these, T_H <= T_D, and Y_H == T_H
  expect_true(all(dat$T_H[idx] <= dat$T_D[idx]))
  expect_true(all(abs(dat$Y_H[idx] - dat$T_H[idx]) < 1e-8))
  # Case 2: non-fatal observed (delta_H == 1) and fatal not observed (delta_D == 0)
  idx2 <- which(dat$delta_H == 1 & dat$delta_D == 0)
  # For these, Y_H should still equal T_H
  expect_true(all(abs(dat$Y_H[idx2] - dat$T_H[idx2]) < 1e-8))
})
test_that("no censoring when lambda_C, beta_A_C, phi_admin are NULL", {
  set.seed(1)
  N <- 1000
  dat <- simulate_dataset(
    N = N,
    lambda_H = 0.1, kappa_H = 1.2, beta_A_H = -0.3,
    lambda_D = 0.08, kappa_D = 1.1, beta_A_D = -0.25,
    incl_latent = TRUE
  )
  # With no censoring, Y_D == T_D and Y_H == min(T_H, T_D)
  expect_true(all(abs(dat$Y_D - dat$T_D) < 1e-8))
  # For those where T_H <= T_D, Y_H == T_H; otherwise Y_H == T_D
  idx_H_first <- which(dat$T_H <= dat$T_D)
  idx_D_first <- which(dat$T_H > dat$T_D)
  expect_true(all(abs(dat$Y_H[idx_H_first] - dat$T_H[idx_H_first]) < 1e-8))
  expect_true(all(abs(dat$Y_H[idx_D_first] - dat$T_D[idx_D_first]) < 1e-8))
})
test_that("administrative censoring truncates follow-up at phi_admin", {
  set.seed(1)
  N <- 1000
  dat <- simulate_dataset(
    N = N,
    lambda_H = 0.1, kappa_H = 1.2, beta_A_H = -0.3,
    lambda_D = 0.08, kappa_D = 1.1, beta_A_D = -0.25,
    phi_admin = 5,
    incl_latent = TRUE
  )
  expect_true(all(dat$Y_D <= 5 + 1e-8))
  expect_true(all(dat$Y_H <= 5 + 1e-8))
})
test_that("random censoring works and yields some censoring", {
  set.seed(1)
  N <- 1000
  dat <- simulate_dataset(
    N = N,
    lambda_H = 0.1, kappa_H = 1.2, beta_A_H = -0.3,
    lambda_D = 0.08, kappa_D = 1.1, beta_A_D = -0.25,
    lambda_C = 0.05, beta_A_C = 0,
    phi_admin = 20
  )
  # We should see some censored observations
  expect_true(any(dat$delta_D == 0))
  expect_true(any(dat$delta_H == 0))
})
test_that("positive beta_A_D increases fatal event times in treatment arm", {
  set.seed(1)
  N <- 5000
  dat <- simulate_dataset(
    N = N,
    lambda_H = 0.1, kappa_H = 1.2, beta_A_H = 0,
    lambda_D = 0.08, kappa_D = 1.1, beta_A_D = 0.5,
    phi_admin = NULL,
    incl_latent = TRUE
  )
  # Look at latent fatal times by treatment arm
  Td_treat <- dat$T_D[dat$A == 1]
  Td_ctrl  <- dat$T_D[dat$A == 0]
  # Treatment arm should have higher median time to fatal event
  expect_gt(median(Td_treat), median(Td_ctrl))
})
test_that("simulate_dataset errors when X_dist Bernoulli but prob_X is NULL", {
  expect_error(
    simulate_dataset(
      N = 100,
      p = 1,
      X_dist = "Bernoulli",
      mu_X = NULL, sd_X = NULL,
      prob_X = NULL,
      lambda_H = 0.1, kappa_H = 1.2, beta_A_H = 0,
      lambda_D = 0.08, kappa_D = 1.1, beta_A_D = 0
    ),
    "dist == 'Bernoulli', but 'prob' is NULL",
    fixed = TRUE
  )
})
### TODO: expect_error is not catching the correct error, fix this
#test_that("simulate_dataset errors when X_dist Normal but mu_X and sd_X are NULL", {
#  expect_error(
#    simulate_dataset(
#      N = 100,
#      p = 1,
#      X_dist = "Normal",
#      mu_X = NULL, sd_X = NULL,
#      lambda_H = 0.1, kappa_H = 1.2, beta_A_H = 0,
#      lambda_D = 0.08, kappa_D = 1.1, beta_A_D = 0
#    ),
#    "mu or 'sd' is NULL",
#    fixed = TRUE
#  )
#})


