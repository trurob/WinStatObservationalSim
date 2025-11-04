test_that("Number of wins, losses, and ties from get_true_WR sum to total pairs", {
  set.seed(2)
  dat <- simulate_dataset(
    N = 400,
    p = 0, q = 0,
    mu_X = numeric(0), sd_X = numeric(0),
    mu_U = numeric(0), sd_U = numeric(0),
    treat_assign = "randomized",
    lambda_H = 0.08, kappa_H = 1, beta_A_H = 0,
    beta_X_H = numeric(0), beta_U_H = numeric(0),
    lambda_D = 0.12, kappa_D = 1, beta_A_D = 0,
    beta_X_D = numeric(0), beta_U_D = numeric(0),
    theta_copula = 1,
    lambda_C = 0, beta_A_C = 0,
    phi_admin = 8
  )
  treatment <- subset(dat, A == 1, select = c(Y_D, delta_D, Y_H, delta_H))
  control <- subset(dat, A == 0, select = c(Y_D, delta_D, Y_H, delta_H))
  result <- get_true_WR(treatment, control)
  expect_equal(result$treated_wins + result$control_wins + result$ties, result$total_pairs)
  expect_gte(result$pi_T, 0)
  expect_gte(result$pi_C, 0)
})
