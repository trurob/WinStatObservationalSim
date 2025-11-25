# This is a modified version of the `win.strategy.default.R` file from the
# WINS package (under terms from the GPLv3 license)

##############################################
# ------------ HELPER FUNCTIONS ------------ #
##############################################
### TODO: move helper functions into their own file

#' Count positions of query times relative to a sorted control time vector
#'
#' @description
#' Given a sorted numeric vector \code{ctrl_sorted} (strictly non-decreasing)
#' and a numeric query vector \code{q}, compute, for each \code{q[i]},
#' how many control times are strictly less than \code{q[i]} (\code{lt}),
#' less than or equal to \code{q[i]} (\code{le}), equal to \code{q[i]} (\code{eq}),
#' and strictly greater than \code{q[i]} (\code{gt}).
#'
#' @param ctrl_sorted A numeric vector of control times in non-decreasing order.
#' @param q A numeric vector of query times.
#'
#' @return A list with numeric vectors \code{lt}, \code{le}, \code{eq}, \code{gt}
#'   of the same length as \code{q}. If \code{length(ctrl_sorted) == 0},
#'   all outputs are zero vectors.
#'
#' @keywords internal
.count_positions <- function(ctrl_sorted, q) {
  m <- length(ctrl_sorted)

  if (m == 0L) {
    zeros <- rep.int(0L, length(q))
    return(list(lt = zeros, le = zeros, eq = zeros, gt = zeros))
  }

  # lower_bound (strictly less): left.open = TRUE
  lt <- findInterval(q, ctrl_sorted, left.open = TRUE)
  # upper_bound (≤): left.open = FALSE
  le <- findInterval(q, ctrl_sorted, left.open = FALSE)

  eq <- le - lt
  gt <- m - le

  list(
    lt = lt,
    le = le,
    eq = eq,
    gt = gt
  )
}

#' Preprocess control arm for death stage
#'
#' @description
#' Split control subjects into those with an observed death and those without.
#' Return the sorted death times and the counts needed by the counting kernel.
#'
#' @param Y_D_c Numeric vector of observed death times in control.
#' @param delta_D_c Integer (0/1) vector of death indicators in control.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{d_ctrl_sorted}: sorted death times among controls with death.
#'     \item \code{m1}: count of controls with death.
#'     \item \code{m0}: count of controls without death.
#'     \item \code{c0_idx}: integer indices of controls without death.
#'   }
#'
#' @keywords internal
.prep_control_death <- function(Y_D_c, delta_D_c) {
  c1_idx <- which(delta_D_c == 1L)
  c0_idx <- which(delta_D_c == 0L)

  list(
    d_ctrl_sorted = sort(Y_D_c[c1_idx]),
    m1 = length(c1_idx),
    m0 = length(c0_idx),
    c0_idx = c0_idx
  )
}

#' Preprocess control arm for hospitalization stage (restricted to no-death)
#'
#' @description
#' Within the subset of controls without observed death, split by hospitalization
#' status and return sorted hospitalization times for those hospitalized.
#'
#' @param Y_H_c0 Numeric vector of observed hospitalization times among controls with no death.
#' @param delta_H_c0 Integer (0/1) vector of hospitalization indicators among controls with no death.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{h_ctrl_sorted}: sorted hospitalization times among controls w/ hospitalization.
#'     \item \code{m0h1}: count of controls (no-death) with hospitalization.
#'     \item \code{m0h0}: count of controls (no-death) without hospitalization.
#'   }
#'
#' @keywords internal
.prep_control_hosp_nodeath <- function(Y_H_c0, delta_H_c0) {
  c0h1_idx <- which(delta_H_c0 == 1L)
  c0h0_idx <- which(delta_H_c0 == 0L)

  list(
    h_ctrl_sorted = sort(Y_H_c0[c0h1_idx]),
    m0h1 = length(c0h1_idx),
    m0h0 = length(c0h0_idx)
  )
}

#' Count wins/losses/ties at the death stage
#'
#' @description
#' Implements the exact all-against-all death comparison by counting, without
#' enumerating pairs. Equal times are counted as ties.
#'
#' @param Y_D_t Numeric vector of observed death times in treated.
#' @param delta_D_t Integer (0/1) vector of death indicators in treated.
#' @param d_ctrl_sorted Sorted control death times (from \code{.prep_control_death()}).
#' @param m1 Number of controls with observed death.
#' @param m0 Number of controls without observed death.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{wins}: treated wins resolved at death stage.
#'     \item \code{losses}: treated losses resolved at death stage.
#'     \item \code{ties}: ties at death stage (equal times).
#'     \item \code{t0_idx}: indices of treated with no observed death (for stage B).
#'   }
#'
#' @keywords internal
.count_stage_death <- function(Y_D_t, delta_D_t, d_ctrl_sorted, m1, m0) {
  t1_idx <- which(delta_D_t == 1L)  # treated w/ death
  t0_idx <- which(delta_D_t == 0L)  # treated w/o death

  wins <- 0.0
  losses <- 0.0
  ties <- 0.0

  # Treated with NO death: win vs all controls who died.
  n0 <- length(t0_idx)
  if (n0 > 0L && m1 > 0L) {
    wins <- wins + as.double(n0) * as.double(m1)
  }

  # Treated with death: compare ordering vs controls with death, and lose vs controls with no death.
  n1 <- length(t1_idx)
  if (n1 > 0L) {
    q <- Y_D_t[t1_idx]
    pos <- .count_positions(d_ctrl_sorted, q)

    # Wins: control's death strictly earlier than treated (lt)
    wins <- wins + sum(as.double(pos$lt))

    # Losses: control's death strictly later than treated (gt) and all controls with no death (m0)
    losses <- losses + sum(as.double(pos$gt)) + as.double(n1) * as.double(m0)

    # Ties at death (equal times)
    ties <- ties + sum(as.double(pos$eq))
  }

  list(
    wins = wins,
    losses = losses,
    ties = ties,
    t0_idx = t0_idx
  )
}

#' Count wins/losses/ties at the hospitalization stage within no-death pairs
#'
#' @description
#' Restrict to pairs where neither side has an observed death, and implement
#' the hospitalization comparison by counting. Equal times are counted as ties.
#'
#' @param Y_H_t0 Numeric vector of treated hospitalization times among no-death treated.
#' @param delta_H_t0 Integer (0/1) vector of treated hospitalization indicators among no-death treated.
#' @param h_ctrl_sorted Sorted hospitalization times among controls with no death and with hospitalization.
#' @param m0h1 Number of controls (no-death) with hospitalization.
#' @param m0h0 Number of controls (no-death) without hospitalization.
#'
#' @return A list with \code{wins}, \code{losses}, \code{ties} resolved at the hospitalization stage.
#'
#' @keywords internal
.count_stage_hosp <- function(Y_H_t0, delta_H_t0, h_ctrl_sorted, m0h1, m0h0) {
  t0h1_idx <- which(delta_H_t0 == 1L)  # treated hospitalized (within no-death)
  t0h0_idx <- which(delta_H_t0 == 0L)  # treated not hospitalized (within no-death)

  wins <- 0.0
  losses <- 0.0
  ties <- 0.0

  # Treated NOT hospitalized: win vs all controls who are hospitalized; tie vs controls not hospitalized.
  n0h0 <- length(t0h0_idx)
  if (n0h0 > 0L) {
    if (m0h1 > 0L) {
      wins <- wins + as.double(n0h0) * as.double(m0h1)
    }
    if (m0h0 > 0L) {
      ties <- ties + as.double(n0h0) * as.double(m0h0)
    }
  }

  # Treated hospitalized: compare ordering vs control hospitalized, and lose vs controls not hospitalized.
  n0h1 <- length(t0h1_idx)
  if (n0h1 > 0L) {
    q <- Y_H_t0[t0h1_idx]
    pos <- .count_positions(h_ctrl_sorted, q)

    # Wins: controls hospitalized strictly earlier than treated (lt)
    wins <- wins + sum(as.double(pos$lt))

    # Losses: Controls hospitalized strictly later (gt) and controls never hospitalized
    losses <- losses + sum(as.double(pos$gt)) + as.double(n0h1) * as.double(m0h0)

    # Ties at hospitalization (equal times)
    ties <- ties + sum(as.double(pos$eq))
  }

  list(
    wins = wins,
    losses = losses,
    ties = ties
  )
}

############################################
# ------------ MAIN FUNCTIONS ------------ #
############################################

#' Compute the "true" win ratio on a superpopulation.
#'
#' @description
#' After splitting the superpopulation into respective treatment and control `data.frame` objects,
#' computes the "population" win probabilities and the respective win ratio.
#'
#' @details
#' Uses the code infrastructure developed in the `WINS` package under the GPLv3 license.
#' Uses an all-against-all pairwise comparison approach against a fatal/non-fatal composite hierarchical outcome to determine winners/losers/and ties. Equal times are counted as ties.
#'
#' @param treated `data.frame`; containing the columns: `Y_D`, `delta_D`, `Y_H`, `delta_H`
#'   for the treated arm.
#' @param control `data frame`; same columns as `treated`
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{treated_wins}, \code{control_wins}, \code{ties}
#'     \item \code{total_pairs}
#'     \item \code{pi_T}, \code{pi_C}
#'     \item \code{WR} = \code{pi_T} / \code{pi_C}
#'   }
#'
#' @export
get_true_WR <- function(treated, control) {

  # Grab relevant data, save into vectors, and grab sizes
  Y_D_t <- as.numeric(treated$Y_D)
  delta_D_t <- as.integer(treated$delta_D)
  Y_H_t <- as.numeric(treated$Y_H)
  delta_H_t <- as.integer(treated$delta_H)
  Y_D_c <- as.numeric(control$Y_D)
  delta_D_c <- as.integer(control$delta_D)
  Y_H_c <- as.numeric(control$Y_H)
  delta_H_c <- as.integer(control$delta_H)
  n_t <- length(Y_D_t)
  n_c <- length(Y_D_c)
  total_pairs <- as.double(n_t) * as.double(n_c)

  #---------------------------
  # 1. COMPARE FATAL
  #---------------------------
  ctrl_death <- .prep_control_death(Y_D_c, delta_D_c)

  death_counts <- .count_stage_death(
    Y_D_t = Y_D_t,
    delta_D_t = delta_D_t,
    d_ctrl_sorted = ctrl_death$d_ctrl_sorted,
    m1 = ctrl_death$m1,
    m0 = ctrl_death$m0
  )

  # Subsets that proceed to hospitalization: no-death on both sides
  t0_idx <- death_counts$t0_idx
  c0_idx <- ctrl_death$c0_idx

  Y_H_t0 <- Y_H_t[t0_idx]
  delta_H_t0 <- delta_H_t[t0_idx]
  Y_H_c0 <- Y_H_c[c0_idx]
  delta_H_c0 <- delta_H_c[c0_idx]

  #---------------------------
  # 2. COMPARE NON-FATAL
  #---------------------------
  ctrl_hosp <- .prep_control_hosp_nodeath(Y_H_c0, delta_H_c0)

  hosp_counts <- .count_stage_hosp(
    Y_H_t0 = Y_H_t0,
    delta_H_t0 = delta_H_t0,
    h_ctrl_sorted = ctrl_hosp$h_ctrl_sorted,
    m0h1 = ctrl_hosp$m0h1,
    m0h0 = ctrl_hosp$m0h0
  )

  #---------------------------
  # 3. COMBINE RESULTS
  #---------------------------

  treated_wins <- death_counts$wins + hosp_counts$wins
  control_wins <- death_counts$losses + hosp_counts$losses
  ties <- death_counts$ties + hosp_counts$ties

  pi_T <- treated_wins / total_pairs
  pi_C <- control_wins / total_pairs
  wr <- pi_T / max(pi_C, .Machine$double.eps)

  list(
    treated_wins = treated_wins,
    control_wins = control_wins,
    ties = ties,
    total_pairs = total_pairs,
    pi_T = pi_T,
    pi_C = pi_C,
    WR = wr
  )
}

#' Monte Carlo approximation of the "true" causal win ratio as defined in Zhang 2021
#'
#' @description
#' Generates a large superpopulation under the specified data model,
#' and computes the empirical population-level win ratio using the observed outcomes.
#'
#' @param N_super integer; size of the superpopulation (default \code{1e6})
#' @param seed integer or \code{NULL}; random seed
#' @param ... Additional parameters passed on to \code{simulate_dataset()}.
#' See documentation of \code{simulate_dataset()} for a list of possible arguments.
#'
#' @details
#' This function computes the hierarchical WR on the observed times \code{Y_*} with indicators
#' \code{delta_*} of a massive simulated dataset. As \code{N_super} grows, this should converge to the analytic truth.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{WR}: Monte Carlo estimate of the true win ratio
#'     \item \code{pi_T}, \code{pi_C}: Population win probabilities
#'     \item \code{treated_wins}, \code{control_wins}, \code{ties}, \code{total_pairs}: Respective counts during the pairwise comparison process
#'     \item \code{N_super}: Superpopulation size used
#'   }
#'
#' @export
true_wr_monte_carlo <- function(N_super = 1e6,
                                seed = NULL,
                                ...) {

  # Set the specified seed
  if (!is.null(seed)) set.seed(seed)

  # Simulated the superpopulation of specified size, given the arguments to `simulate_dataset`
  sim_args <- list(...)
  sim_args$N <- N_super
  superpop <- do.call(simulate_dataset, sim_args)

  # Split the data into treatment and control groups
  treated <- superpop[superpop$A == 1, c("Y_D","delta_D","Y_H","delta_H"), drop = FALSE]
  control <- superpop[superpop$A == 0, c("Y_D","delta_D","Y_H","delta_H"), drop = FALSE]

  # Compute and return population WR
  result <- get_true_WR(treated, control)
  result$N_super <- N_super
  return(result)
}
