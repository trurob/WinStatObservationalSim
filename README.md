# WinStatObservationalSim

**Win Statistics in Observational Studies**

## Overview

`WinStatObservationalSim` is an R package for simulating observational datasets with semi-competing time-to-event outcomes under treatment assignment with potential unmeasured confounding. The package implements data generation mechanisms based on Weibull marginal distributions with Gumbel-Hougaard copula dependence, and provides tools to compute the true causal win ratio via both analytical integration and Monte Carlo approximation.

This package is designed to support simulation studies investigating the performance of win ratio estimators in observational settings, particularly focusing on:

- The impact of unmeasured confounding on bias
- Non-proportional hazards scenarios
- Effects of censoring mechanisms
- Hierarchical composite outcomes with semi-competing risks

The simulation framework builds on theoretical results from Zhang et al. (2021, 2022), Mao (2018, 2019, 2020), and Luo et al. (2015).

## Installation

You can install the development version of `WinStatObservationalSim` from GitHub:

```r
# Install from GitHub using remotes
install.packages("remotes")
remotes::install_github("trurob/WinStatObservationalSim")
```

Or using `devtools`:

```r
# Install from GitHub using devtools
install.packages("devtools")
devtools::install_github("trurob/WinStatObservationalSim")
```

## Main Functions

The package exports four primary functions:

### `simulate_dataset()`

Generates a single dataset of independent subjects under semi-competing risks with:
- Measured covariates $X$ with dimension $p$
- Unmeasured covariates $U$ with dimension $q$
- Treatment assignment $A$ (randomized or confounded)
- Observed semi-competing outcomes: $Y_D$, $\delta_D$, $Y_H$, $\delta_H$

Supports both Normal and Bernoulli distributions for covariates, Weibull marginal distributions for event times with Gumbel-Hougaard copula dependence, and flexible censoring mechanisms (random exponential censoring and administrative censoring).

### `true_wr_analytic()`

Computes the true causal win ratio using numerical integration over the covariate distribution. This function evaluates the integral expressions:

$$
\mathrm{WR}_{\text{causal}} = \frac{\mathbb{E}_{X,U}[A(X,U)+B(X,U)]}{\mathbb{E}_{X,U}[C(X,U)+D(X,U)]}
$$

where $A$, $B$, $C$, $D$ are components derived from the joint distribution of $(T_H, T_D, C)$ under treatment and control (see Mathematical Details below).

### `true_wr_monte_carlo()`

Computes the true causal win ratio via Monte Carlo approximation by generating a large superpopulation and calculating the empirical population-level win ratio. As the superpopulation size increases, this converges to the analytical truth.

### `get_true_WR()`

Given pre-split treatment and control dataframes, performs all-against-all pairwise comparisons to compute the population win ratio. Uses an efficient counting approach to avoid explicit enumeration of all pairs.

## Mathematical Details

### Data Model

For $N$ independent subjects, we simulate a hierarchical composite outcome with two time-to-event components, one fatal $T_{D}$ and one non-fatal $T_{H}$, a treatment indicator variable $A$, $p$ measured covariates $X$, and $q$ unmeasured covariates $U$. 

We assume that covariates can follow either a Normal or a Bernoulli distribution. Measured covariates are parameterized by:

$$
\quad X_{j} \sim \text{Normal}(\mu_{X}, \sigma_{X}^2)\quad\text{or}\quad X_{j}\sim\text{Bernoulli}(\rho_{X})\quad \text{for}\quad j \in \{ 1,\dots,p \}
$$

while unmeasured covariates are parameterized by:

$$
\quad U_{k} \sim \text{Normal}(\mu_{U}, \sigma_{U}^2)\quad\text{or}\quad U_{k}\sim\text{Bernoulli}(\rho_{U})\quad \text{for}\quad k \in \{ 1,\dots,q \}
$$

*Note: For simplicity of notation, the data generating mechanism below uses $p=1$ and $q=1$. Naturally, for multiple covariates, any relevant multiplication below can be replaced by matrix multiplication.*

For treatment assignment, we explore two settings: "randomized" treatment with $A \sim \text{Bernoulli}(0.5)$, and a "confounded" treatment assignment with the logistic model:

$$
\mathbb{P}_{}\left( A=1|X,U \right) = \frac{1}{1+\exp \left\{ -(\alpha_{0} + \alpha_{X}X + \alpha_{U}U) \right\} }
$$

We then assume Weibull marginal distributions for the latent time-to-event components (Cox & Oakes 1984):

$$
T_{D}\sim\text{Weibull}\left( \Lambda_{D}, \kappa_{D}\right), \quad T_{H}\sim\text{Weibull}\left( \Lambda_{H}, \kappa_{H}\right)
$$

shape parameters $\kappa_{D},\kappa_{H}\in(0,\infty)$ and the following proportional hazards model for the Weibull's scale parameters:

$$
\Lambda_{D}=\lambda_{D}\exp \left\{-\left( \beta_{A,D} A+ \beta_{X,D} X + \beta_{U,D} U \right)  \right\}   ,
\quad \Lambda_{H}=\lambda_{H}\exp \left\{-\left( \beta_{A,H} A+ \beta_{X,H} X + \beta_{U,H} U \right)  \right\}
$$

where $\lambda_{D},\lambda_{H}$ are baseline hazard rates. We follow Zhang (2021, 2022) and parameterize the hazard models so that positive coefficients correspond to reduction in event hazards. The marginal survival functions are then given by:

$$
S_{D}(t_{D}|A,X,U)=\exp \left\{-(\Lambda_{D}t_{D})^{\kappa_{D}}  \right\}, 
\quad S_{H}(t_{H}|A,X,U)=\exp \left\{-(\Lambda_{H}t_{H})^{\kappa_{H}}  \right\} 
$$

We model the joint survival function for $(T_{D},T_{H})$ with the Gumbel-Hougaard copula as used in previous win ratio literature (Luo et al., 2015; Oakes, 2016):

$$
S_{D,H}(t_{D},t_{H}|A,X,U)=\exp \left\{ -\left[ (\Lambda_{D}t_{D})^{\kappa_{D}\theta} + (\Lambda_{H}t_{H})^{\kappa_{H}\theta} \right]^{1/\theta}  \right\} 
$$

where $\theta\geq1$ controls the association between $T_{D}$ and $T_{H}$ and $\theta=1$ denotes no association. Kendall's tau would then be $\tau=1-\theta ^{-1}$ (Oakes, 1989).

In our simulation, we don't generate the latent event times $T_{D}$ and $T_{H}$ directly. Instead, we draw the bivariate vector $(V_{D},V_{H})$ from the Gumbel-Hougaard copula (where $V_{D},V_{H}$ are marginally standard uniform random variables), and use the following transformation to derive correlated event times:

$$
T_{D}=  \frac{\left( -\log(V_{D})\right)^{1/\kappa_{D}}}{\Lambda_{D}}  ,\quad T_{H}=\frac{\left( -\log(V_{H})\right)^{1/\kappa_{H}}}{\Lambda_{H}}
$$

Lastly, we introduce a time-to-random-censoring variable $C$ and assume it is independent of $(T_{D},T_{H},X,U)$ conditional on $A$, such that:

$$
C \sim\text{Exponential}(\Lambda_{C})
$$

with rate parameter given by the proportional hazards model:

$$
\Lambda_{C}=\lambda_{C}\exp \left\{ -\beta_{A,C}A \right\} 
$$

and the "censoring" survival function

$$
Q(t|A)=\exp(-\Lambda_{C}t)
$$

Now, let $\varphi$ denote the time until administrative censoring, which we will assume is a fixed known value shared by all participants. The observed time-to-fatal-event, and respective event indicator, are then given by:

$$
Y_{D}=\text{min}(T_{D},C,\varphi),\quad \delta_{D}=\mathbb{1}\{T_{D}\leq\text{min}(C,\varphi)\}
$$

and the observed time-to-non-fatal-event, and respective event indicator, are given by:

$$
Y_{H}=\text{min}(T_{H},T_{D},C,\varphi),\quad \delta_{H}=\mathbb{1}\{T_{H}\leq \text{min}(T_{D},C,\varphi)\}
$$

which characterizes the semi-competing risk paradigm. 

Thus, our observed dataset is:

$$
\mathcal{O}=\left\{ (A_{i},X_{i},Y_{D,i},\delta_{D,i},Y_{H,i},\delta_{H,i}):i=1\dots,N \right\} 
$$

### True Causal Win Ratio

Let $X \in \mathbb{R}^p$ and $U \in \mathbb{R}^q$ denote measured and unmeasured covariates with joint density $f_{X,U}(x,u)$.  
For treatment $a \in \{0,1\}$, define the latent event times $(T_{D,a}, T_{H,a})$ with respective survival functions:

$$
S_{D,a}(t \mid x,u)
  = \exp\!\left\{ -(\Lambda_{D,a}(x,u)\, t)^{\kappa_D} \right\}, \quad S_{H,a}(t \mid x,u)
  = \exp\!\left\{ -(\Lambda_{H,a}(x,u)\, t)^{\kappa_H} \right\},
$$

with respective proportional hazards models:

$$
\Lambda_{D,a}(x,u)
  = \lambda_D \exp\!\left[-(\beta_{A,D}\, a 
          + \beta_{X,D}^\top x + \beta_{U,D}^\top u)\right], \quad \Lambda_{H,a}(x,u)
  = \lambda_H \exp\!\left[-(\beta_{A,H}\, a 
          + \beta_{X,H}^\top x + \beta_{U,H}^\top u)\right].
$$

The joint survival function for $(T_{D,a}, T_{H,a})$ under treatment $a$ is given by the Gumbel–Hougaard copula:

$$
S_{D,H,a}(y_2, y_1 \mid x,u)
  = \exp\!\left\{
     -\left[
       (\Lambda_{D,a}(x,u)\, y_2)^{\kappa_D \theta}
       + (\Lambda_{H,a}(x,u)\, y_1)^{\kappa_H \theta}
     \right]^{1/\theta}
   \right\}
$$

Define the conditional survival of $T_H$ beyond $y_1$ given survival of $T_D$ beyond $y_2$ as

$$
G_a(y_1,y_2 \mid x,u)
  = \frac{S_{D,H,a}(y_2, y_1 \mid x,u)}{S_{D,a}(y_2 \mid x,u)}.
$$

Define the density of $T_D$ at $y_2$ under treatment $a$:

$$
\lambda_{D,a}(y_2 \mid x,u)
 = \frac{\partial}{\partial y_2}
   \left[ 1 - S_{D,a}(y_2 \mid x,u) \right].
$$

Define the hazard for $T_H$ at $y_1$ conditional on $T_D = y_2$ under treatment $a$:

$$
\lambda_{H \mid D, a}(y_1 \mid y_2, x,u)
  = \frac{
      \frac{\partial}{\partial y_1} 
      G_a(y_1,y_2 \mid x,u)
    }{
      G_a(y_1,y_2 \mid x,u)
    }.
$$

Define the censoring survival and hazards under treatment $a$:

$$
Q_a(t)
  = \exp\!\left( -\Lambda_{C,a} t \right),
\qquad
\Lambda_{C,a} = \lambda_C \exp(-\beta_{A,C} a),
$$

$$
\lambda_{C,a}(t) = \Lambda_{C,a}.
$$

Using the notation and derivation established in Zhang 2021, we define the four terms:

$$
\begin{align}
A(x,u)
  & = \int_0^{\infty}
     S_{D,1}(y_2\mid x,u)\,
     Q_1(y_2)\,
     S_{D,0}(y_2\mid x,u)\,
     Q_0(y_2)\,
     \lambda_{D,0}(y_2\mid x,u)\,
   dy_2 \\
 \\
B(x,u)
  & = \int_0^{\infty} \int_0^{y_2}
     G_1(y_1,y_2\mid x,u)\,
     Q_1(y_2)\,
     G_0(y_1,y_2\mid x,u)\,
     Q_0(y_2)\,
     \lambda_{H\mid D,0}(y_1\mid y_2,x,u)\,
     \{\lambda_{C,0}(y_2)+\lambda_{C,1}(y_2)\}
   \, dy_1\, dy_2 \\
 \\
C(x,u)
  & = \int_0^{\infty}
     S_{D,1}(y_2\mid x,u)\,
     Q_1(y_2)\,
     S_{D,0}(y_2\mid x,u)\,
     Q_0(y_2)\,
     \lambda_{D,1}(y_2\mid x,u)\,
   dy_2 \\
 \\
D(x,u)
  & = \int_0^{\infty} \int_0^{y_2}
     G_1(y_1,y_2\mid x,u)\,
     Q_1(y_2)\,
     G_0(y_1,y_2\mid x,u)\,
     Q_0(y_2)\,
     \lambda_{H\mid D,1}(y_1\mid y_2,x,u)\,
     \{\lambda_{C,0}(y_2)+\lambda_{C,1}(y_2)\}
   \, dy_1\, dy_2.
\end{align}
$$

Then the true causal win ratio is defined as:

$$
\mathrm{WR}_{\text{causal}}
  = \frac{
      \mathbb{E}_{X,U}[A(X,U)+B(X,U)]
    }{
      \mathbb{E}_{X,U}[C(X,U)+D(X,U)]
    }
      = \frac{
      \displaystyle
      \int (A(x,u)+B(x,u))\, f_{X,U}(x,u)\, dx\,du
    }{
      \displaystyle
      \int (C(x,u)+D(x,u))\, f_{X,U}(x,u)\, dx\,du
    }
$$

As mentioned by Zhang, when there are a small number of covariates this expression can feasibly be calculated analytically using numerical integration. When the number of covariates is large, the integrals become intractable, and a Monte Carlo approach becomes more feasible.

## References

This package builds on theoretical results and simulation studies from:

- Zhang, D., Tao, Y., Tong, G., Liu, L., Shepherd, B. E., & Shu, D. (2022). Causal inference on win ratio for observational data with dependent subjects. *Statistics in Medicine*, 41(27), 5403–5418.

- Zhang, D. (2021). Inference on win ratio for cluster-randomized semi-competing risk data. *Biometrics*, 77(4), 1433–1445.

- Mao, L. (2020). Sample size formula for general win ratio analysis. *Biometrics*, 76(4), 1386–1396.

- Mao, L. (2019). On the alternative hypotheses for the win ratio. *Biometrics*, 75(1), 347–351.

- Mao, L. (2018). On causal estimation using U-statistics. *Biostatistics*, 19(3), 294–308.

- Luo, X., Tian, H., Mohanty, S., & Tsai, W. Y. (2015). An alternative approach to confidence interval estimation for the win ratio statistic. *Biometrics*, 71(1), 139–145.

- Cox, D. R., & Oakes, D. (1984). *Analysis of Survival Data*. Chapman and Hall.

- Oakes, D. (1989). Bivariate survival models induced by frailties. *Journal of the American Statistical Association*, 84(406), 487–493.

- Oakes, D. (2016). On the win-ratio statistic in clinical trials with multiple types of event. *Biometrika*, 103(3), 742–745.

## License

GPL (>= 3)

## Issues

For bug reports and feature requests, please visit:
https://github.com/trurob/WinStatObservationalSim/issues
