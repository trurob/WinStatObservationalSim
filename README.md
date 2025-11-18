# Simulation Study

This study builds on theoretical results, and simulation studies, found in the following two papers:
- Causal Inference on Win Ratio for Observational Data With Dependent Subjects (Zhang et al., 2022)
- Inference on win ratio for cluster-randomized semi-competing risk data (Zhang et al., 2021)
  
Which, in turn, build on simulation studies and theoretical results found in:
- Sample size formula for general win ratio analysis (Mao 2020)
- On the Alternative Hypotheses for the Win Ratio (Mao 2019)
- On causal estimation using U-statistics (Mao 2018)
- An Alternative Approach to Confidence Interval Estimation for the  Win Ratio Statistic (Luo et al., 2015)
  
We build on previous simulation studies by:
- Modeling marginal time-to-events using the Weibull distribution, which will allow us to study more complex event distributions and non-proportional hazards.
- Studying the effect of unmeasured confounding on bias.
- Using a real observational dataset to calibrate simulation parameters.
  
## Data model

For $N$ independent subjects, we simulate a hierarchical composite outcome with two time-to-event components, one fatal $T_{D}$ and one non-fatal $T_{H}$, a treatment indicator variable $A$, $p$ measured covariates $X$, and $q$ unmeasured covariates $U$.

We assume that covariates can follow either a Normal or a Bernoulli distribution. Measured covariates are parameterized by:

$$
X_{j} \sim \text{Normal}(\mu_{X}, \sigma_{X}^2)\quad\text{or}\quad X_{j}\sim\text{Bernoulli}(\rho_{X})\quad \text{for}\quad j \in \{1,\dots,p\}
$$

while unmeasured covariates are parameterized by:

$$
U_{k} \sim \text{Normal}(\mu_{U}, \sigma_{U}^2)\quad\text{or}\quad U_{k}\sim\text{Bernoulli}(\rho_{U})\quad \text{for}\quad k \in \{1,\dots,q\}
$$

*Note: For simplicity of notation, the data generating mechanism below uses $p=1$ and $q=1$. Naturally, for multiple covariates, any relevant multiplication below can be replaced by matrix multiplication.*

For treatment assignment, we explore two settings: randomized treatment with $A \sim \text{Bernoulli}(0.5)$, and a confounded treatment assignment with the logistic model:

$$
\mathbb{P}(A=1\mid X,U) = \frac{1}{1+\exp\{-(\alpha_{0} + \alpha_{X}X + \alpha_{U}U)\}}
$$

We then assume Weibull marginal distributions for the latent time-to-event components (Cox & Oakes 1984):

$$
T_{D}\sim\text{Weibull}(\Lambda_{D}, \kappa_{D}), \qquad T_{H}\sim\text{Weibull}(\Lambda_{H}, \kappa_{H})
$$

with shape parameters $\kappa_{D},\kappa_{H}\in(0,\infty)$ and proportional hazards scale parameters:

$$
\Lambda_{D}=\lambda_{D}\exp\{-( \beta_{A,D} A+ \beta_{X,D} X + \beta_{U,D} U )\},
$$

$$
\Lambda_{H}=\lambda_{H}\exp\{-( \beta_{A,H} A+ \beta_{X,H} X + \beta_{U,H} U )\}
$$

where $\lambda_{D},\lambda_{H}$ are baseline hazard rates. We follow Zhang (2021, 2022) and parameterize the hazard models so that positive coefficients correspond to event-hazard reduction. The marginal survival functions are:

$$
S_{D}(t_{D}\mid A,X,U)=\exp\{-(\Lambda_{D}t_{D})^{\kappa_{D}}\},
\qquad
S_{H}(t_{H}\mid A,X,U)=\exp\{-(\Lambda_{H}t_{H})^{\kappa_{H}}\}
$$

We model the joint survival function for $(T_{D},T_{H})$ with the Gumbel–Hougaard copula:

$$
S_{D,H}(t_{D},t_{H}\mid A,X,U)
= \exp\!\left(
 -\left[
   (\Lambda_{D}t_{D})^{\kappa_{D}\theta} +
   (\Lambda_{H}t_{H})^{\kappa_{H}\theta}
 \right]^{1/\theta}
\right)
$$

where $\theta\geq1$ controls association, and Kendall’s $\tau = 1 - 1/\theta$ (Oakes, 1989).

In our simulation, we do not generate latent event times $T_{D}$ and $T_{H}$ directly. Instead, we draw $(V_{D},V_{H})$ from the Gumbel–Hougaard copula (with uniform marginals), and transform:

$$
T_{D}= \frac{(-\log V_{D})^{1/\kappa_{D}}}{\Lambda_{D}},
\qquad
T_{H}= \frac{(-\log V_{H})^{1/\kappa_{H}}}{\Lambda_{H}}
$$

We introduce a time-to-random-censoring variable $C$ independently of $(T_{D},T_{H},X,U)$ given $A$:

$$
C \sim \text{Exponential}(\Lambda_{C})
$$

with rate:

$$
\Lambda_{C}=\lambda_{C}\exp\{-\beta_{A,C}A\}
$$

and censoring survival:

$$
Q(t\mid A)=\exp(-\Lambda_{C}t)
$$

Let $\varphi$ denote administrative censoring. The observed fatal-event time and indicator are:

$$
Y_{D}=\min(T_{D},C,\varphi),
\qquad
\delta_{D}=\mathbb{1}\{T_{D}\leq\min(C,\varphi)\}
$$

The observed non-fatal event time and indicator are:

$$
Y_{H}=\min(T_{H},T_{D},C,\varphi),
\qquad
\delta_{H}=\mathbb{1}\{T_{H}\leq \min(T_{D},C,\varphi)\}
$$

