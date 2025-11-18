# Simulation Study

This study builds on theoretical results, and simulation studies, found in the following two papers:
- Causal Inference on Win Ratio for Observational Data With Dependent Subjects (Zhang et al., 2022)
- Inference on win ratio for cluster-randomized semi-competing risk data (Zhang 2021)
  
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

For $N$ independent subjects, we simulate and a hierarchical composite outcome with two time-to-event components, one fatal $T_{D}$ and one non-fatal $T_{H}$, a treatment indicator variable $A$, a measured covariate $X$, and an unmeasured covariate $U$. We assume the covariates are normally distributed:
$$
\quad X \sim \text{Normal}(\mu_{X}, \sigma_{X}^2), \quad U \sim \text{Normal}(\mu_{U},\sigma_{U}^2)
$$
*Note: Naturally, this analysis can be extended to multiple measured and unmeasured covariates by using multivariate normal distributions and replacing any relevant multiplication below with matrix multiplication.*

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
\Lambda_{D}=\lambda_{D}\exp \left\{\beta_{A,D} A+ \beta_{X,D} X + \beta_{U,D} U \right\}   ,
\quad \Lambda_{H}=\lambda_{H}\exp \left\{\beta_{A,H} A+ \beta_{X,H} X + \beta_{U,H} U \right\}
$$
where $\lambda_{D},\lambda_{H}$ are baseline hazard rates. The marginal survival functions are then given by:
$$
S_{D}(t_{D}|A,X,U)=\exp \left\{-(\Lambda_{D}t_{D})^{\kappa_{D}}  \right\}, 
\quad S_{H}(t_{H}|A,X,U)=\exp \left\{-(\Lambda_{H}t_{H})^{\kappa_{H}}  \right\} 
$$
We model the joint survival function for $(T_{D},T_{H})$ with the Gumbel-Hougaard copula as used in previous win ratio literature (Luo et el., 2015; Oakes, 2016):
$$
S_{D,H}(t_{D},t_{H}|A,X,U)=\exp \left\{ -\left[ (\Lambda_{D}t_{D})^{\kappa_{D}\theta} + (\Lambda_{H}t_{H})^{\kappa_{H}\theta} \right]^{1/\theta}  \right\} 
$$
where $\theta\geq1$ controls the association between $T_{D}$ and $T_{H}$ and $\theta=1$ denotes no association. Kendall's tau would then be $\tau=1-\theta ^{-1}$ (Oakes, 1989).

In our simulation, we don't generate the latent event times $T_{D}$ and $T_{H}$ directly. Instead, we draw the bivariate vector $(V_{D},V_{H})$ from the Gumbel-Hougaard copula (where $V_{D},V_{H}$ are marginally standard uniform random variables), and use the following transformation to derive correlated event times:
$$
T_{D}=\Lambda_{D}^{-1}\left( -\log(V_{D})\right)^{1/\kappa_{D}},\quad T_{H}=\Lambda_{H}^{-1}\left( -\log(V_{H})\right)^{1/\kappa_{H}}
$$
Lastly, we introduce a time-to-random-censoring variable $C$ and assume it is independent of $(T_{D},T_{H},X,U)$ conditional on $A$, such that:
$$
C \sim\text{Exponential}(\Lambda_{C})
$$
with scale parameter given by the proportional hazards model:
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
