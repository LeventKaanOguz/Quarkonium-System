"""Model-selection and model-averaging statistics (no ROOT; numpy + scipy).

Every function takes plain numbers, so analyze.py stays a thin driver:

- information criteria: AIC, small-sample AICc, BIC, and a Laplace-approximated
  Bayesian evidence with stated flat priors;
- predictive scores from the exact leave-one-out refits: elpd_loo, pseudo-BMA
  weights and stacking weights;
- Burnham-Anderson model averaging with the unconditional (model + parameter)
  error;
- multiple-testing corrections (Holm, Benjamini-Hochberg);
- Kendall rank correlation between rankings.
"""
import math

import numpy as np
from scipy import stats

LN10 = math.log(10.0)


# ------------------------------------------------------ information criteria --
def aic(chi2, p):
    return chi2 + 2 * p


def aicc(chi2, p, n):
    """AIC with the small-sample correction 2p(p+1)/(n-p-1); inf if n <= p+1."""
    return chi2 + 2 * p + (2 * p * (p + 1) / (n - p - 1) if n > p + 1 else math.inf)


def bic(chi2, p, n):
    return chi2 + p * math.log(n)


def cov_u(names, values, sigma, corr):
    """Covariance of the fit variables u (u = ln theta, except c which is
    linear) from the physical errors and correlations."""
    v = np.asarray(values, float)
    s = np.asarray(sigma, float)
    su = np.array([s[k] if n == "c" else s[k] / v[k] for k, n in enumerate(names)])
    return np.outer(su, su) * np.asarray(corr, float)


def laplace_log_evidence(chi2, names, values, sigma, corr, decades=1.0, c_range=1.0):
    """ln Z (up to a model-independent constant) for a Gaussian likelihood
    exp(-chi2/2) and flat priors: `decades` decades in every log-parameter,
    c_range GeV in c. ln Z = -chi2/2 + p/2 ln 2pi + 1/2 ln det C_u - sum ln W.
    NaN if the covariance is not positive definite (degenerate fit) or the
    posterior is not inside the prior (2 sigma_u > W: Laplace does not apply)."""
    C = cov_u(names, values, sigma, corr)
    sign, logdet = np.linalg.slogdet(C)
    if sign <= 0 or not np.isfinite(logdet):
        return math.nan
    p = len(names)
    W = np.array([c_range if n == "c" else decades * LN10 for n in names])
    # Laplace needs the posterior well inside the prior: +-1 sigma within W
    if np.any(2 * np.sqrt(np.diag(C)) > W):
        return math.nan
    lnW = float(np.sum(np.log(W)))
    return -0.5 * chi2 + 0.5 * p * math.log(2 * math.pi) + 0.5 * logdet - lnW


def ic_weights(values):
    """exp(-Delta/2) weights of a criterion on the -2 ln scale (AIC, AICc, BIC)."""
    x = np.asarray(values, float)
    ok = np.isfinite(x)
    w = np.zeros_like(x)
    w[ok] = np.exp(-0.5 * (x[ok] - x[ok].min()))
    return w / w.sum()


def evidence_weights(lnz):
    """Posterior model probabilities from ln Z with equal prior odds."""
    x = np.asarray(lnz, float)
    ok = np.isfinite(x)
    w = np.zeros_like(x)
    w[ok] = np.exp(x[ok] - x[ok].max())
    return w / w.sum()


# ------------------------------------------------------- predictive scores --
def log_predictive_density(miss, sd):
    """Pointwise log N(miss; 0, sd^2) of held-out data. `sd` is the full
    predictive error from the refit without the data: sqrt((sigma_th w)^2 +
    s_exp^2 + sigma_par^2)."""
    m = np.asarray(miss, float)
    s2 = np.asarray(sd, float) ** 2
    return -0.5 * (np.log(2 * math.pi * s2) + m ** 2 / s2)


def t_coverage(k, nu):
    """P(|T| < k) for T ~ t_nu: the coverage of a +-k predictive interval when
    the scale is estimated with nu degrees of freedom (normal for nu = inf)."""
    return float(stats.norm.cdf(k) - stats.norm.cdf(-k)) if not math.isfinite(nu) \
        else float(stats.t.cdf(k, nu) - stats.t.cdf(-k, nu))


def variance_ratio_test(F, nu1, nu2):
    """Two-sided p-value of F = (chi2_1/nu1)/(chi2_2/nu2) ~ F(nu1, nu2)."""
    c = stats.f.cdf(F, nu1, nu2)
    return float(2 * min(c, 1 - c))


def pseudo_bma_weights(elpd):
    """w_k ~ exp(elpd_k) (Yao et al. 2018, without the Bayesian bootstrap)."""
    return evidence_weights(elpd)


def stacking_weights(logdens, iters=5000, tol=1e-12):
    """Stacking of predictive distributions: maximise sum_i log sum_k w_k p_ik
    over the simplex (EM updates). `logdens`: models x points."""
    L = np.asarray(logdens, float)
    ok = np.all(np.isfinite(L), axis=1)
    P = np.exp(L[ok] - L[ok].max())  # common shift: the argmax is unchanged
    K = P.shape[0]
    w = np.full(K, 1.0 / K)
    for _ in range(iters):
        mix = w @ P
        new = w * (P / mix).mean(axis=1)
        if np.max(np.abs(new - w)) < tol:
            w = new
            break
        w = new
    out = np.zeros(len(L))
    out[ok] = w
    return out


# ---------------------------------------------------------- model averaging --
def model_average(x, sigma, w):
    """Burnham-Anderson: mean = sum w x; unconditional error
    sum w sqrt(sigma^2 + (x - mean)^2); also the pure between-model spread."""
    x, s, w = (np.asarray(a, float) for a in (x, sigma, w))
    ok = np.isfinite(x) & np.isfinite(s) & (w > 0)
    x, s, w = x[ok], s[ok], w[ok] / w[ok].sum()
    mean = float(w @ x)
    uncond = float(w @ np.sqrt(s ** 2 + (x - mean) ** 2))
    spread = float(math.sqrt(w @ (x - mean) ** 2))
    return mean, uncond, spread


# ------------------------------------------------------- multiple testing --
def holm(p):
    """Holm step-down adjusted p-values (family-wise error)."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * p[i]))
        adj[i] = run
    return adj


def benjamini_hochberg(p):
    """Benjamini-Hochberg adjusted p-values (false discovery rate)."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)[::-1]
    adj = np.empty(m)
    run = 1.0
    for k, i in enumerate(order):
        rank = m - k
        run = min(run, p[i] * m / rank)
        adj[i] = run
    return adj


# ------------------------------------------------- residual structure tests --
def runs_test(signs):
    """Exact two-sided Wald-Wolfowitz runs test on a sequence of signs (+1/-1,
    zeros dropped): p-value of a run count at least as far from its mean as the
    observed one, over all arrangements with the same numbers of + and -.
    Few runs = residuals that drift together (underfitting)."""
    from itertools import combinations
    x = [1 if v > 0 else 0 for v in signs if v != 0]
    n, n1 = len(x), sum(x)
    if n1 == 0 or n1 == n:
        return math.nan, 1
    runs = lambda seq: 1 + sum(a != b for a, b in zip(seq, seq[1:]))
    r_obs = runs(x)
    mean = 1 + 2 * n1 * (n - n1) / n
    total = hit = 0
    for pos in combinations(range(n), n1):
        seq = [0] * n
        for i in pos:
            seq[i] = 1
        total += 1
        hit += abs(runs(seq) - mean) >= abs(r_obs - mean) - 1e-12
    return hit / total, r_obs


def normality_p(x):
    """Shapiro-Wilk p-value (NaN for fewer than 3 finite values)."""
    x = np.asarray([v for v in x if v is not None and math.isfinite(v)], float)
    return float(stats.shapiro(x).pvalue) if len(x) >= 3 else math.nan


# ---------------------------------------------------------- rank agreement --
def kendall_tau(a, b):
    """Kendall's tau-b between two score vectors (NaNs dropped pairwise)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return math.nan
    return float(stats.kendalltau(a[ok], b[ok]).statistic)
