"""Checks of analysis/model_selection.py against hand-computed values.
Run: python3 analysis/test_model_selection.py (also registered with CTest)."""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model_selection as ms  # noqa: E402


def close(a, b, tol=1e-12):
    assert np.allclose(a, b, atol=tol, rtol=0), (a, b)


def test_information_criteria():
    close(ms.aic(3.0, 4), 11.0)
    close(ms.aicc(3.0, 4, 8), 11.0 + 40.0 / 3.0)       # 2*4*5/(8-4-1)
    assert ms.aicc(1.0, 4, 5) == math.inf               # n <= p + 1
    close(ms.bic(3.0, 4, 8), 3.0 + 4 * math.log(8))
    close(ms.ic_weights([10.0, 12.0, math.nan]),
          [1 / (1 + math.exp(-1)), math.exp(-1) / (1 + math.exp(-1)), 0.0])


def test_laplace_evidence():
    # one log-parameter theta = 2 +- 0.2 -> sigma_u = 0.1; one decade of prior
    lnz = ms.laplace_log_evidence(4.0, ["m1"], [2.0], [0.2], [[1.0]])
    close(lnz, -2.0 + 0.5 * math.log(2 * math.pi) + math.log(0.1) - math.log(math.log(10)))
    # c is linear: sigma_u = sigma
    lnz_c = ms.laplace_log_evidence(0.0, ["c"], [0.3], [0.05], [[1.0]], c_range=2.0)
    close(lnz_c, 0.5 * math.log(2 * math.pi) + math.log(0.05) - math.log(2.0))
    assert math.isnan(ms.laplace_log_evidence(0.0, ["a", "b"], [1, 1], [1, 1],
                                              [[1, 1], [1, 1]]))   # singular
    # posterior wider than the prior: Laplace refused
    assert math.isnan(ms.laplace_log_evidence(0.0, ["c"], [0.0], [0.8], [[1.0]]))


def test_predictive_weights():
    ld = ms.log_predictive_density([0.0, 0.003], [0.004, 0.004])
    close(ld, [-0.5 * math.log(2 * math.pi * 1.6e-5),
               -0.5 * (math.log(2 * math.pi * 1.6e-5) + 9e-6 / 1.6e-5)])
    close(ms.t_coverage(1.0, math.inf), 0.682689492137, tol=1e-9)
    close(ms.t_coverage(2.0, 3), 0.860674031441, tol=1e-9)
    close(ms.variance_ratio_test(1.0, 9, 9), 1.0)
    assert ms.variance_ratio_test(30.0, 4, 9) < 1e-3
    # one model dominates every point -> stacking puts all weight on it
    L = np.log(np.array([[0.9, 0.9, 0.9], [0.1, 0.1, 0.1]]))
    close(ms.stacking_weights(L), [1.0, 0.0], tol=1e-6)
    # two models each best on a different half -> equal stacking weights
    L = np.log(np.array([[0.9, 0.1], [0.1, 0.9]]))
    close(ms.stacking_weights(L), [0.5, 0.5], tol=1e-6)


def test_model_average():
    mean, unc, spread = ms.model_average([1.0, 3.0], [0.0, 0.0], [0.5, 0.5])
    close([mean, unc, spread], [2.0, 1.0, 1.0])
    mean, unc, _ = ms.model_average([1.0, 1.0], [0.3, 0.4], [0.5, 0.5])
    close([mean, unc], [1.0, 0.35])


def test_multiple_testing():
    close(ms.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
    close(ms.benjamini_hochberg([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])


def test_residual_tests():
    # + + + + - - - -: 2 runs, the most extreme; 2 of the 70 arrangements
    # have 2 runs and 2 have 8 (equally far from the mean 5): p = 4/70
    p, r = ms.runs_test([1, 1, 1, 1, -1, -1, -1, -1])
    assert r == 2
    close(p, 4 / 70)
    p, r = ms.runs_test([1, -1, 1, -1, 1, -1, 1, -1])
    assert r == 8
    close(p, 4 / 70)
    assert ms.runs_test([1, 1, 1])[0] != ms.runs_test([1, 1, 1])[0]  # NaN
    rng = np.random.default_rng(1)
    assert ms.normality_p(rng.normal(size=200)) > 0.01


def test_kendall():
    close(ms.kendall_tau([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
    close(ms.kendall_tau([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)


def test_decay_log_test():
    """The amended decay test: zero at Gamma_NLO = Gamma_PDG, symmetric in
    over- and underprediction, and unbounded for an overprediction (the
    linear pull with sigma_th = f Gamma_NLO is bounded below by -1/f)."""
    from selection import decay_log_chi2

    def system_with(nlo):
        w = {"label": "chi_c0(1P)", "process": "gamma gamma", "in_fit": 1,
             "pdg_keV": 2.0, "pdg_err_keV": 0.1, "nlo_keV": nlo,
             "nlo_band_lo_keV": 0.8 * nlo, "nlo_band_hi_keV": 1.2 * nlo,
             "nlo_err_par_keV": 0.05 * nlo, "v2": 0.3}
        return {"analysis": {"decay_benchmark": {"labels": [w["label"]]}}, "widths": [w]}

    close(decay_log_chi2(system_with(2.0))[0], 0.0)
    _, n, lab, z_over = decay_log_chi2(system_with(2.0 * 13))
    _, _, _, z_under = decay_log_chi2(system_with(2.0 / 13))
    assert n == 1 and lab == "chi_c0(1P)"
    close(z_over, -z_under, 1e-12)
    f = math.hypot(0.2, 0.3, 0.05)                   # relative error of the prediction
    assert z_over < -1 / f                             # beyond the linear bound
    close(z_over, -math.log(13) / math.hypot(f, 0.05), 1e-12)


def test_severity_and_best_compromise():
    # two-sided p -> standard deviations
    close(ms.exceedance_z(0.3173105078629141), 1.0, 1e-9)
    close(ms.exceedance_z(0.0455002638963584), 2.0, 1e-9)
    close(ms.exceedance_z(1.0), 0.0)
    assert ms.exceedance_z(0.0) == math.inf
    # a window: inside -> 0; 0.04 above the edge with sigma 0.01 -> 4
    close(ms.window_z(0.45, 0.40, 0.48, 0.01), 0.0)
    close(ms.window_z(0.52, 0.40, 0.48, 0.01), 4.0)
    close(ms.window_z(0.30, 0.40, 0.48, 0.05), 2.0)
    # the worst failure decides, not the count: one gross failure loses to two
    # marginal ones; ties go to fewer failures, then to the lower AIC
    gross = {"tag": "gross", "severity": 4.1, "n_fail": 1, "aic": 30.0}
    marginal = {"tag": "marginal", "severity": 2.6, "n_fail": 2, "aic": 50.0}
    tie_a = {"tag": "tie_a", "severity": 3.0, "n_fail": 2, "aic": 40.0}
    tie_b = {"tag": "tie_b", "severity": 3.0, "n_fail": 1, "aic": 60.0}
    tie_c = {"tag": "tie_c", "severity": 3.0, "n_fail": 1, "aic": 55.0}
    assert ms.best_compromise([gross, marginal])["tag"] == "marginal"
    assert ms.best_compromise([tie_a, tie_b, tie_c])["tag"] == "tie_c"
    assert ms.best_compromise([]) is None



def test_stale_ppp_guard():
    # a stats run made before a change of the width code is recognized by its widths
    from tables import widths_agree
    w = lambda x: {"widths": [{"label": "J/psi(1S)", "process": "e+ e-", "lo_keV": x}]}
    assert widths_agree(w(11.517), w(11.517 * (1 + 1e-7)))
    assert not widths_agree(w(14.727), w(13.0))        # the S-D sign changed this width
    assert not widths_agree(w(1.0), None)
    assert not widths_agree(w(1.0), {"widths": []})

if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
