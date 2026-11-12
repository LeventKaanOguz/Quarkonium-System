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


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
