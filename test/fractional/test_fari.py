import unittest

import numpy as np

from cami_amber.fractional_metrics import (
    compute_fari,
    compute_fari_from_stats,
    compute_fari_nxn_reference,
    fari_stats_from_memberships,
)


def _onehot(labels, n_classes=None):
    labels = np.asarray(labels)
    k = int(labels.max()) + 1 if n_classes is None else n_classes
    out = np.zeros((len(labels), k))
    out[np.arange(len(labels)), labels] = 1.0
    return out


class TestReferenceFARI(unittest.TestCase):
    """Oracle: line-for-line n×n translation of its-likeli-jeff/FARI R/fari.R
    (commit 9f2e7e8769120a26d8456242ba5ed77bf539f2c2). R was not available;
    the n×n Python translation is the independent oracle for small n.
    Production scoring uses the same equations on compact cross-products.
    """

    def test_stats_match_nxn_soft_vs_soft(self):
        rng = np.random.default_rng(0)
        U = rng.random((12, 3))
        U = U / U.sum(axis=1, keepdims=True)
        V = rng.random((12, 4))
        V = V / V.sum(axis=1, keepdims=True)
        self.assertAlmostEqual(compute_fari(U, V), compute_fari_nxn_reference(U, V), places=12)

    def test_stats_match_nxn_soft_vs_hard(self):
        rng = np.random.default_rng(1)
        U = rng.random((10, 3))
        U = U / U.sum(axis=1, keepdims=True)
        V = _onehot([0, 0, 1, 1, 1, 2, 2, 0, 1, 2])
        self.assertAlmostEqual(compute_fari(U, V), compute_fari_nxn_reference(U, V), places=12)

    def test_stats_match_nxn_hard_vs_hard(self):
        U = _onehot([0, 0, 1, 1, 1, 2])
        V = _onehot([1, 1, 0, 0, 0, 2])
        self.assertAlmostEqual(compute_fari(U, V), compute_fari_nxn_reference(U, V), places=12)

    def test_reflexive_soft(self):
        rng = np.random.default_rng(2)
        U = rng.random((8, 3))
        U = U / U.sum(axis=1, keepdims=True)
        self.assertAlmostEqual(compute_fari(U, U), 1.0, places=12)

    def test_label_permutation(self):
        rng = np.random.default_rng(3)
        U = rng.random((9, 3))
        U = U / U.sum(axis=1, keepdims=True)
        V = _onehot([0, 1, 1, 2, 2, 0, 0, 1, 2])
        Vperm = V[:, [2, 0, 1]]
        self.assertAlmostEqual(compute_fari(U, V), compute_fari(U, Vperm), places=12)

    def test_one_cluster_degenerate(self):
        U = np.ones((6, 1))
        V = np.ones((6, 1))
        self.assertAlmostEqual(compute_fari(U, V), 1.0, places=12)

    def test_weighted_equals_row_replication(self):
        U = np.array([[0.7, 0.3], [1.0, 0.0], [0.2, 0.8]])
        V = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        weights = np.array([2, 3, 1], dtype=float)
        closed = compute_fari(U, V, weights=weights)
        Urep = np.repeat(U, weights.astype(int), axis=0)
        Vrep = np.repeat(V, weights.astype(int), axis=0)
        self.assertAlmostEqual(closed, compute_fari_nxn_reference(Urep, Vrep), places=12)

    def test_wrong_normalization_is_detected(self):
        U = _onehot([0, 0, 1, 1, 1, 2])
        V = _onehot([1, 1, 0, 0, 0, 2])
        n, sA, sB, qA, qB, qAB, tA, tB = fari_stats_from_memberships(U, V)
        correct = compute_fari_from_stats(n, sA, sB, qA, qB, qAB, tA, tB)
        wrong = compute_fari_from_stats(n, sA, sB, qA, qB, qAB / 2.0, tA, tB)
        self.assertGreater(abs(correct - wrong), 1e-6)

    def test_fixed_oracle_hard(self):
        U = _onehot([0, 0, 1, 1, 1, 2])
        V = _onehot([1, 1, 0, 0, 0, 2])
        value = compute_fari_nxn_reference(U, V)
        self.assertAlmostEqual(value, compute_fari(U, V), places=12)
        self.assertTrue(np.isfinite(value))


if __name__ == '__main__':
    unittest.main()
