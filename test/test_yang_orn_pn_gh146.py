"""Tests of analyses/yang_orn_pn_gh146.py's pure functions on small made-up matrices (no model run, no shared directory)."""
import numpy as np
import pandas as pd
import pytest

from analyses import yang_orn_pn_gh146 as Y


def made_up(n_glom=10, n_odor=12, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.gamma(1.5, 20.0, size=(n_glom, n_odor)) - 5.0,
                        index=[f'g{i}' for i in range(n_glom)], columns=Y.ODORS[:n_odor])


def test_inhsc_for_keeps_the_default():
    assert Y.inhsc_for(Y.PN_INHADD_DEFAULT) == pytest.approx(Y.PN_INHSC_DEFAULT)
    assert Y.inhsc_for(1000.0) > Y.PN_INHSC_DEFAULT


def test_kennedy_ss_is_zero_without_odor_and_monotone():
    orn = made_up()
    zero = orn * 0.0
    assert np.allclose(Y.kennedy_ss(zero).to_numpy(), 0.0)
    out = Y.kennedy_ss(orn)
    # within an odor the step is monotone in the glomerulus's own input (one shared denominator)
    for o in orn.columns:
        order = np.argsort(orn[o].to_numpy())
        assert np.all(np.diff(out[o].to_numpy()[order]) >= -1e-12)
    # bounded by the tanh's 200 Hz
    assert out.to_numpy().max() < 200.0


def test_kennedy_ss_divides_by_the_odor_total():
    orn = made_up()
    big = orn.copy()
    big.iloc[1:] = big.iloc[1:] + 100.0        # every other glomerulus much stronger
    a = Y.kennedy_ss(orn).iloc[0]
    b = Y.kennedy_ss(big).iloc[0]
    assert np.all(b.to_numpy() <= a.to_numpy() + 1e-12)


def test_olsen_rectifies_and_normalizes():
    orn = made_up()
    out = Y.olsen(orn, sigma=12.0, m=0.05)
    assert np.all(out.to_numpy()[orn.to_numpy() <= 0] == 0.0)
    assert np.all((out.to_numpy() >= 0) & (out.to_numpy() < 1))
    stronger = Y.olsen(orn, sigma=12.0, m=0.5)
    assert np.all(stronger.to_numpy() <= out.to_numpy() + 1e-12)


def test_signed_power_keeps_sign_and_power_one_is_identity():
    orn = made_up()
    assert np.allclose(Y.signed_power(orn, 1.0).to_numpy(), orn.to_numpy())
    out = Y.signed_power(orn, 0.5).to_numpy()
    assert np.all(np.sign(out) == np.sign(orn.to_numpy()))
    assert np.allclose(np.abs(out), np.sqrt(np.abs(orn.to_numpy())))


def test_fit_power_recovers_a_power():
    orn = made_up(seed=6).abs() + 1.0
    g = Y.signed_power(orn, 1.7)
    p = Y.fit_power(orn, g, list(orn.columns))
    assert p['power'] == pytest.approx(1.7, abs=1e-2) and p['r_train'] == pytest.approx(1.0)


def test_olfsysm_point_gain_factor_one_is_the_ladder():
    a, s, o = Y.olfsysm_point([np.log(1000.0), 0.0, 1.5])
    assert a == pytest.approx(1000.0) and s == pytest.approx(Y.inhsc_for(1000.0)) and o == 1.5
    a, s, o = Y.olfsysm_point([np.log(Y.PN_INHADD_DEFAULT), np.log(2.0), Y.PN_OFFSET_DEFAULT])
    assert s == pytest.approx(2.0 * Y.PN_INHSC_DEFAULT)


def test_metrics_identity_and_scale_invariance():
    g = made_up(seed=1).to_numpy()
    m = Y.metrics(g, g)
    for k in Y.R_METRICS:
        assert m[k] == pytest.approx(1.0)
    m2 = Y.metrics(3.0 * g, g)
    for k in ('r_all', 'rho_all', 'r_pattern', 'r_tuning', 'r_totals', 'r_similarity'):
        assert m2[k] == pytest.approx(1.0)
    assert m2['cv_totals'] == pytest.approx(m['cv_totals'])


def test_metrics_drop_glomeruli_missing_in_the_data():
    g = made_up(seed=2).to_numpy()
    x = g + np.random.default_rng(3).normal(0, 5, g.shape)
    g2 = g.copy()
    g2[0, 3] = np.nan
    m = Y.metrics(x, g2)
    assert m['n_glomeruli'] == g.shape[0] - 1
    assert m['r_all'] == pytest.approx(Y.metrics(x[1:], g[1:])['r_all'])


def test_ceiling_of_identical_flies_is_one_and_of_noise_below():
    g = made_up(seed=4).to_numpy()
    same = np.stack([g] * 7)
    c = Y.ceiling(same).set_index('metric')
    assert np.allclose(c.loc[list(Y.R_METRICS), 'ceiling'], 1.0)
    rng = np.random.default_rng(5)
    noisy = np.stack([g + rng.normal(0, 15, g.shape) for _ in range(7)])
    c2 = Y.ceiling(noisy).set_index('metric')
    assert (c2.loc['r_all', 'ceiling'] < 1.0) and (c2.loc['r_all', 'n_splits'] == 35)


def test_loo_maps_held_out_columns_onto_the_data():
    g = made_up(seed=7)
    calls = []

    def predict(train):
        held = [o for o in g.columns if o not in train]
        assert len(held) == 1
        calls.append(held[0])
        # a different unit in every fold: the map must undo it
        return 10.0 ** len(calls) * g + 3.0, dict(n_train=len(train))

    m, rows = Y.loo(predict, 'x', g, verbose=False)
    assert list(m.columns) == list(g.columns) and len(rows) == len(g.columns) and calls == list(g.columns)
    assert all(r['n_train'] == len(g.columns) - 1 for r in rows)
    assert np.allclose(m.to_numpy(), g.to_numpy())


def test_loo_of_a_fixed_matrix_keeps_its_correlation_close():
    g = made_up(seed=8)
    x = g + np.random.default_rng(9).normal(0, 8, g.shape)
    m, _ = Y.loo(lambda train: (x, {}), 'fixed', g, verbose=False)
    plain = Y.metrics(x.to_numpy(), g.to_numpy())['r_all']
    assert abs(Y.metrics(m.to_numpy(), g.to_numpy())['r_all'] - plain) < 0.02


def test_centered_interval_brackets_the_value():
    lo, hi = Y.centered_interval(0.8, pd.Series([0.70, 0.72, 0.74, 0.76, 0.78]))
    assert lo < 0.8 < hi and (0.8 - lo) == pytest.approx(0.74 - np.percentile([0.70, 0.72, 0.74, 0.76, 0.78], 5))
