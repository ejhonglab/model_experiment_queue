#!/usr/bin/env python3
"""Tests for E90 at the settled pre-odor level (experiments/E90_settled_preodor.py): the settled level is the mean over
-0.25 to 0 s of a made-up trace (the step at -0.25 s in, the odor onset out); E90's default pre-odor window and its
readings' windows are unchanged; the levels read from toy olfsysm variables are the window means of every stage, the
claws formed as E90 forms them, the same for every odor whose pre-odor period is the same; a saved stage dict at the
settled level carries the settled means and its window means unchanged; E82's builder keeps the lab's input by default.
None of these needs the shared directory or a simulation."""

import inspect
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E90_settled_preodor')
E90 = pytest.importorskip('experiments.E90_model_stages')

DT = 0.0005


def _times(t0: float = -0.5, t_end: float = 0.75, dt: float = DT) -> np.ndarray:
    n = int(round((t_end - t0) / dt))
    return t0 + np.arange(n) * dt


def test_settled_level_is_the_mean_over_the_last_quarter_second():
    t = _times()
    rng = np.random.default_rng(0)
    x = np.vstack([8.8 * (1 - np.exp(-(t + 0.5) / 0.06)), 3.0 + 0.01 * t, rng.standard_normal(len(t))])
    inside = (t >= -0.25 - 1e-9) & (t < -1e-9)
    np.testing.assert_allclose(X.settled_level(x, t), x[:, inside].mean(axis=1), rtol=0, atol=1e-12)
    idx = X.pre_window_steps(t, X.SETTLED_PRE_WINDOW, closed_left=True)
    assert np.isclose(t[idx[0]], -0.25) and t[idx[-1]] < 0 and np.isclose(t[idx[-1]], -DT)
    assert len(idx) == int(round(0.25 / DT))
    # a trace that settles by -0.3 s reads its settled value; E90's window reads the start-up too
    y = np.where(t < -0.3, 0.0, 5.0)
    assert X.settled_level(y, t) == 5.0
    assert y[E90._windows(t, E90.MODEL_WINDOW)['pre']].mean() < 5.0


def test_e90_default_window_is_unchanged():
    assert E90.PRE_WINDOW == (-0.5, 0.0)
    assert E90.MODEL_WINDOW == (0.0, 0.5)
    assert X.SETTLED_PRE_WINDOW == (-0.25, 0.0)
    for t0 in (-0.5, -2.0):
        t = _times(t0)
        pre = E90._windows(t, E90.MODEL_WINDOW)['pre']
        np.testing.assert_array_equal(pre, np.flatnonzero((t > -0.5 + 1e-9) & (t < -1e-9)))
        np.testing.assert_array_equal(pre, X.pre_window_steps(t, E90.PRE_WINDOW, closed_left=False))
        np.testing.assert_array_equal(X.window_steps(t)['e90'], pre)
        win = E90._windows(t, E90.MODEL_WINDOW)['win']
        np.testing.assert_array_equal(win, np.flatnonzero((t >= -1e-9) & (t <= 0.5 + 1e-9)))


def _toy_run(n_odors: int = 3, n_gl: int = 2, n_claws: int = 5, dt: float = 0.01):
    """Toy olfsysm variables: two glomeruli with one bouton each, claws on them, APL starting from zero at -0.5 s; the
    odors differ only after 0 s."""
    t = _times(-0.5, 0.75, dt)
    n = len(t)
    rng = np.random.default_rng(1)
    W = pd.DataFrame(np.eye(n_gl)[rng.integers(0, n_gl, n_claws)],
                     index=pd.MultiIndex.from_arrays([np.arange(n_claws), np.arange(n_claws)], names=['claw', H.KC_ID_LEVEL]),
                     columns=pd.MultiIndex.from_arrays([[f'g{i}' for i in range(n_gl)], np.arange(n_gl)],
                                                       names=['glomerulus', 'bouton']))
    apl = 8.0 * (1 - np.exp(-(t + 0.5) / 0.08))
    pn, btn, inh = [], [], []
    for o in range(n_odors):
        p = np.vstack([10.0 + 0.02 * np.sin(t), 6.0 + 0.0 * t]) + np.where(t >= 0, 2.0 + o, 0.0)
        a = apl + np.where(t >= 0, 1.0 + o, 0.0)
        a_prev = np.concatenate([[0.0], a[:-1]])
        b = np.maximum(p - 0.3 * a_prev, 0.0)
        b[:, 0] = np.nan
        pn.append(p)
        btn.append(b)
        inh.append(a[None, :])
    wk = rng.uniform(0.5, 1.5, n_claws)
    rv = SimpleNamespace(orn=SimpleNamespace(sims=[np.zeros((n_gl, n))]), pn=SimpleNamespace(pn_sims=pn, bouton_sims=btn),
                         kc=SimpleNamespace(inh_sims=inh, wAPLKC=wk))
    mp = SimpleNamespace(time_dt=dt, time_pre_start=-2.0, time_start=-0.5, time_end=0.75)
    return mp, rv, W, t, wk


def test_read_levels_on_toy_variables():
    mp, rv, W, t, wk = _toy_run()
    order = np.array([2, 0, 1])
    lv = X.read_levels(mp, rv, W, order, ['a', 'b', 'c'])
    assert lv['spread'] < 1e-12
    for v, idx in X.window_steps(t).items():
        i = int(order[0])
        b = rv.pn.bouton_sims[i][:, idx]
        a = rv.kc.inh_sims[i].ravel()
        np.testing.assert_allclose(lv[v]['pn'], rv.pn.pn_sims[i][:, idx].mean(axis=1), atol=1e-12)
        np.testing.assert_allclose(lv[v]['bouton'], b.mean(axis=1), atol=1e-12)
        claws = np.maximum(W.to_numpy() @ b - np.outer(wk, a[idx - 1]), 0.0).mean(axis=1)
        np.testing.assert_allclose(lv[v]['claw'], claws, atol=1e-12)
        assert np.isclose(lv[v]['apl'], a[idx].mean())
    assert lv['settled']['apl'] > lv['e90']['apl']
    np.testing.assert_allclose(lv['settled']['apl'], X.settled_level(rv.kc.inh_sims[0].ravel(), t))
    assert lv['n_steps'] == {'e90': 49, 'settled': 25}


def test_stage_at_settled_replaces_only_the_pre_odor_means(monkeypatch):
    odors = ['2h', 'eb']
    st = {f'{k}_{w}': pd.DataFrame(np.arange(6.0).reshape(3, 2) + (10.0 if w == 'win' else 1.0), columns=odors)
          for k in X.STAGES for w in ('win', 'pre')}
    st['apl'] = pd.DataFrame({'win': [9.0, 9.5], 'pre': [8.3, 8.3]}, index=odors)
    lv = {c: {'levels': {'e90': {k: np.ones(3) for k in X.STAGES} | {'apl': 8.3},
                         'settled': {k: np.array([0.5, 0.0, 2.0]) for k in X.STAGES} | {'apl': 8.8}}}
          for c in ('on', 'off')}
    monkeypatch.setattr(X, '_saved_stage', lambda part, tag, cc: st)
    monkeypatch.setattr(X, 'load_levels', lambda part, tag: lv)
    assert X.stage_at('B', 'real', 'on', 'e90') is st
    s = X.stage_at('B', 'real', 'acute_off', 'settled')
    for k in X.STAGES:
        assert s[f'{k}_win'] is st[f'{k}_win']
        np.testing.assert_array_equal(s[f'{k}_pre'].to_numpy(), np.repeat(np.array([0.5, 0.0, 2.0])[:, None], 2, axis=1))
        assert s[f'{k}_pre'].index.equals(st[f'{k}_pre'].index) and list(s[f'{k}_pre'].columns) == odors
    assert (s['apl']['pre'] == 8.8).all() and (st['apl']['pre'] == 8.3).all()
    u = E90.stage_units(s, 'claw')
    np.testing.assert_allclose(u.to_numpy(), st['claw_win'].to_numpy() - np.array([0.5, 0.0, 2.0])[:, None])


def test_e82_builder_keeps_the_lab_input_by_default():
    R = pytest.importorskip('experiments.E82_protocol_runs')
    p = inspect.signature(R.build_simulator).parameters
    assert p['orn'].default is None and p['kc_kc_tau'].default == R.KC_KC_TAU
