#!/usr/bin/env python3
"""Tests for E129 (the inhibition that outlasts the block reduced to a static level;
experiments/E129_block_survivor_static.py), on made-up data: the time-mean of the KC-KC input (reading 1, with
olfsysm's filter), the static arms' bases (reading 2), the split (reading 9), the realism of a level, the runner's
identities on a toy simulator with olfsysm's KC equations (the hemibrain matrix at scale 1 reproduces the saved APL-off
run exactly; a static inhibition of zeros changes nothing; the odor order of the static inhibition), the bisection of one
scale with its corner rule, the lanes and claims, and the classes, "differ" rule, gates and outcomes. No model
simulation."""

import types as _types

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E129_block_survivor_static')


def _toy_types(n_g=12, n_apbp=6, n_ab=14):
    return np.array(['g'] * n_g + ["a'b'"] * n_apbp + ['ab'] * n_ab, dtype=object)


def _toy_matrix(types, seed=0):
    rng = np.random.default_rng(seed)
    n = len(types)
    M = rng.poisson(1.0, size=(n, n)).astype(float) * (rng.random((n, n)) < 0.3)
    np.fill_diagonal(M, 0.0)
    return M * (n / M.sum())


# ---------------------------------------------------------------------------
# readings 1, 2, 9, 10

def test_time_mean_input():
    types = _toy_types()
    M = _toy_matrix(types)
    C = np.random.default_rng(1).poisson(2.0, size=(len(types), 5)).astype(float)
    T = X.time_mean_input(C, M, 26.4)
    assert np.allclose(T, 26.4 * M @ C / 0.5)
    assert not X.time_mean_input(np.zeros_like(C), M, 26.4).any()
    with pytest.raises(ValueError):
        X.time_mean_input(C, M[:-1, :-1], 1.0)
    with pytest.raises(ValueError):
        X.time_mean_input(C, M, -1.0)


def test_filter_integrates_to_the_spike_count():
    """Reading 1's premise: olfsysm's forward-Euler filter a(t) = (1 - dt/tau) a(t-1) + spikes(t-1) / tau integrates to
    one per spike (given a tail long enough for it to decay)."""
    dt, tau = 5e-4, 0.05
    rng = np.random.default_rng(0)
    n_steps = int(1.5 / dt)
    spikes = np.zeros(n_steps)
    spikes[rng.choice(int(0.5 / dt), size=7, replace=False)] = 1.0
    a, total = 0.0, 0.0
    for t in range(1, n_steps):
        a = (1 - dt / tau) * a + spikes[t - 1] / tau
        total += a * dt
    assert total == pytest.approx(spikes.sum(), rel=1e-6)


def test_static_bases():
    types = _toy_types()
    Tm = np.random.default_rng(2).gamma(2.0, 10.0, size=(len(types), 6))
    gm = Tm.mean()
    B = {a: X.static_base(a, Tm, types) for a in X.STATIC_ARMS}
    for a, b in B.items():
        assert b.shape == Tm.shape and np.isclose(b.mean(), gm), a
    assert np.array_equal(B['T'], Tm)
    assert np.allclose(B['G'], gm)
    assert np.allclose(B['K'], Tm.mean(axis=1, keepdims=True))
    for t in X.TYPES:
        m = types == t
        assert np.allclose(B['Y'][m], Tm[m].mean())
        # KP: K's values permuted within the type
        assert np.allclose(np.sort(B['KP'][m, 0]), np.sort(B['K'][m, 0]))
    assert not np.allclose(B['KP'], B['K'])
    assert np.array_equal(B['KP'], X.static_base('KP', Tm, types))          # seeded
    for a in ('K', 'Y', 'KP', 'G'):
        assert np.allclose(B[a], B[a][:, :1])                               # constant over the odors
    with pytest.raises(ValueError):
        X.static_base('P', Tm, types)
    with pytest.raises(ValueError):
        X.static_base('T', Tm, types[:-1])


def test_split_kckc():
    types = _toy_types()
    rng = np.random.default_rng(3)
    n, n_od = len(types), 7
    A = rng.gamma(2.0, 5.0, size=(n, n_od))
    s = X.split_kckc(A, types)
    assert np.isclose(s['share_type'] + s['share_kc_within_type'] + s['share_odor'] + s['share_rest'], 1.0)
    assert np.isclose(s['share_kc'], s['share_type'] + s['share_kc_within_type'])
    assert 0 <= s['share_type_x_odor'] <= s['share_rest'] + 1e-12
    # the same per KC / per odor / rest as E95's two-way split
    from experiments.E95_static_apl import two_way_split
    t = two_way_split(A)
    assert np.isclose(s['share_kc'], t['share_kc']) and np.isclose(s['share_odor'], t['share_odor'])
    assert np.isclose(s['share_rest'], t['share_rest'])
    # pure forms
    per_type = np.repeat(np.array([{'g': 1.0, "a'b'": 3.0, 'ab': 7.0}[x] for x in types])[:, None], n_od, axis=1)
    assert np.isclose(X.split_kckc(per_type, types)['share_type'], 1.0)
    per_kc = np.repeat(rng.random(n)[:, None], n_od, axis=1)
    assert np.isclose(X.split_kckc(per_kc, types)['share_kc'], 1.0)
    per_odor = np.repeat(rng.random(n_od)[None, :], n, axis=0)
    assert np.isclose(X.split_kckc(per_odor, types)['share_odor'], 1.0)
    tx = np.zeros((n, n_od))
    tx[types == 'g', 0] = 5.0
    tx[types == 'ab', 1] = 5.0
    sp = X.split_kckc(tx, types)
    assert sp['share_type_x_odor'] > 0.5
    assert np.isclose(sp['mean_g'], tx[types == 'g'].mean())


def test_level_realism():
    types = _toy_types()
    n = len(types)
    thr = np.full(n, 200.0)
    M = np.full((n, 4), 50.0)
    r = X.level_realism(M, thr, types, 100.0)
    assert r['level_mean'] == 50.0 and r['level_over_thr_max'] == 0.25 and r['level_over_fixed_thr_mean'] == 0.5
    assert r['level_above_thr_share'] == 0.0 and not r['level_above_thr']
    M[0, 0] = 300.0
    r = X.level_realism(M, thr, types, 100.0)
    assert r['level_above_thr'] and r['level_above_thr_share'] == pytest.approx(1 / (4 * n))
    assert not X.level_realism(None, thr, types, 100.0)['level_above_thr']


def test_uniform_neuron_and_type_rates_are_e122s():
    E122 = pytest.importorskip('experiments.E122_feedback_neuron')
    types = _toy_types()
    M = _toy_matrix(types)
    assert np.array_equal(X.uniform_neuron(M), E122.arm_array('U', M, types))
    idx = pd.MultiIndex.from_arrays([np.arange(len(types)), types], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    rng = np.random.default_rng(4)
    o = {}
    for c in ('on', 'off'):
        s = pd.DataFrame(rng.poisson(1.0, size=(len(types), 5)).astype(float), index=idx)
        o[c] = ((s >= 2).astype(float), s)
    mine, theirs = X.type_rates(o), E122.type_rates(o)
    for k, v in theirs.items():
        assert mine[k] == pytest.approx(v), k


def test_held_ok():
    assert X.held_ok(4000, 4000) and X.held_ok(4008, 4000) and X.held_ok(3992, 4000)
    assert not X.held_ok(4009, 4000) and not X.held_ok(3990, 4000)


# ---------------------------------------------------------------------------
# the runner on a toy simulator with olfsysm's KC equations

class _Toy:
    """A few KCs with olfsysm's somatic equations: dVm = (-Vm + drive - ext (odor on) - wKCKC @ a) dt / taum, a spike
    above the threshold (Vm reset to 0), a(t) = (1 - dt/tau) a(t-1) + spikes(t-1) / tau; odors in an internal order,
    the outputs' column j the internal odor order[j]; `rv.kc.ext_inh` (KCs x internal odors) and `mp.kc.kc_kc_tau` read
    at every call, as olfsysm does."""

    def __init__(self, types, n_od=4, seed=5, dt=5e-4):
        rng = np.random.default_rng(seed)
        self.types = types
        self.n = len(types)
        self.W = _toy_matrix(types, seed)
        self.drive = rng.gamma(3.0, 60.0, size=(self.n, n_od))
        self.spont = rng.uniform(5, 15, size=self.n)
        self.order = np.array([2, 0, 3, 1])[:n_od]
        self.dt, self.taum = dt, 0.01
        kc = _types.SimpleNamespace(N=self.n, ext_inh_scale=1.0, n_claws_active_to_spike=-1, kc_kc_tau=0.05, ves_p=0.0)
        self.mp = _types.SimpleNamespace(kc=kc, time_stim_start=0.0, time_stim_end=0.5)
        self.rv = _types.SimpleNamespace(kc=_types.SimpleNamespace(ext_inh=np.zeros((0, 0)), thr=None))
        self.idx = pd.MultiIndex.from_arrays([np.arange(100, 100 + self.n), types], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
        self.odors = [f'od{j}' for j in range(n_od)]
        self.n_calls = 0

    def simulate(self, params):
        self.n_calls += 1
        if float(params['wAPLKC_scale']) != 0.0:
            raise AssertionError('the toy is APL-off only')
        dt, tau = self.dt, float(self.mp.kc.kc_kc_tau)
        thr = float(params['fixed_thr']) + 2 * self.spont
        self.rv.kc.thr = thr[:, None]
        Wk = float(params['wKCKC_scale']) * self.W
        ext = np.asarray(self.rv.kc.ext_inh, dtype=float)
        n_od = self.drive.shape[1]
        V = np.zeros((self.n, n_od))
        a = np.zeros((self.n, n_od))
        sp_prev = np.zeros((self.n, n_od))
        counts = np.zeros((self.n, n_od))
        ts = np.arange(-0.1, 0.75, dt)
        for t in ts:
            on = 0.0 <= t < 0.5
            dV = -V + self.spont[:, None] + (self.drive if on else 0.0)
            if ext.size and on:
                dV = dV - ext
            a = (1 - dt / tau) * a + sp_prev / tau
            dV = dV - Wk @ a
            V = V + dV * dt / self.taum
            sp = V > thr[:, None]
            V[sp] = 0.0
            sp_prev = sp.astype(float)
            counts += sp_prev
        s = pd.DataFrame(counts[:, self.order], index=self.idx, columns=self.odors)
        return (s >= 2).astype(float), s


@pytest.fixture
def toy(monkeypatch):
    types = _toy_types()
    m = _Toy(types)
    p = dict(fixed_thr=150.0, wAPLKC_scale=0.0, wKCAPL_scale=1.0, wAPLPN_scale=0.0, wPNAPL_scale=1.0, wKCKC_scale=3.0)
    off = m.simulate(dict(p))
    saved = {'on': off, 'off': off}
    parts = dict(p=p, saved=saved, target_off=X.n_resp(off[0]), size=int(off[0].size), folder='toy', tag='toy',
                 orn=None)

    W0 = m.W.copy()

    def build(w, extra_kws=None, **kw):
        W = W0.copy() if not extra_kws else np.asarray(extra_kws['_wKCKC'], dtype=float)

        def sim(q):
            m.W = W
            return m.simulate(q)

        cells = {k: _types.SimpleNamespace(cell_contents=v) for k, v in
                 (('mp', m.mp), ('rv', m.rv), ('odor_order', m.order), ('wKCKC_unscaled', W.copy()))}
        return sim, dict(cells=cells)

    hb = pd.DataFrame(m.W.copy(), index=pd.Index(m.idx.get_level_values(H.KC_ID_LEVEL)),
                      columns=pd.Index(m.idx.get_level_values(H.KC_ID_LEVEL)))
    monkeypatch.setattr(X.E83, 'wiring_parts', lambda w: parts)
    monkeypatch.setattr(X.E83, 'build_simulator', build)
    monkeypatch.setattr(X, 'hemibrain_matrix', lambda idx: hb.copy())
    monkeypatch.setattr(X, '_preload', lambda: None)
    return m, parts


def test_runner_identities(toy):
    m, parts = toy
    R = X.OffRunner('real', 'HB')
    saved = parts['saved']['off']
    a = R.run(1.0, None)
    # the hemibrain matrix at scale 1: the saved APL-off run exactly
    assert a[1].equals(saved[1]) and a[0].equals(saved[0])
    # a static inhibition of zeros changes nothing
    b = R.run(1.0, np.zeros((m.n, len(m.odors))))
    assert b[1].equals(saved[1]) and b[0].equals(saved[0])
    # the static inhibition is removed after the run
    assert np.asarray(m.rv.kc.ext_inh).size == 0
    # T at factor 1 with KC-KC 0 is read (changes the counts) and its base is the time-mean of the saved counts
    assert np.allclose(R.Tm, X.time_mean_input(saved[1].to_numpy(), m.W, 3.0))
    t = R.run(0.0, R.Tm)
    assert (t[1].to_numpy() != saved[1].to_numpy()).any()
    assert R.thr is not None and R.thr.shape == (m.n,)


def test_runner_odor_order_of_the_static_inhibition(toy):
    """A static inhibition non-zero in one output column changes only that column (the outputs' order mapped to the
    simulator's)."""
    m, parts = toy
    R = X.OffRunner('real', 'HB')
    base = R.run(0.0, None)
    ext = np.zeros((m.n, len(m.odors)))
    ext[:, 1] = 80.0
    out = R.run(0.0, ext)
    d = (out[1].to_numpy() != base[1].to_numpy()).any(axis=0)
    assert d[1] and not d[[0, 2, 3]].any()


def test_runner_kinds(toy):
    m, _ = toy
    S = X.OffRunner('real', 'slow')
    assert m.mp.kc.kc_kc_tau == X.SLOW_TAU_S and S.tau == X.SLOW_TAU_S
    m.mp.kc.kc_kc_tau = 0.05
    with pytest.raises(RuntimeError):
        S.run(1.0)                                                          # tau changed behind its back
    U = X.OffRunner('real', 'U')
    U.run(1.0)
    assert np.allclose(np.asarray(m.W), X.uniform_neuron(_toy_matrix(m.types, 5)))
    assert U.tau == pytest.approx(0.05)
    with pytest.raises(ValueError):
        U.arm_fn('G')                                                       # a static arm on the U runner
    with pytest.raises(ValueError):
        X.OffRunner('real', 'nope')


def test_runner_checks(toy, monkeypatch):
    m, _ = toy
    m.mp.kc.ves_p = 0.1
    with pytest.raises(RuntimeError):
        X.OffRunner('real', 'HB')
    m.mp.kc.ves_p = 0.0
    m.mp.time_stim_end = 1.0
    with pytest.raises(RuntimeError):
        X.OffRunner('real', 'HB')


def test_held_off_bisects_one_scale(toy):
    m, parts = toy
    R = X.OffRunner('real', 'HB')
    # G: one level, bisected to the saved run's APL-off count
    res = X.held_off(R, 'G', 'toy G')
    assert res['factor'] > 0 and np.allclose(res['ext'], res['factor'] * R.Tm.mean())
    assert res['corner'] == (res['at_bound'] or not res['held'])
    assert res['params_off']['wKCKC_scale'] == 0.0 and res['params_off']['wAPLKC_scale'] == 0.0
    if res['held']:
        assert abs(res['n_off'] - parts['target_off']) <= X.RATE_TOL * parts['target_off']
    # Z: not bisected, no inhibition at all
    z = X.held_off(R, 'Z', 'toy Z')
    assert np.isnan(z['factor']) and z['ext'] is None and not z['corner'] and z['n_sims'] == 1
    assert z['n_off'] >= parts['target_off']
    # the factor's corner: a bracket that cannot reach the target
    import experiments.E129_block_survivor_static as mod
    old = mod.SCALE_BRACKET
    try:
        mod.SCALE_BRACKET = (0.999, 1.001)
        c = X.held_off(R, 'G', 'toy G corner')
        assert c['at_bound'] and c['corner']
    finally:
        mod.SCALE_BRACKET = old


# ---------------------------------------------------------------------------
# lanes and claims

def test_planned_lanes():
    s = X.planned('step1a', X.S_WIRINGS)
    t = X.planned('step1a', X.T_WIRINGS)
    assert len(s) == len(t) == 16 and not {X.tag_of(*p) for p in s} & {X.tag_of(*p) for p in t}
    assert s[:2] == [('G', 'real'), ('T', 'real')] and [a for a, _ in s[8:]] == ['P'] * 4 + ['S'] * 4
    assert {X.tag_of(*p) for p in s + t} == {f'S1_{a}_{w}' for a in 'GTPS' for w in X.WIRINGS}
    assert [a for a, _ in X.planned('step1b', ['real'], ('K', 'Y'))] == ['K', 'Y']
    assert X.planned('step2', ['P3']) == [('C50', 'P3'), ('C0', 'P3')]
    assert X.tag_of('Z', 'M2') == 'S1c_Z_M2' and X.tag_of('C0', 'real') == 'S2_C0_real'


def test_claims(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    assert X.claim('S1_G_real', 'lane1')
    assert not X.claim('S1_G_real', 'lane2')
    assert X.claim('S1_G_real', 'lane2', reclaim=True)
    assert X._read_json(tmp_path / 'claims' / 'S1_G_real' / 'claim.json')['lane'] == 'lane2'
    X.release('S1_G_real')
    assert X.claim('S1_G_real', 'lane3')


def test_step_gates_refuse(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    for step in ('step1a', 'step1b', 'step1c', 'step2'):
        ok, why = X.step_allowed(step)
        assert not ok and why
    with pytest.raises(SystemExit):
        X.part_runs('step1a', ['real'])
    with pytest.raises(SystemExit):
        X.part_step3(False)
    with pytest.raises(NotImplementedError):
        X.part_step3(True)
    assert X.step_arms('step1b', {}) == ('K', 'Y')
    assert X.step_arms('step1b', {'kp': {'open': True}}) == ('K', 'Y', 'KP')
    X._write_json(tmp_path / 'identity_h_a.json', {'wirings': {w: {'passes': True} for w in X.S_WIRINGS}})
    assert not X.identity_passed()
    X._write_json(tmp_path / 'identity_h_b.json', {'wirings': {w: {'passes': True} for w in X.T_WIRINGS}})
    assert X.identity_passed()


# ---------------------------------------------------------------------------
# classes, "differ", gates and outcomes

def _rd_runs(values, corner_wirings=None, seed=0):
    """A made-up readings table: per kind, a pattern part per wiring and per population for every target (calibrated),
    and the raw columns E123's classify reads; a runs table with every run in the bands."""
    rng = np.random.default_rng(seed)
    rows, runs = [], []
    for kind, base in values.items():
        for unit in list(X.WIRINGS) + list(X.POPS) + ['all8']:
            r = dict(kind=kind, unit=unit, readout='cal')
            for t in X.targets():
                r[f'{t}_pattern'] = base + rng.normal(0, 1e-4)
            rows.append(r)
            rows.append(dict(kind=kind, unit=unit, readout='raw', mean_effect=0.1, all_loss=0.1))
        for w in X.WIRINGS:
            runs.append(dict(kind=kind, wiring=w, in_bands=True, rate_rho_on=0.5,
                             corner=bool(corner_wirings and kind in corner_wirings and w in corner_wirings[kind]),
                             at_bound_thr=bool(corner_wirings and kind in corner_wirings and w in corner_wirings[kind]),
                             at_bound_apl=False))
    return pd.DataFrame(rows), pd.DataFrame(runs)


def test_classify_and_corners():
    pytest.importorskip('experiments.E123_kckc_form')
    rd, runs = _rd_runs({'family': 0.048, 'E110_A': 0.063, 'G': 0.049, 'T': 0.054, 'P': 0.060, 'S': 0.049},
                        corner_wirings={'S': ('real', 'P1', 'P3')})
    c = {a: X.classify(rd, runs, a) for a in ('G', 'T', 'P', 'S')}
    assert c['G']['cls'] == 'equal' and c['G']['at_least_keeps_gain']
    assert c['T']['cls'] == 'keeps the gain' and c['T']['at_least_keeps_gain']
    assert c['P']['cls'] == 'neither' and not c['P']['at_least_keeps_gain']
    assert c['S']['cls'] == 'corner (no class)' and c['S']['cls_before_corner'] == 'equal'
    assert not c['S']['at_least_keeps_gain']
    assert c['G']['gain_kept']['S4'] == pytest.approx((0.063 - 0.049) / (0.063 - 0.048), abs=0.05)
    g = X.gates(c)
    assert g['step1b']['open'] and g['step1c']['open'] and g['step2']['pair_with'] == ['G']
    assert not g['kp']['open']


def test_differ_is_e126s():
    E126 = pytest.importorskip('experiments.E126_kckc_survivor_ladder')
    # (E126's own arm names, which carry a dose)
    a, b = 'TY25', 'PO25'
    for vals, cw in (({a: 0.050, b: 0.055}, None), ({a: 0.0551, b: 0.0550}, None), ({a: 0.056, b: 0.055}, None),
                     ({a: 0.050, b: 0.055}, {b: ('real', 'P1', 'P2')}), ({a: 0.050, b: 0.055}, {b: ('real', 'P1')})):
        for seed in (3, 4, 5):
            rd, runs = _rd_runs(vals, cw, seed=seed)
            mine, theirs = X.differ(rd, runs, a, b), E126.differ(rd, runs, a, b)
            assert mine['verdict'] == theirs['verdict'], (vals, cw, seed)
    rd, runs = _rd_runs({'a': 0.050, 'b': 0.055})
    assert X.differ(rd, runs, 'a', 'b')['verdict'] == 'a lower'
    assert X.differ(rd, runs[runs['wiring'] != 'M3'], 'a', 'b')['verdict'] == 'incomplete'


def test_outcomes():
    def c(cls):
        return dict(complete=True, cls=cls, at_least_keeps_gain=cls in ('equal', 'keeps the gain'),
                    corner=cls.startswith('corner'))

    o = X.outcomes({'G': c('equal'), 'T': c('equal'), 'P': c('keeps the gain'), 'S': c('neither')},
                   {'S|P': dict(verdict='P lower')})
    assert o['a'] is True and o['b'] is None and o['c'] is None and o['d'] is False and o['e'] == 'P lower'
    assert not o['f']
    o = X.outcomes({'G': c('neither'), 'T': c('neither'), 'P': c('keeps the gain'), 'S': c('corner (no class)'),
                    'K': c('neither'), 'Y': c('equal')}, {'S|P': dict(verdict='no class (corner)')})
    assert o['a'] is False and o['b'] is True and o['c'] is False and o['d'] is True and o['f']
    o = X.outcomes({'G': c('keeps the gain'), 'T': c('equal'), 'K': c('keeps the gain')},
                   {'G|T': dict(verdict='not told apart')})
    assert o['c'] is True and o['f']
