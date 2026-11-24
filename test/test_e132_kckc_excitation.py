#!/usr/bin/env python3
"""Tests for E132 (recurrent excitation among KCs with a pooled inhibition; experiments/E132_kckc_excitation.py), on
made-up data: the refusal of negative doses, the route (the combined matrix d HB + b U is the family's at b 0 and
E129's pooled arm's at d 0, entry for entry), the runner on a toy simulator that reads its KC-KC matrix from a closure
cell as olfsysm's does (the matrix swapped in and restored, 'wKCKC_scale' 1, 'wU_scale' not passed on; the family and
the pooled arm reproduced exactly), the pooled hold (the threshold the fit's, the pooled scale bisected, then the APL
scales), late firing and the parts' inputs, step 0's readings (where the excitation would land, the Perron root and its
uniform share, the slope toward zero), the lanes, claims and archive names, and the classes, "differ" rule and stop
rule on synthetic tables. No model simulation."""

import types as _types

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E132_kckc_excitation')
sparse = pytest.importorskip('scipy.sparse')


def _toy_types(n_g=24, n_apbp=12, n_ab=28):
    return np.array(['g'] * n_g + ["a'b'"] * n_apbp + ['ab'] * n_ab, dtype=object)


def _toy_matrix(types, seed=0):
    rng = np.random.default_rng(seed)
    n = len(types)
    M = rng.poisson(1.0, size=(n, n)).astype(float) * (rng.random((n, n)) < 0.3)
    np.fill_diagonal(M, 0.0)
    return M * (n / M.sum())


# ---------------------------------------------------------------------------
# the dose, the pooled neuron and the route

def test_require_nonnegative():
    X.require_nonnegative(0.0)
    X.require_nonnegative(26.4, 13.0)
    with pytest.raises(X.SignedDoseError):
        X.require_nonnegative(-2.0)
    with pytest.raises(NotImplementedError):
        X.require_nonnegative(-0.5, 10.0)
    for d, b in ((1.0, -1.0), (float('nan'), 0.0), (1.0, float('inf'))):
        with pytest.raises(ValueError):
            X.require_nonnegative(d, b)


def test_uniform_neuron_is_e129s():
    E129 = pytest.importorskip('experiments.E129_block_survivor_static')
    M = _toy_matrix(_toy_types())
    U = X.uniform_neuron(M)
    assert np.array_equal(U, E129.uniform_neuron(M))
    assert U.sum() == pytest.approx(M.sum()) and np.unique(U).size == 1


def test_combined_matrix_is_the_familys_and_e129s():
    """At b 0 the combined matrix is the family's (scale x csc(HB)) entry for entry, with HB's sparsity; at d 0 it is
    E129's pooled arm's (scale x csc(U)); a negative dose is refused."""
    HB = _toy_matrix(_toy_types())
    U = X.uniform_neuron(HB)
    s, b = 26.43217, 0.858171 * 26.43217
    A, F = X.combined_matrix(HB, U, s, 0.0), s * sparse.csc_matrix(HB)
    assert np.array_equal(A.indptr, F.indptr) and np.array_equal(A.indices, F.indices)
    assert np.array_equal(A.data, F.data)
    A, F = X.combined_matrix(HB, U, 0.0, b), b * sparse.csc_matrix(U)
    assert np.array_equal(A.indptr, F.indptr) and np.array_equal(A.indices, F.indices)
    assert np.array_equal(A.data, F.data)
    A = X.combined_matrix(HB, U, 5.0, b).toarray()
    assert np.allclose(A, 5.0 * HB + b * U, rtol=0, atol=0)
    with pytest.raises(X.SignedDoseError):
        X.combined_matrix(HB, U, -2.0, b)
    with pytest.raises(ValueError):
        X.combined_matrix(HB, U[:-1, :-1], 1.0, 1.0)


# ---------------------------------------------------------------------------
# the runner on a toy simulator

class _Toy:
    """A few KCs with olfsysm's somatic equations: dVm = (-Vm + spont + drive (odor on) - wKCKC @ a - APL) dt / taum,
    a spike above the threshold (Vm reset to 0), a(t) = (1 - dt/tau) a(t-1) + spikes(t-1) / tau, APL proportional to
    wAPLKC_scale x the mean of a; the KC-KC matrix read at every call from the closure cell 'wKCKC_unscaled' (CSC) times
    params['wKCKC_scale'], as olfsysm's fixed-parameter simulator does; odors in an internal order."""

    KEYS = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale')

    def __init__(self, types, n_od=4, seed=5, dt=1e-3):
        rng = np.random.default_rng(seed)
        self.types = types
        self.n = len(types)
        self.W = _toy_matrix(types, seed)
        self.cell = _types.SimpleNamespace(cell_contents=sparse.csc_matrix(self.W))
        self.drive = rng.gamma(3.0, 60.0, size=(self.n, n_od))
        self.spont = rng.uniform(5, 15, size=self.n)
        self.order = np.array([2, 0, 3, 1])[:n_od]
        self.dt, self.taum = dt, 0.01
        self.mp = _types.SimpleNamespace(kc=_types.SimpleNamespace(kc_kc_tau=0.05))
        self.rv = _types.SimpleNamespace(kc=_types.SimpleNamespace())
        self.idx = pd.MultiIndex.from_arrays([np.arange(100, 100 + self.n), types], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
        self.odors = [f'od{j} @ -3' for j in range(n_od)]
        self.scales_seen = []

    def simulate(self, params, spike_times=False):
        unknown = set(params) - set(self.KEYS)
        if unknown:
            raise ValueError(f'parameters {unknown}')
        dt, tau = self.dt, float(self.mp.kc.kc_kc_tau)
        self.scales_seen.append(float(params['wKCKC_scale']))
        Wk = (float(params['wKCKC_scale']) * self.cell.cell_contents).toarray()
        thr = float(params['fixed_thr']) + 2 * self.spont
        n_od = self.drive.shape[1]
        V, a, sp_prev = np.zeros((self.n, n_od)), np.zeros((self.n, n_od)), np.zeros((self.n, n_od))
        counts = np.zeros((self.n, n_od))
        spikes = []
        for k, t in enumerate(np.arange(-0.1, 0.75, dt)):
            on = 0.0 <= t < 0.5
            a = (1 - dt / tau) * a + sp_prev / tau
            dV = -V + self.spont[:, None] + (self.drive if on else 0.0) - Wk @ a
            dV = dV - float(params['wAPLKC_scale']) * 40.0 * a.mean(axis=0)[None, :]
            V = V + dV * dt / self.taum
            sp = V > thr[:, None]
            V[sp] = 0.0
            sp_prev = sp.astype(float)
            counts += sp_prev
            if spike_times:
                for i, j in zip(*np.nonzero(sp)):
                    spikes.append((i, int(np.flatnonzero(self.order == j)[0]), k, t))
        s = pd.DataFrame(counts[:, self.order], index=self.idx, columns=self.odors)
        r = (s >= 2).astype(float)
        if not spike_times:
            return r, s
        return r, s, pd.DataFrame(spikes, columns=['kc', 'odor', 'step', 'time_s'])


@pytest.fixture
def toy(monkeypatch):
    E129 = pytest.importorskip('experiments.E129_block_survivor_static')
    pytest.importorskip('pyolfsysm.tuning')
    types = _toy_types()
    m = _Toy(types)
    p = dict(fixed_thr=150.0, wAPLKC_scale=1.0, wKCAPL_scale=1.0, wAPLPN_scale=0.5, wPNAPL_scale=1.0, wKCKC_scale=3.0)
    saved = {'on': m.simulate(dict(p)), 'off': m.simulate(dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))}
    parts = dict(p=p, saved=saved, target_on=X.n_resp(saved['on'][0]), target_off=X.n_resp(saved['off'][0]),
                 size=int(saved['on'][0].size), folder='toy', tag='toy', orn=None)

    def build(w, **kw):
        cells = {'wKCKC_unscaled': m.cell, 'mp': _types.SimpleNamespace(cell_contents=m.mp),
                 'rv': _types.SimpleNamespace(cell_contents=m.rv),
                 'odor_order': _types.SimpleNamespace(cell_contents=m.order)}
        return m.simulate, dict(cells=cells)

    ids = pd.Index(m.idx.get_level_values(H.KC_ID_LEVEL))
    monkeypatch.setattr(X.E83, 'wiring_parts', lambda w: parts)
    monkeypatch.setattr(X.E83, 'build_simulator', build)
    monkeypatch.setattr(E129, 'hemibrain_matrix', lambda idx: pd.DataFrame(m.W.copy(), index=ids, columns=ids))
    monkeypatch.setattr(X, '_preload', lambda: None)
    X._CACHE.pop('runner', None)
    return m, parts


def test_runner_reproduces_the_family_and_the_pooled_arm(toy):
    m, parts = toy
    R = X.Runner('real')
    p = parts['p']
    W0 = m.cell.cell_contents
    for c, q in (('on', R.params(R.s, 0.0)), ('off', X._off(R.params(R.s, 0.0)))):
        assert X._ndiff(R(q), parts['saved'][c]) == dict(counts=0, responses=0)
    assert m.cell.cell_contents is W0
    assert m.scales_seen[-1] == 1.0
    b = 0.7 * R.s
    mine = R(X._off(R.params(0.0, b)))
    m.cell.cell_contents = sparse.csc_matrix(X.uniform_neuron(m.W))
    theirs = m.simulate(dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0, wKCKC_scale=b))
    m.cell.cell_contents = W0
    assert X._ndiff(mine, theirs) == dict(counts=0, responses=0)
    # the pooled neuron changes the run, and the cell is restored after a failing call
    assert X._ndiff(R(X._off(R.params(R.s, 0.3 * R.s))), parts['saved']['off'])['counts'] > 0
    with pytest.raises(X.SignedDoseError):
        R(R.params(-2.0, b))
    assert m.cell.cell_contents is W0
    with pytest.raises(ValueError):
        R(dict(R.params(R.s, 0.0), bogus=1.0))
    assert m.cell.cell_contents is W0


def test_runner_refuses_another_matrix(toy, monkeypatch):
    m, parts = toy
    E129 = pytest.importorskip('experiments.E129_block_survivor_static')
    ids = pd.Index(m.idx.get_level_values(H.KC_ID_LEVEL))
    monkeypatch.setattr(E129, 'hemibrain_matrix', lambda idx: pd.DataFrame(m.W * 1.5, index=ids, columns=ids))
    with pytest.raises(RuntimeError):
        X.Runner('real')


def test_pooled_hold(toy):
    m, parts = toy
    R = X.Runner('real')
    res = X.pooled_hold(R, 0.0, 'toy U0')
    assert res['first_stage'] == 'wU_scale'
    for q in (res['params'], res['params_off']):
        assert q['fixed_thr'] == parts['p']['fixed_thr'] and q['wKCKC_scale'] == 0.0
        assert q['wU_scale'] == pytest.approx(res['pooled_factor'] * R.s)
    assert res['params_off']['wAPLKC_scale'] == 0.0 and res['params_off']['wAPLPN_scale'] == 0.0
    assert res['params']['wAPLPN_scale'] == parts['p']['wAPLPN_scale']
    for k in X.E83.TIED_ON:
        assert res['params'][k] == pytest.approx(parts['p'][k] * res['apl_factor'])
    assert res['held_off'] and res['held_on'] and not res['corner']
    assert X.n_resp(res['off'][0]) == res['n_off']
    hi = X.pooled_hold(R, 2.0, 'toy UI2')
    assert hi['params_off']['wKCKC_scale'] == 2.0
    assert hi['pooled_factor'] < res['pooled_factor']
    with pytest.raises(X.SignedDoseError):
        X.pooled_hold(R, -2.0, 'toy UX2')
    with pytest.raises(X.SignedDoseError):
        X.threshold_hold(R, -2.0, 'toy X2')


def test_run_realism_reproduces_and_reads_late_firing(toy):
    m, parts = toy
    R = X.Runner('real')
    res = X.pooled_hold(R, 2.0, 'toy UI2')
    row, spikes = X.run_realism(R, res, 2.0)
    for c in X.CONDS:
        n = np.zeros(res[c][1].shape)
        np.add.at(n, (spikes[c]['kc'].to_numpy(), spikes[c]['odor'].to_numpy()), 1.0)
        assert np.array_equal(n, res[c][1].to_numpy())
        assert 0.0 <= row[f'late_share_{c}'] <= 1.0 and row[f'n_spikes_{c}'] == len(spikes[c])
        assert row[f'kckc_in_mean_{c}'] > 0 and row[f'pooled_in_mean_{c}'] > 0
    bad = dict(res, on=(res['on'][0], res['on'][1] + 1.0))
    with pytest.raises(RuntimeError):
        X.run_realism(R, bad, 2.0)


def test_held_ok_and_late_share():
    assert X.held_ok(4000, 4000) and X.held_ok(4008, 4000) and not X.held_ok(4009, 4000)
    spk = pd.DataFrame(dict(kc=[0, 1, 2, 3], odor=[0, 0, 1, 1], step=[0, 1, 2, 3], time_s=[0.1, 0.4999, 0.5, 0.6]))
    assert X.late_share(spk) == 0.5
    assert np.isnan(X.late_share(spk.iloc[:0]))


def test_part_inputs():
    HB = _toy_matrix(_toy_types())
    n = pd.DataFrame(np.random.default_rng(2).poisson(2.0, (HB.shape[0], 3)).astype(float))
    r = X.part_inputs(n, HB, 4.0, 10.0, 100.0)
    U = X.uniform_neuron(HB)
    assert r['kckc_in_mean'] == pytest.approx((4.0 * HB @ n.to_numpy() / 0.5).mean())
    assert r['pooled_in_mean'] == pytest.approx((10.0 * U @ n.to_numpy() / 0.5).mean())
    assert r['kckc_in_over_thr'] == pytest.approx(r['kckc_in_mean'] / 100.0)


# ---------------------------------------------------------------------------
# step 0

def test_landing_by_hand():
    """Three KCs, two odors: KC 0 responds to a only, KC 1 to b only, KC 2 to neither; the input is s HB @ n / 0.5."""
    n = pd.DataFrame([[4.0, 0.0], [0.0, 6.0], [0.0, 1.0]], columns=['b @ -2', 'a @ -3'])
    r = (n >= 2).astype(float)
    HB = np.array([[0.0, 1.0, 0.5], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    lp = X.landing(n, r, HB, 3.0)
    assert list(lp[['odor_a', 'odor_b']].iloc[0]) == ['a', 'b']
    I = 3.0 * HB @ n.to_numpy() / 0.5            # columns: b, a
    row = lp.iloc[0]
    # under a (column 1) onto the KCs responding to b and not a (KC 0), over onto neither (KC 2)
    assert row['ratio_ab'] == pytest.approx(I[0, 1] / I[2, 1])
    assert row['ratio_ba'] == pytest.approx(I[1, 0] / I[2, 0])
    assert row['ratio'] == pytest.approx((row['ratio_ab'] + row['ratio_ba']) / 2)
    assert row['n_only_ab'] == 1 and row['n_only_ba'] == 1 and row['n_neither'] == 1


def test_landing_groups():
    pytest.importorskip('experiments.E90_model_stages')
    pairs = pd.DataFrame([('2h', 'IaA', 2.0), ('1-5ol', '1-6ol', 3.0), ('1-6ol', '2h', 5.0), ('B-myr', 'paa', 7.0),
                          ('2-but', 't2h', 11.0)], columns=['odor_a', 'odor_b', 'ratio'])
    g = X.landing_groups(pairs)
    assert g['ester'] == 2.0 and g['alcohol'] == 3.0 and g['ester_alcohol'] == 5.0 and g['other'] == 11.0
    assert g['all120'] == pytest.approx(np.mean([2, 3, 5, 7, 11]))


def test_perron_and_stability():
    types = _toy_types()
    n = len(types)
    U = X.uniform_neuron(np.eye(n))           # total N: every entry 1 / N, mean row sum 1
    r = X.perron(U)
    assert r['perron_root'] == pytest.approx(1.0) and r['uniform_share_right'] == pytest.approx(1.0)
    assert r['uniform_share_left'] == pytest.approx(1.0)
    D = np.diag(np.r_[np.ones(n - 1), 3.0])
    r = X.perron(D)
    assert r['perron_root'] == pytest.approx(3.0) and r['uniform_share_right'] == pytest.approx(1.0 / n)
    st = X.stability(_toy_matrix(types), types)
    assert list(st['block']) == ['all', 'g', 'apbp', 'ab']
    assert list(st['n_kcs']) == [n, 24, 12, 28]
    assert (st['perron_root'] <= st['spectral_radius'] + 1e-9).all() and (st['perron_root'] > 0).all()


def test_slopes():
    rows = []
    for kind, v in (('family', 0.048), ('H50', 0.055), ('E110_A', 0.063)):
        for unit in list(X.WIRINGS) + ['S4', 'T4', 'all8']:
            for ro in X.READOUTS:
                rows.append(dict(kind=kind, unit=unit, readout=ro, all_pattern=v, mean_effect=v, all_loss=v))
    rd = pd.DataFrame(rows)
    scales = {w: 30.0 for w in X.WIRINGS}
    sl = X.slopes(rd, scales).set_index(['unit', 'reading'])
    r = sl.loc[('S4', 'pattern')]
    assert r['slope_hi'] == pytest.approx((0.055 - 0.048) / -15.0)
    assert r['slope_lo'] == pytest.approx((0.063 - 0.055) / -15.0)
    assert r['at_m10'] == pytest.approx(0.063 + 10 * 0.008 / 15.0)
    assert r['at_m2_all'] == pytest.approx(0.063 + 2 * 0.015 / 30.0)
    assert np.isnan(sl.loc[('real', 'max_count_off'), 'at_s'])


# ---------------------------------------------------------------------------
# lanes, claims, archive names

def test_arms_and_tags():
    assert X.tag_of('U0', 'real') == 'S1_U0_real' and X.tag_of('UX5', 'M2') == 'S2_UX5_M2'
    assert X.planned(['P3', 'real'])[:5] == [('U0', 'real'), ('UI2', 'real'), ('UI5', 'real'), ('UI10', 'real'),
                                             ('U0', 'P3')]
    assert len(X.planned()) == 32
    with pytest.raises(X.SignedDoseError):
        X.planned(['real'], ['UX2'])
    for f in (X.part_step2, X.part_step3, X.part_step4):
        with pytest.raises(X.SignedDoseError):
            f()
    assert all(a['dose'] < 0 for k, a in X.ARMS.items() if a['step'] == 2)
    assert all(a['dose'] >= 0 for k, a in X.ARMS.items() if a['step'] == 1)
    assert X.dose0_of('UI5') == 'U0' and X.dose0_of('X5') == 'E110_A' and X.dose0_of('U0') is None


def test_claims(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    assert X.claim('S1_U0_real', 'lane')
    assert not X.claim('S1_U0_real', 'lane')
    assert X.claim('S1_U0_real', 'lane', reclaim=True)
    X.release('S1_U0_real')
    assert X.claim('S1_U0_real', 'lane')


def test_save_load_and_the_archive_names(tmp_path, monkeypatch):
    """model_pair_distances_long.csv's tag is the archive's file name, so the realism bands are scored from it."""
    A = pytest.importorskip('kcq.apl_effect')
    monkeypatch.setattr(X, 'local_dir', lambda part='': tmp_path / 'local' / part)
    monkeypatch.setattr(X, 'archive_dir', lambda: tmp_path / 'archive' / X.EXPERIMENT)
    monkeypatch.setattr(X.H, 'HERE', tmp_path / 'nowhere')
    idx = pd.MultiIndex.from_arrays([[1, 2], ['g', 'ab']], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
    s = pd.DataFrame([[0.0, 3.0], [2.0, 1.0]], index=idx, columns=['a', 'b'])
    spk = pd.DataFrame(dict(kc=[0, 1], odor=[1, 0], step=[3, 4], time_s=[0.1, 0.6]))
    X.save_run('S1_UI2_real', {'on': ((s >= 2).astype(float), s), 'off': ((s >= 2).astype(float), s)},
               {'on': spk, 'off': spk})
    back = X.load_run('S1_UI2_real')
    assert back['off'][1].equals(s) and back['on'][1].equals(s)
    assert (tmp_path / 'archive' / X.EXPERIMENT / 'step1' / 'S1_UI2_real_off_spike_times.parquet').exists()
    tag = X.distances_tag('S1_UI2_real')
    monkeypatch.setattr(A, 'ARCHIVE', tmp_path / 'archive')
    results = tmp_path / 'results' / X.EXPERIMENT
    for name in ('on', 'acute_off'):
        p = A.archive_responses(results, tag, name)
        assert p is not None and p.name == f'S1_UI2_real_{name}_responses.parquet'


# ---------------------------------------------------------------------------
# the classes, the "differ" rule and the stop rule

def _rd_runs(values, corner_wirings=None, seed=0, max_off=None, mean_effect=0.178):
    """A made-up readings table: per kind, a pattern part per wiring and per population for every target (calibrated),
    the mean APL effect under the readout, and the raw columns; a runs table with every run in the bands."""
    rng = np.random.default_rng(seed)
    rows, runs = [], []
    for kind, base in values.items():
        per_w = base if isinstance(base, dict) else {}
        b0 = float(np.mean([v for k, v in per_w.items() if k in X.WIRINGS])) if per_w else float(base)
        for unit in list(X.WIRINGS) + list(X.POPS) + ['all8']:
            v = per_w.get(unit, b0)
            r = dict(kind=kind, unit=unit, readout='cal', mean_effect=mean_effect)
            for t in X.targets():
                r[f'{t}_pattern'] = v + rng.normal(0, 1e-5)
            rows.append(r)
            rows.append(dict(kind=kind, unit=unit, readout='raw', mean_effect=0.12, all_loss=0.1))
        for w in X.WIRINGS:
            c = bool(corner_wirings and kind in corner_wirings and w in corner_wirings[kind])
            runs.append(dict(kind=kind, wiring=w, in_bands=True, rate_rho_on=0.5, corner=c, thr_below_60=False,
                             max_count_off=(max_off or {}).get(kind, 72)))
    return pd.DataFrame(rows), pd.DataFrame(runs)


def test_classify():
    pytest.importorskip('experiments.E123_kckc_form')
    pytest.importorskip('experiments.E108_paa_extrapolation')
    rd, runs = _rd_runs({'family': 0.048, 'E110_A': 0.063, 'U0': 0.058, 'UI2': 0.0495, 'UI5': 0.044, 'UI10': 0.054},
                        corner_wirings={'UI10': ('real', 'P1', 'P3')})
    c = {a: X.classify(rd, runs, a, 0.178) for a in ('U0', 'UI2', 'UI5', 'UI10')}
    assert c['U0']['cls'] == 'costs' and not c['U0']['keeps_gain']
    assert c['UI2']['cls'] == 'equal' and c['UI2']['keeps_gain']
    assert c['UI5']['cls'] == 'below the family'
    assert c['UI10']['cls'] == 'corner (no class)' and c['UI10']['cls_before_corner'] == 'costs'
    rd2, runs2 = _rd_runs({'family': 0.048, 'E110_A': 0.063, 'UI5': 0.044}, mean_effect=0.21)
    assert X.classify(rd2, runs2, 'UI5', 0.178)['cls'] == 'equal'
    assert X.classify(rd, runs[runs['wiring'] != 'M3'], 'U0', 0.178) == dict(arm='U0', complete=False)


def test_differ_is_e129s():
    E129 = pytest.importorskip('experiments.E129_block_survivor_static')
    for vals, cw in (({'a': 0.050, 'b': 0.055}, None), ({'a': 0.0551, 'b': 0.0550}, None),
                     ({'a': 0.056, 'b': 0.055}, None), ({'a': 0.050, 'b': 0.055}, {'b': ('real', 'P1', 'P2')})):
        for seed in (3, 4):
            rd, runs = _rd_runs(vals, cw, seed=seed)
            assert X.differ(rd, runs, 'a', 'b')['verdict'] == E129.differ(rd, runs, 'a', 'b')['verdict']


def test_stop_rule():
    vals = {'family': 0.050, 'E110_A': 0.063, 'U0': 0.058, 'UI5': 0.054, 'UI10': 0.0545,
            'UI2': dict({w: 0.0545 for w in X.WIRINGS}, real=0.060, P1=0.060, S4=0.05725)}
    rd, runs = _rd_runs(vals, max_off={'UI10': 150})
    st = {a: X.stop_rule(rd, runs, a) for a in ('U0', 'UI2', 'UI5', 'UI10')}
    assert st['U0']['passes'] is None and st['U0']['iv_dose0'] is None
    # UI5: within 0.005 of the family on S, below U0 on S and on all of S's wirings
    assert st['UI5']['i_pattern'] and st['UI5']['iv_dose0'] and st['UI5']['passes'] and not st['UI5']['candidate']
    # UI10: the largest count with APL off above twice the family's (150 > 144)
    assert st['UI10']['i_pattern'] and not st['UI10']['ii_max_count'] and not st['UI10']['passes']
    # UI2: lower than U0 on S's population but on only 2 of S's 4 wirings
    assert st['UI2']['dose0']['n_s_wirings_lower'] == 2 and not st['UI2']['iv_dose0'] and not st['UI2']['passes']
    rd, runs = _rd_runs({'family': 0.050, 'E110_A': 0.063, 'U0': 0.060, 'UI5': 0.0551})
    s = X.stop_rule(rd, runs, 'UI5')
    assert not s['i_pattern'] and not s['passes']
