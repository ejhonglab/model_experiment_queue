#!/usr/bin/env python3
"""Tests for E93 (the measured input for the glomeruli that the model's input leaves at zero;
experiments/E93_missing_glomeruli_input.py) and the hooks it added to E10_joint_fit.py (`read_orn_input` with further
glomeruli) and E90_model_stages.py (EXTRA_SUBSETS): the inputs (the 36 rows equal the lab's; the added rows equal their
sources; the unmeasured entries are 0; the controls' permutations; the slopes on made-up tables); part A's departure,
sign agreement, null and rule on synthetic tables; who responds on a made-up claw table; the gate, the equal rule and the
bounds' rule on made-up readings; and part 0's identity on the real wiring (slow, skipped with QUICK=1)."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
X = pytest.importorskip('experiments.E93_missing_glomeruli_input')

ODORS = list(H.APL16_ODORS)
QUICK = os.environ.get('QUICK') == '1'


def _shared() -> bool:
    try:
        return H.shared_dir().exists()
    except Exception:
        return False


def _drosolf() -> bool:
    import importlib.util
    return importlib.util.find_spec('drosolf') is not None


needs_drosolf = pytest.mark.skipif(not _drosolf(), reason="needs drosolf (Hallem's data)")


# ---------------------------------------------------------------------------
# the inputs

@needs_drosolf
def test_inputs_keep_the_lab_rows_exactly():
    lab = X.lab_input()
    assert lab.shape == (36, 16)
    for name in ('I0z', 'F', 'L', 'U', 'V', 'PL1', 'PO1', 'PL2', 'PO2'):
        orn = X.input_of(name)
        assert orn.shape == (46, 16), name
        assert list(orn.index[:36]) == list(lab.index) and list(orn.index[36:]) == list(X.ADDED), name
        assert orn.columns.equals(lab.columns), name
        assert np.array_equal(orn.iloc[:36].to_numpy(), lab.to_numpy()), name
    assert np.array_equal(X.input_of('I0').to_numpy(), lab.to_numpy())
    assert (X.input_of('I0z').iloc[36:].to_numpy() == 0).all()


@needs_drosolf
def test_added_rows_equal_their_sources():
    from pyolfsysm.mb_model import megamat_orn_deltas, validation2_orn_deltas
    lab = X.lab_input()
    ab = [H.odor_abbrev(c) for c in lab.columns]
    meg = [c for c in lab.columns if c[0] == 'megamat']
    val = [c for c in lab.columns if c[0] == 'validation2']
    assert [H.odor_abbrev(c) for c in meg] == list(H.MEGAMAT12) and [H.odor_abbrev(c) for c in val] == list(H.VALIDATION2_4)
    m, v = megamat_orn_deltas(), validation2_orn_deltas()
    hal = X.hallem_table()
    slope = X.hallem_scaling()['slope']
    for kind in ('F', 'U'):
        r = X.added_rows(kind)
        for g in ('VA4', 'VL1'):
            assert np.array_equal(r.loc[g, meg].to_numpy(dtype=float), m.loc[g, meg].to_numpy(dtype=float))
            assert (r.loc[g, val].to_numpy() == 0).all()
        assert np.array_equal(r.loc['VM4', val].to_numpy(dtype=float), v.loc['VM4', val].to_numpy(dtype=float))
        assert (r.loc['VM4', meg].to_numpy() == 0).all()
        assert np.array_equal(r.loc['DL4', val].to_numpy(dtype=float), v.loc['DL4', val].to_numpy(dtype=float))
        for j, o in enumerate(ab):
            f = float(slope[o]) if (kind == 'F' and o in slope.index) else 1.0
            if o in H.MEGAMAT12:
                assert r.loc['DL4'].iloc[j] == pytest.approx(float(hal.loc['DL4', o]) * f, rel=1e-12)
            for g in X.HALLEM_ONLY:
                want = float(hal.loc[g, o]) * f if o in hal.columns else 0.0
                assert r.loc[g].iloc[j] == pytest.approx(want, rel=1e-12, abs=0.0)
    # the bounds: L has the lab's four rows only, V VA4's only
    L, F, V = X.added_rows('L'), X.added_rows('F'), X.added_rows('V')
    for g in X.ADDED:
        if g in X.LAB_MEGAMAT + X.LAB_VALIDATION:
            assert np.array_equal(L.loc[g].to_numpy(), np.nan_to_num(X.lab_rows().loc[g].to_numpy(dtype=float)))
        else:
            assert (L.loc[g].to_numpy() == 0).all()
        assert np.array_equal(V.loc[g].to_numpy(), F.loc[g].to_numpy() if g == 'VA4' else np.zeros(16))
    # euc is not in Hallem's set
    assert 'euc' not in hal.columns


@needs_drosolf
def test_unmeasured_entries_are_zero():
    mm = X.measured_mask()
    assert int(mm.to_numpy().sum()) == 134 and mm.size == 160
    for kind in ('F', 'U', 'L', 'V'):
        r = X.added_rows(kind).to_numpy()
        assert (r[~mm.to_numpy()] == 0).all(), kind
    # the 26 unmeasured: VA4 and VL1 for the four validation odors, VM4 for the 12 megamat odors, the six Hallem rows for euc
    un = {(g, ab) for g in X.ADDED for ab, ok in zip(ODORS, mm.loc[g].to_numpy()) if not ok}
    want = ({(g, o) for g in ('VA4', 'VL1') for o in H.VALIDATION2_4} | {('VM4', o) for o in H.MEGAMAT12}
            | {(g, 'euc') for g in X.HALLEM_ONLY})
    assert un == want
    # measured entries of F are mostly nonzero (Hallem's and the lab's)
    assert (X.added_rows('F').to_numpy()[mm.to_numpy()] != 0).mean() > 0.9


def test_slopes_on_made_up_tables():
    assert X.slope_through_origin([1, 2, 3], [2, 4, 6]) == pytest.approx(2.0)
    assert X.slope_through_origin([1, 2, np.nan, 3], [2, 4, 5, 6]) == pytest.approx(2.0)
    assert X.slope_through_origin([1, -1], [1, 1]) == pytest.approx(0.0)
    assert np.isnan(X.slope_through_origin([0, 0], [1, 2]))
    rng = np.random.default_rng(0)
    gl = [f'g{i}' for i in range(6)]
    hal = pd.DataFrame(rng.uniform(-10, 100, (6, 3)), index=gl, columns=['a', 'b', 'c'])
    k = {'a': 0.5, 'b': 2.0, 'c': 1.25}
    lab = pd.DataFrame({o: hal[o].iloc[:4] * k[o] for o in k}).reindex(gl[:4] + ['x1', 'x2']).fillna(7.0)
    sl = X.hallem_scaling(lab=lab, hal=hal)
    for o, v in k.items():
        assert sl.loc[o, 'slope'] == pytest.approx(v, rel=1e-12)
        assert sl.loc[o, 'r'] == pytest.approx(1.0) and sl.loc[o, 'n_glomeruli'] == 4
    x = np.concatenate([hal[o].iloc[:4].to_numpy() for o in k])
    y = np.concatenate([hal[o].iloc[:4].to_numpy() * k[o] for o in k])
    assert sl.loc['pooled', 'slope'] == pytest.approx(float(x @ y / (x @ x)))


def test_permuted_controls():
    cols = pd.MultiIndex.from_tuples([('megamat' if o in H.MEGAMAT12 else 'validation2', f'{o} @ -3') for o in ODORS],
                                     names=['panel', 'odor'])
    rows = pd.DataFrame(np.arange(160, dtype=float).reshape(10, 16), index=pd.Index(X.ADDED, name='glomerulus'), columns=cols)
    for seed in (1, 2, 3):
        p = X.permute_labels(rows, seed)
        assert list(p.index) == list(rows.index)
        moved = [int(np.flatnonzero((rows.to_numpy() == p.iloc[i].to_numpy()).all(axis=1))[0]) for i in range(10)]
        assert sorted(moved) == list(range(10)) and all(m != i for i, m in enumerate(moved))
        q = X.permute_odors(rows, seed)
        for i in range(10):
            a, b = rows.iloc[i].to_numpy(), q.iloc[i].to_numpy()
            assert sorted(a[:12]) == sorted(b[:12]) and sorted(a[12:]) == sorted(b[12:])
        assert not np.array_equal(q.to_numpy(), rows.to_numpy())
    assert np.array_equal(X.permute_labels(rows, 1).to_numpy(), X.permute_labels(rows, 1).to_numpy())


# ---------------------------------------------------------------------------
# the E10 and E90 hooks

def test_read_orn_input_with_further_glomeruli(tmp_path):
    from experiments import E10_joint_fit as J
    lab = X.lab_input()
    ab = H.abbrev_columns(lab)
    ab.index.name = 'glomerulus'
    p36 = tmp_path / 'lab.csv'
    ab.to_csv(p36)
    out36 = J.read_orn_input(str(p36), lab)
    assert out36.equals(lab.astype(float))
    extra = pd.DataFrame(np.arange(32, dtype=float).reshape(2, 16) / 3.0, index=['VA4', 'VL1'], columns=ab.columns)
    p38 = tmp_path / 'more.csv'
    pd.concat([ab, extra]).rename_axis('glomerulus').to_csv(p38)
    out = J.read_orn_input(str(p38), lab)
    assert list(out.index) == list(lab.index) + ['VA4', 'VL1'] and out.columns.equals(lab.columns)
    assert np.array_equal(out.iloc[:36].to_numpy(), lab.to_numpy())
    assert np.allclose(out.iloc[36:].to_numpy(), extra.to_numpy(), rtol=0, atol=1e-12)
    bad = pd.concat([extra, ab]).rename_axis('glomerulus')
    bad.to_csv(tmp_path / 'bad.csv')
    with pytest.raises(ValueError):
        J.read_orn_input(str(tmp_path / 'bad.csv'), lab)
    pd.concat([ab, ab.iloc[:1]]).rename_axis('glomerulus').to_csv(tmp_path / 'dup.csv')
    with pytest.raises(ValueError):
        J.read_orn_input(str(tmp_path / 'dup.csv'), lab)


def test_e90_extra_subsets():
    from experiments import E90_model_stages as E90
    units = pd.Index(['VA4', 'DL1', 'DA1'], name='glomerulus')
    meta = dict(claw_glomerulus=np.array(['VA4', 'DL1', 'DA1', 'VA4']), bouton_glomerulus=np.array(['VA4', 'DA1']))
    try:
        E90.EXTRA_SUBSETS['input46'] = ['VA4', 'DL1']
        assert E90.subset_mask('claw', None, meta, 'input46').tolist() == [True, True, False, True]
        assert E90.subset_mask('bouton', None, meta, 'input46').tolist() == [True, False]
        assert E90.subset_mask('pn', units, None, 'input46').tolist() == [True, True, False]
    finally:
        E90.EXTRA_SUBSETS.pop('input46', None)
    with pytest.raises(ValueError):
        E90.subset_mask('claw', None, meta, 'input46')
    assert E90.subset_mask('claw', pd.Index(range(4)), meta, 'all').all()


# ---------------------------------------------------------------------------
# part A

def test_departure_on_synthetic_tables():
    rng = np.random.default_rng(1)
    od, gl = [f'o{i}' for i in range(6)], [f'g{i}' for i in range(5)]
    O = pd.DataFrame(rng.uniform(0.1, 1.0, (6, 5)), index=od, columns=gl)
    P = O * np.array([0.5, 1.0, 2.0, 3.0, 0.25])                    # a gain per glomerulus
    assert np.allclose(X.departure(P, O, 'glom').to_numpy(), 0.0)
    assert not np.allclose(X.departure(P, O, 'odor').to_numpy(), 0.0)
    Q = O.mul(np.array([0.5, 1, 2, 3, 4, 5]), axis=0)               # a gain per odor
    assert np.allclose(X.departure(Q, O, 'odor').to_numpy(), 0.0)
    assert np.allclose(X.departure(O * 7.0, O, 'global').to_numpy(), 0.0)
    # a dendrite response above the ORNs' at one entry, within its glomerulus: a positive departure there
    P2 = O.copy()
    P2.iloc[2, 3] *= 3.0
    d = X.departure(P2, O, 'glom')
    assert d.iloc[2, 3] > 0.5 and (d.drop(index='o2').iloc[:, 3] < 0).all()
    # negative values are clipped at 0 before the floor
    Pn = O.copy()
    Pn.iloc[0, 0] = -5.0
    assert np.isfinite(X.departure(Pn, O, 'glom').iloc[0, 0])
    with pytest.raises(ValueError):
        X.departure(P, O, 'bogus')


def test_sign_agreement_and_null():
    f = np.array([1.0, -1.0, 2.0, -0.5, 0.3, -0.2])
    assert X.sign_agreement(np.array([2.0, -1.0, 1.0, -3.0, 0.1, -9.0]), f) == 1.0
    assert X.sign_agreement(np.array([-2.0, 1.0, 1.0, -3.0, 0.0, -9.0]), f) == pytest.approx(3 / 6)   # a zero disagrees
    assert np.isnan(X.sign_agreement(np.array([]), np.array([])))
    # 30 entries in 5 odors, departures that agree everywhere: the null (permuted within odor) agrees less
    rng = np.random.default_rng(2)
    odors = np.repeat([f'o{i}' for i in range(5)], 6)
    fac = rng.normal(0, 1, 30)
    dep = np.sign(fac) * rng.uniform(0.1, 1.0, 30)
    mask = np.ones(30, dtype=bool)
    nul = X.null_agreement(pd.Series(dep), odors, fac, mask, n=300, seed=0)
    assert X.sign_agreement(dep[mask], fac[mask]) == 1.0
    assert nul.shape == (300,) and nul.max() <= 1.0 and nul.mean() < 0.8
    # permutation within odor keeps each odor's departures
    nul2 = X.null_agreement(pd.Series(np.where(odors == 'o0', 1.0, -1.0)), odors, np.where(odors == 'o0', 1.0, -1.0), mask, n=50)
    assert np.allclose(nul2, 1.0)


def test_a2_rule_on_made_up_agreements():
    def r(a, lo=0.42, hi=0.62):
        return dict(agree=a, null_p5=lo, null_p95=hi)
    assert X.a2_rule([r(0.7), r(0.65), r(0.5)]) == "in the list's direction"
    assert X.a2_rule([r(0.7), r(0.5), r(0.5)]) == 'neither'
    assert X.a2_rule([r(0.3), r(0.4), r(0.5)]) == 'opposite'
    assert X.a2_rule([r(0.3), r(0.7), r(0.5)]) == 'neither'
    assert X.a2_rule([r(0.62), r(0.42), r(np.nan)]) == 'neither'          # at the percentiles, not beyond
    assert X.a2_rule([r(0.7), r(0.3), r(0.7)]) == "in the list's direction"


def test_who_responds_on_a_made_up_claw_table(monkeypatch):
    lab = list(X.lab_input().index)
    gl = lab[:3] + ['VA4', 'DA4l', 'DA1']
    M = pd.DataFrame([[2, 1, 1, 0, 0, 0],      # no zero-filled claws
                      [1, 1, 0, 1, 0, 1],      # VA4 + DA1: share zero-filled 0.5, of the ten 0.25
                      [0, 0, 0, 1, 1, 0],      # VA4 + DA4l: share 1.0
                      [3, 0, 0, 0, 0, 1]],     # DA1: 0.25
                     index=pd.Index([10, 11, 12, 13], name=H.KC_ID_LEVEL), columns=gl)
    cls = X.kc_classes(M)
    assert cls.share_zero.tolist() == [0.0, 0.5, 1.0, 0.25]
    assert cls.share_ten.tolist() == [0.0, 0.25, 1.0, 0.0]
    assert cls.va4.tolist() == [False, True, True, False] and cls.da4l_or_da3.tolist() == [False, False, True, False]
    assert cls.bin_zero.tolist() == ['none', 'above 0.4', 'above 0.4', '0.2 to 0.4']
    idx = pd.MultiIndex.from_tuples([(k, 'g') for k in M.index], names=[H.KC_ID_LEVEL, 'kc_type'])
    on = pd.DataFrame(0.0, index=idx, columns=ODORS)
    off = on.copy()
    off.iloc[0, :4] = 1.0
    off.iloc[1, 0] = 1.0
    on.iloc[0, 0] = 1.0
    o = {'on': (on, on * 3), 'off': (off, off * 3)}
    t = X.kc_table(o, M)
    assert t.n_odors_off.tolist() == [4.0, 1.0, 0.0, 0.0] and t.n_odors_on.tolist() == [1.0, 0.0, 0.0, 0.0]
    rows = pd.DataFrame(X.who_responds_rows(t, unit='toy'))
    va4 = rows[(rows.kind == 'VA4_claw')].set_index('cls')
    assert va4.loc['True', 'n_kcs'] == 2 and va4.loc['True', 'rate_off'] == pytest.approx(1 / 32)
    assert va4.loc['False', 'rate_2h_off'] == pytest.approx(0.5) and va4.loc['True', 'rate_2h_off'] == pytest.approx(0.5)
    z = rows[rows.kind == 'zero_filled_share'].set_index('cls')
    assert z.loc['none', 'odors_per_kc_off'] == 4.0 and z.loc['above 0.4', 'silent_share_off'] == 0.5


def test_input_shares_and_distances_on_made_up_tables():
    cols = pd.MultiIndex.from_tuples([('megamat' if o in H.MEGAMAT12 else 'validation2', f'{o} @ -3') for o in ODORS],
                                     names=['panel', 'odor'])
    base = pd.DataFrame(10.0, index=['g1', 'g2'], columns=cols)
    full = pd.concat([base, pd.DataFrame(-5.0, index=['g3'], columns=cols), pd.DataFrame(20.0, index=['g4'], columns=cols)])
    claws = pd.Series({'g1': 1.0, 'g2': 3.0, 'g3': 2.0, 'g4': 1.0, 'g5': 9.0})
    sh = X.input_shares(full, base, claws)
    assert sh.index.tolist() == ODORS and np.allclose(sh.to_numpy(), 20.0 / 60.0)
    rng = np.random.default_rng(3)
    t = pd.DataFrame(rng.normal(size=(6, 16)), index=[f'g{i}' for i in range(6)], columns=cols)
    d = X.weighted_distances(t, None)
    assert len(d) == 66
    a, b = d.index[0]
    ab = H.abbrev_columns(t)
    assert d.iloc[0] == pytest.approx(1.0 - np.corrcoef(ab[a], ab[b])[0, 1])
    w = X.weighted_distances(t, pd.Series(2.0, index=t.index))    # a constant weight is no weight
    assert np.allclose(w.to_numpy(), d.to_numpy())


def test_half_agreement():
    rng = np.random.default_rng(4)
    prof = rng.normal(size=12)
    flies = prof[None, :] + 0.01 * rng.normal(size=(8, 12))
    assert X.half_agreement(flies, n_splits=50) > 0.99
    noise = rng.normal(size=(8, 12))
    assert X.half_agreement(noise, n_splits=200) < X.half_agreement(flies, n_splits=50) - 0.5
    assert np.isnan(X.half_agreement(flies[:1]))


# ---------------------------------------------------------------------------
# the gate, the equal rule and the bounds on made-up readings

TGT = ['z8', 'z6'] + [f'half{h}' for h in X.HALVES]


def _readings(pat: dict, raw: dict = None, mean: dict = None, pat66: dict = None, singles: dict = None) -> pd.DataFrame:
    """Made-up readings: per input a pattern part (the same against every target) on S, T, all8 and each wiring."""
    rows = []
    for name, p in pat.items():
        for unit in ('S', 'T', 'all8') + tuple(X.ALL_WIRINGS):
            pu = (singles or {}).get((name, unit), p)
            base = {f'{t}_pattern': pu for t in TGT}
            base.update({f'{t}_pattern66': (pat66 or {}).get(name, p) for t in TGT})
            rows.append(dict(input=name, unit=unit, readout='cal', mean_effect=(mean or {}).get(name, 0.178), **base))
            rows.append(dict(input=name, unit=unit, readout='raw', z8_loss=(raw or {}).get(name, 0.073), **base))
    return pd.DataFrame(rows)


def _runs(names, flagged=()) -> pd.DataFrame:
    rows = [dict(input=n, wiring=w, flagged=(n, w) in flagged, tag=f'{n}_{w}') for n in names for w in X.ALL_WIRINGS]
    return pd.DataFrame(rows)


def test_gate_on_made_up_readings(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    names = ['I0', 'F', 'PL1', 'PO1']
    rd = _readings({'I0': 0.050, 'F': 0.040, 'PL1': 0.049, 'PO1': 0.048})
    g = X.gate93(rd, _runs(names), flies_mean=0.178)
    assert g['complete'] and g['gain_all'] == pytest.approx(0.010) and g['n_cells_lower'] == 16
    assert all(g[k] for k in ('i', 'ii', 'iii', 'iv', 'v', 'vi')) and g['passes']
    # the controls gain as much: (v) fails
    g = X.gate93(_readings({'I0': 0.050, 'F': 0.040, 'PL1': 0.044, 'PO1': 0.048}), _runs(names), flies_mean=0.178)
    assert not g['v'] and not g['passes'] and g['controls']['PL'] == dict(n_draws=1, gains=[pytest.approx(0.006)],
                                                                          mean_gain=pytest.approx(0.006))
    # a flagged run of F anywhere, a raw-count loss above the slack, a mean effect off: (iii) fails
    assert not X.gate93(rd, _runs(names, flagged={('F', 'real')}), flies_mean=0.178)['iii']
    assert not X.gate93(_readings({'I0': 0.05, 'F': 0.04, 'PL1': 0.05, 'PO1': 0.05}, raw={'F': 0.079}), _runs(names),
                        flies_mean=0.178)['iii']
    assert not X.gate93(_readings({'I0': 0.05, 'F': 0.04, 'PL1': 0.05, 'PO1': 0.05}, mean={'F': 0.165}), _runs(names),
                        flies_mean=0.178)['iii']
    # lower on the 120 pairs, not on the 66 megamat pairs: (vi) fails, outcome (d)
    g = X.gate93(_readings({'I0': 0.05, 'F': 0.04, 'PL1': 0.05, 'PO1': 0.05}, pat66={'F': 0.06}), _runs(names), flies_mean=0.178)
    assert g['passes_i_to_v'] and not g['vi'] and not g['passes'] and g['passes_120_not_66']
    # a gain of 0.002 is below (i); only five of T's 16 cells lower fails (ii)
    assert not X.gate93(_readings({'I0': 0.050, 'F': 0.048, 'PL1': 0.05, 'PO1': 0.05}), _runs(names), flies_mean=0.178)['i']
    singles = {('F', w): 0.06 for w in X.T_WIRINGS[:3]}
    g = X.gate93(_readings({'I0': 0.050, 'F': 0.040, 'PL1': 0.05, 'PO1': 0.05}, singles=singles), _runs(names), flies_mean=0.178)
    assert g['n_cells_lower'] == 4 and not g['ii']
    # (iv): a gain of 0.004 from 0.050 is a squared reduction of 15.4 %, from 0.060 of 12.9 %
    assert X.gate93(_readings({'I0': 0.050, 'F': 0.046, 'PL1': 0.05, 'PO1': 0.05}), _runs(names), flies_mean=0.178)['iv']
    assert not X.gate93(_readings({'I0': 0.060, 'F': 0.056, 'PL1': 0.06, 'PO1': 0.06}), _runs(names), flies_mean=0.178)['iv']
    # without F's T runs the gate is incomplete
    g = X.gate93(rd, _runs(['I0']), flies_mean=0.178)
    assert not g['complete'] and not g['passes']


def test_equal_rule_on_made_up_readings():
    names = ['I0', 'F']
    eq = X.equal93(_readings({'I0': 0.050, 'F': 0.052}), _runs(names))
    assert eq['equal'] and not eq['better_everywhere'] and eq['max_abs_d_pattern'] == pytest.approx(0.002)
    assert not X.equal93(_readings({'I0': 0.050, 'F': 0.054}), _runs(names))['equal']
    assert X.equal93(_readings({'I0': 0.050, 'F': 0.048}), _runs(names))['better_everywhere']
    assert not X.equal93(_readings({'I0': 0.05, 'F': 0.05}, raw={'F': 0.079}), _runs(names))['equal']
    assert not X.equal93(_readings({'I0': 0.05, 'F': 0.05}, mean={'F': 0.19}), _runs(names))['equal']
    assert not X.equal93(_readings({'I0': 0.05, 'F': 0.05}), _runs(names, flagged={('F', 'M3')}))['equal']
    assert not X.equal93(_readings({'I0': 0.05, 'F': 0.05}), _runs(['I0']))['equal']


def test_bounds_rule_on_made_up_readings(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    rd = _readings({'I0': 0.050, 'F': 0.049, 'L': 0.0495, 'U': 0.046, 'V': 0.050})
    ch = X.bound_choices(rd)
    assert [c['bound'] for c in ch] == ['U'] * 4 and all(c['to_T'] for c in ch)
    assert X.bound_choices(_readings({'I0': 0.05, 'F': 0.049, 'L': 0.049, 'U': 0.048, 'V': 0.05}))[0]['to_T'] is False
    # U read on T against the other half: passes the gate and is below F there
    b = X.bounds93(rd, _runs(['I0', 'F', 'L', 'U', 'V']), flies_mean=0.178)
    assert all(x['read_on_T'] and x['passes'] and x['preferred_to_F'] for x in b)
    # without U's T runs nothing is read
    runs = _runs(['I0', 'F', 'L', 'V'])
    runs = pd.concat([runs, pd.DataFrame([dict(input='U', wiring=w, flagged=False, tag=f'U_{w}') for w in X.S_WIRINGS])])
    assert not any(x['read_on_T'] for x in X.bounds93(rd, runs, flies_mean=0.178))


# ---------------------------------------------------------------------------
# part 0 (slow)

@pytest.mark.skipif(QUICK or not _shared() or not _drosolf(), reason='slow (builds the real wiring\'s simulator with I0z); '
                    'skipped with QUICK=1 or without the shared directory')
def test_identity_on_the_real_wiring(tmp_path, monkeypatch):
    """I0 with the ten rows written out as zeros: the same wiring and APL weights (E85's build raises otherwise), the
    breadth recomputed for it equal to the lab input's, and the saved APL-on and APL-off runs reproduced exactly."""
    from experiments import E83_cross_wiring_fit as E83
    monkeypatch.setattr(X, 'local_dir', lambda part='': tmp_path)
    P = E83.wiring_parts('real')
    (sim, info), tag = X.build93('real', X.input_of('I0z'), 'I0z')
    assert info['W'].equals(P['W'])
    assert info['breadth'].index.equals(P['breadth'].index)
    assert np.array_equal(info['breadth'].to_numpy(dtype=float), P['breadth'].to_numpy(dtype=float))
    on = sim(dict(P['p']))
    off = sim(dict(P['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))
    assert E83._same({'on': on, 'off': off}, P['saved']) == (0, 0)


def test_point_provenance_added_to_a_parts_record(tmp_path, monkeypatch):
    """A point's provenance is added to provenance_B.json whether or not the part's dated record (written at its start,
    without 'points') is already there."""
    X = pytest.importorskip('experiments.E93_missing_glomeruli_input')
    import json
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    (tmp_path / 'provenance_B.json').write_text(json.dumps(dict(part='B', runs=[dict(time='t')])))
    X._provenance_add('B', 'B_F_real', dict(a=1))
    X._provenance_add('B', 'B_F_P1', dict(a=2))
    d = json.loads((tmp_path / 'provenance_B.json').read_text())
    assert d['runs'] == [dict(time='t')] and d['points'] == {'B_F_real': dict(a=1), 'B_F_P1': dict(a=2)}
