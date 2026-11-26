"""Tests of experiments/E144_a1_malecns_calyces.py (the weaker antennal-lobe normalization and the rise of the
inhibition among KCs on the male CNS's two calyces): the design's points and names, the settings' starts and targets,
the builder's arguments against E133's, the step-0 matrix statistics on planted matrices, the read's bookkeeping and
the lock. No model is run."""
import json
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from experiments import E144_a1_malecns_calyces as X


FIT = {c: dict(params=dict(fixed_thr=100.0 + i, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.4,
                           wPNAPL_scale=1e-3, wKCKC_scale=10.0 * (i + 1)),
               target_on=200, target_off=440, size=3000, n_kcs=1500, name=f'fit_{c}')
       for i, c in enumerate(X.CALYCES)}


@pytest.fixture
def fits(monkeypatch):
    monkeypatch.setattr(X, 'calyx_fit', lambda c: FIT[c])
    monkeypatch.setattr(X, 'hemibrain_real_kckc', lambda: 26.43)
    return FIT


def test_points_and_names():
    pts = X.planned_points(['1'])
    assert len(pts) == 148 and len(set(pts)) == 148
    # the families' megamat reads first: L's and R's real wirings, their permutations, then H's
    assert [(p.calyx, p.wiring) for p in pts[:9]] == [('L', 'real'), ('R', 'real'), ('L', 'P1'), ('R', 'P1'),
                                                       ('L', 'P2'), ('R', 'P2'), ('L', 'P3'), ('R', 'P3'), ('H', 'real')]
    assert all(p.setting == 'family' and p.panel == 'megamat' for p in pts[:9])
    n = Counter((p.calyx, p.setting) for p in pts)
    for c in X.MALE:
        assert all(n[(c, s)] == 16 for s in ('family', 'B', 'A1')) and all(n[(c, s)] == 4 for s in ('F', 'D1.8', 'D5.8'))
    assert n[('L', 'Dabs')] == 4 and n[('R', 'Dabs')] == 0 and n[('H', 'Dabs')] == 0
    assert all(n[('H', s)] == 4 for s in ('family', 'F', 'B', 'D1.8', 'A1', 'D5.8'))
    assert sum(v for (c, s), v in n.items() if c in X.MALE and s in ('family', 'B', 'A1')) == 96
    tags = [X.point_tag(p) for p in pts]
    assert all(t.isascii() and all(ch.isalnum() or ch in '._-' for ch in t) for t in tags)
    assert X.point_tag(X.Point('L', 'A1', 'P2', 'megamat')) == 'L_A1_megamat_P2'
    assert all(X.parse_tag(t) == p for t, p in zip(tags, pts))
    for bad in (X.Point('H', 'A1', 'P1', 'megamat'), X.Point('R', 'Dabs', 'real', 'apl16'),
                X.Point('L', 'A1', 'P4', 'apl16'), X.Point('L', 'A1', 'real', 'kmix')):
        with pytest.raises(ValueError):
            X.point_tag(bad)
    assert X.seed_index('L', 'P3') == 3 and X.seed_index('H', 'real') == 0
    assert X.planned_points(['2'], calyces=['R']) == [X.Point('R', s, w, p) for s in ('D1.8', 'D5.8') for p in X.PANELS
                                                     for w in ('P1', 'P2', 'P3')]
    with pytest.raises(ValueError):
        X.planned_points(['2'], calyces=['H'])


def test_step2_calyces(monkeypatch, tmp_path):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    with pytest.raises(RuntimeError):
        X.step2_calyces()
    (tmp_path / 'classes.json').write_text(json.dumps(dict(calyces=dict(L=dict(a1_holds=True), R=dict(a1_holds=False)))))
    assert X.step2_calyces() == ('R',) and len(X.planned_points(['2'])) == 24


def test_settings(fits):
    assert X.SETTINGS['D1.8'].kckc == pytest.approx(1.8) and X.SETTINGS['D5.8'].kckc == pytest.approx(5.832)
    p, kws = X.setting_params('R', 'family')
    assert p == FIT['R']['params'] and kws == {}
    p, kws = X.setting_params('R', 'F')
    assert p == FIT['R']['params'] and kws['pn_inhadd'] == pytest.approx(1000.0)
    from experiments.E77_input_screen import inhsc_for
    assert kws['pn_inhsc'] == pytest.approx(inhsc_for(1000.0))
    p, _ = X.setting_params('L', 'A1')
    assert p['wAPLPN_scale'] == pytest.approx(0.4 * 0.155) and p['wKCKC_scale'] == pytest.approx(20.0 * 3.24)
    assert p['fixed_thr'] == FIT['L']['params']['fixed_thr'] and p['wAPLKC_scale'] == 0.1
    p, _ = X.setting_params('L', 'Dabs')
    assert p['wKCKC_scale'] == pytest.approx(3.24 * 26.43) and X.dabs_scale() == pytest.approx(85.63, abs=0.01)
    assert p['wAPLPN_scale'] == pytest.approx(0.4 * 0.155)
    with pytest.raises(ValueError):
        X.setting_params('R', 'Dabs')
    assert FIT['L']['params']['wKCKC_scale'] == 20.0      # the fit's own parameters are not changed in place


def test_read_targets(fits, monkeypatch):
    from experiments import E134_population_baseline as B
    monkeypatch.setattr(B, 'held_rate', lambda panel, lk=None: dict(megamat=0.0915, mixtures=0.054)[panel])
    lk = dict(populations={'kc_soma__validation2__wild_type': dict(response_rate=dict(mean_over_flies=0.05))})
    t = X.read_targets('L', 'megamat', 1500, 17, lk)
    assert t[1] == round(0.0915 * 1500 * 17) and t[0] == round(0.0915 * 2.2 * 1500 * 17)
    assert t[2] == pytest.approx(0.0915 * 2.2) and t[3] == pytest.approx(0.0915)
    assert X.read_targets('L', 'validation2', 1500, 22, lk)[1] == round(0.05 * 1500 * 22)
    assert X.read_targets('L', 'mixtures', 1500, 16, lk)[3] == pytest.approx(0.054)
    assert X.read_targets('L', 'apl16', 1500, 16, lk)[:2] == (440 * 8, 200 * 8)
    with pytest.raises(ValueError):
        X.read_targets('L', 'kmix', 1500, 8, lk)


def _fake_tables(kckc):
    W = pd.DataFrame(np.eye(2))
    return dict(orn=pd.DataFrame(np.ones((2, 2))), wnorm=dict(wAPLKC=pd.Series([1.0, 2.0]), wKCAPL='kcapl',
                                                              wAPLPN='aplpn', wPNAPL='pnapl'),
                breadth=pd.Series([3.0, 4.0]), W=W, kckc=kckc, calyx='malecns-left', meta={})


def test_build_simulator_is_e133s(monkeypatch):
    """With no input and no kwargs given, the simulator is built with E133's arguments exactly; with them, the input and
    the normalization's terms replace or join the connectome's kwargs and nothing else changes."""
    import pyolfsysm.mb_model as M
    from experiments import E133_calyx_transplants as E133
    from experiments import E10_joint_fit as J
    calls = []

    def fake(orn, **kw):
        calls.append(dict(orn=orn, **kw))
        mp, rv, odor_order = 'mp', 'rv', [0]

        def sim(p):
            return mp, rv, odor_order
        return sim, None
    monkeypatch.setattr(M, 'make_fixed_param_kc_simulator', fake)
    monkeypatch.setattr(J, 'breadth_weights', lambda w, b, k: ('bw', k))
    K = pd.DataFrame(np.ones((2, 2)))
    T = _fake_tables(K)
    monkeypatch.setattr(E133, 'calyx_tables', lambda c, seed: dict(T, seed=seed))
    p = dict(wKCKC_scale=7.7, fixed_thr=150.0)
    E133.build_simulator(T, p)
    X.build_simulator('L', 'P2', p)
    a, b = calls
    assert set(a) == set(b) and all(a[k] is b[k] or a[k] == b[k] for k in a if k != 'initial_params')
    assert a['initial_params'] == b['initial_params'] and b['_wKCKC'] is K and b['_wAPLKC'] == ('bw', 2.0)
    orn = pd.DataFrame(np.zeros((2, 3)))
    X.build_simulator('L', 'real', p, orn=orn, extra_kws=dict(pn_inhadd=1000.0, pn_inhsc=544.4))
    c = calls[-1]
    assert c['orn'] is orn and c['pn_inhadd'] == 1000.0 and c['pn_inhsc'] == 544.4
    assert {k for k in c if k not in a} == {'pn_inhadd', 'pn_inhsc'}
    assert all(c[k] is a[k] or c[k] == a[k] for k in a if k not in ('orn', 'initial_params'))
    # the hemibrain's matrix is olfsysm's own: no _wKCKC
    monkeypatch.setattr(E133, 'calyx_tables', lambda c, seed: dict(_fake_tables(None), calyx='hemibrain-argmax'))
    X.build_simulator('H', 'real', p)
    assert '_wKCKC' not in calls[-1]


def test_matrix_table_planted():
    """Four KCs of two types: KCs 0 and 1 share glomerulus A, 2 and 3 share B; all the weight on the sharing pairs,
    reciprocal, within type."""
    cols = pd.MultiIndex.from_tuples([('A', 1, 1), ('B', 2, 2)], names=['glomerulus', 'pn_id', 'bouton_id'])
    idx = pd.MultiIndex.from_tuples([(10, 'g'), (11, 'g'), (12, 'ab'), (13, 'ab')], names=[X.H.KC_ID_LEVEL, 'kc_type'])
    W = pd.DataFrame([[1, 0], [1, 0], [0, 1], [0, 1]], index=idx, columns=cols, dtype=float)
    ids = [10, 11, 12, 13]
    A = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float)
    M = pd.DataFrame(A, index=ids, columns=ids)
    r = X.matrix_table('planted', M, W, np.array(['g', 'g', 'ab', 'ab'], dtype=object))
    assert r['shared_glomerulus_pair_share'] == pytest.approx(4 / 12)
    assert r['shared_glomerulus_weight_share'] == 1.0 and r['shared_glomerulus_ratio'] == pytest.approx(3.0)
    assert r['reciprocity_weighted'] == 1.0 and r['reciprocal_share'] == 1.0 and r['row_sum_cv'] == 0.0
    assert r['within_type_weight_share'] == 1.0 and r['block_g_from_g'] == 0.5 and r['block_g_from_ab'] == 0.0
    assert r['row_sum_gini'] == pytest.approx(0.0, abs=1e-12)
    B2 = pd.DataFrame(np.array([[0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=float), index=ids,
                      columns=ids)
    r2 = X.matrix_table('one', B2, W, np.array(['g', 'g', 'ab', 'ab'], dtype=object))
    assert r2['shared_glomerulus_weight_share'] == 0.0 and r2['reciprocity_weighted'] == 0.0
    assert r2['row_sum_gini'] == pytest.approx(0.75) and r2['block_g_from_ab'] == 1.0


def test_read_point_bookkeeping(fits, monkeypatch, tmp_path):
    """A read's row, run info and provenance; a done tag is skipped; a failed read releases its claim."""
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    root, loc = tmp_path / 'res', tmp_path / 'loc'
    monkeypatch.setattr(X, 'out_dir', lambda: root)
    monkeypatch.setattr(R, 'require_lock', lambda: {})
    r = pd.DataFrame(np.array([[1.0, 0.0], [0.0, 0.0]]))
    s = pd.DataFrame(np.array([[3.0, 1.0], [0.0, 0.0]]))
    hold = dict(params=dict(FIT['L']['params'], fixed_thr=120.0), params_off=dict(FIT['L']['params'], wAPLKC_scale=0.0),
                on=(r, s), off=(r, s), thr_factor=1.2, apl_factor=0.9, at_bound_thr=False, at_bound_apl=False,
                n_sims=7, kws={}, history=[])
    monkeypatch.setattr(X, 'hold_point', lambda pt, lk=None: dict(hold=hold, start=FIT['L']['params'],
                                                                  kws=X.pn_terms(), orn=pd.DataFrame(np.ones((2, 2))),
                                                                  targets=(4, 2, 0.2, 0.1), n_kcs=2))

    class Built:
        meta, orn, info = {}, pd.DataFrame(np.ones((2, 2))), dict(model='C2')
    monkeypatch.setattr(X, 'build', lambda *a, **k: Built())
    monkeypatch.setattr(B, 'run_stages', lambda pt, b, p, poff, **k: {c: dict(responses=r, spike_counts=s)
                                                                       for c in ('on', 'off')})
    saved = []
    monkeypatch.setattr(R, 'save_step3_run', lambda tag, out, meta, orn, info, dirs: saved.append((tag, info, dirs)))
    monkeypatch.setattr(X.H, 'model_provenance', lambda **k: dict(prov=True))
    pt = X.Point('L', 'A1', 'P2', 'megamat')
    row = X.read_point(pt, 'viki-1', dirs=[loc])
    assert row['tag'] == 'L_A1_megamat_P2' and row['half'] == 'held' and row['variant'] == 'L_A1'
    assert row['rate_on'] == 0.25 and row['n_on'] == 1 and row['pn_inhadd'] == pytest.approx(1000.0)
    assert row['fixed_thr'] == 120.0 and row['off_wAPLKC_scale'] == 0.0 and row['start_fixed_thr'] == 101.0
    tag, info, dirs = saved[-1]
    assert info['wiring'] == 'L-P2' and info['calyx_wiring'] == 'P2' and info['seed_index'] == 2 and dirs == [loc]
    assert info['kckc_factor'] == 3.24 and info['setting'] == 'A1'
    t = pd.read_csv(root / f'reads_{X._host()}_viki-1.csv')
    assert list(t.tag) == ['L_A1_megamat_P2']
    prov = json.loads((root / f'provenance_read_{X._host()}_viki-1.json').read_text())
    assert prov['points']['L_A1_megamat_P2']['model'] == dict(prov=True)
    n = len(saved)
    assert X.read_point(pt, 'viki-1', done={'L_A1_megamat_P2'}) is None and len(saved) == n
    # a held claim skips; a failure releases the claim
    assert X.read_point(X.Point('L', 'B', 'real', 'apl16'), 'x', dirs=[loc]) is not None
    (root / 'claims' / 'L_B_megamat_real').mkdir(parents=True)
    assert X.read_point(X.Point('L', 'B', 'real', 'megamat'), 'x', dirs=[loc]) is None
    monkeypatch.setattr(B, 'run_stages', lambda *a, **k: {c: dict(responses=r, spike_counts=s * 2) for c in ('on', 'off')})
    with pytest.raises(RuntimeError):
        X.read_point(X.Point('R', 'B', 'real', 'megamat'), 'x', dirs=[loc])
    assert not (root / 'claims' / 'R_B_megamat_real').exists()


def test_lock_roundtrip(monkeypatch, tmp_path):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path / 'res')
    monkeypatch.setattr(X, 'committed_lock_path', lambda: tmp_path / 'E144_lock.json')
    d = dict(experiment=X.EXPERIMENT, definitions=dict(rise=3.24), amendments=[])
    monkeypatch.setattr(X, 'lock_dict', lambda: json.loads(json.dumps(d)))
    w = X.write_lock()
    assert w['definitions'] == dict(rise=3.24) and 'written' in w and X.require_lock()['experiment'] == X.EXPERIMENT
    assert X.write_lock() == w
    d['definitions'] = dict(rise=3.0)
    with pytest.raises(RuntimeError):
        X.require_lock()
    (tmp_path / 'res' / 'claims' / 'x').mkdir(parents=True)
    with pytest.raises(RuntimeError):
        X.write_lock()


def test_lanes(monkeypatch, tmp_path):
    monkeypatch.setattr(X, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(X, 'require_lock', lambda: {})
    monkeypatch.setattr(X, 'clean_commits', lambda: {})
    monkeypatch.setattr(X, 'still_clean', lambda cc: None)
    with pytest.raises(RuntimeError):
        X.part_read('x')
    (tmp_path / 'step0.json').write_text(json.dumps(dict(identity_passes=True, glomeruli_pass=True)))
    got = []
    monkeypatch.setattr(X, 'read_point', lambda pt, lane, reclaim, done: got.append(pt))
    X.part_read('x', reverse=True)
    assert len(got) == 148 and got[0] == X.Point('L', 'Dabs', 'real', 'apl16') and got[-1].setting == 'family'


def test_write_variants(monkeypatch, tmp_path):
    arch, res = tmp_path / 'arch', tmp_path / 'res'
    (arch / 'read').mkdir(parents=True)
    res.mkdir()
    monkeypatch.setattr(X, 'out_dir', lambda: res)
    on = dict(fixed_thr=110.0, wAPLKC_scale=0.15, wKCAPL_scale=8e-5, wAPLPN_scale=0.044, wPNAPL_scale=3e-4,
              wKCKC_scale=25.0)
    (arch / 'read' / 'L_A1_megamat_P2_run.json').write_text(json.dumps(dict(
        model='C2', wiring='L-P2', panel='megamat', reading='read', tag='L_A1_megamat_P2', variant='L_A1', params=on,
        params_off=dict(on, wAPLKC_scale=0.0, wAPLPN_scale=0.0), start=dict(on, fixed_thr=157.2),
        pn_kws=dict(pn_inhadd=1000.0, pn_inhsc=544.4), fit='malecns-left', kckc_factor=3.24, aplpn_factor=0.155)))
    df = X.write_variants(arch)
    r = df.iloc[0]
    assert len(df) == 1 and (res / 'variants.csv').exists() and r.calyx == 'L' and r.wiring_seed == 2
    assert r.off_wAPLKC_scale == 0.0 and r.start_fixed_thr == 157.2 and r.pn_inhadd == 1000.0 and r.kind.isascii()


def test_score_classes_and_price():
    """E144's classify on a calyx's four wirings against its family (E136's rule with the wirings given), a single
    wiring's bar the spread, the price against the family, and the outcomes' letters."""
    from experiments import E136_stage_population_fit as F
    from experiments import E144_score as S
    ws = ('real', 'P1', 'P2', 'P3')
    st = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
    spreads = dict(stage=0.02, **{s: 0.05 for s in st}, **{'/'.join(k): 0.05 for k in F.SCORED_PANELS})
    fam = {w: dict(stage=1.0, **{s: 1.0 for s in st}) for w in ws}
    a1 = {w: dict(stage=0.9, **{s: 0.98 for s in st}) for w in ws}
    pv = {'L_family': {k: {w: 1.0 for w in ws} for k in F.SCORED_PANELS},
          'L_A1': {k: {w: 1.01 for w in ws} for k in F.SCORED_PANELS}}
    rec = S.classify_on(F, 'L_A1', 'L_family', ws, st, {'L_family': fam, 'L_A1': a1}, pv, spreads, {}, None)
    assert rec['cls'] == 'holds' and rec['megamat']['change'] == pytest.approx(-0.1) and rec['bar'] == pytest.approx(0.02)
    pv['L_A1'][F.HOLDS_PANELS[0]] = {w: 1.2 for w in ws}
    assert S.classify_on(F, 'L_A1', 'L_family', ws, st, {'L_family': fam, 'L_A1': a1}, pv, spreads, {},
                         None)['cls'] == 'megamat only'
    flagged = S.classify_on(F, 'L_A1', 'L_family', ws, st, {'L_family': fam, 'L_A1': a1}, pv, spreads,
                            {'L_A1': {'P2': 'rate_on'}}, None)
    assert flagged['cls'] == 'no class' and 'P2' in flagged['reasons'][0]
    one = S.classify_on(F, 'L_D1.8', 'L_family', ['real'], st, {'L_family': fam, 'L_D1.8': {'real': a1['real']}},
                        {'L_family': pv['L_family'], 'L_D1.8': {k: {'real': 1.0} for k in F.SCORED_PANELS}}, spreads,
                        {}, None)
    assert one['cls'] == 'holds' and one['bar'] == pytest.approx(0.02)
    parts = {'L_A1': {w: {'z8_pattern': 0.07, 'halfA_pattern': 0.07, 'halfB_pattern': 0.07} for w in ws},
             'L_family': {w: {'z8_pattern': 0.06, 'halfA_pattern': 0.06, 'halfB_pattern': 0.06} for w in ws}}
    from experiments import E139_normalization_ladder as E139
    if set(E139.PRICE_TARGETS) <= {'z8', 'halfA', 'halfB'}:
        pr = S.price_on(parts, 'L_A1', 'L_family', ws)
        assert pr['price'] is True and pr['price_value'] == pytest.approx(0.01)
    inside = {s: dict(inside=True) for s in S.FAMILY_INSIDE_STAGES}
    cal = {'L': dict(a1_holds=True, a1_megamat_falls=True, b_holds=False, family_inside=inside),
           'R': dict(a1_holds=False, a1_megamat_falls=True, b_holds=False,
                     family_inside=dict(inside, pn_bouton=dict(inside=False)))}
    o = S.outcome_letters(cal, {'A1': dict(class_differs=False)})
    assert o == dict(a=False, b=False, c=True, d=False, e=True, f=False)
    cal['R']['a1_holds'] = None
    assert S.outcome_letters(cal, {})['c'] is None
