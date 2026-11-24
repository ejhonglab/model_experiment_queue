#!/usr/bin/env python3
"""Tests for E133 (experiments/E133_calyx_transplants.py; model_experiment_queue): the pair groups and their counts,
the grouped errors' identity, the glomerulus and responder readings on made-up tables, step 1's class and the guard on
made-up values, the tags, plans, lanes and claims, the step-2 gate, and that steps 0 and 1 never load the builders. No
simulations. Needs the shared directory (E10 reads the data-target registry at import); skipped without it."""

import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
E = pytest.importorskip('experiments.E133_calyx_transplants')


def _pairs(odors=None):
    od = list(H.APL16_ODORS if odors is None else odors)
    return pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(od) for b in od[i + 1:]], names=['odor_a', 'odor_b'])


# ---------------------------------------------------------------------------
# the pair groups

def test_pair_group_counts():
    """On the 16 odors with 2-but, paa and B-myr: 120 pairs, 42 with any of the three, 78 with none, 39 with exactly
    one, 81 the rest, 3 among the three; the groups partition as they should."""
    m = E.pair_groups(_pairs())
    assert {g: int(v.sum()) for g, v in m.items()} == dict(p120=120, p42=42, p78=78, p39=39, p81=81, p3=3)
    assert np.array_equal(m['p42'], m['p39'] | m['p3']) and np.array_equal(m['p81'], m['p78'] | m['p3'])
    assert not (m['p42'] & m['p78']).any() and not (m['p39'] & m['p81']).any()
    # canonical order (sorted names) gives the same counts
    assert int(E.pair_groups(E._canonical_pairs())['p39'].sum()) == 39
    assert len(E._canonical_pairs()) == 120


def test_grouped_errors_add_up():
    """The pattern part squared is the pair-count-weighted mean of each split's two groups' squares."""
    rng = np.random.default_rng(0)
    err = pd.Series(rng.normal(size=120), index=_pairs())
    g = E.grouped_errors(err)
    assert g['p120'] == pytest.approx(float(np.sqrt(np.mean(err.to_numpy() ** 2))))
    for a, b, na, nb in (('p42', 'p78', 42, 78), ('p39', 'p81', 39, 81)):
        assert 120 * g['p120'] ** 2 == pytest.approx(na * g[a] ** 2 + nb * g[b] ** 2)


def test_group_patterns_match_pattern_error():
    """group_patterns' 120-pair value is E83's pattern error's rms; the data's effect on the 3 pairs is reported."""
    from experiments import E83_cross_wiring_fit as E83
    rng = np.random.default_rng(1)
    P = _pairs()
    off = pd.Series(rng.uniform(0.4, 0.8, 120), index=P)
    on = off + rng.uniform(0.0, 0.3, 120)
    target = pd.Series(rng.uniform(0.0, 0.3, 120), index=P)
    row = E.group_patterns({'on': on, 'off': off}, {'all': target})
    want = float(np.sqrt(np.mean(E83.pattern_error({'on': on, 'off': off}, target).to_numpy() ** 2)))
    assert row['all_pattern_p120'] == pytest.approx(want)
    assert sum(k.startswith('p3_effect_model_') for k in row) == 3
    assert sum(k.startswith('all_p3_effect_data_') for k in row) == 3
    assert row['p3_effect_model_2-but_paa'] == pytest.approx(float(on[('2-but', 'paa')] - off[('2-but', 'paa')]))


# ---------------------------------------------------------------------------
# step 0's model-table readings on made-up tables

def _tiny_wiring():
    rows = [(1, 0, 'g'), (1, 1, 'g'), (1, 2, 'g'), (2, 0, 'ab'), (2, 1, 'ab'), (3, 0, "a'b'")]
    idx = pd.MultiIndex.from_tuples(rows, names=[H.KC_ID_LEVEL, 'claw_id', H.KC_TYPE_LEVEL])
    cols = pd.MultiIndex.from_tuples([('DM1', 10, 0), ('DM1', 10, 1), ('VA2', 11, 0)],
                                     names=['glomerulus', 'pn_id', 'bouton_id'])
    # claws: kc1 DM1, DM1, VA2; kc2 DM1, VA2; kc3 VA2
    v = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0], [0, 0, 1], [0, 0, 1]], dtype=float)
    W = pd.DataFrame(v, index=idx, columns=cols)
    counts = dict(apl_kc=pd.Series([1., 2., 3., 4., 5., 6.], index=idx), kc_apl=pd.Series([1.] * 6, index=idx),
                  apl_pn=pd.Series([10., 20., 30.], index=cols), pn_apl=pd.Series([1., 1., 2.], index=cols))
    return W, counts


def test_odor_glomeruli_share_rule():
    orn = pd.DataFrame({'2h': [10.0, 1.0, 0.95, -5.0], 'paa': [0.0, 0.0, 3.0, 1.0]},
                       index=pd.Index(['A', 'B', 'C', 'D'], name='glomerulus'))
    g = E.odor_glomeruli(orn, odors=('2h', 'paa'))
    # 2h: positive sum 11.95, a tenth 1.195: only A; paa: sum 4, a tenth 0.4: C and D
    assert list(g['2h'].index) == ['A']
    assert list(g['paa'].index) == ['C', 'D']


def test_glomerulus_tables_on_made_up_wiring():
    W, counts = _tiny_wiring()
    rows, by = E.glomerulus_tables('hemibrain-argmax', {'paa': pd.Series([5.0], index=['DM1'])}, W=W, counts=counts)
    r = rows[0]
    assert (r['n_claws'], r['n_claws_g'], r['n_claws_ab'], r['n_claws_apbp']) == (3, 2, 1, 0)
    assert (r['kcs_2plus'], r['kcs_2plus_g'], r['kcs_1plus']) == (1, 1, 2)
    assert r['n_boutons'] == 2 and r['apl_kc'] == 1 + 2 + 4 and r['apl_pn'] == 30.0
    assert r['share_apl_pn'] == pytest.approx(30 / 60) and r['apl_pn_per_bouton_rel'] == pytest.approx(15 / 20)
    assert by[0]['glomerulus'] == 'DM1' and by[0]['n_claws'] == 3


# ---------------------------------------------------------------------------
# step 1's class and the guard

def test_class_step1_rules():
    hemi = [0.050, 0.052, 0.054, 0.056, 0.058, 0.060, 0.062]          # mean 0.056, median 0.056
    kept = [x + 0.008 for x in hemi]                                 # gap +0.008; 6 of 7 above 0.056
    assert E.class_step1(kept, hemi)['cls'] == 'kept'
    redealt = [x + 0.004 for x in hemi]                              # gap +0.004
    assert E.class_step1(redealt, hemi)['cls'] == 're-dealt'
    mixed = [x + 0.0055 for x in hemi]
    assert E.class_step1(mixed, hemi)['cls'] == 'mixed'
    # a gap of 0.007 or more with only 4 of 7 above the hemibrain's median is "mixed"
    lopsided = [0.050, 0.051, 0.052, 0.1, 0.1, 0.1, 0.1]
    c = E.class_step1(lopsided, hemi)
    assert c['gap'] >= E.KEPT_GAP and c['left_above_hemi_median'] == 4 and c['cls'] == 'mixed'
    # no class when a run lacks one; incomplete without seven
    nc = E.class_step1(kept, hemi, all_classed=False)
    assert nc['cls'] is None and nc['would_be'] == 'kept'
    assert E.class_step1(kept[:6] + [np.nan], hemi)['cls'] == 'incomplete'


def _made_up_step1(left_shift=0.008, no_class_tag=None):
    rows, runs = [], []
    for c, shift in (('hemibrain-argmax', 0.0), ('malecns-left', left_shift), ('malecns-right', 0.003)):
        for s in (0,) + E.SEEDS1:
            tag = E.tag_of(E.Point(1, c, s, None))
            base = 0.050 + 0.002 * s + shift
            r = dict(tag=tag, calyx=c, seed=s, readout='cal')
            for t in E.TARGETS_READ:
                r[f'{t}_pattern'] = base
                for g in E.GROUPS:
                    r[f'{t}_pattern_{g}'] = base
            rows.append(r)
            rows.append(dict(r, readout='raw'))
            runs.append(dict(tag=tag, calyx=c, seed=s, no_class=(tag == no_class_tag), corner=False, guard_max_count=False,
                             realism_flag='', **{b: True for b in E.BAND_COLS}))
    return pd.DataFrame(rows), pd.DataFrame(runs)


def test_verdict1_on_made_up_readings():
    rd, runs = _made_up_step1()
    v = E.verdict1(rd, runs)
    assert v['class'] == 'kept' and v['all']['gap'] == pytest.approx(0.008)
    assert v['all']['right']['R_minus_H'] == pytest.approx(0.003)
    assert v['all']['own']['malecns-left'] == pytest.approx(0.058)
    assert set(v['all']['groups']) == set(E.GROUPS)
    rd, runs = _made_up_step1(no_class_tag='S1_malecns-left_seed3')
    v = E.verdict1(rd, runs)
    assert v['class'] is None and v['all']['would_be'] == 'kept' and v['no_class_runs'] == ['S1_malecns-left_seed3']
    # a run of the right without a class does not take step 1's class away
    rd, runs = _made_up_step1(left_shift=0.002, no_class_tag='S1_malecns-right_seed3')
    assert E.verdict1(rd, runs)['class'] == 're-dealt'


def test_guard_flags():
    assert not E.guard_flags(80, 160, (40, 80))['guard_max_count']        # exactly twice is not above
    assert E.guard_flags(81, 10, (40, 80))['guard_max_count']
    assert E.guard_flags(10, 161, (40, 80))['guard_max_count']


def test_rel():
    assert E._rel(1.0620204552640171, 1.062020455264017) <= E.REL_TOL
    assert E._rel(0.0, 0.0) == 0.0 and E._rel(1.0, 1.1) == pytest.approx(0.1 / 1.1)


def test_beside_right():
    b = E.beside_right([0.06] * 7, [0.05] * 7, [0.07] * 7)
    assert b['R_minus_H'] == pytest.approx(0.01) and b['R_minus_L'] == pytest.approx(-0.01)
    assert b['right_above_hemi_median'] == 7


# ---------------------------------------------------------------------------
# tags, plans, lanes, claims, the gate

def test_tags_round_trip():
    for pt in (E.Point(1, 'malecns-left', 3, None), E.Point(2, 'hemibrain-argmax', 0, 'ALL'), E.Point(2, 'malecns-left', 2, 'C')):
        assert E.parse_tag(E.tag_of(pt)) == pt
    assert E.tag_of(E.Point(1, 'malecns-left', 1, None)) == 'S1_malecns-left_seed1'
    with pytest.raises(ValueError):
        E.tag_of(E.Point(1, 'malecns-left', 1, 'C'))


def test_plans():
    p1 = E.planned_step1()
    assert len(p1) == 21 and len(set(p1)) == 21
    # the left's and the hemibrain's 14 first, the right's 7 last
    assert {pt.calyx for pt in p1[:14]} == {'hemibrain-argmax', 'malecns-left'}
    assert [pt.calyx for pt in p1[14:]] == ['malecns-right'] * 7 and [pt.seed for pt in p1[14:]] == list(range(1, 8))
    assert p1[:2] == [E.Point(1, 'hemibrain-argmax', 1, None), E.Point(1, 'malecns-left', 1, None)]
    r1 = E.planned_step1(('malecns-left', 'hemibrain-argmax', 'malecns-right'), reverse=True)
    assert set(r1) == set(p1) and r1[0] == E.Point(1, 'malecns-left', 7, None) and r1[14] == E.Point(1, 'malecns-right', 7, None)
    assert [pt.calyx for pt in E.planned_step1(('malecns-right',))] == ['malecns-right'] * 7
    with pytest.raises(ValueError):
        E.planned_step1(seeds=(0,))
    p2 = E.planned_step2()
    assert len(p2) == 48 and {pt.arm for pt in p2} == {'C', 'S', 'A', 'B', 'K', 'ALL'}
    assert {pt.seed for pt in p2} == {0, 1, 2, 3} and {pt.calyx for pt in p2} == set(E.HOSTS)
    assert E.lane_name(('malecns-left', 'hemibrain-argmax')) == 'ml-hb'
    assert E.lane_name(('malecns-left',), (1,), True) == 'ml-s1-rev'


def test_claims(tmp_path):
    assert E.claim('S1_malecns-left_seed1', 'lane-a', root=tmp_path)
    assert not E.claim('S1_malecns-left_seed1', 'lane-b', root=tmp_path)
    assert E.claim('S1_malecns-left_seed1', 'lane-b', reclaim=True, root=tmp_path)
    rec = json.loads((tmp_path / 'claims' / 'S1_malecns-left_seed1' / 'claim.json').read_text())
    assert rec['lane'] == 'lane-b'
    E.release('S1_malecns-left_seed1', root=tmp_path)
    assert not (tmp_path / 'claims' / 'S1_malecns-left_seed1').exists()


def test_step2_gate(tmp_path, monkeypatch):
    monkeypatch.setattr(E, 'out_dir', lambda: tmp_path)
    assert 'step 1' in E.step2_gate()
    for cls, ok in (('re-dealt', False), ('kept', True), ('mixed', True), (None, False)):
        (tmp_path / 'verdict1.json').write_text(json.dumps({'class': cls}))
        assert (E.step2_gate() is None) == ok


def test_steps_0_and_1_do_not_load_the_builders():
    """E108's `_still_clean` stops a lane when a file it loaded changes: steps 0 and 1 (the module, its preloads) must
    not load experiments/E133_builders.py, so the builders can be committed while step-1 lanes run."""
    code = ('import sys; import experiments.E133_calyx_transplants as m; m._preload(); '
            "print('experiments.E133_builders' in sys.modules)")
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=str(H.HERE), timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip().splitlines()[-1] == 'False'
