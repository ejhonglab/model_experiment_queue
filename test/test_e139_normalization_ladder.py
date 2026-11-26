"""Tests of experiments/E139_normalization_ladder.py (E139, the weaker antennal-lobe normalization without the calyx
refit, along its strength): the points' coordinates in E136's A1, the planned reads and a lane's order, the read through
E136's `read_at`, the variants' runs for E136's score (E139's reads, E102's saved ladder for F on the APL16 input), the
price on the APL effect, the holds, the rung named on S, A1's evaluations binned by rung and the step's signature. No
model is run."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from experiments import E136_stage_population_fit as F
from experiments import E139_normalization_ladder as L
from kcq import stage_readings as R


def test_point_z():
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    z = L.point_z('B', 1000.0)
    assert z['aplpn'] == pytest.approx(np.log(0.155)) and z['kckc'] == 0.0 and z['inhsc'] == 0.0
    assert PN_INHADD_DEFAULT * np.exp(z['inhadd']) == pytest.approx(1000.0)
    kw = F.pn_norm_kws(z)
    from experiments.E77_input_screen import inhsc_for
    assert kw['pn_inhadd'] == pytest.approx(1000.0) and kw['pn_inhsc'] == pytest.approx(inhsc_for(1000.0))
    assert L.point_z('K', 314.0)['kckc'] == pytest.approx(np.log(3.24))


def test_planned_points():
    p1, p2 = L.planned_points(['1']), L.planned_points(['2'])
    assert len(p1) == 96 and len(p2) == 96
    assert {(p.setting, p.rung) for p in p2} == {('B', 314.0), ('B', 1000.0), ('K', 1000.0)}
    assert {p.panel for p in p1} == {'megamat', 'validation2', 'mixtures'} and 'apl16' in {p.panel for p in p2}
    assert len({L.point_tag(p) for p in p1 + p2}) == 192
    assert L.point_tag(L.Point('F', 314.0, 'megamat', 'P3')) == 'F_inhadd314_megamat_P3'
    assert p1[0] == L.Point('F', 100.0, 'megamat', F.WIRINGS[0])


def test_variants_and_tags():
    v = L.variants()
    assert v == ['F_inhadd100', 'F_inhadd314', 'F_inhadd1000', 'F_inhadd3141', 'B_inhadd314', 'B_inhadd1000',
                 'K_inhadd1000']
    for name in v:
        s, r = L.parse_variant(name)
        assert L.variant(s, r) == name
        pt = L.Point(s, r, 'validation2', 'M2')
        assert L.point_tag(pt) == f'{name}_validation2_M2'
    for bad in ('A1', 'family', 'X_inhadd10', 'F314'):
        with pytest.raises(ValueError):
            L.parse_variant(bad)
    assert L.variants(['1']) == v[:4]


def test_read_point_calls_read_at(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(L, 'out_dir', lambda: tmp_path / 'res')
    monkeypatch.setattr(L, 'read_dirs', lambda: [tmp_path / 'loc', tmp_path / 'arch'])
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='search', rung=None, outcome='e'))

    def fake(tag, arm, z, wiring, panel, half, lane, root, dirs, reclaim=False, done=None, extra=None, label='',
             script=''):
        seen.update(tag=tag, arm=arm, z=z, wiring=wiring, panel=panel, half=half, lane=lane, root=root, dirs=dirs,
                    reclaim=reclaim, extra=extra, label=label, script=script)
        return 'row'
    monkeypatch.setattr(F, 'read_at', fake)
    pt = L.Point('B', 1000.0, 'apl16', 'M1')
    assert L.read_point(pt, 'disco-smoke') == 'row'
    assert seen == dict(tag='B_inhadd1000_apl16_M1', arm='A1', z=L.point_z('B', 1000.0), wiring='M1', panel='apl16',
                        half=None, lane='disco-smoke', root=tmp_path / 'res', dirs=[tmp_path / 'loc', tmp_path / 'arch'],
                        reclaim=False, extra=dict(variant='B_inhadd1000', setting='B', rung=1000.0), label='E139',
                        script=L.SCRIPT)
    assert L.read_point(pt, 'x', done={'B_inhadd1000_apl16_M1'}) is None
    with pytest.raises(ValueError):
        L.read_point(L.Point('B', 1000.0, 'kmix', 'M1'), 'x')
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='table', rung='R1', outcome='a'))
    with pytest.raises(RuntimeError, match='mode'):
        L.read_point(pt, 'x')


def test_part_read_order(monkeypatch):
    got = []
    monkeypatch.setattr(L, 'read_point', lambda pt, lane, reclaim, done: got.append(pt))
    monkeypatch.setattr(L, 'reads_table', lambda: pd.DataFrame(columns=['tag']))
    L.part_read('x', ['1', '2'])
    assert got == L.planned_points(['1', '2'])
    got.clear()
    L.part_read('x', ['2'], reverse=True)
    assert got == L.planned_points(['2'])[::-1]
    with pytest.raises(NotImplementedError, match='beside'):
        L.part_read('x', ['1', 'beside'])


def _write_e102(base: Path, step: str, wiring: str, prefix: str, n: int = 6, odors=('2h @ -3', 'paa @ -4')):
    rng = np.random.default_rng(0)
    base.mkdir(parents=True, exist_ok=True)
    for c in ('on', 'acute_off'):
        s = pd.DataFrame(rng.poisson(2.0, (n, len(odors))).astype(float), columns=list(odors))
        s.to_parquet(base / f'{prefix}{wiring}_held_{c}_spike_counts.parquet')
        (s > 1.5).astype(float).to_parquet(base / f'{prefix}{wiring}_held_{c}_responses.parquet')
    (base / f'{prefix}{wiring}_point.json').write_text(json.dumps(dict(n_on=3, target_on=3, n_off=5, target_off=5,
                                                                       at_bound_thr=False, at_bound_apl=False)))


def test_load_variant_runs(monkeypatch, tmp_path):
    """F on the APL16 input is E102's saved point (the local layout first, then the archive's 'B_<step>_' names), every
    other point E139's read; the family's identity point is E102's 'identity' folder; through E136's registry."""
    from kcq import harness as H
    here, shared = tmp_path / 'here', tmp_path / 'shared'
    monkeypatch.setattr(H, 'HERE', here)
    monkeypatch.setattr(H, 'shared_dir', lambda: shared)
    monkeypatch.setattr(H, 'apl16_orn_deltas', lambda: pd.DataFrame({'2h @ -3': [1.0], 'paa @ -4': [2.0]}, index=['DM1']))
    assert L.load_e102_run(314.0, 'P3') is None
    _write_e102(shared / 'data' / 'model_outputs' / L.E102 / 'B' / 'inhadd314', 'inhadd314', 'P3', 'B_inhadd314_')
    r = L.load_e102_run(314.0, 'P3')
    assert r.source == 'E102:B' and r.tag == 'B_inhadd314_P3_held' and r.input_panel == 'apl16'
    assert r.stages() == ('orn', 'kc_soma') and r.info['n_on'] == 3 and r.info['e102_step'] == 'inhadd314'
    _write_e102(here / 'outputs' / L.E102 / 'B' / 'identity', 'identity', 'P3', '')
    assert L.load_e102_run(None, 'P3').tag == 'B_identity_P3_held'
    seen = []
    monkeypatch.setattr(L, 'read_dirs', lambda: ['d'])

    def step3(tag, dirs=None):
        seen.append((tag, dirs))
        raise FileNotFoundError(tag)
    monkeypatch.setattr(R, 'load_step3_run', step3)
    assert L.load_variant_run(F.ReadPoint('F_inhadd314', 'apl16', 'P3', None)).tag == 'B_inhadd314_P3_held'
    assert L.load_variant_run(F.ReadPoint('B_inhadd314', 'apl16', 'P3', None)) is None
    assert L.load_variant_run(F.ReadPoint('F_inhadd314', 'megamat', 'P3', 'held')) is None
    assert seen == [('B_inhadd314_apl16_P3', ['d']), ('F_inhadd314_megamat_P3', ['d'])]
    monkeypatch.setattr(F, 'VARIANT_LOADERS', {})
    L.register_variants(F, ['F_inhadd314'])
    assert F.load_read(F.ReadPoint('F_inhadd314', 'apl16', 'P3', None)).tag == 'B_inhadd314_P3_held'
    with pytest.raises(ValueError):
        L.register_variants(F, ['A1'])


def _pairs(n_odors=5):
    od = [f'o{i}' for i in range(n_odors)]
    idx = pd.MultiIndex.from_tuples([(a, b) for i, a in enumerate(od) for b in od[i + 1:]], names=['odor_a', 'odor_b'])
    return od, idx


def test_apl_effect_parts(monkeypatch):
    """With every pair defined, the price's statistics are E83's `effect_scores`; a pair a silent odor leaves undefined
    is left out on both sides."""
    from experiments import E83_cross_wiring_fit as E83
    rng = np.random.default_rng(1)
    od, idx = _pairs()
    d = {'on': pd.Series(rng.uniform(0.5, 0.9, len(idx)), index=idx),
         'off': pd.Series(rng.uniform(0.3, 0.7, len(idx)), index=idx)}
    tg = {t: pd.Series(rng.normal(0.18, 0.05, len(idx)), index=idx) for t in L.PRICE_TARGETS}
    monkeypatch.setattr(E83, 'pair_distances', lambda o, readout, K=None: {c: E83._cs(v) for c, v in o['d'].items()})
    monkeypatch.setattr(E83, 'targets', lambda: tg)
    p = L.apl_effect_parts({'d': d})
    sc = E83.effect_scores({c: E83._cs(v) for c, v in d.items()}, tg)
    for t in L.PRICE_TARGETS:
        for k in ('loss', 'pattern', 'rho'):
            assert p[f'{t}_{k}'] == pytest.approx(sc[f'{t}_{k}'])
        assert p[f'{t}_amplitude'] ** 2 + p[f'{t}_shape'] ** 2 == pytest.approx(p[f'{t}_pattern'] ** 2, rel=1e-6)
    assert p['n_pairs'] == len(idx) and p['mean_effect'] == pytest.approx(sc['mean_effect'])
    d2 = {c: v.copy() for c, v in d.items()}
    silent = [i for i, (a, b) in enumerate(idx) if 'o0' in (a, b)]       # a silent odor leaves its pairs undefined
    d2['on'].iloc[silent] = np.nan
    p2 = L.apl_effect_parts({'d': d2})
    assert p2['n_pairs'] == len(idx) - len(silent) and np.isfinite(p2['z8_pattern'])
    keep = d['on'].index[[i for i in range(len(idx)) if i not in silent]]
    sc2 = E83.effect_scores({c: E83._cs(v.reindex(keep)) for c, v in d.items()}, {t: v.reindex(keep) for t, v in tg.items()})
    assert p2['z8_pattern'] == pytest.approx(sc2['z8_pattern'])


def test_price_of():
    fam = {w: {f'{t}_pattern': 0.05 for t in L.PRICE_TARGETS} for w in F.WIRINGS}
    up = {w: {f'{t}_pattern': 0.05 + (0.004 if (t == 'halfB' and w in F.T_WIRINGS) else 0.0) for t in L.PRICE_TARGETS}
          for w in F.WIRINGS}
    rec = L.price_of({'family': fam, 'V': up}, 'V', 'family')
    assert rec['complete'] and rec['price'] and rec['price_value'] == pytest.approx(0.004)
    assert rec['by_target']['halfB']['change_8'] == pytest.approx(0.002) and not rec['price_8']
    small = {w: {f'{t}_pattern': 0.052 for t in L.PRICE_TARGETS} for w in F.WIRINGS}
    rec = L.price_of({'family': fam, 'V': small}, 'V', 'family')
    assert rec['price'] is False and rec['price_value'] == pytest.approx(0.002)
    del small['M3']
    rec = L.price_of({'family': fam, 'V': small}, 'V', 'family')
    assert not rec['complete'] and rec['price'] is None and rec['missing'] == ['M3']


def test_holds_failed():
    import types
    reads = pd.DataFrame([dict(tag='F_inhadd100_megamat_P3', variant='F_inhadd100', rate_on=0.0915, rate_off=0.2,
                               frac_on_target=0.0915, frac_off_target=0.2, at_bound_thr=False, at_bound_apl=False),
                          dict(tag='F_inhadd100_megamat_P4', variant='F_inhadd100', rate_on=0.0915, rate_off=0.21,
                               frac_on_target=0.0915, frac_off_target=0.2, at_bound_thr=False, at_bound_apl=False),
                          dict(tag='F_inhadd100_megamat_M2', variant='F_inhadd100', rate_on=0.0915, rate_off=0.2,
                               frac_on_target=0.0915, frac_off_target=0.2, at_bound_thr=False, at_bound_apl=True)])
    run_ok = types.SimpleNamespace(source='E102:B', tag='ok', info=dict(n_on=3, target_on=3, n_off=5, target_off=5))
    run_bad = types.SimpleNamespace(source='E102:B', tag='bad', info=dict(n_on=3, target_on=4, n_off=5, target_off=5))
    bad = L.holds_failed('F_inhadd100', reads, {'P3': run_ok, 'P4': run_bad}, 0.005)
    assert bad == ['F_inhadd100_megamat_P4', 'F_inhadd100_megamat_M2', 'bad']
    assert L.holds_failed('B_inhadd314', reads, {}, 0.005) == []


def test_name_rungs():
    mega_fit = {L.variant('F', r): {w: dict(stage=1.0 - (0.1 if r == 314.0 else 0.0)) for w in F.WIRINGS}
                for r in L.RUNGS}
    mega_fit['K_inhadd1000'] = {w: dict(stage=0.9) for w in F.WIRINGS}
    out = L.name_rungs(mega_fit, L.variants(), F.S_WIRINGS)
    assert out['F']['named'] == 314.0 and out['F']['values']['314'] == pytest.approx(0.9) and set(out) == {'F', 'K'}


def test_rung_bins_and_a1_by_rung(tmp_path):
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    b = L.rung_bins(np.array([PN_INHADD_DEFAULT, 50.0, 60.0, 300.0, 1017.3, 2580.0, 9000.0]))
    assert list(b) == [PN_INHADD_DEFAULT, PN_INHADD_DEFAULT, 100.0, 314.0, 1000.0, 3141.0, 3141.0]
    zi = np.log(np.array([PN_INHADD_DEFAULT, 307.0, 290.0, 1017.0]) / PN_INHADD_DEFAULT)
    hist = pd.DataFrame(dict(eval=[0, 1, 2, 3], z_aplpn=[0.0, -1.0, -0.5, -1.86], z_kckc=[0.0, 0.5, 0.2, 1.17],
                             z_inhadd=zi, z_inhsc=[0.0, 0.1, 0.0, 0.09]))
    hist.to_csv(tmp_path / 'fit_history_A1.csv', index=False)
    rows = []
    for i in range(4):
        for half, add in (('fit', 0.0), ('held', 0.05)):
            v = 1.0 - 0.03 * i + add + (0.02 if i == 1 else 0.0)
            rows.append(dict(arm='A1', eval=i, measure='strength', half=half, stage_error=v, stage_error_downstream=v,
                             pn_dendrite_error=v, pn_bouton_error=v, kc_claw_error=v, kc_soma_error=v))
            rows.append(dict(rows[-1], measure='all', stage_error=9.0))
    pd.DataFrame(rows).to_csv(tmp_path / 'rescore_A1.csv', index=False)
    t = L.a1_by_rung(tmp_path)
    best = t[(t.statistic == 'best')].set_index(['bin_rung', 'half'])
    assert best.loc[(314.0, 'fit'), 'eval'] == 2 and best.loc[(314.0, 'held'), 'eval'] == 2
    assert best.loc[(314.0, 'fit'), 'n_evals'] == 2 and best.loc[(1000.0, 'held'), 'stage_error'] == pytest.approx(0.96)
    assert best.loc[(314.0, 'fit'), 'aplpn_factor'] == pytest.approx(np.exp(-0.5))
    med = t[(t.statistic == 'median') & (t.bin_rung == 314.0) & (t.half == 'fit')].iloc[0]
    assert med.stage_error == pytest.approx(np.median([0.99, 0.94]))
    assert set(t.columns) >= {'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma', 'stage_error_downstream'}


def test_signature():
    """A proportional step has no signature; a divisive step by the odor's total gives the strong odors a negative one."""
    rng = np.random.default_rng(2)
    x = rng.uniform(0.0, 2.0, (6, 8))
    s = rng.uniform(0.5, 2.0, 6)
    assert np.allclose(L._signature(x, s[:, None] * x), 0.0, atol=1e-12)
    tot = x.sum(axis=0)
    y = s[:, None] * x / (1.0 + tot / tot.mean())
    sig = L._signature(x, y)
    from scipy import stats
    assert stats.spearmanr(sig, tot)[0] < -0.5
    assert np.allclose(L._signature(x, 3.0 * y), sig)       # scale-free
