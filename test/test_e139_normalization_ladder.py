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


@pytest.fixture(autouse=True)
def no_simulation(monkeypatch):
    """No test here runs the model: E83's simulator, E77's hold, E102's extended block and E136's fit raise unless a
    test plants its own."""
    from experiments import E77_input_screen as S
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E102_antennal_lobe_normalization as E

    def refuse(*a, **k):
        raise RuntimeError('a test reached a model simulation')
    for mod, name in ((E83, 'build_simulator'), (S, 'hold_rates'), (E, 'extended_block'), (E, 'run_extended'),
                      (F, 'fit_arm'), (F, 'determined'), (F, 'rescore')):
        monkeypatch.setattr(mod, name, refuse)


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
    monkeypatch.setattr(L, 'read_beside', lambda *a, **k: pytest.fail('no beside read without --step beside'))
    monkeypatch.setattr(L, 'reads_table', lambda: pd.DataFrame(columns=['tag']))
    L.part_read('x', ['1', '2'])
    assert got == L.planned_points(['1', '2'])
    got.clear()
    L.part_read('x', ['2'], reverse=True)
    assert got == L.planned_points(['2'])[::-1]


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
    # a count within the reads' relative tolerance holds (E102's 4068 against 4067); one beyond it, at a bound or
    # without counts fails
    run_near = types.SimpleNamespace(source='E102:B', tag='near', info=dict(n_on=1852, target_on=1852, n_off=4068,
                                                                           target_off=4067))
    run_far = types.SimpleNamespace(source='E102:B', tag='far', info=dict(n_on=1852, target_on=1852, n_off=4100,
                                                                         target_off=4067))
    run_bound = types.SimpleNamespace(source='E102:B', tag='bound', info=dict(n_on=3, target_on=3, n_off=5,
                                                                             target_off=5, at_bound_apl=True))
    run_none = types.SimpleNamespace(source='E102:B', tag='none', info={})
    run_read = types.SimpleNamespace(source='E139', tag='read', info=dict(n_on=1, target_on=9))
    bad = L.holds_failed('F_inhadd100', reads, {'P3': run_ok, 'P4': run_bad, 'M2': run_near, 'M3': run_far,
                                                'P1': run_bound, 'P2': run_none, 'real': run_read}, 0.005)
    assert bad == ['F_inhadd100_megamat_P4', 'F_inhadd100_megamat_M2', 'bad', 'far', 'bound', 'none']
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


def test_write_variants(monkeypatch, tmp_path):
    from experiments.E77_input_screen import inhsc_for
    arch, res = tmp_path / 'arch', tmp_path / 'res'
    for sub in ('read', 'fit_A0_inhadd314', 'extended/A0_inhadd314/B/inhadd314', 'extended/A1/B/inhadd1017.3'):
        (arch / sub).mkdir(parents=True)
    res.mkdir()
    monkeypatch.setattr(L, 'out_dir', lambda: res)
    on = dict(fixed_thr=30.0, wAPLKC_scale=1.0, wKCAPL_scale=6e-4, wAPLPN_scale=0.06, wPNAPL_scale=3e-3, wKCKC_scale=40.0)
    off = dict(on, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    d = dict(model='C2', wiring='P3', panel='megamat', reading='read', arm='A1', z=L.point_z('B', 314.0), params=on,
             params_off=off, tag='B_inhadd314_megamat_P3', variant='B_inhadd314',
             setting='B', rung=314.0, pn_inhadd=314.0, pn_inhsc=420.0, host='disco', lane='1', time='t')
    (arch / 'read' / 'B_inhadd314_megamat_P3_run.json').write_text(json.dumps(d))
    # step 3's read and its fit's evaluation record neither the setting's PN terms nor (the evaluation) the setting
    z0 = dict(aplpn=-1.65, kckc=0.0)
    (arch / 'read' / 'A0_inhadd314_validation2_M1_run.json').write_text(json.dumps(dict(
        model='C2', wiring='M1', panel='validation2', reading='read', arm='A0_inhadd314', z=z0, params=on,
        params_off=off, tag='A0_inhadd314_validation2_M1', variant='A0_inhadd314', setting='A0', rung=314.0)))
    (arch / 'fit_A0_inhadd314' / 'A0_inhadd314_e58_M1_run.json').write_text(json.dumps(dict(
        model='C2', wiring='M1', panel='megamat', reading='fit', arm='A0_inhadd314', z=dict(z0, kckc=0.41), params=on,
        params_off=off, tag='A0_inhadd314_e58_M1', eval=58)))
    # the beside step's runs: <wiring>_point.json, E102's part-B layout, the PN terms in pn_kws
    for sub, arm, rung, pn in (('extended/A0_inhadd314/B/inhadd314', 'A0_inhadd314', 314.0, None),
                               ('extended/A1/B/inhadd1017.3', 'A1', float('nan'), dict(pn_inhadd=1017.3, pn_inhsc=598.9))):
        (arch / sub / 'M1_point.json').write_text(json.dumps(dict(
            tag=L.beside_tag(arm, 'M1'), variant=arm, setting='A0' if arm != 'A1' else 'A1', rung=rung, arm=arm,
            wiring='M1', z=z0, params=on, params_off=off, **({'pn_kws': pn} if pn else {}))))
    df = L.write_variants(arch).set_index(['subdir', 'tag'])
    assert len(df) == 5 and (res / 'variants.csv').exists()
    b = df.loc[('read', 'B_inhadd314_megamat_P3')]
    assert b.variant == 'B_inhadd314' and b.off_wAPLKC_scale == 0.0 and b.z_aplpn == pytest.approx(np.log(0.155))
    assert b.pn_inhadd == 314.0 and b.pn_inhsc == 420.0 and 'APL->bouton x0.155' in b.kind
    for key in (('read', 'A0_inhadd314_validation2_M1'), ('fit_A0_inhadd314', 'A0_inhadd314_e58_M1'),
                ('extended/A0_inhadd314/B/inhadd314', 'M1')):
        r = df.loc[key]
        assert r.pn_inhadd == pytest.approx(314.0) and r.pn_inhsc == pytest.approx(inhsc_for(314.0)), key
        assert r.setting == 'A0' and r.rung == 314.0 and "step 3's refit" in r.kind, key
    f = df.loc[('fit_A0_inhadd314', 'A0_inhadd314_e58_M1')]
    assert f.reading == 'fit' and f.eval == 58 and f.part == 'fit_A0_inhadd314' and 'evaluation 58' in f.kind
    a = df.loc[('extended/A1/B/inhadd1017.3', 'M1')]
    assert a.panel == 'extended' and a.reading == 'extended' and a.part == 'extended' and a.pn_inhadd == 1017.3
    assert "E136's A1" in a.kind and "extended input" in a.kind
    assert L._truthy(True) and L._truthy('True') and not L._truthy(False) and not L._truthy('False')
    assert not L._truthy(np.nan)


# ---------------------------------------------------------------------------
# the beside step and step 3

def test_beside_points_and_names():
    bp = L.beside_points()
    assert len(bp) == 16 and bp[0] == ('B_inhadd1000', F.WIRINGS[0]) and bp[-1] == ('A1', F.WIRINGS[-1])
    assert L.beside_tag('A1', 'P3') == 'A1_extended_P3'
    assert L.fit_name(314.0) == 'A0_inhadd314' and L.fit_rung('A0_inhadd314') == 314.0
    assert L.fit_rung('B_inhadd314') is None and L.fit_rung('A0') is None and L.fit_rung('A1') is None
    assert (L._setting_of('A0_inhadd314'), L._rung_of('A0_inhadd314')) == ('A0', 314.0)
    assert L._setting_of('A1') == 'A1' and np.isnan(L._rung_of('A1'))
    arm = L.fit_arm_def(314.0)
    from experiments.E77_input_screen import inhsc_for
    assert arm.name == 'A0_inhadd314' and arm.searched == ('aplpn', 'kckc')
    assert arm.kws['pn_inhadd'] == pytest.approx(314.0) and arm.kws['pn_inhsc'] == pytest.approx(inhsc_for(314.0))
    assert L.beside_model('B_inhadd1000') == ('A1', L.point_z('B', 1000.0))


def _beside_fakes(monkeypatch, tmp_path, ok=True):
    from experiments import E98_apl_response_strength as E98
    from experiments import E102_antennal_lobe_normalization as E
    from kcq import harness as H
    got = {}
    monkeypatch.setattr(L, 'out_dir', lambda: (tmp_path / 'res').mkdir(parents=True, exist_ok=True) or tmp_path / 'res')
    monkeypatch.setattr(L, 'local_dir', lambda part='': (tmp_path / 'loc' / part).mkdir(parents=True, exist_ok=True)
                        or tmp_path / 'loc' / part)
    monkeypatch.setattr(L, 'archive_dir', lambda: tmp_path / 'arch')
    monkeypatch.setattr(F, 'require_e135', lambda: dict(mode='search', rung=None, outcome='e'))
    monkeypatch.setattr(R, 'require_lock', lambda: {})
    monkeypatch.setattr(F, '_commits', lambda: {'queue_commit': 'x'})
    monkeypatch.setattr(L, '_same_as_read', lambda name, wiring, hold: dict(tag='t', params_equal=True,
                                                                             counts_equal=True))
    monkeypatch.setattr(H, 'model_provenance', lambda orn_deltas=None: {})
    monkeypatch.setattr(E98, 'b1_input', lambda: (pd.DataFrame({'x': [1.0]}), pd.DataFrame()))
    r = pd.DataFrame(np.r_[np.ones(3), np.zeros(17)].reshape(10, 2), columns=['2h @ -3', 'paa @ -4'])
    p = dict(fixed_thr=30.0, wAPLKC_scale=1.0, wKCAPL_scale=6e-4, wAPLPN_scale=0.06, wPNAPL_scale=3e-3,
             wKCKC_scale=40.0)

    def hold_at(arm, z, wiring, panel, label='', lock=None):
        got['hold_at'] = dict(arm=arm, z=dict(z), wiring=wiring, panel=panel)
        kws = F.pn_norm_kws(z)
        hold = dict(params=dict(p), params_off=dict(p, wAPLKC_scale=0.0), on=(r, r), off=(r, r), thr_factor=1.1,
                    apl_factor=0.9, n_on=3, n_off=3, at_bound_thr=False, at_bound_apl=False, n_sims=9)
        return dict(hold=hold, kws=kws, orn=None, breadth='B', claw_thr=None, targets=(3, 3, 0.15, 0.15), p0=dict(p))
    monkeypatch.setattr(F, 'hold_at', hold_at)

    def extended_block(rec, d, wiring, step, held, breadth, commits, prefix=None, **ext):
        got['ext'] = dict(step=step, breadth=breadth, build=ext.get('build'), label=ext.get('label'))
        (d / 'b1_stimuli.csv').write_text('label\n')
        for c in ('on', 'acute_off'):
            (d / f'{wiring}_{c}_apl.csv').write_text('stimulus,win\n')
            (d / f'{wiring}_{c}_spike_counts.parquet').write_text('x')
        pd.DataFrame([dict(wiring=wiring, cond=c, pn_inhadd=step, pn_inhsc=1.0) for c in ('on', 'acute_off')]).to_csv(
            d / f'{wiring}_runs.csv', index=False)
        rec.update(extended_matches_held=ok, a0_identity_passes=True)
        return ok
    monkeypatch.setattr(E, 'extended_block', extended_block)
    return got


def test_read_beside(monkeypatch, tmp_path):
    """A beside read on planted model calls: held as E136 holds a point (hold_at on the APL16 input at the variant's
    coordinates), the extended block run at the held parameters with E136's breadth into E102's layout, the PN terms
    of runs.csv the variant's, the files copied to the archive with the point record last, the row in E139's reads;
    a failed check raises and releases the claim."""
    import json as _json
    got = _beside_fakes(monkeypatch, tmp_path)
    row = L.read_beside('A1', 'P3', 'lane', arm='A1', z=dict(aplpn=np.log(0.155), kckc=np.log(3.24),
                                                             inhadd=np.log(1017.3 / 31.4088), inhsc=np.log(1.0936)))
    assert got['hold_at']['panel'] == 'apl16' and got['hold_at']['arm'] == 'A1' and got['ext']['breadth'] == 'B'
    assert got['ext']['step'] == pytest.approx(1017.3) and callable(got['ext']['build'])
    d = tmp_path / 'loc' / 'extended' / 'A1' / 'B' / 'inhadd1017.3'
    a = tmp_path / 'arch' / 'extended' / 'A1' / 'B' / 'inhadd1017.3'
    rc = pd.read_csv(d / 'P3_runs.csv')
    assert np.allclose(rc.pn_inhsc, F.pn_norm_kws(got['hold_at']['z'])['pn_inhsc'])
    assert (rc.variant == 'A1').all()
    assert {p.name for p in a.iterdir()} == {p.name for p in d.iterdir()} and (a / 'P3_point.json').exists()
    rec = _json.loads((a / 'P3_point.json').read_text())
    assert rec['complete'] and rec['operating_point'].startswith("the wiring's own fit's APL16 counts")
    assert row['step'] == 'beside' and row['panel'] == 'extended' and row['same_as_apl16_read'] is True
    tab = pd.read_csv(tmp_path / 'res' / f'reads_{L._host()}_lane.csv')
    assert list(tab.tag) == ['A1_extended_P3'] and tab.frac_on_target.iloc[0] == 0.15
    assert L.read_beside('A1', 'P3', 'lane', done={'A1_extended_P3'}) is None
    assert L.read_beside('A1', 'P3', 'lane', arm='A1', z={}) is None          # claimed
    _beside_fakes(monkeypatch, tmp_path, ok=False)
    with pytest.raises(RuntimeError, match='check'):
        L.read_beside('B_inhadd1000', 'M2', 'lane')
    assert not (tmp_path / 'res' / 'claims' / 'B_inhadd1000_extended_M2').exists()
    assert not (tmp_path / 'arch' / 'extended' / 'B_inhadd1000').exists()


def test_part_read_with_beside(monkeypatch):
    got = []
    monkeypatch.setattr(L, 'read_point', lambda pt, lane, reclaim, done: got.append(('point', pt)))
    monkeypatch.setattr(L, 'read_beside', lambda name, w, lane, reclaim, done: got.append(('beside', (name, w))))
    monkeypatch.setattr(L, 'reads_table', lambda: pd.DataFrame(columns=['tag']))
    L.part_read('x', ['beside'])
    assert got == [('beside', bp) for bp in L.beside_points()]
    got.clear()
    L.part_read('x', ['2', 'beside'], reverse=True)
    assert got == ([('point', p) for p in L.planned_points(['2'])] + [('beside', bp) for bp in L.beside_points()])[::-1]


def test_beside_readings(monkeypatch, tmp_path):
    """The beside readings: E102's step readings of each variant's folder (the archive's), summarized as E102's
    b_steps rows, the low-concentration class against E102's identity point, E102's own rows beside."""
    from experiments import E102_antennal_lobe_normalization as E
    from kcq import harness as H
    monkeypatch.setattr(L, 'archive_dir', lambda: tmp_path / 'arch')
    monkeypatch.setattr(L, 'local_dir', lambda part='': tmp_path / 'loc' / part)
    rd = tmp_path / 'e102'
    rd.mkdir()
    monkeypatch.setattr(H, 'results_dir', lambda exp: rd)
    pd.DataFrame([dict(step='identity', n_wirings=8, v2='flatter', v2_held='flatter', low_on=0.0136, low_off=0.04,
                       low_flies_on=0.05, low_flies_off=0.11, low_class_on='about the same'),
                  dict(step='inhadd1000', n_wirings=8, v2='inside', v2_held='inside', low_on=0.0048, low_off=0.02,
                       low_flies_on=0.05, low_flies_off=0.11, low_class_on='further')]).to_csv(rd / 'b_steps.csv',
                                                                                              index=False)
    d = tmp_path / 'arch' / 'extended' / 'B_inhadd1000' / 'B' / 'inhadd1000'
    d.mkdir(parents=True)
    for w in ('P3', 'P4'):
        (d / f'{w}_point.json').write_text('{}')
    seen = {}

    def step_readings(tag, root=None, wirings=E.ALL_WIRINGS):
        seen.update(tag=tag, root=root, wirings=list(wirings))
        comps = pd.DataFrame([dict(comparison='V2, paa', reading_set='carrying', category_linear='flatter',
                                   model_linear_min=0.7, model_linear_max=0.8)])
        low = pd.DataFrame([dict(stimulus=E.LOW_CARRYING, wiring=w, rate_on=v, rate_off=2 * v, off_over_on=2.0,
                                 flies_on=0.05, flies_off=0.11) for w, v in (('P3', 0.004), ('P4', 0.006))])
        return dict(step=tag, wirings=list(wirings), e98=dict(V2='inside', V2_held='inside', V1='inside'),
                    e98_comparisons=comps, e98_wirings=pd.DataFrame([dict(wiring='P3', v=1)]),
                    e100=dict(W1='inside'), low=low, e101=pd.DataFrame([dict(ratio_id='r1', quantity='APL', held='x')]),
                    pn=dict(cv_totals=0.1, dist_rms_diff=0.2))
    monkeypatch.setattr(E, 'step_readings', step_readings)
    s, w = L.beside_readings(['B_inhadd1000', 'A1'])
    assert seen == dict(tag='inhadd1000', root=tmp_path / 'arch' / 'extended' / 'B_inhadd1000', wirings=['P3', 'P4'])
    b = s.set_index('variant')
    assert b.loc['B_inhadd1000', 'n_wirings'] == 2 and b.loc['B_inhadd1000', 'v2_held'] == 'inside'
    assert b.loc['B_inhadd1000', 'low_on'] == pytest.approx(0.005) and b.loc['B_inhadd1000', 'x2_r1_APL'] == 'x'
    assert b.loc['B_inhadd1000', 'low_class_on'] == E.low_conc_class(0.005, 0.05, 0.0136)
    assert b.loc['A1', 'n_wirings'] == 0 and {'family', 'F_inhadd1000 (E102)'} <= set(b.index)
    assert set(w.table) == {'e98', 'low'}


def test_part_fit(monkeypatch, tmp_path):
    """Step 3 refuses a rung the score has not named unless forced; it runs E136's fit, single-factor returns and
    rescore of the arm with the normalization fixed, in E139's folders, then the 28 reads through read_at at the
    classed point and the 8 beside reads; E136's arms, paths and measure are as before after it."""
    import json as _json
    res = tmp_path / 'res'
    res.mkdir()
    monkeypatch.setattr(L, 'out_dir', lambda: res)
    monkeypatch.setattr(L, 'local_dir', lambda part='': tmp_path / 'loc' / part)
    monkeypatch.setattr(L, 'archive_dir', lambda: tmp_path / 'arch')
    monkeypatch.setattr(L, 'read_dirs', lambda: [tmp_path / 'loc' / 'read', tmp_path / 'arch' / 'read'])
    with pytest.raises(RuntimeError, match='not a step-3 candidate'):
        L.part_fit(314.0)
    (res / 'classes.json').write_text(_json.dumps(dict(step3_candidates=[dict(variant='F_inhadd314', rung=314.0)])))
    assert L.candidate_rungs() == [314.0]
    calls, reads, besides = [], [], []
    z = dict(aplpn=-0.5, kckc=0.2)

    def fit_arm(name, maxfev=None, **kw):
        assert F.out_dir() == res and F.local_dir('fit_x') == tmp_path / 'loc' / 'fit_x'
        assert F.archive_dir() == tmp_path / 'arch' and F.MEASURE == L.MEASURE
        assert F.arm_kws(name)['pn_inhadd'] == pytest.approx(L.fit_rung(name)) and maxfev == F.MAXFEV['A0']
        calls.append(('fit', name))
        (res / f'fit_{name}.json').write_text(_json.dumps(dict(best_eval=3, best_z=z, n_evals=40)))
        return {}

    def determined(name):
        calls.append(('determined', name))
        fit = _json.loads((res / f'fit_{name}.json').read_text())
        fit.update(determined={}, classed_z=z)
        (res / f'fit_{name}.json').write_text(_json.dumps(fit))
        return fit
    monkeypatch.setattr(F, 'fit_arm', fit_arm)
    monkeypatch.setattr(F, 'determined', determined)
    monkeypatch.setattr(F, 'history_evals', lambda name: [0, 1, 2, 3])
    monkeypatch.setattr(F, 'rescore', lambda name: calls.append(('rescore', name)))

    def read_at(tag, arm, zz, wiring, panel, half, lane, root, dirs, reclaim=False, done=None, extra=None, **kw):
        reads.append(dict(tag=tag, arm=arm, z=zz(), root=root, dirs=dirs, extra=extra, half=half))
    monkeypatch.setattr(F, 'read_at', read_at)
    monkeypatch.setattr(L, 'read_beside', lambda name, w, lane, reclaim, done, arm=None, z=None:
                        besides.append((name, w, arm, dict(z))))
    monkeypatch.setattr(L, 'reads_table', lambda: pd.DataFrame(columns=['tag']))
    out = L.part_fit(314.0, lane='hal-fit')
    assert calls == [('fit', 'A0_inhadd314'), ('determined', 'A0_inhadd314'), ('rescore', 'A0_inhadd314')]
    assert len(reads) == 28 and all(r['arm'] == 'A0_inhadd314' and r['z'] == z and r['root'] == res for r in reads)
    assert reads[0]['tag'] == 'A0_inhadd314_megamat_P3' and reads[0]['half'] == 'held'
    assert reads[0]['extra'] == dict(variant='A0_inhadd314', setting='A0', rung=314.0, step='3')
    assert len(besides) == 8 and besides[0] == ('A0_inhadd314', F.WIRINGS[0], 'A0_inhadd314', z)
    assert out['classed_z'] == z and out['candidate'] and not out['forced']
    assert 'A0_inhadd314' not in F.ARMS and F.PATHS == {} and F.MEASURE == 'all'
    # a second call skips the done parts; a forced rung that is not a candidate runs
    calls.clear()
    pd.DataFrame(dict(eval=[0, 1, 2, 3])).to_csv(res / 'rescore_A0_inhadd314.csv', index=False)
    L.part_fit(314.0, parts=('fit', 'determined', 'rescore'))
    assert calls == []
    L.part_fit(100.0, force=True, parts=('fit',))
    assert calls == [('fit', 'A0_inhadd100')] and L.fitted_names() == ['A0_inhadd100', 'A0_inhadd314']


def test_load_fit_run(monkeypatch, tmp_path):
    import json as _json
    res = tmp_path / 'res'
    res.mkdir()
    monkeypatch.setattr(L, 'out_dir', lambda: res)
    monkeypatch.setattr(L, 'local_dir', lambda part='': tmp_path / 'loc' / part)
    monkeypatch.setattr(L, 'archive_dir', lambda: tmp_path / 'arch')
    monkeypatch.setattr(L, 'read_dirs', lambda: ['reads'])
    (res / 'fit_A0_inhadd314.json').write_text(_json.dumps(dict(best_eval=7)))
    seen = []
    monkeypatch.setattr(R, 'load_step3_run', lambda tag, dirs=None: seen.append((tag, dirs)) or tag)
    assert L.load_fit_run(F.ReadPoint('A0_inhadd314', 'megamat', 'P1', None)) == 'A0_inhadd314_e7_P1'
    assert seen[-1][1] == [tmp_path / 'loc' / 'fit_A0_inhadd314', tmp_path / 'arch' / 'fit_A0_inhadd314']
    assert L.load_fit_run(F.ReadPoint('A0_inhadd314', 'megamat', 'P3', 'held')) == 'A0_inhadd314_megamat_P3'
    assert L.load_fit_run(F.ReadPoint('A0_inhadd314', 'apl16', 'real', None)) == 'A0_inhadd314_apl16_real'
    assert seen[-1][1] == ['reads']
