#!/usr/bin/env python3
"""Tests for the three-calyx loader of E100's part D in E79_connectome_maps.py (the hemibrain's right calyx and the
male CNS's left and right, from Pratyush's tables, with no ROI rule): on made-up tables, the one rule (only the
counted claws, the boutons that carry them, their release-site weights and centroids, APL's synapses with those claws
and boutons, the sides), the arguments passed to E97's reader, the outline and posterior edge, the official-ROI
outline beside it, the neighbours' axis check passing and failing, and the guard against held paths; and one smoke
check on the real tables, skipped when they or E97's reader are not available."""

import inspect
import os
import re

import numpy as np
import pandas as pd
import pytest

E = pytest.importorskip('experiments.E79_connectome_maps')


# ---------------------------------------------------------------------------
# made-up tables in the forms the loader reads

def _vox(um):
    """um -> 8 nm voxels."""
    return np.asarray(um, float) * 1000.0 / E.VOXEL_NM


def made_up_tables():
    """Two brains. The hemibrain (one side): KC 1 with claws 0 (from PN 10, bouton 0, 6 synapses, counted) and 1 (from
    PN 10, bouton 1, 4 synapses, not counted); KC 2 with claw -1 (PN 11, bouton 0, counted). PN 10's bouton 0 has 3
    release sites, bouton 1 has 2, PN 11's bouton 0 one. APL synapses: two with KC 1 claw 0 (one each way), one with KC
    1 claw 1, one with KC 1 off the claws (NaN), one with KC 2 claw -1, one with PN 10 bouton 0, one with PN 10 bouton 1.
    The male CNS: KC 5 on the left (PN 20, bouton 0) and KC 6 on the right (PN 21, bouton 0), one APL synapse each."""
    hb_pn_kc = pd.DataFrame({
        'kc_id': [1] * 10 + [2] * 5, 'pn_id': [10] * 10 + [11] * 5,
        'claw_id': [0.0] * 6 + [1.0] * 4 + [-1.0] * 5, 'bouton_id': [0] * 6 + [1] * 4 + [0] * 5,
        'roi_pre': ['CA(R)'] * 5 + ['SLP(R)'] + ['CA(R)'] * 9, 'roi_post': ['CA(R)'] * 15, 'side': ['R'] * 15,
        'x_post': _vox([1.0] * 15), 'y_post': _vox([2.0] * 15), 'z_post': _vox([3.0] * 15)})
    hb_release = pd.DataFrame({'pn_id': [10, 10, 10, 10, 10, 11], 'bouton': [0, 0, 0, 1, 1, 0],
                               'x': _vox([0, 3, 6, 50, 50, 9]), 'y': _vox([0, 0, 0, 0, 0, 9]),
                               'z': _vox([0, 0, 0, 0, 0, 9])})

    def apl(rows, side):
        df = pd.DataFrame(rows, columns=['partner_id', 'kind', 'apl_input', 'claw_id', 'bouton', 'roi_pre'])
        n = len(df)
        df['apl_id'] = 99
        df['apl_unit'] = 0.0
        df['side'] = side
        df['roi_post'] = df['roi_pre']
        for i, c in enumerate(('x', 'y', 'z')):
            df[f'{c}_pre'] = _vox(np.arange(n) * 2.0 + i)
            df[f'{c}_post'] = _vox(np.arange(n) * 2.0 + i + 2.0)
        df['confidence_pre'] = df['confidence_post'] = 0.9
        return df
    hb_apl = apl([(1, 'KC', False, 0.0, np.nan, 'CA(R)'), (1, 'KC', True, 0.0, np.nan, 'CA(R)'),
                  (1, 'KC', False, 1.0, np.nan, 'CA(R)'), (1, 'KC', False, np.nan, np.nan, 'gL(R)'),
                  (2, 'KC', True, -1.0, np.nan, 'CA(R)'), (10, 'PN', False, np.nan, 0.0, 'SLP(R)'),
                  (10, 'PN', True, np.nan, 1.0, 'CA(R)')], 'R')
    mc_pn_kc = pd.DataFrame({'kc_id': [5] * 5 + [6] * 5, 'pn_id': [20] * 5 + [21] * 5, 'claw_id': [0.0] * 10,
                             'bouton_id': [0] * 10, 'roi_pre': ['CA(L)'] * 5 + ['CA(R)'] * 5,
                             'roi_post': ['CA(L)'] * 5 + ['CA(R)'] * 5, 'side': ['L'] * 5 + ['R'] * 5,
                             'x_post': _vox([300.0] * 5 + [100.0] * 5), 'y_post': _vox([10.0] * 10),
                             'z_post': _vox([10.0] * 10)})
    mc_release = pd.DataFrame({'pn_id': [20, 21], 'bouton': [0, 0], 'x': _vox([300, 100]), 'y': _vox([10, 10]),
                               'z': _vox([10, 10])})
    mc_apl = pd.concat([apl([(5, 'KC', False, 0.0, np.nan, 'CA(L)')], 'L'),
                        apl([(6, 'KC', True, 0.0, np.nan, 'CentralBrain-unspecified')], 'R')], ignore_index=True)
    return {'hemibrain': dict(pn_kc=hb_pn_kc, release=hb_release, apl=hb_apl),
            'malecns': dict(pn_kc=mc_pn_kc, release=mc_release, apl=mc_apl)}


def made_up_claws(with_e97_extras: bool = True):
    """The counted claws as E97's counted_claws gives them (KC 1's 4-synapse claw 1 is not counted)."""
    c = pd.DataFrame({'mb': ['hemibrain', 'hemibrain', 'malecns_L', 'malecns_R'], 'kc_id': [1, 2, 5, 6],
                      'claw_id': [0, -1, 0, 0], 'pn_id': [10, 11, 20, 21], 'glomerulus': ['DM1', 'DL1', 'VA7m', 'D'],
                      'n_synapses': [6, 5, 5, 5]})
    if with_e97_extras:
        c['bouton_ids'] = [(0,), (0,), (0,), (0,)]
        c['main_bouton_id'] = [0, 0, 0, 0]
        c['x_post'] = _vox([1.0, 1.0, 300.0, 100.0])
        c['y_post'] = _vox([2.0, 2.0, 10.0, 10.0])
        c['z_post'] = _vox([3.0, 3.0, 10.0, 10.0])
    return c


# ---------------------------------------------------------------------------
# the one rule

@pytest.mark.parametrize('extras', [True, False])
def test_the_one_rule_on_made_up_tables(extras):
    cal = E._assemble_calyces(made_up_claws(extras), made_up_tables())
    assert list(cal) == list(E.PRAT_CALYCES)
    hb = cal['hemibrain_R']
    # only the counted claws, in the (R, D, A) frame: hemibrain (-x, -z, +y)
    assert sorted(zip(hb['claws']['kc_id'], hb['claws']['claw_id'])) == [('1', 0), ('2', -1)]
    c1 = hb['claws'].set_index('kc_id').loc['1']
    assert (c1['R'], c1['D'], c1['A']) == pytest.approx((-1.0, -3.0, 2.0))
    assert c1['n_synapses'] == c1['n_PN_synapses'] == 6
    # one of its six synapses has a presynaptic end outside the ROI
    assert c1['n_synapses_in_roi'] == 5 and c1['roi_share'] == pytest.approx(5 / 6)
    # only the boutons that carry a counted claw: PN 10's bouton 1 carries only the uncounted claw
    b = hb['boutons'].set_index('bouton_id')
    assert sorted(b.index) == ['10_0', '11_0']
    # the weight is the number of release sites, the position their centroid
    assert b.loc['10_0', 'n_tbars'] == 3 and b.loc['11_0', 'n_tbars'] == 1
    assert b.loc['10_0', 'R'] == pytest.approx(-3.0) and b.loc['10_0', 'glomerulus'] == 'DM1'
    assert b.loc['10_0', 'n_counted_claws'] == 1
    # APL's synapses with those claws (both directions, claw -1 matched by (KC, -1)) and boutons, nothing else
    a = hb['apl']
    got = sorted(zip(a['kind'], a['partner_id'], a['apl_input'], a['claw_id'].astype('float').fillna(-9),
                     a['bouton'].astype('float').fillna(-9)))
    assert got == [('KC', '1', False, 0.0, -9.0), ('KC', '1', True, 0.0, -9.0), ('KC', '2', True, -1.0, -9.0),
                   ('PN', '10', False, -9.0, 0.0)]
    assert a.loc[a['kind'] == 'PN', 'in_roi'].tolist() == [False]
    # the male CNS: each claw in its KC's side's calyx, with its boutons and APL synapses
    for calyx, kc, in_roi in (('malecns_L', '5', True), ('malecns_R', '6', False)):
        p = cal[calyx]
        assert p['claws']['kc_id'].tolist() == [kc] and len(p['boutons']) == 1
        assert p['apl']['partner_id'].tolist() == [kc] and p['apl']['in_roi'].tolist() == [in_roi]
    # male CNS (-x, -y, -z): the left KC at x = 300 um lies at R = -300
    assert cal['malecns_L']['claws']['R'].iloc[0] == pytest.approx(-300.0)


def test_apl_position_is_the_midpoint():
    cal = E._assemble_calyces(made_up_claws(), made_up_tables())
    a = cal['malecns_L']['apl'].iloc[0]
    # pre (0, 1, 2) and post (2, 3, 4) um in x, y, z -> midpoint (1, 2, 3) -> (R, D, A) = (-1, -2, -3)
    assert (a['R'], a['D'], a['A']) == pytest.approx((-1.0, -2.0, -3.0))


def test_a_claw_counted_twice_is_refused():
    c = pd.concat([made_up_claws(), made_up_claws().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError):
        E._assemble_calyces(c, made_up_tables())


def test_claws_come_from_e97_with_no_roi_rule(monkeypatch):
    """With claws=None the loader asks E97's counted_claws for both brains with roi 'none' (and passes syn_min and
    min_confidence); a made-up reader applying the syn_min rule leaves the 4-synapse claw out."""
    calls = []

    class FakeE97:
        @staticmethod
        def counted_claws(brain, *, roi=None, syn_min=5, min_confidence=None):
            calls.append((brain, roi, syn_min, min_confidence))
            s = made_up_tables()[brain]['pn_kc']
            c = s.groupby(['kc_id', 'claw_id', 'pn_id', 'side']).size().rename('n_synapses').reset_index()
            c = c[c['n_synapses'] >= syn_min]
            c['mb'] = 'hemibrain' if brain == 'hemibrain' else 'malecns_' + c['side']
            c['glomerulus'] = 'DM1'
            return c

    monkeypatch.setattr(E, '_e97', lambda: FakeE97)
    monkeypatch.setattr(E, '_brain_tables', lambda b: made_up_tables()[b])
    cal = E.pratyush_calyces(syn_min=5, min_confidence=0.5)
    assert calls == [('hemibrain', 'none', 5, 0.5), ('malecns', 'none', 5, 0.5)]
    assert sorted(cal['hemibrain_R']['claws']['claw_id']) == [-1, 0]
    # at syn_min 4 the 4-synapse claw counts, and its bouton and APL synapse come with it
    cal = E.pratyush_calyces(syn_min=4)
    assert sorted(cal['hemibrain_R']['claws']['claw_id']) == [-1, 0, 1]
    assert '10_1' in set(cal['hemibrain_R']['boutons']['bouton_id'])
    assert (cal['hemibrain_R']['apl']['claw_id'] == 1).sum() == 1


def test_summary_counts():
    cal = E._assemble_calyces(made_up_claws(), made_up_tables())
    s = E.pratyush_calyx_summary(cal).set_index('calyx')
    assert s.loc['hemibrain_R', 'n_claws'] == 2 and s.loc['hemibrain_R', 'n_boutons'] == 2
    assert s.loc['hemibrain_R', 'n_tbars'] == 4 and s.loc['hemibrain_R', 'n_apl'] == 4
    assert s.loc['hemibrain_R', 'n_apl_to_kc'] == 1 and s.loc['hemibrain_R', 'n_kc_to_apl'] == 2
    assert s.loc['hemibrain_R', 'n_apl_to_pn'] == 1 and s.loc['hemibrain_R', 'n_pn_to_apl'] == 0
    assert s.loc['hemibrain_R', 'n_apl_on_claw_minus1'] == 1
    assert s.loc['hemibrain_R', 'share_claws_in_roi'] == 1.0


# ---------------------------------------------------------------------------
# the outline and the posterior edge

def _cylinder(n=20000, radius=20.0, a0=10.0, a1=40.0, seed=0):
    rng = np.random.default_rng(seed)
    r = radius * np.sqrt(rng.random(n))
    t = 2 * np.pi * rng.random(n)
    return pd.DataFrame({'R': r * np.cos(t), 'D': r * np.sin(t), 'A': rng.uniform(a0, a1, n)})


def test_calyx_outlines_on_made_up_points():
    apl = _cylinder()
    o = E.calyx_outlines(apl)
    # the posterior edge is the 0.5th percentile of A
    assert o['posterior_edge_A'] == pytest.approx(10.0 + 0.005 * 30.0, abs=0.1)
    # the whole outline is about the disk's area (pi r^2, the 97.5 % rule shaves a little)
    assert 0.85 * np.pi * 400 < o['whole'].sum() < 1.1 * np.pi * 400
    f = o['frame']
    iy, ix = int(round(f['D1'] / E.PX_UM)), int(round(-f['R0'] / E.PX_UM))
    assert o['whole'][iy, ix] and not o['whole'][0, 0]
    # slabs every 3 um from the edge, depths measured from it
    depths = sorted(o['slabs'])
    assert depths[0] == pytest.approx(E.SLAB_THICK_UM / 2) and np.allclose(np.diff(depths), E.SLAB_STEP_UM)
    # a slab too sparse for an outline is left out
    sparse = pd.concat([apl, pd.DataFrame({'R': [0.0] * 5, 'D': [0.0] * 5, 'A': [80.0] * 5})], ignore_index=True)
    assert max(E.calyx_outlines(sparse)['slabs']) < 40.0


def test_official_roi_outline_on_made_up_points():
    apl = _cylinder()
    # the official ROI misses the posterior 6 um
    apl['in_roi'] = apl['A'] >= 16.0
    claws = pd.DataFrame({'n_synapses_in_roi': [10, 2, 7, 0], 'A': [30.0, 12.0, 20.0, 11.0], 'roi_share': np.nan})
    cal = {'malecns_L': dict(apl=apl, claws=claws, boutons=None)}
    s, outs = E.official_roi_outline(cal)
    r = s.iloc[0]
    assert r['share_apl_in_roi'] == pytest.approx(0.8, abs=0.02)
    assert r['share_claws_in_roi'] == 0.5
    assert r['roi_minus_all_edge_um'] == pytest.approx(6.0, abs=0.3)
    assert r['apl_A_outside_minus_inside_um'] < -10
    assert 0.95 < r['whole_iou'] <= 1.0
    assert outs['malecns_L']['posterior_edge_A'] == pytest.approx(r['posterior_edge_roi_um'])
    assert r['n_slabs_roi'] < r['n_slabs_all']


# ---------------------------------------------------------------------------
# the neighbours' axis check

def _neighbour_setup(flip_dorsal: bool = False, n_ped: int = 30):
    """A right calyx (centre at R = 50) and a left one (R = -50) in one brain, the hemibrain's right one beside, and
    the neighbours where they lie in a fly: PED ventral and anterior, SLP dorsal and anterior, PLP ventral and lateral."""
    def calyx(r):
        return {'apl': pd.DataFrame({'R': [r - 1, r + 1], 'D': [0.0, 0.0], 'A': [0.0, 0.0]})}
    calyces = {'hemibrain_R': calyx(50.0), 'malecns_L': calyx(-50.0), 'malecns_R': calyx(50.0)}
    sd = -1.0 if flip_dorsal else 1.0
    syn = {}
    for brain, sides in (('hemibrain', ('R',)), ('malecns', ('L', 'R'))):
        rows = []
        for side in sides:
            lat = 1.0 if side == 'R' else -1.0
            c = 50.0 * lat
            for nb, (dr, dd, da), n in (('PED', (5, -20, 10), n_ped), ('SLP', (10, 10, 10), 40),
                                        ('PLP', (20, -20, 0), 40)):
                rows += [dict(roi=f'{nb}({side})', R=c + lat * dr, D=sd * dd, A=da)] * n
        syn[brain] = pd.DataFrame(rows)
    return calyces, syn


def test_axis_check_passes():
    calyces, syn = _neighbour_setup()
    t = E._neighbour_rows(calyces, syn)
    assert set(t['status']) == {'pass', 'n/a'} and t['calyx_passed'].all()
    assert t.loc[(t.calyx == 'hemibrain_R') & (t.check == 'right mushroom body at larger R'), 'status'].item() == 'n/a'
    row = t[(t.calyx == 'malecns_L') & (t.check == 'PLP')].iloc[0]
    assert row['d_lateral_um'] == pytest.approx(20.0) and row['d_R_um'] == pytest.approx(-20.0)
    E._raise_if_failed(t, strict=True)


def test_axis_check_fails_clearly(monkeypatch):
    calyces, syn = _neighbour_setup(flip_dorsal=True)
    t = E._neighbour_rows(calyces, syn)
    failed = t[t.status == 'fail']
    assert set(failed['check']) == {'PED', 'SLP', 'PLP'} and not t['calyx_passed'].any()
    with pytest.raises(E.AxisCheckError, match='FAILED') as ex:
        E._raise_if_failed(t, strict=True)
    assert ex.value.table is t
    with pytest.warns(UserWarning, match='part D stops'):
        E._raise_if_failed(t, strict=False)
    # through the public function
    monkeypatch.setattr(E, '_neighbour_synapses', lambda b: syn[b])
    with pytest.raises(E.AxisCheckError):
        E.neighbour_axis_checks(calyces)
    # left and right swapped: the right calyx at smaller R fails the two-sides check
    calyces, syn = _neighbour_setup()
    swapped = {'malecns_L': calyces['malecns_R'], 'malecns_R': calyces['malecns_L']}
    t = E._neighbour_rows(swapped, syn)
    assert (t.loc[t.check == 'right mushroom body at larger R', 'status'] == 'fail').all()


def test_axis_check_fails_on_too_few_synapses():
    calyces, syn = _neighbour_setup(n_ped=5)
    t = E._neighbour_rows(calyces, syn, min_synapses=20)
    ped = t[t.check == 'PED']
    assert (ped['status'] == 'fail').all() and ped['note'].str.contains('fewer than').all()


# ---------------------------------------------------------------------------
# the held paths

NEW_FUNCTIONS = ('_e97', '_prat_dirs', '_prat_files', '_held_path', '_read_prat', '_to_rda_voxels', '_side_letter',
                 '_release_sites', '_harmonize_apl', '_brain_tables', '_calyx_roi', '_assemble_calyces',
                 'pratyush_calyces', '_claws_in_roi', 'pratyush_calyx_summary', 'calyx_outlines',
                 '_neighbour_synapses', '_expect', '_neighbour_rows', 'neighbour_axis_checks', '_raise_if_failed',
                 'official_roi_outline', 'prat_inputs', 'write_provenance_prat', '_with_caveats', 'run_pratyush')


def test_new_code_paths_name_no_held_data():
    """The new functions' code names neither the held tables (Yang's own Claude's project) nor FlyWire, FAFB or BANC,
    nor the readers and paths of the held tables (CLAUDE.md, "Connectome data"), and the files they open are not
    held."""
    banned = [re.escape(m.lower()) for m in E.HELD_PATH_MARKERS]
    banned += [r'\byang_apl\b', r'(?<![\w.])_apl_table\(', r'(?<![\w.])load_calyces\(', r'(?<![\w.])axis_checks\(',
               r'(?<![\w.])_prat_boutons\(', r'\bprat_dir\b', r'\bprat_pn_']
    for name in NEW_FUNCTIONS:
        src = inspect.getsource(getattr(E, name)).lower()
        code = '\n'.join(ln for ln in src.splitlines() if not ln.strip().startswith('#'))
        for pat in banned:
            assert not re.search(pat, code), (name, pat)
    for f in list(E.HEMIBRAIN_FILES.values()) + list(E.MALECNS_FILES.values()):
        assert not E._held_path(f), f


def test_read_prat_refuses_held_and_outside_paths(tmp_path):
    root = tmp_path / 'from_pratyush' / '2025-12-05'
    root.mkdir(parents=True)
    pd.DataFrame({'a': [1, 2]}).to_parquet(root / 'ok.parquet')
    pd.DataFrame({'a': [3]}).to_csv(root / 'ok.csv', index=False)
    assert E._read_prat(root / 'ok.parquet', roots=[root])['a'].tolist() == [1, 2]
    assert E._read_prat(root / 'ok.csv', roots=[root])['a'].tolist() == [3]
    # outside the roots, also by '..'
    other = tmp_path / 'from_pratyush' / '2024-09-13'
    other.mkdir()
    pd.DataFrame({'a': [1]}).to_parquet(other / 'x.parquet')
    for p in (other / 'x.parquet', root / '..' / '2024-09-13' / 'x.parquet', tmp_path / 'y.parquet'):
        with pytest.raises(PermissionError):
            E._read_prat(p, roots=[root])
    # a held marker inside an allowed root
    held = root / 'claude_APL_modeling_YZ4TC'
    held.mkdir()
    pd.DataFrame({'a': [1]}).to_parquet(held / 'x.parquet')
    with pytest.raises(PermissionError):
        E._read_prat(held / 'x.parquet', roots=[root])
    with pytest.raises(PermissionError):
        E._read_prat(root / 'flywire_x.parquet', roots=[root])


# ---------------------------------------------------------------------------
# the real tables (skipped without them or without E97's reader)

def _real_e97():
    try:
        E97 = E._e97()
        dirs = E._prat_dirs()
    except Exception as ex:  # noqa: BLE001
        pytest.skip(f'E97 not importable: {ex!r}')
    if not all(Path_exists(d) for d in dirs.values()):
        pytest.skip('Pratyush\'s tables not on this host')
    try:
        E97.counted_claws('hemibrain', roi='none')
    except NotImplementedError:
        pytest.skip("E97's counted_claws is not implemented yet")
    return E97


def Path_exists(d) -> bool:
    return os.path.isdir(str(d))


def test_real_three_calyces_smoke():
    """The three calyces' counts beside the planner's probes (notes/E97_E100_replan_notes_2026-11-17.md), printed for
    the record; the axis check must pass."""
    _real_e97()
    cal = E.pratyush_calyces()
    s = E.pratyush_calyx_summary(cal)
    pd.set_option('display.width', 200)
    print(s.T.to_string())
    gl = ['VM5d', 'DM4', 'VM2', 'DC3', 'VM5v', 'VM7d']
    print(pd.DataFrame({c: p['boutons'].groupby('glomerulus').size().reindex(gl) for c, p in cal.items()}).to_string())
    roi, _ = E.official_roi_outline(cal)
    print(roi.T.to_string())
    ax = E.neighbour_axis_checks(cal, strict=True)
    print(ax.to_string())
    assert list(cal) == list(E.PRAT_CALYCES)
    for c, p in cal.items():
        assert len(p['claws']) > 5000 and len(p['boutons']) > 200 and len(p['apl']) > 10000, c
        assert (p['boutons']['n_tbars'] > 0).all() and p['boutons']['glomerulus'].notna().all()
    assert ax['calyx_passed'].all()


# ---------------------------------------------------------------------------
# the readers of the held tables are switched off (CLAUDE.md, "Connectome data (Tom, 2026-11-17)")

def test_held_readers_refuse_naming_the_rule(monkeypatch, tmp_path):
    """load_calyces, _apl_table, axis_checks and write_provenance raise HeldTablesError (a PermissionError) naming the
    rule, before opening anything; main() refuses without --calyces pratyush; run_pratyush never writes to OUT."""
    import builtins
    opened = []
    real_open = builtins.open

    def spy(path, *a, **k):
        opened.append(str(path))
        return real_open(path, *a, **k)
    monkeypatch.setattr(builtins, 'open', spy)
    monkeypatch.setattr(pd, 'read_csv', lambda p, *a, **k: opened.append(str(p)) or pd.DataFrame())
    for call in (lambda: E.load_calyces(), lambda: E._apl_table('hemibrain'), lambda: E.axis_checks(),
                 lambda: E.write_provenance({})):
        with pytest.raises(E.HeldTablesError, match='Connectome data'):
            call()
    assert not any('claude_apl_modeling' in p.lower() for p in opened), opened
    assert issubclass(E.HeldTablesError, PermissionError)
    monkeypatch.setattr('sys.argv', ['E79_connectome_maps.py', '--analyses', '4'])
    with pytest.raises(E.HeldTablesError):
        E.main()
    with pytest.raises(ValueError, match='withheld'):
        E.run_pratyush(E.OUT)


def test_prat_inputs_are_the_allowed_tables():
    try:
        got = E.prat_inputs()
    except Exception as ex:  # noqa: BLE001
        pytest.skip(f"E97's folders not available: {ex!r}")
    assert {k.split('/')[0] for k in got} == {'hemibrain', 'malecns'}
    assert {k.split('/')[1] for k in got} == set(E.HEMIBRAIN_FILES) | {'pn_kc'}
    for v in got.values():
        assert not E._held_path(v)
        assert str(v.parent).endswith(('from_pratyush/2025-12-05', 'data/Male-CNS_Tom')), v
