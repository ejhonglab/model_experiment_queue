"""Tests of reports/yang_model_scheme.py on small made-up tables and source snippets (no shared directory needed)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reports import yang_model_scheme as S


def write_tables(shared: Path, *, a1_boutons=5, a1_claws=15, stats_kcs=6, target_rate=0.05, real_row=True):
    """The four tables of S.TABLES under `shared`: two glomeruli of the wiring (one with ORN input), two left out."""
    for p in S.TABLES.values():
        (shared / p).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dict(glomerulus=['A', 'B', 'C', 'D'], model_n_pn=[2.0, 1.0, np.nan, 0.0],
        model_n_boutons=[3.0, 2.0, np.nan, 0.0], model_n_claws=[10.0, 5.0, np.nan, 0.0],
        driven=[True, False, True, True])).to_csv(shared / S.TABLES['glomerulus_substrate'], index=False)
    pd.DataFrame(dict(wiring=['real', 'P1'] if real_row else ['P1', 'P2'], n_boutons=[a1_boutons, 5],
        n_claws=[a1_claws, 15], n_kcs=[6, 6])).to_csv(shared / S.TABLES['a1_path'], index=False)
    pd.DataFrame(dict(tag=['bouton_base', 'bouton_coincidence-2'], rate_on=[0.0516, 0.03], rate_off=[0.0978, 0.09],
        target_rate=[target_rate, 0.05], n_kcs=[stats_kcs, 6])).to_csv(shared / S.TABLES['stats_variants'],
        index=False)
    rows = [dict(tag='bouton_base', condition=c, odor_a=a, odor_b=b, distance=0.5)
        for c in ('on', 'off') for a, b in (('x', 'y'), ('x', 'w'), ('y', 'w'))]
    rows.append(dict(tag='bouton_coincidence-2', condition='on', odor_a='x', odor_b='z', distance=0.4))
    pd.DataFrame(rows).to_csv(shared / S.TABLES['pair_distances'], index=False)


def test_load_counts(tmp_path):
    write_tables(tmp_path)
    c = S.load_counts(tmp_path)
    assert (c['n_glomeruli'], c['n_driven'], c['n_undriven'], c['n_pn']) == (2, 1, 1, 3)
    assert (c['n_boutons'], c['n_claws'], c['n_kcs'], c['n_odors']) == (5, 15, 6, 3)
    assert c['claws_per_kc'] == 2.5 and c['claws_per_bouton'] == 3
    assert (c['rate_on'], c['rate_off'], c['target_rate']) == (0.0516, 0.0978, 0.05)


@pytest.mark.parametrize('kw, match', [(dict(a1_boutons=6), 'disagree'), (dict(a1_claws=16), 'disagree'),
    (dict(stats_kcs=7), 'KCs'), (dict(real_row=False), "wiring == 'real'")])
def test_load_counts_fails_loudly(tmp_path, kw, match):
    write_tables(tmp_path, **kw)
    with pytest.raises(ValueError, match=match):
        S.load_counts(tmp_path)


CPP = """ModelParams const DEFAULT_PARAMS = []() {
    p.pn.taum       = 0.01;
    p.kc.apl_taum              = 0.05;
    p.kc.taum                  = 0.01;
    std::string msg = cat("p.kc.apl_taum=", p.kc.apl_taum);
"""
E4 = """CLAWS = dict(one_row_per_claw=True, prat_claws=True,
    n_spikes_for_response=2
)
VARIANTS = {'1spk': (dict(CLAWS, n_spikes_for_response=None), False)}
"""


def snippet_roots(tmp_path, cpp=CPP):
    olf, queue = tmp_path / 'olfsysm', tmp_path / 'queue'
    (olf / S.OLFSYSM_CPP).parent.mkdir(parents=True)
    (olf / S.OLFSYSM_CPP).write_text(cpp)
    (queue / 'experiments').mkdir(parents=True)
    (queue / 'experiments' / 'E4_substrate.py').write_text(E4)
    return {'olfsysm': olf, 'queue': queue}


def snippet_constants():
    claws = next(c for c in S.CONSTANTS if c.name == 'n_spikes_for_response')
    return [S._default('p.pn.taum', 0.01, 2), S._default('p.kc.apl_taum', 0.05, 3), S._default('p.kc.taum', 0.01, 4),
        claws]


def test_constants_match_the_source(tmp_path):
    out = S.check_constants(snippet_constants(), snippet_roots(tmp_path))
    assert [r['status'] for r in out] == ['checked'] * 4
    # the field in a string (line 5) and the variant overriding the CLAWS dict (line 4 of E4) are not matched
    assert [r['line'] for r in out] == [2, 3, 4, 2]


def test_constants_mismatch_fails(tmp_path):
    with pytest.raises(ValueError, match='p.pn.taum is 0.02'):
        S.check_constants(snippet_constants(), snippet_roots(tmp_path, CPP.replace('0.01;', '0.02;', 1)))
    roots = snippet_roots(tmp_path / 'b', CPP + '    p.pn.taum = 0.01;\n')
    with pytest.raises(ValueError, match='2 matches'):
        S.check_constants(snippet_constants(), roots)


def test_constants_without_the_source(tmp_path):
    roots = snippet_roots(tmp_path)
    for olf in (None, tmp_path / 'nowhere'):
        out = S.check_constants(snippet_constants(), dict(roots, olfsysm=olf))
        assert [r['status'] for r in out] == ['source not found'] * 3 + ['checked']
        assert all(r['line'] is None for r in out[:3])


# sim_PN_layer's lines, and the scalar copy of the PN dynamics olfsysm 7528f6d added (sim_PN_glom_shifted), which repeats
# its literals
PN_CPP = """\
void sim_PN_layer() {
        dPNdt +=
            200.0*((orn_delta.array()+p.pn.offset)*p.pn.tanhsc/200.0*inh_PN).matrix(
                ).unaryExpr<double(*)(double)>(&tanh);
        inh_PN = p.pn.inhsc/(p.pn.inhadd+0.25*inhA(t)+0.75*inhB(t));
        pn_t.col(t) = pn_t.col(t-1) + dPNdt*p.time.dt/p.pn.taum;
}
void sim_PN_glom_shifted() {
        dPNdt += 200.0*std::tanh((orn_delta+p.pn.offset)*p.pn.tanhsc/200.0*inh_PN);
        inh_PN = p.pn.inhsc/(p.pn.inhadd+0.25*inhA(t)+0.75*inhB(t));
        double v = out(t-1) + dPNdt*p.time.dt/p.pn.taum;
}
"""


def test_pn_literals_are_sim_pn_layers(tmp_path):
    """The PN literals are matched once each, in sim_PN_layer's lines, not in the scalar copy."""
    pn = [c for c in S.CONSTANTS if c.name.startswith('PN ')]
    assert len(pn) == 4
    out = S.check_constants(pn, snippet_roots(tmp_path, PN_CPP))
    assert [r['status'] for r in out] == ['checked'] * 4
    assert [r['line'] for r in out] == [3, 3, 5, 5]


def test_table_matches_the_installed_simulator():
    """CONSTANTS against the olfsysm checkout and this repository's experiment scripts (skipped without olfsysm)."""
    olf = S.olfsysm_source_dir()
    if olf is None:
        pytest.skip('no olfsysm source on this machine')
    out = S.check_constants(S.CONSTANTS, {'olfsysm': olf, 'queue': S.REPO})
    assert all(r['status'] == 'checked' for r in out)
    assert len({c.name for c in S.CONSTANTS}) == len(S.CONSTANTS)


def made_up_counts(**kw):
    c = dict(n_glomeruli=52, n_driven=35, n_undriven=17, n_pn=108, n_boutons=401, n_claws=9512, n_kcs=1690,
        n_odors=16, claws_per_kc=5.6, claws_per_bouton=24, rate_on=0.0507, rate_off=0.1012, target_rate=0.05)
    c.update(kw)
    return c


def test_figure_values():
    v = S.figure_values(made_up_counts(), S.constant_values())
    assert (v['n_claws'], v['n_kcs'], v['claws_per_kc'], v['claws_per_bouton']) == ('9,512', '1,690', '5.6', '24')
    assert (v['rate_on'], v['rate_off']) == ('0.051', '0.101')
    assert (v['tau_orn'], v['tau_apl'], v['odor_s'], v['window_from'], v['window_to']) == ('10', '50', '0.5',
        '−0.5', '0.75')
    assert (v['pct_pre_apl'], v['pct_target'], v['n_spikes'], v['thr_factor']) == ('10', '5', '2', '2')
    assert (v['tanhsc'], v['offset'], v['inhsc'], v['inhadd']) == ('5.34', '2.94', '368.7', '31.4')
    assert (v['w_a'], v['w_b']) == (r'\dfrac{1}{4}', r'\dfrac{3}{4}')


def test_draw_writes_a_png(tmp_path):
    out = tmp_path / 'scheme.png'
    size = S.draw(S.figure_values(made_up_counts(), S.constant_values()), out)
    assert size == (round(S.FIG_W * S.DPI), round(S.FIG_H * S.DPI))
    data = out.read_bytes()
    assert data[:8] == b'\x89PNG\r\n\x1a\n' and len(data) > 100_000


def test_layout_check_sees_overlaps():
    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    with matplotlib.rc_context(S.RC):
        fig = Figure(figsize=(4, 3), dpi=100)
        FigureCanvasAgg(fig)
        cv = S._Canvas(fig)
        cv.box(0.5, 0.5, 2.0, 1.5, fill='w', edge='k')
        cv.text(0.6, 1.0, 'inside', 10, va='center')
        assert S._layout_problems(cv) == []
        cv.text(1.8, 1.0, 'across the edge', 10, va='center')
        cv.text(0.7, 1.0, 'over another', 10, va='center')
        cv.arrow((2.5, 2.0), (3.5, 2.0))
        cv.text(2.7, 2.0, 'on an arrow', 10, va='center')
        found = '\n'.join(S._layout_problems(cv))
    assert 'box edge' in found and 'texts overlap' in found and "arrow: 'on an arrow'" in found


def test_tokens_keep_math_whole():
    assert S._tokens(r'$\theta_k = \theta + 2\times$ (the KC) $a$:') == [r'$\theta_k = \theta + 2\times$', '(the',
        'KC)', '$a$:']


def test_run_without_the_simulator_source(tmp_path, capsys):
    shared, out = tmp_path / 'shared', tmp_path / 'out'
    write_tables(shared)
    prov = S.run(shared, out, roots={'olfsysm': None, 'queue': S.REPO})
    assert (out / S.PNG_NAME).stat().st_size > 0
    saved = json.loads((out / S.PROV_NAME).read_text())
    assert saved['constants_checked'] is False and saved['olfsysm_source'] is None
    unchecked = {c.name for c in S.CONSTANTS if c.root == 'olfsysm'}
    assert set(saved['constants_unchecked']) == unchecked
    assert {r['name'] for r in saved['constants'] if r['status'] == 'checked'} == {'TARGET', 'n_spikes_for_response'}
    assert set(saved['inputs']) == set(S.TABLES) and all(len(i['sha256']) == 64 for i in saved['inputs'].values())
    assert saved['counts']['n_claws'] == 15 and saved['values_on_figure']['n_claws'] == '15'
    assert saved['queue_commit'] and prov['figure']['sha256'] == saved['figure']['sha256']
    assert 'not checked' in capsys.readouterr().out


def test_run_refuses_a_run_tuned_elsewhere(tmp_path):
    write_tables(tmp_path, target_rate=0.068)
    with pytest.raises(ValueError, match='tuned to 0.068'):
        S.run(tmp_path, tmp_path / 'out', roots={'olfsysm': None, 'queue': S.REPO})
