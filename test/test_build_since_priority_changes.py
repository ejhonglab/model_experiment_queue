"""Tests of the digest pages' scripts (reports/build_since_priority_changes.py and build_since_yang_new_data.py, both
build_summary.page_main) on a scratch template, with a scratch shared directory (its data-target registry and
experiment index) and every output in tmp_path, never in the shared directory or in outputs/_artifacts."""
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
# each digest page's name: its script is reports/build_<name>.py
NAMES = ['since_priority_changes', 'since_yang_new_data']

TEMPLATE = """# Since the priority changes

last model experiment: {{meta:last_experiment}}, edited: {{meta:generated}}

[TOC]

%%% hover
APL effect: the change in odor-odor distance when APL output is blocked
%%%

## What we learned {#learned}

The APL effect in the KC data is {{data:apl_effect_mean ~ +.3f}} (90 % interval {{data:apl_effect_mean.lo ~ .2f}}
to {{data:apl_effect_mean.hi ~ .2f}}).

??? How it was measured
| statistic | APL effect{span=2} |
|---|---|---|
| KC data | {{data:apl_effect_mean ~ .2f}} | {{data:apl_effect_mean.hi ~ .2f}} |
???
"""


def scratch_shared(tmp_path):
    shared = tmp_path / 'shared'
    (shared / 'results' / 'data_targets').mkdir(parents=True)
    (shared / 'results' / 'data_targets' / 'data_targets.json').write_text(json.dumps(dict(apl_effect_mean=dict(
        value=0.12345, ci=[0.05, 0.2], definition='mean APL effect over odor pairs', function='apl_effect_mean',
        commit='abcdef0123'))))
    (shared / 'experiment_priorities.md').write_text(
        '| E1 | the baseline | completed | hal | report | report | 11-04 | 11-04 09:30 (#1) |\n')
    return shared


def build(tmp_path, template, name):
    shared = scratch_shared(tmp_path)
    t = tmp_path / 'template.md'
    if template is not None:
        t.write_text(template)
    out_md, out_html = tmp_path / 'md' / 'page.md', tmp_path / 'html' / 'index.html'
    script = REPO / 'reports' / f'build_{name}.py'
    r = subprocess.run([sys.executable, str(script), '--template', str(t), '--out-md', str(out_md),
        '--out-html', str(out_html)], env=dict(os.environ, KC_REMAPPING_SHARED=str(shared)), cwd=str(REPO),
        capture_output=True, text=True)
    return r, shared, t, out_md, out_html


@pytest.mark.parametrize('name', NAMES)
def test_build(tmp_path, name):
    r, shared, t, out_md, out_html = build(tmp_path, TEMPLATE, name)
    assert r.returncode == 0, r.stdout + r.stderr
    md = out_md.read_text()
    assert '{{' not in md and '⟦pv:' not in md
    assert 'is +0.123 (90 % interval 0.05' in md and '| KC data | 0.12 | 0.20 |' in md
    assert f'last model experiment: {datetime.now().year}-11-04 09:30, edited: ' in md
    side = json.loads(out_md.with_suffix('.values.json').read_text())
    assert set(side) == {'template', 'template_commit', 'generated', 'n_values', 'lint', 'values'}
    assert side['template'] == str(t.resolve()) and side['n_values'] == 7 and side['lint'] == []
    assert [v['key'] for v in side['values']][:2] == ['meta:last_experiment', 'meta:generated']
    # the page: source tooltips, the template's markup converted
    page = out_html.read_text()
    assert '<title>Since the priority changes</title>' in page
    assert 'KC / ORN data: apl_effect_mean' in page and '<details' in page and 'id="learned"' in page
    assert (out_html.parent / f'{name}_filled.md').read_text().count('⟦pv:') == 7
    # nothing written to the (scratch) shared directory
    assert not (shared / 'add_to_reports').exists()


@pytest.mark.parametrize('name', NAMES)
def test_experiment_key_links(tmp_path, name):
    """The experiment key's IDs link to the experiments' pages saved in the shared directory's reports/ (relative
    links); an experiment without a page keeps a plain ID, and without reports/ nothing is linked."""
    key = '\n??? Experiment key: the ID of each experiment\n- E1: the baseline (completed)\n- E10: the joint fit (running)\n???\n'
    shared = scratch_shared(tmp_path / 'first')
    (shared / 'reports').mkdir()
    (shared / 'reports' / 'E1_baseline.html').write_text('<body>E1</body>')
    t = tmp_path / 'template.md'
    t.write_text(TEMPLATE + key)
    out_html = tmp_path / 'html' / 'index.html'
    r = subprocess.run([sys.executable, str(REPO / 'reports' / f'build_{name}.py'), '--template', str(t), '--out-md',
        str(tmp_path / 'md' / 'page.md'), '--out-html', str(out_html)],
        env=dict(os.environ, KC_REMAPPING_SHARED=str(shared)), cwd=str(REPO), capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert 'experiment key: 1 IDs linked to their pages; no page for E10' in r.stdout
    page = out_html.read_text()
    assert '<li><a href="E1_baseline.html">E1</a>: the baseline (completed)</li>' in page
    assert '<li>E10: the joint fit (running)</li>' in page
    # no reports/ folder: the page is built, without links
    r, shared, t, out_md, out_html = build(tmp_path / 'second', TEMPLATE + key, name)
    assert r.returncode == 0, r.stdout + r.stderr
    assert '<li>E1: the baseline (completed)</li>' in out_html.read_text()


@pytest.mark.parametrize('name', NAMES)
def test_unresolved_placeholder_fails(tmp_path, name):
    r, shared, t, out_md, out_html = build(tmp_path, TEMPLATE + '\nAnd {{data:no_such_target ~ .2f}}.\n', name)
    assert r.returncode != 0
    assert 'unresolved placeholders' in r.stderr and 'data:no_such_target' in r.stderr
    assert not out_md.exists() and not out_html.exists() and not out_html.parent.exists()


@pytest.mark.parametrize('name', NAMES)
def test_missing_template_fails(tmp_path, name):
    r, shared, t, out_md, out_html = build(tmp_path, None, name)
    assert r.returncode != 0 and f'no template at {t}' in r.stderr
    assert not out_md.exists() and not out_html.parent.exists()


def test_default_paths(tmp_path):
    """Each script's defaults follow its name (the template beside it, the page under outputs/_artifacts/)."""
    shared = scratch_shared(tmp_path)
    for name in NAMES:
        r = subprocess.run([sys.executable, str(REPO / 'reports' / f'build_{name}.py'), '--help'], cwd=str(REPO),
            env=dict(os.environ, KC_REMAPPING_SHARED=str(shared)), capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        help_text = ' '.join(r.stdout.split())
        assert f'reports/{name}_template.md' in help_text and f'outputs/_artifacts/{name}/index.html' in help_text
        assert f'add_to_reports/{name}.md' in help_text
