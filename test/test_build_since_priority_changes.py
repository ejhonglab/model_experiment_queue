"""Tests of reports/build_since_priority_changes.py on a scratch template, with a scratch shared directory (its
data-target registry and experiment index) and every output in tmp_path, never in the shared directory or in
outputs/_artifacts."""
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / 'reports' / 'build_since_priority_changes.py'

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


def build(tmp_path, template):
    shared = scratch_shared(tmp_path)
    t = tmp_path / 'template.md'
    t.write_text(template)
    out_md, out_html = tmp_path / 'md' / 'page.md', tmp_path / 'html' / 'index.html'
    r = subprocess.run([sys.executable, str(SCRIPT), '--template', str(t), '--out-md', str(out_md),
        '--out-html', str(out_html)], env=dict(os.environ, KC_REMAPPING_SHARED=str(shared)), cwd=str(REPO),
        capture_output=True, text=True)
    return r, shared, t, out_md, out_html


def test_build(tmp_path):
    r, shared, t, out_md, out_html = build(tmp_path, TEMPLATE)
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
    assert (out_html.parent / 'since_priority_changes_filled.md').read_text().count('⟦pv:') == 7
    # nothing written to the (scratch) shared directory
    assert not (shared / 'add_to_reports').exists()


def test_unresolved_placeholder_fails(tmp_path):
    r, shared, t, out_md, out_html = build(tmp_path, TEMPLATE + '\nAnd {{data:no_such_target ~ .2f}}.\n')
    assert r.returncode != 0
    assert 'unresolved placeholders' in r.stderr and 'data:no_such_target' in r.stderr
    assert not out_md.exists() and not out_html.exists() and not out_html.parent.exists()
