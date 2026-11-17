"""Tests of reports/yang_short_report.py (a short hand-written report for Yang, filled and converted as her version of
the summary is) on small markdown files in tmp_path, with the placeholder resolver replaced by a fake one (as in
test_yang_summary.py) and the shared directory and outputs/_artifacts pointed at tmp_path: nothing here reads or writes
the shared directory or the repository's outputs."""
import json
import os
import re
import sys
import tempfile
from pathlib import Path

import pytest

from kcq import harness as H

try:
    H.shared_dir()
except FileNotFoundError:   # report_values looks the shared directory up on import; give it an empty one
    os.environ['KC_REMAPPING_SHARED'] = tempfile.mkdtemp(prefix='yang_short_report_test_shared_')

from reports import report_values as R   # noqa: E402
from reports import yang_short_report as S   # noqa: E402

DATA = {'a_dist': 0.123, 'b_dist': 0.5, 'rho': 0.65, 'count': 92}
TABLE = {'E1/t.csv:tag=x:col_a': 0.3}
REAL_RESOLVE = R.resolve
PNG = b'\x89PNG\r\n\x1a\n' + b'\x00' * 16   # embedded as it is (the converter does not decode it)


def fake_resolve(key):
    """data: and table keys from the dicts above; lit:, diff: and meta:generated as report_values resolves them
    (diff: calls back into this function through the patched module global)."""
    key = key.strip()
    if key.startswith('data:'):
        name = key[5:]
        if name not in DATA:
            raise KeyError(f'{key}: not in the data-target registry')
        return R.Resolved(key, DATA[name], dict(kind='data', name=name, part=None, definition=f'def of {name}',
            function='f', commit='abcdef0'))
    if key in TABLE:
        return R.Resolved(key, TABLE[key], dict(kind='table', file='results/E1/t.csv', row='tag=x', column='col_a',
            sha256='0', commit=None))
    if key.startswith(('lit:', 'diff:', 'meta:generated')):
        return REAL_RESOLVE(key)
    raise KeyError(f'{key}: no table row / column matches')


@pytest.fixture(autouse=True)
def world(monkeypatch, tmp_path):
    """The fake resolver, a shared directory with one figure and outputs/_artifacts, all in tmp_path."""
    monkeypatch.setattr(R, 'resolve', fake_resolve)
    shared = tmp_path / 'shared'
    (shared / 'results').mkdir(parents=True)
    (shared / 'results' / 'fig.png').write_bytes(PNG)
    monkeypatch.setattr(S, '_shared', lambda: shared)
    monkeypatch.setattr(S, 'ARTIFACTS', tmp_path / 'artifacts')
    return shared


SAMPLE = """# Inhibition in **short**

Built {{meta:generated}}.

[TOC]

## The data at {{data:rho ~ .2f}}

The mean distance is {{data:a_dist ~ .2f}}, a correlation of {{inv:data:a_dist ~ .2f}}; APL changes it by
{{neg:data:rho ~ +.2f}}, and a hand-entered {{lit:0.05 # Yang, by hand}}.

| statistic | control flies |
|---|---|
| mean pair correlation | {{inv:data:b_dist ~ .2f}} |
| a table value | {{E1/t.csv:tag=x:col_a ~ .1f}} |

![The figure at {{data:rho ~ .2f}}](results/fig.png)

??? More detail

- the difference {{diff:data:a_dist -- data:b_dist ~ +.3f}}

???

%%% hover
TNT: the toxin, in {{data:count ~ d}} flies
%%%

TNT blocks APL output.
"""


def write_md(tmp_path, text, name='short.md'):
    p = tmp_path / name
    p.write_text(text)
    return p


def out_dir(name):
    return S.ARTIFACTS / name


def pv(text, tip):
    """A filled value as the page shows it."""
    return f'<span class="pv" title="{tip}">{text}</span>'


def test_build_fills_and_converts(tmp_path, world):
    md = write_md(tmp_path, SAMPLE)
    rep = S.build(md)
    assert rep['ok'] and rep['name'] == 'short_yang' and rep['title'] == 'Inhibition in short'
    out = out_dir('short_yang')
    page = (out / 'index.html').read_text()
    assert '<title>Inhibition in short</title>' in page and '<h1 id="inhibition-in-short">Inhibition in <b>short</b></h1>' in page
    a_tip = 'KC / ORN data: a_dist: def of a_dist [kcq/data_targets.py abcdef0]'
    # a plain value, inv: and neg: with the main page's value, a literal, a table cell, a difference
    assert pv('0.12', a_tip) in page
    assert pv('0.88', 'shown as 1 − distance; the main page gives 0.12: ' + a_tip) in page
    assert pv('−0.65', 'sign reversed (a change of a correlation); the main page gives +0.65: KC / ORN data: rho: '
        'def of rho [kcq/data_targets.py abcdef0]') in page
    assert pv('0.05', 'hand-entered: Yang, by hand') in page
    assert pv('0.50', 'shown as 1 − distance; the main page gives 0.50: KC / ORN data: b_dist: def of b_dist '
        '[kcq/data_targets.py abcdef0]') in page
    assert '<td class="num">' + pv('0.3', 'results/E1/t.csv, row tag=x, column col_a') + '</td>' in page
    assert re.search(r'title="difference: \(KC / ORN data: a_dist: [^"]*\) minus \(KC / ORN data: b_dist: [^"]*\)">'
        r'−0\.377</span>', page)
    # meta:generated in the byline, with its source
    assert re.search(r'<p class="meta">Built <span class="pv" title="time this page was built">'
        r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}</span>\.</p>', page)
    # the heading with a value keeps the anchor of its plain text, in the table of contents too
    assert '<h2 id="the-data-at-0-65">The data at ' in page and '<a href="#the-data-at-0-65">The data at 0.65</a>' in page
    # values in attributes keep their text only: the figure's alt text, the hover definition
    assert 'alt="The figure at 0.65"' in page and '<figcaption>The figure at <span class="pv"' in page
    assert '<span class="def" title="the toxin, in 92 flies">TNT</span> blocks APL output.' in page
    assert '<details class="fold"><summary>More detail</summary>' in page and 'src="data:image/png;base64,' in page
    assert '⟦' not in page and '{{' not in page
    plain = (out / 'short_yang.md').read_text().split('\n')
    assert plain[8:10] == ['The mean distance is 0.12, a correlation of 0.88; APL changes it by',
        '−0.65, and a hand-entered 0.05.']
    filled = (out / 'short_yang_filled.md').read_text()
    assert '## The data at ⟦pv:0.65§' in filled and filled.split('\n')[6].endswith(' {#the-data-at-0-65}')
    side = json.loads((out / 'short_yang.values.json').read_text())
    assert side['n_values'] == 11 and side['n_inv'] == 2 and side['n_neg'] == 1
    inv = next(e for e in side['values'] if e.get('transform') == 'inv')
    assert inv['line'] == 9 and inv['main_text'] == '0.12' and inv['text'] == '0.88' and abs(inv['value'] - 0.877) < 1e-12
    assert rep['values'] == dict(filled=11, plain=8, inv=2, neg=1,
        by_source=dict(data=7, difference=1, literal=1, meta=1, table=1))
    assert json.loads((out / 'build_report.json').read_text()) == rep
    assert rep['tool'] == S.TOOL and rep['lint'] == [] and rep['fill_lint'] == []
    assert rep['inputs']['markdown_git'] == 'outside the repository' and rep['inputs']['shared'] == str(world)
    assert [(f['line'], f['path']) for f in rep['inputs']['figures']] == [(17, 'results/fig.png')]
    assert sorted(p.name for p in out.iterdir()) == ['build_report.json', 'index.html', 'short_yang.md',
        'short_yang.values.json', 'short_yang_filled.md']


def test_unresolved_placeholder_fails_and_removes_the_page(tmp_path):
    md = write_md(tmp_path, SAMPLE)
    S.build(md, 'page')
    out = out_dir('page')
    assert (out / 'index.html').exists()
    text = SAMPLE + '\nA missing value {{data:nope ~ .2f}} and {{inv:lit:abc # s}}.\n'
    n = text.split('\n').index('A missing value {{data:nope ~ .2f}} and {{inv:lit:abc # s}}.') + 1
    md.write_text(text)
    with pytest.raises(S.BuildError) as e:
        S.build(md, 'page')
    msg = str(e.value)
    assert msg.startswith('unresolved placeholders:') and f'line {n}: data:nope' in msg
    assert f'line {n}: inv:lit:abc # s: inv: needs a number' in msg
    # nothing of the earlier build is left but the report of the failure
    assert sorted(p.name for p in out.iterdir()) == ['build_report.json']
    rep = json.loads((out / 'build_report.json').read_text())
    assert rep['ok'] is False and rep['failed'] == 'unresolved placeholders' and rep['tool'] == S.TOOL
    assert any(f'line {n}: data:nope' in x for x in rep['errors']) and len(rep['removed']) == 4
    # and the next good build replaces it
    md.write_text(SAMPLE)
    assert S.build(md, 'page')['ok'] and (out / 'index.html').exists()


@pytest.mark.parametrize('extra, what, needle', [
    ('The **unclosed bold here.', 'the page does not pass md_to_artifact.lint', 'bold marker left in text'),
    ('A value with `an unclosed backtick {{data:rho ~ .2f}} here.', 'the page does not pass md_to_artifact.lint',
        'backtick left in text'),
    ('![a missing figure](results/missing.png)', 'missing figures', 'not found: results/missing.png'),
    ('A marker ⟦pv:1§x⟧ written by hand.', 'characters the filled page uses as markers', 'contains ⟦ or ⟧'),
])
def test_failures_name_the_line(tmp_path, extra, what, needle):
    text = SAMPLE + '\n' + extra + '\n'
    n = text.split('\n').index(extra) + 1
    md = write_md(tmp_path, text)
    out = out_dir('short_yang')
    out.mkdir(parents=True)
    (out / 'index.html').write_text('stale')
    (out / 'build_report.json').write_text(json.dumps(dict(tool=S.TOOL)))
    with pytest.raises(S.BuildError) as e:
        S.build(md)
    assert e.value.what == what and any(p.startswith(f'line {n}: ') and needle in p for p in e.value.problems), \
        e.value.problems
    assert not (out / 'index.html').exists()
    assert json.loads((out / 'build_report.json').read_text())['failed'] == what


def test_footnote_lint_names_the_heading_and_a_title_is_required(tmp_path):
    text = '# T\n\nText.\n\n## Footnotes {#footnotes}\n\n### A note {#fn-a}\n\nNot cited.\n'
    with pytest.raises(S.BuildError) as e:
        S.build(write_md(tmp_path, text), 'fn')
    assert e.value.problems == ['line 7: footnote #fn-a is not linked from the text']
    with pytest.raises(S.BuildError) as e:
        S.build(write_md(tmp_path, '## Only a section\n\nText.\n'), 'notitle')
    assert e.value.what == 'no title' and not (out_dir('notitle') / 'index.html').exists()


@pytest.mark.parametrize('name', ['summary', 'summary_yang', 'summary_all_experiments', 'Summary_Yang', 'SUMMARY',
    'summary-yang', 'summary.x', '', 'a/b', '../x', '.hidden', 'naïve', 'two words'])
def test_refused_names(tmp_path, name):
    with pytest.raises(S.Refused):
        S.check_name(name)
    with pytest.raises(S.Refused):
        S.build(write_md(tmp_path, SAMPLE), name)
    assert not S.ARTIFACTS.exists()


def test_names(tmp_path):
    assert S.default_name(S.DEFAULT_MARKDOWN) == 'inhibition_short_yang'
    assert S.default_name(Path('x/notes_yang.md')) == 'notes_yang'
    for ok in ('inhibition_short_yang', 'scratch_short_report_check', 'a-b.c', 'E1_short'):
        assert S.check_name(ok) == ok
    # a name derived from the file is checked as well
    with pytest.raises(S.Refused):
        S.build(write_md(tmp_path, SAMPLE, 'summary_notes.md'))
    # a folder another tool wrote is refused and left alone
    other = out_dir('questions_for_yang')
    other.mkdir(parents=True)
    (other / 'index.html').write_text('theirs')
    (other / 'build_report.json').write_text('{"generated": "x"}')
    with pytest.raises(S.Refused, match='another page'):
        S.build(write_md(tmp_path, SAMPLE), 'questions_for_yang')
    assert (other / 'index.html').read_text() == 'theirs'


def test_command_line(tmp_path, monkeypatch, capsys):
    md = write_md(tmp_path, SAMPLE)
    monkeypatch.setattr(sys, 'argv', ['yang_short_report.py', 'build', str(md), '--name', 'cli_page'])
    S.main()
    lines = capsys.readouterr().out.strip().split('\n')
    assert lines[0] == 'values: 11 filled (2 inv, 1 neg; by source: 7 data, 1 difference, 1 literal, 1 meta, 1 table)'
    wrote = [x for x in lines if x.startswith('wrote ')]
    assert len(wrote) == 5 and any('index.html' in x and '"Inhibition in short"' in x for x in wrote)
    assert lines[-1] == 'lint clean'
    # refused names exit with the usage error, before anything is written
    monkeypatch.setattr(sys, 'argv', ['yang_short_report.py', 'build', str(md), '--name', 'summary_yang'])
    with pytest.raises(SystemExit) as e:
        S.main()
    assert e.value.code == 2 and 'summary pages' in capsys.readouterr().err
    assert not out_dir('summary_yang').exists()
    # a failed build exits non-zero with the line named
    md.write_text(SAMPLE + '\nBroken {{data:nope}}.\n')
    monkeypatch.setattr(sys, 'argv', ['yang_short_report.py', 'build', str(md), '--name', 'cli_page'])
    with pytest.raises(SystemExit) as e:
        S.main()
    assert e.value.code not in (0, None) and 'line 31: data:nope' in str(e.value.code)
    assert 'no index.html written' in str(e.value.code)
    out = capsys.readouterr().out
    assert 'removed ' in out and 'index.html (from an earlier build)' in out
    assert not (out_dir('cli_page') / 'index.html').exists()
