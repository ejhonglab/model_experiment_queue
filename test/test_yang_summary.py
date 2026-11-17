"""Tests of reports/yang_summary.py (Yang's version of the summary page) on tiny in-memory templates, with the
placeholder resolver replaced by a fake one and the shared directory pointed at tmp_path: nothing here reads or writes
the shared directory, the repository's translation memory or outputs/_artifacts."""
import json
import os
import re
import tempfile

import pytest

from kcq import harness as H

try:
    H.shared_dir()
except FileNotFoundError:   # report_values looks the shared directory up on import; give it an empty one
    os.environ['KC_REMAPPING_SHARED'] = tempfile.mkdtemp(prefix='yang_summary_test_shared_')

from reports import md_to_artifact as M   # noqa: E402
from reports import report_values as R   # noqa: E402
from reports import yang_summary as Y   # noqa: E402

DATA = {'a_dist': 0.123, 'b_dist': 0.5, 'tiny': 0.001, 'rho': 0.65, 'count': 92}
TABLE = {'E1/t.csv:tag=x;model=a_rather_long_model_name_here:col_a,col_b': 0.3}
REAL_RESOLVE = R.resolve


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
def fake_world(monkeypatch, tmp_path):
    monkeypatch.setattr(R, 'resolve', fake_resolve)
    shared = tmp_path / 'shared'
    shared.mkdir()
    monkeypatch.setattr(Y, '_shared', lambda: shared)
    return shared


def paths(tmp_path, template):
    (tmp_path / 'data').mkdir(exist_ok=True)
    t = tmp_path / 'template.md'
    t.write_text(template)
    return Y.Paths(template=t, data=tmp_path / 'data', out=tmp_path / 'out')


def check(template, n, answer):
    """(template line, errors) of an answer to line n of a template."""
    lines = Y.analyze(template)
    out, errs, _ = Y.expand(answer, lines[n - 1].text, Y.Filler())
    if out is not None:
        errs = Y.structure_errors(lines[n - 1], out)
    return out, errs


# ---------------------------------------------------------------- pre-pass, compact form, tokens

def test_prepass():
    line = 'ρ on odor-odor distances, and the distance {{lit:3 # on odor-odor distances}} on odor-odor distances'
    assert Y.prepass(line) == ('ρ on odor-odor correlations, and the distance {{lit:3 # on odor-odor distances}} '
        'on odor-odor correlations')
    assert Y.prepass('the mean distance on odor pairs') == 'the mean distance on odor pairs'
    assert not Y.watch_hit('ρ on odor-odor distances {{data:a_dist}}')
    assert Y.watch_hit('pairs far apart') and Y.watch_hit('{{lit:1 # x}} sits above the fit')


def test_compact_form():
    line = ('| mean pair distance | {{data:a_dist ~ .2f}} | {{E1/t.csv:tag=x;model=a_rather_long_model_name_here:col_a,col_b ~ .1f}} '
        '| {{lit:0.05 # src}} | {{diff:data:a_dist -- data:b_dist ~ +.3f}} | {{data:a_dist.lo ~ .2f}} |')
    c = Y.compact(line, Y.Filler())
    assert c.startswith('| mean pair distance | ⟨1|0.12|a_dist⟩ | ⟨2|0.3|…')
    assert ':col_a⟩' in c and len(re.search(r'⟨2\|0\.3\|([^:]*):', c).group(1)) == 30
    assert '| ⟨3|0.05|⟩ | ⟨4|−0.377|diff(a_dist − b_dist)⟩ | ⟨5|??|a_dist.lo⟩ |' in c
    assert Y.compact('edited: {{meta:generated}}', Y.Filler()).endswith('|meta:generated⟩')
    assert Y._tok('a|b⟩c⟨d') == 'a/b/c/d'


def test_token_round_trip():
    src = ('The pairs are {{data:a_dist ~ .2f}} apart, from {{data:b_dist ~ .2f}} to {{lit:0.88 # src x}} '
        '({{diff:data:a_dist -- data:b_dist ~ +.3f}}; ρ {{data:rho ~ .2f}}, n {{data:count ~ .0f}}).')
    f = Y.Filler()
    c = Y.compact(src, f)
    toks = re.findall(r'⟨[^⟩]*⟩', c)
    # reordered, inv, neg, lit, the source token copied verbatim and the short form
    answer = (f'The pairs have a correlation of ⟨1|inv⟩, from ⟨3|lit|0.12⟩ to ⟨2|inv⟩ (⟨4|neg⟩; ρ {toks[4]}, '
        'n ⟨6⟩).')
    out, errs, warns = Y.expand(answer, src, f)
    assert errs == [] and warns == []
    assert out == ('The pairs have a correlation of {{inv:data:a_dist ~ .2f}}, from {{lit:0.12 # Yang\'s version, '
        'rewritten from "0.88"; src x}} to {{inv:data:b_dist ~ .2f}} ({{neg:diff:data:a_dist -- data:b_dist ~ +.3f}}; '
        'ρ {{data:rho ~ .2f}}, n {{data:count ~ .0f}}).')
    assert Y.Filler().line(out, 't', tooltips=False) == 'The pairs have a correlation of 0.88, from 0.12 to 0.50 (+0.377; ρ 0.65, n 92).'
    # the source line in compact form unchanged expands back to the source
    assert Y.expand(c, src, f)[0] == src
    # a kept meta: token from an earlier export (another build time) is still the source token
    assert Y.expand('edited: ⟨1|2020-01-01 00:00|meta:generated⟩', 'edited: {{meta:generated}}', f)[0] == 'edited: {{meta:generated}}'


TEMPLATE_V = """# Title

## Pairs far apart

## Section distance {#sec}

- **Lead:** the distance is {{data:a_dist ~ .2f}} and {{meta:generated}}.

1. the distance

The distance [is fitted](#fn-isotonic) with `code` and !button[m](results/x.png) {{data:missing}}.

| | ORNs | KCs {span=2} |
|---|---|---|---|
| a distance | {{data:a_dist ~ .2f}} | b | c |

??? The figure with distance

![caption distance](results/a.png)

%%% hover
VCF: Odor distances
%%%
"""


@pytest.mark.parametrize('n, answer, reason', [
    (7, '- **Lead:** the correlation is ⟨1|inv⟩.', 'token ⟨2⟩ is missing'),
    (7, '- **Lead:** ⟨1|inv⟩ ⟨1⟩ ⟨2⟩.', 'token ⟨1⟩ is used 2 times'),
    (7, '- **Lead:** ⟨1|inv⟩ ⟨2⟩ ⟨3⟩.', 'the source line has 2 placeholder(s)'),
    (7, '- **Lead:** ⟨a|inv⟩ ⟨1⟩ ⟨2⟩.', 'does not start with a placeholder number'),
    (7, '- **Lead:** ⟨1⟩ ⟨2|inv⟩.', 'needs a number'),
    (7, '- **Lead:** ⟨1|lit|0.9⟩ ⟨2⟩.', 'lit only replaces a hand-entered value'),
    (7, '- **Lead:** ⟨1|0.88|a_dist⟩ ⟨2⟩.', 'changes the value'),
    (7, '- **Lead:** ⟨1⟩ ⟨2⟩ {{data:b_dist}}.', 'contains {{ }}'),
    (7, '- **Lead:** ⟨1⟩ ⟨2⟩ ⟪x⟫.', 'contains ⟪ ⟫'),
    (11, 'The correlation [is fitted](#fn-isotonic) with `code` and !button[m](results/x.png) ⟨1|inv⟩.', 'does not resolve'),
    (3, '### Pairs decorrelated', 'heading level'),
    (3, '## Pairs decorrelated {#pairs}', 'heading id changed'),
    (5, '## Section correlation {#other}', 'heading id changed'),
    (5, '## Section correlation', 'heading id changed'),
    (7, '* **Lead:** the correlation is ⟨1|inv⟩ and ⟨2⟩.', 'list marker'),
    (9, '2. the correlation', 'list marker'),
    (9, 'the correlation', 'the source is a list line, the answer a para line'),
    (11, '- The correlation [is fitted](#fn-isotonic) with `code` and !button[m](results/x.png) ⟨1⟩.', 'the answer a list line'),
    (11, 'The correlation [is fitted](#fn-other) with `code` and !button[m](results/x.png) ⟨1⟩.', 'link targets changed'),
    (11, 'The correlation [is fitted](#fn-isotonic) with `code` and !button[n](results/x.png) ⟨1⟩.', '!button'),
    (11, 'The **correlation [is fitted](#fn-isotonic) with `code` and !button[m](results/x.png) ⟨1⟩.', '** occurs an odd'),
    (11, 'The correlation [is fitted](#fn-isotonic) with `code and !button[m](results/x.png) ⟨1⟩.', 'backticks'),
    (13, '| | ORNs {span=2} | KCs |', 'span'),
    (15, '| a correlation | ⟨1|inv⟩ | b |', '3 cells instead of 4'),
    (14, '|---|---|---|', 'cannot be converted'),
    (17, 'The figure with correlation', 'the source is a fold line'),
    (19, '![caption correlation](results/b.png)', 'figure path changed'),
    (19, '![caption correlation](results/summary_yang/missing.png)', 'figure path changed'),
    (22, 'VCF2: Odor correlations', 'hover term'),
    (21, '%%% hover x', 'cannot be converted'),
])
def test_validation_failures(n, answer, reason):
    out, errs = check(TEMPLATE_V, n, answer)
    assert any(reason in e for e in errs), errs


def test_validation_passes(fake_world):
    assert check(TEMPLATE_V, 3, '## Pairs decorrelated')[1] == []
    assert check(TEMPLATE_V, 5, '## Section correlation {#sec}')[1] == []
    assert check(TEMPLATE_V, 13, '| | ORNs | KCs {span=2} |')[1] == []
    assert check(TEMPLATE_V, 17, '??? The figure with correlation')[1] == []
    assert check(TEMPLATE_V, 22, 'VCF: Odor correlations (the natural-source distances reversed)')[1] == []
    (fake_world / 'results' / 'summary_yang').mkdir(parents=True)
    (fake_world / 'results' / 'summary_yang' / 'c.png').write_bytes(b'png')
    assert check(TEMPLATE_V, 19, '![caption correlation](results/summary_yang/c.png)')[1] == []
    # lit text: none, or with # or |, is refused
    src = 'x {{lit:0.88 # src}} apart'
    for bad in ('x ⟨1|lit|⟩ y', 'x ⟨1|lit|0.1 # mine⟩ y', 'x ⟨1|lit|0.1|0.2⟩ y'):
        assert any('lit|new text' in e for e in Y.expand(bad, src, Y.Filler())[1])
    assert Y.expand('x ⟨1|inv⟩ y', 'x {{data:count ~ .0f}} y', Y.Filler())[2]   # inv on 92: a warning


def test_distance_left():
    assert Y.distance_left('the mean distance is ⟨1|inv⟩')
    assert not Y.distance_left('the total variation distance and ⟨1|distance hint⟩ and 1 - distance')
    assert not Y.distance_left('drawn as distances (distance = 1 − correlation)')


# ---------------------------------------------------------------- import, memory, status, export

TEMPLATE_I = """# Title

## One

The pairs are {{data:a_dist ~ .2f}} apart.

A context line.

The distance of one recording is {{data:b_dist ~ .2f}}.

## Two

| statistic | value |
|---|---|
| mean distance | {{data:a_dist ~ .2f}} |
| ρ | {{data:rho ~ .2f}} |
"""


def test_import_keep_memory_and_status(tmp_path):
    p = paths(tmp_path, TEMPLATE_I)
    ans = tmp_path / 'answers.out.md'
    ans.write_text('header text the agent copied\n\n@@ 5 *\nThe pairs have a correlation of ⟨1|inv⟩.\n\n'
        '@@ 9 KEEP\n@@ 15\n| mean distance | ⟨1|inv⟩ |\n@@ 16\n| ρ | ⟨1|inv⟩ | extra |\n@@ 7\n'
        'one line\nand another\n@@ 99\nx\n')
    n, failures, warnings = Y.import_answers([ans], p, by='test agent')
    assert n == 3
    reasons = {x['line']: x['reason'] for x in failures}
    assert set(reasons) == {16, 7, 99}
    assert '3 cells instead of 2' in reasons[16] and 'exactly one line' in reasons[7] and 'not in the template' in reasons[99]
    assert any('line 15' in w and 'distanc' in w for w in warnings)
    fail_md = (p.convert / 'import_failures.md').read_text()
    assert '@@ 16 *\n| ρ | ⟨1|0.65|rho⟩ |' in fail_md and 'answer given:\n\n| ρ | ⟨1|inv⟩ | extra |' in fail_md
    mem = [json.loads(t) for t in p.memory.read_text().splitlines()]
    assert [e['line'] for e in mem] == [5, 9, 15] and list(mem[0]) == list(Y.MEMORY_KEYS)
    assert mem[0]['out'] == 'The pairs have a correlation of {{inv:data:a_dist ~ .2f}}.' and not mem[0]['keep']
    assert mem[1]['keep'] and mem[1]['out'] == mem[1]['src'] and mem[0]['by'] == 'test agent'
    assert mem[0]['sha'] == Y.sha('The pairs are {{data:a_dist ~ .2f}} apart.')
    assert re.fullmatch(r'\d{4}-\d\d-\d\d \d\d:\d\d', mem[0]['when'])
    s = Y.status(p)
    # lines 13 and 16 are missing: rows of a table with a starred row
    assert (s['converted'], s['keep'], s['missing'], s['stale']) == (2, 1, 2, [])
    assert list(s['by_section'].values()) == [[13, 16]]
    # re-importing replaces the entry; the template changing a line leaves its entry stale
    ans.write_text('@@ 5\nThe pairs are decorrelated, ⟨1|inv⟩.\n@@ 16 KEEP\n')
    n, failures, warnings = Y.import_answers([ans], p)
    assert n == 2 and failures == [] and not (p.convert / 'import_failures.md').exists()
    mem = {e['line']: e for e in map(json.loads, p.memory.read_text().splitlines())}
    assert mem[5]['out'] == 'The pairs are decorrelated, {{inv:data:a_dist ~ .2f}}.' and mem[5]['by'] == 'unspecified'
    p.template.write_text(TEMPLATE_I.replace('apart.', 'apart!'))
    assert Y.status(p)['stale'] == [5]


def test_import_warns_on_unanswered_starred_lines(tmp_path):
    p = paths(tmp_path, TEMPLATE_I)
    Y.export(p)
    ans = p.convert / 'chunk_01.out.md'
    ans.write_text('@@ 5\nThe pairs have a correlation of ⟨1|inv⟩.\n')
    n, failures, warnings = Y.import_answers([ans], p)
    assert n == 1 and failures == []
    assert any('chunk_01.out.md: 4 starred line(s) of lines 1-17 got no answer: 9, 13, 15, 16' in w for w in warnings), warnings
    assert 'got no answer' in (p.convert / 'import_warnings.md').read_text()


def test_export_chunks(tmp_path):
    body = '\n\n'.join(f'## Section {k}\n\nA plain line {k} of some length to fill the chunk {"x" * 300}.'
        + ('\n\nThe pairs are {{data:a_dist ~ .2f}} apart.' if k in (1, 5) else '') for k in range(1, 7))
    p = paths(tmp_path, '# Title\n\n' + body + '\n')
    p.convert.mkdir(parents=True)
    (p.convert / 'chunk_09.md').write_text('old')
    (p.convert / 'chunk_09.out.md').write_text('kept')
    rows = Y.export(p, chunk_kb=0.7)
    assert [r['written'] for r in rows].count(True) == 2 and len(rows) > 2
    names = sorted(x.name for x in p.convert.glob('chunk_*.md'))
    assert 'chunk_09.out.md' in names and 'chunk_09.md' not in names
    first = (p.convert / 'chunk_01.md').read_text().split('\n')
    assert first[0].startswith(f'Chunk 1 of {len(rows)}: lines 1-') and 'conversion_rules.md' in first[1]
    assert '@@ 7 *\nThe pairs are ⟨1|0.12|a_dist⟩ apart.' in (p.convert / 'chunk_01.md').read_text()
    later = [r for r in rows if r['written']][1]
    # numbering by position: a later chunk keeps its number though the ones before it are not written
    assert later['chunk'] != 'chunk_02' and (p.convert / f'{later["chunk"]}.md').exists()
    all_rows = Y.export(p, all_chunks=True, chunk_kb=0.7)
    assert [r['chunk'] for r in all_rows] == [r['chunk'] for r in rows] and all(r['written'] for r in all_rows)
    # a line in the memory is shown converted, filled, with '='
    Y.write_memory(p.memory, {Y.sha('The pairs are {{data:a_dist ~ .2f}} apart.'): dict(sha=Y.sha(
        'The pairs are {{data:a_dist ~ .2f}} apart.'), line=7, src='The pairs are {{data:a_dist ~ .2f}} apart.',
        out='The pairs have a correlation of {{inv:data:a_dist ~ .2f}}.', keep=False, by='x', when='t')}, Y.analyze(p.template.read_text()))
    rows = Y.export(p, chunk_kb=0.7)
    assert not any(r['written'] for r in rows)
    Y.export(p, all_chunks=True, chunk_kb=0.7)
    assert '= \nThe pairs have a correlation of 0.88.' not in (p.convert / 'chunk_01.md').read_text()
    assert '@@ 7 =\nThe pairs have a correlation of 0.88.' in (p.convert / 'chunk_01.md').read_text()


def test_parse_answers_format():
    e = Y.parse_answers('@@ 3 * KEEP\n@@ 4 KEEP\nx\n@@ x\ny\n@@ 5 maybe\nz\n@@ 6\n')
    assert [x['error'] for x in e][0] is None and e[0]['keep']
    assert 'KEEP takes no converted line' in e[1]['error'] and 'not a line number' in e[2]['error']
    assert 'unexpected text' in e[3]['error'] and 'no converted line' in e[4]['error']


# ---------------------------------------------------------------- values

def test_inv_neg_formatting():
    f = Y.Filler()
    assert f.line('{{inv:data:a_dist ~ .2f}} {{neg:data:a_dist ~ +.3f}} {{neg:data:a_dist ~ .2f}} {{neg:data:tiny ~ .2f}}',
        't', tooltips=False) == '0.88 −0.123 −0.12 0.00'
    assert f.line('{{inv:lit:0.05 # s}} {{neg:lit:+0.120 # s}} {{inv:lit:5 # s}} {{neg:lit:−0.03 # s}} {{inv:lit:0.05 # s ~ .3f}}',
        't', tooltips=False) == '0.95 −0.120 −4 0.03 0.950'
    assert f.line('{{neg:diff:data:a_dist -- data:b_dist ~ +.3f}}', 't', tooltips=False) == '+0.377'
    tip = f.line('{{inv:data:a_dist ~ .2f}}', 't', tooltips=True)
    assert tip.startswith('⟦pv:0.88§shown as 1 − distance; the main page gives 0.12: KC / ORN data: a_dist: def of a_dist')
    tip = f.line('{{neg:data:a_dist ~ +.3f}}', 't', tooltips=True)
    assert tip.startswith('⟦pv:−0.123§sign reversed (a change of a correlation); the main page gives +0.123: ')
    f.check()
    side = []
    f2 = Y.Filler()
    f2.line('{{inv:data:a_dist ~ .2f}} {{data:b_dist ~ .1f}}', 't', tooltips=False, record=dict(line=1))
    assert f2.side[0]['transform'] == 'inv' and f2.side[0]['main_text'] == '0.12' and abs(f2.side[0]['value'] - 0.877) < 1e-12
    assert 'transform' not in f2.side[1] and f2.side[1]['line'] == 1
    f3 = Y.Filler()
    assert f3.line('a {{inv:lit:abc # s}} b {{inv:meta:generated}} {{data:nope}}', 'line 7', tooltips=False) == 'a ?? b ?? ??'
    with pytest.raises(KeyError) as e:
        f3.check()
    msg = e.value.args[0]
    assert 'line 7: inv:lit:abc # s: inv: needs a number' in msg and 'line 7: inv:meta:generated' in msg and 'line 7: data:nope' in msg
    del side


# ---------------------------------------------------------------- Yang's lines

TEMPLATE_Y = """# Title

## Pairs far apart

### ρ on odor-odor distances

## Kept heading

The pairs are far apart {{data:a_dist ~ .2f}}.

- the distance item

| statistic | distance |
|---|---|
| a | {{data:a_dist ~ .2f}} |

??? The figure with distance

![caption distance](results/a.png)

%%% hover
VCF: Odor distances
%%%

ρ on odor-odor distances {{data:rho ~ .2f}}.

Converted already: distance.
"""


def test_yang_lines_heading_ids_and_unreviewed():
    lines = Y.analyze(TEMPLATE_Y)
    conv = 'Converted already: distance.'
    mem = {Y.sha('## Pairs far apart'): dict(out='## Pairs decorrelated', keep=False),
        Y.sha(conv): dict(out='Converted already: correlation.', keep=False)}
    out, info = Y.yang_lines(lines, mem, Y.starred(lines))
    assert out[2] == '## Pairs decorrelated {#pairs-far-apart}'
    assert out[4] == '### ρ on odor-odor correlations {#on-odor-odor-distances}' == '### ρ on odor-odor correlations {#' + M._slug('ρ on odor-odor distances') + '}'
    assert out[6] == '## Kept heading'
    assert out[8] == '⟪ydist⟫The pairs are far apart {{data:a_dist ~ .2f}}.'
    assert out[10] == '- ⟪ydist⟫the distance item'
    assert out[12] == '| statistic | distance |' and out[13] == '|---|---|' and out[14] == '| ⟪ydist⟫a | {{data:a_dist ~ .2f}} |'
    assert out[16] == '??? ⟪ydist⟫The figure with distance' and out[18] == '![caption distance](results/a.png)'
    assert out[21] == 'VCF: Odor distances'
    assert out[24] == 'ρ on odor-odor correlations {{data:rho ~ .2f}}.'
    assert out[26] == 'Converted already: correlation.'
    assert info['converted'] == [3, 27] and info['unreviewed'] == [9, 11, 13, 15, 17, 22]
    assert info['marked'] == [9, 11, 15, 17]
    assert Y.keep_anchor('## A {#x}', '## B {#x}') == '## B {#x}' and Y.keep_anchor('## A', '## A') == '## A'


# ---------------------------------------------------------------- questions

QUESTIONS = """# Yang's questions

Free text, ignored.

## q-one
asked: 2026-11-16
sentence: The mean correlation is 0.88 in flies.
sentence_main: The mean distance is 0.12 in flies.
question: Why?

Because {{inv:data:a_dist ~ .2f}}.

- a list

## q-two
sentence_main: the earlier correction (too low)
question: Too low?

| statistic | x |
|---|---|
| a | {{data:a_dist ~ .2f}} |
"""


def test_parse_questions():
    qs = Y.parse_questions(QUESTIONS)
    assert [q['id'] for q in qs] == ['q-one', 'q-two']
    assert qs[0]['asked'] == '2026-11-16' and qs[0]['sentence'].startswith('The mean correlation')
    assert qs[0]['answer'] == 'Because {{inv:data:a_dist ~ .2f}}.\n\n- a list' and qs[1]['sentence'] == ''
    assert qs[1]['answer'].startswith('| statistic | x |')


@pytest.mark.parametrize('text, reason', [
    ('## q1\nquestion: x\n\nanswer', 'needs sentence: or sentence_main:'),
    ('## q1\nsentence: s\n\nanswer', 'needs question:'),
    ('## q1\nsentence: s\nquestion: x\n\n', 'the answer is empty'),
    ('## q1\nsentence: s\nquestion: x\n\na\n## q1\nsentence: s\nquestion: y\n\nb', 'duplicate id'),
    ('## Q 1\nsentence: s\nquestion: x\n\na', 'the id must be'),
    ('## q1\nsentence: s\nquestoin: x\n\na', 'expected one of'),
    ('## q1\nsentence: s\nsentence: t\nquestion: x\n\na', 'sentence given twice'),
])
def test_parse_questions_errors(text, reason):
    with pytest.raises(Y.InputError) as e:
        Y.parse_questions(text)
    assert reason in str(e.value)


def q(qid, sentence='', sentence_main=''):
    return dict(id=qid, asked='', sentence=sentence, sentence_main=sentence_main, question='Why?', answer='a')


def wrapped(doc, questions, main=lambda i: None):
    warnings = []
    anchors = Y.locate(questions, doc, main, warnings)
    lines = [t for t, _ in doc]
    groups = Y.wrap(lines, anchors, [x['id'] for x in questions])
    return lines, groups, anchors, warnings


def test_sentence_matching():
    doc = [('The mean is ⟦pv:0.88§tip [x]⟧ in flies. Next sentence.', 'para'),
        ('- **Lead:** first sentence. The second ⟦pv:0.5§t⟧ sentence here.', 'list'),
        ('| the same by the earlier correction (too low) | ⟦pv:0.79§t⟧ |', 'table'),
        ('The control matrix [is fitted with the best monotone](#fn-isotonic) function of it.', 'para'),
        ('**The remapping in control flies.** KCs push pairs apart.', 'para'),
        ('The odor–odor value is ⟦pv:−0.03§t⟧ here.', 'para')]
    qs = [q('changed-numbers', 'The mean is 0.12 in flies.'), q('in-list', 'The second 0.7 sentence here.'),
        q('in-cell', 'the same by the earlier correction (too low)'), q('across-cells', 'too low) 0.79'),
        q('across-link', 'matrix is fitted with the best'), q('across-bold', 'control flies. KCs push'),
        q('value-end', 'The mean is 0.88'), q('minus', 'The odor-odor value is -0.05 here.')]
    lines, groups, anchors, warnings = wrapped(doc, qs)
    ids = {g[0]: k for k, g in enumerate(groups)}
    assert lines[0] == f'⟪yq:{ids["changed-numbers"]}⟫The mean is ⟦pv:0.88§tip [x]⟧ in flies.⟪/yq⟫ Next sentence.'
    assert anchors['value-end'].start == anchors['changed-numbers'].start   # overlapping spans share one wrapper
    assert ['changed-numbers', 'value-end'] in groups
    assert lines[1] == f'- **Lead:** first sentence. ⟪yq:{ids["in-list"]}⟫The second ⟦pv:0.5§t⟧ sentence here.⟪/yq⟫'
    assert lines[2] == f'| ⟪yq:{ids["in-cell"]}⟫the same by the earlier correction (too low)⟪/yq⟫ | ⟦pv:0.79§t⟧ |'
    assert 'across-cells' not in anchors and any('across-cells' in w and 'unanchored' in w for w in warnings)
    assert lines[3] == (f'The control ⟪yq:{ids["across-link"]}⟫matrix [is fitted with the best monotone](#fn-isotonic)⟪/yq⟫'
        ' function of it.')
    assert lines[4] == f'⟪yq:{ids["across-bold"]}⟫**The remapping in control flies.** KCs push⟪/yq⟫ pairs apart.'
    assert lines[5] == f'⟪yq:{ids["minus"]}⟫The odor–odor value is ⟦pv:−0.03§t⟧ here.⟪/yq⟫'
    assert all(a.how == 'sentence' for a in anchors.values())


def test_sentence_fallbacks():
    doc = [('Text before.', 'para'),
        ('Totally reworded correlation line ⟦pv:0.88§t⟧.', 'para'),
        ('| a | ⟪ydist⟫the correlation |', 'table'),
        ('ρ on odor-odor correlations is high.', 'para'),
        ('Twice here.', 'para'), ('Twice here.', 'para')]
    main = {1: ('The original distance line ⟦pv:0.12§t⟧.', 'para'), 2: ('| a | the old cell wording |', 'table')}
    qs = [q('whole', 'Not there at all.', 'The original distance line 0.5.'),
        q('whole-cell', sentence_main='the old cell wording'),
        q('prepassed', sentence_main='ρ on odor-odor distances is high.'),
        q('gone', 'Nowhere.', 'Nowhere either.'), q('twice', 'Twice here.')]
    lines, groups, anchors, warnings = wrapped(doc, qs, main=lambda i: main.get(i))
    ids = {g[0]: k for k, g in enumerate(groups)}
    assert anchors['whole'].how == 'whole' and lines[1] == f'⟪yq:{ids["whole"]}⟫Totally reworded correlation line ⟦pv:0.88§t⟧.⟪/yq⟫'
    assert lines[2] == f'| a | ⟪ydist⟫⟪yq:{ids["whole-cell"]}⟫the correlation⟪/yq⟫ |'
    assert anchors['prepassed'].how == 'sentence_main'
    assert 'gone' not in anchors and anchors['twice'].doc == 4
    assert sum('whole passage marked' in w for w in warnings) == 2
    assert any('twice' in w and 'occurs in 2 lines' in w for w in warnings)
    block = Y.questions_block(qs, anchors)
    assert block[0] == "??? Yang's questions (5)" and block[-2:] == ['???', '']
    assert '- [Why?](#yq-whole)' in block and '- Why? (the sentence is no longer on the page)' in block
    assert Y._cut('word ' * 60).endswith(' …') and len(Y._cut('word ' * 60)) <= 142


# ---------------------------------------------------------------- the page

def test_postprocess(tmp_path):
    (tmp_path / 'f.png').write_bytes(b'\x89PNG fake')
    md = ('# Title\n\nbyline\n\n[TOC]\n\n## Head ⟪yq:0⟫distance⟪/yq⟫ {#h}\n\n⟪ydist⟫Some ⟪yq:1⟫sentence **bold**⟪/yq⟫ text.\n\n'
        '- item ⟪yq:2⟫with [a link](#h)⟪/yq⟫\n\n| statistic | x |\n|---|---|\n| ⟪yq:3⟫cell⟪/yq⟫ | 1 |\n\n'
        '![a ⟪yq:4⟫caption⟪/yq⟫](f.png)\n')
    page = M.convert_markdown(md, tmp_path, 'T')
    groups = [['h1'], ['p1', 'p2'], ['l1'], ['c1'], ['f1']]
    answers = {i: f'<div class="yq-q">q {i}</div><div class="yq-a"><p>a</p></div>' for g in groups for i in g}
    out, warnings = Y.postprocess(page, groups, answers)
    assert '⟪' not in out and '⟫' not in out and warnings == []
    assert '<a id="yq-p2"></a><span class="yq" id="yq-p1" data-yq="p1 p2" tabindex="0" role="button" aria-haspopup="dialog">sentence <b>bold</b></span>' in out
    assert '<h2 id="h">Head <span class="yq" id="yq-h1"' in out
    assert re.search(r'<nav class="toc">.*?</nav>', out, re.S).group(0).count('yq') == 0
    assert '<span class="ydist" title="This passage has not been converted yet' in out and '>distances</span> Some' in out
    assert 'alt="a caption"' in out and '<figcaption>a <span class="yq" id="yq-f1"' in out
    assert out.count('<template id="yqa-') == 6 and out.index('<template id="yqa-h1">') < out.index('</main>')
    assert 'div.yq-pop' not in out and '.yq-pop {' in out and "Yang's question" in out
    assert M.lint(out) == []
    # a marked passage that crosses an element is widened to its block
    bad = '<html><head><style></style></head><body><main><p>Start ⟪yq:0⟫a <b>bold⟪/yq⟫ text</b> end.</p></main></body></html>'
    out, warnings = Y.postprocess(bad, [['u1']], {})
    assert '<p><span class="yq" id="yq-u1" data-yq="u1" tabindex="0" role="button" aria-haspopup="dialog">Start a <b>bold text</b> end.</span></p>' in out
    assert warnings == ['u1: the sentence crossed an HTML element; the whole p is marked']
    assert Y.balanced('a <b>x</b><br> <span class="pv" title="t">1</span>') and not Y.balanced('a <b>x')


def test_build_end_to_end(tmp_path, fake_world):
    template = ('# Title\n\nlast model experiment: x, edited: {{meta:generated}}\n\n[TOC]\n\n## Data\n\n'
        'The mean distance is {{data:a_dist ~ .2f}} in flies. Another sentence.\n\nSee [the note](#fn-x).\n\n'
        '- **Noise floor:** the noise floor is {{data:b_dist ~ .2f}} here.\n\n'
        '| statistic | control flies |\n|---|---|\n| mean pair distance | {{data:a_dist ~ .2f}} |\n'
        '| the earlier correction (too low) | {{data:b_dist ~ .2f}} |\n\n## Footnotes {#footnotes}\n\n### Note {#fn-x}\n\n'
        'Text.\n')
    p = paths(tmp_path, template)
    (p.data / 'intro.md').write_text("**Yang's version.** Correlation = 1 − distance, {{inv:data:a_dist ~ .2f}}.\n")
    (p.data / 'questions.md').write_text('# Q\n\n## q-mean\nasked: 2026-11-16\nsentence: The mean correlation is 0.9 in flies.\n'
        'question: Which mean?\n\nThe mean over {{inv:data:a_dist ~ .2f}} pairs.\n\n### A heading in an answer\n\nmore.\n\n'
        '## q-low\nsentence_main: the earlier correction (too low)\nquestion: Too low?\n\nYes.\n\n'
        '## q-gone\nsentence: Not on the page.\nquestion: Where?\n\nGone.\n')
    lines = Y.analyze(template)
    src = 'The mean distance is {{data:a_dist ~ .2f}} in flies. Another sentence.'
    Y.write_memory(p.memory, {Y.sha(src): dict(sha=Y.sha(src), line=9, src=src,
        out='The mean correlation is {{inv:data:a_dist ~ .2f}} in flies. Another sentence.', keep=False, by='t', when='t')}, lines)
    ok, rep = Y.build(p, skip_figure=True, shared=fake_world)
    assert ok, rep['lint']
    page = (p.out / 'index.html').read_text()
    assert '<title>Model experiments (Yang&#x27;s version)</title>' in page or "<title>Model experiments (Yang's version)</title>" in page
    assert 'id="yq-q-mean"' in page and 'id="yq-q-low"' in page and 'id="yq-q-gone"' not in page
    assert '<template id="yqa-q-mean">' in page and '<template id="yqa-q-gone">' not in page
    assert '<b>A heading in an answer</b>' in page and 'Yang, 2026-11-16' in page
    assert page.index("Yang's version.") < page.index('<nav class="toc">') and page.index("Yang's questions (3)") < page.index('<nav class="toc">')
    assert 'Where? (the sentence is no longer on the page)' in page and 'href="#yq-q-mean"' in page
    assert page.count('class="ydist"') == 3   # the list item and the two table body rows; not the header row
    plain = (p.out / 'summary_yang.md').read_text()
    assert '⟪' not in plain and '⟦' not in plain and 'The mean correlation is 0.88 in flies.' in plain
    side = json.loads((p.out / 'summary_yang.values.json').read_text())
    assert side['n_inv'] == 3 and side['n_values'] == 7   # 5 in the template, 1 in the intro, 1 in an answer
    vals = [v for v in side['values'] if v.get('template_line') == 9]
    assert vals[0]['line'] == plain.split('\n').index('The mean correlation is 0.88 in flies. Another sentence.') + 1
    rq = rep['questions']
    assert rq['anchored'] == ['q-mean', 'q-low'] and rq['unanchored'] == ['q-gone'] and rq['whole_passage'] == []
    assert rep['lines']['converted'] == 1 and rep['lines']['unreviewed_lines'] == [13, 15, 17, 18]
    assert rep['values'] == dict(filled=7, inv=3, neg=0)
    assert json.loads((p.out / 'build_report.json').read_text())['questions']['unanchored'] == ['q-gone']
    # an unresolved placeholder in the memory fails the build, naming every key
    Y.write_memory(p.memory, {Y.sha(src): dict(sha=Y.sha(src), line=9, src=src,
        out='The mean correlation is {{inv:data:nope ~ .2f}} and {{neg:lit:abc # s}}.', keep=False, by='t', when='t')}, lines)
    with pytest.raises(KeyError) as e:
        Y.build(p, skip_figure=True, shared=fake_world)
    assert 'line 9: inv:data:nope' in e.value.args[0] and 'line 9: neg:lit:abc # s' in e.value.args[0]


# ---------------------------------------------------------------- carry

# the template before the split (one line per case) and after it; the memory holds entries for the old lines
OLD = {
    'para': 'The pairs are {{data:a_dist ~ .2f}} apart. They move closer with APL off. The ρ is {{data:rho ~ .2f}}.',
    'list': 'Three things follow. First, the pairs move apart. Second, the noise floor rises. Third, nothing else.',
    'sub': '- **Lead:** the distance is {{data:a_dist ~ .2f}}. It is far from {{data:b_dist ~ .2f}} and {{lit:0.88 # src}}.',
    'bold': '**The block is large. It sits above the curve. More text here, not about the block.**',
    'lead': '**The lead is long. It says the pairs move apart.** The rest is plain. It ends here.',
    'renum': 'A plain first sentence. The distances run from {{data:a_dist ~ .2f}} to {{data:b_dist ~ .2f}} ({{data:count ~ .0f}} pairs).',
    'nomatch': 'The pairs are far apart. The noise floor is low.',
    'inside': 'The pairs are far apart and the noise floor is low, which is all.',
    'changed': 'The old wording of a line about the distance.',
}
OUT = {
    'para': 'The pairs have a correlation of {{inv:data:a_dist ~ .2f}}. They become more correlated with APL off. '
        'The ρ is {{data:rho ~ .2f}}.',
    'list': None,   # KEEP
    'sub': '- **Lead:** the correlation is {{inv:data:a_dist ~ .2f}}. It is weakly correlated with '
        '{{inv:data:b_dist ~ .2f}} and {{lit:0.12 # Yang\'s version, rewritten from "0.88"; src}}.',
    'bold': '**The block is large. It sits below the curve. More text here, not about the block.**',
    'lead': '**The lead is long. It says the pairs are decorrelated.** The rest is plain. It ends here.',
    'renum': 'A plain first sentence. The correlations run from {{inv:data:b_dist ~ .2f}} to {{inv:data:a_dist ~ .2f}} '
        '({{data:count ~ .0f}} pairs).',
    'nomatch': 'The pairs are decorrelated and the noise ceiling is high.',
    'inside': 'The pairs are decorrelated and the noise ceiling is high, which is all.',
    'changed': 'The old wording of a line about the correlation.',
}
NEW = """# Title

## Data

The pairs are {{data:a_dist ~ .2f}} apart. They move closer with APL off.

The ρ is {{data:rho ~ .2f}}.

Three things follow.

- First, the pairs move apart.
- Second, the noise floor rises.
- Third, nothing else.

- **Lead:** the distance is {{data:a_dist ~ .2f}}.
  - It is far from {{data:b_dist ~ .2f}} and {{lit:0.88 # src}}.

**The block is large.** It sits above the curve.

More text here, not about the block.

**The lead is long.** It says the pairs move apart.

The rest is plain. It ends here.

A plain first sentence.

The distances run from {{data:a_dist ~ .2f}} to {{data:b_dist ~ .2f}} ({{data:count ~ .0f}} pairs).

The pairs are far apart.

The noise floor is low.

The pairs are far apart

and the noise floor is low, which is all.

The new wording of a line about the distance.
"""


def carry_world(tmp_path, cases=tuple(OLD)):
    p = paths(tmp_path, NEW)
    mem = {}
    for c in cases:
        src, out = OLD[c], OUT[c] or OLD[c]
        mem[Y.sha(src)] = dict(sha=Y.sha(src), line=1, src=src, out=out, keep=OUT[c] is None, by='a converter',
            when='2026-11-16 15:34')
    Y.write_memory(p.memory, mem, Y.analyze(NEW))
    return p


def new_entries(p):
    """sha -> entry of the memory's entries made by carry."""
    return {e['sha']: e for e in map(json.loads, p.memory.read_text().splitlines()) if Y.CARRIED.match(e['by'])}


def line_no(text):
    return NEW.split('\n').index(text) + 1


def test_find_split_forms():
    lines = Y.analyze(NEW)
    # a paragraph split in two, turned into a list, a list item into sub-bullets, a moved bold
    for case, first, n in (('para', 'The pairs are', 2), ('list', 'Three things', 4), ('sub', '- **Lead:**', 2),
            ('bold', '**The block', 2), ('lead', '**The lead', 2)):
        splits, why, _ = Y.find_split(lines, OLD[case])
        assert len(splits) == 1 and why == '', (case, why)
        assert len(splits[0].idx) == n and lines[splits[0].idx[0]].text.startswith(first)
    assert Y.find_split(lines, OLD['bold'])[0][0].moved == (0, len('**The block is large.'), 0)
    # a long bold lead whose closing ** moved to the end of its first sentence: the line's first bold span
    assert Y.find_split(lines, OLD['lead'])[0][0].moved == (0, len('**The lead is long.'), 0)
    assert Y.find_split(Y.analyze('A **b.** C. **D e.** F g.\n\nH.\n'), 'A **b.** C. **D e. F g.** H.')[0][0].moved[2] == 1
    # a change that is not a split, and kinds that are never split
    assert 'not a pure split' in Y.find_split(lines, OLD['changed'])[1]
    assert 'only paragraphs and list items' in Y.find_split(lines, '| a | b |')[1]
    # the same text with its bold dropped is not a pure split; nor is a paragraph that became one list item, or a list
    # item that became a paragraph
    assert '** changed' in Y.find_split(Y.analyze('A b. C d.\n\nE f.\n'), '**A b. C d. E f.**')[1]
    assert 'list marker' in Y.find_split(Y.analyze('- A b. C d.\n'), 'A b. C d.')[1]
    assert 'not a pure split' in Y.find_split(Y.analyze('A b. C d.\n'), '- A b. C d.')[1]
    # a sibling list item is not a sub-bullet
    assert not Y.find_split(Y.analyze('- A b.\n- C d.\n'), '- A b. C d.')[0]
    # a wholly bold paragraph whose closing ** moved, without a split
    one = Y.find_split(Y.analyze('**A b.** C d.\n'), '**A b. C d.**')[0]
    assert len(one) == 1 and one[0].moved == (0, len('**A b.'), 0)
    # a bold lead's closing ** moved earlier in one line; a closing ** may only move earlier, inside its span
    assert Y.find_split(Y.analyze('**A b.** C d. E f.\n'), '**A b. C d.** E f.')[0][0].moved == (0, len('**A b.'), 0)
    assert '** changed' in Y.find_split(Y.analyze('**A b. C d. E f.**\n'), '**A b. C d.** E f.')[1]


def test_carry_paragraph_list_sub_bullets_and_bold(tmp_path):
    p = carry_world(tmp_path, ('para', 'list', 'sub', 'bold', 'lead'))
    r = Y.carry(p)
    assert r['failures'] == [] and r['stale'] == 5 and len(r['carried']) == 5
    got = {e['src']: e for e in new_entries(p).values()}
    # a paragraph split in two: the converted sentences go with their piece, the unstarred unchanged piece needs none
    e = got['The pairs are {{data:a_dist ~ .2f}} apart. They move closer with APL off.']
    assert e['out'] == 'The pairs have a correlation of {{inv:data:a_dist ~ .2f}}. They become more correlated with APL off.'
    assert not e['keep'] and e['line'] == line_no(e['src'])
    assert e['by'] == f'carried over from {Y.sha(OLD["para"])} by yang_summary.py carry (from an entry by a converter, 2026-11-16 15:34)'
    assert re.fullmatch(r'\d{4}-\d\d-\d\d \d\d:\d\d', e['when']) and list(e) == list(Y.MEMORY_KEYS)
    assert 'The ρ is {{data:rho ~ .2f}}.' not in got
    # a paragraph turned into a list (KEEP): KEEP entries for the starred pieces only
    assert {s for s, e in got.items() if e['keep']} == {'- First, the pairs move apart.', '- Second, the noise floor rises.'}
    assert all(got[s]['out'] == s for s in ('- First, the pairs move apart.', '- Second, the noise floor rises.'))
    # a list item turned into sub-bullets: each piece keeps its own marker and indentation, a rewritten lit: kept
    assert got['- **Lead:** the distance is {{data:a_dist ~ .2f}}.']['out'] == '- **Lead:** the correlation is {{inv:data:a_dist ~ .2f}}.'
    assert got['  - It is far from {{data:b_dist ~ .2f}} and {{lit:0.88 # src}}.']['out'] == (
        '  - It is weakly correlated with {{inv:data:b_dist ~ .2f}} and {{lit:0.12 # Yang\'s version, rewritten from '
        '"0.88"; src}}.')
    # a moved bold: the converted line's closing ** moves to the end of the same sentence
    assert got['**The block is large.** It sits above the curve.']['out'] == '**The block is large.** It sits below the curve.'
    assert 'More text here, not about the block.' not in got
    assert got['**The lead is long.** It says the pairs move apart.']['out'] == '**The lead is long.** It says the pairs are decorrelated.'
    assert 'The rest is plain. It ends here.' not in got
    # the page shows the carried pieces; the old entries stay in the memory, stale, and status counts both
    lines = Y.analyze(NEW)
    out, info = Y.yang_lines(lines, Y.load_memory(p.memory), Y.starred(lines))
    assert out[line_no('The pairs are {{data:a_dist ~ .2f}} apart. They move closer with APL off.') - 1].startswith('The pairs have')
    assert not any(Y.YDIST in t for t in out[:line_no('More text here, not about the block.')])
    s = Y.status(p)
    assert s['carried'] == 7 and len(s['stale']) == 5 and len(s['stale_carried']) == 5   # the KEEP paragraph's too
    assert not (p.convert / 'carry_failures.md').exists()
    # carrying again changes nothing
    before = p.memory.read_text()
    r2 = Y.carry(p)
    assert r2['entries'] == [] and r2['failures'] == [] and p.memory.read_text() == before
    assert all(a in ('has this entry already', 'no entry needed') for c in r2['carried'] for _, a in c['actions'])


def test_carry_renumbers_placeholders(tmp_path):
    p = carry_world(tmp_path, ('renum',))
    lines = Y.analyze(NEW)
    piece = 'The distances run from {{data:a_dist ~ .2f}} to {{data:b_dist ~ .2f}} ({{data:count ~ .0f}} pairs).'
    # in the source the placeholders are 1-3 of the line and so in the piece, which is the whole second sentence;
    # the converted piece, its range swapped, is the answer ⟨2|inv⟩ ... ⟨1|inv⟩ ... ⟨3⟩ against the piece
    conv = 'The correlations run from {{inv:data:b_dist ~ .2f}} to {{inv:data:a_dist ~ .2f}} ({{data:count ~ .0f}} pairs).'
    assert Y.as_answer(conv, piece) == ('The correlations run from ⟨2|inv⟩ to ⟨1|inv⟩ (⟨3⟩ pairs).', [])
    # with a placeholder before it in the line, the numbers shift by one: ⟨2⟩-⟨4⟩ of the line are ⟨1⟩-⟨3⟩ of the piece
    line = 'Value {{data:rho ~ .2f}}. ' + piece
    c = Y.compact(line, Y.Filler())
    assert '⟨2|0.12|a_dist⟩' in c and Y.compact(piece, Y.Filler()).startswith('The distances run from ⟨1|0.12|a_dist⟩')
    assert Y.as_answer('{{inv:data:nope ~ .2f}}', piece)[1]
    r = Y.carry(p)
    assert r['failures'] == []
    got = {e['src']: e['out'] for e in new_entries(p).values()}
    assert got == {piece: conv}
    assert Y.Filler().line(conv, 't', tooltips=False) == 'The correlations run from 0.50 to 0.88 (92 pairs).'
    del lines


def test_carry_failures(tmp_path):
    p = carry_world(tmp_path, ('nomatch', 'inside', 'changed', 'list'))
    r = Y.carry(p)
    reasons = {x['src']: x['reason'] for x in r['failures']}
    assert set(reasons) == {OLD['nomatch'], OLD['inside'], OLD['changed']}
    # sentence boundaries that do not match, a split inside a sentence, a change that is not a split
    assert 'has 1 sentences where its source has 2' in reasons[OLD['nomatch']]
    assert 'inside a sentence' in reasons[OLD['inside']]
    assert 'not a pure split' in reasons[OLD['changed']]
    # none of them gets entries; the KEEP paragraph is carried (a KEEP split inside a sentence would be too)
    assert {e['src'] for e in new_entries(p).values()} == {'- First, the pairs move apart.', '- Second, the noise floor rises.'}
    md = (p.convert / 'carry_failures.md').read_text()
    assert f'## {Y.sha(OLD["nomatch"])} (recorded line 1)' in md
    assert f'@@ {line_no("The pairs are far apart.")}\nThe pairs are far apart.\n@@ {line_no("The noise floor is low.")}' in md
    assert f'converted:\n\n{OUT["changed"]}' in md
    # for a change that is not a split, the most similar line is shown
    assert f'@@ {line_no("The new wording of a line about the distance.")}\nThe new wording' in md
    # a converted piece that fails import's checks is not carried either
    (tmp_path / 'q').mkdir()
    q = carry_world(tmp_path / 'q', ('para',))
    mem = Y.load_memory(q.memory)
    mem[Y.sha(OLD['para'])]['out'] = OUT['para'].replace('{{inv:data:a_dist ~ .2f}}. They', 'nothing. They {{inv:data:a_dist ~ .2f}}', 1).replace(
        'The ρ is {{data:rho ~ .2f}}.', 'The ρ {{inv:data:a_dist ~ .2f}} is {{data:rho ~ .2f}}.')
    Y.write_memory(q.memory, mem, Y.analyze(NEW))
    r = Y.carry(q)
    assert len(r['failures']) == 1 and 'none of the placeholders of its piece' in r['failures'][0]['reason']
    assert new_entries(q) == {}


def test_carry_dry_run_and_command(tmp_path, monkeypatch, capsys):
    p = carry_world(tmp_path)
    before = p.memory.read_text()
    r = Y.carry(p, dry_run=True)
    assert r['entries'] and r['failures'] and p.memory.read_text() == before and not p.convert.exists()
    monkeypatch.setattr(Y, 'Paths', lambda: p)
    monkeypatch.setattr('sys.argv', ['yang_summary.py', 'carry', '--dry-run'])
    with pytest.raises(SystemExit) as e:
        Y.main()
    assert e.value.code == 1 and p.memory.read_text() == before and not p.convert.exists()
    out = capsys.readouterr().out
    assert 'WOULD CARRY' in out and 'NOT CARRIED' in out and 'would be written' in out
    monkeypatch.setattr('sys.argv', ['yang_summary.py', 'carry'])
    with pytest.raises(SystemExit):
        Y.main()
    assert p.memory.read_text() != before and (p.convert / 'carry_failures.md').exists()
    monkeypatch.setattr('sys.argv', ['yang_summary.py', 'status'])
    Y.main()
    assert 'carried over by carry' in capsys.readouterr().out


def test_carry_warns_on_questions_across_a_split(tmp_path):
    p = carry_world(tmp_path, ('para',))
    (p.data / 'questions.md').write_text(
        '# Q\n\n## q-within\nsentence_main: They move closer with APL off.\nquestion: Why?\n\nBecause.\n\n'
        '## q-across\nsentence_main: closer with APL off. The ρ is 0.7\nquestion: Across?\n\nYes.\n\n'
        '## q-yang\nsentence: more correlated with APL off. The ρ\nquestion: Yang side?\n\nYes.\n')
    r = Y.carry(p, dry_run=True)
    warned = [w for w in r['warnings'] if 'question' in w]
    assert len(warned) == 2 and any('q-across' in w for w in warned) and any('q-yang' in w for w in warned)
    assert not any('q-within' in w for w in warned)
