#!/usr/bin/env python3
"""Yang's version of the all-experiments summary page: every odor-odor distance given as a correlation
(correlation = 1 - distance), and every sentence Yang has asked about underlined, opening her question and its
answer in a floating window when clicked.

The page is built from the summary's template (reports/summary_template.md) and a translation memory of
converted template lines, reports/yang_summary/correlation_lines.jsonl (one JSON object per line: `sha`, the sha1
of the source template line and the lookup key; `line`; `src`, the source line; `out`, the converted template line;
`keep`; `by`; `when`). The rewording is done line by line by language-model agents and by hand, following
reports/yang_summary/conversion_rules.md; this script only moves lines out and back in and checks them.

- export [--all] [--chunk-kb N]: writes outputs/_artifacts/summary_yang/convert/chunk_NN.md, every non-blank
  template line as `@@ <line> <flag>` and the line in compact form (below). Flag `*`: to be answered (after the
  pre-pass it matches WATCH outside placeholders, or it is a row of a table with such a row); `=`: in the memory
  (shown converted, as plain text); none: context. Chunks are cut at `## ` / `### ` headings, numbered by position
  in the template; without --all a chunk with nothing starred is not written.
- import [FILE ...] [--by TEXT]: reads answers (`@@ <line>` and the converted compact line, or `@@ <line> KEEP`;
  default: every convert/chunk_*.out.md), checks each against its template line and writes the valid ones to the
  memory; failures go to convert/import_failures.md (exit status 1), warnings to convert/import_warnings.md.
- status: lines converted / KEEP / starred and not in the memory, stale memory entries, the missing lines by section.
- carry [--dry-run]: carries memory entries over a pure split of their template line (a paragraph or list item broken
  between sentences into paragraphs, a bullet list or sub-bullets; see `carry`): each piece gets an entry of its own,
  `by` "carried over from <old sha> by yang_summary.py carry"; what cannot be carried goes to convert/carry_failures.md
  (exit status 1). --dry-run writes nothing and prints what it would do.
- build [--skip-figure]: draws the background figure in correlations (results/summary_yang/ in the shared
  directory), builds Yang's lines (memory, else the pre-pass, unreviewed starred lines marked "distances"), inserts
  reports/yang_summary/intro.md and the list of Yang's questions before the table of contents, fills, marks the
  sentences of reports/yang_summary/questions.md and writes outputs/_artifacts/summary_yang/ (index.html,
  summary_yang.md with its .values.json sidecar, build_report.json). Nothing is written to the shared directory but
  the figure.

Placeholders: besides the template's `{{key ~ spec}}`, Yang's lines may use `{{inv:KEY ~ spec}}` (1 minus the value
of KEY, a distance shown as a correlation) and `{{neg:KEY ~ spec}}` (minus the value, a signed change of a distance
shown as a change of a correlation), computed from the unrounded value with the same spec (a `lit:` value without a
spec keeps its own number of decimals); their tooltips say so and give the main page's value.

Compact form of a line: each placeholder becomes `⟨k|VALUE|HINT⟩` (k its 1-based index in the line, VALUE as the
main page prints it, HINT what it is). An answer uses every k once, as the source token or `⟨k⟩` (kept), `⟨k|inv⟩`,
`⟨k|neg⟩`, or `⟨k|lit|new text⟩` (a `lit:` placeholder only).

questions.md: free text, then one entry per `## <id>` ([a-z0-9-]+), with `asked:`, `sentence:` (Yang's wording)
and/or `sentence_main:` (the main page's), `question:` header lines up to the first blank line, then the answer
(markdown, placeholders allowed) up to the next `## ` line. A sentence is found in the page's visible text with every
number a wildcard; failing that, in the main page's text of a line (the whole line, or table cell, is then marked).
Usage: ./venv/bin/python reports/yang_summary.py {export [--all] [--chunk-kb N] | import [FILE ...] [--by TEXT] |
       status | carry [--dry-run] | build [--skip-figure]}
"""
from __future__ import annotations

import argparse
import hashlib
import html as html_mod
import json
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from kcq import harness as H
from reports import fill_report as F
from reports import md_to_artifact as M
from reports import report_values as R

TEMPLATE = H.HERE / 'reports' / 'summary_template.md'
DATA_DIR = H.HERE / 'reports' / 'yang_summary'
OUT_DIR = H.HERE / 'outputs' / '_artifacts' / 'summary_yang'
MEMORY_NAME = 'correlation_lines.jsonl'
FIGURE = 'results/summary_yang/summary_background_correlation.png'   # relative to the shared directory
FIGURE_PREFIX = 'results/summary_yang/'   # a converted figure line may point only here
TITLE = "Model experiments (Yang's version)"

# lines to be converted: after the pre-pass, outside placeholders
WATCH = re.compile(r'distanc|apart|farther|closer|far from|close to|too close|too far|near-maximal|noise floor|inflat|'
    r'cosine|dissimilar|similar|separat|(above|below) the (curve|trend|monotone|fit|diagonal)|compress|decorrelat', re.I)
# phrases in which "distance" may stay in a converted line (no warning)
DIST_ALLOWED = ('total variation distance', 'distance of one recording', 'natural-source distance', 'chemical distance',
    'Pearson correlation distance', 'drawn as distance', 'drawn in distance', '1 − distance', 'one minus the distance',
    'distance = 1 − correlation')
# the deterministic pre-pass on lines with no memory entry (outside placeholders)
PREPASS = (('on odor-odor distances', 'on odor-odor correlations'),)
# characters the compact form and the page's markers use; none may occur in the template
RESERVED = ('⟦', '⟧', '⟨', '⟩', '⟪', '⟫', '§')

PH = F.PH
TOKEN = re.compile(r'⟨([^⟨⟩]*)⟩')
HEADING = re.compile(r'^(#{1,6}) (.*)$')
HEAD_ID = re.compile(r'\s*\{#([A-Za-z0-9_-]+)\}\s*$')
LIST = re.compile(r'^(\s*)(\d+\. |[-*] )')
FIG = M._FIG_RE
STRUCTURAL = re.compile(r'^(%%% hover|%%%|%%% html \S+|\[TOC\]|\?\?\?|---)\s*$')
SEP_CELL = re.compile(r':?-{2,}:?')
SPAN = re.compile(r'\{span=(\d+)\}\s*$')
LINK = re.compile(r'\]\(([^)\s]*)\)')
BUTTON = re.compile(r'!button\[[^\]]*\]\([^)]*\)')
QID = re.compile(r'^[a-z0-9-]+$')
YDIST = '⟪ydist⟫'
MARK = re.compile(r'⟪[^⟪⟫]*⟫')


class InputError(ValueError):
    """questions.md, intro.md or the template cannot be used as they are (the message says why)."""


def _shared() -> Path:
    """The shared directory (a function, so that tests can point it elsewhere)."""
    return H.shared_dir()


@dataclass(frozen=True)
class Paths:
    """Where one run reads and writes (the defaults are the repository's; tests use their own)."""
    template: Path = TEMPLATE
    data: Path = DATA_DIR
    out: Path = OUT_DIR

    @property
    def memory(self) -> Path:
        return self.data / MEMORY_NAME

    @property
    def convert(self) -> Path:
        return self.out / 'convert'


# ---------------------------------------------------------------- template lines

@dataclass
class Line:
    """One line of a template (or of intro.md) with what the conversion needs to know about it."""
    n: int
    text: str
    kind: str   # blank, structural, heading, hover, fold, figure, table, separator, list, para
    section: str = ''
    table: Optional[int] = None
    header_row: bool = False


def classify(text: str, in_hover: bool = False) -> str:
    """The kind of a line: blank, structural, hover (a `term: definition` line of a `%%% hover` block), heading,
    fold (`??? title`), figure, separator (a table's `|---|` row), table, list or para."""
    s = text.strip()
    if not s:
        return 'blank'
    if STRUCTURAL.match(s):
        return 'structural'
    if in_hover:
        return 'hover'
    if HEADING.match(s):
        return 'heading'
    if s.startswith('??? '):
        return 'fold'
    if FIG.match(s):
        return 'figure'
    if s.startswith('|'):
        return 'separator' if all(SEP_CELL.fullmatch(c.strip()) for c in s.strip('|').split('|')) else 'table'
    if LIST.match(text):
        return 'list'
    return 'para'


def analyze(text: str) -> List[Line]:
    """Every line of a template with its kind, section heading and table membership (index i is line i + 1)."""
    out: List[Line] = []
    in_hover, section, table, prev_table = False, '', -1, False
    for n, t in enumerate(text.split('\n'), 1):
        s = t.strip()
        if s == '%%% hover':
            kind, in_hover = 'structural', True
        elif in_hover and s == '%%%':
            kind, in_hover = 'structural', False
        else:
            kind = classify(t, in_hover)
        if kind == 'heading' and re.match(r'^#{1,3} ', s):
            section = s
        is_table = kind in ('table', 'separator')
        if is_table and not prev_table:
            table += 1
        prev_table = is_table
        out.append(Line(n, t, kind, section or s, table if is_table else None))
    # header rows: before a table's first separator (its first row if it has none), as md_to_artifact.py reads them
    by_table: Dict[int, List[Line]] = {}
    for ln in out:
        if ln.table is not None:
            by_table.setdefault(ln.table, []).append(ln)
    for rows in by_table.values():
        sep = next((k for k, r in enumerate(rows) if r.kind == 'separator'), None)
        for k, r in enumerate(rows):
            r.header_row = r.kind == 'table' and (k < sep if sep else k == 0)
    return out


def sha(line: str) -> str:
    return hashlib.sha1(line.encode('utf-8')).hexdigest()


def _outside(line: str, fn: Callable[[str], str]) -> str:
    """`line` with `fn` applied to the text between placeholders."""
    out, pos = [], 0
    for m in PH.finditer(line):
        out += [fn(line[pos:m.start()]), m.group(0)]
        pos = m.end()
    return ''.join(out + [fn(line[pos:])])


def prepass(line: str) -> str:
    """The deterministic pre-pass for a template line with no memory entry (outside placeholders only): "on
    odor-odor distances" becomes "on odor-odor correlations", nothing else."""
    def sub(t: str) -> str:
        for a, b in PREPASS:
            t = t.replace(a, b)
        return t
    return _outside(line, sub)


def watch_hit(line: str) -> bool:
    """Whether a line, after the pre-pass and outside placeholders, matches WATCH."""
    return bool(WATCH.search(PH.sub(' ', prepass(line))))


def starred(lines: Sequence[Line]) -> Set[int]:
    """Line numbers to be answered: WATCH matches (structural lines and separators never), and every row but the
    separators of a table in which any row matches."""
    hit = {ln.n for ln in lines if ln.kind not in ('blank', 'structural', 'separator') and watch_hit(ln.text)}
    tables = {ln.table for ln in lines if ln.n in hit and ln.table is not None}
    return hit | {ln.n for ln in lines if ln.table in tables and ln.kind == 'table'}


def sections(lines: Sequence[Line]) -> List[Tuple[int, int]]:
    """(first, last + 1) indices of the template's sections: a new one starts at every `## ` or `### ` heading."""
    starts = [0] + [k for k, ln in enumerate(lines) if k and re.match(r'^#{2,3} ', ln.text)]
    return list(zip(starts, starts[1:] + [len(lines)]))


# ---------------------------------------------------------------- values

def split_transform(key: str) -> Tuple[Optional[str], str]:
    """('inv' / 'neg' / None, the key it applies to)."""
    for t in ('inv', 'neg'):
        if key.startswith(t + ':'):
            return t, key[len(t) + 1:].strip()
    return None, key


def as_number(value) -> float:
    """A resolved value as a float (a literal's text parsed, Unicode minus allowed); ValueError if it is none."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{value!r} is not a number')
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip().replace('−', '-'))
        except ValueError:
            raise ValueError(f'{value!r} is not a number') from None
    raise ValueError(f'{value!r} is not a number')


def own_spec(text: str) -> str:
    """The format spec that prints a number as `text` does: its number of decimals, and a '+' if it has one."""
    m = re.fullmatch(r'([+\-−]?)\d*(?:\.(\d*))?', text.strip())
    if not m or not re.search(r'\d', text):
        return 'g'
    return ('+' if m.group(1) == '+' else '') + f'.{len(m.group(2) or "")}f'


def _plain(v):
    return v.item() if hasattr(v, 'item') else v


class Values:
    """Placeholder keys resolved once per run (report_values.resolve, looked up at call time)."""

    def __init__(self):
        self._cache: Dict[str, object] = {}

    def get(self, key: str) -> R.Resolved:
        key = key.strip()
        if key not in self._cache:
            try:
                self._cache[key] = R.resolve(key)
            except Exception as e:   # noqa: BLE001 (reported with the line)
                self._cache[key] = e
        v = self._cache[key]
        if isinstance(v, Exception):
            raise v
        return v


def _msg(e: Exception) -> str:
    return str(e.args[0]) if isinstance(e, KeyError) and e.args else str(e)


class Filler:
    """Fills `{{key ~ spec}}` placeholders, including `inv:` and `neg:`, as fill_report.fill does (tooltips as
    `⟦pv:value§tip⟧`), collecting a sidecar entry per value when asked and every unresolved key."""

    def __init__(self, values: Optional[Values] = None):
        self.values = values or Values()
        self.side: List[dict] = []
        self.errors: List[str] = []

    def value(self, body: str) -> Tuple[str, str, dict]:
        """(text, tooltip, sidecar entry) of one placeholder body; raises on a key that does not resolve or a
        transform of a non-numeric value."""
        key, spec = F.split_placeholder(body)
        transform, base = split_transform(key)
        res = self.values.get(base)
        main = R.format_value(res.value, spec)
        tip = R.tooltip(res)
        entry = dict(key=key, spec=spec, text=main, value=_plain(res.value), source=res.source)
        if transform is None:
            return main, tip, entry
        try:
            x = as_number(res.value)
        except ValueError as e:
            raise ValueError(f'{transform}: needs a number, the value is {e}') from None
        new = 1.0 - x if transform == 'inv' else -x
        text = R.format_value(new, spec or own_spec(main))
        lead = (f'shown as 1 − distance; the main page gives {main}: ' if transform == 'inv'
            else f'sign reversed (a change of a correlation); the main page gives {main}: ')
        entry.update(text=text, value=new, transform=transform, main_text=main, main_value=_plain(res.value))
        return text, lead + tip, entry

    def line(self, line: str, where: str, *, tooltips: bool, record: Optional[dict] = None) -> str:
        """`line` filled; an unresolved placeholder becomes '??' and is listed in `errors` (with `where`).

        Args:
            record: if given, a sidecar entry is kept for every value, with these fields added (line numbers).
        """
        def sub(m):
            try:
                text, tip, entry = self.value(m.group(1))
            except Exception as e:   # noqa: BLE001
                self.errors.append(f'{where}: {F.split_placeholder(m.group(1))[0]}: {_msg(e)}')
                return '??'
            if record is not None:
                self.side.append(dict(record, **entry))
            if not tooltips:
                return text
            tip = tip.replace('§', ' ').replace('⟧', ')').replace('|', '/').replace('"', "'")
            return f'⟦pv:{text}§{tip}⟧'
        return PH.sub(sub, line)

    def check(self) -> None:
        """Raises KeyError listing every unresolved placeholder so far (as fill_report.fill does)."""
        if self.errors:
            raise KeyError('unresolved placeholders:\n' + '\n'.join(self.errors))


# ---------------------------------------------------------------- compact form

def hint(key: str) -> str:
    """A short description of a placeholder key for the compact form."""
    key = key.strip()
    t, base = split_transform(key)
    if t:
        return f'{t}({hint(base)})'
    if key.startswith('data:'):
        return key[5:]
    if key.startswith('lit:'):
        return ''
    if key.startswith('meta:'):
        return key
    if key.startswith('diff:'):
        a, _, b = key[5:].partition(' -- ')
        return f'diff({hint(a)} − {hint(b)})'
    parts = key.split(':', 2)
    if len(parts) < 3:
        return key[-40:]
    row, col = parts[1].strip(), parts[2].split(',')[0].strip()
    return f"{row if len(row) <= 30 else '…' + row[-29:]}:{col}"


def _tok(s: str) -> str:
    return s.replace('|', '/').replace('⟩', '/').replace('⟨', '/')


def compact(line: str, filler: Filler) -> str:
    """`line` with each placeholder replaced by `⟨k|VALUE|HINT⟩` (VALUE '??' if it does not resolve)."""
    k = 0

    def sub(m):
        nonlocal k
        k += 1
        key, _ = F.split_placeholder(m.group(1))
        try:
            text = filler.value(m.group(1))[0]
        except Exception:   # noqa: BLE001
            text = '??'
        return f'⟨{k}|{_tok(text)}|{_tok(hint(key))}⟩'
    return PH.sub(sub, line)


def _norm_value(s: str) -> str:
    return re.sub(r'\s+', '', s).replace('−', '-')


def expand(answer: str, src: str, filler: Filler) -> Tuple[Optional[str], List[str], List[str]]:
    """The template line an answer in compact form stands for.

    Args:
        answer: the converted line, with tokens.
        src: the source template line.
        filler: resolves the source line's placeholders (for the numeric check and the kept tokens' values).

    Returns:
        (the template line or None, errors, warnings).
    """
    phs = list(PH.finditer(src))
    info = []
    for m in phs:
        key, spec = F.split_placeholder(m.group(1))
        try:
            text, _, entry = filler.value(m.group(1))
            value = entry['value']
        except Exception:   # noqa: BLE001
            text, value = None, None
        info.append(dict(m=m, key=key, spec=spec, text=text, value=value))
    errors, warnings, seen = [], [], Counter()

    def sub(m):
        inside = m.group(1)
        fields = inside.split('|')
        if not fields[0].strip().isdigit():
            errors.append(f'token ⟨{inside}⟩ does not start with a placeholder number')
            return ''
        k = int(fields[0])
        if not 1 <= k <= len(phs):
            errors.append(f'token ⟨{inside}⟩: the source line has {len(phs)} placeholder(s)')
            return ''
        seen[k] += 1
        p = info[k - 1]
        op = fields[1].strip() if len(fields) > 1 else ''
        spec = f' ~ {p["spec"]}' if p['spec'] else ''
        if op in ('inv', 'neg'):
            if p['text'] is None:
                errors.append(f'token ⟨{inside}⟩: {p["key"]} does not resolve')
                return ''
            try:
                x = as_number(p['value'])
            except ValueError:
                errors.append(f'token ⟨{inside}⟩: {op} needs a number, {p["key"]} is {p["text"]!r}')
                return ''
            if op == 'inv' and not -0.5 <= x <= 2.0:
                warnings.append(f'⟨{k}|inv⟩ on {p["text"]}, outside the range of a distance: check that it is one')
            return '{{' + f'{op}:{p["key"]}{spec}' + '}}'
        if op == 'lit':
            if not p['key'].startswith('lit:'):
                errors.append(f'token ⟨{inside}⟩: lit only replaces a hand-entered value (lit:), not {p["key"]}')
                return ''
            new = fields[2].strip() if len(fields) == 3 else ''
            if len(fields) != 3 or not new or '#' in new or '{' in new or '}' in new:
                errors.append(f'token ⟨{inside}⟩: write ⟨{k}|lit|new text⟩, the text non-empty and without | # {{ }}')
                return ''
            old, _, note = p['key'][4:].partition('#')
            return '{{' + f'lit:{new} # Yang\'s version, rewritten from "{old.strip()}"; {note.strip()}{spec}' + '}}'
        # a kept token must show the value the page gives (a value typed by hand is refused); meta: values (the build
        # time) change from one run to the next, so theirs is not compared
        if len(fields) > 1 and not p['key'].startswith('meta:') and (
                p['text'] is None or _norm_value(fields[1]) != _norm_value(_tok(p['text']))):
            errors.append(f'token ⟨{inside}⟩ changes the value (the main page gives {p["text"]}): write ⟨{k}|inv⟩ or '
                f'⟨{k}|neg⟩ to convert it, ⟨{k}⟩ to keep it, and never type a value a token carries')
            return ''
        return p['m'].group(0)

    out = TOKEN.sub(sub, answer)
    rest = TOKEN.sub(' ', answer)
    for k in range(1, len(phs) + 1):
        if seen[k] == 0:
            errors.append(f'token ⟨{k}⟩ is missing')
        elif seen[k] > 1:
            errors.append(f'token ⟨{k}⟩ is used {seen[k]} times')
    bad = [c for c in ('{{', '}}', '⟦', '⟧', '⟪', '⟫', '⟨', '⟩') if c in rest]
    if bad:
        errors.append(f'the text outside tokens contains {" ".join(bad)}')
    return (None if errors else out), errors, warnings


def _cells(text: str) -> List[str]:
    return [c.strip() for c in text.strip().strip('|').split('|')]


def _links(text: str, kind: str) -> Tuple[Counter, Counter]:
    t = PH.sub(' ', text)
    if kind == 'figure':
        t = FIG.match(t.strip()).group(1)
    return Counter(LINK.findall(t)), Counter(BUTTON.findall(t))


def structure_errors(src: Line, out: str) -> List[str]:
    """What a converted template line changed of its source's structure (empty list = nothing)."""
    if src.kind in ('structural', 'separator'):
        return [] if out == src.text else [f'a {src.kind} line cannot be converted']
    kind = classify(out, src.kind == 'hover')
    if kind != src.kind:
        return [f'the source is a {src.kind} line, the answer a {kind} line']
    errs = []
    s, o = src.text.strip(), out.strip()
    if kind == 'heading':
        hs, ho = HEADING.match(s), HEADING.match(o)
        if len(hs.group(1)) != len(ho.group(1)):
            errs.append('the heading level changed')
        ids, ido = HEAD_ID.search(s), HEAD_ID.search(o)
        if (ids and ids.group(1)) != (ido and ido.group(1)):
            errs.append(f'the heading id changed (it must stay {"{#" + ids.group(1) + "}" if ids else "absent"})')
    elif kind == 'list':
        if LIST.match(src.text).group(0) != LIST.match(out).group(0):
            errs.append('the list marker or its indentation changed')
    elif kind == 'table':
        cs, co = _cells(src.text), _cells(out)
        if len(cs) != len(co):
            errs.append(f'the row has {len(co)} cells instead of {len(cs)}')
        elif [SPAN.search(c) and SPAN.search(c).group(1) for c in cs] != [SPAN.search(c) and SPAN.search(c).group(1) for c in co]:
            errs.append('a {span=N} marker moved or changed')
    elif kind == 'hover':
        if ':' not in out or out.split(':', 1)[0] != src.text.split(':', 1)[0]:
            errs.append('the hover term before the colon changed')
    elif kind == 'figure':
        ps, po = FIG.match(s).group(2), FIG.match(o).group(2)
        if po != ps and not (po.startswith(FIGURE_PREFIX) and (_shared() / po).is_file()):
            errs.append(f'the figure path changed ({po}); only a figure that exists under {FIGURE_PREFIX} in the shared '
                'directory may replace it')
    ls, bs = _links(src.text, src.kind)
    lo, bo = _links(out, kind)
    if ls != lo:
        errs.append('the link targets changed')
    if bs != bo:
        errs.append('a !button[...](...) changed')
    bare = PH.sub(' ', out)
    if bare.count('**') % 2:
        errs.append('** occurs an odd number of times')
    if bare.count('`') % 2:
        errs.append('backticks occur an odd number of times')
    return errs


def _dash(s: str) -> str:
    return re.sub('[‐‑‒–—−]', '-', s).lower()


def distance_left(answer: str) -> bool:
    """Whether a converted line still says "distanc", outside tokens and the allowed phrases."""
    t = _dash(TOKEN.sub(' ', answer))
    for p in DIST_ALLOWED:
        t = t.replace(_dash(p), ' ')
    return 'distanc' in t


# ---------------------------------------------------------------- translation memory

MEMORY_KEYS = ('sha', 'line', 'src', 'out', 'keep', 'by', 'when')


def load_memory(path: Path) -> Dict[str, dict]:
    """sha -> entry of the translation memory (empty if the file does not exist)."""
    if not path.exists():
        return {}
    mem = {}
    for k, t in enumerate(path.read_text().split('\n'), 1):
        if t.strip():
            try:
                e = json.loads(t)
            except json.JSONDecodeError as ex:
                raise ValueError(f'{path}:{k}: not JSON ({ex})') from None
            mem[e['sha']] = e
    return mem


def write_memory(path: Path, mem: Dict[str, dict], lines: Sequence[Line]) -> None:
    """Writes the memory sorted by the entries' line numbers in the current template (stale ones by their recorded
    line), one object per line with a fixed key order."""
    first = {}
    for ln in lines:
        first.setdefault(sha(ln.text), ln.n)
    entries = sorted(mem.values(), key=lambda e: (first.get(e['sha'], e['line']), e['sha']))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(''.join(json.dumps({k: e[k] for k in MEMORY_KEYS}, ensure_ascii=False) + '\n' for e in entries))


# ---------------------------------------------------------------- export

def _chunk_name(k: int, total: int) -> str:
    return f'chunk_{k:0{max(2, len(str(total)))}d}'


def chunk_ranges(lines: Sequence[Line], compact_lines: Dict[int, str], star: Set[int], limit: int) -> List[Tuple[int, int]]:
    """(first, last + 1) indices of each chunk: whole sections, packed up to `limit` bytes of the template-only
    content (every flag as with an empty memory), so the numbering depends on the template alone."""
    out: List[Tuple[int, int]] = []
    cur, size = None, 0
    for a, b in sections(lines):
        n = sum(len(_entry(lines[k].n, '*' if lines[k].n in star else '', compact_lines[lines[k].n]).encode())
            for k in range(a, b) if lines[k].kind != 'blank')
        if cur is not None and size + n <= limit:
            cur, size = (cur[0], b), size + n
            continue
        if cur is not None:
            out.append(cur)
        cur, size = (a, b), n
    if cur is not None:
        out.append(cur)
    return out


def _entry(n: int, flag: str, text: str) -> str:
    return f'@@ {n}{" " + flag if flag else ""}\n{text}\n'


def check_reserved(text: str, what: str) -> None:
    bad = [c for c in RESERVED if c in text]
    assert not bad, f'{what} contains {" ".join(bad)}, which the compact form and the page markers use'


def export(paths: Paths = Paths(), *, all_chunks: bool = False, chunk_kb: float = 70) -> List[dict]:
    """Writes the chunks for conversion (see the module docstring) and returns one row per chunk."""
    template = paths.template.read_text()
    check_reserved(template, str(paths.template))
    lines = analyze(template)
    star = starred(lines)
    mem = load_memory(paths.memory)
    filler = Filler()
    comp = {ln.n: compact(prepass(ln.text), filler) for ln in lines if ln.kind != 'blank'}
    ranges = chunk_ranges(lines, comp, star, int(chunk_kb * 1024))
    paths.convert.mkdir(parents=True, exist_ok=True)
    for p in paths.convert.glob('chunk_*.md'):
        if not p.name.endswith('.out.md'):
            p.unlink()
    rel = _rel(paths.template)
    rows = []
    for k, (a, b) in enumerate(ranges, 1):
        name = _chunk_name(k, len(ranges))
        body, n_star = [], 0
        for ln in lines[a:b]:
            if ln.kind == 'blank':
                continue
            hit = mem.get(sha(ln.text))
            if hit:
                body.append(_entry(ln.n, '=', filler.line(hit['out'], '', tooltips=False)))
            elif ln.n in star:
                body.append(_entry(ln.n, '*', comp[ln.n])); n_star += 1
            else:
                body.append(_entry(ln.n, '', comp[ln.n]))
        head = (f'Chunk {k} of {len(ranges)}: lines {lines[a].n}-{lines[b - 1].n} of {rel}, for Yang\'s version of the '
            f'summary (every odor-odor distance given as a correlation).\n'
            f'The rules are in reports/yang_summary/conversion_rules.md. Answer every line flagged * in {name}.out.md '
            f'("@@ <line>" and the converted line on the next line, or "@@ <line> KEEP"); = is converted already, no '
            f'flag is context.\n\n')
        text = head + ''.join(body)
        written = n_star > 0 or all_chunks
        if written:
            (paths.convert / f'{name}.md').write_text(text)
        rows.append(dict(chunk=name, lines=f'{lines[a].n}-{lines[b - 1].n}', bytes=len(text.encode()), starred=n_star,
            written=written))
    return rows


def _rel(p: Path) -> str:
    try:
        return str(p.resolve().relative_to(H.HERE))
    except ValueError:
        return str(p)


# ---------------------------------------------------------------- import

ENTRY = re.compile(r'^@@\s+(\S+)(.*)$')


def parse_answers(text: str) -> List[dict]:
    """The entries of an answer file: dicts with `at` (the file's line number), `line` (template line number or
    None), `keep`, `answer` (the converted compact line) and `error` (a format problem, else None). Text before the
    first `@@` is ignored."""
    out, cur = [], None
    for k, t in enumerate(text.split('\n'), 1):
        m = ENTRY.match(t)
        if m:
            flags = m.group(2).split()
            keep = any(f.upper() == 'KEEP' for f in flags)
            cur = dict(at=k, line=int(m.group(1)) if m.group(1).isdigit() else None, keep=keep, body=[], error=None)
            if cur['line'] is None:
                cur['error'] = f'"@@ {m.group(1)}" is not a line number'
            elif any(f.upper() != 'KEEP' and f not in ('*', '=') for f in flags):
                cur['error'] = f'unexpected text after the line number: {m.group(2).strip()!r}'
            out.append(cur)
        elif cur is not None and t.strip():
            cur['body'].append(t.rstrip())
    for e in out:
        body = e.pop('body')
        e['answer'] = body[0] if body else ''
        if e['error']:
            continue
        if e['keep'] and body:
            e['error'] = 'KEEP takes no converted line'
        elif not e['keep'] and not body:
            e['error'] = 'no converted line (write "@@ <line> KEEP" for a line that needs no change)'
        elif len(body) > 1:
            e['error'] = f'the answer has {len(body)} lines; a converted line is exactly one line'
    return out


def import_answers(files: Sequence[Path], paths: Paths = Paths(), *, by: str = 'unspecified') -> Tuple[int, List[dict], List[str]]:
    """Checks answer files against the template and writes the valid entries to the memory.

    Returns:
        (number of entries imported, failures as dicts with file, at, line, reason, source, answer, warnings).
    """
    template = paths.template.read_text()
    lines = analyze(template)
    star = starred(lines)
    mem = load_memory(paths.memory)
    filler = Filler()
    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    comp_all = {}
    ranges = None
    failures, warnings, imported, answered = [], [], 0, {}

    def source_compact(n):
        if n not in comp_all:
            comp_all[n] = compact(prepass(lines[n - 1].text), filler)
        return comp_all[n]

    for f in files:
        f = Path(f)
        entries = parse_answers(f.read_text())
        for e in entries:
            n = e['line']
            ln = lines[n - 1] if n and 1 <= n <= len(lines) else None
            reasons = [e['error']] if e['error'] else []
            if not reasons and ln is None:
                reasons = [f'line {n} is not in the template ({len(lines)} lines)']
            elif not reasons and ln.kind == 'blank':
                reasons = [f'line {n} is blank in the template']
            out, warn = None, []
            if not reasons:
                if e['keep']:
                    out = ln.text
                else:
                    out, reasons, warn = expand(e['answer'], ln.text, filler)
                    if out is not None:
                        reasons = structure_errors(ln, out)
            if reasons:
                failures.append(dict(file=f.name, at=e['at'], line=n, reason='; '.join(reasons),
                    source=source_compact(n) if ln is not None and ln.kind != 'blank' else '', answer=e['answer'],
                    flag='*' if n in star else ''))
                continue
            if n in answered:
                warnings.append(f'{f.name}, line {e["at"]}: line {n} answered again ({answered[n]}); this answer is used')
            answered[n] = f'{f.name}, line {e["at"]}'
            warnings += [f'{f.name}, line {e["at"]} (template line {n}): {w}' for w in warn]
            if not e['keep'] and distance_left(e['answer']):
                warnings.append(f'{f.name}, line {e["at"]} (template line {n}): still says "distanc..." outside the '
                    'allowed phrases')
            mem[sha(ln.text)] = dict(sha=sha(ln.text), line=n, src=ln.text, out=out, keep=e['keep'] or out == ln.text,
                by=by, when=now)
            imported += 1
        m = re.match(r'^chunk_(\d+)\.out\.md$', f.name)
        sib = f.with_name(f.name[:-len('.out.md')] + '.md') if m else None
        head = re.search(r'lines (\d+)-(\d+) of', sib.read_text().split('\n', 1)[0]) if sib and sib.exists() else None
        if head:
            lo, hi = int(head.group(1)), int(head.group(2))
        elif m:
            if ranges is None:
                comp = {ln.n: compact(prepass(ln.text), filler) for ln in lines if ln.kind != 'blank'}
                ranges = chunk_ranges(lines, comp, star, 70 * 1024)
                comp_all.update(comp)
            k = int(m.group(1))
            lo, hi = (lines[ranges[k - 1][0]].n, lines[ranges[k - 1][1] - 1].n) if 1 <= k <= len(ranges) else (0, -1)
        else:
            nums = [e['line'] for e in entries if e['line']]
            lo, hi = (min(nums), max(nums)) if nums else (0, -1)
        missing = [n for n in sorted(star) if lo <= n <= hi and n not in answered and sha(lines[n - 1].text) not in mem
            and not any(x['line'] == n and x['file'] == f.name for x in failures)]
        if missing:
            warnings.append(f'{f.name}: {len(missing)} starred line(s) of lines {lo}-{hi} got no answer: '
                + ', '.join(map(str, missing)))
    write_memory(paths.memory, mem, lines)
    paths.convert.mkdir(parents=True, exist_ok=True)
    _write_report(paths.convert / 'import_failures.md', failures and _failures_md(failures, now))
    _write_report(paths.convert / 'import_warnings.md', warnings and (f'# Import warnings ({now})\n\n'
        + ''.join(f'- {w}\n' for w in warnings)))
    return imported, failures, warnings


def _write_report(path: Path, text) -> None:
    if text:
        path.write_text(text)
    elif path.exists():
        path.unlink()


def _failures_md(failures: List[dict], now: str) -> str:
    out = [f'# Answers that did not import ({now})\n\nEach entry: the answer file and line, the reason, the source '
        'line as the chunk showed it and the answer given. Answer again in the chunk format (rules: '
        'reports/yang_summary/conversion_rules.md).\n']
    for x in failures:
        out.append(f'\n## line {x["line"]} ({x["file"]}, line {x["at"]})\n\nreason: {x["reason"]}\n\nsource:\n\n'
            + (f'@@ {x["line"]}{" " + x["flag"] if x["flag"] else ""}\n{x["source"]}\n' if x['source'] else '(none)\n')
            + f'\nanswer given:\n\n{x["answer"] or "(none)"}\n')
    return ''.join(out)


# ---------------------------------------------------------------- status

def status(paths: Paths = Paths()) -> dict:
    """Counts of the conversion's progress (see the module docstring)."""
    lines = analyze(paths.template.read_text())
    star = starred(lines)
    mem = load_memory(paths.memory)
    shas = {sha(ln.text) for ln in lines}
    conv = [ln.n for ln in lines if ln.kind != 'blank' and sha(ln.text) in mem and not mem[sha(ln.text)]['keep']]
    keep = [ln.n for ln in lines if ln.kind != 'blank' and sha(ln.text) in mem and mem[sha(ln.text)]['keep']]
    missing = [ln for ln in lines if ln.n in star and sha(ln.text) not in mem]
    by_section: Dict[str, List[int]] = {}
    for ln in missing:
        by_section.setdefault(ln.section, []).append(ln.n)
    # entries made by `carry` (by their `by`), and the stale entries they were carried from
    origin = {m.group(1) for e in mem.values() for m in [CARRIED.match(e['by'])] if m}
    return dict(lines=sum(ln.kind != 'blank' for ln in lines), starred=len(star), converted=len(conv), keep=len(keep),
        missing=len(missing), stale=sorted(e['line'] for s, e in mem.items() if s not in shas), by_section=by_section,
        carried=sum(s in shas and bool(CARRIED.match(e['by'])) for s, e in mem.items()),
        stale_carried=sorted(e['line'] for s, e in mem.items() if s not in shas and s in origin))


# ---------------------------------------------------------------- carry

CARRIED = re.compile(r'^carried over from ([0-9a-f]{40}) by yang_summary\.py carry')
_ORIGIN = re.compile(r' \((from an entry by .*)\)$', re.S)
# where a sentence ends: . ! or ? (also : and ; when read loosely), with any closing brackets, quotes, * or _ after it,
# before whitespace
_END = {False: re.compile(r'[.!?][)\]"\'’”*_]*(?=\s)'), True: re.compile(r'[.!?:;][)\]"\'’”*_]*(?=\s)')}
_CODE = re.compile(r'`[^`]*`')
_LIT_REWRITE = re.compile(r'^lit:(.*?) # Yang\'s version, rewritten from "(.*?)";\s?(.*)$', re.S)


def _mask(text: str) -> str:
    """`text` with placeholders and inline code blanked out (same length), so that no sentence end is found inside."""
    for rx in (PH, _CODE):
        text = rx.sub(lambda m: 'x' * len(m.group(0)), text)
    return text


def sentence_spans(text: str, loose: bool = False) -> List[Tuple[int, int]]:
    """(start, end) of each sentence of a line's text (the text after its list marker), by the ends in `_END`; the
    whitespace between sentences belongs to none."""
    end = len(text.rstrip())
    spans, start = [], len(text) - len(text.lstrip())
    for m in _END[loose].finditer(_mask(text), 0, end):
        if m.end() > start:
            spans.append((start, m.end()))
        start = m.end()
        while start < end and text[start].isspace():
            start += 1
    if start < end:
        spans.append((start, end))
    return spans


def _parts(text: str, kind: str) -> Tuple[str, str]:
    """(prefix, content) of a paragraph or list line: its indentation and list marker, and the rest, stripped."""
    if kind == 'list':
        m = LIST.match(text)
        return m.group(0), text[m.end():].strip()
    return text[:len(text) - len(text.lstrip())], text.strip()


def _nb(s: str) -> str:
    return s.replace('**', '')


def _offsets(contents: Sequence[str], body: str) -> Optional[List[int]]:
    """Where each of `contents` starts in `body` if `body` is `contents` in order, separated by whitespace and by
    nothing else; else None."""
    pos, starts = 0, []
    for k, c in enumerate(contents):
        if not c or not body.startswith(c, pos):
            return None
        starts.append(pos)
        pos += len(c)
        if k < len(contents) - 1:
            if pos >= len(body) or not body[pos].isspace():
                return None
            while pos < len(body) and body[pos].isspace():
                pos += 1
    return starts if pos == len(body) else None


def _joins(contents: Sequence[str], body: str) -> bool:
    """Whether `body` is `contents` in order, separated by whitespace and by nothing else."""
    return _offsets(contents, body) is not None


def _bold_marks(text: str) -> List[int]:
    """Offsets of the `**` of a text, outside placeholders and inline code."""
    return [m.start() for m in re.finditer(r'\*\*', _mask(text))]


@dataclass
class Split:
    """A run of template lines that is a pure split of a stale entry's source line."""
    idx: List[int]                        # indices into the template's lines of the pieces (blank lines skipped)
    contents: List[str]                   # the pieces' contents (after the list marker), the moved ** taken out
    body: str                             # the source's content as the pieces give it back (its closing ** dropped if moved)
    # (piece, offset in its content, which bold span of the line: 0 the first) of a closing ** moved earlier, else None
    moved: Optional[Tuple[int, int, int]] = None


class _Source:
    """The kind, prefix and content of a stale entry's source line, and the rule for which lines may be its pieces."""

    def __init__(self, src: str):
        self.kind = classify(src)
        self.ok = self.kind in ('para', 'list')
        if self.ok:
            self.prefix, self.body = _parts(src, self.kind)
            self.indent = len(LIST.match(src).group(1)) if self.kind == 'list' else len(self.prefix)

    def fits(self, ln: Line, prev: Optional[Line]) -> bool:
        """Whether a template line may be the piece after `prev` (None: the first piece).

        Of a paragraph: a paragraph at its indentation, or a "- " / "* " item of a well-formed list, nested to any
        depth: the first item after a paragraph piece (or first of all) at the paragraph's indentation, every other
        item at most two spaces deeper than the item before it and not shallower than the paragraph. Of a list item:
        the item itself (same marker) first, then "- " / "* " items indented deeper than it.
        """
        if ln.kind not in ('para', 'list'):
            return False
        pre, _ = _parts(ln.text, ln.kind)
        if ln.kind == 'list':
            m = LIST.match(ln.text)
            indent, bullet = len(m.group(1)), m.group(2) in ('- ', '* ')
        if self.kind == 'para':
            if ln.kind == 'para':
                return pre == self.prefix
            if not bullet or not ln.text.startswith(self.prefix):
                return False
            if prev is None or prev.kind == 'para':
                return indent == self.indent
            return self.indent <= indent <= len(LIST.match(prev.text).group(1)) + 2
        if prev is None:
            return ln.kind == 'list' and pre == self.prefix
        return ln.kind == 'list' and indent > self.indent and bullet


def _runs(lines: Sequence[Line], s: _Source) -> List[List[int]]:
    """The runs of consecutive lines (blank lines between them allowed) that may be pieces of `s` and whose contents,
    with every ** dropped, give its content back (likewise): the candidates, checked exactly by `_split`."""
    target = _nb(s.body).strip()
    out = []
    for k, ln in enumerate(lines):
        if not s.fits(ln, None) or not target.startswith(_nb(_parts(ln.text, ln.kind)[1]).strip() or '\0'):
            continue
        run, pos = [], 0
        for i in range(k, len(lines)):
            ln2 = lines[i]
            if ln2.kind == 'blank':
                continue
            c = _nb(_parts(ln2.text, ln2.kind)[1]).strip()
            if not s.fits(ln2, lines[run[-1]] if run else None) or not c or not target.startswith(c, pos):
                break
            run.append(i)
            pos += len(c)
            if pos == len(target):
                out.append(run)
                break
            if not target[pos].isspace():
                break
            while pos < len(target) and target[pos].isspace():
                pos += 1
    return out


def _split(lines: Sequence[Line], run: List[int], s: _Source) -> Tuple[Optional[Split], str]:
    """(the Split, '') if a candidate run is a pure split of `s` exactly, else (None, why not)."""
    contents = [_parts(lines[i].text, lines[i].kind)[1] for i in run]
    odd = [lines[i].n for i, c in zip(run, contents) if PH.sub(' ', c).count('**') % 2]
    if odd:
        return None, f'** occurs an odd number of times in line(s) {", ".join(map(str, odd))}'
    if _joins(contents, s.body):
        if len(run) == 1:
            return None, f'line {lines[run[0]].n} differs from the source only in its list marker or indentation'
        return Split(list(run), contents, s.body), ''
    # a bold span (a whole bold paragraph, or a long bold lead) whose closing ** moved earlier, to the end of a word
    b = s.body
    marks = _bold_marks(b)
    for k in range(0, len(marks) - 1, 2):
        o, c_ = marks[k], marks[k + 1]
        bw = b[:c_] + b[c_ + 2:]
        for p, c in enumerate(contents):
            for j in _bold_marks(c):
                if j + 2 < len(c) and not c[j + 2].isspace():
                    continue
                cs = contents[:p] + [c[:j] + c[j + 2:]] + contents[p + 1:]
                starts = _offsets(cs, bw)
                if starts is not None and o < starts[p] + j < c_:
                    if len(run) == 1 and _parts(lines[run[0]].text, lines[run[0]].kind)[0] != s.prefix:
                        return None, f'line {lines[run[0]].n} changes the list marker or indentation'
                    return Split(list(run), cs, bw, (p, j, k // 2)), ''
    return None, ('the lines join back to the source only with its ** changed beyond the closing ** of a bold span '
        'moved earlier')


def find_split(lines: Sequence[Line], src: str) -> Tuple[List[Split], str, List[int]]:
    """The pure splits of a stale source line in the template.

    Returns:
        (the splits found (usually one), why none was found if none was, the line indices to show with that reason).
    """
    s = _Source(src)
    if not s.ok:
        return [], f'a {s.kind} line: only paragraphs and list items are split', []
    runs = _runs(lines, s)
    found, why, shown = [], '', []
    for run in runs:
        sp, reason = _split(lines, run, s)
        if sp:
            found.append(sp)
        elif not why:
            why, shown = reason, run
    if not found and not runs:
        why = 'no run of consecutive template lines joins back to the source line (not a pure split)'
    return found, why, shown


def _out_pieces(sp: Split, s: _Source, out: str) -> Tuple[Optional[List[str]], str]:
    """The contents of a converted line's pieces: its text cut at the sentence boundaries at which the source was
    split (with the moved ** at the end of the same sentence), or (None, why not)."""
    kind = classify(out)
    if kind != s.kind:
        return None, f'the converted line is a {kind} line, its source a {s.kind} line'
    _, ob = _parts(out, kind)
    if sp.moved is not None:
        marks = _bold_marks(ob)
        k = 2 * sp.moved[2] + 1
        if len(marks) <= k:
            return None, (f'the closing ** of bold span {sp.moved[2] + 1} of the source moved, but the converted line '
                'has no such span')
        ob = ob[:marks[k]] + ob[marks[k] + 2:]
    for loose in (False, True):
        sent = [sp.body[a:b] for a, b in sentence_spans(sp.body, loose)]
        per = [[c[a:b] for a, b in sentence_spans(c, loose)] for c in sp.contents]
        if sum(per, []) == sent:
            break
    else:
        return None, 'a split falls inside a sentence (by the sentence ends . ! ? : ;)'
    ends = '. ! ? : ;' if loose else '. ! ?'
    os_ = sentence_spans(ob, loose)
    if len(os_) != len(sent):
        return None, (f'the converted line has {len(os_)} sentences where its source has {len(sent)} (by the sentence '
            f'ends {ends}), so the pieces cannot be matched')
    first = [sum(len(x) for x in per[:p]) for p in range(len(per))]
    if sp.moved is not None:
        p, j, _ = sp.moved
        t = [b for _, b in sentence_spans(sp.contents[p], loose)]
        if j not in t:
            return None, 'the moved closing ** does not end a sentence'
        at = os_[first[p] + t.index(j)][1]
        ob = ob[:at] + '**' + ob[at:]
        os_ = [(a + 2 * (a >= at), b + 2 * (b >= at)) for a, b in os_]
    return [ob[os_[f][0]:os_[f + len(x) - 1][1]] for f, x in zip(first, per)], ''


def as_answer(out: str, src: str) -> Tuple[str, List[str]]:
    """The answer, in compact form against the template line `src`, that `expand` turns back into the converted line
    `out` (kept placeholders as ⟨k⟩, inv: and neg: as ⟨k|inv⟩ and ⟨k|neg⟩, a rewritten lit: as ⟨k|lit|text⟩, each k
    the placeholder's number in `src`), and the placeholders of `out` that match none of `src`'s."""
    phs = list(PH.finditer(src))
    keys = [F.split_placeholder(m.group(1)) for m in phs]
    used: Set[int] = set()
    errors: List[str] = []

    def sub(m):
        key, spec = F.split_placeholder(m.group(1))
        t, base = split_transform(key)
        lit = _LIT_REWRITE.match(key)
        for k, (pm, (pk, ps)) in enumerate(zip(phs, keys)):
            if k in used:
                continue
            if m.group(0) == pm.group(0):
                tok = f'⟨{k + 1}⟩'
            elif t and base == pk and spec == ps:
                tok = f'⟨{k + 1}|{t}⟩'
            elif lit and pk.startswith('lit:') and spec == ps and [x.strip() for x in pk[4:].partition('#')[::2]] == [
                    lit.group(2).strip(), lit.group(3).strip()]:
                tok = f'⟨{k + 1}|lit|{lit.group(1).strip()}⟩'
            else:
                continue
            used.add(k)
            return tok
        errors.append(f'{m.group(0)} is none of the placeholders of its piece of the template')
        return m.group(0)
    return PH.sub(sub, out), errors


def _questions_split(questions: Sequence[dict], before: Tuple[str, str], after: Sequence[Tuple[str, str]],
        filler: Filler, main: bool) -> List[str]:
    """The ids of the questions whose sentence is found in a line (`before`: text, kind) but in none of its pieces
    (`after`), as `locate` finds sentences (filled values, every number a wildcard)."""
    field_ = 'sentence_main' if main else 'sentence'

    def lt(text, kind):
        return LineText(filler.line(text, '', tooltips=False), kind)
    whole, pieces = lt(*before), [lt(*a) for a in after]
    out = []
    for q in questions:
        sent = norm_sentence(q[field_]) if q[field_] else ''
        if sent and whole.find(sent) and not any(p.find(sent) for p in pieces):
            out.append(q['id'])
    return out


def carry(paths: Paths = Paths(), *, dry_run: bool = False) -> dict:
    """Carries stale memory entries over a pure split of their source line.

    A pure split breaks a paragraph or a list item only between sentences and changes no word: a blank line between
    two sentences, a run of sentences turned into a bullet list ("- " before each; its items may nest, each at most two
    spaces deeper than the one before, the list starting at the paragraph's indentation), a list item broken into
    sub-bullets ("- " items indented deeper than it), and, where a whole paragraph or a long
    lead of a line was bold, its closing ** moved earlier, to the end of its first sentence (or another sentence of
    the span). For every memory entry whose `src` is no longer a line of the template, the runs of consecutive
    template lines (blank lines between them allowed) that are a pure split of it are found: their contents, list
    markers and indentation dropped, joined by whitespace, give `src`'s content back exactly, allowing for the moved
    **. Then:

    - a KEEP entry gives a KEEP entry to each piece that export would star (`starred`: WATCH after the pre-pass), and
      to an unstarred piece only if the pre-pass would change it (so that the page shows each piece as it showed the
      line); the other pieces need none;
    - a converted entry gives converted pieces when its `out` has the same number of sentences as `src` (sentence ends
      . ! ?, else also : and ;) and the pieces cover the same sentences of each: `out` is cut where `src` was cut, the
      pieces' own list markers and indentation put back, and the closing ** of the same bold span of `out` moved to the
      end of the same sentence. Each piece is then checked as `import` checks an answer: it is turned into an answer
      in compact form against its template line (`as_answer`; its placeholders renumbered per piece) and run through
      `expand` and `structure_errors`. A piece whose text equals its template line is treated as KEEP above.

    A new entry has the memory's usual fields: `line` the piece's line, `src` the piece, `by` "carried over from <old
    sha> by yang_summary.py carry (from an entry by <the old entry's by>, <its when>)", `when` now. A piece that has an
    entry already keeps it (a different wording is reported as a warning). The stale entries are left in place, since
    the memory's lines have fixed fields (`MEMORY_KEYS`) and `write_memory` would drop a new one; `status` counts the
    carried entries by their `by`. Running carry again carries nothing twice.

    What cannot be carried (a stale entry with no pure-split run, a split inside a sentence, sentences that do not
    match, a piece that fails the checks, two stale entries giving one piece two wordings) is listed in
    convert/carry_failures.md with the old sha, the reason and the template lines.

    Yang's questions (questions.md) need nothing: `locate` finds a question's sentence within one line (one table cell)
    of the page, and a pure split keeps every sentence whole, in a line of its own or with its neighbours, so a
    sentence is found in its piece as it was in the line, with its numbers; a moved ** is not part of the visible text
    it is matched on. Only a question whose quoted text runs across a split point would be lost; carry checks every
    carried line for this and warns. None of questions.md's 9 questions on 2026-11-17 quotes across a sentence end
    (. ! ?); three quote across a colon or a semicolon, which a split between sentences does not cut.

    Args:
        paths: the template, the data folder (memory, questions.md) and the output folder (convert/carry_failures.md).
        dry_run: write nothing (neither the memory nor the report); the result says what would be done.

    Returns:
        dict: `stale` (number of stale entries), `carried` (one dict per stale entry carried: sha, line, splits, the
        actions per piece), `entries` (the new memory entries), `failures` (sha, line, reason, src, out, the template
        lines shown) and `warnings`.
    """
    template = paths.template.read_text()
    lines = analyze(template)
    star = starred(lines)
    mem = load_memory(paths.memory)
    shas = {sha(ln.text) for ln in lines}
    filler = Filler()
    q_path = paths.data / 'questions.md'
    warnings: List[str] = []
    try:
        questions = parse_questions(q_path.read_text()) if q_path.exists() else []
    except InputError as e:
        questions = []
        warnings.append(f'questions.md not checked: {e}')
    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    new: Dict[str, dict] = {}
    carried, failures = [], []
    stale = sorted((e for s, e in mem.items() if s not in shas), key=lambda e: (e['line'], e['sha']))
    for e in stale:
        src, keep = e['src'], e['keep']
        s = _Source(src)
        splits, why, shown = find_split(lines, src)

        def fail(reason, idx=None):
            idx = shown if idx is None else idx
            failures.append(dict(sha=e['sha'], line=e['line'], reason=reason, src=src, out=None if keep else e['out'],
                lines=[(lines[i].n, lines[i].text) for i in (idx or _near(lines, src))]))
        if not splits:
            fail(why)
            continue
        origin = _ORIGIN.search(e['by']) if CARRIED.match(e['by']) else None
        by = (f'carried over from {e["sha"]} by yang_summary.py carry ('
            + (origin.group(1) if origin else f'from an entry by {e["by"]}, {e["when"]}') + ')')
        plan, reason = [], ''
        for sp in splits:
            if keep:
                outs = [lines[i].text for i in sp.idx]
            else:
                pieces, reason = _out_pieces(sp, s, e['out'])
                if pieces is None:
                    fail(reason, sp.idx)
                    break
                outs = []
                for i, c in zip(sp.idx, pieces):
                    ln = lines[i]
                    piece = _parts(ln.text, ln.kind)[0] + c
                    answer, errs = as_answer(piece, ln.text)
                    out2, errs2, warn = (None, [], []) if errs else expand(answer, ln.text, filler)
                    errs = errs or errs2 or structure_errors(ln, out2)
                    if errs:
                        reason = f'line {ln.n}: ' + '; '.join(errs) + f' (the piece: {piece})'
                        break
                    warnings += [f'{e["sha"][:8]}, line {ln.n}: {w}' for w in warn]
                    if out2 != piece:
                        warnings.append(f'{e["sha"][:8]}, line {ln.n}: the piece was rebuilt from its template line as '
                            f'{out2} (cut from the converted line: {piece})')
                    outs.append(out2)
                if reason:
                    fail(reason, sp.idx)
                    break
            # Yang's questions: still found in a piece wherever they were found in the line
            yang_src = src if keep else e['out']
            for main_side, before, after in ((True, src, [lines[i].text for i in sp.idx]), (False, yang_src, outs)):
                lost = _questions_split(questions, (before, s.kind),
                    [(t, lines[i].kind) for t, i in zip(after, sp.idx)], filler, main_side)
                warnings += [f'{e["sha"][:8]}: question {qid}: its sentence runs across the split at lines '
                    f'{lines[sp.idx[0]].n}-{lines[sp.idx[-1]].n} and will no longer be found' for qid in lost]
            plan.append((sp, outs))
        if reason:
            continue
        actions, pending, conflict = [], {}, None
        for sp, outs in plan:
            for i, o in zip(sp.idx, outs):
                ln = lines[i]
                h = sha(ln.text)
                if o == ln.text and ln.n not in star and prepass(ln.text) == ln.text:
                    actions.append((ln.n, 'no entry needed'))
                elif h in mem:
                    same = mem[h]['out'] == o
                    actions.append((ln.n, 'has this entry already' if same else 'has an entry of its own (kept)'))
                    if not same:
                        warnings.append(f'{e["sha"][:8]}: line {ln.n} has an entry of its own, which is kept; the '
                            f'carried wording would be: {o}')
                elif (h in new and new[h]['out'] != o) or (h in pending and pending[h]['out'] != o):
                    conflict = conflict or (ln.n, sp.idx)
                else:
                    pending.setdefault(h, dict(sha=h, line=ln.n, src=ln.text, out=o, keep=o == ln.text, by=by,
                        when=now))
                    actions.append((ln.n, 'KEEP' if o == ln.text else 'converted'))
        if conflict:
            fail(f'line {conflict[0]} is a piece of another stale entry too, which gives it another wording', conflict[1])
            continue
        for h, x in pending.items():
            new.setdefault(h, x)
        carried.append(dict(sha=e['sha'], line=e['line'], keep=keep, actions=actions,
            runs=[[lines[i].n for i in sp.idx] for sp, _ in plan]))
    if not dry_run:
        if new:
            write_memory(paths.memory, {**mem, **new}, lines)
        paths.convert.mkdir(parents=True, exist_ok=True)
        _write_report(paths.convert / 'carry_failures.md', failures and _carry_failures_md(failures, now))
    return dict(stale=len(stale), carried=carried, entries=list(new.values()), failures=failures, warnings=warnings)


def _words(text: str) -> Set[str]:
    return set(re.findall(r'\w+', PH.sub(' ', text).lower()))


def _near(lines: Sequence[Line], src: str) -> List[int]:
    """Indices of the template line most like `src` (by the overlap of their sets of words) and of the paragraph and
    list lines around it made mostly of `src`'s words: where a line that changed, or was split with a change, probably
    went."""
    a = _words(src)
    if not a:
        return []
    best, score = None, 0.0
    for k, ln in enumerate(lines):
        b = _words(ln.text) if ln.kind != 'blank' else set()
        if b and len(a & b) / len(a | b) > score:
            best, score = k, len(a & b) / len(a | b)
    if best is None or score < 0.3:
        return []
    out = [best]
    for step in (-1, 1):
        k, n = best + step, 0
        while 0 <= k < len(lines) and lines[k].kind in ('blank', 'para', 'list') and n < 8:
            if lines[k].kind != 'blank':
                b = _words(lines[k].text)
                if not b or len(a & b) / len(b) < 0.8:
                    break
                out.append(k)
                n += 1
            k += step
    return sorted(out)


def _carry_failures_md(failures: List[dict], now: str) -> str:
    out = [f'# Memory entries that were not carried over ({now})\n\nEach entry: the stale memory entry (its sha and the '
        'line it was recorded at), why it was not carried, its source line, its converted line (or KEEP) and the '
        'template lines it was matched against (or, if no run was found, the most similar line and those after it). '
        'Convert the lines anew (`export`, then `import`), or fix the split.\n']
    for x in failures:
        out.append(f'\n## {x["sha"]} (recorded line {x["line"]})\n\nreason: {x["reason"]}\n\nsource:\n\n{x["src"]}\n\n'
            f'converted:\n\n{x["out"] if x["out"] is not None else "KEEP"}\n\ntemplate lines:\n\n'
            + (''.join(f'@@ {n}\n{t}\n' for n, t in x['lines']) or '(none found)\n'))
    return ''.join(out)


# ---------------------------------------------------------------- questions

FIELDS = ('asked', 'sentence', 'sentence_main', 'question')


def parse_questions(text: str) -> List[dict]:
    """The entries of questions.md (see the module docstring), in file order.

    Raises:
        InputError: listing every problem (a bad or duplicate id, an unknown header field, no sentence, no question,
            an empty answer).
    """
    lines = text.split('\n')
    starts = [k for k, t in enumerate(lines) if t.startswith('## ')]
    entries, errors, ids = [], [], set()
    for a, b in zip(starts, starts[1:] + [len(lines)]):
        qid = lines[a][3:].strip()
        e = dict(id=qid, line=a + 1, **{f: '' for f in FIELDS})
        where = f'questions.md line {a + 1} ({qid})'
        if not QID.match(qid):
            errors.append(f'{where}: the id must be [a-z0-9-]+')
        if qid in ids:
            errors.append(f'{where}: duplicate id')
        ids.add(qid)
        k, seen = a + 1, set()
        while k < b and lines[k].strip():
            m = re.match(r'^([a-z_]+):\s?(.*)$', lines[k])
            if not m or m.group(1) not in FIELDS:
                errors.append(f'questions.md line {k + 1}: expected one of {", ".join(f + ":" for f in FIELDS)} up to '
                    'the first blank line')
            elif m.group(1) in seen:
                errors.append(f'questions.md line {k + 1}: {m.group(1)} given twice')
            else:
                seen.add(m.group(1))
                e[m.group(1)] = m.group(2).strip()
            k += 1
        e['answer'] = '\n'.join(lines[k:b]).strip('\n')
        if not e['sentence'] and not e['sentence_main']:
            errors.append(f'{where}: needs sentence: or sentence_main:')
        if not e['question']:
            errors.append(f'{where}: needs question:')
        if not e['answer'].strip():
            errors.append(f'{where}: the answer is empty')
        bad = [c for c in ('⟪', '⟫', '⟦', '⟧') if any(c in e[f] for f in FIELDS + ('answer',))]
        if bad:
            errors.append(f'{where}: contains {" ".join(bad)}')
        entries.append(e)
    if errors:
        raise InputError('questions.md:\n' + '\n'.join(errors))
    return entries


# ---------------------------------------------------------------- sentence matching

_VIS = re.compile(r'⟦pv:(.*?)§.*?⟧|⟪[^⟪⟫]*⟫|!button\[([^\]]*)\]\([^)\s]*\)|\[([^\]]+)\]\(([^)\s]+)\)|\*\*|`([^`]+)`'
    r'|<br\s*/?>', re.I)
_CHARMAP = {'−': '-', '–': '-', '—': '-', '‐': '-', '‑': '-', '‒': '-', '‘': "'", '’': "'", '“': '"', '”': '"',
    '×': 'x'}
_NUM = re.compile(r'[+-]?\d[\d.,]*%?')


def _emit(raw: str, lo: int, hi: int, chars: list, starts: list, ends: list, atoms: list) -> None:
    """The visible characters of raw[lo:hi] with their raw offsets; `atoms` gets the spans no marker may cut
    (a value with its tooltip, a link, a button, inline code, a ** pair)."""
    k, bold = lo, None
    for m in _VIS.finditer(raw, lo, hi):
        for p in range(k, m.start()):
            chars.append(raw[p]); starts.append(p); ends.append(p + 1)
        tok = m.group(0)
        if tok.startswith('⟦'):
            for p in range(m.start(1), m.end(1)):
                chars.append(raw[p]); starts.append(p); ends.append(p + 1)
            atoms.append((m.start(), m.end()))
        elif tok.startswith('⟪'):
            pass
        elif tok.startswith('!'):
            for p in range(m.start(2), m.end(2)):
                chars.append(raw[p]); starts.append(p); ends.append(p + 1)
            atoms.append((m.start(), m.end()))
        elif tok.startswith('['):
            _emit(raw, m.start(3), m.end(3), chars, starts, ends, atoms)
            atoms.append((m.start(), m.end()))
        elif tok == '**':
            if bold is None:
                bold = m.start()
            else:
                atoms.append((bold, m.end())); bold = None
        elif tok.startswith('`'):
            for p in range(m.start(5), m.end(5)):
                chars.append(raw[p]); starts.append(p); ends.append(p + 1)
            atoms.append((m.start(), m.end()))
        else:   # <br>
            chars.append(' '); starts.append(m.start()); ends.append(m.end())
        k = m.end()
    for p in range(k, hi):
        chars.append(raw[p]); starts.append(p); ends.append(p + 1)


@dataclass
class Seg:
    """The normalized visible text of one part of a line (a table cell, or the text after a list marker) and the map
    of each of its characters back to raw offsets."""
    lo: int
    hi: int
    norm: str
    starts: List[int]
    ends: List[int]
    atoms: List[Tuple[int, int]]


def visible(raw: str, lo: int = 0, hi: Optional[int] = None) -> Seg:
    """raw[lo:hi] normalized for matching: tooltip markup reduced to its value, `**`, backticks and link syntax to the
    visible text, page markers dropped, whitespace collapsed, dashes and the Unicode minus as '-', curly quotes
    straight, lower case, and every number (`[+-]?\\d[\\d.,]*%?`) replaced by '#'."""
    hi = len(raw) if hi is None else hi
    chars, starts, ends, atoms = [], [], [], []
    _emit(raw, lo, hi, chars, starts, ends, atoms)
    c2, s2, e2 = [], [], []
    for c, s, e in zip(chars, starts, ends):
        c = _CHARMAP.get(c, c)
        low = c.lower()
        c = low if len(low) == 1 else c
        if c.isspace():
            if not c2 or c2[-1] == ' ':
                continue
            c = ' '
        c2.append(c); s2.append(s); e2.append(e)
    while c2 and c2[-1] == ' ':
        c2.pop(); s2.pop(); e2.pop()
    text = ''.join(c2)
    c3, s3, e3, pos = [], [], [], 0
    for m in _NUM.finditer(text):
        c3 += c2[pos:m.start()]; s3 += s2[pos:m.start()]; e3 += e2[pos:m.start()]
        c3.append('#'); s3.append(s2[m.start()]); e3.append(e2[m.end() - 1])
        pos = m.end()
    c3 += c2[pos:]; s3 += s2[pos:]; e3 += e2[pos:]
    return Seg(lo, hi, ''.join(c3), s3, e3, atoms)


def norm_sentence(s: str) -> str:
    return visible(s.strip()).norm


def _ranges(text: str, kind: str) -> List[Tuple[int, int]]:
    """Raw (lo, hi) of the matchable parts of a line: its text after a list marker, `??? ` or the heading's `#`s (and
    before a heading id), a figure's caption, each table cell; none for other kinds."""
    if kind == 'para':
        return [(len(text) - len(text.lstrip()), len(text))]
    if kind == 'list':
        return [(LIST.match(text).end(), len(text))]
    if kind == 'fold':
        return [(text.index('??? ') + 4, len(text))]
    if kind == 'heading':
        m = HEADING.match(text)
        idm = HEAD_ID.search(text)
        return [(m.start(2), idm.start() if idm else len(text))]
    if kind == 'figure':
        m = FIG.match(text)
        return [(m.start(1), m.end(1))] if m else []
    if kind == 'table':
        s = len(text) - len(text.lstrip())
        e = len(text.rstrip())
        lo = s + 1 if text[s:s + 1] == '|' else s
        hi = e - 1 if e - 1 >= lo and text[e - 1] == '|' else e
        bars = [lo - 1] + [k for k in range(lo, hi) if text[k] == '|'] + [hi]
        return [(a + 1, b) for a, b in zip(bars, bars[1:])]
    return []


def _trim(text: str, lo: int, hi: int) -> Tuple[int, int]:
    while lo < hi and text[lo].isspace():
        lo += 1
    while hi > lo and text[hi - 1].isspace():
        hi -= 1
    if text.startswith(YDIST, lo):
        lo += len(YDIST)
    return lo, hi


def _snap(s: int, e: int, atoms: Sequence[Tuple[int, int]]) -> Tuple[int, int]:
    """Moves a span's start left and its end right out of any construct a marker may not cut."""
    changed = True
    while changed:
        changed = False
        for a, b in atoms:
            if a < s < b:
                s, changed = a, True
            if a < e < b:
                e, changed = b, True
    return s, e


class LineText:
    """The segments of one page line, normalized once for all questions."""

    def __init__(self, text: str, kind: str):
        self.text, self.kind = text, kind
        self.segs = [visible(text, a, b) for a, b in _ranges(text, kind)]

    def find(self, sent: str) -> Optional[Tuple[int, int, int]]:
        """(raw start, raw end, segment index) of the first match of a normalized sentence, or None."""
        if not sent:
            return None
        for ci, sg in enumerate(self.segs):
            j = sg.norm.find(sent)
            if j >= 0:
                s, e = _snap(sg.starts[j], sg.ends[j + len(sent) - 1], sg.atoms)
                return s, e, ci
        return None

    def whole(self, ci: int) -> Optional[Tuple[int, int]]:
        """The raw span of segment `ci` (its trimmed whole text, after a "distances" marker)."""
        r = _ranges(self.text, self.kind)
        if not r:
            return None
        a, b = r[min(ci, len(r) - 1)]
        a, b = _trim(self.text, a, b)
        return (a, b) if b > a else None


@dataclass
class Anchor:
    """Where a question is marked: a line of the document, a raw span in it, and how it was found."""
    doc: int
    start: int
    end: int
    how: str   # 'sentence' (Yang's wording), 'sentence_main' (the main page's, on Yang's line), 'whole'


def locate(questions: Sequence[dict], doc: Sequence[Tuple[str, str]], main: Callable[[int], Optional[Tuple[str, str]]],
        warnings: List[str]) -> Dict[str, Anchor]:
    """Finds each question's sentence in the page.

    Args:
        questions: parsed questions.md entries.
        doc: (filled line with tooltip markup, kind) of each line of Yang's page.
        main: for a document line index, the main page's filled text of its template line and its kind, if it differs
            from Yang's (None otherwise: nothing new to search).
        warnings: appended to.

    Returns:
        id -> Anchor for the questions found.
    """
    texts = [LineText(t, k) for t, k in doc]
    main_texts: Dict[int, Optional[LineText]] = {}
    out = {}

    def hits(sent, pool):
        return [(i, h) for i, lt in pool for h in [lt.find(sent)] if h]

    for q in questions:
        tries = []
        if q['sentence']:
            tries.append(('sentence', norm_sentence(q['sentence'])))
        if q['sentence_main']:
            tries.append(('sentence_main', norm_sentence(q['sentence_main'])))
            tries.append(('sentence_main', norm_sentence(prepass(q['sentence_main']))))
        found = None
        for how, sent in tries:
            hs = hits(sent, enumerate(texts))
            if hs:
                found = (how, hs)
                break
        if found:
            how, hs = found
            i, (s, e, _) = hs[0]
            out[q['id']] = Anchor(i, s, e, how)
            if len(hs) > 1:
                warnings.append(f'{q["id"]}: the sentence occurs in {len(hs)} lines; the first is marked')
            continue
        # the main page's text of each line that Yang's version changed
        for i in range(len(doc)):
            if i not in main_texts:
                m = main(i)
                main_texts[i] = LineText(*m) if m else None
        pool = [(i, lt) for i, lt in main_texts.items() if lt is not None]
        hs = []
        for _, sent in [t for t in tries if t[0] == 'sentence_main'][:1] + [t for t in tries if t[0] == 'sentence']:
            hs = hits(sent, pool)
            if hs:
                break
        if hs:
            i, (_, _, ci) = hs[0]
            w = texts[i].whole(ci)
            if w:
                out[q['id']] = Anchor(i, w[0], w[1], 'whole')
                warnings.append(f'{q["id"]}: sentence reworded since the question was asked: whole passage marked')
                if len(hs) > 1:
                    warnings.append(f'{q["id"]}: the sentence occurs in {len(hs)} lines of the main page; the first is used')
                continue
        warnings.append(f'{q["id"]}: the sentence is not on the page (unanchored)')
    return out


def wrap(doc_lines: List[str], anchors: Dict[str, Anchor], order: Sequence[str]) -> List[List[str]]:
    """Inserts `⟪yq:<n>⟫` ... `⟪/yq⟫` around each anchored span, in place; questions on overlapping spans of a line
    share one wrapper. Returns the ids of wrapper n (in questions.md order). The wrapper carries an index, not the
    ids, so that no hover term of the page can match inside it."""
    by_line: Dict[int, List[Tuple[int, int, List[str]]]] = {}
    for qid in order:
        if qid in anchors:
            a = anchors[qid]
            by_line.setdefault(a.doc, []).append((a.start, a.end, [qid]))
    groups: List[List[str]] = []
    for i in sorted(by_line):
        spans = sorted(by_line[i], key=lambda x: (x[0], -x[1]))
        merged = [list(spans[0])]
        for s, e, ids in spans[1:]:
            if s < merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e); merged[-1][2] += ids
            else:
                merged.append([s, e, ids])
        t = doc_lines[i]
        for s, e, ids in reversed(merged):
            groups.append(sorted(ids, key=list(order).index))
            g = len(groups) - 1
            t = t[:s] + f'⟪yq:{g}⟫' + t[s:e] + '⟪/yq⟫' + t[e:]
        doc_lines[i] = t
    return groups


# ---------------------------------------------------------------- Yang's lines

def yang_lines(lines: Sequence[Line], mem: Dict[str, dict], star: Set[int]) -> Tuple[List[str], dict]:
    """Template lines -> Yang's lines (with `⟪ydist⟫` before the text of unreviewed lines).

    A memory hit gives its `out`; a heading without an explicit id whose text changed gets the main page's anchor
    (`{#<slug of the original text>}`); a line with no hit gets the pre-pass, and if it is still starred (and is not a
    heading, a figure, a hover definition or structural) it is unreviewed: marked at the start of its text, after a
    list marker, inside the first cell of a table body row (not of a header row, where it would change how
    md_to_artifact.py lays the table out).

    Returns:
        (lines, info with converted / keep / unreviewed / marked line numbers).
    """
    out, info = [], dict(converted=[], keep=[], unreviewed=[], marked=[])
    for ln in lines:
        if ln.kind == 'blank':
            out.append(ln.text)
            continue
        hit = mem.get(sha(ln.text))
        if hit:
            t = hit['out']
            info['keep' if hit['keep'] else 'converted'].append(ln.n)
        else:
            t = prepass(ln.text)
            if ln.n in star and ln.kind not in ('heading', 'figure', 'structural', 'separator'):
                info['unreviewed'].append(ln.n)
                if ln.kind != 'hover' and not ln.header_row:
                    t = mark_unreviewed(t, ln.kind)
                    info['marked'].append(ln.n)
        if ln.kind == 'heading':
            t = keep_anchor(ln.text, t)
        out.append(t)
    return out, info


def keep_anchor(src: str, out: str) -> str:
    """A converted heading with the main page's anchor: `{#<slug of the source's text>}` appended when the source had
    no explicit id and the text changed."""
    if HEAD_ID.search(src) or HEAD_ID.search(out):
        return out
    hs, ho = HEADING.match(src.strip()), HEADING.match(out.strip())
    if not hs or not ho or hs.group(2) == ho.group(2):
        return out
    return out.rstrip() + ' {#' + M._slug(hs.group(2)) + '}'


def mark_unreviewed(t: str, kind: str) -> str:
    if kind == 'list':
        m = LIST.match(t)
        return t[:m.end()] + YDIST + t[m.end():]
    if kind == 'table':
        return re.sub(r'^(\s*\|\s*)', lambda m: m.group(1) + YDIST, t, count=1)
    if kind == 'fold':
        k = t.index('??? ') + 4
        return t[:k] + YDIST + t[k:]
    k = len(t) - len(t.lstrip())
    return t[:k] + YDIST + t[k:]


# ---------------------------------------------------------------- the page

PILL = ('<span class="ydist" title="This passage has not been converted yet: its distances are as on the main page '
    '(correlation = 1 − distance).">distances</span> ')

YQ_CSS = """
  span.ydist { display: inline-block; font-size: 0.72em; line-height: 1.45; padding: 0 6px; margin-right: 1px; border: 1px solid var(--line); border-radius: 999px; background: var(--soft); color: var(--muted); font-weight: 400; vertical-align: 0.12em; cursor: help; white-space: nowrap; }
  span.yq { text-decoration-line: underline; text-decoration-style: solid; text-decoration-color: var(--accent); text-decoration-thickness: 2px; text-underline-offset: 3px; text-decoration-skip-ink: none; cursor: pointer; border-radius: 2px; -webkit-box-decoration-break: clone; box-decoration-break: clone; }
  span.yq:hover, span.yq:focus-visible, span.yq.yq-open { background: rgba(42, 120, 214, 0.13); background: color-mix(in srgb, var(--accent) 15%, transparent); }
  span.yq:focus { outline: none; }
  span.yq:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
  .yq-pop { position: fixed; z-index: 900; display: none; flex-direction: column; width: min(560px, calc(100vw - 24px)); max-height: 70vh; box-sizing: border-box; background: var(--bg); color: var(--ink); border: 1px solid var(--line); border-radius: 10px; box-shadow: 0 12px 32px rgba(0, 0, 0, 0.24), 0 2px 6px rgba(0, 0, 0, 0.14); font-size: 0.94em; line-height: 1.5; }
  .yq-pop.open { display: flex; }
  .yq-head { display: flex; align-items: center; gap: 8px; padding: 6px 6px 6px 14px; border-bottom: 1px solid var(--line); background: var(--soft); border-radius: 10px 10px 0 0; cursor: move; user-select: none; -webkit-user-select: none; touch-action: none; }
  .yq-title { flex: 1; font-weight: 600; color: var(--muted); font-size: 0.9em; }
  .yq-close { font: inherit; font-size: 20px; line-height: 1; width: 30px; height: 30px; padding: 0; border: 1px solid transparent; border-radius: 6px; background: transparent; color: var(--muted); cursor: pointer; }
  .yq-close:hover, .yq-close:focus-visible { border-color: var(--line); color: var(--ink); background: var(--bg); }
  .yq-body { overflow-y: auto; padding: 2px 16px 14px; overscroll-behavior: contain; }
  .yq-q { font-weight: 600; margin: 12px 0 8px; }
  .yq-meta { color: var(--muted); font-weight: 400; font-size: 0.86em; margin-top: 2px; }
  .yq-a table { font-size: 0.8em; }
  .yq-a ol, .yq-a ul { padding-left: 22px; }
  hr.yq-sep { border: 0; border-top: 1px solid var(--line); margin: 14px 0 4px; }
"""

# The floating window of Yang's questions: one div.yq-pop, filled from the <template id="yqa-<id>"> of each id of the
# clicked span.yq; opened by a click, Enter or Space on a marked sentence, by a link to #yq-<id> or by the page's hash;
# closed by Esc, the x or a click outside it and outside any marked sentence. Figures inside an answer open in the
# page's zoom overlay (z-index 1000, above the window's 900).
YQ_JS = """<script>
(function () {
  if (!document.querySelector('span.yq')) return;
  var pop = null, body = null, title = null, closeBtn = null, current = null;
  function vw() { return document.documentElement.clientWidth || window.innerWidth; }
  function vh() { return document.documentElement.clientHeight || window.innerHeight; }
  function clamp(v, lo, hi) { return Math.max(lo, Math.min(v, hi)); }
  function make() {
    pop = document.createElement('div');
    pop.className = 'yq-pop';
    pop.setAttribute('role', 'dialog'); pop.setAttribute('aria-modal', 'false'); pop.setAttribute('aria-labelledby', 'yq-pop-title');
    pop.innerHTML = '<div class="yq-head"><span class="yq-title" id="yq-pop-title"></span>' +
      '<button type="button" class="yq-close" aria-label="close">\\u00d7</button></div><div class="yq-body"></div>';
    document.body.appendChild(pop);
    body = pop.querySelector('.yq-body'); title = pop.querySelector('.yq-title'); closeBtn = pop.querySelector('.yq-close');
    closeBtn.addEventListener('click', function (e) { e.stopPropagation(); close(); });
    body.addEventListener('click', function (e) { var z = e.target.closest('img.zoomable, button.zoombtn'); if (z) zoom(z); });
    body.addEventListener('keydown', function (e) {
      var z = e.target.closest && e.target.closest('img.zoomable');
      if (z && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); zoom(z); }
    });
    drag(pop.querySelector('.yq-head'));
  }
  function drag(head) {
    var sx = 0, sy = 0, ox = 0, oy = 0, pid = null;
    head.addEventListener('pointerdown', function (e) {
      if (e.button !== 0 || e.target.closest('.yq-close')) return;
      var r = pop.getBoundingClientRect();
      sx = e.clientX; sy = e.clientY; ox = r.left; oy = r.top; pid = e.pointerId;
      try { head.setPointerCapture(pid); } catch (x) {}
      e.preventDefault();
    });
    head.addEventListener('pointermove', function (e) {
      if (pid !== e.pointerId) return;
      pop.style.left = clamp(ox + e.clientX - sx, 0, Math.max(0, vw() - pop.offsetWidth)) + 'px';
      pop.style.top = clamp(oy + e.clientY - sy, 0, Math.max(0, vh() - pop.offsetHeight)) + 'px';
    });
    function end(e) { if (pid === e.pointerId) { try { head.releasePointerCapture(pid); } catch (x) {} pid = null; } }
    head.addEventListener('pointerup', end); head.addEventListener('pointercancel', end);
  }
  function fill(span) {
    var ids = (span.getAttribute('data-yq') || '').split(/\\s+/).filter(Boolean);
    body.innerHTML = '';
    ids.forEach(function (id, k) {
      var t = document.getElementById('yqa-' + id);
      if (!t) return;
      if (k) { var hr = document.createElement('hr'); hr.className = 'yq-sep'; body.appendChild(hr); }
      body.appendChild(t.content ? t.content.cloneNode(true) : document.createRange().createContextualFragment(t.innerHTML));
    });
    title.textContent = ids.length > 1 ? "Yang's questions (" + ids.length + ')' : "Yang's question";
    body.scrollTop = 0;
  }
  function place(span) {
    // below the sentence's last line box, else above its first; if neither side has room for the whole window,
    // the larger side, the window shortened to fit it (unless that leaves too little, then clamped over the text)
    var rects = span.getClientRects();
    var r = rects.length ? rects[rects.length - 1] : span.getBoundingClientRect();
    var r0 = rects.length ? rects[0] : r;
    pop.style.maxHeight = '';
    var w = pop.offsetWidth, h = pop.offsetHeight, W = vw(), Hh = vh();
    var below = Hh - r.bottom - 14, above = r0.top - 14, top;
    if (h <= below) top = r.bottom + 6;
    else if (h <= above) top = r0.top - 6 - h;
    else if (Math.max(below, above) >= 160) {
      pop.style.maxHeight = Math.max(below, above) + 'px'; h = pop.offsetHeight;
      top = below >= above ? r.bottom + 6 : r0.top - 6 - h;
    } else top = r.bottom + 6;
    pop.style.left = clamp(r.left, 12, Math.max(12, W - w - 12)) + 'px';
    pop.style.top = clamp(top, 8, Math.max(8, Hh - h - 8)) + 'px';
  }
  function focusQuiet(el) { try { el.focus({ preventScroll: true }); } catch (x) { el.focus(); } }
  function open(span) {
    if (!pop) make();
    if (current && current !== span) current.classList.remove('yq-open');
    current = span; span.classList.add('yq-open');
    fill(span);
    pop.style.left = '0px'; pop.style.top = '0px';
    pop.classList.add('open');
    place(span);
    focusQuiet(closeBtn);
  }
  function close() {
    if (!pop || !pop.classList.contains('open')) return;
    pop.classList.remove('open');
    var s = current; current = null;
    if (s) { s.classList.remove('yq-open'); focusQuiet(s); }
  }
  function overlay() {
    var ov = document.querySelector('.zoom-overlay');
    if (ov) return ov;
    ov = document.createElement('div');
    ov.className = 'zoom-overlay'; ov.setAttribute('role', 'dialog'); ov.setAttribute('aria-modal', 'true');
    ov.innerHTML = '<button class="zoom-close" aria-label="close">\\u00d7</button><img alt="">';
    document.body.appendChild(ov);
    function hide() {
      if (document.fullscreenElement === ov && document.exitFullscreen) document.exitFullscreen().catch(function () {});
      ov.classList.remove('open'); ov.querySelector('img').removeAttribute('src');
    }
    ov.addEventListener('click', hide);
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape' && ov.classList.contains('open')) hide(); });
    return ov;
  }
  function zoom(el) {
    var ov = overlay(), big = ov.querySelector('img');
    big.src = el.getAttribute('data-src') || el.src; big.alt = el.alt || el.title || '';
    ov.classList.add('open');
    if (ov.requestFullscreen) ov.requestFullscreen().catch(function () {});
  }
  function find(id) {
    var el = document.getElementById('yq-' + id);
    if (el && el.matches('span.yq')) return el;
    var all = document.querySelectorAll('span.yq');
    for (var k = 0; k < all.length; k++) {
      if ((' ' + all[k].getAttribute('data-yq') + ' ').indexOf(' ' + id + ' ') >= 0) return all[k];
    }
    return null;
  }
  function go(id) {
    var s = find(id);
    if (!s) return;
    for (var d = s.parentElement; d; d = d.parentElement) if (d.tagName === 'DETAILS') d.open = true;
    s.scrollIntoView({ block: 'center' });
    open(s);
    // the browser may still scroll to the fragment itself: centre the sentence again and follow it
    setTimeout(function () { if (current === s) { s.scrollIntoView({ block: 'center' }); place(s); } }, 0);
  }
  function fromHash() {
    var h = location.hash || '';
    try { h = decodeURIComponent(h); } catch (x) {}
    if (h.indexOf('#yq-') === 0) go(h.slice(4));
  }
  document.addEventListener('click', function (e) {
    var t = e.target;
    if (!t || !t.closest || t.closest('.zoom-overlay')) return;
    var a = t.closest('a[href^="#yq-"]');
    if (a) {
      e.preventDefault();
      var h = a.getAttribute('href');
      if (location.hash !== h) { try { history.pushState(null, '', h); } catch (x) {} }
      go(h.slice(4));
      return;
    }
    if (pop && pop.contains(t)) return;
    var s = t.closest('span.yq');
    if (s) {
      if (t.closest('a[href]')) return;
      var sel = window.getSelection && window.getSelection();
      if (sel && !sel.isCollapsed && s.contains(sel.anchorNode)) return;
      open(s);
      return;
    }
    close();
  });
  document.addEventListener('keydown', function (e) {
    var t = e.target;
    if (t && t.matches && t.matches('span.yq') && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); open(t); }
  });
  window.addEventListener('keydown', function (e) {
    if (e.key !== 'Escape' || !pop || !pop.classList.contains('open')) return;
    var ov = document.querySelector('.zoom-overlay.open');
    if (ov) return;
    close();
  }, true);
  window.addEventListener('resize', function () {
    if (!pop || !pop.classList.contains('open')) return;
    pop.style.left = clamp(pop.getBoundingClientRect().left, 0, Math.max(0, vw() - pop.offsetWidth)) + 'px';
    pop.style.top = clamp(pop.getBoundingClientRect().top, 0, Math.max(0, vh() - pop.offsetHeight)) + 'px';
  });
  window.addEventListener('hashchange', fromHash);
  window.addEventListener('load', function () {
    if (current && pop.classList.contains('open') && (location.hash || '').indexOf('#yq-') === 0) { current.scrollIntoView({ block: 'center' }); place(current); }
  });
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', fromHash); else fromHash();
})();
</script>"""

DEF_CSS = """
  .def-tip { position: fixed; z-index: 950; display: none; max-width: min(420px, calc(100vw - 24px)); box-sizing: border-box; padding: 8px 12px; background: var(--bg); color: var(--ink); border: 1px solid var(--line); border-left: 3px solid var(--accent); border-radius: 8px; box-shadow: 0 8px 24px rgba(0, 0, 0, 0.22), 0 1px 4px rgba(0, 0, 0, 0.12); font-size: 0.86em; font-weight: 400; line-height: 1.45; white-space: normal; text-align: left; }
  .def-tip.open { display: block; }
  .def-term { display: block; font-weight: 600; margin-bottom: 2px; }
  span.def.def-on { background: rgba(127, 127, 127, 0.16); border-radius: 2px; }
"""

# The definitions (span.def, the dotted terms; span.ydist, the "distances" pill): their text, in the title attribute,
# shows in a small box under the term, instead of the browser's own tooltip, which appears late or not at all (never on
# a touch screen). A mouse resting on a term shows it until the mouse leaves the term and the box; a click or a tap
# keeps it until Esc, another click or a click elsewhere. A click on a term inside one of Yang's sentences shows the
# definition and does not open her question; the rest of the sentence opens it. Terms inside her window work too (their
# titles are read when first pointed at, so the answers' templates need nothing done to them).
DEF_JS = """<script>
(function () {
  var SEL = 'span.def, span.ydist', tip = null, cur = null, pinned = false, timer = null;
  function vw() { return document.documentElement.clientWidth || window.innerWidth; }
  function vh() { return document.documentElement.clientHeight || window.innerHeight; }
  function text(el) {
    // the definition moves from title to data-def the first time the term is pointed at, so the browser's own
    // tooltip does not show as well
    var t = el.getAttribute('title');
    if (t !== null) { el.setAttribute('data-def', t); el.removeAttribute('title'); el.setAttribute('aria-description', t); }
    return el.getAttribute('data-def') || '';
  }
  function make() {
    tip = document.createElement('div'); tip.className = 'def-tip'; tip.setAttribute('role', 'tooltip');
    tip.addEventListener('pointerleave', function (e) { if (e.pointerType === 'mouse' && !pinned) later(hide); });
    tip.addEventListener('pointerenter', function () { clearTimeout(timer); });
    document.body.appendChild(tip);
  }
  function place() {
    if (!cur) return;
    var rects = cur.getClientRects(), r = rects.length ? rects[rects.length - 1] : cur.getBoundingClientRect();
    var r0 = rects.length ? rects[0] : r;
    if (r.bottom < 0 || r0.top > vh()) { hide(); return; }
    var w = tip.offsetWidth, h = tip.offsetHeight, top = r.bottom + 6;
    if (top + h > vh() - 8 && r0.top - 6 - h >= 8) top = r0.top - 6 - h;
    tip.style.left = Math.max(12, Math.min(r.left, vw() - w - 12)) + 'px';
    tip.style.top = Math.max(8, Math.min(top, vh() - h - 8)) + 'px';
  }
  function show(el, pin) {
    var t = text(el);
    if (!t) return;
    if (!tip) make();
    clearTimeout(timer);
    if (cur && cur !== el) cur.classList.remove('def-on');
    cur = el; pinned = !!pin; el.classList.add('def-on');
    tip.innerHTML = '';
    var term = document.createElement('span'); term.className = 'def-term'; term.textContent = el.textContent;
    tip.appendChild(term); tip.appendChild(document.createTextNode(t));
    tip.style.left = '0px'; tip.style.top = '0px'; tip.classList.add('open');
    place();
  }
  function hide() {
    clearTimeout(timer);
    if (tip) tip.classList.remove('open');
    if (cur) cur.classList.remove('def-on');
    cur = null; pinned = false;
  }
  function later(f) { clearTimeout(timer); timer = setTimeout(f, 200); }
  function term(t) { return t && t.closest ? t.closest(SEL) : null; }
  document.addEventListener('pointerover', function (e) {
    var el = term(e.target);
    if (!el) return;
    text(el);
    if (e.pointerType !== 'mouse' || pinned) return;
    clearTimeout(timer); timer = setTimeout(function () { show(el, false); }, 250);
  });
  document.addEventListener('pointerout', function (e) {
    var el = term(e.target);
    if (!el || e.pointerType !== 'mouse' || pinned) return;
    if (e.relatedTarget && (el.contains(e.relatedTarget) || (tip && tip.contains(e.relatedTarget)))) return;
    if (cur === el) later(hide); else clearTimeout(timer);
  });
  // capture phase: a click on a term shows its definition and goes no further (it neither opens nor closes Yang's
  // window); a click elsewhere closes a kept definition and carries on as usual
  document.addEventListener('click', function (e) {
    var el = term(e.target);
    if (el) {
      if (el.closest('a[href]')) return;
      var sel = window.getSelection && window.getSelection();
      if (sel && !sel.isCollapsed && el.contains(sel.anchorNode)) return;
      e.stopPropagation(); e.preventDefault();
      if (cur === el && pinned) hide(); else show(el, true);
      return;
    }
    if (tip && tip.contains(e.target)) { e.stopPropagation(); return; }
    if (cur) hide();
  }, true);
  // registered before Yang's window's own Esc handler (this script comes first in the page), so that Esc closes the
  // definition alone when one is shown
  window.addEventListener('keydown', function (e) {
    if (e.key === 'Escape' && cur && tip.classList.contains('open')) { hide(); e.stopImmediatePropagation(); }
  }, true);
  window.addEventListener('scroll', function () { if (cur) place(); }, true);
  window.addEventListener('resize', function () { if (cur) place(); });
})();
</script>"""

_VOID = {'br', 'img', 'hr', 'input', 'meta', 'link', 'wbr', 'source', 'col', 'area', 'base', 'embed', 'track'}
_TAG = re.compile(r'<(/?)([a-zA-Z][a-zA-Z0-9]*)\b[^>]*?(/?)>')
_BLOCKS = ('p', 'li', 'td', 'th', 'figcaption', 'summary', 'h1', 'h2', 'h3')


def balanced(fragment: str) -> bool:
    """Whether every element opened in an HTML fragment is closed in it, and nothing else is closed."""
    stack = []
    for m in _TAG.finditer(fragment):
        name = m.group(2).lower()
        if name in _VOID or m.group(3):
            continue
        if m.group(1):
            if not stack or stack[-1] != name:
                return False
            stack.pop()
        else:
            stack.append(name)
    return not stack


def _enclosing(page: str, s: int, e: int) -> Optional[Tuple[int, int, str]]:
    """(content start, content end, tag) of the innermost p / li / td / th / figcaption (or summary, heading) holding
    page[s:e]."""
    opens = list(re.finditer(r'<(' + '|'.join(_BLOCKS) + r')\b[^>]*>', page[:s]))
    for m in reversed(opens):
        name, depth = m.group(1), 1
        for t in re.finditer(r'<(/?)' + name + r'\b[^>]*>', page[m.end():]):
            depth += -1 if t.group(1) else 1
            if depth == 0:
                close = m.end() + t.start()
                if close >= e:
                    return m.end(), close, name
                break
    return None


def postprocess(page: str, groups: Sequence[Sequence[str]], answers: Dict[str, str]) -> Tuple[str, List[str]]:
    """The converted page with the pills, the marked sentences, the answers' templates, their CSS and script.

    Args:
        page: md_to_artifact.convert_markdown's output, with `⟪ydist⟫` and `⟪yq:<n>⟫` ... `⟪/yq⟫` markers.
        groups: the ids of wrapper n.
        answers: id -> the inner HTML of its `<template id="yqa-<id>">`.

    Returns:
        (page, warnings).
    """
    warnings = []
    # markers inside a tag (a figure's alt text) or in the table of contents (a heading's copy) are dropped
    page = re.sub(r'<[^<>]*>', lambda m: MARK.sub('', m.group(0)), page)
    page = re.sub(r'<nav class="toc">.*?</nav>', lambda m: MARK.sub('', m.group(0)), page, flags=re.S)
    page = page.replace(YDIST, PILL)
    for g, ids in enumerate(groups):
        start = f'⟪yq:{g}⟫'
        s = page.find(start)
        if s < 0:
            warnings.append(f'{", ".join(ids)}: the marked sentence was lost in the conversion')
            continue
        e = page.find('⟪/yq⟫', s)
        inner = page[s + len(start):e]
        tag = (''.join(f'<a id="yq-{i}"></a>' for i in ids[1:]) + f'<span class="yq" id="yq-{ids[0]}" '
            f'data-yq="{" ".join(ids)}" tabindex="0" role="button" aria-haspopup="dialog">')
        if balanced(inner):
            page = page[:s] + tag + inner + '</span>' + page[e + len('⟪/yq⟫'):]
            continue
        page = page[:s] + inner + page[e + len('⟪/yq⟫'):]
        enc = _enclosing(page, s, s + len(inner))
        if enc is None:
            raise RuntimeError(f'{", ".join(ids)}: the marked sentence crosses elements and has no enclosing block')
        a, b, name = enc
        page = page[:a] + tag + page[a:b] + '</span>' + page[b:]
        warnings.append(f'{", ".join(ids)}: the sentence crossed an HTML element; the whole {name} is marked')
    tpl = ''.join(f'<template id="yqa-{i}">{h}</template>\n' for i, h in answers.items())
    k = page.rfind('</main>')
    page = page[:k] + tpl + page[k:]
    page = page.replace('</style>', YQ_CSS + DEF_CSS + '</style>', 1)
    k = page.rfind('</body>')
    page = page[:k] + DEF_JS + '\n' + YQ_JS + '\n' + page[k:]
    assert '⟪' not in page and '⟫' not in page, 'a page marker survived the conversion'
    return page, warnings


def answer_html(q: dict, filler: Filler, figdir: Optional[Path]) -> str:
    """The inner HTML of a question's template: the question with "Yang, <asked>", and the answer filled and
    converted (its headings as bold lines, so that the page's footnote check does not read them)."""
    md = '\n'.join(re.sub(r'^#{1,6} (.*)$', r'**\1**', t) for t in q['answer'].split('\n'))
    md = '\n'.join(filler.line(t, f'questions.md {q["id"]}', tooltips=True, record=dict(line=None, template_line=None,
        question=q['id'])) for t in md.split('\n'))
    meta = 'Yang' + (f', {html_mod.escape(q["asked"])}' if q['asked'] else '')
    return (f'<div class="yq-q">{M.inline(q["question"])}<div class="yq-meta">{meta}</div></div>'
        f'<div class="yq-a">{M.convert_body(md, figdir, byline=False)}</div>')


def _cut(s: str, n: int = 140) -> str:
    s = s.replace('[', '(').replace(']', ')')
    if len(s) <= n:
        return s
    return s[:n].rsplit(' ', 1)[0].rstrip(',;:') + ' …'


def questions_block(questions: Sequence[dict], anchors: Dict[str, Anchor]) -> List[str]:
    """The collapsed list of Yang's questions inserted before the table of contents."""
    if not questions:
        return []
    out = [f"??? Yang's questions ({len(questions)})"]
    for q in questions:
        text = _cut(q['question'])
        out.append(f'- [{text}](#yq-{q["id"]})' if q['id'] in anchors else f'- {text} (the sentence is no longer on the page)')
    return out + ['???', '']


# ---------------------------------------------------------------- build

def background_figure(shared: Path) -> str:
    """Draws the background figure in correlations at <shared>/results/summary_yang/, unless the file is newer than
    the data-target registry and the scripts that draw it. Returns 'drawn' or 'kept'."""
    out = shared / FIGURE
    deps = [shared / 'results' / 'data_targets' / 'data_targets.json', Path(__file__),
        H.HERE / 'reports' / 'summary_background_figure.py']
    if out.exists() and all(out.stat().st_mtime > p.stat().st_mtime for p in deps if p.exists()):
        print(f'figure kept: {FIGURE} is newer than the data-target registry and the scripts')
        return 'kept'
    from reports import summary_background_figure as SB
    SB.draw_similarity(out)
    print(f'figure drawn: {FIGURE}')
    return 'drawn'


def build(paths: Paths = Paths(), *, skip_figure: bool = False, shared: Optional[Path] = None) -> Tuple[bool, dict]:
    """Builds Yang's page (steps in the module docstring).

    Args:
        paths: the template, the data folder (memory, intro.md, questions.md) and the output folder.
        skip_figure: do not draw the background figure in correlations.
        shared: the shared directory (figures are read from it; default H.shared_dir()).

    Returns:
        (whether md_to_artifact.lint found nothing, the build report).

    Raises:
        KeyError: a placeholder did not resolve (all are listed; nothing is written).
        InputError: questions.md or the template cannot be used.
    """
    shared = shared or _shared()
    figure = 'skipped' if skip_figure else background_figure(shared)
    template = paths.template.read_text()
    check_reserved(template, str(paths.template))
    lines = analyze(template)
    star = starred(lines)
    mem = load_memory(paths.memory)
    ylines, info = yang_lines(lines, mem, star)
    intro_path, q_path = paths.data / 'intro.md', paths.data / 'questions.md'
    intro_text = intro_path.read_text().strip('\n') if intro_path.exists() else ''
    check_reserved(intro_text, str(intro_path))
    intro = analyze(intro_text) if intro_text else []
    questions = parse_questions(q_path.read_text()) if q_path.exists() else []
    toc = next((k for k, ln in enumerate(lines) if ln.text.strip() == '[TOC]'), None)
    if toc is None:
        raise InputError('the template has no [TOC] line to put the introduction before')
    # the document: (Yang's line, kind, template line number or None, where for messages); the introduction goes
    # right before the [TOC] line, and the list of questions between the two once the sentences are found
    doc = [(t, ln.kind, ln.n, f'line {ln.n}') for t, ln in zip(ylines[:toc], lines[:toc])]
    doc += [(ln.text, ln.kind, None, f'intro.md line {ln.n}') for ln in intro] + ([('', 'blank', None, '')] if intro else [])
    k_toc = len(doc)
    doc += [(t, ln.kind, ln.n, f'line {ln.n}') for t, ln in zip(ylines[toc:], lines[toc:])]
    filler = Filler()
    tipped = [filler.line(t, where, tooltips=True) for t, _, _, where in doc]
    filler.check()
    main_filler = Filler(filler.values)

    def main_text(i):
        n = doc[i][2]
        if n is None or MARK.sub('', doc[i][0]) == lines[n - 1].text:
            return None
        return main_filler.line(lines[n - 1].text, '', tooltips=True), lines[n - 1].kind

    warnings: List[str] = []
    anchors = locate(questions, [(t, d[1]) for t, d in zip(tipped, doc)], main_text, warnings)
    # a marked heading without an explicit id keeps the main page's anchor (the marker would change its slug)
    for a in anchors.values():
        t = tipped[a.doc]
        if doc[a.doc][1] == 'heading' and not HEAD_ID.search(t):
            tipped[a.doc] = t.rstrip() + ' {#' + M._slug(HEADING.match(t).group(2)) + '}'
    groups = wrap(tipped, anchors, [q['id'] for q in questions])
    block = questions_block(questions, anchors)
    md = '\n'.join(tipped[:k_toc] + block + tipped[k_toc:])
    # the plain markdown (no tooltips, no markers) and its sidecar, one entry per value, line numbers of that file
    plain = []
    for i, (t, _, n, where) in enumerate(doc):
        at = i + 1 + (len(block) if i >= k_toc else 0)
        plain.append(filler.line(MARK.sub('', t), where, tooltips=False, record=dict(line=at, template_line=n)))
    plain = plain[:k_toc] + block + plain[k_toc:]
    answers = {q['id']: answer_html(q, filler, shared) for q in questions if q['id'] in anchors}
    filler.check()
    fill_lint = F.lint('\n'.join(MARK.sub('', d[0]) for d in doc))
    page = M.convert_markdown(md, shared, TITLE)
    page, pw = postprocess(page, groups, answers)
    warnings += pw
    probs = M.lint(page)
    paths.out.mkdir(parents=True, exist_ok=True)
    (paths.out / 'summary_yang_filled.md').write_text(md)
    (paths.out / 'index.html').write_text(page)
    plain_path = paths.out / 'summary_yang.md'
    plain_path.write_text('\n'.join(plain))
    side = filler.side
    n_inv = sum(e.get('transform') == 'inv' for e in side)
    n_neg = sum(e.get('transform') == 'neg' for e in side)
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=H.HERE, capture_output=True, text=True).stdout.strip()
    mem_text = paths.memory.read_text() if paths.memory.exists() else ''
    plain_path.with_suffix('.values.json').write_text(json.dumps(dict(template=_rel(paths.template),
        template_commit=commit, memory=_rel(paths.memory), memory_sha1=hashlib.sha1(mem_text.encode()).hexdigest(),
        generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), n_values=len(side), n_inv=n_inv,
        n_neg=n_neg, lint=fill_lint, values=side), indent=1, ensure_ascii=False, default=str) + '\n')
    shas = {sha(ln.text) for ln in lines}
    how = {qid: a.how for qid, a in anchors.items()}
    report = dict(generated=datetime.now().strftime('%Y-%m-%d %H:%M'), template=_rel(paths.template),
        figure=figure, lines=dict(non_blank=sum(ln.kind != 'blank' for ln in lines), starred=len(star),
            converted=len(info['converted']), keep=len(info['keep']), unreviewed=len(info['unreviewed']),
            unreviewed_lines=info['unreviewed'], marked_distances=len(info['marked']),
            stale_memory_entries=sorted(e['line'] for s, e in mem.items() if s not in shas)),
        questions=dict(total=len(questions), anchored=[q for q, h in how.items() if h != 'whole'],
            whole_passage=[q for q, h in how.items() if h == 'whole'],
            unanchored=[q['id'] for q in questions if q['id'] not in anchors],
            shared_spans=[g for g in groups if len(g) > 1], found_by=how),
        values=dict(filled=len(side), inv=n_inv, neg=n_neg), fill_lint=fill_lint, warnings=warnings, lint=probs,
        page=_rel(paths.out / 'index.html'), markdown=_rel(plain_path), bytes=len(page.encode()))
    (paths.out / 'build_report.json').write_text(json.dumps(report, indent=1, ensure_ascii=False) + '\n')
    return not probs, report


# ---------------------------------------------------------------- command line

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('export', help='write the chunks for conversion')
    p.add_argument('--all', action='store_true', help='also chunks with nothing to answer')
    p.add_argument('--chunk-kb', type=float, default=70, help='size limit of a chunk (default %(default)s)')
    p = sub.add_parser('import', help='import converted lines into the memory')
    p.add_argument('files', nargs='*', type=Path, help='default: every convert/chunk_*.out.md')
    p.add_argument('--by', default='unspecified', help='who converted them (recorded with each entry)')
    sub.add_parser('status', help='progress of the conversion')
    p = sub.add_parser('carry', help='carry memory entries over pure splits of their template lines')
    p.add_argument('--dry-run', action='store_true', help='write nothing, print what would be done')
    p = sub.add_parser('build', help="build Yang's page")
    p.add_argument('--skip-figure', action='store_true')
    args = ap.parse_args()
    paths = Paths()
    if args.cmd == 'export':
        rows = export(paths, all_chunks=args.all, chunk_kb=args.chunk_kb)
        print(f'{"chunk":<10} {"lines":<12} {"bytes":>8} {"starred":>8}')
        for r in rows:
            print(f'{r["chunk"]:<10} {r["lines"]:<12} {r["bytes"]:>8} {r["starred"]:>8}' + ('' if r['written'] else '  (not written)'))
        print(f'{sum(r["written"] for r in rows)} of {len(rows)} chunks written to {_rel(paths.convert)}; '
            f'{sum(r["starred"] for r in rows)} starred lines')
    elif args.cmd == 'import':
        files = args.files or sorted(paths.convert.glob('chunk_*.out.md'))
        if not files:
            sys.exit(f'no answer files in {_rel(paths.convert)}')
        n, failures, warnings = import_answers(files, paths, by=args.by)
        for x in failures:
            print(f'FAILED {x["file"]}, line {x["at"]} (template line {x["line"]}): {x["reason"]}')
        for w in warnings:
            print(f'WARNING {w}')
        print(f'{n} entries imported into {_rel(paths.memory)}; {len(failures)} failed'
            + (f' (convert/import_failures.md)' if failures else '') + f'; {len(warnings)} warnings')
        if failures:
            sys.exit(1)
    elif args.cmd == 'status':
        s = status(paths)
        rest = sorted((Counter(s['stale']) - Counter(s['stale_carried'])).elements())
        print(f'{s["lines"]} non-blank template lines, {s["starred"]} starred; memory: {s["converted"]} converted, '
            f'{s["keep"]} KEEP' + (f' ({s["carried"]} of them carried over by carry)' if s['carried'] else '')
            + f'; {s["missing"]} starred lines not in the memory; {len(s["stale"])} stale entries'
            + (f', {len(s["stale_carried"])} of them carried over' if s['stale_carried'] else '')
            + (f' ({"others: " if s["stale_carried"] else ""}recorded lines {", ".join(map(str, rest))})' if rest else ''))
        for sec, ns in s['by_section'].items():
            print(f'  {sec[:90]}: {", ".join(map(str, ns))}')
    elif args.cmd == 'carry':
        r = carry(paths, dry_run=args.dry_run)
        would = 'would be ' if args.dry_run else ''
        for c in r['carried']:
            print(f'{"WOULD CARRY" if args.dry_run else "CARRIED"} {c["sha"][:8]} (recorded line {c["line"]}, '
                f'{"KEEP" if c["keep"] else "converted"}) -> ' + ', '.join(f'line {n}: {a}' for n, a in c['actions']))
        for x in r['failures']:
            print(f'NOT CARRIED {x["sha"][:8]} (recorded line {x["line"]}): {x["reason"]}')
        for w in r['warnings']:
            print(f'WARNING {w}')
        print(f'{r["stale"]} stale entries: {len(r["carried"])} carried over, {len(r["failures"])} not; '
            f'{len(r["entries"])} new entries {would}written to {_rel(paths.memory)}'
            + (f'; failures {would}listed in {_rel(paths.convert / "carry_failures.md")}' if r['failures'] else ''))
        if r['failures']:
            sys.exit(1)
    else:
        try:
            ok, rep = build(paths, skip_figure=args.skip_figure)
        except (KeyError, InputError) as e:
            sys.exit(_msg(e))
        ln, q = rep['lines'], rep['questions']
        print(f'lines: {ln["converted"]} converted, {ln["keep"]} KEEP, {ln["unreviewed"]} unreviewed '
            f'({ln["marked_distances"]} tagged "distances"); {len(ln["stale_memory_entries"])} stale memory entries')
        print(f'questions: {q["total"]}: {len(q["anchored"])} anchored, {len(q["whole_passage"])} whole passage, '
            f'{len(q["unanchored"])} unanchored' + (f' ({", ".join(q["unanchored"])})' if q['unanchored'] else ''))
        print(f'values: {rep["values"]["filled"]} filled ({rep["values"]["inv"]} inv, {rep["values"]["neg"]} neg); '
            f'figure {rep["figure"]}')
        for w in rep['warnings']:
            print('WARNING', w)
        for p_ in rep['fill_lint']:
            print('LINT', p_)
        print(f'page -> {rep["page"]} ({rep["bytes"] // 1024} KB); markdown -> {rep["markdown"]}; report -> '
            f'{_rel(paths.out / "build_report.json")}')
        if not ok:
            print('LINT PROBLEMS:\n' + '\n'.join(' - ' + p_ for p_ in rep['lint']))
            sys.exit(1)
        print('lint clean')


if __name__ == '__main__':
    main()
