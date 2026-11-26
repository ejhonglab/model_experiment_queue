#!/usr/bin/env python3
"""Resolve report placeholders to values with recorded provenance.

A placeholder in report markdown is `{{key ~ spec}}` (spec: a Python format spec such as `.2f`
or `+.3f`, or `w` for a count in words, "two" for 2; negative numbers are written with a Unicode
minus; `|` is not used because it separates table cells). Keys:

- `data:<name>` or `data:<name>.lo` / `.hi`: an entry of the data-target registry
  (results/data_targets/data_targets.json in the shared directory, written by
  data_targets.py: value, optional 90 % interval, definition, function, inputs, commit).
- `<EXP>[/<table.csv>]:<row>:<column>`: one cell of an experiment's results table
  (results/<EXP>_*/<table>.csv, or results/<EXP>/ when <EXP> names a folder; without a table
  name every csv with a `tag` or `variant` column is searched, and a key matching more than one
  table is an error). `<row>` is a tag (matched against `tag`, else `variant`), or
  `col=value[;col=value]` to select by other columns (`col!=value` excludes rows; `*` selects every row, for
  `min:` / `max:`). `<column>` may list
  alternatives separated by commas; the first one the table has is used.
- `min:<table key>` / `max:<table key>`: the smallest / largest value of the column over every row the
  table key selects (one table), e.g. the lowest wiring's value without naming the wiring; the source
  records the row it came from.
- `lit:<text> # <source>`: a hand-entered value, allowed only with its source after `#`.
- `diff:<key> -- <key>`: the first key's value minus the second's (each any key above, e.g. a
  model's APL-on minus its APL-off statistic), computed from the unrounded values; both sources
  are recorded.
- `meta:generated` (the time the page is built) and `meta:last_experiment` (the latest
  completion time in the experiment index of experiment_priorities.md), as YYYY-MM-DD HH:MM.
- `json:<file>:<path>`: a value in a JSON file of the shared directory (`<file>` relative to it, e.g.
  results/E134_population_baseline/lock.json), `<path>` its keys joined by `.` (a list index is a number), or
  `len(<path>)` for the length of a list or dict there.
- `const:<file>:<NAME>[.<index>]`: a module-level constant of this repository's code (`<file>` relative to the
  repository, e.g. kcq/apl_effect.py), read from the file without importing it: the right-hand side of the
  file's one `NAME = ...` assignment, a literal or arithmetic on literals; `.<index>` picks an item of a tuple or
  list. The source records the file's sha256, its last commit and whether it has uncommitted changes.

Every resolution returns the value and a source record (file, row, column, sha256 of the
file, the commit of the experiment's scripts recorded in the folder's provenance if any), which
fill_report.py writes to the page's sidecar. That commit is `queue_commit` (this repository) in
provenance written since the scripts moved out of al_analysis (2026-11-11), and
`al_analysis_commit` (an al_analysis commit, where they lived) in older provenance.
"""
from __future__ import annotations

import glob
import hashlib
import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from kcq import harness as H

RESULTS = H.shared_dir() / 'results'
DATA_TARGETS = RESULTS / 'data_targets' / 'data_targets.json'


@dataclass
class Resolved:
    key: str
    value: Any
    source: Dict[str, Any] = field(default_factory=dict)


@lru_cache(maxsize=None)
def _sha256(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def _registry() -> Dict[str, Any]:
    return json.loads(DATA_TARGETS.read_text())


@lru_cache(maxsize=None)
def _table(path: str) -> pd.DataFrame:
    """A results table; leading lines that start with '#' (a caveat line, e.g. E101's "# Yang has not checked ...")
    are skipped, and '#' anywhere else is read as data."""
    n = 0
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if not line.startswith('#'):
                break
            n += 1
    return pd.read_csv(path, skiprows=n)


def _folders(exp: str) -> List[Path]:
    if (RESULTS / exp).is_dir():
        return [RESULTS / exp]
    return [Path(p) for p in sorted(glob.glob(str(RESULTS / f'{exp}_*'))) if Path(p).is_dir()]


@lru_cache(maxsize=None)
def _commit(folder: str) -> Optional[str]:
    for p in sorted(Path(folder).glob('provenance*.json')):
        try:
            d = json.loads(p.read_text())
        except Exception:
            continue
        # queue_commit first: newer provenance also has al_analysis_commit, the model code's commit
        for k in ('queue_commit', 'al_analysis_commit', 'commit'):
            if isinstance(d.get(k), str) and d[k]:
                return d[k]
        mp = d.get('model_provenance') or {}
        if isinstance(mp, dict):
            for repo in ('queue', 'al_analysis'):
                info = mp.get(repo) if isinstance(mp.get(repo), dict) else {}
                c = info.get('git_hash') or info.get('commit')
                if c:
                    return c
    return None


def _select(df: pd.DataFrame, row: str) -> pd.DataFrame:
    if row == '*':
        return df
    if '=' in row:
        m = pd.Series(True, index=df.index)
        for cond in row.split(';'):
            col, val = cond.split('=', 1)
            neg = col.endswith('!')
            col, val = col.rstrip('!').strip(), val.strip()
            if col not in df.columns:
                return df.iloc[0:0]
            m &= (df[col].astype(str) != val) if neg else (df[col].astype(str) == val)
        return df[m]
    for col in ('tag', 'variant'):
        if col in df.columns:
            return df[df[col].astype(str) == row]
    return df[df.iloc[:, 0].astype(str) == row]


def _resolve_table(key: str) -> Resolved:
    head, row, cols = key.split(':', 2)
    exp, _, table = head.partition('/')
    folders = _folders(exp)
    if not folders:
        raise KeyError(f'{key}: no results folder for {exp!r}')
    paths = [f / table for f in folders] if table else [Path(p) for f in folders for p in sorted(glob.glob(str(f / '*.csv')))]
    candidates = [c.strip() for c in cols.split(',')]
    hits = []
    for p in paths:
        if not p.exists():
            continue
        df = _table(str(p))
        if not table and not ({'tag', 'variant'} & set(df.columns)) and '=' not in row:
            continue
        col = next((c for c in candidates if c in df.columns), None)
        if col is None:
            continue
        sel = _select(df, row)
        if len(sel) == 0:
            continue
        if len(sel) > 1:
            raise KeyError(f'{key}: {len(sel)} rows match in {p}')
        hits.append((p, col, sel.iloc[0][col]))
    if not hits:
        raise KeyError(f'{key}: no table row / column matches')
    if len(hits) > 1:
        vals = {round(float(v), 9) if isinstance(v, (int, float, np.floating)) else v for _, _, v in hits}
        if len(vals) > 1:
            raise KeyError(f'{key}: matches in {len(hits)} tables with different values: {[str(h[0]) for h in hits]}')
    p, col, v = hits[0]
    v = v.item() if hasattr(v, 'item') else v
    return Resolved(key, v, dict(kind='table', file=str(p.relative_to(H.shared_dir())), row=row, column=col,
        sha256=_sha256(str(p)), commit=_commit(str(p.parent))))


def _resolve_extreme(key: str) -> Resolved:
    """`min:<table key>` / `max:<table key>` (module docstring)."""
    how, inner = key[:3], key[4:]
    head, row, cols = inner.split(':', 2)
    exp, _, table = head.partition('/')
    if not table:
        raise KeyError(f'{key}: name the table (<EXP>/<table.csv>)')
    paths = [f / table for f in _folders(exp) if (f / table).exists()]
    if len(paths) != 1:
        raise KeyError(f'{key}: {len(paths)} tables {table} for {exp!r}')
    df = _table(str(paths[0]))
    col = next((c.strip() for c in cols.split(',') if c.strip() in df.columns), None)
    sel = _select(df, row) if col else df.iloc[0:0]
    vals = pd.to_numeric(sel[col], errors='coerce').dropna() if col else pd.Series(dtype=float)
    if vals.empty:
        raise KeyError(f'{key}: no numeric value selected')
    i = vals.idxmin() if how == 'min' else vals.idxmax()
    v = vals[i]
    v = v.item() if hasattr(v, 'item') else v
    picked = {c: str(sel.loc[i, c]) for c in sel.columns if c in ('tag', 'variant', 'wiring', 'arm', 'group', 'fly')}
    return Resolved(key, v, dict(kind='table', file=str(paths[0].relative_to(H.shared_dir())), row=row, column=col,
        extreme=how, n_rows=int(len(vals)), picked=picked, sha256=_sha256(str(paths[0])), commit=_commit(str(paths[0].parent))))


def _last_experiment() -> Resolved:
    """Latest completion time in the experiment index (rows `| E<n> | ... | MM-DD HH:MM (#k) |`,
    or `(staged)` in place of `(#k)` before a Fable agent places the entry in the Completed list)."""
    import re
    from datetime import datetime
    path = H.shared_dir() / 'experiment_priorities.md'
    year = datetime.now().year
    best = None
    for line in path.read_text().split('\n'):
        m = re.match(r'^\| (E\d+) \|.*\| completed \|.*\| (\d{2})-(\d{2}) (\d{2}):(\d{2}) \((?:#\d+|staged)\) \|\s*$', line)
        if m:
            t = datetime(year, int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5)))
            if best is None or t > best[0]:
                best = (t, m.group(1))
    if best is None:
        raise KeyError('meta:last_experiment: no completed experiment in the index')
    return Resolved('meta:last_experiment', best[0].strftime('%Y-%m-%d %H:%M'), dict(kind='meta',
        description=f'completion time of the latest completed experiment ({best[1]}) in experiment_priorities.md\'s index',
        file='experiment_priorities.md', sha256=_sha256(str(path))))


def _json_path(obj: Any, path: str) -> Any:
    """The value at a dotted `path` in `obj`. A dict key may itself contain dots (a variant named 'D1.8_inhadd1000'):
    where a part is not a key, it is joined with the parts after it, shortest first, until one is."""
    parts = path.split('.') if path else []
    i = 0
    while i < len(parts):
        part = parts[i]
        if isinstance(obj, list):
            obj = obj[int(part)]
        elif isinstance(obj, dict) and part in obj:
            obj = obj[part]
        else:
            j = next((j for j in range(i + 2, len(parts) + 1)
                      if isinstance(obj, dict) and '.'.join(parts[i:j]) in obj), None)
            if j is None:
                raise KeyError(part)
            obj = obj['.'.join(parts[i:j])]
            i = j
            continue
        i += 1
    return obj


def _resolve_json(key: str) -> Resolved:
    """`json:<file>:<path>` or `json:<file>:len(<path>)` (module docstring)."""
    _, file, path = key.split(':', 2)
    f = H.shared_dir() / file
    if not f.is_file():
        raise KeyError(f'{key}: no file {file} in the shared directory')
    length = path.startswith('len(') and path.endswith(')')
    inner = path[4:-1] if length else path
    try:
        v = _json_path(json.loads(f.read_text()), inner)
    except (KeyError, IndexError, ValueError) as e:
        raise KeyError(f'{key}: {inner!r} not in {file} ({e})') from None
    if length:
        v = len(v)
    elif isinstance(v, (dict, list)):
        raise KeyError(f'{key}: {inner!r} is a {type(v).__name__}, not a value (len(...) gives its length)')
    return Resolved(key, v, dict(kind='json', file=file, path=path, sha256=_sha256(str(f)), commit=_commit(str(f.parent))))


def _literal_expr(node) -> Any:
    """The value of a constant expression: literals, tuples and lists of them, and + - * / on them."""
    import ast
    import operator
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float, str, bool)):
        return node.value
    if isinstance(node, (ast.Tuple, ast.List)):
        return tuple(_literal_expr(e) for e in node.elts)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        v = _literal_expr(node.operand)
        return -v if isinstance(node.op, ast.USub) else v
    ops = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}
    if isinstance(node, ast.BinOp) and type(node.op) in ops:
        return ops[type(node.op)](_literal_expr(node.left), _literal_expr(node.right))
    raise ValueError(f'not a constant expression: {ast.dump(node)[:80]}')


def _resolve_const(key: str) -> Resolved:
    """`const:<file>:<NAME>[.<index>]` (module docstring)."""
    import ast
    import subprocess
    _, file, name = key.split(':', 2)
    name, _, index = name.partition('.')
    f = H.HERE / file
    if not f.is_file() or f.suffix != '.py':
        raise KeyError(f'{key}: no Python file {file} in the repository')
    nodes = [n for n in ast.parse(f.read_text()).body if isinstance(n, (ast.Assign, ast.AnnAssign))
             and any(isinstance(tg, ast.Name) and tg.id == name for tg in (n.targets if isinstance(n, ast.Assign) else [n.target]))]
    if len(nodes) != 1:
        raise KeyError(f'{key}: {len(nodes)} module-level assignments of {name} in {file}')
    try:
        v = _literal_expr(nodes[0].value)
        if index:
            v = v[int(index)]
    except (ValueError, IndexError, TypeError) as e:
        raise KeyError(f'{key}: {e}') from None
    git = lambda *a: subprocess.run(['git', *a], cwd=H.HERE, capture_output=True, text=True).stdout.strip()
    return Resolved(key, v, dict(kind='const', file=file, name=name + (f'[{index}]' if index else ''),
        line=nodes[0].lineno, sha256=_sha256(str(f)), commit=git('log', '-1', '--format=%H', '--', file) or None,
        dirty=bool(git('status', '--porcelain', '--', file))))


def resolve(key: str) -> Resolved:
    """The value of a placeholder key, with its source."""
    key = key.strip()
    if key == 'meta:generated':
        from datetime import datetime
        return Resolved(key, datetime.now().strftime('%Y-%m-%d %H:%M'), dict(kind='meta', description='time this page was built'))
    if key == 'meta:last_experiment':
        return _last_experiment()
    if key.startswith('diff:'):
        a, sep, b = key[5:].partition(' -- ')
        if not sep or not a.strip() or not b.strip():
            raise KeyError(f'{key}: a difference is written diff:<key> -- <key>')
        ra, rb = resolve(a), resolve(b)
        return Resolved(key, float(ra.value) - float(rb.value), dict(kind='difference',
            a=dict(key=ra.key, source=ra.source, tooltip=tooltip(ra)),
            b=dict(key=rb.key, source=rb.source, tooltip=tooltip(rb))))
    if key.startswith('lit:'):
        text, sep, src = key[4:].partition('#')
        if not sep or not src.strip():
            raise KeyError(f'{key}: a literal needs its source after #')
        return Resolved(key, text.strip(), dict(kind='literal', source=src.strip()))
    if key.startswith(('min:', 'max:')):
        return _resolve_extreme(key)
    if key.startswith('json:'):
        return _resolve_json(key)
    if key.startswith('const:'):
        return _resolve_const(key)
    if key.startswith('data:'):
        name = key[5:]
        part = None
        if name.endswith('.lo') or name.endswith('.hi'):
            name, part = name[:-3], name[-2:]
        reg = _registry()
        if name not in reg:
            raise KeyError(f'{key}: not in the data-target registry')
        e = reg[name]
        v = e['value'] if part is None else e['ci'][0 if part == 'lo' else 1]
        return Resolved(key, v, dict(kind='data', name=name, part=part, definition=e.get('definition'),
            function=e.get('function'), commit=e.get('commit'), file=str(DATA_TARGETS.relative_to(H.shared_dir())),
            sha256=_sha256(str(DATA_TARGETS))))
    return _resolve_table(key)


_WORDS = ('zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve',
          'thirteen', 'fourteen', 'fifteen', 'sixteen', 'seventeen', 'eighteen', 'nineteen', 'twenty')


def format_value(value: Any, spec: Optional[str]) -> str:
    """`value` formatted with `spec` (Python format spec), Unicode minus for negatives; a value
    that rounds to zero carries no minus sign; strings (literals) without a spec are returned as
    written. The spec `w` writes a count in words: the value rounded to an integer, spelled out
    from zero to twenty ("two flies") and in digits outside that range."""
    if isinstance(value, str):
        if not spec:
            return value
        try:
            value = float(value)
        except ValueError:
            return value
    if spec == 'w':
        n = int(round(float(value)))
        return _WORDS[n] if 0 <= n < len(_WORDS) else str(n).replace('-', '−')
    if spec:
        s = format(float(value) if not spec.endswith('d') else int(round(float(value))), spec)
    else:
        s = str(value)
    if s.startswith('-') and not any(c in '123456789' for c in s.split('e')[0].split('E')[0]):
        s = s[1:]
    return s.replace('-', '−')


def tooltip(res: Resolved) -> str:
    """One-line description of a value's source, for the page's hover text."""
    s = res.source
    if s.get('kind') == 'literal':
        return f"hand-entered: {s['source']}"
    if s.get('kind') == 'meta':
        return s['description']
    if s.get('kind') == 'difference':
        return f"difference: ({s['a']['tooltip']}) minus ({s['b']['tooltip']})"
    if s.get('kind') == 'json':
        return f"{s['file']}, {s['path']}" + (f" [commit {s['commit'][:7]}]" if s.get('commit') else '')
    if s.get('kind') == 'const':
        return (f"{s['file']}, {s['name']} (line {s['line']})" + (f" [commit {s['commit'][:7]}]" if s.get('commit') else '')
                + (' (uncommitted changes)' if s.get('dirty') else ''))
    if s.get('kind') == 'data':
        part = f" ({'lower' if s['part'] == 'lo' else 'upper'} end of the 90 % interval)" if s.get('part') else ''
        return f"KC / ORN data: {s['name']}{part}: {s.get('definition') or ''} [kcq/data_targets.py {str(s.get('commit') or '')[:7]}]"
    ext = (f", the {'lowest' if s['extreme'] == 'min' else 'highest'} of {s['n_rows']} rows"
           + (f" ({', '.join(f'{k} {v}' for k, v in s['picked'].items())})" if s.get('picked') else '')) if s.get('extreme') else ''
    return f"{s['file']}, row {s['row']}, column {s['column']}{ext}" + (f" [commit {s['commit'][:7]}]" if s.get('commit') else '')
