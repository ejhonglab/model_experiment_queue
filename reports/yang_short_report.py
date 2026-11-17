#!/usr/bin/env python3
"""A short report for Yang on one topic, from a hand-written markdown file, filled and converted as her version of the
summary page is (reports/yang_summary.py).

The markdown (default reports/yang_summary/inhibition_short.md) is a whole page: a `# ` title line, then anything
reports/md_to_artifact.py converts ([TOC], tables, `??? title` ... `???` collapsed blocks, `%%% hover` blocks, figures
`![caption](results/...png)` relative to the shared directory). Its placeholders are those of the answers in
reports/yang_summary/questions.md: `{{key ~ spec}}` with the keys of reports/report_values.py (`data:<name>`,
`<EXP>[/<table.csv>]:<row>:<column>`, `lit:<text> # <source>`, `diff:<key> -- <key>`, `meta:generated`,
`meta:last_experiment`), and `{{inv:KEY ~ spec}}` (1 minus the value: a distance shown as a correlation) and
`{{neg:KEY ~ spec}}` (minus the value). They are filled by yang_summary.Filler, as her answers are: every value's
hover text names its source, and an `inv:` or `neg:` value's also gives the main page's value. `meta:` keys resolve as
on her page (report_values.resolve, once per build).

- build [MARKDOWN] [--name NAME]: writes outputs/_artifacts/<NAME>/ (NAME default: the markdown's file name without
  .md, plus `_yang`): index.html (the page, titled by the markdown's `# ` line), <NAME>_filled.md (the markdown that
  was converted, each value with its source), <NAME>.md (the filled markdown as plain text) with its sidecar
  <NAME>.values.json (every value: line, key, spec, text, value, source) and build_report.json (inputs, values filled
  by kind, lint). An unresolved placeholder, a missing figure or a problem found by md_to_artifact.lint fails the
  build (exit status 1) with the line named: the outputs of an earlier build are removed, index.html included, and
  build_report.json records the failure. Bare numbers in "KC data" rows (fill_report.lint) are warnings, as on her
  page. NAME `summary` or `summary_...` is refused (in any case, and `summary-...` and `summary....` too): those names
  belong to the summary pages (summary_yang is hers). So is a folder of outputs/_artifacts/ that another tool wrote.
  Nothing is written to the shared directory.

Usage: ./venv/bin/python reports/yang_short_report.py build [MARKDOWN] [--name NAME]
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import html as html_mod
import json
import os
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from kcq import harness as H
from reports import fill_report as F
from reports import md_to_artifact as M
from reports import yang_summary as Y

DEFAULT_MARKDOWN = H.HERE / 'reports' / 'yang_summary' / 'inhibition_short.md'
ARTIFACTS = H.HERE / 'outputs' / '_artifacts'
TOOL = 'reports/yang_short_report.py'   # recorded in build_report.json: the folders this tool may rewrite
NAME_OK = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*')   # ASCII only, one path component (CLAUDE.md, "File and folder names")
NAME_TAKEN = re.compile(r'summary($|[_.-])', re.I)   # the summary pages' names (summary, summary_yang, summary_...)
MARKERS = ('⟦', '⟧')   # fill_report's value markers, `⟦pv:value§source⟧`
H123 = re.compile(r'^#{1,3} ')   # the headings md_to_artifact.py converts (deeper ones stay paragraphs)
TAGS = re.compile(r'<[^<>]*>')
QUOTED = re.compile(r": ('(?:[^'\\]|\\.)*'|\"(?:[^\"\\]|\\.)*\")$")   # the text md_to_artifact.lint quotes, as repr


class Refused(ValueError):
    """A page name or output folder this tool may not use (nothing has been written or removed)."""


class BuildError(RuntimeError):
    """The page cannot be built from the markdown as it is.

    Attributes:
        what: what failed (e.g. 'unresolved placeholders').
        problems: one message per problem, naming its line of the markdown where it can be found.
        report: the build report written for the failure (set by `build`).
    """

    def __init__(self, what: str, problems: Sequence[str]):
        self.what, self.problems, self.report = what, list(problems), {}
        super().__init__(f'{what}:\n' + '\n'.join('  ' + p for p in self.problems))


def _shared() -> Path:
    """The shared directory (a function, so that tests can point it elsewhere)."""
    return H.shared_dir()


# ---------------------------------------------------------------- names and folders

def default_name(markdown: Path) -> str:
    """The page name of a markdown file: its file name without .md, plus `_yang` unless it ends so already."""
    stem = Path(markdown).stem
    return stem if stem.endswith('_yang') else stem + '_yang'


def check_name(name: str) -> str:
    """`name` if this tool may write outputs/_artifacts/<name>/; raises Refused otherwise."""
    if not NAME_OK.fullmatch(name or ''):
        raise Refused(f'page name {name!r}: use ASCII letters, digits, "_", "-" and "." only, starting with a letter or '
            'a digit (pass --name)')
    if NAME_TAKEN.match(name):
        raise Refused(f'page name {name!r}: the names summary and summary_... belong to the all-experiments summary pages '
            "(summary_yang is Yang's version of it); pass another --name")
    return name


def check_folder(out: Path) -> None:
    """Raises Refused if `out` exists and holds anything but this tool's outputs (by its build_report.json)."""
    if out.exists() and not out.is_dir():
        raise Refused(f'{Y._rel(out)} exists and is not a folder')
    if not out.is_dir() or not any(out.iterdir()):
        return
    try:
        tool = json.loads((out / 'build_report.json').read_text()).get('tool')
    except (OSError, ValueError, AttributeError):
        tool = None
    if tool != TOOL:
        raise Refused(f'{Y._rel(out)} holds another page\'s outputs (no build_report.json written by {TOOL}); '
            'pass another --name')


def outputs(out: Path, name: str) -> Dict[str, Path]:
    """The files a build writes in `out`."""
    return dict(page=out / 'index.html', filled_markdown=out / f'{name}_filled.md', markdown=out / f'{name}.md',
        values=out / f'{name}.values.json', report=out / 'build_report.json')


# ---------------------------------------------------------------- the page

def page_title(line: str) -> str:
    """The page title from a filled `# ` line (no tooltips): its text as the converter renders it, as plain text."""
    text = Y.HEAD_ID.sub('', line.strip()[2:])
    return html_mod.unescape(re.sub(r'<[^>]+>', '', M.inline(text))).strip()


def figure_refs(lines: Sequence[Y.Line], plain: Sequence[str]) -> List[Tuple[int, str]]:
    """(line, path) of every file the converter embeds: figure lines, buttons in table rows, `%%% html` blocks."""
    refs = []
    for ln in lines:
        t = plain[ln.n - 1].strip()
        if ln.kind == 'figure':
            refs += [(ln.n, m.group(2)) for m in [M._FIG_RE.match(t)] if m]
        elif ln.kind == 'table':
            refs += [(ln.n, m.group(2)) for m in M._BUTTON_RE.finditer(t)]
        elif ln.kind == 'structural':
            refs += [(ln.n, m.group(1)) for m in [re.match(r'^%%% html (\S+)$', t)] if m]
    return refs


def _head_id(line: str) -> Optional[str]:
    """The anchor md_to_artifact.py gives a `#`-`###` heading line (its `{#id}`, else the slug of its text)."""
    s = line.strip()
    if not H123.match(s):
        return None
    m = Y.HEAD_ID.search(s)
    return m.group(1) if m else M._slug(s.split(' ', 1)[1])


def lint_lines(msg: str, plain: Sequence[str]) -> List[int]:
    """The markdown lines (1-based) an md_to_artifact.lint message is about, as far as they can be found: the line
    holding the start of the text it quotes (unique, or one of at most three), or a footnote's heading."""
    m = re.match(r'footnote #(\S+) is not linked', msg)
    if m:
        return [n for n, t in enumerate(plain, 1) if _head_id(t) == m.group(1)]
    if msg.startswith('table-of-contents marker'):
        snippet = '[TOC]'
    else:
        q = QUOTED.search(msg)
        if not q:
            return []
        try:
            snippet = html_mod.unescape(ast.literal_eval(q.group(1))).strip()
        except (ValueError, SyntaxError):
            return []
    # the quoted text starts a text node of the page: a piece of one line (a paragraph's lines are joined)
    for k in range(len(snippet), min(len(snippet), 8) - 1, -1):
        hits = [n for n, t in enumerate(plain, 1) if snippet[:k] in t]
        if hits:
            return hits if len(hits) <= 3 else []
    return []


def _at(lines: Sequence[int]) -> str:
    return f'line {lines[0]}: ' if len(lines) == 1 else (f'lines {", ".join(map(str, lines))}: ' if lines else '')


def make(text: str, shared: Path) -> dict:
    """Fills and converts a short report (the build without its files).

    Args:
        text: the markdown.
        shared: the shared directory (figure paths are relative to it).

    Returns:
        dict(page, filled (the markdown converted), plain (filled, no tooltips), side (sidecar entries), title,
        figures, fill_lint).

    Raises:
        BuildError: value markers in the text, no title line, an unresolved placeholder, a missing figure or a lint
            problem.
    """
    src = text.split('\n')
    bad = [f'line {n}: contains {" or ".join(c for c in MARKERS if c in t)}, which marks filled values'
        for n, t in enumerate(src, 1) if any(c in t for c in MARKERS)]
    if bad:
        raise BuildError('characters the filled page uses as markers', bad)
    lines = Y.analyze(text)
    title_at = next((ln.n for ln in lines if ln.kind == 'heading' and ln.text.strip().startswith('# ')), None)
    if title_at is None:
        raise BuildError('no title', ['the markdown has no "# " title line'])
    filler = Y.Filler()
    filled = [filler.line(t, f'line {n}', tooltips=True) for n, t in enumerate(src, 1)]
    if filler.errors:
        raise BuildError('unresolved placeholders', filler.errors)
    plain = [filler.line(t, f'line {n}', tooltips=False, record=dict(line=n)) for n, t in enumerate(src, 1)]
    title = page_title(plain[title_at - 1])
    if not title:
        raise BuildError('no title', [f'line {title_at}: the title line is empty'])
    # a heading with a value keeps the anchor of its plain text (the value's source would be in its slug)
    for ln in lines:
        k = ln.n - 1
        if ln.kind == 'heading' and H123.match(src[k].strip()) and filled[k] != plain[k] and not Y.HEAD_ID.search(src[k]):
            filled[k] = filled[k].rstrip() + ' {#' + _head_id(plain[k]) + '}'
    figures, missing = [], []
    for n, ref in figure_refs(lines, plain):
        p = Path(ref) if Path(ref).is_absolute() else shared / ref
        if p.is_file():
            figures.append(dict(line=n, path=ref, sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        else:
            missing.append(f'line {n}: not found: {ref} (looked for {p})')
    if missing:
        raise BuildError('missing figures', missing)
    try:
        page = M.convert_markdown('\n'.join(filled), shared, html_mod.escape(title, quote=False))
    except FileNotFoundError as e:   # a file the checks above did not see
        raise BuildError('missing figures', [str(e)]) from None
    # values inside attributes (a figure's alt text, a hover definition's title) keep their text only
    page = TAGS.sub(lambda m: M.PV.sub(r'\1', m.group(0)), page)
    probs = M.lint(page)
    probs += [f'a filled-value marker was left in the page: {t[:70]!r}' for t in re.findall(r'⟦pv:[^⟧]*⟧?', page)]
    if probs:
        raise BuildError('the page does not pass md_to_artifact.lint', [_at(lint_lines(p, plain)) + p for p in probs])
    return dict(page=page, filled=filled, plain=plain, side=filler.side, title=title, figures=figures,
        fill_lint=F.lint(text))


# ---------------------------------------------------------------- build

def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(['git', *args], cwd=H.HERE, capture_output=True, text=True)


def git_state(path: Path) -> str:
    """'committed', 'modified', 'untracked', 'ignored' or 'outside the repository'."""
    r = _git('status', '--porcelain', '--ignored', '--', str(Path(path).resolve()))
    if r.returncode:
        return 'outside the repository'
    return {'': 'committed', '??': 'untracked', '!!': 'ignored'}.get(r.stdout[:2], 'modified')


def value_counts(side: Sequence[dict]) -> dict:
    """Values filled: in all, plain / inv / neg, and by kind of source (data, table, literal, difference, meta)."""
    t = Counter(e.get('transform') for e in side)
    return dict(filled=len(side), plain=t[None], inv=t['inv'], neg=t['neg'],
        by_source=dict(sorted(Counter(str((e.get('source') or {}).get('kind')) for e in side).items())))


def _dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, ensure_ascii=False, default=str) + '\n')


def build(markdown: Path = DEFAULT_MARKDOWN, name: Optional[str] = None, *, shared: Optional[Path] = None,
        out_root: Optional[Path] = None) -> dict:
    """Builds a short report (steps in the module docstring).

    Args:
        markdown: the report's markdown.
        name: the page's folder under `out_root` (default: `default_name(markdown)`).
        shared: the shared directory (default H.shared_dir()); only read.
        out_root: where the page's folder goes (default outputs/_artifacts/).

    Returns:
        The build report (also written to build_report.json).

    Raises:
        Refused: the name or its folder may not be used (nothing is written or removed).
        BuildError: the page cannot be built (the outputs of an earlier build are removed; build_report.json records
            the failure, and the exception's `report` holds it).
    """
    markdown = Path(markdown)
    name = check_name(name if name is not None else default_name(markdown))
    out = Path(out_root or ARTIFACTS) / name
    check_folder(out)
    text = markdown.read_text()
    shared = Path(shared or _shared())
    files = outputs(out, name)
    report = dict(tool=TOOL, ok=False, name=name, title=None, generated=datetime.now().strftime('%Y-%m-%d %H:%M'),
        inputs=dict(markdown=Y._rel(markdown), markdown_sha1=hashlib.sha1(text.encode()).hexdigest(),
            markdown_git=git_state(markdown), commit=_git('rev-parse', 'HEAD').stdout.strip(), shared=str(shared)))
    # nothing of an earlier build survives a failed one (index.html goes first), so that a broken page never looks
    # current; the report, rewritten at once as "building", keeps the folder this tool's if the build is interrupted
    removed = [Y._rel(p) for k, p in files.items() if k != 'report' and p.exists()]
    for k, p in files.items():
        if k != 'report':
            p.unlink(missing_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    _dump(files['report'], dict(report, state='building'))
    try:
        r = make(text, shared)
    except BuildError as e:
        report.update(failed=e.what, errors=e.problems, removed=removed, outputs=dict(report=Y._rel(files['report'])))
        _dump(files['report'], report)
        e.report = report
        raise
    side, page = r['side'], r['page']
    counts = value_counts(side)
    files['filled_markdown'].write_text('\n'.join(r['filled']))
    files['markdown'].write_text('\n'.join(r['plain']))
    _dump(files['values'], dict(markdown=report['inputs']['markdown'], markdown_sha1=report['inputs']['markdown_sha1'],
        commit=report['inputs']['commit'], generated=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        n_values=len(side), n_inv=counts['inv'], n_neg=counts['neg'], lint=r['fill_lint'], values=side))
    tmp = out / '.index.html.tmp'
    tmp.write_text(page)
    os.replace(tmp, files['page'])
    report.update(ok=True, title=r['title'], inputs=dict(report['inputs'], figures=r['figures']), values=counts,
        fill_lint=r['fill_lint'], lint=[], outputs={k: Y._rel(p) for k, p in files.items()}, bytes=len(page.encode()))
    _dump(files['report'], report)
    return report


# ---------------------------------------------------------------- command line

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build', help='fill and convert a short report into outputs/_artifacts/<NAME>/')
    b.add_argument('markdown', nargs='?', type=Path, default=DEFAULT_MARKDOWN,
        help=f'the report (default: {Y._rel(DEFAULT_MARKDOWN)})')
    b.add_argument('--name', default=None, help="the page's folder under outputs/_artifacts/ (default: the markdown's "
        'file name without .md, plus _yang); never summary or summary_...')
    args = ap.parse_args()
    if not args.markdown.is_file():
        ap.error(f'{args.markdown}: no such file')
    try:
        rep = build(args.markdown, args.name)
    except Refused as e:
        ap.error(str(e))
    except BuildError as e:
        for p in e.report.get('removed', []):
            print(f'removed {p} (from an earlier build)')
        print(f'wrote {e.report["outputs"]["report"]} (the failure)')
        sys.exit(f'FAILED: {Y._rel(args.markdown)}: {e}\nno index.html written')
    v = rep['values']
    print(f'values: {v["filled"]} filled ({v["inv"]} inv, {v["neg"]} neg; by source: '
        + (', '.join(f'{n} {k}' for k, n in v['by_source'].items()) or 'none') + ')')
    for w in rep['fill_lint']:
        print('WARNING (bare number, fill_report.lint)', w)
    o = rep['outputs']
    for k in ('filled_markdown', 'markdown', 'values'):
        print(f'wrote {o[k]}')
    print(f'wrote {o["page"]} ({rep["bytes"] // 1024} KB): "{rep["title"]}"')
    print(f'wrote {o["report"]}')
    print('lint clean')


if __name__ == '__main__':
    main()
