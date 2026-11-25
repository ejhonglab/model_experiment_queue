#!/usr/bin/env python3
"""Build the all-experiments summary report, "Summary", from its template (reports/summary_template.md).

Since 2026-11-25 the Summary is organized along the stage target (CLAUDE.md, "The main goal for new experiments");
the summary organized around the APL effect until then is the APL16 Summary (reports/apl16_summary_template.md,
built by build_apl16_summary.py on request).

1. scores every saved model variant's per-pair APL effect and draws the matrix figures of both summaries'
   buttons (kcq/apl_effect.py), regenerates the APL16 Summary's progress plot (summary_progress_plot.py),
   draws the Summary's stage figures (summary_stage_figures.py) and regenerates its stage-target progress plot
   (summary_stage_progress_plot.py; the template inserts it once the fit of all megamat stages is scored);
2. fills the template from the data-target registry and the experiments' results tables
   (fill_report.py; fails on any unresolved placeholder, prints lint for bare numbers in KC-data
   rows and the best-models table);
3. writes the plain filled markdown (no tooltip markup) to the shared directory's
   add_to_reports/summary_all_experiments.md (generated: edit the template, never
   that file), and the sidecar of every value's source next to it
   (summary_all_experiments.values.json);
4. builds the page with source tooltips (md_to_artifact.py) at
   outputs/_artifacts/summary/index.html, the IDs of its experiment key linked to the experiments' own
   pages in the shared directory's reports/ (relative links: they work in the copy saved there).
save_report_copy.py stays manual; the page is not published as an artifact (CLAUDE.md, "Report copies").
Steps 2-4 are fill_and_build(), which the digest pages (build_since_priority_changes.py,
build_since_yang_new_data.py) use for their own templates, through page_main().
Usage: ./venv/bin/python reports/build_summary.py [--skip-progress-plot] [--skip-stage-figures]
       [--skip-stage-progress-plot] [--data-targets]
       (--data-targets first rebuilds the registry, kcq/data_targets.py)
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from kcq import harness as H
from reports import fill_report as F

TEMPLATE = H.HERE / 'reports' / 'summary_template.md'
OUT_MD = H.shared_dir() / 'add_to_reports' / 'summary_all_experiments.md'
OUT_HTML = H.HERE / 'outputs' / '_artifacts' / 'summary' / 'index.html'


def fill_and_build(template_path: Path, out_md: Path, out_html: Path, filled_name: str) -> bool:
    """Fills a report template and builds its page (steps 2-4 of this module's docstring), for the
    summary's template or another page's with the same syntax.

    Args:
        template_path: the template.
        out_md: where the plain filled markdown goes; its sidecar of every value's source is written
            next to it (`.values.json`).
        out_html: where the page with source tooltips goes (md_to_artifact.py, figures relative to
            the shared directory).
        filled_name: file name of the filled markdown with tooltip markup, written next to
            `out_html` as md_to_artifact.py's input.

    Returns:
        False if md_to_artifact.py failed or found markdown left in the page.

    Raises:
        KeyError: a placeholder did not resolve (before anything is written).
    """
    template = template_path.read_text()
    probs = F.lint(template)
    for p in probs:
        print('LINT', p)
    plain, side = F.fill(template, tooltips=False)
    tipped, _ = F.fill(template, tooltips=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(plain)
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=H.HERE, capture_output=True, text=True).stdout.strip()
    try:
        template_name = str(template_path.resolve().relative_to(H.HERE))
    except ValueError:
        template_name = str(template_path.resolve())
    out_md.with_suffix('.values.json').write_text(json.dumps(dict(template=template_name, template_commit=commit,
        generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), n_values=len(side), lint=probs, values=side),
        indent=1, ensure_ascii=False, default=str) + '\n')
    out_html.parent.mkdir(parents=True, exist_ok=True)
    tipped_md = out_html.parent / filled_name
    tipped_md.write_text(tipped)
    # the experiment key's IDs link to the experiments' own pages, saved next to this page's copy in reports/
    pages_dir = H.shared_dir() / 'reports'
    r = subprocess.run([sys.executable, str(H.HERE / 'reports' / 'md_to_artifact.py'), str(tipped_md), str(out_html), '--figdir', str(H.shared_dir())]
        + (['--experiment-pages', str(pages_dir)] if pages_dir.is_dir() else []),
        capture_output=True, text=True)
    print('\n'.join(l for l in (r.stdout + r.stderr).split('\n') if l and 'Warning' not in l and 'warnings.warn' not in l)[-600:])
    print(f'{len(side)} values filled; {len(probs)} lint problems; markdown -> {out_md}; page -> {out_html}')
    return r.returncode == 0 and 'LINT PROBLEMS' not in r.stdout


def page_main(name: str, doc: str) -> None:
    """The command line of a page built from its own template with fill_and_build() alone (the digest pages:
    build_since_priority_changes.py, build_since_yang_new_data.py), reading only what main() already produced.

    Options `--template`, `--out-md` and `--out-html` override the defaults below; exits non-zero on a missing
    template, an unresolved placeholder (nothing written) or a failed page build.

    Args:
        name: the page's name: its template is reports/<name>_template.md, its plain filled markdown the shared
            directory's add_to_reports/<name>.md (sidecar <name>.values.json next to it), its page
            outputs/_artifacts/<name>/index.html (beside <name>_filled.md, the markdown with tooltip markup).
        doc: the calling script's docstring; its first paragraph is the command's description.
    """
    ap = argparse.ArgumentParser(description=doc.split('\n\n')[0])
    ap.add_argument('--template', type=Path, default=H.HERE / 'reports' / f'{name}_template.md',
        help='default: %(default)s')
    ap.add_argument('--out-md', type=Path, default=None,
        help=f"the plain filled markdown, its .values.json sidecar next to it (default: the shared directory's "
        f'add_to_reports/{name}.md)')
    ap.add_argument('--out-html', type=Path, default=H.HERE / 'outputs' / '_artifacts' / name / 'index.html',
        help='the page (default: %(default)s)')
    args = ap.parse_args()
    if not args.template.is_file():
        sys.exit(f'no template at {args.template}')
    out_md = args.out_md or H.shared_dir() / 'add_to_reports' / f'{name}.md'
    try:
        ok = fill_and_build(args.template, out_md, args.out_html, f'{name}_filled.md')
    except KeyError as e:
        sys.exit(e.args[0])
    if not ok:
        sys.exit(1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--skip-progress-plot', action='store_true')
    ap.add_argument('--data-targets', action='store_true')
    ap.add_argument('--skip-stage-figures', action='store_true',
        help="don't redraw the Summary's stage figures (summary_stage_figures.py)")
    ap.add_argument('--skip-stage-progress-plot', action='store_true',
        help="don't regenerate the Summary's stage-target progress plot (summary_stage_progress_plot.py)")
    args = ap.parse_args()
    py = sys.executable
    if args.data_targets:
        subprocess.run([py, str(H.HERE / 'kcq' / 'data_targets.py')], check=True)
    if not args.skip_progress_plot:
        subprocess.run([py, str(H.HERE / 'kcq' / 'apl_effect.py')], check=True)
        subprocess.run([py, str(H.HERE / 'reports' / 'summary_progress_plot.py')], check=True)
    if not args.skip_stage_figures:
        subprocess.run([py, str(H.HERE / 'reports' / 'summary_stage_figures.py')], check=True)
    if not args.skip_stage_progress_plot:
        subprocess.run([py, str(H.HERE / 'reports' / 'summary_stage_progress_plot.py')], check=True)
    if not fill_and_build(TEMPLATE, OUT_MD, OUT_HTML, 'summary_filled.md'):
        sys.exit(1)


if __name__ == '__main__':
    main()
