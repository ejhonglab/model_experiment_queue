#!/usr/bin/env python3
"""Build the "since the priority changes" page from its template (reports/since_priority_changes_template.md).

What was learned since the last two changes of the main goal: a snapshot written at Tom's request on 2026-11-13,
rebuilt only on request (hosts do not keep it current at each completion, unlike the summary). It reads only what
build_summary.py already produced (the data-target registry, the experiments' results tables, figures in the shared
directory): it does not run kcq/apl_effect.py or the progress plot, and does not touch the summary's outputs.

1. fills the template from the data-target registry and the experiments' results tables (fill_report.py; fails on
   any unresolved placeholder before writing anything, prints lint for bare numbers in KC-data rows);
2. writes the plain filled markdown (no tooltip markup) to the shared directory's
   add_to_reports/since_priority_changes.md (generated: edit the template, never that file), and the sidecar of
   every value's source next to it (since_priority_changes.values.json, the same fields as the summary's);
3. builds the page with source tooltips (md_to_artifact.py) at outputs/_artifacts/since_priority_changes/index.html.
The template's syntax is the summary's (`{{...}}` placeholders, `{{meta:generated}}` and `{{meta:last_experiment}}`
included; `%%% hover` blocks, `[TOC]`, `{#id}` anchors, `??? ... ???` folds, `{span=N}` table headers). Then save the
copy the queue monitor shows (Report tab, "Since the priority changes"; not published as an artifact):
    ./venv/bin/python reports/save_report_copy.py outputs/_artifacts/since_priority_changes/index.html \\
        since_priority_changes --current current-since-priority-changes.html
Usage: ./venv/bin/python reports/build_since_priority_changes.py [--template T] [--out-md MD] [--out-html HTML]
(the command line is build_summary.page_main's, which build_since_yang_new_data.py shares)
"""
from __future__ import annotations

from reports.build_summary import page_main

NAME = 'since_priority_changes'

if __name__ == '__main__':
    page_main(NAME, __doc__)
