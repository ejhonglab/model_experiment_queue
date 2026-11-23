#!/usr/bin/env python3
"""Build "The last 1.5 weeks" page from its template (reports/recent_weeks_template.md).

A concise summary of the experiments completed from 2026-11-12 to 2026-11-23 and of the two then running, written at
Tom's request on 2026-11-23 and rebuilt only on request (hosts do not keep it current at each completion, unlike the
summary). It is built as the other digest pages are (build_summary.page_main: the same template syntax and steps,
reading only what build_summary.py already produced), to the shared directory's add_to_reports/recent_weeks.md
(generated: edit the template, never that file) with its sidecar recent_weeks.values.json, and the page
outputs/_artifacts/recent_weeks/index.html. Then save the copy the queue monitor shows (Report tab, "The last 1.5
weeks"; not published as an artifact):
    ./venv/bin/python reports/save_report_copy.py outputs/_artifacts/recent_weeks/index.html \\
        recent_weeks --current current-recent-weeks.html
Usage: ./venv/bin/python reports/build_recent_weeks.py [--template T] [--out-md MD] [--out-html HTML]
"""
from __future__ import annotations

from reports.build_summary import page_main

NAME = 'recent_weeks'

if __name__ == '__main__':
    page_main(NAME, __doc__)
