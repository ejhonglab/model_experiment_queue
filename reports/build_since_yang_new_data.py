#!/usr/bin/env python3
"""Build the "since Yang's new data" page from its template (reports/since_yang_new_data_template.md).

What was learned from Yang's new data, by type of data, with what remains to be checked: a snapshot written at Tom's
request on 2026-11-16, rebuilt only on request (hosts do not keep it current at each completion, unlike the summary).
It is built as build_since_priority_changes.py builds its page (build_summary.page_main: the same template syntax and
steps, reading only what build_summary.py already produced), to the shared directory's
add_to_reports/since_yang_new_data.md (generated: edit the template, never that file) with its sidecar
since_yang_new_data.values.json, and the page outputs/_artifacts/since_yang_new_data/index.html. Then save the copy
the queue monitor shows (Report tab, "Since Yang's new data"; not published as an artifact):
    ./venv/bin/python reports/save_report_copy.py outputs/_artifacts/since_yang_new_data/index.html \\
        since_yang_new_data --current current-since-yang-new-data.html
Usage: ./venv/bin/python reports/build_since_yang_new_data.py [--template T] [--out-md MD] [--out-html HTML]
"""
from __future__ import annotations

from reports.build_summary import page_main

NAME = 'since_yang_new_data'

if __name__ == '__main__':
    page_main(NAME, __doc__)
