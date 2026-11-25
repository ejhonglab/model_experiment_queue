#!/usr/bin/env python3
"""Build the APL16 Summary from its template (reports/apl16_summary_template.md).

The APL16 Summary is the all-experiments summary as it was organized around the APL effect on Yang's APL16 panel, the
main target until 2026-11-24. On 2026-11-25 Tom asked for a new summary along the stage target, titled "Summary"
(reports/summary_template.md, built by build_summary.py), and for this one to be renamed "APL16 Summary". It is kept as
it stood then and rebuilt only on request: new experiments are folded into the Summary. It is built as the digest pages
are (build_summary.page_main: the same template syntax and steps, reading only what build_summary.py already produced:
the APL-effect scores and matrices, the progress plot), to the shared directory's add_to_reports/apl16_summary.md
(generated: edit the template, never that file) with its sidecar apl16_summary.values.json, and the page
outputs/_artifacts/apl16_summary/index.html. Then save the copy the queue monitor shows (Report tab, "APL16 Summary";
not published as an artifact):
    ./venv/bin/python reports/save_report_copy.py outputs/_artifacts/apl16_summary/index.html apl16_summary
Usage: ./venv/bin/python reports/build_apl16_summary.py [--template T] [--out-md MD] [--out-html HTML]
"""
from __future__ import annotations

from reports.build_summary import page_main

NAME = 'apl16_summary'

if __name__ == '__main__':
    page_main(NAME, __doc__)
