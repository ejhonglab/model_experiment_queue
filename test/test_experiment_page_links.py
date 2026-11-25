"""Tests of the links from experiment IDs to the experiments' own pages (E<n>_<short>.html in the shared directory's
reports/): md_to_artifact.py's experiment_pages and experiment-key links (the summary and the digest pages) and
build_report_artifact.py's experiment list (the report page), on scratch folders and on the real templates; nothing
is written outside tmp_path."""
import re
import subprocess
import sys
from pathlib import Path

import pytest

from reports import build_report_artifact as B
from reports import md_to_artifact as M

REPO = Path(__file__).resolve().parents[1]

KEY_MD = """# A page

last model experiment: 2026-11-16 12:00

??? Experiment key: the ID of each experiment and the name used for it in this page
- E1: the APL knockout baseline (completed)
- E10: the joint fit to both matrices (completed)
- E97: the glomeruli's shares of the claws in other connectomes (running; on hold)
- E100: the flies' APL at the weak odors (running)
???

## Results

- E1: a bullet outside the key, not linked
- E10 and E1 in prose.

??? Other fold
- E1: not a key either
???
"""


def pages_folder(tmp_path, names, readme_rows=()):
    d = tmp_path / 'reports'
    d.mkdir()
    for n in names:
        (d / n).write_text('<body>page</body>')
    if readme_rows is not None:
        (d / 'README.md').write_text('# Report copies\n\n| name | title | artifact | copied | host |\n|---|---|---|---|---|\n'
            + ''.join(f'| {r} | {r} | local | 2026-11-16 12:00 | hal |\n' for r in readme_rows))
    return d


def test_experiment_pages_unique_match(tmp_path):
    d = pages_folder(tmp_path, ['E1_apl_knockout.html', 'E10_joint_fit.html', 'E100_apl_weak_odors.html',
        # not experiment pages: other pages, other extensions, uppercase or missing short names, a folder
        'summary_all_experiments.html', 'current-summary.html', 'E2_local_apl.zip', 'E3_Drive.html', 'E4.html',
        'E5_.html.bak', 'note_e6_x.html', 'XE7_x.html'])
    (d / 'E8_folder.html').mkdir()
    assert M.experiment_pages(d) == {1: 'E1_apl_knockout.html', 10: 'E10_joint_fit.html',
        100: 'E100_apl_weak_odors.html'}


def test_experiment_pages_e1_is_not_e10(tmp_path):
    """E1's page is never taken for E10 or E100 (or the other way round), with or without E1 having a page."""
    d = pages_folder(tmp_path, ['E10_joint_fit.html', 'E100_apl_weak_odors.html'])
    assert M.experiment_pages(d) == {10: 'E10_joint_fit.html', 100: 'E100_apl_weak_odors.html'}
    assert M.link_experiment_id('E1: the baseline', M.experiment_pages(d)) == 'E1: the baseline'
    assert M.link_experiment_id('E10: the joint fit', {1: 'E1_a.html'}) == 'E10: the joint fit'


def test_experiment_pages_duplicates(tmp_path):
    # two pages for one experiment: the one in README.md wins
    d = pages_folder(tmp_path, ['E5_kc_heterogeneity.html', 'E5_old_name.html', 'E6_readout_window.html'],
        readme_rows=['E5_kc_heterogeneity', 'E6_readout_window'])
    assert M.experiment_pages(d) == {5: 'E5_kc_heterogeneity.html', 6: 'E6_readout_window.html'}
    # neither listed, or both listed: refused, naming the files
    for rows in ([], ['E5_kc_heterogeneity', 'E5_old_name']):
        (d / 'README.md').write_text(''.join(f'| {r} | t | local | x | hal |\n' for r in rows))
        with pytest.raises(ValueError, match='E5 has several pages.*E5_kc_heterogeneity.html, E5_old_name.html'):
            M.experiment_pages(d)
    (d / 'README.md').unlink()
    with pytest.raises(ValueError, match='none of them is listed'):
        M.experiment_pages(d)


def test_link_experiment_id():
    pages = {1: 'E1_apl_knockout.html', 10: 'E10_joint_fit.html'}
    assert M.link_experiment_id('E1: the baseline (completed)', pages) == \
        '<a href="E1_apl_knockout.html">E1</a>: the baseline (completed)'
    assert M.link_experiment_id('<b>E10, joint fit.</b> More.', pages) == \
        '<b><a href="E10_joint_fit.html">E10</a>, joint fit.</b> More.'
    # no page, no ID at the start, or not followed by ':' / ','
    for item in ('E2: local APL units', 'the E1: baseline', 'E1 and E10', 'E10s: x', '<i>E1:</i> x'):
        assert M.link_experiment_id(item, pages) == item


def test_experiment_key_links(tmp_path):
    d = pages_folder(tmp_path, ['E1_apl_knockout.html', 'E10_joint_fit.html'])
    page = M.convert_markdown(KEY_MD, tmp_path, experiment_pages=M.experiment_pages(d))
    assert M.lint(page) == []
    key = re.search(r'<details class="fold"><summary>Experiment key.*?</details>', page, flags=re.S).group(0)
    assert '<li><a href="E1_apl_knockout.html">E1</a>: the APL knockout baseline (completed)</li>' in key
    assert '<li><a href="E10_joint_fit.html">E10</a>: the joint fit' in key
    # experiments without a page keep plain IDs
    assert '<li>E97: the glomeruli' in key and '<li>E100: the flies' in key
    assert M.experiment_key_links(page) == ([1, 10], [97, 100])
    # only the experiment key's items are linked: not other bullets, prose or other folds
    assert page.count('href="E1_apl_knockout.html"') == 1 and page.count('href="E10_joint_fit.html"') == 1
    assert '<li>E1: a bullet outside the key' in page and '<li>E1: not a key either</li>' in page
    # without pages (the default, every other page): no links, and the same page otherwise
    plain = M.convert_markdown(KEY_MD, tmp_path)
    assert 'href="E' not in plain and M.experiment_key_links(plain) == ([], [1, 10, 97, 100])
    assert re.sub(r'<a href="[^"]*">(E\d+)</a>', r'\1', page) == plain


def test_md_to_artifact_command_line(tmp_path):
    d = pages_folder(tmp_path, ['E10_joint_fit.html'])
    md, out = tmp_path / 'page.md', tmp_path / 'page.html'
    md.write_text(KEY_MD)
    r = subprocess.run([sys.executable, str(REPO / 'reports' / 'md_to_artifact.py'), str(md), str(out),
        '--experiment-pages', str(d)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert 'experiment key: 1 IDs linked to their pages; no page for E1, E97, E100' in r.stdout
    assert '<li><a href="E10_joint_fit.html">E10</a>:' in out.read_text()
    # several pages for one experiment, none in README.md: the build fails
    (d / 'E10_old.html').write_text('x')
    r = subprocess.run([sys.executable, str(REPO / 'reports' / 'md_to_artifact.py'), str(md), str(out),
        '--experiment-pages', str(d)], capture_output=True, text=True)
    assert r.returncode != 0 and 'E10 has several pages' in r.stderr


REPORT_MD = """# The report

## Candidate mechanisms, ranked

1. **E1, not an experiment item here (mechanism 1).** Left alone.

## In silico experiments (olfsysm + al\\_analysis)

1. **E1, APL knockout (mechanism 2 vs 1).** Acute.
2. **E10, joint fit.** Extend.

*Added later: E97 onwards.*

3. **E97, the glomeruli's shares (mechanism 1).** Goal.

## Model changes likely to move distances toward the data

1. **E10, also not linked here.** Text.
"""


def test_report_experiment_list():
    pages = {1: 'E1_apl_knockout.html', 10: 'E10_joint_fit.html'}
    page = B.convert(REPORT_MD, pages)
    assert '<li><b><a href="E1_apl_knockout.html">E1</a>, APL knockout (mechanism 2 vs 1).</b> Acute.</li>' in page
    assert '<li><b><a href="E10_joint_fit.html">E10</a>, joint fit.</b> Extend.</li>' in page
    assert '<li><b>E97, the glomeruli' in page   # no page
    # only in the In silico experiments section
    assert page.count('href="E1_apl_knockout.html"') == 1 and page.count('href="E10_joint_fit.html"') == 1
    assert '<li><b>E1, not an experiment item here' in page and '<li><b>E10, also not linked here.' in page
    # without pages: as before
    assert 'href="E' not in B.convert(REPORT_MD)


@pytest.mark.parametrize('name', ['summary_template.md', 'apl16_summary_template.md', 'since_priority_changes_template.md',
    'since_yang_new_data_template.md'])
def test_templates_keys_all_linkable(tmp_path, name):
    """Every item of the real templates' experiment keys has the `E<n>: ` form the links need, so each gets its
    link once its page exists."""
    md = (REPO / 'reports' / name).read_text()
    key = re.search(r'^(\?\?\? Experiment key.*?\n)(.*?)^\?\?\?$', md, flags=re.S | re.M)
    items = [ln for ln in key.group(2).split('\n') if ln.strip()]
    ids = [int(m.group(1)) for ln in items for m in [re.match(r'- E(\d+): ', ln)] if m]
    assert len(ids) == len(items) > 0
    page = M.convert_markdown('# T\n\n' + key.group(0) + '\n', tmp_path,
        experiment_pages={n: f'E{n}_x.html' for n in ids})
    assert M.experiment_key_links(page) == (ids, [])


def test_report_list_all_linkable():
    """Every item of the real report's experiment list starts with `**E<n>, `, so each gets its link once its page
    exists."""
    md = (REPO / 'reports' / 'apl_orn_to_kc_remapping_report.md').read_text()
    section = md.split('\n## In silico experiments', 1)[1].split('\n## ', 1)[0]
    ids = [int(n) for n in re.findall(r'^\d+\. \*\*E(\d+), ', section, flags=re.M)]
    assert len(ids) == len(re.findall(r'^\d+\. ', section, flags=re.M)) > 100
    page = B.convert(md, {n: f'E{n}_x.html' for n in ids})
    assert [int(n) for n in re.findall(r'<li><b><a href="E\d+_x\.html">E(\d+)</a>, ', page)] == ids
