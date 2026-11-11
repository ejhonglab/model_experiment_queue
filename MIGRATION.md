# Moving the queue out of al_analysis (cutover checklist)

Until 2026-11-11 the queue lived in al_analysis, as `scripts/claude_pair_remapping_experiments`
(below, "the old folder"). This repository holds its history (a `git subtree split` of al_analysis
16f552e, 87c871f here) and the changes for living on its own. Al_analysis keeps its history: commits
cited before the move are al_analysis commits, and the folder is deleted there in a new commit.
This checklist is for the pause after the round that is running from the old folder. Do the steps
in the order a, c, b, d, e, f. The host step (c) comes before al_analysis's removal commit (b), so
the queue monitor and the hooks never point at deleted files.
`migration/cutover_hal.sh` does what can be scripted; run its phases from this repository on hal.

## Before the pause

- Every host commits and pushes its changes to the old folder in al_analysis, including hal's two
  commits that were not on origin at 12:45 (14fa8e6, 25f58af) and its uncommitted
  E83_cross_wiring_fit.py. `git -C ~/src/al_analysis status --short scripts/claude_pair_remapping_experiments`
  is empty on hal, atlas and viki.
- No model run, chat watcher or session runs in or from the old folder on any host. The
  `host` phase checks this on hal and join_queue.sh checks it on viki.

## (a) Merge the folder's last commits

On hal, from this repository: `migration/cutover_hal.sh merge`. It fetches al_analysis's main,
splits the old folder's history again in a scratch clone (the split is deterministic: done again at
16f552e it gives 87c871f exactly), and merges the commits made since. If main was rewritten in the
meantime (from_remy/ and apl_pair_analysis/ removed, then force pushed by Tom), the raw split must not be
merged: it would bring the old history and both folders back. The script then finds 87c871f's
rewritten counterpart by author, date and subject. It replays the new commits onto it without those
two folders (a commit that only touched them is dropped) and merges the result. Both paths were
tested on a scratch clone, with filter-branch standing in for the rewrite. The phase then runs
`migration/copy_al_analysis_tests.sh` again, for tests added to or changed in al_analysis's test/
since 699e470, and commits them. Conflicts stop it; resolve and commit by hand. Then run the tests
(below) and `git push origin HEAD:main`.

## (c) hal itself

`migration/cutover_hal.sh host`, once nothing runs in the old folder (stop the chat watcher with
TaskStop; the session that runs the script may stay):
1. moves `outputs/`, `from_yang/` and `.fable_consent.local` from the old folder here (`mv`: both
   are on /home; entry by entry where both exist, keeping duplicates for a person to merge);
2. makes `./venv` (python3.14) with `pip install -e ../al_analysis[dev]`, then `../hong2p[all]`, then
   `../olfsysm`. al_analysis comes first because it pins the other two to their GitHub versions.
   olfsysm is installed only when no process has the extension in ~/src/olfsysm loaded, since `pip
   install -e` rebuilds it in place. Otherwise the phase lists those processes; run
   `migration/cutover_hal.sh olfsysm` once they end. Until then fit_mb_model refuses to run
   (olfsysm not clean);
3. relinks the queue monitor:
   `systemctl --user stop/disable claude-queue-monitor` (disable also removes the link to the old
   folder's unit file), `systemctl --user link ~/src/model_experiment_queue/monitor-server/claude-queue-monitor.service`,
   `daemon-reload`, `enable --now`;
4. points the 12 hook commands in `~/src/dotfiles/claude_settings.json` (= `~/.claude/settings.json`,
   lines 30, 39, 51, 63, 75, 87, 99, 111, 123, 135, 147, 159) here: in each,
   `$HOME/src/al_analysis/scripts/claude_pair_remapping_experiments/` becomes
   `$HOME/src/model_experiment_queue/`, so line 30 reads
   `"command": "f=\"$HOME/src/model_experiment_queue/session_title_hook.sh\"; [ -x \"$f\" ] && exec \"$f\"; cat >/dev/null; exit 0",`
   and the other 11 read
   `"command": "f=\"$HOME/src/model_experiment_queue/monitor-server/presence_hook.py\"; [ -f \"$f\" ] && exec timeout -k 2 15 python3 -I \"$f\"; cat >/dev/null; exit 0",`
   Running sessions reload the settings. The change is left uncommitted: dotfiles is Tom's and not
   one of the queue's four repositories, so Tom commits and pushes it (and atlas pulls it).

Then start hal's queue session(s) in ~/src/model_experiment_queue (accept the workspace-trust
prompt; the session title is unchanged, "al_analysis modeling - hal") and restart the chat watcher
from there (`./chat_watch.sh`). presence_hook.py counts sessions in both folders, so the monitor
keeps showing sessions still in the old one.

Not scripted, for Tom to decide: Claude Code keeps auto-memory per project directory. Sessions in
the old folder used `~/.claude/projects/-home-tom-src-al-analysis/memory/` (the repository root).
Sessions here get `~/.claude/projects/-home-tom-src-model-experiment-queue/memory/`, which starts
empty. Copy the KC-remapping entries over if they should carry on (most are already in CLAUDE.md's
Gotchas).

## (b) al_analysis's removal commit

By hal's session in ~/src/al_analysis, after (c); `migration/cutover_hal.sh al-analysis` prints the
commands with the current list of tests:
- `git pull --rebase`; `git rm -r scripts/claude_pair_remapping_experiments`; `git rm` the queue's
  tests in test/ (every test_*.py naming the folder: the 18 copied here at 699e470, and any added
  later and copied by (a)). `test/test_model_provenance.py` stays (it tests mb_model);
- README.md:163: "The experiment harness (`scripts/claude_pair_remapping_experiments`)" becomes "The
  KC-remapping experiment harness (ejhonglab/model_experiment_queue since 2026-11-11)";
- al_analysis/mb_model.py, `olfsysm_checkout_path`'s docstring: "(see
  scripts/claude_pair_remapping_experiments/join_queue.sh)" becomes "(see join_queue.sh in
  ejhonglab/model_experiment_queue)";
- .gitignore:56, the comment "(the experiments folder's .claude/settings.json is shared)" may become
  "(a checkout's .claude/settings.json is shared)"; the two patterns under it stay;
- `QUICK=1 ./venv/bin/python -m pytest test/test_model_provenance.py`, commit, push.

The old folder keeps only ignored files afterwards (caches, and anything (c) left because it existed
on both sides). Delete it when nothing runs there.

## (d) Tom

- After his force push of the rewritten main (from_remy/ and apl_pair_analysis/ removed), and only
  then: `./github_rules.sh` from this repository. It creates the three rulesets the new repository
  lacks ("no deleting or force push" on main, nobody exempt; non-default branches and tags for
  admins only) and makes its Actions token read-only. Run it before the force push and the force
  push is refused. `./github_rules.sh --check` shows the state at any time.
- Commit and push the dotfiles change from (c), so atlas gets it with a pull.
- After viki re-runs join_queue.sh (e): `./approve_push_keys.sh viki`, comparing the one new
  fingerprint (model_experiment_queue) with the one Yang sends. Its three other keys are already in
  place.
- Start atlas's session(s) in ~/src/model_experiment_queue once atlas has done (e).

## (e) atlas and viki (chat.md text)

> For atlas (Tom's host, ~/src layout), once your runs from the old folder have finished: commit
> and push anything left in the old folder of al_analysis; stop the chat watcher; then
> `git clone https://github.com/ejhonglab/model_experiment_queue ~/src/model_experiment_queue`,
> `cd ~/src/model_experiment_queue`, `python3.14 -m venv venv`, `./venv/bin/pip install -e
> ../al_analysis[dev]`, `./venv/bin/pip install -e ../hong2p[all]`, and when no process has olfsysm
> loaded `./venv/bin/pip install -e ../olfsysm`. Move `outputs/`, `from_yang/` (if present) and
> `.fable_consent.local` from ~/src/al_analysis/scripts/claude_pair_remapping_experiments with `mv`
> (same filesystem). Pull ~/src/dotfiles once Tom has pushed the hook change, then pull al_analysis
> after its removal commit. Tom starts atlas's sessions in ~/src/model_experiment_queue.
>
> For viki (Yang's machine, ~/src/kc_queue), once your runs from the old folder have finished:
> commit and push anything left in the old folder, quit the session, then in a terminal
> `git clone https://github.com/ejhonglab/model_experiment_queue ~/src/kc_queue/model_experiment_queue`
> and `~/src/kc_queue/model_experiment_queue/join_queue.sh`. It moves outputs/, from_yang/ and
> .fable_consent.local over (and refuses while anything still runs in the old folder), builds the
> new venv in ~/src/kc_queue/model_experiment_queue/venv, and makes one new key. Send Tom its
> fingerprint. Then `git -C ~/src/kc_queue/al_analysis pull` and start the session with
> `cd ~/src/kc_queue/model_experiment_queue && claude` (accept the trust prompt). Delete
> ~/src/kc_queue/al_analysis/venv once the new session works.

## (f) The shared directory's live documents

Update the path mentions; the old ones are history once the move is done. Line numbers are from
2026-11-11 12:45.
- experiment_priorities.md:3: "Scheme: `al_analysis/scripts/claude_pair_remapping_experiments/CLAUDE.md`"
  becomes "Scheme: CLAUDE.md of ejhonglab/model_experiment_queue";
- consent/fable/README.md:5;
- hal_local_data.md:23, 29, 71 and atlas_local_data.md:40, 59: outputs/ and from_yang/ are now
  under ~/src/model_experiment_queue; viki_local_data.md:18: under
  ~/src/kc_queue/model_experiment_queue;
- notes/remy-summary.md:4 (from_remy/ is in the shared directory now) and notes/yang-summary.md:7
  (outputs/_scratch/yang/ moved with outputs/);
- E7_confound_controls.md:5 ("both in `scripts/claude_pair_remapping_experiments/`");
- the generated data READMEs naming the old folder: data/yang_apl16_kc_noise_units/README.md,
  README_unnormalized.md and olfdb/README_olfdb*.md, data/yang_weak_odor_kc/README.md,
  data/yang_apl16_kc_traces/README.md. Regenerate them with their scripts at the next run (the
  generators now name this repository), or edit the path text by hand.
  data/yang_apl16_kc_olfdb/README.md is hand-written: edit it.
Provenance JSONs and data/model_outputs archives are records of the past: leave them.

## Tests

From this repository's root, with the new venv (or, before it exists, al_analysis's venv as
`/home/tom/src/al_analysis/venv/bin/python`):
```
QUICK=1 OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -m pytest
./venv/bin/python -m pytest monitor-server
```
and al_analysis's `QUICK=1 ./venv/bin/python -m pytest test/test_model_provenance.py` from its
root.
