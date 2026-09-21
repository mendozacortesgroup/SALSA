# Changelog

All notable changes to SALSA (Substitution Approximation evoLutionary Search and
Ab-initio) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-09-15

First public release: the SALSA pipeline as used in the 2023 paper
([arXiv:2310.00118](https://arxiv.org/abs/2310.00118)), cleaned for
distribution. Covers stages 1-4 - ionic substitution, property approximation,
USPEX structure prediction and CRYSTAL setup. Stage 5, hybrid-DFT convex-hull
phase stability, follows in v2.0.0.

Every behavioural change below was verified by running the code rather than by
reading it, including against real SLURM on MSU HPCC where the scheduler's own
behaviour is what is under test.

### Added

- `setup_env.sh`, which writes a marked `SALSA_DIR` (and optional
  `VASP_PP_PATH`) block to `~/.bashrc`. Re-running replaces that block instead
  of appending a duplicate. `--print` previews without writing; `RC_FILE`
  redirects the target file.
- `requirements.txt` with pinned versions, replacing a prose `pip install` line
  that listed the standard-library module `zipfile` and omitted two packages the
  code actually imports. Requires Python >= 3.10.
- `LICENSE` (MIT) and `CITATION.cff`, with the 2023 paper as the preferred
  citation.
- `CHANGELOG.md`, `.gitignore` and `.gitattributes`.
- `3. USPEX_setup/POTCARs/README.md` documenting the pseudopotential policy.

### Fixed

- **`checkStatus_local.py` treated any scheduler failure as "job finished".**
  A bare `except: return True` meant that if `squeue` was missing or erroring,
  USPEX was told the job had completed and advanced past results that were
  never produced. It now raises instead. Note the distinction that matters
  here: once a job ages past SLURM's `MinJobAge`, `squeue -j <id>` exits
  non-zero with `Invalid job id specified`, and that *is* the ordinary
  completion signal - it is reported as finished, while any other failure
  raises. Both branches were exercised against SLURM 26.05.3.
- **`compile_POTCARs_locally.sh` copied nothing while reporting success.**
  Unquoted variables word-split on the space in `3. USPEX_setup`, so every
  `cp` failed silently, all 118 elements printed "Copied", and the script
  exited 0 - while `mkdir -p` created stray directories relative to the
  caller's working directory. Variables are quoted, a failed copy is now fatal,
  and the script reports counts and exits non-zero when it finds nothing.
- **Stale `Scripts/` directory layout.** Seven call sites across
  `setup_project.sh`, `SALSA_WorkFlow_Manager.py`, `print_bandgap_CRYSTAL.sh`
  and `instantiate_Step{1,2}.sh` referenced a `Scripts/` directory that does
  not exist in the repository, so project setup and the CRYSTAL helpers failed
  outright. They now point at the real locations (`data/`, `banner.txt`,
  `4. CRYSTAL_setup/CRYSTAL_Scripts/`).
- **`$SALSA_CIR` typo** in `instantiate_Step{1,2}.sh` left the basis-set path
  empty, and the path it built referenced a `ReferenceFiles/` layout that no
  longer exists.
- **`sed -i 's/CALCULATION_NAME/$calcname/g'`** was single-quoted, so it wrote
  the literal text `$calcname` into every generated submission script instead
  of the calculation name.
- **Hardcoded absolute path** `/mnt/home/staf6068/SALSA/...` in
  `submitJob_local.py`, replaced by `SALSA_DIR` with a path-relative fallback.
- **`#!/usr/bin/env python`** in the two CRYSTAL helper scripts, which fails on
  systems that ship only `python3`. Both are now `python3` and executable, and
  their callers invoke them through `python3` explicitly.
- **README** corrected throughout: the Python range said 3.6-3.9 when the pinned
  dependencies require >= 3.10; the manual `export SALSA_DIR = "..."` example
  was not valid bash; and the directory tree listed files that no longer exist.

### Removed

- 89 licensed VASP POTCAR files, purged from the full history with
  `git filter-repo`. These are proprietary and cannot be redistributed; supply
  your own and point `VASP_PP_PATH` at them. They are needed only for stage 3.
- `Projects/`, 626 files and roughly 90 MB of `CHGCAR`, `WAVECAR` and further
  POTCAR data belonging to a specific study rather than to the tool. The
  directory is now created at runtime by `setup_project.sh` and is gitignored.
- `.DS_Store` and `__pycache__` artefacts.

Together these reduced the tracked tree from 825 to 138 files and the packed
repository from about 15 MB to 1.1 MB.
