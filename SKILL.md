---
name: msms-lipid-annotation
description: Annotate on-tissue lipid MS/MS spectra from Thermo .raw files (averaged spectrum per precursor window) with diagnostic-fragment rules, record the user's decision, and plot the annotated spectra in Origin. Use when the user asks to annotate, check or plot lipid MS/MS data.
---

# Lipid MS/MS annotation

Version 0.7.0 (2026-11-22). Runs on Windows. The plot step needs OriginPro on the same computer.

The scripts do all calculations (masses, matching, levels). You (the LLM) run the steps, read
`results.json`, explain the result to the user, and record the user's decision.
**Do not calculate m/z values yourself. Use only the values in results.json.**

All paths below are relative to this skill folder.

## Folder
```
config.yaml                 settings (tolerances, threshold, adducts, plot)
rules/positive.yaml         fragment rules, positive mode
rules/negative.yaml         fragment rules, negative mode
rules/references.md         literature for the rules
rules/rules_review.md       rule notes: scoring, rules not included, example m/z values
data/precursor_db.csv       precursor m/z database (LIPID MAPS COMP_DB + rules)
scripts/check_env.py        step 0
scripts/annotate.py         step 2 (uses average.py, match.py, report.py, lipidmass.py)
scripts/set_status.py       step 4
scripts/plot_origin.py      step 5
scripts/build_precursor_db.py  only after a rule file changes
```

## Step 0. Check the environment
1. Run `python scripts/check_env.py` and read the JSON report.
2. Tell the user in short form: Python version, missing packages, Origin found or not.
3. If packages are missing, show the `install_command` from the report and **ask the user before you
   install anything**. If the user prefers a new environment, ask which tool (venv or conda) and its name.
4. `ready_for_annotation` must be true before step 2. If `ready_for_plot` is false, steps 1–4 still work;
   tell the user that the plot step needs Origin on Windows.

## Step 1. Get the input
Ask for the .raw file(s) or a folder of .raw files, if the user did not give them.
Tell the user the current adducts and tolerances from `config.yaml` in one line. Change `config.yaml` only if the user asks.

## Step 2. Run the annotation
`python scripts/annotate.py <file.raw or folder>`
For each .raw file the script writes `<rawfile>_annotation/results.json` and `results.xlsx` next to the file.

What the script does for each MS/MS window (one window = same polarity, activation and collision energy,
and precursor m/z within the tolerance; this also merges repeated selections in data-dependent files):
1. Averages all scans of the window and centroids the result with the Thermo averager (needs `fisher_py`). Without `fisher_py` (no .NET, e.g. Linux), the script reads the file with `native-fisher-py` and uses its own average.
2. Finds the precursor m/z: MS1 peak, else precursor peak in MS2, else filter value (2 decimals; wider tolerance).
3. Finds candidates in `data/precursor_db.csv` (polarity and adducts from config). No candidate → `no_match`.
   Species of classes (or adducts) without rules are reported by accurate mass only ("no rules").
4. Matches the rule fragments. Levels:
   - `high`: ≥2 diagnostic, or 1 diagnostic + ≥1 supporting → MS/MS confirmed
   - `medium`: 1 diagnostic → MS/MS confirmed
   - `low`: supporting only → MS1 annotation; MS/MS does not confirm the class
   - `ms1_only`: no fragment → MS1 annotation
   - TG, DG, MG (`class_rule: chains`): `high` if ≥2 chain fragments fit one chain combination
     (1 is enough if all chains of that combination are the same)
5. Chain combination (e.g. 16:0_18:1) is given only if fragments for all its chains are present.

## Step 3. Explain the result
Read `results.json`. For each window give one short line or a table row:
window, precursor m/z (and source), level, top annotation(s), the diagnostic/supporting fragments found,
chain combinations, other candidates, and large unexplained peaks if there are any.

Rules for your explanation:
- If two or more candidates have the same best level, report all of them. Do not choose one.
- If several chain combinations are complete, write that fragments for each were detected
  (e.g. "chain fragments for 16:0_18:1 and 18:0_16:1 detected"). Do **not** say which one is more abundant.
- Ether species (O-, P-) get species level only.
- Candidates marked "accurate mass only; no rules" are not MS/MS annotations. Say so.
- For `low` and `ms1_only`, say clearly that the MS/MS spectrum does not confirm the class.
- If the precursor m/z came from the filter (2 decimals), say so: the precursor match is less certain.
- Mark windows with large unexplained peaks for the user (possible mixture or a class not in the rules).

## Step 4. Record the user's decision
Ask the user to confirm, change or reject each window (they can confirm all at once).
- `python scripts/set_status.py <annotation_folder> <window_id> confirmed`
- `python scripts/set_status.py <annotation_folder> <window_id> changed --annotation "..." --comment "..."`
- `python scripts/set_status.py <annotation_folder> <window_id> rejected --comment "..."`
- `python scripts/set_status.py <annotation_folder> all confirmed` (all windows still pending)
The script rewrites `results.json` and `results.xlsx`.

## Step 5. Plot (Windows + Origin)
`python scripts/plot_origin.py <annotation_folder>` plots confirmed and changed windows
(`--all` plots all windows, `--windows id1,id2` plots selected windows,
`--confirmed-level` plots windows with level high or medium). Output: `plots/<window_id>.opju` and `plots/<window_id>.tif`.

The script controls Origin in the background (Origin does not open on the screen).
Plot design: sticks, x = m/z, y = relative intensity (%) from 0 to 110. Colors (Okabe–Ito):
unannotated gray, class fragments blue, chain fragments vermillion; ▼ marks the precursor. No legend.
Labels only on matched fragments, two lines: observed m/z (4 decimals) and the assignment
("-183.066" = neutral loss, "C5H15NO4P+" = fragment ion, "FA 18:1" = chain ion, "-FA 18:1" = chain loss).
Header above the frame: annotation, precursor m/z, activation and collision energy. Arial, axis text 9 pt,
labels 7 pt, TIF 600 dpi.
If the script cannot start Origin, tell the user that Origin must be installed on this Windows computer.

## Step 6. Finish
Tell the user where the folder is and how many windows were confirmed, changed and rejected.

## When the rules change
Edit `rules/*.yaml` (add the source to `rules/references.md` and the change to `rules/rules_review.md`),
then run `python scripts/build_precursor_db.py` to rebuild `data/precursor_db.csv`.
Options: `--include-oxidized` (oxidized species), `--no-ether` (remove ether species).
The script keeps only species that the chain range in the rule files can make.

## Use as a subskill
This folder is self-contained (own config, rules, data and scripts; relative paths only).
A parent skill can call the scripts above from the command line.
