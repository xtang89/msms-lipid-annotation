# msms-lipid-annotation

A Claude skill that annotates lipid MS/MS spectra from Thermo `.raw` files with diagnostic-fragment rules,
records your decision for each spectrum, and plots the annotated spectra in Origin.

## What it does
1. Reads a Thermo `.raw` file and groups the MS/MS scans into windows
   (same polarity, activation, collision energy and precursor m/z).
2. Averages all scans of each window and centroids the result.
3. Finds the precursor m/z (MS1 peak, else precursor peak in MS2, else filter value).
4. Finds candidate lipids in a precursor database made from LIPID MAPS COMP_DB.
5. Matches the fragments in the rule files and gives a level for each candidate:

   | Level | Condition | Meaning |
   |---|---|---|
   | high | ≥ 2 diagnostic fragments, or 1 diagnostic + ≥ 1 supporting | MS/MS confirmed |
   | medium | 1 diagnostic fragment | MS/MS confirmed |
   | low | supporting fragments only | MS1 annotation |
   | ms1_only | no fragment | MS1 annotation |

   TG, DG and MG are confirmed by chain fragments. A chain combination (for example PC 34:1 (16:0_18:1))
   is given only when fragments for all its chains are present. Isobaric candidates at the same level are all reported.
6. Claude explains the result. You confirm, change or reject each annotation. The decision is saved.
7. Origin plots (`.opju` and `.tif`) of the confirmed spectra.

The scripts do all mass calculations. Claude does not calculate m/z values.

## Lipid classes with rules
- Positive: PC, PE, SM, PS, LPC, LPE ([M+H]+, and Li/Na/K adducts where confirmed); TG, DG, MG (Na, NH4, Li, K where confirmed).
- Negative: PC and LPC (formate and acetate adducts), PE, LPE, PI, LPI, PS, LPS, PG ([M-H]-); PGE2/PGD2/PGH2, PGF2a.
- Other classes in LIPID MAPS COMP_DB get an accurate-mass annotation only.

See `rules/rules_review.md` for each rule and `rules/references.md` for the literature.

## Requirements
- Windows with Python 3.9 or later (Anaconda is fine).
- Python packages: `pip install -r requirements.txt`
  (numpy, pyyaml, openpyxl, fisher_py, native-fisher-py, originpro).
- OriginPro (2021 or later) for the plot step only. Steps 1–4 work without Origin.
- Thermo `.raw` files. Other formats are not supported.

## Install
**Claude Code (recommended; all steps work, including Origin plots)**
```
git clone https://github.com/xtang89/msms-lipid-annotation.git %USERPROFILE%\.claude\skills\msms-lipid-annotation
pip install -r %USERPROFILE%\.claude\skills\msms-lipid-annotation\requirements.txt
```
Claude Code must run on the same Windows computer as Origin.

**Claude app (claude.ai or desktop chat)**
Download the repository as a ZIP file and upload it as a skill in the Claude app settings.
In the chat, the scripts run in a Linux workspace. Thus the annotation works, but the Origin plot step does not.

## Use
Ask Claude, for example:
> Annotate the MS/MS data in C:\data\brain_pos.raw and C:\data\brain_neg.raw

Claude checks the environment (and asks before it installs packages), runs the annotation, explains the result,
records your decisions and makes the plots.

You can also run the scripts yourself:
```
python scripts/check_env.py
python scripts/annotate.py C:\data\brain_pos.raw C:\data\brain_neg.raw
python scripts/set_status.py C:\data\brain_pos_annotation all confirmed
python scripts/plot_origin.py C:\data\brain_pos_annotation
```

## Output
For each `.raw` file, a folder `<rawfile>_annotation/` next to the file:
- `results.json`: all windows, candidates, fragments and decisions (for Claude)
- `results.xlsx`: the same result as tables (for people)
- `plots/<window_id>.opju` and `plots/<window_id>.tif`

## Settings
`config.yaml`: mass tolerance (10 ppm), intensity threshold (1 % of base peak), adducts, plot size and fonts.
After you change a rule file, run `python scripts/build_precursor_db.py` to make the precursor database again.

## Known limits
- Scans are averaged with the Thermo averager through `fisher_py` (needs .NET; Windows has it).
  Without .NET (for example the Linux workspace of the Claude app), the script reads the file with `native-fisher-py`
  and uses its own averaging (peaks within the tolerance are combined; intensity = sum / number of scans).
- Ether lipids (O-, P-) are annotated at species level only.
- Isobaric species (for example PC 34:1 [M+Na]+ and PC 36:4 [M+H]+) are both reported when the fragments agree with both.
- The rules are for low-energy CID/HCD spectra.

## Data and license
Code: MIT license (see `LICENSE`).
`COMP_DB.tsv` and `data/precursor_db.csv` are derived from the LIPID MAPS® Computationally-generated Bulk Lipids database
(COMP_DB), licensed under CC BY 4.0. Please cite LIPID MAPS (https://www.lipidmaps.org) when you use them:
Conroy MJ et al. LIPID MAPS: update to databases and tools for the lipidomics community. Nucleic Acids Res 2024, 52, D1677–D1682.
