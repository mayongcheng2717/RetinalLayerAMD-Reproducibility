# RetinalLayerAMD-Reproducibility

Release: **v1.0.0, 27 October 2026**

## Scope

This repository preserves the actual returned results, consolidated numerical analyses, final parameters, random-seed records, software environments, figure source data, plotting scripts, final figures, supplementary files, and validation evidence supporting **Manuscript v7**. The `_v6` suffixes retained in several code filenames are historical execution filenames; they do not identify the current manuscript version, and no new scientific analysis is implied by this release.

The archive is an execution record and source-data package, not a fresh-download one-click installer. Large public GWAS, eQTL, pQTL, retinal-atlas, and Visium raw downloads are not redistributed. Their identifiers, source URLs, local audits, and available checksums are supplied. The submission DOCX is distributed separately in the journal submission package.

A public repository URL is assigned when the authors upload this release to GitHub. Insert that same URL at both marked positions in the manuscript's **Availability of data and materials** statement. **No DOI has been assigned.**

## Directory guide

- `Code`: exact analysis and build scripts, the portable archive checker, the documented Figure 1 text-correction script, and an idempotent supplementary-workbook README updater. Original execution-path constants are preserved for provenance; adjust root constants before rerunning elsewhere.
- `Config`: final parameters, seed maps, trait definitions, panel annotation, and the data-source catalog. The superseded MR-PRESSO configuration with 3,000 simulations and seed 42 is not a final parameter file.
- `Environment`: local, full54/scDRS, robust-MVMR, and finalization environment records.
- `Inputs`: compact legacy inputs, the legacy supplementary workbook used for Figure 6C and 6D, and retained Figure 2 source files.
- `Results/Returned`: returned MR, MVMR, SMR/colocalization, Visium, sample-size, environment, pQTL, scDRS, and diagnostic outputs.
- `Results/Consolidated`: independently reconciled SMR, MVMR, overlap, and reverse-MR outputs and manifests.
- `Source_Data`: 37 CSV files supplying every finalized quantitative figure panel, plus the Figure 1 reconstruction prompt.
- `Figures`: ten final white-background PNGs; nine vector PDFs and nine SVGs for Figures 2–6 and Supplementary Figures S1–S4. Figure 1 is raster-only.
- `Figures/Source`: the untouched author-supplied Figure 1, the corrected native-resolution image, and the documented provenance record.
- `Additional_Files`: the final supplementary workbook and four-page supplementary-figure PDF.
- `Provenance`: authoritative-input checksums, figure-to-source mapping, key numeric spot checks, evidence tiers, three-level data provenance, and image-preparation records.
- `QA`: independent numerical QA, current v7 manuscript checks, workbook verification, figure QA, citation checks, outstanding analyses, and release status.

## Reproduction and verification

The scientific execution order is recorded by the scripts themselves: the authoritative local continuation pipeline, consolidated reconciliation, robust MVMR, independent numerical QA, supplementary-table construction, and figure construction. Exact scientific scripts retain original absolute paths so the executed code is not silently rewritten for packaging.

Run integrity verification from the extracted repository root:

```bash
python3 Code/Verify_Archive_Integrity.py
```

The checker uses the script location to find the repository, verifies every listed SHA256, detects missing or unlisted files, and checks for non-ASCII paths and empty directories. `MANIFEST.csv` and `SHA256SUMS.txt` are the two inventory files excluded from the checksum list to avoid circular hashes.

The historical table builder retains its original pre-repository README string as part of the execution record. After rerunning it, use `Code/Update_Supplementary_Workbook_Repository_Note.py` to restore the current wording without changing formulas, caches, or styles. See `Provenance/Historical_Code_Notes.json`.

`Code/Correct_Figure_1_Sample_Size_Label.py` reproduces the documented Figure 1 text-only correction. It requires NumPy, OpenCV, Pillow, and Liberation Sans. Set `FIGURE1_FONT_PATH` if the font is installed elsewhere. The exact font hash and pixel-edit bounds are recorded in `Provenance/Figure_1_Correction.json`.

## Figure 1 correction

The author-supplied artwork labeled the combined six-disease totals as “early AMD: 54,082 cases, 1,306,001 controls.” Those totals describe all six disease datasets. Winkler early AMD (GCST010723) contains 14,034 cases and 91,214 controls. The final figure therefore reads “6 traits; 54,082 cases, 1,306,001 controls.” Only the erroneous attribution was removed; the combined numerical totals were not changed. The untouched original, corrected native image, final 300-dpi PNG, and hashes are retained.

## Interpretation boundaries

- The RPE early-AMD result is an rs10737680-dependent and estimator-inconsistent sensitivity finding.
- INL and OPL are not reported as two proven independent causal effects.
- Colocalization supports shared association signals under the specified single-causal-variant model. It does not establish layer-thickness mediation of disease.
- The TNXB indirect effect has P=0.084807; mediation is not confirmed.
- scDRS sensitivity analyses remain within one six-donor atlas.
- Visium analyses are descriptive because disease status is completely confounded with donor identity and matched-random tests were null. ENSG00000258311 was not detected under the available alias.
- No non-overlapping early AMD replication, retinal-specific QTL validation, independent-donor spatial validation, or functional perturbation is claimed.
- Outstanding analyses are listed in `QA/Outstanding_Analysis_Status.csv`; their absence is not counted as a failed completed-QA check.

## Figure editing

PDF is the submission-first vector format when available. SVG text is outlined to avoid missing fonts and is no longer editable as text. Appearance edits are permitted, but data elements such as point positions, error bars, axis ranges, heatmap values, and bar heights must only be changed through the analysis and plotting code.

## License and third-party terms

`LICENSE` applies the MIT License to original project code. Third-party software, public datasets, derived source data, retained artwork, and fonts remain subject to their original licenses and terms. The license does not grant rights to restricted underlying data and does not replace source-specific citation requirements.
