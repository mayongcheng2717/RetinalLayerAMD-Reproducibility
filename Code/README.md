# Code Records

The `_v6` filename suffixes are preserved historical execution names. The current release supports Manuscript v7; no additional scientific analysis is claimed from repackaging.

- `Review_Revision_Continuation_Pipeline.py`: authoritative user-local continuation pipeline, copied byte-for-byte.
- `Consolidated_Analysis_v6.py`: SMR global correction, early-AMD MVMR reconciliation, overlap scenarios, and unified reverse MR.
- `Robust_MVMR_v6.py`: MVMR-Median and MVMR-Egger sensitivity analyses.
- `Independent_Numerical_QA_v6.py`: 149 independent numerical checks using separate implementations.
- `Build_Figures_v6.py`: final figure and per-figure source-data construction, including direct vector-PDF-to-SVG conversion for retained Figure 2.
- `Build_Supplementary_Tables_v6.py`: final 50-sheet workbook construction.
- `Finalize_Supplementary_Tables_v6.py`: post-recalculation workbook-view restoration.
- `Core_Analysis_Requirements.txt`: original core dependency record returned with the local pipeline.
- `Verify_Archive_Integrity.py`: portable integrity checker for this archive.

The scientific scripts are preserved verbatim, including original absolute path constants. This avoids presenting edited code as the exact executed record. To rerun elsewhere, copy the relevant script and change only its root paths after reviewing the required inputs. Package integrity verification is separately portable.

- `Correct_Figure_1_Sample_Size_Label.py`: portable reproduction of the documented Figure 1 text-only correction; no numerical totals are altered.
- `Figure_1_Correction_Requirements.txt`: dependencies for that image-correction utility.

- `Update_Supplementary_Workbook_Repository_Note.py`: idempotent README-only workbook patch for use after rerunning the historical table builder. It leaves formulas, cached values, and style parts untouched.
