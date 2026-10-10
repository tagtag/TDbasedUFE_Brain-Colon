# Validation performed during assembly

Date: 2026-10-09.

## Completed checks

The supplied reference files were recalculated in the available Python/SciPy environment, independently of the new R driver:

- 75 unique mouse identifiers; the score and behavioral files use the same individuals.
- 16,080 unique background gene IDs and 1,577 unique selected gene IDs.
- All 36 ordinary Pearson correlations, P-values, and BH-adjusted P-values match the archived analysis.
- All 36 correlations after categorical age adjustment, with 69 degrees of freedom, match the archived analysis.
- All 32 ENS contingency tables and Fisher/BH results match. The 29 Table S4 sets were rebuilt from the original workbook using `coefD > 0` and `padjD < 0.05`; the three small published-list tests used the archived marker memberships.
- All 12 Boruta comparisons (six sources at two run-count thresholds) match. The union contains 506 genes, overlaps the tensor set at 64 genes, and leaves 1,513 genes.
- All 6,132 post-Boruta enrichment tests were recalculated from the frozen GMT files with the documented library-specific backgrounds and matched the archived results.
- The highlighted Table 7 overlaps are 31/31, 43/45, 28/29, 32/32, and 23/24.
- The archive contains 228 original Enrichr library exports, not 218.

Numeric comparisons used relative tolerance `1e-8` (with small absolute tolerance where appropriate). The archive integrity check verifies every packaged file against its source before delivery; the R checksum script verifies the reference snapshot after unpacking.

## Not executed here

An R runtime was unavailable in the authoring environment. R installation was blocked or unavailable in this environment, and the attempt was stopped. Consequently:

- No successful R parse or R execution is claimed.
- Neither the new R reference driver nor the full download/HOSVD workflow was executed here.
- The complete approximately 630 MB source archive and Ensembl GTF were not downloaded again solely to test this package.
- The original runtime package versions, HOSVD object, sigma-search object, and reusable mapping CSV were not supplied. The original R tensor and feature-selection functions are included; the mapping scope is reconstructed from the original mapping script.

The code was reviewed for path consistency, input column names, tensor ordering, age-adjustment degrees of freedom, enrichment backgrounds, and reference comparison logic. This review is not a substitute for a completed R run.

## Checks to run in R

```sh
Rscript scripts/install_dependencies.R
Rscript scripts/verify_bundle.R
Rscript run_analysis.R reference
sh run_sequential.sh full
```

`reference` mode explicitly recalculates downstream statistics. `full` mode rebuilds the molecular tensor and checks exact background/selected gene identity and numeric molecular scores before a successful final validation. Mismatches stop the pipeline and leave diagnostics; they never trigger automatic parameter tuning. If available, supply the original transcript map at `data/transcript_to_gene_map.csv`.

The included single-job GitHub Actions workflow performs an R syntax check, checksum verification, and reference-mode calculation after publication. It does not certify a full HOSVD rerun; that requires the `full` command and its separate validation report.

Retain `results/full/validation.csv` and `results/full/sessionInfo.txt` with the eventual manuscript release. A successful future run should be documented as such, without attributing that validation to the package-assembly session.
