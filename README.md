# Cognitive aging: hippocampal and colonic multi-omics tensor analysis

Reproducibility materials for **Unsupervised tensor integration of hippocampal and colonic multi-omics reveals enteric neural signatures in cognitive aging**.

This repository provides an R workflow for the common-gene tensor analysis, behavioral correlations, enteric marker enrichment, and comparison with published Boruta-selected features. All analysis stages run sequentially.

## Quick start

Install R (4.3 or newer recommended), unzip this repository, and open a terminal in its root directory:

```sh
Rscript scripts/install_dependencies.R
Rscript scripts/verify_bundle.R
Rscript run_analysis.R reference
```

`reference` mode requires `readxl`; the statistical calculations otherwise use base R. It uses the archived molecular scores and selected genes to **recalculate** all 36 ordinary correlations, all 36 age-adjusted correlations, all 32 ENS marker tests, Boruta overlaps, and all 6,132 Boruta-excluded enrichment tests. It does not rerun HOSVD. Tables 1 and 2 are reconstructed from archived Enrichr exports, not represented as newly calculated tests.

For the full reconstruction from the public molecular matrices:

```sh
sh run_sequential.sh full
```

The shell launcher sets single-thread BLAS/OpenMP limits before R starts. It invokes `Rscript run_analysis.R full`. On Windows, use R with a single-thread BLAS or set the corresponding environment variables before starting R, then run:

```powershell
$env:OMP_NUM_THREADS="1"
$env:OPENBLAS_NUM_THREADS="1"
$env:MKL_NUM_THREADS="1"
Rscript run_analysis.R full
```

`RhpcBLASctl`, when installed, also limits threads from R. No script starts workers or parallel jobs. Dependency installation uses `Ncpus=1`.

## Run modes

| Command | Purpose | Network access |
|---|---|---|
| `Rscript run_analysis.R reference` | Recalculate downstream statistics from frozen inputs | None |
| `Rscript run_analysis.R download` | Download/extract original inputs and the pinned Ensembl GTF | Yes, for missing files |
| `Rscript run_analysis.R full` | Download if needed, rebuild tensor and HOSVD, select genes, run downstream statistics | Yes, for missing files |
| `Rscript scripts/refresh_enrichr.R` | Optional new Enrichr submission after a full run | Yes |

Run commands from the repository root. Results are written to separate `results/reference/`, `results/full/`, or `results/enrichr_live/` directories. Existing downloads are reused; there are no automatic download retries. A failed download leaves a `.part` file and never becomes a valid cached input. Full analysis recomputes its outputs; it does not silently substitute reference scores when reconstruction fails.

## Data acquisition

### Original cognitive-aging data

Flor et al., *Multi-omics analysis highlights the link of aging-related cognitive decline with systemic inflammation and alterations of tissue-maintenance*:

- Study: https://doi.org/10.1101/2025.07.13.662751
- Data: https://doi.org/10.5281/zenodo.15849341
- Download: https://zenodo.org/records/15849341/files/ScriptsAndFilesCognitiveAgingManuscript.zip?download=1

The downloader fetches the approximately 630 MB ZIP and extracts only the required processed matrices, behavioral table, original annotation, Boruta summary CSVs, and source data license. It does not rerun raw RNA-seq or methylation processing, or the original Boruta model fitting. Source matrices are the published VST RNA values and mean promoter methylation values.

| Tissue/modality | Source filename |
|---|---|
| Hippocampal RNA | `varianceStabilizedRNASeq_Brain.tsv` |
| Colonic RNA | `varianceStabilizedRNASeq_Colon.tsv` |
| Hippocampal methylation | `brainMeanPromoterMethylationNoMissingValues.tsv` |
| Colonic methylation | `intestineMeanPromoterMethylationNoMissingValues.tsv` |
| Behavior | `df_allCognitionData_2025-02-17.tsv` |

To work offline, place the original ZIP at `data/ScriptsAndFilesCognitiveAgingManuscript.zip` and the following GTF at `data/Mus_musculus.GRCm39.116.gtf.gz` before running `full`.

### Ensembl annotation

The full workflow downloads the release-116 mouse GTF:

https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz

Transcript-to-gene conversion for tensor construction and gene-to-symbol conversion for Enrichr use this pinned annotation. ENS marker matching and the Boruta comparison use the **original study annotation**, `Ensembl2GeneSymbol.csv`, matching the manuscript analyses. These two annotation uses are deliberately distinct. Ensembl identifiers, not symbols, are the counting unit for ENS tests.

### Enteric marker data

Drokhlyansky et al., *The human and mouse enteric nervous system at single-cell resolution*, Cell (2020):

- https://doi.org/10.1016/j.cell.2020.08.003
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8358727/

The supplied `reference/ens/Drokhlyansky_Table_S4.xlsx` is the supplementary Table S4 used for the analysis. It is included so that the workflow does not depend on changing supplementary-download URLs. The R code reads the `Mouse Colon All` sheet and recalculates marker sets directly from it.

### Enrichr snapshots

`reference/enrichr/` contains the original 228 library exports and submission metadata. The saved run submitted 1,574 unique symbols representing 1,577 gene IDs. Tables 1 and 2 preserve the saved export order. The archive contains **228** successful library exports; it should not be described as 218 libraries.

The two frozen GMT libraries used for Table 7 are in `reference/gene_sets/`:

- `GO_Biological_Process_2026.gmt`
- `PanglaoDB_Augmented_2021.gmt`

Their source endpoint is `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LIBRARY_NAME`. The default pipeline uses the frozen files because a current download may change membership and P-values even when the library name is unchanged. The optional live Enrichr script submits the selected symbols sequentially to the same 228 library names and writes a separate output directory; current results are not guaranteed to match the manuscript.

## Analysis specification

1. Trim whitespace from sample identifiers. Match all four molecular matrices and the behavioral table. Require finite values in all nine behavioral columns and apply the recorded exclusion of `30M_16N`.
2. There are 81 common individuals and 78 eligible individuals. Use the exact frozen 75-individual publication map, with 15 individuals per age group and **no replacement or duplicate mice**. The original map was generated with seed `20261005` by sampling sorted eligible IDs without replacement. The default workflow uses the archived map to avoid RNG/locale-dependent changes in replicate ordering. Replicate positions do not denote longitudinally repeated mice across ages.
3. Reconstruct the scope of the reused transcript-to-gene mapping: the original mapping file retained uniquely mapped transcripts from the finite methylation transcript intersection across the two tissues. Read Ensembl release 116, restrict it to that same transcript universe, and take the arithmetic mean of mapped promoter values within each gene. Retain the finite common gene intersection across all four layers: 16,080 genes, ordered by gene ID. If the original `overlap_1121/transcript_to_gene_map.csv` is available, placing it at `data/transcript_to_gene_map.csv` uses it directly. This is an annotation-universe restriction; no preliminary SVD compression is performed in this workflow.
4. For each tissue, modality, and individual, subtract the mean across those genes and divide by their sample standard deviation (`sample_zscore`, denominator `N-1`). No preliminary per-modality SVD compression is used.
5. Construct a `2 x 2 x 5 x 15 x 16080` tensor, with tissue, modality, age, replicate, and gene modes. Perform HOSVD with ranks `c(2,2,5,15,300)`. The gene mode retains all possible nonzero directions; the other-mode product is 300, not 300 independent mice.
6. Inspect the core slice at `(1,1,2,1,:)`. Select the gene component with the largest absolute core entry, expected to be component 5. This choice does not use cognitive scores. Save the core slice and mode-factor plots for inspection.
7. Apply the original mean-centered Gaussian null and 100-bin histogram-flatness selection procedure to the unit-norm gene loading vector. Use BH `<= 0.01`, an initial log-scale sigma range of `sd(vector) * c(0.01,100)`, 101 coarse points plus the initial SD, two 61-point refinement rounds, and a minimum background fraction of 0.5. Select 1,577 gene IDs. Full settings, all tested genes, search evaluations, and diagnostic plots are saved. The code stops on a gene-set mismatch rather than changing thresholds to force agreement.
8. Project each normalized molecular profile onto the selected gene loading vector. SVD signs are arbitrary; score orientation is matched to the archived **molecular** scores, without using behavioral outcomes. Save that sign choice. Figures may have an equivalent global sign reversal in other implementations.
9. Compute two-sided Pearson tests for 4 molecular scores x 9 behavioral measures (36 tests), with BH correction across all 36. For age adjustment, regress each score on a five-level categorical age factor and correlate residuals. Use `df = 75 - 5 - 1 = 69` for the partial-correlation test and apply BH separately across the 36 adjusted tests.
10. For Table S4, retain `coefD > 0` and `padjD < 0.05`, using all 29 cell categories. Add the published neuronal, glial, and union signature lists. Match case-sensitive symbols to the original study annotation and the 16,080-gene background. Test all 32 sets with one-sided Fisher tests and a joint BH correction. Keep the literal unresolved `Cacna2d` entry; do not replace it with a guessed gene.
11. Reconstruct the published Boruta gene sets from the six relevant feature-count summaries (five individual molecular layers and the combined analysis). The broadest set is their union at count >= 1, including liver. A single missing transcript mapping, `ENSMUST00000239425 -> ENSMUSG00000071984` (`Fndc1`), uses the archived Ensembl lookup. The lookup response is included, not fetched dynamically.
12. Remove the broadest Boruta set from the tensor-selected genes: 1,577 - 64 = 1,513. For each GMT library, remove Boruta symbols from the union of annotated genes, restrict the query and each term to that background, and calculate upper-tail hypergeometric P-values. Apply BH separately to all 5,954 GO terms and all 178 PanglaoDB terms, including zero-overlap terms.

The ENS analysis uses the measured common-gene background. Table 7 uses library-specific annotation backgrounds after Boruta exclusion. These backgrounds must not be interchanged. Table 7's overlap denominator is the number of selected term-associated genes **before exclusion**, not the term's total size.

## Expected results and manuscript outputs

| Output | Expected result / role |
|---|---|
| `Table1_Orphanet_Augmented_2021.csv` | Top 10 archived Enrichr rows |
| `Table2_GTEx_Aging_Signatures_2021.csv` | Top 10 archived Enrichr rows |
| `Table3_correlations.csv` | Significant correlations at BH < 0.05 |
| `Table4_ENS.csv` | Neuron 1, Neuron 2, glia, fibroblast |
| `Table5_source_samples.csv` (full mode) | Hippocampal RNA: 82; other three matrices and behavior: 83 each |
| `Table6_sample_selection.csv` (full mode) | Common 81; eligible 78; selected 75 |
| `Table7_TD_only_enrichment.csv` | Five highlighted synaptic / enteric terms |
| `tensor/` and `selection/` (full mode) | Tensor, HOSVD, all gene tests, sigma optimization diagnostics |
| `Figure1_mode_factors.pdf` (full mode) | Tissue/modality/age/replicate factors |
| `Figure2_gene_loadings.pdf` (full mode) | Gene loadings with selected genes highlighted |
| `cognition_scatterplots.pdf` | Overall and age-adjusted molecular–cognition scatterplots |
| `validation.csv` | Numeric comparisons with the archived analysis |
| `sessionInfo.txt` | Software/platform versions from your actual run |

The plots reproduce the numerical analysis; graphic styling is not intended to be pixel-identical to the manuscript.

- Colonic methylation vs overall cognition: `r = 0.5888819697`, BH `P = 9.8713734e-7`.
- No age-adjusted test passes BH < 0.05.
- Neuron 1: 38/76, fold enrichment 5.0983, BH `P = 1.6862e-17`.
- Glia: 31/68, fold enrichment 4.6484, BH `P = 3.9233e-13`.
- Boruta union: 506 genes; overlap: 64; tensor-only: 1,513.
- Boruta-excluded GO background/query: 15,287 / 730; PanglaoDB: 6,479 / 525.

`reference/` files are archived inputs and expected results, not outputs of a new full run. Each run writes freshly calculated results outside that directory. R's Fisher conditional odds-ratio estimate can differ from the original Python sample odds ratio; both are exported, while the P-values and BH values use equivalent tests.

## Validation status and reproducibility limits

See [docs/VALIDATION.md](docs/VALIDATION.md) for the checks performed when this package was assembled. The assembled R workflow has not been executed end-to-end in the authoring environment, which lacks an R runtime. Accordingly, this repository is not labeled as a verified full rerun. The original R construction/selection functions are retained, and the full run checks the background, selected gene identities, molecular scores, and downstream statistics against the archived analysis.

The original complete R session information, HOSVD object, sigma-search object, and reusable transcript-mapping CSV were not available in the supplied files. The mapping scope is reconstructed from the original scripts; an original mapping CSV can instead be supplied as described above. No historical package versions or optimized sigma are invented. Numerical linear algebra and a discrete histogram objective can produce small platform-dependent differences; any mismatch is reported, not silently corrected. Record the generated session information after a successful full run and retain the results with the manuscript release.

This analysis supports age-associated molecular and marker signatures. It does not establish an age-independent association, enteric specificity relative to the CNS, or causal gut-to-brain effects.

## Repository layout

- `run_analysis.R`: single entry point.
- `R/`: input preparation, tensor construction, feature selection, statistics, and validation.
- `scripts/`: explicit dependency installation, reference checksums, optional live Enrichr run.
- `reference/`: frozen sample map, gene sets, study annotation, marker spreadsheet, original exports, and numerical reference results.
- `data/`: downloaded files; excluded from Git.
- `results/`: generated output; excluded from Git.
- `docs/`: provenance and validation notes.

## Publishing on GitHub

Create an empty GitHub repository, then run these commands from this directory, replacing `YOUR_ACCOUNT` and `YOUR_REPOSITORY`:

```sh
git init
git add .
git commit -m "Add cognitive-aging tensor reproducibility workflow"
git branch -M main
git remote add origin https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
```

The `.gitignore` excludes downloaded ZIP/GTF files and generated results. The included reference files are intentionally tracked. For a publication release, attach the successful-run `sessionInfo.txt` and `validation.csv`, choose a code license, and archive the tagged repository with a DOI. See [DATA_SOURCES.md](DATA_SOURCES.md) for third-party data attribution. This package does not assign a new license to the original datasets or supplementary material.
