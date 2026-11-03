# NMD3 functional ribosome supply — figure reproduction

Reproducible **R analysis driver** for the study *Restricted functional ribosome supply uncouples protein expression from transcript abundance in advanced prostate cancer* (figure layout dated **October 30, 2026**).

The script executes the computational panels in their **manuscript/PPT order**: **Figures 1–6, then Supplementary Figures S1–S8**. Individual plot PDFs are written and closed when each analysis finishes. Wet-lab blots, polysome traces, IHC images, hand-drawn schematics, and other panels that were not made from R inputs are **not regenerated**.

> **Important availability note:** The R code is shareable, but a full reproduction requires both public datasets and the authors' in-house RNA-seq/DIA inputs, as well as a few locally assembled public-cohort derivative files. This repository package does not include the original source data or in-house experimental inputs. Downloading the public study archives alone does not create the authors' derived `.xlsx` or `.RData` objects.

## Repository contents

```text
NMD3_Figures_1_to_S8.R       # Full analysis, Main followed by Supplementary
README.md                    # This guide
.gitignore                   # Excludes data and output
 data/                        # Create the input tree below (not version-controlled)
 results/                     # Created automatically for each run (not version-controlled)
```

## 1. Software

Use R with Bioconductor and CRAN packages required at the beginning of the script. The driver checks package availability before starting analysis.

- **Bioconductor (common):** `DESeq2`, `apeglm`, `limma`, `GEOquery`, `Biobase`, `GSVA`, `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `enrichplot`, `cBioPortalData`, `fgsea`.
- **CRAN (common):** `ashr`, `readxl`, `statmod`, `dplyr`, `tibble`, `tidyr`, `ggplot2`, `ggrepel`, `ggpubr`, `patchwork`, `pheatmap`, `gridExtra`, `msigdbr`, `janitor`, `readr`, `data.table`, `stringr`, `survival`, `survminer`, `VennDiagram`, `purrr`, `Seurat`, `forcats`, `ggforce`, `nlme`, `scales`.

Some packages above are distributed through Bioconductor rather than CRAN depending on the R/Bioconductor release. The **actual authoritative package list is `.nmd3_required_packages` in the R script**. Install the packages with `BiocManager::install()` or `install.packages()` as appropriate. Do not use an unrelated historical `.RData` or `.RDS` to substitute for missing input files.

**Versioning:** GSEA gene-set membership can vary with `msigdbr`/MSigDB releases. The analysis originated with human MSigDB 2025.1.Hs gene sets. To reproduce exact enrichment results, use the matching package/database versions; the script writes `R_SESSION_INFO.txt` at the end. No dependency lockfile is provided because this script has not been run with all original data in the distribution environment.

## 2. Input directory layout

Create this structure underneath the repository root. All **paths are relative** to `data/` and refer to **input filenames expected by the code**, not to old lab computers or old R-output directories.

```text
data/
  in_house/
    Count_Matrix.txt
    LCMSMSRNA-Seq/03_DataAnalysis/192_0002_XPO1_Report.xlsx
  TCGA_PRAD/
    gdc_sample_sheet.2024-08-29.tsv
    counts/                    # directory with the original gene-count files
  GSE35988/
    GSE35988-GPL6480_series_matrix.txt.gz
  prad_su2c_2019/
    ARPCvNEPCvDNPC_revision.RData   # analysis-ready, author-prepared derivative
    expression_with_clinical.xlsx   # analysis-ready, author-prepared derivative
    data_clinical_sample.txt
  nepc_wcm_2016/
    data_mrna_seq_fpkm.txt
    data_clinical_sample.txt
  public_proteomics/
    10780432ccr181229-sup-200506_2_supp_4875649_pbgdn5.xlsx
  clinical.project-tcga-prad.2024-09-06/
    clinical.tsv
  prad_mskcc/
    data_mrna_agilent_microarray.txt
    data_clinical_sample.txt
  GSE4016_series_matrix.txt
  references/
    Savitski_halflives.csv    # optional: auto-downloaded if absent
```

**In-house input requirements:** `Count_Matrix.txt` is the original gene-by-sample raw RNA count matrix, with sample IDs that match the code's PC3 contrasts; the DIA workbook must have a sheet named **`Proteins`** containing `Gene` and the experiment's normalized abundance columns (the script validates a required subset). The repository does **not** contain these two files. Any distribution requires an appropriate public release or permission from their owners.

**Public-cohort derived inputs:** `ARPCvNEPCvDNPC_revision.RData` and `expression_with_clinical.xlsx` are specifically named, author-prepared analysis-ready objects from the prior workflow. They are **not** guaranteed to be files distributed by cBioPortal. Users must obtain them from the study's associated source-data release or reproduce their preprocessing before the full script can run. Do not rename unrelated cBioPortal downloads to satisfy these filenames.

The GDC sample sheet is the manifest originally used to match files/samples; its historical date is part of the expected filename. A fresh GDC download may use different names/workflows or updated expression quantification. Verify the assay, sample mapping, and gene identifiers before substituting files.

## 3. Public sources and download locations

| Public resource used | Original study / content | Official study/download location | Expected local file(s) / important note |
|---|---|---|---|
| **TCGA-PRAD** (GDC) | Primary prostate tumors, tumor–normal matched RNA, clinical data | [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/projects/TCGA-PRAD) and [GDC Repository](https://portal.gdc.cancer.gov/) | `TCGA_PRAD/counts/`, `TCGA_PRAD/gdc_sample_sheet.2024-08-29.tsv`, `clinical.project-tcga-prad.2024-09-06/clinical.tsv`. Reproduce original sample sheet and count workflow where possible. |
| **TCGA Prostate Adenocarcinoma** (cBioPortal) | Expression/clinical API requests for patient-level correlations | [TCGA Firehose Legacy/TCGA study browsing](https://www.cbioportal.org/) and [cBioPortal API](https://www.cbioportal.org/api/swagger-ui/index.html) | Queried programmatically in portions of the script; study identifiers/profile availability can change. |
| **GSE35988**, **GPL6480** | Tomlins et al., prostate tissue microarray data | [GEO GSE35988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35988) | `GSE35988/GSE35988-GPL6480_series_matrix.txt.gz` (select **GPL6480** series matrix). |
| **GSE4016**, **GPL3341** | Gene expression in androgen-sensitive and androgen-insensitive prostate cancer cell lines | [GEO GSE4016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4016) | `GSE4016_series_matrix.txt` (decompress GEO's series matrix if necessary); GPL annotation is retrieved with `GEOquery`. |
| **GSE21032** (MSKCC / Taylor) | Expression in prostate tumors, used in FigS1G alongside GSE4016 | [GEO GSE21032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032) | The R script uses the `prad_mskcc` expression table for this panel; GSE21032 is a related original public accession, not an extra local input filename. |
| **SU2C/PCF 2019** (`prad_su2c_2019`) | Abida et al., advanced prostate cancer expression/clinical outcomes | [cBioPortal study](https://www.cbioportal.org/study/summary?id=prad_su2c_2019) · [cBioPortal Datahub](https://github.com/cBioPortal/datahub/tree/master/public/prad_su2c_2019) | `data_clinical_sample.txt` is public; `expression_with_clinical.xlsx` and `ARPCvNEPCvDNPC_revision.RData` are locally prepared derivatives and must be reconstructed/obtained separately. |
| **WCM NEPC 2016** (`nepc_wcm_2016`) | Beltran et al., CRPC-adenocarcinoma vs neuroendocrine CRPC | [cBioPortal study](https://www.cbioportal.org/study/summary?id=nepc_wcm_2016) · [Datahub](https://github.com/cBioPortal/datahub/tree/master/public/nepc_wcm_2016) | `data_mrna_seq_fpkm.txt`, `data_clinical_sample.txt`; the analysis aggregates repeat samples at the patient level for specific plots. |
| **MSKCC prostate 2010** (`prad_mskcc`) | Taylor et al., integrative genomic profiling | [cBioPortal study](https://www.cbioportal.org/study/summary?id=prad_mskcc) · [Datahub](https://github.com/cBioPortal/datahub/tree/master/public/prad_mskcc) | `data_mrna_agilent_microarray.txt`, `data_clinical_sample.txt`; check original profile name and patient mapping. |
| **Clinical Cancer Research 2018 prostate proteomics** | *The Proteome of Prostate Cancer Bone Metastasis Reveals Heterogeneity with Prognostic Implications* | [AACR Clinical Cancer Research article and supplementary data](https://aacrjournals.org/clincancerres/article/24/21/5433/281650/The-Proteome-of-Prostate-Cancer-Bone-Metastasis) | Download Table S1 workbook with exact name `10780432ccr181229-sup-200506_2_supp_4875649_pbgdn5.xlsx`. |
| **Savitski half-life reference** | Processed protein half-life measurements, Leduc/Slavov deposition | [Zenodo 14827610](https://zenodo.org/records/14827610) · [CSV download](https://zenodo.org/records/14827610/files/Savitski_halflives.csv?download=1) | `references/Savitski_halflives.csv`. If absent, the S7F part attempts to download the same CSV to its current output folder. |
| **MSigDB** | Hallmark, GO:BP, Reactome, KEGG legacy and other gene sets | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) · [2025.1.Hs release files](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/) | Read by `msigdbr` in the R analysis. Be aware that dynamically loaded gene sets can change with package/database versions. |

For cBioPortal, use its **Datasets** section or Datahub to obtain the correctly named studies. [Official cBioPortal download guidance](https://docs.cbioportal.org/downloads/) documents both procedures. The exact filenames used above are based on the analysis code, not a promise that all such files are directly downloadable under those names.

### In-house dataset availability

The study's original PC3 nuclear/cytoplasmic RNA-seq and DIA analyses require the in-house count matrix and DIA workbook. Accession IDs associated with the project include **GEO GSE336778** (RNA-seq) and **PRIDE PXD080083** (proteomics), but the analysis driver still expects the authors' specific *processed* matrix/workbook schema. Repository deposition and public availability should be checked before publishing download instructions; raw repositories do **not** automatically recreate the exact intermediate `.txt` and `.xlsx` files.

## 4. Run

From the repository's top directory in RStudio:

```r
source("NMD3_Figures_1_to_S8.R")
```

Or in a shell (with the working directory set to the repository):

```bash
Rscript NMD3_Figures_1_to_S8.R
```

On startup the script creates `results/NMD3_ALL_FIGURES_YYYYMMDD_HHMMSS/`, prints `[RUN STARTED]`, and validates inputs and packages. Each computational Figure PDF is saved to `FIGURES_NEW_ORDER/` **immediately when generated**; CSV statistical outputs are written within the run directory. A failure will stop the script and leave successfully completed PDFs intact.

Optional environment variables (the *script itself* contains no personal absolute paths):

```r
Sys.setenv(NMD3_DATA_ROOT = file.path(getwd(), "data"))
Sys.setenv(NMD3_OUTPUT_ROOT = file.path(getwd(), "results"))
source("NMD3_Figures_1_to_S8.R")
```

`NMD3_PROJECT_ROOT` can instead set a different project root. `NMD3_HALF_LIFE_FILE` can point to a separately downloaded half-life CSV. Without overrides, all inputs are read from this repository's `data/`, and outputs go into this repository's `results/`. The script does **not** use an old `.RData`, `.RDS`, old figure directory, or hard-coded user account location.

## 5. Outputs and figure-specific notes

- **Main Figures:** the computational panels of Fig1–Fig6, in order. Figure 5 and Fig6J/6K have no R-source plot.
- **Supplementary Figures:** computational panels of FigS1–FigS8, in order. **FigS4H is excluded** because it is hand-drawn; experimental panels are not reconstructed.
- **Fig2F and FigS7D** use **all eligible RPL/RPS genes** satisfying the original panel-level data-quality filters, rather than an author-selected short gene list.
- **FigS2A** contains **three complete GSEA plots side by side**: KEGG Prostate Cancer, Reactome Translation, and Ribosomal Large Subunit Export from Nucleus. Results are calculated from the existing TCGA NMD3-high vs NMD3-low analysis; no NES/FDR values are hard-coded into the plot.
- **FigS6B** (siNMD3 nuclear/cytoplasmic RNA volcano) uses `padj < 0.05` and `abs(log2FC) >= 1` to classify colored points. Its label candidates are ranked among these qualifying significant points. **No special inclusion or exclusion is applied to PEG10**: it is treated like every other gene, and receives a label only if it qualifies under the same ranking and significance rules.
- **Fig6I** uses independently measured finite KPNB1/CCT5 protein values, with paired samples only for Normal–Primary paired comparisons, and no arbitrary `n >= 3` group requirement.
- **No extra PDF figures:** non-indexed PDF requests are rejected by the output helper. Additional *intermediate* CSVs and sensitivity statistics may still be written to subfolders; they are not extra PPT panels.

Expected tracking files (subject to completion) include `PDF_SAVE_LOG.csv`, `ANALYSIS_EXECUTION.csv`, `FIGURE_PANEL_INDEX.csv`, `RUN_STATUS.txt`, `R_SESSION_INFO.txt`, `ALL_OUTPUT_FILES.csv`, and `NMD3_COMPLETE_ANALYSIS_WORKSPACE.RData`. The `.RData` is a generated **output**, not a required input. The session/workspace files can be very large and should not be added to Git.

## 6. Reproducibility and limitations

- Analysis thresholds, tests, patient-level grouping and modeling were carried forward from the supplied final script. The only figure-specific changes in this release are **FigS2A layout/complete plot rendering** and **FigS6B label behavior**, together with **portable file-path resolution**.
- The analysis depends on specific original upstream data formats and cohort versions. An independently downloaded or updated public cohort can legitimately differ in gene coverage and per-comparison sample counts. The script does not fabricate missing observations or replace unavailable data.
- This is a **GitHub-ready source distribution**, **not a claim of end-to-end re-execution**. Without the authors' original in-house/derived files and compatible R packages, a full live run cannot be verified here. Please report an exact error message and input schema if a module fails; debugging should target the affected module instead of changing biological definitions or rerunning completed analyses unnecessarily.
- Check permissions, patient-data terms, and repository release status before uploading any source data. The `.gitignore` intentionally excludes `data/` and `results/`.
