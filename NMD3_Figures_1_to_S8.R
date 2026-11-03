# NMD3 and functional ribosome supply in advanced prostate cancer
# Figure-reproduction analysis: Main Figures 1–6, followed by Supplementary S1–S8.
# Reorganized from the 2026-10-30 author analysis. Statistical procedures retained.
# Place required inputs in data/ (see README.md). Output goes to results/.
# No user-specific paths or dependency on previous results/workspaces.
# Note: requires the in-house RNA-seq and DIA input datasets.

nmd3_run_main <- function() {
.nmd3_previous_dir <- getwd()
on.exit(setwd(.nmd3_previous_dir), add=TRUE)
# Project directory defaults to the location of this script when sourced.
.nmd3_script_paths <- vapply(sys.frames(), function(frame) {
  p <- frame$ofile
  if (is.null(p) || !is.character(p) || !length(p)) "" else p[1L]
}, character(1))
.nmd3_script_paths <- .nmd3_script_paths[nzchar(.nmd3_script_paths)]
.nmd3_script_dir <- if (length(.nmd3_script_paths)) {
  dirname(normalizePath(tail(.nmd3_script_paths, 1L), mustWork=TRUE))
} else normalizePath(getwd(), mustWork=TRUE)
.nmd3_project_root <- path.expand(Sys.getenv(
  "NMD3_PROJECT_ROOT", unset=.nmd3_script_dir))
.nmd3_data_root <- path.expand(Sys.getenv(
  "NMD3_DATA_ROOT", unset=file.path(.nmd3_project_root, "data")))
.nmd3_output_base <- path.expand(Sys.getenv(
  "NMD3_OUTPUT_ROOT", unset=file.path(.nmd3_project_root, "results")))
.nmd3_run_root <- file.path(.nmd3_output_base,
  paste0("NMD3_ALL_FIGURES_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (dir.exists(.nmd3_run_root)) .nmd3_run_root <- paste0(.nmd3_run_root,"_",Sys.getpid())
# Create the actual output folder BEFORE package checks and file discovery.
# This makes the run visible immediately in Finder, including if setup fails.
dir.create(.nmd3_run_root, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(.nmd3_run_root)) stop("Cannot create output directory: ", .nmd3_run_root)
dir.create(file.path(.nmd3_run_root, "FIGURES_NEW_ORDER"), recursive=TRUE, showWarnings=FALSE)
writeLines(c(
  paste0("STARTED: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("OUTPUT_ROOT: ", .nmd3_run_root),
  "Status: setup started; this file does NOT imply PDFs have been created."
), file.path(.nmd3_run_root, "STARTED_SETUP.txt"))
cat("\n[RUN STARTED] ", .nmd3_run_root, "\n", sep="")
cat("[PDF FOLDER] ", file.path(.nmd3_run_root, "FIGURES_NEW_ORDER"), "\n", sep="")
flush.console()
.nmd3_required_packages <- c("DESeq2","apeglm","ashr","readxl","limma","statmod",
 "dplyr","tibble","tidyr","ggplot2","ggrepel","ggpubr","patchwork","pheatmap",
 "gridExtra","GEOquery","Biobase","GSVA","msigdbr","clusterProfiler",
 "org.Hs.eg.db","AnnotationDbi","enrichplot","janitor","readr","data.table",
 "stringr","survival","survminer","VennDiagram","purrr","cBioPortalData",
 "Seurat","fgsea","forcats","ggforce","nlme","scales")
.nmd3_missing_packages <- .nmd3_required_packages[!vapply(.nmd3_required_packages,
 requireNamespace,logical(1),quietly=TRUE)]
if(length(.nmd3_missing_packages)) stop("Missing R packages: ",
 paste(.nmd3_missing_packages,collapse=", "),". Install these before running.")
options(timeout=max(600,getOption("timeout")))
# missing cohort data, refit from saved graphics, or download unpublished data.
#
# Inputs are resolved exclusively from <project>/data/ (or NMD3_DATA_ROOT).
# The working output folder contains symlinks into that source directory, not copies.
# This deliberately avoids scanning Desktop/Downloads or a previous R workspace.
if (!dir.exists(.nmd3_data_root)) {
  stop("Data folder not found: ", .nmd3_data_root,
       ". Create the documented data/ tree or set NMD3_DATA_ROOT.")
}
.nmd3_data_root <- normalizePath(.nmd3_data_root, mustWork=TRUE)
.nmd3_roots <- .nmd3_data_root

# Relative locations expected by the existing manuscript code. These strings
# are the ACTUAL input expectations of that code, not download URLs.
.nmd3_expected <- c(
  "data/TCGA_PRAD/gdc_sample_sheet.2024-08-29.tsv",
  "data/TCGA_PRAD/counts",
  "data/GSE35988/GSE35988-GPL6480_series_matrix.txt.gz",
  "data/prad_su2c_2019/ARPCvNEPCvDNPC_revision.RData",
  "data/nepc_wcm_2016/data_mrna_seq_fpkm.txt",
  "data/nepc_wcm_2016/data_clinical_sample.txt",
  "data/public_proteomics/10780432ccr181229-sup-200506_2_supp_4875649_pbgdn5.xlsx",
  "data/clinical.project-tcga-prad.2024-09-06/clinical.tsv",
  "data/in_house/Count_Matrix.txt",
  "data/in_house/LCMSMSRNA-Seq/03_DataAnalysis/192_0002_XPO1_Report.xlsx"
)

# All file discovery happens ONCE before the analyses. Symlink located inputs
# into the output area's own data/ tree, keeping originals unchanged.
.nmd3_is_present <- function(p) {
  if (!file.exists(p)) return(FALSE)
  if (basename(p) == "counts" && dir.exists(p)) {
    return(length(list.files(p, recursive=TRUE, all.files=FALSE)) > 0L)
  }
  TRUE
}
.nmd3_relative <- sub("^data/", "", .nmd3_expected)
.nmd3_found <- rep(NA_character_,length(.nmd3_expected))
names(.nmd3_found) <- .nmd3_expected

# Resolve exact input filenames; no automated search of private machine folders.
for (.i in seq_along(.nmd3_expected)) {
  .candidate <- file.path(.nmd3_data_root, .nmd3_relative[.i])
  if (.nmd3_is_present(.candidate)) .nmd3_found[.i] <- .candidate
}

# -----------------------------------------------------------------------------
# Validate the original in-house DIA workbook in data/in_house/.
# Do not search private folders or substitute a different proteomics experiment.
# -----------------------------------------------------------------------------
.nmd3_dia_key <- paste0(
  "data/in_house/LCMSMSRNA-Seq/03_DataAnalysis/",
  "192_0002_XPO1_Report.xlsx"
)
.nmd3_dia_required_columns <- c(
  "Gene", "Norm_DMSO_Cyto_D1C", "Norm_DMSO_Cyto_D2C", "Norm_DMSO_Cyto_D3C",
  "Norm_XPO1inhibitor_Cyto_X1C", "Norm_XPO1inhibitor_Cyto_X2C",
  "Norm_XPO1inhibitor_Cyto_X3C"
)
.nmd3_valid_dia <- function(p) {
  if (length(p) != 1L || is.na(p) || !nzchar(p) || !file.exists(p) || dir.exists(p))
    return(FALSE)
  .sheets <- tryCatch(readxl::excel_sheets(p), error=function(e) character(0))
  if (!("Proteins" %in% .sheets)) return(FALSE)
  .cols <- tryCatch(names(readxl::read_excel(p, sheet="Proteins", n_max=0,
                           .name_repair="minimal")), error=function(e) character(0))
  all(.nmd3_dia_required_columns %in% .cols)
}
.nmd3_dia_real_path <- file.path(.nmd3_data_root,
  "in_house/LCMSMSRNA-Seq/03_DataAnalysis/192_0002_XPO1_Report.xlsx")
if (!.nmd3_valid_dia(.nmd3_dia_real_path)) {
  stop("Required in-house DIA workbook missing or has unexpected columns: ",
       .nmd3_dia_real_path, ". See README.md for input schema.")
}
.nmd3_dia_real_path <- normalizePath(.nmd3_dia_real_path, mustWork=TRUE)
.nmd3_found[.nmd3_dia_key] <- .nmd3_dia_real_path
cat("\nVerified in-house DIA workbook.\n")

# Fail-safe: all results saved in an isolated folder so a previous project's
# data/ symlink cannot accidentally cause writes into your original data.
dir.create(.nmd3_run_root,recursive=TRUE,showWarnings=FALSE)
setwd(.nmd3_run_root)
dir.create("data",recursive=TRUE,showWarnings=FALSE)
dir.create("results",recursive=TRUE,showWarnings=FALSE)
dir.create("extra_sensitivity",recursive=TRUE,showWarnings=FALSE)
.nmd3_link <- function(target, original) {
  # A previous failed run can leave a broken or outdated symlink.
  # Remove symlinks only; never delete a real user data file.
  .existing_link <- suppressWarnings(Sys.readlink(target))
  if (length(.existing_link) == 1L && !is.na(.existing_link) &&
      nzchar(.existing_link)) {
    .old_real <- if (file.exists(target)) normalizePath(target, mustWork=TRUE) else ""
    .new_real <- if (file.exists(original)) normalizePath(original, mustWork=TRUE) else ""
    if (!identical(.old_real, .new_real)) {
      if (unlink(target) != 0L) stop("Could not clear stale data link: ", target)
    }
  }
  if (.nmd3_is_present(target)) return(invisible(TRUE))
  if (!.nmd3_is_present(original)) return(invisible(FALSE))
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  status <- suppressWarnings(file.symlink(original,target))
  if(!isTRUE(status)) stop("Cannot link ",original," to ",target,
    ". Check permissions; no input files have been copied.")
  invisible(TRUE)
}
for(.i in which(!is.na(.nmd3_found)))
  .nmd3_link(.nmd3_expected[.i],.nmd3_found[.i])

if (!.nmd3_valid_dia(.nmd3_dia_key)) {
  stop("The DIA symlink does not resolve to a valid Proteins workbook at: ",
       file.path(.nmd3_run_root, .nmd3_dia_key))
}
.nmd3_dia_link_real <- normalizePath(.nmd3_dia_key, mustWork=TRUE)
if (!identical(.nmd3_dia_link_real, .nmd3_dia_real_path)) {
  stop("DIA source/link mismatch: expected ", .nmd3_dia_real_path,
       " but linked ", .nmd3_dia_link_real)
}

.nmd3_present <- vapply(.nmd3_expected,.nmd3_is_present,logical(1))
.nmd3_inventory <- data.frame(
  Required_input=.nmd3_expected,
  Found=.nmd3_present,
  Original_file=ifelse(is.na(.nmd3_found),"",.nmd3_found),
  stringsAsFactors=FALSE
)
utils::write.csv(.nmd3_inventory,"INPUT_INVENTORY.csv",row.names=FALSE)
cat("\nNMD3 input discovery: ",sum(.nmd3_present),"/",
    length(.nmd3_present)," input paths found.\n",sep="")
.nmd3_missing <- .nmd3_expected[!.nmd3_present]
if(length(.nmd3_missing)) {
  cat("\nStill missing actual source data (see INPUT_INVENTORY.csv):\n")
  print(.nmd3_missing)
  stop("Analysis not started: ",length(.nmd3_missing),
    " input paths still missing. Public source data have to be present ",
    "under NMD3_DATA_ROOT; populate the documented data/ tree. Inventory: ",
    file.path(.nmd3_run_root,"INPUT_INVENTORY.csv"))
}
rm(list=intersect(c(".nmd3_candidate_files", ".nmd3_present", ".nmd3_found"), ls(all.names=TRUE)))

# Shared MSigDB calls are cached so analyses do not separately download/re-read
# the same MSigDB gene set multiple times.
.nmd3_msig_cache <- new.env(parent=emptyenv())
nmd3_cached_msig <- function(...) {
  .args <- list(...)
  .key <- paste(capture.output(dput(.args)),collapse="")
  if(!exists(.key,envir=.nmd3_msig_cache,inherits=FALSE)) {
    assign(.key,do.call(msigdbr::msigdbr,.args),envir=.nmd3_msig_cache)
  }
  get(.key,envir=.nmd3_msig_cache,inherits=FALSE)
}

# Work around `AnnotationDbi::select` / `dplyr::select` masking.
suppressPackageStartupMessages(library(dplyr))
select <- dplyr::select
set.seed(123)
cat("\nNMD3 analysis starting (Main Figures 1–6, then Supplementary S1–S8). Output: ",.nmd3_run_root,"\n",sep="")


# Raw data are read once. Module environments prevent cross-cohort object overwrites.
prot_raw <- readxl::read_excel("data/in_house/LCMSMSRNA-Seq/03_DataAnalysis/192_0002_XPO1_Report.xlsx",sheet="Proteins")
expr_matrix <- read.table("data/in_house/Count_Matrix.txt",header=TRUE,sep="\t",row.names=1,check.names=FALSE)
output_dir <- .nmd3_run_root
select <- dplyr::select
msigdbr <- nmd3_cached_msig

.nmd3_panels <- data.frame(do.call(rbind,list(
c("Fig1A","","Pre-60S/NMD3 export experimental schematic","non-R schematic"),
c("Fig1B","Fig1B_siNMD3_Nuclear_Retention.pdf","siNMD3 protein nuclear-retention volcano (left); right EIF6 WB experimental","partial"),
c("Fig1C","","siNMD3 polysome profiles","experimental"),
c("Fig1D","","siNMD3 puromycin immunoblot and densitometry","experimental"),
c("Fig1E","","NMD3 rescue WB and cell viability","experimental"),
c("Fig1F","","LNCaP NMD3 overexpression, puromycin, proliferation","experimental"),
c("Fig2A","Fig2A_KPT330_Nuclear_vs_Cytoplasmic_RNA.pdf","KPT-330 nuclear versus cytoplasmic RNA log2FC","computational"),
c("Fig2B","Fig2B_KPT330_Protein_Nuclear_Retention.pdf","KPT-330 nuclear/cytoplasmic protein redistribution","computational"),
c("Fig2C","","KPT-330 polysome profiles","experimental"),
c("Fig2D","","KPT-330 puromycin pulse labeling","experimental"),
c("Fig2E","Fig2E_KPT330_Nuclear_Cytoplasmic_Protein_Volcano.pdf","KPT-330 nuclear and cytoplasmic DIA volcano","computational"),
c("Fig2F","Fig2F_RPL_RPS_RNA_Protein_Heatmaps.pdf","ALL eligible RPL/RPS, KPT-330 RNA and protein heatmaps","computational"),
c("Fig3A-left","Fig3A_left_KPT330_RNA_Protein_Quadrants.pdf","KPT-330 paired RNA/protein quadrants","computational"),
c("Fig3A-right","Fig3A_right_siNMD3_RNA_Protein_Quadrants.pdf","siNMD3 paired RNA/protein quadrants","computational"),
c("Fig3B","Fig3B_Perturbed_Protein_RNA_Summary.pdf","Perturbed protein/RNA counts from the FigS6D gene sets","computational"),
c("Fig3C-upper","Fig3C_upper_Central95_RNA_Protein.pdf","Central 95 percent quadrant regions","computational"),
c("Fig3C-lower","Fig3C_lower_RNA_Protein_Spearman.pdf","Spearman protein/RNA correlations","computational"),
c("Fig3D","Fig3D_CAT_RNA_Protein.pdf","Concordance-at-the-top (CAT) RNA/protein analysis","computational"),
c("Fig3E","Fig3E_4Groups.pdf","KPT-330 four-group nuclear-retention sensitivity","computational"),
c("Fig3F","Fig3F_siNMD3_4Groups.pdf","siNMD3 four-group nuclear-retention sensitivity","computational"),
c("Fig4A","Fig4A_RNA_Stable_Hallmark.pdf","RNA-stable Hallmark pathways siNMD3","computational"),
c("Fig4B","Fig4B_Protein_RNA_Residual_Pathways.pdf","Hallmark pathway protein-RNA residuals","computational"),
c("Fig4C","Fig4C_Proteostasis_GSEA.pdf","Proteostasis/turnover protein enrichment","computational"),
c("Fig4D","Fig4D_MYC_Targets_RNA_log2FC.pdf","MYC targets versus other RNA log2FC","computational"),
c("Fig4E","Fig4E_MYC_Targets_Protein_RNA_Residual.pdf","MYC targets versus other protein-RNA residual","computational"),
c("Fig4F","","KPNB1/CCT5 mRNA and protein validation","experimental"),
c("Fig5A","","siNMD3 PEG10 immunoblot in PC3 and DU145","experimental"),
c("Fig5B","","PEG10 after siXPO1/KPT-330/Portimine A","experimental"),
c("Fig5C","","PEG10 and ENO2 qRT-PCR","experimental"),
c("Fig5D","","NMD3 WT versus 3Lmut PEG10 rescue WB","experimental"),
c("Fig5E","","NMD3 NES mutant rescue of cell proliferation","experimental"),
c("Fig5F","","siPEG10 proliferation and knockdown Western blot","experimental"),
c("Fig6A","Fig6A_Proteomics_Ribosome_Export_Score_7_Pairs.pdf","7 matched normal/primary protein ribosome-export score","computational"),
c("Fig6B","Fig6B_Proteomics_Ribosome_Export_Score_Three_Groups.pdf","Normal/primary/metastatic protein export-score violin","computational"),
c("Fig6C","Fig6C_Ribosome_Export_Survival.pdf","TCGA ribosome-export score survival","computational"),
c("Fig6D","Fig6D_WCM_SU2C_PEG10_ENO2_Heatmaps.pdf","PEG10/ENO2 WCM and SU2C expression heatmaps","computational"),
c("Fig6E","Fig6E_GSE35988_PEG10_Primary_vs_Metastatic.pdf","GSE35988 PEG10 primary/metastatic expression","computational"),
c("Fig6F","Fig6F_GSE35988_ENO2_Primary_vs_Metastatic.pdf","GSE35988 ENO2 primary/metastatic expression","computational"),
c("Fig6G","Fig6G_level_NE_score.pdf","WCM NEPC neuroendocrine score","computational"),
c("Fig6H","Fig6H_level_PEG10.pdf","WCM NEPC PEG10 patient-level expression","computational"),
c("Fig6I","Fig6I_KPNB1_CCT5_Three_Group_Protein_Combined.pdf","Public proteomics KPNB1 and CCT5","computational"),
c("Fig6J","","Clinical IHC images (manuscript Fig6J)","experimental"),
c("Fig6K","","Mechanistic schematic (manuscript Fig6K)","non-R schematic")

)),stringsAsFactors=FALSE)
names(.nmd3_panels) <- c("Panel","PDF","Description","Coverage")

.nmd3_order_folder <- file.path(.nmd3_run_root,"FIGURES_NEW_ORDER")
dir.create(.nmd3_order_folder,recursive=TRUE,showWarnings=FALSE)
# Only indexed PPT panels are generated; no extra-PDF directory is created.
.nmd3_run_log <- file.path(.nmd3_run_root,"RUN_STATUS.txt")
.nmd3_pdf_log <- file.path(.nmd3_run_root,"PDF_SAVE_LOG.csv")
.nmd3_panels$Status <- ifelse(nzchar(.nmd3_panels$PDF),"PENDING","NOT_R_INPUT")
.nmd3_panels$Path <- ""
.nmd3_state <- new.env(parent=emptyenv())
.nmd3_state$modules <- list()
.nmd3_state$devices <- list()
.nmd3_state$last_panel <- 0L
.nmd3_state$current <- "setup"
.nmd3_state$saved <- character()
.nmd3_state$steps <- data.frame(Module=character(),Status=character(),Elapsed_seconds=numeric())
.nmd3_context <- environment()
nmd3_write_status <- function() {
 utils::write.csv(.nmd3_panels,file.path(.nmd3_run_root,"FIGURE_PANEL_INDEX.csv"),row.names=FALSE)
 utils::write.csv(.nmd3_state$steps,file.path(.nmd3_run_root,"ANALYSIS_EXECUTION.csv"),row.names=FALSE)
}
nmd3_log <- function(...) {
 message(paste0(...))
 cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S")," ",paste0(...),"\n",file=.nmd3_run_log,append=TRUE,sep="")
 flush.console()
}
nmd3_target <- function(filename) {
 stopifnot(is.character(filename),length(filename)==1L,nzchar(filename))
 name <- basename(filename)
 if (!grepl("\\.pdf$",name,ignore.case=TRUE)) stop("Expected PDF filename: ",filename)
 i <- match(name,.nmd3_panels$PDF)
 if(!is.na(i)) {
  file.path(.nmd3_order_folder,name)
 } else {
  stop("Non-PPT/extra PDF requested: ", name,
       ". This final run only writes indexed Figure 1–6 and S1–S8 panels.")
 }
}
nmd3_check_pdf <- function(path) {
 if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size<200)
  stop("PDF missing or empty: ",path)
 con<-file(path,"rb");on.exit(close(con))
 if (!identical(rawToChar(readBin(con,"raw",n=5)),"%PDF-")) stop("Invalid PDF header: ",path)
}
nmd3_record_pdf <- function(path) {
 nmd3_check_pdf(path)
 i<-match(basename(path),.nmd3_panels$PDF)
 if(!is.na(i)) {
  if(i<.nmd3_state$last_panel)stop("PDF order regression at ",.nmd3_panels$Panel[i])
  .nmd3_state$last_panel<-i
  .nmd3_panels$Status[i]<<-"SAVED"
  .nmd3_panels$Path[i]<<-path
 }
 row<-data.frame(Time=format(Sys.time(),"%Y-%m-%d %H:%M:%S"),Module=.nmd3_state$current,
  Panel=if(is.na(i))"EXTRA" else .nmd3_panels$Panel[i],PDF=path,Bytes=file.info(path)$size)
 utils::write.table(row,.nmd3_pdf_log,append=file.exists(.nmd3_pdf_log),sep=",",
  row.names=FALSE,col.names=!file.exists(.nmd3_pdf_log),quote=TRUE)
 .nmd3_state$saved<-c(.nmd3_state$saved,path)
 nmd3_write_status()
 nmd3_log("[PDF SAVED] ",path)
 invisible(path)
}
nmd3_save_pdf <- function(filename,plot=ggplot2::last_plot(),...) {
 target<-nmd3_target(filename)
 if(target %in% .nmd3_state$saved)
   stop("This PDF has already been saved; refusing silent overwrite: ",target)
 nmd3_log("[PDF] ",basename(target))
 args<-list(...)
 # cairo_pdf is used where available, with a native PDF fallback.
 if(!isTRUE(capabilities("cairo")) && !is.null(args$device) && identical(args$device,grDevices::cairo_pdf))
  args$device<-grDevices::pdf
 if(inherits(plot,"ggplot")) {
  # Legends outside the plotting region for survival plots.
  if(grepl("Survival|survival",basename(target)))plot<-plot+ggplot2::theme(legend.position="top")
 }
 do.call(ggplot2::ggsave,c(list(filename=target,plot=plot),args))
 nmd3_record_pdf(target)
 invisible(target)
}
nmd3_open_pdf <- function(filename,width=7,height=7,...,cairo=FALSE) {
 target<-nmd3_target(filename)
 if(target %in% .nmd3_state$saved)stop("Duplicate direct PDF device: ",target)
 if(cairo && isTRUE(capabilities("cairo")))grDevices::cairo_pdf(filename=target,width=width,height=height,...)
 else grDevices::pdf(file=target,width=width,height=height,...)
 .nmd3_state$devices[[as.character(grDevices::dev.cur())]]<-target
 invisible(target)
}
nmd3_pdf <- function(file,width=7,height=7,...)nmd3_open_pdf(file,width,height,...,cairo=FALSE)
nmd3_cairo_pdf <- function(filename,width=7,height=7,...)nmd3_open_pdf(filename,width,height,...,cairo=TRUE)
nmd3_dev_off <- function(which=grDevices::dev.cur()) {
 key<-as.character(which);target<-.nmd3_state$devices[[key]]
 ans<-grDevices::dev.off(which=which)
 .nmd3_state$devices[[key]]<-NULL
 if(!is.null(target))nmd3_record_pdf(target)
 invisible(ans)
}
nmd3_step <- function(name,dependencies,code) {
 if(name %in% names(.nmd3_state$modules))stop("Duplicate module: ",name)
 absent<-setdiff(dependencies,names(.nmd3_state$modules))
 if(length(absent))stop("Missing upstream modules for ",name,": ",paste(absent,collapse=", "))
 env<-new.env(parent=.nmd3_context)
 for(d in dependencies)list2env(as.list.environment(.nmd3_state$modules[[d]],all.names=TRUE),envir=env)
 .nmd3_state$current<-name
 nmd3_log("[START] ",name)
 start<-proc.time()[[3]]
 tryCatch(eval(code,envir=env),error=function(e){
  # Close only devices opened by this module. Completed PDFs stay intact.
  for(key in names(.nmd3_state$devices))try(grDevices::dev.off(as.integer(key)),silent=TRUE)
  .nmd3_state$devices<-list()
  .nmd3_state$steps<-rbind(.nmd3_state$steps,data.frame(Module=name,Status="FAILED",Elapsed_seconds=proc.time()[[3]]-start))
  nmd3_write_status();nmd3_log("[FAILED] ",name,": ",conditionMessage(e))
  stop("Module ",name," failed: ",conditionMessage(e),"\nAlready saved PDFs: ",.nmd3_order_folder,call.=FALSE)
 })
 .nmd3_state$modules[[name]]<-env
 .nmd3_state$steps<-rbind(.nmd3_state$steps,data.frame(Module=name,Status="DONE",Elapsed_seconds=proc.time()[[3]]-start))
 nmd3_write_status();nmd3_log("[DONE] ",name)
 invisible(NULL)
}
nmd3_finish <- function() {
 missing<-.nmd3_panels$Panel[.nmd3_panels$Status=="PENDING"]
 nmd3_write_status()
 writeLines(capture.output(utils::sessionInfo()),file.path(.nmd3_run_root,"R_SESSION_INFO.txt"))
 writeLines(c(
 "MAIN AND SUPPLEMENTARY: follows the 2026-10-30 figure order (Fig1–6 then FigS1–S8).",
 "Figure 5 and Figure 6J/K are wet-lab or schematic panels and have no raw R plot here.",
 "Fig3B uses the original 16-set gene threshold calculations also used for supplementary FigS6D; no supplementary PDF is generated.",
 "Fig2F includes all detected RPL/RPS genes passing original panel-specific QC.",
 "Fig6I takes matrices, metadata, paired sample design and plot styles directly from the completed Fig6A module.",
 "Fig6I uses finite KPNB1/CCT5 values per protein, without a fixed >=3-sample minimum per group.",
 "All original statistical tests and significance thresholds in the remaining main-figure analyses are retained."
 ),file.path(.nmd3_run_root,"SOURCE_NOTES.txt"))
 if(length(missing))stop("Run ended with missing computational panels: ",paste(missing,collapse=", "),call.=FALSE)
 nmd3_log("[COMPLETE] ",sum(.nmd3_panels$Status=="SAVED")," panel PDFs saved to ",.nmd3_order_folder)
}
nmd3_write_status()
nmd3_log("Output folder: ",.nmd3_order_folder)

# -----------------------------------------------------------------------------
# FINAL RUN BACKUP: save every live object (not just the PDFs) at function exit.
# This captures the module environments, complete data objects and plots so a
# completed run can be reopened without recreating upstream analysis objects.
# It runs on normal completion and on a later module failure. It never alters
# statistics and never overwrites files from an earlier timestamped run.
# -----------------------------------------------------------------------------
.nmd3_save_all_outputs <- function() {
  tryCatch({
    nmd3_write_status()
    .workspace_file <- file.path(.nmd3_run_root,
                                 "NMD3_COMPLETE_ANALYSIS_WORKSPACE.RData")
    .object_names <- ls(envir=.nmd3_context, all.names=TRUE)
    base::save(list=.object_names, file=.workspace_file,
               envir=.nmd3_context, compress="gzip")
    nmd3_log("[WORKSPACE SAVED] ",.workspace_file)
  }, error=function(e) {
    message("[WORKSPACE BACKUP FAILED] ",conditionMessage(e))
    try(cat("WORKSPACE BACKUP FAILED: ",conditionMessage(e),"\n",
            file=file.path(.nmd3_run_root,"WORKSPACE_BACKUP_ERROR.txt")),silent=TRUE)
  })

  # List only actual output files, not linked input datasets in data/.
  tryCatch({
    .entries <- list.files(.nmd3_run_root,full.names=TRUE,recursive=FALSE,
                           all.files=FALSE)
    .entries <- .entries[basename(.entries)!="data"]
    .files <- unlist(lapply(.entries,function(p) {
      if(dir.exists(p))list.files(p,full.names=TRUE,recursive=TRUE,
                                  include.dirs=FALSE,all.files=FALSE)
      else p
    }),use.names=FALSE)
    .files <- unique(.files[file.exists(.files) & !dir.exists(.files)])
    .info <- file.info(.files)
    .manifest <- data.frame(
      Relative_path=substring(.files,nchar(.nmd3_run_root)+2L),
      Bytes=.info$size,
      Modified=as.character(.info$mtime),
      stringsAsFactors=FALSE
    )
    utils::write.csv(.manifest,file.path(.nmd3_run_root,"ALL_OUTPUT_FILES.csv"),
                     row.names=FALSE)
    nmd3_log("[OUTPUT MANIFEST] ",nrow(.manifest)," saved files indexed")
  }, error=function(e) {
    message("[OUTPUT MANIFEST FAILED] ",conditionMessage(e))
  })
}
on.exit(.nmd3_save_all_outputs(),add=TRUE)

# Save the exact source R script when executed with source() in RStudio.
if (length(.nmd3_script_paths)>0L && file.exists(tail(.nmd3_script_paths,1L))) {
  file.copy(tail(.nmd3_script_paths,1L),
            file.path(.nmd3_run_root,"EXECUTED_FINAL_SCRIPT.R"),
            overwrite=TRUE)
}

# =============================================================================
# Fig1B
# =============================================================================
nmd3_step("Fig1B", character(), quote({
# Fig1B — analysis below
# Nuclear retention after NMD3 knockdown
#
# One limma interaction analysis and one display version:
#   - RPL proteins: one yellow color
#   - RPS proteins: one green color
#   - Colors do not depend on the P value
#   - No P = 0.05 horizontal line
#
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

stopifnot(
  exists("prot_raw"),
  "Gene" %in% colnames(prot_raw)
)

# ------------------------------------------------------------
# 1. Gene definitions
# ------------------------------------------------------------

is_RPS_gene_5C <- function(x) {
  grepl(
    "^RPS[0-9]+[A-Z]?$",
    toupper(x),
    perl = TRUE
  )
}

split_symbols_5C <- function(x) {
  toupper(
    trimws(
      unlist(
        strsplit(
          as.character(x),
          ";"
        )
      )
    )
  )
}

rm_genes_5C <- c(
  "RPS6KA1",
  "RPS6KA3",
  "RPS6KA4",
  "RPS6KA5",
  "RPS6KB1",
  "RPS6KB2",
  "RPS6KC1",
  "RPS19BP1"
)

# ------------------------------------------------------------
# 2. Select the 12 siCTRL/siNMD3 samples
# ------------------------------------------------------------

sample_cols_5C <- grep(
  paste0(
    "^Norm_siRNAx(Ctrl|NMD3xKD)_",
    "(Cyto|Nuc)_[A-Z][123][CN]$"
  ),
  colnames(prot_raw),
  value = TRUE
)

stopifnot(
  length(sample_cols_5C) == 12L,
  !anyDuplicated(sample_cols_5C)
)

meta_5C <- data.frame(
  Sample = sample_cols_5C,
  Treatment = ifelse(
    grepl(
      "siRNAxNMD3xKD",
      sample_cols_5C,
      ignore.case = TRUE
    ),
    "siNMD3",
    "siCTRL"
  ),
  Compartment = ifelse(
    grepl(
      "_Nuc_",
      sample_cols_5C,
      ignore.case = TRUE
    ),
    "Nucleus",
    "Cytoplasm"
  ),
  Replicate = sub(
    ".*_[A-Z]([123])[CN]$",
    "R\\1",
    sample_cols_5C
  ),
  stringsAsFactors = FALSE,
  row.names = sample_cols_5C
)

meta_5C$Treatment <- factor(
  meta_5C$Treatment,
  levels = c("siCTRL", "siNMD3")
)

meta_5C$Compartment <- factor(
  meta_5C$Compartment,
  levels = c("Cytoplasm", "Nucleus")
)

meta_5C$Replicate <- factor(
  meta_5C$Replicate,
  levels = c("R1", "R2", "R3")
)

sample_table_5C <- table(
  meta_5C$Treatment,
  meta_5C$Compartment,
  meta_5C$Replicate
)

cat("\nFigure 5C sample annotation:\n")
print(meta_5C)

cat("\nTreatment x compartment x replicate:\n")
print(sample_table_5C)

stopifnot(
  nrow(meta_5C) == 12L,
  all(sample_table_5C == 1L)
)

# ------------------------------------------------------------
# 3. Protein matrix
# ------------------------------------------------------------

expr_raw_5C <- prot_raw %>%
  dplyr::filter(
    !is.na(Gene),
    trimws(as.character(Gene)) != ""
  ) %>%
  dplyr::select(
    Gene,
    dplyr::all_of(sample_cols_5C)
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols_5C),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

bad_gene_5C <- vapply(
  expr_raw_5C$Gene,
  function(x) {
    any(
      split_symbols_5C(x) %in%
        rm_genes_5C
    )
  },
  logical(1)
)

expr_raw_5C <- expr_raw_5C[
  !bad_gene_5C,
  ,
  drop = FALSE
]

# Aggregate duplicated gene entries once.

expr_5C_df <- expr_raw_5C %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_cols_5C),
      ~ if (all(is.na(.x))) {
        NA_real_
      } else {
        mean(.x, na.rm = TRUE)
      }
    ),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(expr_5C_df) <- as.character(
  expr_5C_df$Gene
)

expr_5C <- as.matrix(
  expr_5C_df[
    ,
    sample_cols_5C,
    drop = FALSE
  ]
)

storage.mode(expr_5C) <- "double"

# ------------------------------------------------------------
# 4. Require >=2 observations in all four conditions
# ------------------------------------------------------------

enough_obs_5C <- function(
    compartment_name,
    treatment_name
) {
  idx <- which(
    meta_5C$Compartment == compartment_name &
      meta_5C$Treatment == treatment_name
  )
  
  rowSums(
    !is.na(
      expr_5C[
        ,
        idx,
        drop = FALSE
      ]
    )
  ) >= 2L
}

keep_5C <- (
  enough_obs_5C("Cytoplasm", "siCTRL") &
    enough_obs_5C("Cytoplasm", "siNMD3") &
    enough_obs_5C("Nucleus", "siCTRL") &
    enough_obs_5C("Nucleus", "siNMD3")
)

keep_5C[is.na(keep_5C)] <- FALSE

expr_5C <- expr_5C[
  keep_5C,
  ,
  drop = FALSE
]

stopifnot(
  nrow(expr_5C) > 0L
)

# ------------------------------------------------------------
# 5. Limma interaction model
#
# effect =
# (Nucleus - Cytoplasm)siNMD3 -
# (Nucleus - Cytoplasm)siCTRL
#
# Positive effect = increased nuclear retention after siNMD3.
# ------------------------------------------------------------

design_5C <- model.matrix(
  ~ Replicate +
    Compartment * Treatment,
  data = meta_5C
)

colnames(design_5C) <- make.names(
  colnames(design_5C)
)

stopifnot(
  qr(design_5C)$rank ==
    ncol(design_5C)
)

coef_interaction_5C <- grep(
  paste0(
    "CompartmentNucleus.*",
    "TreatmentsiNMD3"
  ),
  colnames(design_5C),
  value = TRUE
)

stopifnot(
  length(coef_interaction_5C) == 1L
)

fit_5C <- limma::lmFit(
  expr_5C,
  design_5C
)

fit_5C <- limma::eBayes(
  fit_5C,
  trend = TRUE,
  robust = TRUE
)

tt_5C <- limma::topTable(
  fit_5C,
  coef = coef_interaction_5C,
  number = Inf,
  sort.by = "P"
)

# ------------------------------------------------------------
# 6. Result table and RPL/RPS classification
# ------------------------------------------------------------

res_5C <- data.frame(
  Gene = rownames(tt_5C),
  tt_5C,
  row.names = NULL
) %>%
  dplyr::mutate(
    effect = logFC,
    fdr = adj.P.Val,
    neglogP = -log10(
      pmax(
        P.Value,
        .Machine$double.xmin
      )
    ),
    is_RPL = stringr::str_detect(
      Gene,
      "(?i)^RPL"
    ),
    is_RPS = is_RPS_gene_5C(Gene),
    Family = dplyr::case_when(
      is_RPL ~ "RPL",
      is_RPS ~ "RPS",
      TRUE ~ "Other"
    ),
    
    # Mirror the x-axis so increased nuclear retention is on the left.
    x_plot = -effect
  )

# ------------------------------------------------------------
# 7. RPL/RPS directional statistics
# ------------------------------------------------------------

rpl_5C <- res_5C %>%
  dplyr::filter(is_RPL)

rps_5C <- res_5C %>%
  dplyr::filter(is_RPS)

stopifnot(
  nrow(rpl_5C) > 0L,
  nrow(rps_5C) > 0L
)

binom_rpl_5C <- stats::binom.test(
  x = sum(
    rpl_5C$effect > 0,
    na.rm = TRUE
  ),
  n = nrow(rpl_5C),
  p = 0.5,
  alternative = "greater"
)

binom_rps_5C <- stats::binom.test(
  x = sum(
    rps_5C$effect > 0,
    na.rm = TRUE
  ),
  n = nrow(rps_5C),
  p = 0.5,
  alternative = "two.sided"
)

family_stats_5C <- data.frame(
  Family = c("RPL", "RPS"),
  Total = c(
    nrow(rpl_5C),
    nrow(rps_5C)
  ),
  Nuclear_direction = c(
    sum(rpl_5C$effect > 0),
    sum(rps_5C$effect > 0)
  ),
  Cytoplasmic_direction = c(
    sum(rpl_5C$effect < 0),
    sum(rps_5C$effect < 0)
  ),
  Nominal_nuclear_P_lt_0.05 = c(
    sum(
      rpl_5C$effect > 0 &
        rpl_5C$P.Value < 0.05
    ),
    sum(
      rps_5C$effect > 0 &
        rps_5C$P.Value < 0.05
    )
  ),
  Nominal_cytoplasmic_P_lt_0.05 = c(
    sum(
      rpl_5C$effect < 0 &
        rpl_5C$P.Value < 0.05
    ),
    sum(
      rps_5C$effect < 0 &
        rps_5C$P.Value < 0.05
    )
  ),
  Binomial_P = c(
    binom_rpl_5C$p.value,
    binom_rps_5C$p.value
  )
)

cat("\nFigure 5C RPL/RPS directional statistics:\n")
print(family_stats_5C)

# ------------------------------------------------------------
# 8. Figure 5C
#
# Colors identify protein family only and do not depend on P.
# ------------------------------------------------------------

p_5C <- ggplot(
  res_5C,
  aes(
    x = x_plot,
    y = neglogP
  )
) +
  geom_vline(
    xintercept = 0,
    color = "gray30",
    linewidth = 0.4
  ) +
  geom_point(
    data = subset(
      res_5C,
      Family == "Other"
    ),
    color = "gray75",
    alpha = 0.3,
    size = 1.5
  ) +
  geom_point(
    data = subset(
      res_5C,
      Family == "RPL"
    ),
    aes(fill = "RPL"),
    shape = 21,
    size = 3.4,
    stroke = 0.7,
    color = "black",
    alpha = 0.9
  ) +
  geom_point(
    data = subset(
      res_5C,
      Family == "RPS"
    ),
    aes(fill = "RPS"),
    shape = 21,
    size = 3.4,
    stroke = 0.7,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "RPL" = "#FDB863",
      "RPS" = "#66C2A5"
    ),
    breaks = c("RPL", "RPS"),
    name = NULL,
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 3.4,
        color = "black",
        stroke = 0.7,
        alpha = 1
      )
    )
  ) +
  labs(
    x = expression(
      paste(
        "Mirrored ",
        -Delta,
        Delta,
        " log"[2],
        "FC  (left = nuclear retention)"
      )
    ),
    y = expression(
      -log[10]~"P-value"
    )
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10)
  )

# ------------------------------------------------------------
# 9. Add XPO1 and EIF6 labels
# ------------------------------------------------------------

targets_5C <- c(
  "XPO1",
  "EIF6"
)

target_idx_5C <- which(
  vapply(
    res_5C$Gene,
    function(x) {
      any(
        split_symbols_5C(x) %in%
          targets_5C
      )
    },
    logical(1)
  )
)

if (length(target_idx_5C) > 0L) {
  target_df_5C <- res_5C[
    target_idx_5C,
    c(
      "Gene",
      "x_plot",
      "neglogP"
    ),
    drop = FALSE
  ]
  
  target_df_5C$Label <- vapply(
    target_df_5C$Gene,
    function(x) {
      matched <- split_symbols_5C(x)
      matched <- matched[
        matched %in% targets_5C
      ]
      
      if (length(matched) == 0L) {
        as.character(x)
      } else {
        matched[1]
      }
    },
    character(1)
  )
  
  p_5C <- p_5C +
    geom_point(
      data = target_df_5C,
      aes(
        x = x_plot,
        y = neglogP
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 5,
      stroke = 1.2,
      fill = "white",
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = target_df_5C,
      aes(
        x = x_plot,
        y = neglogP,
        label = Label
      ),
      inherit.aes = FALSE,
      direction = "x",
      hjust = 1,
      nudge_x = -2,
      box.padding = 0.6,
      point.padding = 0.4,
      min.segment.length = 0,
      segment.size = 0.5,
      segment.color = "black",
      max.overlaps = Inf,
      force_pull = 2,
      force = 1
    )
}


nmd3_save_pdf("Fig1B_siNMD3_Nuclear_Retention.pdf", p_5C, width=7, height=5, device=grDevices::cairo_pdf)

}))

# =============================================================================
# Fig2A
# =============================================================================
nmd3_step("Fig2A", character(), quote({
# Fig2A — analysis below
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Read raw counts
# -----------------------------------------------------------------------------

count_file <- "data/in_house/Count_Matrix.txt"

stopifnot(file.exists(count_file))

# Shared count table was loaded once.

sample_columns <- grep(
  "^[DX][0-9]+[CN]$",
  colnames(expr_matrix),
  value = TRUE
)

count_matrix <- as.matrix(
  expr_matrix[, sample_columns, drop = FALSE]
)

stopifnot(
  length(sample_columns) > 0L,
  !anyDuplicated(sample_columns),
  !anyDuplicated(rownames(count_matrix)),
  !anyNA(count_matrix),
  all(is.finite(count_matrix)),
  all(count_matrix >= 0),
  all(count_matrix == round(count_matrix))
)

storage.mode(count_matrix) <- "integer"

# -----------------------------------------------------------------------------
# 2. Sample information
# -----------------------------------------------------------------------------

sample_info <- data.frame(
  Sample = colnames(count_matrix),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    Group = dplyr::recode(
      substr(Sample, 1, 1),
      "D" = "DMSO",
      "X" = "KPT-330"
    ),
    Type = ifelse(
      grepl("N$", Sample),
      "Nucleus",
      "Cytoplasm"
    ),
    Replicate = sub(
      "^[DX]([0-9]+)[CN]$",
      "\\1",
      Sample
    )
  )

design_check <- sample_info %>%
  dplyr::count(Type, Group, name = "n")

print(design_check)

stopifnot(
  nrow(design_check) == 4L,
  all(design_check$n >= 2L),
  setequal(sample_info$Group, c("DMSO", "KPT-330")),
  setequal(sample_info$Type, c("Cytoplasm", "Nucleus")),
  identical(sample_info$Sample, colnames(count_matrix))
)

# -----------------------------------------------------------------------------
# 3. DESeq2 analysis for each compartment
# -----------------------------------------------------------------------------

run_deseq2_by_type <- function(type_name) {
  
  si <- sample_info %>%
    dplyr::filter(Type == type_name)
  
  counts <- count_matrix[
    ,
    si$Sample,
    drop = FALSE
  ]
  
  coldata <- data.frame(
    condition = factor(
      si$Group,
      levels = c("DMSO", "KPT-330")
    ),
    row.names = si$Sample
  )
  
  stopifnot(
    identical(colnames(counts), rownames(coldata)),
    nlevels(coldata$condition) == 2L,
    all(table(coldata$condition) >= 2L)
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )
  
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )
  
  coefficient <- "condition_KPT.330_vs_DMSO"
  
  stopifnot(
    coefficient %in% DESeq2::resultsNames(dds)
  )
  
  res_shr <- DESeq2::lfcShrink(
    dds,
    coef = coefficient,
    type = "apeglm"
  )
  
  result <- as.data.frame(res_shr) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::transmute(
      gene = toupper(as.character(gene)),
      baseMean = as.numeric(baseMean),
      log2FC_RNA = as.numeric(log2FoldChange),
      padj_RNA = as.numeric(padj)
    ) %>%
    dplyr::filter(
      is.finite(log2FC_RNA)
    )
  
  stopifnot(
    nrow(result) > 0L,
    !anyDuplicated(result$gene),
    all(is.finite(result$log2FC_RNA))
  )
  
  result
}

# -----------------------------------------------------------------------------
# 4. Cytoplasmic and nuclear changes
# -----------------------------------------------------------------------------

rna_cyto_df <- run_deseq2_by_type("Cytoplasm") %>%
  dplyr::rename(
    baseMean_cyto = baseMean,
    log2FC_RNA_cyto = log2FC_RNA,
    padj_RNA_cyto = padj_RNA
  )

rna_nuc_df <- run_deseq2_by_type("Nucleus") %>%
  dplyr::rename(
    baseMean_nuc = baseMean,
    log2FC_RNA_nuc = log2FC_RNA,
    padj_RNA_nuc = padj_RNA
  )

rna_cytonuc <- dplyr::inner_join(
  rna_cyto_df,
  rna_nuc_df,
  by = "gene"
) %>%
  dplyr::filter(
    is.finite(log2FC_RNA_cyto),
    is.finite(log2FC_RNA_nuc)
  )

stopifnot(
  nrow(rna_cytonuc) > 0L,
  !anyDuplicated(rna_cytonuc$gene),
  all(is.finite(rna_cytonuc$log2FC_RNA_cyto)),
  all(is.finite(rna_cytonuc$log2FC_RNA_nuc))
)

# -----------------------------------------------------------------------------
# 5. Spearman correlation across all matched genes
# -----------------------------------------------------------------------------

cor_all <- stats::cor.test(
  x = rna_cytonuc$log2FC_RNA_cyto,
  y = rna_cytonuc$log2FC_RNA_nuc,
  method = "spearman",
  alternative = "two.sided",
  exact = FALSE
)

rho <- unname(cor_all$estimate)
p_value <- cor_all$p.value

p_text <- if (p_value < 2.2e-16) {
  "p < 2.2e-16"
} else {
  sprintf("p = %.2e", p_value)
}

correlation_statistics <- data.frame(
  Analysis = "Cytoplasmic vs nuclear RNA log2FC",
  Method = "Two-sided Spearman correlation",
  n_matched_genes = nrow(rna_cytonuc),
  Spearman_rho = rho,
  P_value = p_value,
  stringsAsFactors = FALSE
)

print(cor_all)
print(correlation_statistics)

# -----------------------------------------------------------------------------
# 6. Figure 2F
# -----------------------------------------------------------------------------

fig2f <- ggplot2::ggplot(
  rna_cytonuc,
  ggplot2::aes(
    x = log2FC_RNA_cyto,
    y = log2FC_RNA_nuc
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    color = "black",
    size = 1.4,
    alpha = 0.45
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "red",
    linewidth = 0.9
  ) +
  ggplot2::annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1.3,
    label = paste0(
      "Spearman \u03c1 = ",
      sprintf("%.2f", rho),
      "\n",
      p_text
    ),
    size = 5
  ) +
  ggplot2::labs(
    x = expression(
      paste(
        "Cytoplasmic RNA ",
        log[2],
        "FC (KPT-330 vs DMSO)"
      )
    ),
    y = expression(
      paste(
        "Nuclear RNA ",
        log[2],
        "FC (KPT-330 vs DMSO)"
      )
    )
  ) +
  ggplot2::coord_fixed() +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(
      face = "bold",
      color = "black"
    ),
    axis.text = ggplot2::element_text(
      color = "black"
    ),
    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.2
    )
  )


# -----------------------------------------------------------------------------
# 7. Save outputs
# -----------------------------------------------------------------------------

nmd3_save_pdf(
  filename = "Fig2A_KPT330_Nuclear_vs_Cytoplasmic_RNA.pdf",
  plot = fig2f,
  width = 6,
  height = 6,
  units = "in",
  device = grDevices::cairo_pdf
)

utils::write.csv(
  rna_cytonuc,
  "Fig2F_Cytoplasm_vs_Nuclear_RNA_source_data.csv",
  row.names = FALSE
)

utils::write.csv(
  correlation_statistics,
  "Fig2F_Cytoplasm_vs_Nuclear_RNA_statistics.csv",
  row.names = FALSE
)




}))

# =============================================================================
# Fig2B
# =============================================================================
nmd3_step("Fig2B", character(), quote({
# Fig2B — analysis below
# DeltaDelta = (N-C)DMSO - (N-C)KPT
# Negative = Nuclear Retention; positive = Cytoplasmic Shift
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Sample annotation and protein matrix
# -----------------------------------------------------------------------------

expr_raw <- prot_raw %>%
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::select(Gene, dplyr::starts_with("Norm_")) %>%
  dplyr::mutate(
    Gene = trimws(as.character(Gene)),
    dplyr::across(-Gene, ~ suppressWarnings(as.numeric(.x)))
  )

samples <- colnames(expr_raw)[-1]
meta_prot <- data.frame(
  Sample = samples,
  Treatment = dplyr::case_when(
    stringr::str_detect(
      samples,
      stringr::regex("KPT[-]?330|XPO1inhibitor", ignore_case = TRUE)
    ) ~ "KPT",
    stringr::str_detect(
      samples,
      stringr::regex("DMSO", ignore_case = TRUE)
    ) ~ "DMSO",
    TRUE ~ NA_character_
  ),
  Compartment = dplyr::case_when(
    stringr::str_detect(
      samples,
      stringr::regex("Nuc", ignore_case = TRUE)
    ) ~ "Nucleus",
    stringr::str_detect(
      samples,
      stringr::regex("Cyto", ignore_case = TRUE)
    ) ~ "Cytoplasm",
    TRUE ~ NA_character_
  ),
  Replicate = dplyr::case_when(
    stringr::str_detect(samples, "D1|X1") ~ "R1",
    stringr::str_detect(samples, "D2|X2") ~ "R2",
    stringr::str_detect(samples, "D3|X3") ~ "R3",
    TRUE ~ NA_character_
  ),
  stringsAsFactors = FALSE,
  row.names = samples
)

keep_samples <- with(
  meta_prot,
  !is.na(Treatment) & !is.na(Compartment) & !is.na(Replicate)
)
meta_prot <- meta_prot[keep_samples, , drop = FALSE]
expr_raw <- expr_raw[, c("Gene", rownames(meta_prot)), drop = FALSE]

cell_counts <- meta_prot %>%
  dplyr::count(Treatment, Compartment, name = "n")
replicate_counts <- meta_prot %>%
  dplyr::count(Treatment, Compartment, Replicate, name = "n")

stopifnot(
  nrow(meta_prot) == 12L,
  !anyDuplicated(meta_prot$Sample),
  !anyNA(meta_prot),
  nrow(cell_counts) == 4L,
  all(cell_counts$n == 3L),
  nrow(replicate_counts) == 12L,
  all(replicate_counts$n == 1L),
  identical(colnames(expr_raw)[-1], rownames(meta_prot))
)

num_cols <- setdiff(colnames(expr_raw), "Gene")
expr_mat_df <- expr_raw %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(num_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

expr_mat_df <- as.data.frame(expr_mat_df)
rownames(expr_mat_df) <- make.unique(expr_mat_df$Gene)
expr_mat <- as.matrix(expr_mat_df[, num_cols, drop = FALSE])
storage.mode(expr_mat) <- "numeric"

enough_obs <- function(matrix, metadata) {
  ok_by_cell <- function(compartment, treatment) {
    index <- which(
      metadata$Compartment == compartment &
        metadata$Treatment == treatment
    )
    rowSums(!is.na(matrix[, index, drop = FALSE])) >= 2L
  }
  
  ok_by_cell("Cytoplasm", "DMSO") &
    ok_by_cell("Cytoplasm", "KPT") &
    ok_by_cell("Nucleus", "DMSO") &
    ok_by_cell("Nucleus", "KPT")
}

expr_mat <- expr_mat[enough_obs(expr_mat, meta_prot), , drop = FALSE]

stopifnot(
  nrow(expr_mat) > 0L,
  ncol(expr_mat) == 12L,
  identical(colnames(expr_mat), rownames(meta_prot))
)

# -----------------------------------------------------------------------------
# 2. Original limma interaction model and significance definition
# -----------------------------------------------------------------------------

is_RPL <- function(x) grepl("^RPL", toupper(x))
is_RPS <- function(x) grepl("^RPS", toupper(x))

cat("RPL proteins entering Fig. 2I:", sum(is_RPL(rownames(expr_mat))), "\n")
cat("RPS proteins entering Fig. 2I:", sum(is_RPS(rownames(expr_mat))), "\n")

Comp <- factor(meta_prot$Compartment, levels = c("Cytoplasm", "Nucleus"))
Trt <- factor(meta_prot$Treatment, levels = c("DMSO", "KPT"))
Rep <- factor(meta_prot$Replicate)

design_dd <- model.matrix(~ Rep + Comp * Trt)
rownames(design_dd) <- rownames(meta_prot)
colnames(design_dd) <- make.names(colnames(design_dd))

stopifnot(
  nrow(design_dd) == ncol(expr_mat),
  identical(rownames(design_dd), colnames(expr_mat)),
  "CompNucleus.TrtKPT" %in% colnames(design_dd)
)

fit0 <- limma::lmFit(expr_mat, design_dd)
contrast_dd <- limma::makeContrasts(
  DeltaDelta = -CompNucleus.TrtKPT,
  levels = design_dd
)
fit_dd <- limma::contrasts.fit(fit0, contrast_dd)

# Preserve the original eBayes analysis. Do not replace with treat().
fit_dd <- limma::eBayes(fit_dd)
tt_dd <- limma::topTable(
  fit_dd,
  coef = "DeltaDelta",
  number = Inf,
  sort.by = "none"
)
res_dd <- data.frame(Gene = rownames(tt_dd), tt_dd, row.names = NULL)

fc_cut <- 1
alpha <- 0.05

res_dd_labeled <- res_dd %>%
  dplyr::mutate(
    neglog10FDR = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    family = dplyr::case_when(
      is_RPL(Gene) ~ "RPL",
      is_RPS(Gene) ~ "RPS",
      TRUE ~ "Other"
    ),
    sig = dplyr::case_when(
      adj.P.Val < alpha & logFC <= -fc_cut ~ "Nuclear Retention",
      adj.P.Val < alpha & logFC >= fc_cut ~ "Cytoplasmic Shift",
      TRUE ~ "NS"
    ),
    is_sig = sig != "NS"
  )

retained_ribosomal <- res_dd_labeled %>%
  dplyr::filter(
    sig == "Nuclear Retention",
    family %in% c("RPL", "RPS")
  )
nuclear_retention_genes <- res_dd_labeled %>%
  dplyr::filter(sig == "Nuclear Retention")

cat("RPL in final Fig. 2I dataset:", sum(res_dd_labeled$family == "RPL"), "\n")
cat("RPS in final Fig. 2I dataset:", sum(res_dd_labeled$family == "RPS"), "\n")
cat("Total Nuclear Retention genes:", nrow(nuclear_retention_genes), "\n")
print(retained_ribosomal %>% dplyr::count(family))
print(
  retained_ribosomal %>%
    dplyr::select(Gene, family, logFC, adj.P.Val) %>%
    dplyr::arrange(family, logFC)
)

# -----------------------------------------------------------------------------
# 3. Figure 2I
# -----------------------------------------------------------------------------

df_nonsig <- res_dd_labeled %>% dplyr::filter(!is_sig)
df_sig_other <- res_dd_labeled %>% dplyr::filter(is_sig, family == "Other")
df_rpl_sig <- res_dd_labeled %>% dplyr::filter(is_sig, family == "RPL")
df_rps_sig <- res_dd_labeled %>% dplyr::filter(is_sig, family == "RPS")

p_2I <- ggplot2::ggplot(
  res_dd_labeled,
  ggplot2::aes(x = logFC, y = neglog10FDR)
) +
  ggplot2::geom_hline(
    yintercept = -log10(alpha),
    linetype = "dotted",
    color = "gray70"
  ) +
  ggplot2::geom_vline(
    xintercept = c(-fc_cut, fc_cut),
    linetype = "dotted",
    color = "gray70"
  ) +
  ggplot2::geom_point(
    data = df_nonsig,
    color = "gray80",
    size = 1.2,
    alpha = 0.25
  ) +
  ggplot2::geom_point(
    data = df_sig_other,
    ggplot2::aes(color = sig),
    shape = 16,
    size = 1.6,
    alpha = 0.9
  ) +
  ggplot2::geom_point(
    data = df_rpl_sig,
    ggplot2::aes(shape = "RPL", fill = sig),
    size = 3.6,
    color = "black",
    stroke = 0.6
  ) +
  ggplot2::geom_point(
    data = df_rps_sig,
    ggplot2::aes(shape = "RPS", color = sig),
    size = 3.6,
    stroke = 1.1
  ) +
  ggplot2::scale_color_manual(
    name = "Shift",
    values = c(
      "Nuclear Retention" = "#D73027",
      "Cytoplasmic Shift" = "#4575B4",
      "NS" = "gray80"
    )
  ) +
  ggplot2::scale_fill_manual(
    name = "Shift",
    values = c(
      "Nuclear Retention" = "#D73027",
      "Cytoplasmic Shift" = "#4575B4",
      "NS" = "gray80"
    )
  ) +
  ggplot2::scale_shape_manual(
    name = "Ribosomal family",
    values = c("RPL" = 24, "RPS" = 7),
    labels = c("RPL", "RPS")
  ) +
  ggplot2::labs(
    x = expression(
      Delta * Delta ~ log[2] * FC ~
        "[(N-C)"[DMSO] - "(N-C)"[KPT] * "]"
    ),
    y = expression(-log[10] ~ FDR)
  ) +
  ggplot2::theme_classic(base_size = 12)


nmd3_save_pdf(
  filename = file.path(
    output_dir,
    "Fig2B_KPT330_Protein_Nuclear_Retention.pdf"
  ),
  plot = p_2I,
  width = 8,
  height = 6,
  device = grDevices::cairo_pdf
)

utils::write.csv(
  res_dd_labeled,
  file.path(output_dir, "Fig2I_all_protein_DeltaDelta_statistics.csv"),
  row.names = FALSE
)
utils::write.csv(
  retained_ribosomal,
  file.path(output_dir, "Fig2I_retained_RPL_RPS.csv"),
  row.names = FALSE
)
utils::write.csv(
  meta_prot,
  file.path(output_dir, "Fig2I_protein_sample_metadata.csv"),
  row.names = FALSE
)





}))

# =============================================================================
# Fig2E
# =============================================================================
nmd3_step("Fig2E", character(), quote({
# Fig2E — analysis below
# Proteomic volcano plots in cytoplasm and nucleus
# KPT-330 vs DMSO
#
# Cytoplasm:
#   RPL/RPS proteins decreased after KPT-330 are highlighted
#
# Nucleus:
#   RPL/RPS proteins increased after KPT-330 are highlighted
#
# Reuses prot_raw generated in Figure 2
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

stopifnot(
  exists("prot_raw"),
  "Gene" %in% colnames(prot_raw)
)

# ------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------

alpha_S3B <- 0.05
fc_threshold_S3B <- 0
y_cap_S3B <- 14
number_labels_S3B <- 10

# Automatic RPL/RPS definitions used in the previous analyses

is_RPL_S3B <- function(x) {
  grepl(
    "^RPL[0-9]+[A-Z]?$",
    toupper(x),
    perl = TRUE
  )
}

is_RPS_S3B <- function(x) {
  grepl(
    "^RPS[0-9]+[A-Z]?$",
    toupper(x),
    perl = TRUE
  )
}

# ------------------------------------------------------------
# 2. Prepare proteomics result for one compartment
# ------------------------------------------------------------

prepare_proteomics_S3B <- function(
    compartment = c(
      "Cytoplasm",
      "Nucleus"
    )
) {
  
  compartment <- match.arg(
    compartment
  )
  
  if (compartment == "Cytoplasm") {
    
    log2fc_column <-
      "Log2FC_XPO1inhibitor_Cyto.DMSO_Cyto"
    
    pvalue_column <-
      "Tukey_XPO1inhibitor_Cyto.DMSO_Cyto"
    
  } else {
    
    log2fc_column <-
      "Log2FC_XPO1inhibitor_Nuc.DMSO_Nuc"
    
    pvalue_column <-
      "Tukey_XPO1inhibitor_Nuc.DMSO_Nuc"
  }
  
  stopifnot(
    all(
      c(
        log2fc_column,
        pvalue_column
      ) %in% colnames(prot_raw)
    )
  )
  
  result <- prot_raw %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != ""
    ) %>%
    dplyr::transmute(
      gene = toupper(
        as.character(Gene)
      ),
      log2FC = suppressWarnings(
        as.numeric(
          .data[[log2fc_column]]
        )
      ),
      pvalue = suppressWarnings(
        as.numeric(
          .data[[pvalue_column]]
        )
      )
    ) %>%
    dplyr::filter(
      is.finite(log2FC),
      is.finite(pvalue)
    ) %>%
    dplyr::mutate(
      pvalue = pmax(
        pvalue,
        .Machine$double.xmin
      ),
      minus_log10_pvalue = pmin(
        -log10(pvalue),
        y_cap_S3B
      ),
      Direction = dplyr::case_when(
        pvalue < alpha_S3B &
          log2FC > fc_threshold_S3B ~
          "KPT-330 increased",
        
        pvalue < alpha_S3B &
          log2FC < -fc_threshold_S3B ~
          "KPT-330 decreased",
        
        TRUE ~ "NS"
      ),
      Family = dplyr::case_when(
        is_RPL_S3B(gene) ~ "RPL",
        is_RPS_S3B(gene) ~ "RPS",
        TRUE ~ "Other"
      ),
      Compartment = compartment
    )
  
  result
}

cytoplasm_data_S3B <- prepare_proteomics_S3B(
  "Cytoplasm"
)

nucleus_data_S3B <- prepare_proteomics_S3B(
  "Nucleus"
)

# ------------------------------------------------------------
# 3. Identify RPL/RPS proteins highlighted in each panel
# ------------------------------------------------------------

cytoplasm_ribosomal_S3B <- cytoplasm_data_S3B %>%
  dplyr::filter(
    Family %in% c(
      "RPL",
      "RPS"
    ),
    log2FC < 0,
    pvalue < alpha_S3B
  )

nucleus_ribosomal_S3B <- nucleus_data_S3B %>%
  dplyr::filter(
    Family %in% c(
      "RPL",
      "RPS"
    ),
    log2FC > 0,
    pvalue < alpha_S3B
  )

cat("\nFigure S3B cytoplasmic RPL/RPS:\n")

print(
  table(
    cytoplasm_ribosomal_S3B$Family
  )
)

cat("Cytoplasmic RPL genes:\n")

print(
  sort(
    unique(
      cytoplasm_ribosomal_S3B$gene[
        cytoplasm_ribosomal_S3B$Family == "RPL"
      ]
    )
  )
)

cat("Cytoplasmic RPS genes:\n")

print(
  sort(
    unique(
      cytoplasm_ribosomal_S3B$gene[
        cytoplasm_ribosomal_S3B$Family == "RPS"
      ]
    )
  )
)

cat("\nFigure S3B nuclear RPL/RPS:\n")

print(
  table(
    nucleus_ribosomal_S3B$Family
  )
)

cat("Nuclear RPL genes:\n")

print(
  sort(
    unique(
      nucleus_ribosomal_S3B$gene[
        nucleus_ribosomal_S3B$Family == "RPL"
      ]
    )
  )
)

cat("Nuclear RPS genes:\n")

print(
  sort(
    unique(
      nucleus_ribosomal_S3B$gene[
        nucleus_ribosomal_S3B$Family == "RPS"
      ]
    )
  )
)

# ------------------------------------------------------------
# 4. Volcano plot function
# ------------------------------------------------------------

make_proteomics_volcano_S3B <- function(
    plot_data,
    ribosomal_data,
    compartment_name,
    show_y_axis = TRUE
) {
  
  top_increased <- plot_data %>%
    dplyr::filter(
      Direction == "KPT-330 increased"
    ) %>%
    dplyr::arrange(
      pvalue,
      dplyr::desc(abs(log2FC))
    ) %>%
    dplyr::slice_head(
      n = number_labels_S3B
    )
  
  top_decreased <- plot_data %>%
    dplyr::filter(
      Direction == "KPT-330 decreased"
    ) %>%
    dplyr::arrange(
      pvalue,
      dplyr::desc(abs(log2FC))
    ) %>%
    dplyr::slice_head(
      n = number_labels_S3B
    )
  
  label_data <- dplyr::bind_rows(
    top_increased,
    top_decreased
  )
  
  cat(
    "\n",
    compartment_name,
    " top labeled proteins:\n",
    sep = ""
  )
  
  print(
    label_data %>%
      dplyr::select(
        gene,
        log2FC,
        pvalue,
        Direction
      )
  )
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = log2FC,
      y = minus_log10_pvalue
    )
  ) +
    
    # Threshold lines
    
    ggplot2::geom_hline(
      yintercept = -log10(
        alpha_S3B
      ),
      linetype = "dotted",
      color = "gray30",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "gray30",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        -1,
        1
      ),
      linetype = "dashed",
      color = "gray60",
      linewidth = 0.4
    ) +
    
    # Nonsignificant proteins
    
    ggplot2::geom_point(
      data = plot_data %>%
        dplyr::filter(
          Direction == "NS"
        ),
      color = "gray80",
      size = 2.5,
      alpha = 0.80
    ) +
    
    # KPT-330-decreased proteins
    
    ggplot2::geom_point(
      data = plot_data %>%
        dplyr::filter(
          Direction ==
            "KPT-330 decreased"
        ),
      color = "#D73027",
      size = 2.5,
      alpha = 0.80
    ) +
    
    # KPT-330-increased proteins
    
    ggplot2::geom_point(
      data = plot_data %>%
        dplyr::filter(
          Direction ==
            "KPT-330 increased"
        ),
      color = "#4575B4",
      size = 2.5,
      alpha = 0.80
    ) +
    
    # RPL/RPS highlights
    
    ggplot2::geom_point(
      data = ribosomal_data,
      ggplot2::aes(
        fill = Family
      ),
      shape = 21,
      color = "black",
      size = 3.4,
      stroke = 0.6,
      alpha = 1
    ) +
    
    # Top labels
    
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(
        label = gene
      ),
      color = "black",
      size = 3.4,
      fontface = "italic",
      max.overlaps = Inf,
      force = 2,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray40",
      min.segment.length = 0,
      seed = 123,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "RPL" = "gold",
        "RPS" = "springgreen3"
      ),
      breaks = c(
        "RPL",
        "RPS"
      ),
      name = "Ribosomal family"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        -6,
        6
      ),
      breaks = seq(
        -6,
        6,
        by = 3
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(
        0,
        y_cap_S3B
      ),
      breaks = seq(
        0,
        y_cap_S3B,
        by = 2
      ),
      labels = function(x) {
        ifelse(
          x == y_cap_S3B,
          paste0(">", y_cap_S3B),
          x
        )
      },
      expand = ggplot2::expansion(
        mult = c(
          0.01,
          0.02
        )
      )
    ) +
    ggplot2::labs(
      title = compartment_name,
      x = expression(
        log[2]~"Fold change (KPT-330 vs DMSO)"
      ),
      y = if (show_y_axis) {
        expression(
          -log[10]~"Tukey p-value"
        )
      } else {
        NULL
      }
    ) +
    ggplot2::theme_classic(
      base_size = 12
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 10,
        color = "black"
      ),
      legend.position = "top",
      panel.grid = ggplot2::element_blank()
    )
  
  p
}

# ------------------------------------------------------------
# 5. Create cytoplasm and nucleus panels
# ------------------------------------------------------------

p_S3B_cytoplasm <- make_proteomics_volcano_S3B(
  plot_data = cytoplasm_data_S3B,
  ribosomal_data = cytoplasm_ribosomal_S3B,
  compartment_name = "Cytoplasm",
  show_y_axis = TRUE
)

p_S3B_nucleus <- make_proteomics_volcano_S3B(
  plot_data = nucleus_data_S3B,
  ribosomal_data = nucleus_ribosomal_S3B,
  compartment_name = "Nucleus",
  show_y_axis = FALSE
)

# ------------------------------------------------------------
# 6. Combine Figure S3B
# ------------------------------------------------------------

p_S3B <- (
  p_S3B_cytoplasm +
    p_S3B_nucleus
) +
  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "top"
  )


# ------------------------------------------------------------
# 7. Save Figure S3B
# ------------------------------------------------------------

nmd3_save_pdf(
  filename = "Fig2E_KPT330_Nuclear_Cytoplasmic_Protein_Volcano.pdf",
  plot = p_S3B,
  width = 14,
  height = 7,
  units = "in",
  device = grDevices::cairo_pdf
)

# ------------------------------------------------------------
# 8. Save plotted RPL/RPS proteins
# ------------------------------------------------------------

ribosomal_results_S3B <- dplyr::bind_rows(
  cytoplasm_ribosomal_S3B,
  nucleus_ribosomal_S3B
)

utils::write.csv(
  ribosomal_results_S3B,
  file = "FigureS3B_Highlighted_RPL_RPS_Proteins.csv",
  row.names = FALSE
)

}))

# =============================================================================
# Fig2F
# =============================================================================
nmd3_step("Fig2F", character(), quote({
# =============================================================================
# PART 4: ALL RPL/RPS HEATMAPS IN NATURAL NUMERIC ORDER
# =============================================================================


# =============================================================================
# NMD3 IN-HOUSE: ALL DETECTED RPL/RPS HEATMAPS — NATURAL NUMERIC ORDER
# =============================================================================
# Source: NMD3_inhouse_RNAseq_DIA_analysis(2).R, Figures S3D / S6E,
# and the expanded All-RPs heatmap sections.
#
# 1) KPT-330 vs DMSO: cytoplasmic RNA + protein (side-by-side)
# 2) siNMD3 vs siCtrl: nuclear + cytoplasmic protein (stacked)
#
# Does not rerun DESeq2, limma, or the original 16,000-line analysis.
#
# Original preprocessing retained:
#   KPT RNA: log2(TPM + 1); KPT protein: normalized intensities
#   siNMD3 protein: conditional log2(x + 1) if compartment median > 50
#   Duplicated proteomics Gene: mean per sample (ignoring NA)
#   Complete 6-sample QC and SD > 0 per heatmap (all 12 for siNMD3)
#   Row-wise z-score, clipped to [-2, 2]
#
# Only new gene arrangement: ALL eligible RPLs in natural-number order,
# followed by ALL eligible RPSs in natural-number order.
# Examples: RPLP0, RPLP1, RPLP2, RPL3, RPL4, ..., RPL10, RPL10A,
#           RPL11, ... / RPS2, RPS3, RPS3A, RPS4X, ...
# RPSA (no numeric designation), if present, is placed after numbered RPSs.
# No legend or sample/gene-family annotations.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(pheatmap)
  library(grid)
  library(gridExtra)
})

# Both inputs are guaranteed by the main run; do not reload or restore anything.
stopifnot(exists("expr_matrix"), exists("prot_raw"))

stopifnot(
  "Length" %in% colnames(expr_matrix),
  "Gene" %in% colnames(prot_raw)
)

rp_out <- file.path(getwd(), "RPL_RPS_Natural_Numeric_Order")
dir.create(rp_out, recursive = TRUE, showWarnings = FALSE)

cat("RNA input: ", nrow(expr_matrix), " genes; ", ncol(expr_matrix), " columns\n", sep = "")
cat("Protein input: ", nrow(prot_raw), " rows; ", ncol(prot_raw), " columns\n", sep = "")
cat("Output: ", normalizePath(rp_out), "\n", sep = "")

# ============================================================================
# 1. Gene-family definitions: same regex as the expanded All-RPs script
# ============================================================================

rp_is_rpl <- function(x) {
  x <- toupper(as.character(x))
  grepl("^RPL(P[0-2]|[0-9]+(A|B|AL|L[0-9]*)?)$", x)
}

rp_is_rps <- function(x) {
  x <- toupper(as.character(x))
  grepl("^RPS(A|[0-9]+[A-Z]?|4Y[12])$", x)
}

# Natural numeric sorting instead of lexical sorting / the old hand-picked list.
# For example: RPLP0 -> RPLP1 -> RPLP2 -> RPL3 -> ... -> RPL10 -> RPL10A;
# RPS2 -> RPS3 -> RPS3A -> RPS4X -> ... -> RPS10 -> ... -> RPSA.
# The prefix families are kept separate, and non-numbered RPSA comes last.
rp_sort_numeric <- function(genes) {
  genes <- sort(unique(toupper(as.character(genes))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) return(character())
  
  group <- ifelse(rp_is_rpl(genes), 1L,
                  ifelse(rp_is_rps(genes), 2L, 3L))
  digit_end <- regexpr("[0-9]+", genes)
  digit_length <- attr(digit_end, "match.length")
  has_digit <- digit_end > 0L
  number_text <- rep(NA_character_, length(genes))
  number_text[has_digit] <- substr(
    genes[has_digit],
    digit_end[has_digit],
    digit_end[has_digit] + digit_length[has_digit] - 1L
  )
  number <- suppressWarnings(as.numeric(number_text))
  number[!is.finite(number)] <- Inf
  
  # For equal numbers (e.g. RPL10 and RPL10A), place the unsuffixed
  # form first, then suffix variants alphabetically.
  suffix <- ifelse(
    digit_end > 0,
    substring(genes, digit_end + digit_length),
    genes
  )
  suffix <- ifelse(grepl("^RPLP[0-2]$", genes), "", suffix)
  
  genes[order(group, number, nchar(suffix), suffix, genes, na.last = TRUE)]
}

rp_valid <- function(x) {
  all(is.finite(x)) && is.finite(stats::sd(x)) && stats::sd(x) > 0
}

rp_zscore <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  z <- t(base::scale(t(mat), center = TRUE, scale = TRUE))
  stopifnot(
    is.matrix(z),
    identical(dim(z), dim(mat)),
    identical(rownames(z), rownames(mat)),
    identical(colnames(z), colnames(mat)),
    all(is.finite(z))
  )
  # Subset assignment preserves row names, column names and 2D shape.
  z[z >  2] <-  2
  z[z < -2] <- -2
  stopifnot(is.matrix(z), identical(dim(z), dim(mat)))
  z
}

rp_colors <- grDevices::colorRampPalette(
  c("#4575B4", "#F7F7F7", "#D73027")
)(100)
rp_breaks <- seq(-2, 2, length.out = 101)

# Clean, deduplicate and aggregate all original protein sample intensities.
# Same first-Gene-before-semicolon and per-sample mean rules as the source.
rp_protein_matrix <- function(prot_df, sample_cols) {
  stopifnot(all(sample_cols %in% colnames(prot_df)))
  tab <- prot_df %>%
    dplyr::mutate(gene = toupper(trimws(sub(";.*$", "", as.character(Gene))))) %>%
    dplyr::filter(!is.na(gene), nzchar(gene)) %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(sample_cols),
      ~ suppressWarnings(as.numeric(.x))
    )) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dplyr::across(
      dplyr::all_of(sample_cols),
      ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
    ), .groups = "drop")
  
  out <- as.matrix(tab[, sample_cols, drop = FALSE])
  rownames(out) <- tab$gene
  storage.mode(out) <- "double"
  stopifnot(!anyDuplicated(rownames(out)))
  out
}

rp_save_matrix <- function(mat, filename) {
  utils::write.csv(
    data.frame(Gene_or_sample = rownames(mat), mat, check.names = FALSE),
    file.path(rp_out, filename), row.names = FALSE
  )
}

rp_save_genelist <- function(genes, filename) {
  utils::write.csv(
    data.frame(Order = seq_along(genes), Gene = genes,
               Family = ifelse(rp_is_rpl(genes), "RPL", "RPS")),
    file.path(rp_out, filename), row.names = FALSE
  )
}

# ============================================================================
# 2. XPO1i vs DMSO — cytoplasmic RNA and protein heatmaps
# ============================================================================

rp_xpo_rna_cols <- c("D1C", "D2C", "D3C", "X1C", "X2C", "X3C")
rp_xpo_prot_cols <- c(
  paste0("Norm_DMSO_Cyto_D", 1:3, "C"),
  paste0("Norm_XPO1inhibitor_Cyto_X", 1:3, "C")
)
stopifnot(
  all(rp_xpo_rna_cols %in% colnames(expr_matrix)),
  all(rp_xpo_prot_cols %in% colnames(prot_raw))
)

# Original RNA normalization: TPM on all genes with valid gene length,
# followed by log2(TPM+1), with the six cytoplasmic columns.
rp_gene_length_kb <- suppressWarnings(as.numeric(expr_matrix$Length)) / 1000
rp_length_ok <- is.finite(rp_gene_length_kb) & rp_gene_length_kb > 0
rp_xpo_counts <- as.matrix(expr_matrix[rp_length_ok, rp_xpo_rna_cols, drop = FALSE])
storage.mode(rp_xpo_counts) <- "double"
stopifnot(all(is.finite(rp_xpo_counts)), all(rp_xpo_counts >= 0))

rp_xpo_rpk <- sweep(rp_xpo_counts, 1, rp_gene_length_kb[rp_length_ok], "/")
rp_xpo_library <- colSums(rp_xpo_rpk)
stopifnot(all(is.finite(rp_xpo_library)), all(rp_xpo_library > 0))
rp_xpo_tpm <- sweep(rp_xpo_rpk, 2, rp_xpo_library, "/") * 1e6
rp_xpo_rna_log <- log2(rp_xpo_tpm + 1)

# Original XPO1i heatmap used normalized protein intensities as-is.
rp_xpo_protein <- rp_protein_matrix(prot_raw, rp_xpo_prot_cols)

rp_xpo_candidates <- base::intersect(
  rownames(rp_xpo_rna_log), rownames(rp_xpo_protein)
)
rp_xpo_candidates <- rp_xpo_candidates[
  rp_is_rpl(rp_xpo_candidates) | rp_is_rps(rp_xpo_candidates)
]
rp_xpo_keep <- rp_xpo_candidates[vapply(rp_xpo_candidates, function(g) {
  rp_valid(rp_xpo_rna_log[g, ]) && rp_valid(rp_xpo_protein[g, ])
}, logical(1))]
rp_xpo_genes <- rp_sort_numeric(rp_xpo_keep)

stopifnot(length(rp_xpo_genes) > 1L,
          any(rp_is_rpl(rp_xpo_genes)), any(rp_is_rps(rp_xpo_genes)))

rp_xpo_rna_z <- rp_zscore(rp_xpo_rna_log[rp_xpo_genes, , drop = FALSE])
rp_xpo_protein_z <- rp_zscore(rp_xpo_protein[rp_xpo_genes, , drop = FALSE])
stopifnot(identical(rownames(rp_xpo_rna_z), rownames(rp_xpo_protein_z)))

rp_xpo_gap <- sum(rp_is_rpl(rp_xpo_genes))

rp_draw_xpo <- function(mat, panel_title) {
  stopifnot(is.matrix(mat), ncol(mat) == 6L,
            is.numeric(rp_xpo_gap), length(rp_xpo_gap) == 1L,
            rp_xpo_gap >= 1L, rp_xpo_gap < nrow(mat))
  message(panel_title, ": ", nrow(mat), " RPL/RPS genes x ", ncol(mat),
          " samples; RPL/RPS boundary after row ", rp_xpo_gap)
  pheatmap::pheatmap(
    mat = mat,
    color = rp_colors,
    breaks = rp_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 7.5,
    border_color = "grey90",
    gaps_row = rp_xpo_gap,
    gaps_col = 3, # valid now: exactly six sample columns, group boundary after column 3
    legend = FALSE,
    annotation_legend = FALSE,
    main = panel_title,
    silent = TRUE
  )
}

rp_xpo_h_rna <- rp_draw_xpo(rp_xpo_rna_z, "RNA")
rp_xpo_h_prot <- rp_draw_xpo(rp_xpo_protein_z, "Protein")

# Align the two halves.
rp_xpo_shared_widths <- grid::unit.pmax(
  rp_xpo_h_rna$gtable$widths,
  rp_xpo_h_prot$gtable$widths
)
rp_xpo_h_rna$gtable$widths <- rp_xpo_shared_widths
rp_xpo_h_prot$gtable$widths <- rp_xpo_shared_widths

rp_xpo_plot <- gridExtra::arrangeGrob(
  rp_xpo_h_rna$gtable, rp_xpo_h_prot$gtable,
  ncol = 2, widths = c(1, 1)
)


rp_xpo_pdf <- file.path(rp_out, "Fig2F_RPL_RPS_RNA_Protein_Heatmaps.pdf")
nmd3_cairo_pdf(
  rp_xpo_pdf,
  width = 12,
  height = max(9, length(rp_xpo_genes) * 0.18)
)
grid::grid.newpage()
grid::grid.draw(rp_xpo_plot)
nmd3_dev_off()

rp_save_genelist(rp_xpo_genes, "XPO1i_RPL_RPS_GeneOrder.csv")
rp_save_matrix(rp_xpo_rna_z, "XPO1i_RNA_RowZ.csv")
rp_save_matrix(rp_xpo_protein_z, "XPO1i_Protein_RowZ.csv")

cat("\nXPO1i: ", sum(rp_is_rpl(rp_xpo_genes)), " RPL; ",
    sum(rp_is_rps(rp_xpo_genes)), " RPS\n", sep = "")
cat("XPO1i excluded after QC: ",
    paste(setdiff(rp_xpo_candidates, rp_xpo_keep), collapse = ", "), "\n", sep = "")

}))

# =============================================================================
# Fig3A
# =============================================================================
nmd3_step("Fig3A", character(), quote({
# Fig3A — analysis below
# RNA-protein four-quadrant analysis
# siNMD3 vs siCTRL and KPT-330 vs DMSO
# Cytoplasmic fraction
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(limma)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

alpha_thr <- 0.05

# ============================================================
# 1. Rebuild the complete RNA count matrix
# expr_matrix and prot_raw were loaded in Figure 2
# ============================================================

sample_cols <- grep(
  "^[CNDX][0-9]+[CN]$",
  colnames(expr_matrix),
  value = TRUE
)

count_matrix_all <- expr_matrix[
  ,
  sample_cols,
  drop = FALSE
]

# ============================================================
# 2. RNA differential-expression function
# ============================================================

run_deseq2_cyto <- function(
    count_mat,
    samples,
    groups,
    coef_name
) {
  
  counts <- as.matrix(
    count_mat[
      ,
      samples,
      drop = FALSE
    ]
  )
  
  storage.mode(counts) <- "integer"
  
  coldata <- data.frame(
    row.names = samples,
    condition = groups
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )
  
  dds <- DESeq(
    dds,
    quiet = TRUE
  )
  
  # Raw P values
  
  res_raw <- results(
    dds,
    name = coef_name
  )
  
  # apeglm-shrunken log2FC
  
  res_shrunk <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )
  
  data.frame(
    gene = toupper(rownames(res_raw)),
    log2FC_RNA = res_shrunk$log2FoldChange,
    pval_RNA = res_raw$pvalue,
    padj_RNA = res_raw$padj
  ) %>%
    mutate(
      pval_RNA = coalesce(pval_RNA, 1),
      padj_RNA = coalesce(padj_RNA, 1)
    ) %>%
    filter(
      is.finite(log2FC_RNA)
    )
}

# ============================================================
# 3. Protein differential-expression function
# ============================================================

run_limma_3v3 <- function(
    prot_df,
    col_ctrl,
    col_treat
) {
  
  stopifnot(
    length(col_ctrl) == 3,
    length(col_treat) == 3
  )
  
  expr_cols <- c(
    col_ctrl,
    col_treat
  )
  
  tab <- prot_df %>%
    filter(!is.na(Gene)) %>%
    mutate(
      across(
        all_of(expr_cols),
        ~ suppressWarnings(as.numeric(.x))
      )
    ) %>%
    group_by(Gene) %>%
    summarise(
      across(
        all_of(expr_cols),
        ~ if (all(is.na(.x))) {
          NA_real_
        } else {
          median(.x, na.rm = TRUE)
        }
      ),
      .groups = "drop"
    ) %>%
    filter(
      rowSums(
        is.na(
          across(all_of(expr_cols))
        )
      ) == 0
    )
  
  X <- as.matrix(
    tab[
      ,
      expr_cols,
      drop = FALSE
    ]
  )
  
  rownames(X) <- tab$Gene
  
  # Preserve original analysis
  
  if (median(X, na.rm = TRUE) > 50) {
    X <- log2(X + 1)
  }
  
  group <- factor(
    c(
      rep("CTRL", 3),
      rep("TREAT", 3)
    ),
    levels = c(
      "CTRL",
      "TREAT"
    )
  )
  
  design <- model.matrix(
    ~ group
  )
  
  fit <- lmFit(
    X,
    design
  )
  
  contrast <- makeContrasts(
    TreatVsCtrl = groupTREAT,
    levels = design
  )
  
  fit2 <- contrasts.fit(
    fit,
    contrast
  )
  
  fit2 <- eBayes(
    fit2,
    robust = TRUE,
    trend = TRUE
  )
  
  tt <- topTable(
    fit2,
    coef = "TreatVsCtrl",
    number = Inf,
    sort.by = "none"
  )
  
  data.frame(
    gene = toupper(rownames(tt)),
    log2FC_Prot = tt$logFC,
    pval_Prot = tt$P.Value,
    padj_Prot = tt$adj.P.Val
  ) %>%
    mutate(
      pval_Prot = coalesce(pval_Prot, 1),
      padj_Prot = coalesce(padj_Prot, 1)
    )
}

# ============================================================
# 4. siNMD3 vs siCTRL RNA
# ============================================================

si_samples <- grep(
  "^[CN][0-9]+C$",
  colnames(count_matrix_all),
  value = TRUE
)

si_group <- factor(
  ifelse(
    grepl("^C", si_samples),
    "siCTRL",
    "siNMD3"
  ),
  levels = c(
    "siCTRL",
    "siNMD3"
  )
)

print(si_samples)
print(table(si_group))

stopifnot(
  sum(si_group == "siCTRL") == 3,
  sum(si_group == "siNMD3") == 3
)

rna_si <- run_deseq2_cyto(
  count_mat = count_matrix_all,
  samples = si_samples,
  groups = si_group,
  coef_name = "condition_siNMD3_vs_siCTRL"
)

# ============================================================
# 5. siNMD3 vs siCTRL protein
# ============================================================

si_ctrl_cols <- grep(
  "^Norm_siRNAxCtrl_Cyto_",
  colnames(prot_raw),
  value = TRUE
)

si_nmd3_cols <- grep(
  "^Norm_siRNAxNMD3xKD_Cyto_",
  colnames(prot_raw),
  value = TRUE
)

stopifnot(
  length(si_ctrl_cols) == 3,
  length(si_nmd3_cols) == 3
)

prot_si <- run_limma_3v3(
  prot_df = prot_raw,
  col_ctrl = si_ctrl_cols,
  col_treat = si_nmd3_cols
)

multi_si <- inner_join(
  rna_si,
  prot_si,
  by = "gene"
) %>%
  filter(
    is.finite(log2FC_RNA),
    is.finite(log2FC_Prot)
  )

# ============================================================
# 6. KPT-330 vs DMSO RNA
# ============================================================

kpt_samples <- grep(
  "^[DX][0-9]+C$",
  colnames(count_matrix_all),
  value = TRUE
)

kpt_group <- factor(
  ifelse(
    grepl("^D", kpt_samples),
    "DMSO",
    "KPT.330"
  ),
  levels = c(
    "DMSO",
    "KPT.330"
  )
)

print(kpt_samples)
print(table(kpt_group))

stopifnot(
  sum(kpt_group == "DMSO") == 3,
  sum(kpt_group == "KPT.330") == 3
)

rna_kpt <- run_deseq2_cyto(
  count_mat = count_matrix_all,
  samples = kpt_samples,
  groups = kpt_group,
  coef_name = "condition_KPT.330_vs_DMSO"
)

# ============================================================
# 7. KPT-330 vs DMSO protein
# ============================================================

dms_cols <- c(
  "Norm_DMSO_Cyto_D1C",
  "Norm_DMSO_Cyto_D2C",
  "Norm_DMSO_Cyto_D3C"
)

kpt_cols <- c(
  "Norm_XPO1inhibitor_Cyto_X1C",
  "Norm_XPO1inhibitor_Cyto_X2C",
  "Norm_XPO1inhibitor_Cyto_X3C"
)

stopifnot(
  all(dms_cols %in% colnames(prot_raw)),
  all(kpt_cols %in% colnames(prot_raw))
)

prot_kpt <- run_limma_3v3(
  prot_df = prot_raw,
  col_ctrl = dms_cols,
  col_treat = kpt_cols
)

multi_kpt <- inner_join(
  rna_kpt,
  prot_kpt,
  by = "gene"
) %>%
  filter(
    is.finite(log2FC_RNA),
    is.finite(log2FC_Prot)
  )

# ============================================================
# 8. Common RNA-protein gene universe
# ============================================================

common_genes <- intersect(
  multi_si$gene,
  multi_kpt$gene
)

multi_si_common <- multi_si %>%
  filter(gene %in% common_genes)

multi_kpt_common <- multi_kpt %>%
  filter(gene %in% common_genes)

cat(
  "Common RNA-protein universe:",
  length(common_genes),
  "\n"
)

# ============================================================
# 9. Four-quadrant classification
#
# Original criterion:
# RNA raw P < 0.05
# Protein raw P < 0.05
# No additional log2FC threshold
# ============================================================

get_quadrant <- function(df) {
  
  df %>%
    filter(
      pval_RNA < alpha_thr,
      pval_Prot < alpha_thr
    ) %>%
    mutate(
      Quadrant = case_when(
        
        log2FC_RNA > 0 &
          log2FC_Prot > 0 ~ "High-High",
        
        log2FC_RNA < 0 &
          log2FC_Prot < 0 ~ "Low-Low",
        
        log2FC_RNA > 0 &
          log2FC_Prot < 0 ~ "High-Low",
        
        log2FC_RNA < 0 &
          log2FC_Prot > 0 ~ "Low-High",
        
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Quadrant))
}

quad_si <- get_quadrant(
  multi_si_common
)

quad_kpt <- get_quadrant(
  multi_kpt_common
)

# ============================================================
# 10. Quadrant summary
# ============================================================

summary_si <- quad_si %>%
  dplyr::count(
    Quadrant,
    name = "n"
  ) %>%
  dplyr::mutate(
    pct = n / sum(n),
    Contrast = "siNMD3 vs siCTRL"
  )

summary_kpt <- quad_kpt %>%
  dplyr::count(
    Quadrant,
    name = "n"
  ) %>%
  dplyr::mutate(
    pct = n / sum(n),
    Contrast = "KPT-330 vs DMSO"
  )

quadrant_summary <- dplyr::bind_rows(
  summary_si,
  summary_kpt
) %>%
  tidyr::complete(
    Contrast,
    Quadrant,
    fill = list(
      n = 0,
      pct = 0
    )
  )

print(quadrant_summary)

cat(
  "siNMD3 significant RNA-protein genes:",
  nrow(quad_si),
  "\n"
)

cat(
  "KPT-330 significant RNA-protein genes:",
  nrow(quad_kpt),
  "\n"
)

# ============================================================
# 11. Plot
# ============================================================

quad_cols <- c(
  "High-High" = "#E41A1C",
  "Low-Low"   = "#377EB8",
  "High-Low"  = "#FF7F00",
  "Low-High"  = "#4DAF4A"
)

# siNMD3 vs siCTRL

p_4A_si <- ggplot(
  quad_si,
  aes(
    x = log2FC_RNA,
    y = log2FC_Prot,
    color = Quadrant
  )
) +
  geom_point(
    size = 1.6,
    alpha = 0.65
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = quad_cols
  ) +
  labs(
    title = "siNMD3 vs siCTRL",
    x = "RNA log2FC",
    y = "Protein log2FC"
  ) +
  theme_classic(base_size = 11)

# KPT-330 vs DMSO

p_4A_kpt <- ggplot(
  quad_kpt,
  aes(
    x = log2FC_RNA,
    y = log2FC_Prot,
    color = Quadrant
  )
) +
  geom_point(
    size = 1.6,
    alpha = 0.65
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = quad_cols
  ) +
  labs(
    title = "KPT-330 vs DMSO",
    x = "RNA log2FC",
    y = "Protein log2FC"
  ) +
  theme_classic(base_size = 11)



# ============================================================
# 12. Save PDFs
# ============================================================

nmd3_save_pdf(
  "Fig3A_left_KPT330_RNA_Protein_Quadrants.pdf",
  p_4A_kpt,
  width = 7,
  height = 5
)

nmd3_save_pdf(
  "Fig3A_right_siNMD3_RNA_Protein_Quadrants.pdf",
  p_4A_si,
  width = 7,
  height = 5
)


}))

# =============================================================================
# Prepare_Fig3B_S6
# =============================================================================
nmd3_step("Prepare_Fig3B_S6", character(), quote({
# Prepare_Fig3B_S6 — analysis below
# siNMD3 versus siCTRL in nucleus and cytoplasm
# S5B: RNA-seq
# S5C: Proteomics
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ------------------------------------------------------------
# 0. Parameters
# ------------------------------------------------------------

rna_fdr_threshold_S5BC <- 0.05
rna_lfc_threshold_S5BC <- 1.0
rna_label_number_S5BC <- 5

protein_p_threshold_S5BC <- 0.05
protein_label_number_S5BC <- 10

volcano_colors_S5BC <- c(
  "siNMD3 increased" = "#4575B4",
  "siNMD3 decreased" = "#D73027",
  "NS" = "gray80"
)

# ============================================================
# Prepare_Fig3B_S6 — analysis below
# ============================================================

# ------------------------------------------------------------
# 1. Reuse the RNA count matrix
# ------------------------------------------------------------

stopifnot(
  exists("expr_matrix"),
  "Length" %in% colnames(expr_matrix)
)

rna_counts_all_S5B <- expr_matrix[
  ,
  setdiff(
    colnames(expr_matrix),
    "Length"
  ),
  drop = FALSE
]

# C prefix: siCTRL
# N prefix: siNMD3
# C suffix: Cytoplasm
# N suffix: Nucleus

rna_sample_info_S5B <- data.frame(
  Sample = colnames(rna_counts_all_S5B),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(
    grepl(
      "^[CN][0-9]+[CN]$",
      Sample
    )
  ) %>%
  dplyr::mutate(
    Prefix = substr(
      Sample,
      1,
      1
    ),
    Suffix = substr(
      Sample,
      nchar(Sample),
      nchar(Sample)
    ),
    Condition = dplyr::case_when(
      Prefix == "C" ~ "siCTRL",
      Prefix == "N" ~ "siNMD3",
      TRUE ~ NA_character_
    ),
    Compartment = dplyr::case_when(
      Suffix == "C" ~ "Cytoplasm",
      Suffix == "N" ~ "Nucleus",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(
    !is.na(Condition),
    !is.na(Compartment)
  )

cat("RNA-seq samples:\n")

print(
  table(
    rna_sample_info_S5B$Condition,
    rna_sample_info_S5B$Compartment
  )
)

stopifnot(
  all(
    table(
      rna_sample_info_S5B$Condition,
      rna_sample_info_S5B$Compartment
    ) >= 2
  )
)

# ------------------------------------------------------------
# 2. DESeq2 function
# Positive log2FC = higher in siNMD3
# ------------------------------------------------------------

run_deseq_S5B <- function(
    compartment_name
) {
  
  sample_table <- rna_sample_info_S5B %>%
    dplyr::filter(
      Compartment == compartment_name
    )
  
  sample_names <- sample_table$Sample
  
  count_matrix <- as.matrix(
    rna_counts_all_S5B[
      ,
      sample_names,
      drop = FALSE
    ]
  )
  
  storage.mode(count_matrix) <- "integer"
  
  coldata <- data.frame(
    condition = factor(
      sample_table$Condition,
      levels = c(
        "siCTRL",
        "siNMD3"
      )
    ),
    row.names = sample_names
  )
  
  stopifnot(
    identical(
      colnames(count_matrix),
      rownames(coldata)
    )
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
  )
  
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )
  
  result <- DESeq2::results(
    dds,
    contrast = c(
      "condition",
      "siNMD3",
      "siCTRL"
    )
  )
  
  result_table <- as.data.frame(
    result
  ) %>%
    tibble::rownames_to_column(
      "gene"
    ) %>%
    dplyr::transmute(
      gene = toupper(
        as.character(gene)
      ),
      baseMean = as.numeric(baseMean),
      log2FC = as.numeric(log2FoldChange),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj),
      Compartment = compartment_name
    ) %>%
    dplyr::mutate(
      log2FC = ifelse(
        is.finite(log2FC),
        log2FC,
        0
      ),
      pvalue = ifelse(
        is.finite(pvalue),
        pvalue,
        1
      ),
      padj = ifelse(
        is.finite(padj),
        padj,
        1
      ),
      padj_plot = pmax(
        padj,
        .Machine$double.xmin
      ),
      neg_log10_padj = -log10(
        padj_plot
      ),
      Direction = dplyr::case_when(
        padj < rna_fdr_threshold_S5BC &
          log2FC >= rna_lfc_threshold_S5BC ~
          "siNMD3 increased",
        
        padj < rna_fdr_threshold_S5BC &
          log2FC <= -rna_lfc_threshold_S5BC ~
          "siNMD3 decreased",
        
        TRUE ~ "NS"
      ),
      Direction = factor(
        Direction,
        levels = c(
          "siNMD3 increased",
          "siNMD3 decreased",
          "NS"
        )
      )
    )
  
  result_table
}

rna_nucleus_S5B <- run_deseq_S5B(
  "Nucleus"
)

rna_cytoplasm_S5B <- run_deseq_S5B(
  "Cytoplasm"
)

# Prepare_Fig3B_S6 — analysis below
# Nucleus-cytoplasm overlap of treatment-responsive genes
#
#
# Definitions:
# "promoted"   = expression decreases after inhibition/KD
# "suppressed" = expression increases after inhibition/KD
#
# RNA:
# abs(log2FC) > 1 and DESeq2 BH-FDR < 0.05
#
# Protein:
# direction of log2FC and Tukey P < 0.05
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(patchwork)
})

# ------------------------------------------------------------
# 0. Required upstream objects
# ------------------------------------------------------------

stopifnot(
  exists("expr_matrix"),
  exists("rna_nucleus_S5B"),
  exists("rna_cytoplasm_S5B")
)

if (exists("prot_raw")) {
  
  protein_source_S5D <- prot_raw
  
} else if (exists("protein_data")) {
  
  protein_source_S5D <- protein_data
  
} else {
  
  stop(
    paste0(
      "Neither prot_raw nor protein_data exists. ",
      "Load the original proteomics report first."
    )
  )
}

stopifnot(
  "Gene" %in%
    colnames(protein_source_S5D)
)

rna_lfc_threshold_S5D <- 1
rna_fdr_threshold_S5D <- 0.05
protein_p_threshold_S5D <- 0.05


# ============================================================
# 1. NMD3 RNA sets
#
# Reuse the completed Figure S5B DESeq2 results.
# Positive log2FC = higher in siNMD3 than siCTRL.
# ============================================================

required_rna_columns_S5D <- c(
  "gene",
  "log2FC",
  "padj"
)

stopifnot(
  all(
    required_rna_columns_S5D %in%
      colnames(rna_nucleus_S5B)
  ),
  all(
    required_rna_columns_S5D %in%
      colnames(rna_cytoplasm_S5B)
  )
)

extract_rna_set_S5D <- function(
    result_table,
    direction
) {
  
  clean_table <- result_table %>%
    dplyr::transmute(
      gene = toupper(
        trimws(
          as.character(gene)
        )
      ),
      log2FC = suppressWarnings(
        as.numeric(log2FC)
      ),
      padj = suppressWarnings(
        as.numeric(padj)
      )
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      is.finite(log2FC),
      is.finite(padj),
      padj < rna_fdr_threshold_S5D
    )
  
  if (identical(direction, "promoted")) {
    
    selected <- clean_table %>%
      dplyr::filter(
        log2FC <
          -rna_lfc_threshold_S5D
      )
    
  } else if (
    identical(direction, "suppressed")
  ) {
    
    selected <- clean_table %>%
      dplyr::filter(
        log2FC >
          rna_lfc_threshold_S5D
      )
    
  } else {
    
    stop(
      "direction must be promoted or suppressed."
    )
  }
  
  sort(
    unique(
      selected$gene
    )
  )
}

RNA_NMD3_promoted_Nucleus_S5D <-
  extract_rna_set_S5D(
    rna_nucleus_S5B,
    "promoted"
  )

RNA_NMD3_promoted_Cytoplasm_S5D <-
  extract_rna_set_S5D(
    rna_cytoplasm_S5B,
    "promoted"
  )

RNA_NMD3_suppressed_Nucleus_S5D <-
  extract_rna_set_S5D(
    rna_nucleus_S5B,
    "suppressed"
  )

RNA_NMD3_suppressed_Cytoplasm_S5D <-
  extract_rna_set_S5D(
    rna_cytoplasm_S5B,
    "suppressed"
  )


# ============================================================
# 2. XPO1 RNA results
#
# Reconstruct the Figure 2 comparison separately in nucleus
# and cytoplasm.
#
# Positive log2FC = higher in KPT-330 than DMSO.
# ============================================================

run_xpo1_rna_S5D <- function(
    compartment_name
) {
  
  suffix <- if (
    identical(
      compartment_name,
      "Nucleus"
    )
  ) {
    "N"
  } else if (
    identical(
      compartment_name,
      "Cytoplasm"
    )
  ) {
    "C"
  } else {
    stop(
      "Unknown compartment."
    )
  }
  
  sample_pattern <- paste0(
    "^[DX][0-9]+",
    suffix,
    "$"
  )
  
  samples <- grep(
    sample_pattern,
    colnames(expr_matrix),
    value = TRUE
  )
  
  stopifnot(
    length(samples) >= 4
  )
  
  condition <- factor(
    ifelse(
      substr(
        samples,
        1,
        1
      ) == "D",
      "DMSO",
      "KPT-330"
    ),
    levels = c(
      "DMSO",
      "KPT-330"
    )
  )
  
  stopifnot(
    all(
      table(condition) >= 2
    )
  )
  
  count_matrix <- as.matrix(
    expr_matrix[
      ,
      samples,
      drop = FALSE
    ]
  )
  
  storage.mode(
    count_matrix
  ) <- "integer"
  
  stopifnot(
    all(
      is.finite(count_matrix)
    ),
    all(
      count_matrix >= 0
    )
  )
  
  coldata <- data.frame(
    condition = condition,
    row.names = samples
  )
  
  stopifnot(
    identical(
      colnames(count_matrix),
      rownames(coldata)
    )
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
  )
  
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )
  
  # Preserve the original Figure S5D definition:
  # use the standard, unshrunken DESeq2 log2 fold change.
  result_unshrunken <- DESeq2::results(
    dds,
    contrast = c(
      "condition",
      "KPT-330",
      "DMSO"
    )
  )
  
  as.data.frame(
    result_unshrunken
  ) %>%
    tibble::rownames_to_column(
      "gene"
    ) %>%
    dplyr::transmute(
      gene = toupper(
        trimws(
          as.character(gene)
        )
      ),
      log2FC = suppressWarnings(
        as.numeric(log2FoldChange)
      ),
      padj = suppressWarnings(
        as.numeric(padj)
      ),
      Compartment = compartment_name
    )
}

rna_xpo1_nucleus_S5D <-
  run_xpo1_rna_S5D(
    "Nucleus"
  )

rna_xpo1_cytoplasm_S5D <-
  run_xpo1_rna_S5D(
    "Cytoplasm"
  )

RNA_XPO1_promoted_Nucleus_S5D <-
  extract_rna_set_S5D(
    rna_xpo1_nucleus_S5D,
    "promoted"
  )

RNA_XPO1_promoted_Cytoplasm_S5D <-
  extract_rna_set_S5D(
    rna_xpo1_cytoplasm_S5D,
    "promoted"
  )

RNA_XPO1_suppressed_Nucleus_S5D <-
  extract_rna_set_S5D(
    rna_xpo1_nucleus_S5D,
    "suppressed"
  )

RNA_XPO1_suppressed_Cytoplasm_S5D <-
  extract_rna_set_S5D(
    rna_xpo1_cytoplasm_S5D,
    "suppressed"
  )


# ============================================================
# 3. Protein comparisons
#
# Only same-compartment treatment-control comparisons are used.
# Positive log2FC = higher after inhibition/KD.
# ============================================================

protein_columns_S5D <- list(
  
  NMD3_Nucleus = c(
    log2FC =
      "Log2FC_siRNAxNMD3xKD_Nuc.siRNAxCtrl_Nuc",
    pvalue =
      "Tukey_siRNAxNMD3xKD_Nuc.siRNAxCtrl_Nuc"
  ),
  
  NMD3_Cytoplasm = c(
    log2FC =
      "Log2FC_siRNAxNMD3xKD_Cyto.siRNAxCtrl_Cyto",
    pvalue =
      "Tukey_siRNAxNMD3xKD_Cyto.siRNAxCtrl_Cyto"
  ),
  
  XPO1_Nucleus = c(
    log2FC =
      "Log2FC_XPO1inhibitor_Nuc.DMSO_Nuc",
    pvalue =
      "Tukey_XPO1inhibitor_Nuc.DMSO_Nuc"
  ),
  
  XPO1_Cytoplasm = c(
    log2FC =
      "Log2FC_XPO1inhibitor_Cyto.DMSO_Cyto",
    pvalue =
      "Tukey_XPO1inhibitor_Cyto.DMSO_Cyto"
  )
)

required_protein_columns_S5D <- unique(
  unname(
    unlist(
      protein_columns_S5D
    )
  )
)

stopifnot(
  all(
    required_protein_columns_S5D %in%
      colnames(protein_source_S5D)
  )
)

extract_protein_result_S5D <- function(
    comparison_name,
    compartment_name
) {
  
  comparison_columns <-
    protein_columns_S5D[[comparison_name]]
  
  protein_source_S5D %>%
    dplyr::transmute(
      gene = toupper(
        trimws(
          as.character(Gene)
        )
      ),
      log2FC = suppressWarnings(
        as.numeric(
          .data[[
            comparison_columns[["log2FC"]]
          ]]
        )
      ),
      pvalue = suppressWarnings(
        as.numeric(
          .data[[
            comparison_columns[["pvalue"]]
          ]]
        )
      ),
      Comparison = comparison_name,
      Compartment = compartment_name
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      is.finite(log2FC),
      is.finite(pvalue),
      pvalue >= 0,
      pvalue <= 1
    )
}

extract_protein_set_S5D <- function(
    result_table,
    direction
) {
  
  significant_table <- result_table %>%
    dplyr::filter(
      pvalue <
        protein_p_threshold_S5D
    )
  
  if (identical(direction, "promoted")) {
    
    selected <- significant_table %>%
      dplyr::filter(
        log2FC < 0
      )
    
  } else if (
    identical(direction, "suppressed")
  ) {
    
    selected <- significant_table %>%
      dplyr::filter(
        log2FC > 0
      )
    
  } else {
    
    stop(
      "direction must be promoted or suppressed."
    )
  }
  
  sort(
    unique(
      selected$gene
    )
  )
}

protein_nmd3_nucleus_S5D <-
  extract_protein_result_S5D(
    "NMD3_Nucleus",
    "Nucleus"
  )

protein_nmd3_cytoplasm_S5D <-
  extract_protein_result_S5D(
    "NMD3_Cytoplasm",
    "Cytoplasm"
  )

protein_xpo1_nucleus_S5D <-
  extract_protein_result_S5D(
    "XPO1_Nucleus",
    "Nucleus"
  )

protein_xpo1_cytoplasm_S5D <-
  extract_protein_result_S5D(
    "XPO1_Cytoplasm",
    "Cytoplasm"
  )

Protein_NMD3_promoted_Nucleus_S5D <-
  extract_protein_set_S5D(
    protein_nmd3_nucleus_S5D,
    "promoted"
  )

Protein_NMD3_promoted_Cytoplasm_S5D <-
  extract_protein_set_S5D(
    protein_nmd3_cytoplasm_S5D,
    "promoted"
  )

Protein_NMD3_suppressed_Nucleus_S5D <-
  extract_protein_set_S5D(
    protein_nmd3_nucleus_S5D,
    "suppressed"
  )

Protein_NMD3_suppressed_Cytoplasm_S5D <-
  extract_protein_set_S5D(
    protein_nmd3_cytoplasm_S5D,
    "suppressed"
  )

Protein_XPO1_promoted_Nucleus_S5D <-
  extract_protein_set_S5D(
    protein_xpo1_nucleus_S5D,
    "promoted"
  )

Protein_XPO1_promoted_Cytoplasm_S5D <-
  extract_protein_set_S5D(
    protein_xpo1_cytoplasm_S5D,
    "promoted"
  )

Protein_XPO1_suppressed_Nucleus_S5D <-
  extract_protein_set_S5D(
    protein_xpo1_nucleus_S5D,
    "suppressed"
  )

Protein_XPO1_suppressed_Cytoplasm_S5D <-
  extract_protein_set_S5D(
    protein_xpo1_cytoplasm_S5D,
    "suppressed"
  )


# ============================================================
# 4. Store the 16 gene sets
# ============================================================

sets_S5D <- list(
  
  RNA_NMD3_promoted_Nucleus =
    RNA_NMD3_promoted_Nucleus_S5D,
  
  RNA_NMD3_promoted_Cytoplasm =
    RNA_NMD3_promoted_Cytoplasm_S5D,
  
  RNA_NMD3_suppressed_Nucleus =
    RNA_NMD3_suppressed_Nucleus_S5D,
  
  RNA_NMD3_suppressed_Cytoplasm =
    RNA_NMD3_suppressed_Cytoplasm_S5D,
  
  RNA_XPO1_promoted_Nucleus =
    RNA_XPO1_promoted_Nucleus_S5D,
  
  RNA_XPO1_promoted_Cytoplasm =
    RNA_XPO1_promoted_Cytoplasm_S5D,
  
  RNA_XPO1_suppressed_Nucleus =
    RNA_XPO1_suppressed_Nucleus_S5D,
  
  RNA_XPO1_suppressed_Cytoplasm =
    RNA_XPO1_suppressed_Cytoplasm_S5D,
  
  Protein_NMD3_promoted_Nucleus =
    Protein_NMD3_promoted_Nucleus_S5D,
  
  Protein_NMD3_promoted_Cytoplasm =
    Protein_NMD3_promoted_Cytoplasm_S5D,
  
  Protein_NMD3_suppressed_Nucleus =
    Protein_NMD3_suppressed_Nucleus_S5D,
  
  Protein_NMD3_suppressed_Cytoplasm =
    Protein_NMD3_suppressed_Cytoplasm_S5D,
  
  Protein_XPO1_promoted_Nucleus =
    Protein_XPO1_promoted_Nucleus_S5D,
  
  Protein_XPO1_promoted_Cytoplasm =
    Protein_XPO1_promoted_Cytoplasm_S5D,
  
  Protein_XPO1_suppressed_Nucleus =
    Protein_XPO1_suppressed_Nucleus_S5D,
  
  Protein_XPO1_suppressed_Cytoplasm =
    Protein_XPO1_suppressed_Cytoplasm_S5D
)

stopifnot(
  length(sets_S5D) == 16
)


# ============================================================
# 5. Two-set proportional Venn plotting function
# ============================================================

make_venn_S5D <- function(
    nucleus_set,
    cytoplasm_set,
    panel_title,
    panel_subtitle
) {
  
  nucleus_set <- unique(
    as.character(nucleus_set)
  )
  
  cytoplasm_set <- unique(
    as.character(cytoplasm_set)
  )
  
  overlap <- length(
    intersect(
      nucleus_set,
      cytoplasm_set
    )
  )
  
  regulator_label <- sub(
    "^(RNA-seq|Proteomics):[[:space:]]*",
    "",
    panel_title
  )
  
  venn_grobs <- VennDiagram::draw.pairwise.venn(
    area1 = length(nucleus_set),
    area2 = length(cytoplasm_set),
    cross.area = overlap,
    category = c(
      paste0(
        "Nucleus\n",
        regulator_label
      ),
      paste0(
        "Cytoplasm\n",
        regulator_label
      )
    ),
    scaled = TRUE,
    inverted = FALSE,
    ext.text = TRUE,
    rotation.degree = 0,
    col = c(
      "#0072C6",
      "#D00000"
    ),
    fill = c(
      "transparent",
      "transparent"
    ),
    alpha = c(
      0,
      0
    ),
    lwd = c(
      2.2,
      2.2
    ),
    lty = c(
      "solid",
      "solid"
    ),
    cex = 1.35,
    fontface = "bold",
    cat.col = c(
      "#0072C6",
      "#D00000"
    ),
    cat.cex = 1.05,
    cat.fontface = "bold",
    cat.pos = c(
      -20,
      20
    ),
    cat.dist = c(
      0.055,
      0.055
    ),
    margin = 0.08,
    ind = FALSE
  )
  
  venn_tree <- grid::gTree(
    children = do.call(
      grid::gList,
      venn_grobs
    )
  )
  
  patchwork::wrap_elements(
    full = venn_tree
  ) +
    patchwork::plot_annotation(
      subtitle = panel_subtitle,
      theme = ggplot2::theme(
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = 9
        )
      )
    )
}


}))

# =============================================================================
# Fig3B
# =============================================================================
nmd3_step("Fig3B", c("Prepare_Fig3B_S6"), quote({

# Same sixteen thresholded gene sets as FigS6D; each gene is counted once.
set_union <- function(prefix) unique(unlist(sets_S5D[startsWith(names(sets_S5D), prefix)], use.names=FALSE))
summary_3B <- expand.grid(Assay=c("Protein","RNA"), Treatment=c("XPO1","NMD3"), stringsAsFactors=FALSE)
summary_3B$N <- mapply(function(a,t) length(set_union(paste0(a,"_",t,"_"))), summary_3B$Assay, summary_3B$Treatment)
summary_3B$KPT_reference <- summary_3B$N[match(paste(summary_3B$Assay,"XPO1"), paste(summary_3B$Assay,summary_3B$Treatment))]
stopifnot(all(summary_3B$KPT_reference>0))
summary_3B$Percent_of_KPT <- 100 * summary_3B$N / summary_3B$KPT_reference
summary_3B$Assay <- factor(summary_3B$Assay, levels=c("Protein","RNA"))
summary_3B$Treatment <- factor(summary_3B$Treatment,levels=c("XPO1","NMD3"),labels=c("KPT-330","siNMD3"))
p_3B <- ggplot2::ggplot(summary_3B,ggplot2::aes(Assay,Percent_of_KPT,fill=Treatment)) +
 ggplot2::geom_col(position=ggplot2::position_dodge(.75),width=.65) +
 ggplot2::geom_text(ggplot2::aes(label=sprintf("n = %d\n(%.1f%%)",N,Percent_of_KPT)),
 position=ggplot2::position_dodge(.75),vjust=-.25,size=3.8) +
 ggplot2::scale_fill_manual(values=c("KPT-330"="#D73027","siNMD3"="#4575B4")) +
 ggplot2::scale_y_continuous(limits=c(0,120),breaks=c(0,25,50,75,100)) +
 ggplot2::labs(x=NULL,y="Perturbed genes (% of KPT-330 count)",fill=NULL) +
 ggplot2::theme_classic(base_size=13) + ggplot2::theme(legend.position="bottom")
nmd3_save_pdf("Fig3B_Perturbed_Protein_RNA_Summary.pdf",p_3B,width=6,height=5)
utils::write.csv(summary_3B,"Fig3B_Counts_Percentages.csv",row.names=FALSE)

}))

# =============================================================================
# Fig3C
# =============================================================================
nmd3_step("Fig3C", c("Fig3A"), quote({
# Fig3C — analysis below
# Central 95% RNA-protein changes
# KPT-330 vs DMSO and siNMD3 vs siCTRL
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# 1. Calculate central 95% boundaries
# ------------------------------------------------------------

get_95_bounds <- function(df) {
  
  c(
    xmin = unname(
      quantile(
        df$log2FC_RNA,
        0.025,
        na.rm = TRUE
      )
    ),
    
    xmax = unname(
      quantile(
        df$log2FC_RNA,
        0.975,
        na.rm = TRUE
      )
    ),
    
    ymin = unname(
      quantile(
        df$log2FC_Prot,
        0.025,
        na.rm = TRUE
      )
    ),
    
    ymax = unname(
      quantile(
        df$log2FC_Prot,
        0.975,
        na.rm = TRUE
      )
    )
  )
}

kpt_95 <- get_95_bounds(multi_kpt)
si_95  <- get_95_bounds(multi_si)

cat("\nKPT-330 vs DMSO central 95%:\n")
print(round(kpt_95, 2))

cat("\nsiNMD3 vs siCTRL central 95%:\n")
print(round(si_95, 2))

# ------------------------------------------------------------
# 2. Plotting function
# ------------------------------------------------------------

plot_center95 <- function(
    df,
    bounds,
    box_color,
    x_label,
    y_label,
    x_lim,
    y_lim = c(-3, 3)
) {
  
  ggplot(
    df,
    aes(
      x = log2FC_RNA,
      y = log2FC_Prot
    )
  ) +
    
    # Points
    geom_point(
      size = 1.6,
      alpha = 0.55,
      color = "grey20"
    ) +
    
    # Zero lines
    geom_vline(
      xintercept = 0,
      color = "black",
      linewidth = 0.6
    ) +
    
    geom_hline(
      yintercept = 0,
      color = "black",
      linewidth = 0.6
    ) +
    
    # 95% rectangle
    annotate(
      "rect",
      xmin = bounds["xmin"],
      xmax = bounds["xmax"],
      ymin = bounds["ymin"],
      ymax = bounds["ymax"],
      fill = NA,
      color = box_color,
      linewidth = 1
    ) +
    
    # Vertical lines to bottom X axis
    annotate(
      "segment",
      x = bounds["xmin"],
      xend = bounds["xmin"],
      y = y_lim[1],
      yend = bounds["ymin"],
      color = box_color,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    annotate(
      "segment",
      x = bounds["xmax"],
      xend = bounds["xmax"],
      y = y_lim[1],
      yend = bounds["ymin"],
      color = box_color,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    # Horizontal lines to left Y axis
    annotate(
      "segment",
      x = x_lim[1],
      xend = bounds["xmin"],
      y = bounds["ymin"],
      yend = bounds["ymin"],
      color = box_color,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    annotate(
      "segment",
      x = x_lim[1],
      xend = bounds["xmin"],
      y = bounds["ymax"],
      yend = bounds["ymax"],
      color = box_color,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    # X-axis boundary labels
    annotate(
      "text",
      x = bounds["xmin"],
      y = y_lim[1],
      label = sprintf("%.2f", bounds["xmin"]),
      color = box_color,
      size = 5,
      angle = 45,
      hjust = 1,
      vjust = 1.2
    ) +
    
    annotate(
      "text",
      x = bounds["xmax"],
      y = y_lim[1],
      label = sprintf("%.2f", bounds["xmax"]),
      color = box_color,
      size = 5,
      angle = 45,
      hjust = 1,
      vjust = 1.2
    ) +
    
    # Y-axis boundary labels
    annotate(
      "text",
      x = x_lim[1],
      y = bounds["ymin"],
      label = sprintf("%.2f", bounds["ymin"]),
      color = box_color,
      size = 5,
      hjust = 1.15
    ) +
    
    annotate(
      "text",
      x = x_lim[1],
      y = bounds["ymax"],
      label = sprintf("%.2f", bounds["ymax"]),
      color = box_color,
      size = 5,
      hjust = 1.15
    ) +
    
    labs(
      x = x_label,
      y = y_label
    ) +
    
    coord_cartesian(
      xlim = x_lim,
      ylim = y_lim,
      clip = "off",
      expand = FALSE
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      axis.line = element_blank(),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      axis.text = element_text(
        color = "black"
      ),
      
      plot.margin = margin(
        t = 10,
        r = 10,
        b = 35,
        l = 45
      )
    )
}
# ------------------------------------------------------------
# 3. KPT-330 vs DMSO
# ------------------------------------------------------------

# KPT-330

p_4C_kpt <- plot_center95(
  df = multi_kpt,
  bounds = kpt_95,
  box_color = "red",
  x_label = "RNA log2FC (KPT-330 vs DMSO)",
  y_label = "Protein log2FC (KPT-330 vs DMSO)",
  x_lim = c(-3, 3),
  y_lim = c(-3, 3)
)

# siNMD3

p_4C_si <- plot_center95(
  df = multi_si,
  bounds = si_95,
  box_color = "#2448E8",
  x_label = "RNA log2FC (siNMD3 vs siCTRL)",
  y_label = "Protein log2FC (siNMD3 vs siCTRL)",
  x_lim = c(-1.8, 1.8),
  y_lim = c(-3, 3)
)

p_4C <- p_4C_kpt + p_4C_si


# New PPT Fig.3C: write the already calculated top row once.
nmd3_save_pdf(
  "Fig3C_upper_Central95_RNA_Protein.pdf",
  plot = p_4C, width = 10, height = 5,
  device = grDevices::cairo_pdf
)








# ============================================================
# Fig3C — analysis below
# RNA-protein concordance
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# 1. Spearman correlations
# ------------------------------------------------------------

cor_kpt <- suppressWarnings(
  cor.test(
    multi_kpt$log2FC_RNA,
    multi_kpt$log2FC_Prot,
    method = "spearman",
    exact = FALSE
  )
)

cor_si <- suppressWarnings(
  cor.test(
    multi_si$log2FC_RNA,
    multi_si$log2FC_Prot,
    method = "spearman",
    exact = FALSE
  )
)

label_kpt <- sprintf(
  "rho = %.3f, p = %.2e\nn = %d",
  unname(cor_kpt$estimate),
  cor_kpt$p.value,
  nrow(multi_kpt)
)

label_si <- sprintf(
  "rho = %.3f, p = %.2e\nn = %d",
  unname(cor_si$estimate),
  cor_si$p.value,
  nrow(multi_si)
)

cat(label_kpt, "\n")
cat(label_si, "\n")

# ------------------------------------------------------------
# 2. KPT-330 vs DMSO
# ------------------------------------------------------------

p_4C_cor_kpt <- ggplot(
  multi_kpt,
  aes(
    x = log2FC_RNA,
    y = log2FC_Prot
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.5,
    color = "grey30"
  ) +
  geom_smooth(
    method = "lm",
    color = "red",
    fill = "grey70",
    linewidth = 0.8,
    se = TRUE,
    alpha = 0.3
  ) +
  annotate(
    "text",
    x = min(multi_kpt$log2FC_RNA) +
      0.05 * diff(range(multi_kpt$log2FC_RNA)),
    y = min(multi_kpt$log2FC_Prot) +
      0.08 * diff(range(multi_kpt$log2FC_Prot)),
    label = label_kpt,
    hjust = 0,
    vjust = 0,
    size = 5
  ) +
  labs(
    x = "RNA log2FC (KPT-330 vs DMSO)",
    y = "Protein log2FC (KPT-330 vs DMSO)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  )

# ------------------------------------------------------------
# 3. siNMD3 vs siCTRL
# ------------------------------------------------------------

p_4C_cor_si <- ggplot(
  multi_si,
  aes(
    x = log2FC_RNA,
    y = log2FC_Prot
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.5,
    color = "grey30"
  ) +
  geom_smooth(
    method = "lm",
    color = "#2456FF",
    fill = "grey70",
    linewidth = 0.8,
    se = TRUE,
    alpha = 0.3
  ) +
  annotate(
    "text",
    x = min(multi_si$log2FC_RNA) +
      0.05 * diff(range(multi_si$log2FC_RNA)),
    y = max(multi_si$log2FC_Prot) -
      0.05 * diff(range(multi_si$log2FC_Prot)),
    label = label_si,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  labs(
    x = "RNA log2FC (siNMD3 vs siCTRL)",
    y = "Protein log2FC (siNMD3 vs siCTRL)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  )

# ------------------------------------------------------------
# 4. Combine
# ------------------------------------------------------------

p_4C_lower <- p_4C_cor_kpt + p_4C_cor_si


nmd3_save_pdf(
  "Fig3C_lower_RNA_Protein_Spearman.pdf",
  p_4C_lower,
  width = 10,
  height = 5,
  device = cairo_pdf
)



}))

# =============================================================================
# Fig3D
# =============================================================================
nmd3_step("Fig3D", character(), quote({
# Fig3D — analysis below
# RNA-protein Concordance At the Top (CAT)
#
# Top:    KPT-330 vs DMSO, cytoplasm
# Bottom: siNMD3 vs siCTRL, cytoplasm
#
# This script is standalone. It rebuilds both RNA and protein
# comparisons from the original count matrix and proteomics report.
# It preserves the submitted analysis:
#   RNA:    DESeq2 with apeglm shrinkage
#   Protein: limma, 3 controls vs 3 treated samples
#   Ranking: abs(log2FC) * -log10(BH-adjusted P)
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ashr)
  library(readxl)
  library(limma)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

# ------------------------------------------------------------
# 0. Shared input objects and CAT parameters
# ------------------------------------------------------------

kstar_frac_S5E <- 0.20
max_k_frac_S5E <- 0.50
tiny_S5E <- 1e-300

# Reuse the in-house RNA and proteomics tables already loaded upstream.
stopifnot(
  exists("expr_matrix", inherits = TRUE),
  exists("prot_raw", inherits = TRUE),
  "Gene" %in% colnames(prot_raw)
)

# ------------------------------------------------------------
# 1. Reuse original data
# ------------------------------------------------------------

expr_matrix_S5E <- expr_matrix

protein_raw_S5E <- as.data.frame(
  prot_raw,
  check.names = FALSE
)

stopifnot(
  "Gene" %in% colnames(protein_raw_S5E),
  !anyDuplicated(rownames(expr_matrix_S5E))
)

# ------------------------------------------------------------
# 2. Exact cytoplasmic RNA sample definitions
# ------------------------------------------------------------

xpo1_ctrl_rna_S5E <- c(
  "D1C", "D2C", "D3C"
)

xpo1_treat_rna_S5E <- c(
  "X1C", "X2C", "X3C"
)

nmd3_ctrl_rna_S5E <- c(
  "C1C", "C2C", "C3C"
)

nmd3_treat_rna_S5E <- c(
  "N1C", "N2C", "N3C"
)

all_rna_samples_S5E <- c(
  xpo1_ctrl_rna_S5E,
  xpo1_treat_rna_S5E,
  nmd3_ctrl_rna_S5E,
  nmd3_treat_rna_S5E
)

stopifnot(
  all(
    all_rna_samples_S5E %in%
      colnames(expr_matrix_S5E)
  )
)

# ------------------------------------------------------------
# 3. DESeq2 helper
#
# Positive log2FC is treated minus control.
# ------------------------------------------------------------

shrink_lfc_S5E <- function(
    dds,
    coefficient_name
) {
  tryCatch(
    DESeq2::lfcShrink(
      dds,
      coef = coefficient_name,
      type = "apeglm"
    ),
    error = function(e) {
      message(
        "apeglm failed; using ashr: ",
        conditionMessage(e)
      )
      
      DESeq2::lfcShrink(
        dds,
        coef = coefficient_name,
        type = "ashr"
      )
    }
  )
}

run_deseq2_cytoplasm_S5E <- function(
    count_data,
    control_samples,
    treated_samples,
    control_label,
    treated_label
) {
  stopifnot(
    length(control_samples) == 3,
    length(treated_samples) == 3,
    all(
      c(
        control_samples,
        treated_samples
      ) %in% colnames(count_data)
    )
  )
  
  samples <- c(
    control_samples,
    treated_samples
  )
  
  counts <- as.matrix(
    count_data[
      ,
      samples,
      drop = FALSE
    ]
  )
  
  storage.mode(counts) <- "integer"
  
  stopifnot(
    identical(
      colnames(counts),
      samples
    ),
    all(is.finite(counts)),
    all(counts >= 0)
  )
  
  sample_group <- factor(
    c(
      rep(control_label, 3),
      rep(treated_label, 3)
    ),
    levels = c(
      control_label,
      treated_label
    )
  )
  
  coldata <- data.frame(
    condition = sample_group,
    row.names = samples
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )
  
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )
  
  coefficient_names <- setdiff(
    DESeq2::resultsNames(dds),
    "Intercept"
  )
  
  stopifnot(
    length(coefficient_names) == 1
  )
  
  coefficient_name <- coefficient_names[1]
  
  result_shrunken <- shrink_lfc_S5E(
    dds,
    coefficient_name
  )
  
  result_raw <- DESeq2::results(
    dds,
    name = coefficient_name
  )
  
  shrunken_df <- as.data.frame(
    result_shrunken
  ) %>%
    tibble::rownames_to_column(
      "gene"
    )
  
  raw_df <- as.data.frame(
    result_raw
  ) %>%
    tibble::rownames_to_column(
      "gene"
    )
  
  rna_result <- shrunken_df %>%
    dplyr::transmute(
      gene = toupper(
        trimws(
          as.character(gene)
        )
      ),
      log2FC_RNA = dplyr::coalesce(
        suppressWarnings(
          as.numeric(log2FoldChange)
        ),
        0
      )
    ) %>%
    dplyr::left_join(
      raw_df %>%
        dplyr::transmute(
          gene = toupper(
            trimws(
              as.character(gene)
            )
          ),
          pval_RNA = dplyr::coalesce(
            suppressWarnings(
              as.numeric(pvalue)
            ),
            1
          ),
          padj_RNA = dplyr::coalesce(
            suppressWarnings(
              as.numeric(padj)
            ),
            1
          )
        ),
      by = "gene"
    ) %>%
    dplyr::filter(
      gene != ""
    ) %>%
    dplyr::distinct(
      gene,
      .keep_all = TRUE
    )
  
  rna_result
}

# ------------------------------------------------------------
# 4. Run both RNA comparisons
# ------------------------------------------------------------

rna_kpt_S5E <- run_deseq2_cytoplasm_S5E(
  count_data = expr_matrix_S5E,
  control_samples = xpo1_ctrl_rna_S5E,
  treated_samples = xpo1_treat_rna_S5E,
  control_label = "DMSO",
  treated_label = "KPT330"
)

rna_si_S5E <- run_deseq2_cytoplasm_S5E(
  count_data = expr_matrix_S5E,
  control_samples = nmd3_ctrl_rna_S5E,
  treated_samples = nmd3_treat_rna_S5E,
  control_label = "siCTRL",
  treated_label = "siNMD3"
)

# ------------------------------------------------------------
# 5. Exact cytoplasmic proteomics sample definitions
# ------------------------------------------------------------

xpo1_ctrl_protein_S5E <- c(
  "Norm_DMSO_Cyto_D1C",
  "Norm_DMSO_Cyto_D2C",
  "Norm_DMSO_Cyto_D3C"
)

xpo1_treat_protein_S5E <- c(
  "Norm_XPO1inhibitor_Cyto_X1C",
  "Norm_XPO1inhibitor_Cyto_X2C",
  "Norm_XPO1inhibitor_Cyto_X3C"
)

nmd3_ctrl_protein_S5E <- grep(
  "^Norm_siRNAxCtrl_Cyto_",
  colnames(protein_raw_S5E),
  value = TRUE
)

nmd3_treat_protein_S5E <- grep(
  "^Norm_siRNAxNMD3xKD_Cyto_",
  colnames(protein_raw_S5E),
  value = TRUE
)

nmd3_ctrl_protein_S5E <- sort(
  nmd3_ctrl_protein_S5E
)

nmd3_treat_protein_S5E <- sort(
  nmd3_treat_protein_S5E
)

stopifnot(
  all(
    c(
      xpo1_ctrl_protein_S5E,
      xpo1_treat_protein_S5E
    ) %in% colnames(protein_raw_S5E)
  ),
  length(nmd3_ctrl_protein_S5E) == 3,
  length(nmd3_treat_protein_S5E) == 3
)

cat("KPT-330 protein control samples:\n")
print(xpo1_ctrl_protein_S5E)

cat("KPT-330 protein treated samples:\n")
print(xpo1_treat_protein_S5E)

cat("siNMD3 protein control samples:\n")
print(nmd3_ctrl_protein_S5E)

cat("siNMD3 protein treated samples:\n")
print(nmd3_treat_protein_S5E)

# ------------------------------------------------------------
# 6. limma 3-vs-3 helper
#
# Positive log2FC is treated minus control.
# ------------------------------------------------------------

run_limma_3v3_S5E <- function(
    protein_data,
    control_columns,
    treated_columns
) {
  stopifnot(
    length(control_columns) == 3,
    length(treated_columns) == 3
  )
  
  expression_columns <- c(
    control_columns,
    treated_columns
  )
  
  protein_table <- protein_data %>%
    dplyr::filter(
      !is.na(Gene),
      trimws(as.character(Gene)) != ""
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(expression_columns),
        ~ suppressWarnings(
          as.numeric(.)
        )
      )
    ) %>%
    dplyr::mutate(
      Gene = toupper(
        trimws(
          as.character(Gene)
        )
      )
    ) %>%
    dplyr::group_by(
      Gene
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(expression_columns),
        ~ if (all(is.na(.))) {
          NA_real_
        } else {
          stats::median(
            .,
            na.rm = TRUE
          )
        }
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      rowSums(
        is.na(
          dplyr::across(
            dplyr::all_of(expression_columns)
          )
        )
      ) == 0
    )
  
  expression_matrix <- as.matrix(
    protein_table[
      ,
      expression_columns,
      drop = FALSE
    ]
  )
  
  rownames(expression_matrix) <-
    protein_table$Gene
  
  storage.mode(expression_matrix) <-
    "numeric"
  
  stopifnot(
    all(is.finite(expression_matrix))
  )
  
  if (
    stats::median(
      expression_matrix,
      na.rm = TRUE
    ) > 50
  ) {
    expression_matrix <- log2(
      expression_matrix + 1
    )
  }
  
  treatment_group <- factor(
    c(
      rep("CTRL", 3),
      rep("TREAT", 3)
    ),
    levels = c(
      "CTRL",
      "TREAT"
    )
  )
  
  design <- stats::model.matrix(
    ~ treatment_group
  )
  
  fit <- limma::lmFit(
    expression_matrix,
    design
  )
  
  contrast_matrix <- limma::makeContrasts(
    TreatVsCtrl = treatment_groupTREAT,
    levels = design
  )
  
  fit_contrast <- limma::contrasts.fit(
    fit,
    contrast_matrix
  )
  
  fit_contrast <- limma::eBayes(
    fit_contrast,
    robust = TRUE,
    trend = TRUE
  )
  
  limma::topTable(
    fit_contrast,
    coef = "TreatVsCtrl",
    number = Inf,
    sort.by = "none"
  ) %>%
    tibble::rownames_to_column(
      "gene"
    ) %>%
    dplyr::transmute(
      gene = toupper(
        trimws(
          as.character(gene)
        )
      ),
      log2FC_Prot = suppressWarnings(
        as.numeric(logFC)
      ),
      pval_Prot = suppressWarnings(
        as.numeric(P.Value)
      ),
      padj_Prot = stats::p.adjust(
        pval_Prot,
        method = "BH"
      )
    ) %>%
    dplyr::distinct(
      gene,
      .keep_all = TRUE
    )
}

# ------------------------------------------------------------
# 7. Run both protein comparisons
# ------------------------------------------------------------

prot_kpt_S5E <- run_limma_3v3_S5E(
  protein_data = protein_raw_S5E,
  control_columns = xpo1_ctrl_protein_S5E,
  treated_columns = xpo1_treat_protein_S5E
)

prot_si_S5E <- run_limma_3v3_S5E(
  protein_data = protein_raw_S5E,
  control_columns = nmd3_ctrl_protein_S5E,
  treated_columns = nmd3_treat_protein_S5E
)

# ------------------------------------------------------------
# 8. Match RNA and protein independently for each comparison
#
# Do not restrict either comparison to a shared gene universe.
# ------------------------------------------------------------

make_multiomics_input_S5E <- function(
    rna_result,
    protein_result
) {
  rna_result %>%
    dplyr::inner_join(
      protein_result,
      by = "gene"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(
          pval_RNA,
          padj_RNA,
          pval_Prot,
          padj_Prot
        ),
        ~ dplyr::coalesce(
          as.numeric(.),
          1
        )
      )
    ) %>%
    dplyr::filter(
      is.finite(log2FC_RNA),
      is.finite(log2FC_Prot)
    ) %>%
    dplyr::distinct(
      gene,
      .keep_all = TRUE
    )
}

multi_kpt_S5E <- make_multiomics_input_S5E(
  rna_kpt_S5E,
  prot_kpt_S5E
)

multi_si_S5E <- make_multiomics_input_S5E(
  rna_si_S5E,
  prot_si_S5E
)

cat(
  "\nKPT-330 matched RNA-protein genes:",
  nrow(multi_kpt_S5E),
  "\n"
)

cat(
  "siNMD3 matched RNA-protein genes:",
  nrow(multi_si_S5E),
  "\n"
)

# ------------------------------------------------------------
# 9. CAT calculation
# ------------------------------------------------------------

compute_CAT_S5E <- function(
    matched_data,
    kstar_fraction = 0.20,
    max_k_fraction = 0.50
) {
  clean_data <- matched_data %>%
    dplyr::transmute(
      gene = toupper(
        trimws(
          as.character(gene)
        )
      ),
      log2FC_RNA = as.numeric(log2FC_RNA),
      padj_RNA = as.numeric(padj_RNA),
      log2FC_Prot = as.numeric(log2FC_Prot),
      padj_Prot = as.numeric(padj_Prot)
    ) %>%
    dplyr::filter(
      gene != "",
      is.finite(log2FC_RNA),
      is.finite(log2FC_Prot)
    ) %>%
    dplyr::mutate(
      padj_RNA = dplyr::if_else(
        is.finite(padj_RNA) &
          padj_RNA >= 0 &
          padj_RNA <= 1,
        pmax(
          padj_RNA,
          tiny_S5E
        ),
        1
      ),
      padj_Prot = dplyr::if_else(
        is.finite(padj_Prot) &
          padj_Prot >= 0 &
          padj_Prot <= 1,
        pmax(
          padj_Prot,
          tiny_S5E
        ),
        1
      ),
      RNA_rank_score =
        abs(log2FC_RNA) *
        -log10(padj_RNA),
      Protein_rank_score =
        abs(log2FC_Prot) *
        -log10(padj_Prot)
    ) %>%
    dplyr::group_by(
      gene
    ) %>%
    dplyr::slice_max(
      order_by =
        RNA_rank_score +
        Protein_rank_score,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
  
  N <- nrow(clean_data)
  
  stopifnot(N >= 60)
  
  ranked_rna <- clean_data %>%
    dplyr::arrange(
      dplyr::desc(RNA_rank_score),
      gene
    ) %>%
    dplyr::pull(
      gene
    )
  
  ranked_protein <- clean_data %>%
    dplyr::arrange(
      dplyr::desc(Protein_rank_score),
      gene
    ) %>%
    dplyr::pull(
      gene
    )
  
  k_step <- max(
    25,
    floor(N / 60)
  )
  
  max_k <- min(
    floor(
      N * max_k_fraction
    ),
    5000
  )
  
  K_values <- seq(
    from = k_step,
    to = max_k,
    by = k_step
  )
  
  cat_curve <- purrr::map_dfr(
    K_values,
    function(k) {
      top_rna <- ranked_rna[
        seq_len(k)
      ]
      
      top_protein <- ranked_protein[
        seq_len(k)
      ]
      
      overlap_n <- length(
        intersect(
          top_rna,
          top_protein
        )
      )
      
      expected_overlap <- k^2 / N
      
      overlap_variance <-
        k *
        (k / N) *
        (1 - k / N) *
        ((N - k) / (N - 1))
      
      se_proportion <-
        sqrt(overlap_variance) / k
      
      tibble::tibble(
        K = k,
        overlap_n = overlap_n,
        observed = overlap_n / k,
        expected = k / N,
        null_lower = pmax(
          0,
          expected_overlap / k -
            1.96 * se_proportion
        ),
        null_upper = pmin(
          1,
          expected_overlap / k +
            1.96 * se_proportion
        )
      )
    }
  )
  
  delta_auc <- mean(
    cat_curve$observed -
      cat_curve$expected
  )
  
  target_k <- floor(
    N * kstar_fraction
  )
  
  Kstar <- cat_curve$K[
    which.min(
      abs(
        cat_curve$K -
          target_k
      )
    )
  ]
  
  overlap_star <- cat_curve$overlap_n[
    cat_curve$K == Kstar
  ]
  
  p_star <- stats::phyper(
    q = overlap_star - 1,
    m = Kstar,
    n = N - Kstar,
    k = Kstar,
    lower.tail = FALSE
  )
  
  statistics <- tibble::tibble(
    N = N,
    AUC = delta_auc,
    Kstar = Kstar,
    overlap_at_Kstar = overlap_star,
    p_at_Kstar = p_star
  )
  
  list(
    curve = cat_curve,
    statistics = statistics,
    ranked_data = clean_data
  )
}

# ------------------------------------------------------------
# 10. Run CAT for both comparisons
# ------------------------------------------------------------

cat_xpo1_S5E <- compute_CAT_S5E(
  multi_kpt_S5E,
  kstar_fraction = kstar_frac_S5E,
  max_k_fraction = max_k_frac_S5E
)

cat_nmd3_S5E <- compute_CAT_S5E(
  multi_si_S5E,
  kstar_fraction = kstar_frac_S5E,
  max_k_fraction = max_k_frac_S5E
)

cat("\nXPO1 CAT statistics:\n")
print(cat_xpo1_S5E$statistics)

cat("\nNMD3 CAT statistics:\n")
print(cat_nmd3_S5E$statistics)

# ------------------------------------------------------------
# 11. CAT plot helper
# ------------------------------------------------------------

make_CAT_plot_S5E <- function(
    cat_result,
    comparison_title,
    show_x_title = TRUE
) {
  panel_statistics <- cat_result$statistics
  cat_curve <- cat_result$curve
  
  p_value_text <- if (
    is.finite(
      panel_statistics$p_at_Kstar[1]
    ) &&
    panel_statistics$p_at_Kstar[1] > 0
  ) {
    formatC(
      panel_statistics$p_at_Kstar[1],
      format = "e",
      digits = 2
    )
  } else {
    "< machine precision"
  }
  
  panel_title <- paste0(
    "CAT - (",
    comparison_title,
    ")\n",
    "All genes (N = ",
    panel_statistics$N[1],
    ")  AUC = ",
    sprintf(
      "%.3f",
      panel_statistics$AUC[1]
    ),
    "\nK* = ",
    panel_statistics$Kstar[1],
    "   p = ",
    p_value_text
  )
  
  ggplot2::ggplot(
    cat_curve,
    ggplot2::aes(
      x = K,
      y = observed
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = null_lower,
        ymax = null_upper
      ),
      fill = "grey85",
      color = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = expected
      ),
      linetype = "dashed",
      color = "black",
      linewidth = 0.6
    ) +
    ggplot2::geom_line(
      color = "black",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      color = "black",
      size = 1.6
    ) +
    ggplot2::labs(
      title = panel_title,
      x = if (show_x_title) {
        "Top-K ranked genes"
      } else {
        NULL
      },
      y = paste0(
        "Concordance\n",
        "(Proportion of overlapping genes)"
      )
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::comma
    ) +
    ggplot2::scale_y_continuous(
      limits = c(
        0,
        NA
      ),
      expand = ggplot2::expansion(
        mult = c(
          0.02,
          0.08
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 12
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 12,
        lineheight = 1.05
      ),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        size = 12
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        size = 11
      ),
      axis.text = ggplot2::element_text(
        color = "black",
        size = 10
      ),
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      ),
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      ),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(
        8,
        8,
        8,
        8
      )
    )
}

# ------------------------------------------------------------
# 12. Draw Figure S5E
# ------------------------------------------------------------

p_CAT_xpo1_S5E <- make_CAT_plot_S5E(
  cat_xpo1_S5E,
  comparison_title = "KPT-330 vs Ctrl",
  show_x_title = FALSE
)

p_CAT_nmd3_S5E <- make_CAT_plot_S5E(
  cat_nmd3_S5E,
  comparison_title = "siNMD3 vs siCtrl",
  show_x_title = TRUE
)

p_S5E <- p_CAT_xpo1_S5E /
  p_CAT_nmd3_S5E +
  patchwork::plot_layout(
    heights = c(
      1,
      1
    )
  )


# ------------------------------------------------------------
# 13. Save Figure S5E
# ------------------------------------------------------------

nmd3_save_pdf(
  filename = "Fig3D_CAT_RNA_Protein.pdf",
  plot = p_S5E,
  width = 6,
  height = 10,
  units = "in",
  device = grDevices::cairo_pdf
)

# Optional PNG export intentionally omitted: the PPT workflow requires PDFs only.

# ------------------------------------------------------------
# 14. Export source data and calculated results
# ------------------------------------------------------------

CAT_statistics_S5E <- dplyr::bind_rows(
  cat_xpo1_S5E$statistics %>%
    dplyr::mutate(
      Comparison = "KPT-330 vs DMSO"
    ),
  cat_nmd3_S5E$statistics %>%
    dplyr::mutate(
      Comparison = "siNMD3 vs siCTRL"
    )
) %>%
  dplyr::select(
    Comparison,
    dplyr::everything()
  )

print(CAT_statistics_S5E)

readr::write_csv(
  CAT_statistics_S5E,
  "FigureS5E_CAT_statistics.csv"
)

readr::write_csv(
  cat_xpo1_S5E$curve %>%
    dplyr::mutate(
      Comparison = "KPT-330 vs DMSO"
    ),
  "FigureS5E_CAT_XPO1_curve.csv"
)

readr::write_csv(
  cat_nmd3_S5E$curve %>%
    dplyr::mutate(
      Comparison = "siNMD3 vs siCTRL"
    ),
  "FigureS5E_CAT_NMD3_curve.csv"
)

readr::write_csv(
  cat_xpo1_S5E$ranked_data,
  "FigureS5E_CAT_XPO1_ranked_genes.csv"
)

readr::write_csv(
  cat_nmd3_S5E$ranked_data,
  "FigureS5E_CAT_NMD3_ranked_genes.csv"
)

cat(
  "\nFigure S5E completed.\n",
  "KPT-330 N = ",
  cat_xpo1_S5E$statistics$N,
  "; siNMD3 N = ",
  cat_nmd3_S5E$statistics$N,
  "\n",
  sep = ""
)



}))

# =============================================================================
# Fig3E_F
# =============================================================================
nmd3_step("Fig3E_F", character(), quote({
# Fig3E_F — analysis below
# RPL vs RPS nuclear-retention shift
#
# DeltaDelta =
# (N-C)Treatment - (N-C)Control
#
# Positive value = increased nuclear retention
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

# prot_raw was loaded previously

# ============================================================
# 1. Protein DeltaDelta function
# ============================================================

run_protein_dd <- function(
    prot_df,
    ctrl_pattern,
    treat_pattern
) {
  
  ctrl_cols <- grep(
    ctrl_pattern,
    colnames(prot_df),
    value = TRUE
  )
  
  treat_cols <- grep(
    treat_pattern,
    colnames(prot_df),
    value = TRUE
  )
  
  sample_cols <- c(
    ctrl_cols,
    treat_cols
  )
  
  stopifnot(
    length(ctrl_cols) == 6,
    length(treat_cols) == 6
  )
  
  # Sample annotation
  
  meta <- data.frame(
    Sample = sample_cols,
    
    Treatment = ifelse(
      sample_cols %in% treat_cols,
      "Treat",
      "Ctrl"
    ),
    
    Compartment = ifelse(
      grepl("_Nuc_", sample_cols),
      "Nucleus",
      "Cytoplasm"
    ),
    
    Replicate = sub(
      ".*_[A-Z]([123])[CN]$",
      "\\1",
      sample_cols
    ),
    
    stringsAsFactors = FALSE,
    row.names = sample_cols
  )
  
  meta$Treatment <- factor(
    meta$Treatment,
    levels = c("Ctrl", "Treat")
  )
  
  meta$Compartment <- factor(
    meta$Compartment,
    levels = c("Cytoplasm", "Nucleus")
  )
  
  meta$Replicate <- factor(
    meta$Replicate
  )
  
  # ----------------------------------------------------------
  # Protein expression matrix
  # ----------------------------------------------------------
  
  expr_df <- prot_df %>%
    dplyr::filter(!is.na(Gene)) %>%
    dplyr::select(
      Gene,
      dplyr::all_of(sample_cols)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(sample_cols),
        ~ suppressWarnings(as.numeric(.x))
      )
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(sample_cols),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  expr_df <- as.data.frame(expr_df)
  
  rownames(expr_df) <- make.unique(
    expr_df$Gene
  )
  
  expr_mat <- as.matrix(
    expr_df[
      ,
      sample_cols,
      drop = FALSE
    ]
  )
  
  # ----------------------------------------------------------
  # Require >=2 measurements in each of the four conditions
  # ----------------------------------------------------------
  
  enough_obs <- function(comp, trt) {
    
    idx <- which(
      meta$Compartment == comp &
        meta$Treatment == trt
    )
    
    rowSums(
      !is.na(
        expr_mat[
          ,
          idx,
          drop = FALSE
        ]
      )
    ) >= 2
  }
  
  keep <- (
    enough_obs("Cytoplasm", "Ctrl") &
      enough_obs("Cytoplasm", "Treat") &
      enough_obs("Nucleus", "Ctrl") &
      enough_obs("Nucleus", "Treat")
  )
  
  expr_mat <- expr_mat[
    keep,
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # Limma interaction model
  # ----------------------------------------------------------
  
  design <- model.matrix(
    ~ Replicate +
      Compartment * Treatment,
    data = meta
  )
  
  colnames(design) <- make.names(
    colnames(design)
  )
  
  interaction_coef <- grep(
    "CompartmentNucleus.*TreatmentTreat",
    colnames(design),
    value = TRUE
  )
  
  stopifnot(
    length(interaction_coef) == 1
  )
  
  fit <- limma::lmFit(
    expr_mat,
    design
  )
  
  fit <- limma::eBayes(
    fit
  )
  
  tt <- limma::topTable(
    fit,
    coef = interaction_coef,
    number = Inf,
    sort.by = "none"
  )
  
  data.frame(
    Gene = rownames(tt),
    logFC = tt$logFC,
    P.Value = tt$P.Value,
    adj.P.Val = tt$adj.P.Val,
    row.names = NULL
  )
}

# ============================================================
# 2. KPT-330 vs DMSO
# ============================================================

res_dd_kpt <- run_protein_dd(
  prot_df = prot_raw,
  
  ctrl_pattern =
    "^Norm_DMSO_(Cyto|Nuc)_",
  
  treat_pattern =
    "^Norm_XPO1inhibitor_(Cyto|Nuc)_"
)

# ============================================================
# 3. siNMD3 vs siCTRL
# ============================================================

res_dd_si <- run_protein_dd(
  prot_df = prot_raw,
  
  ctrl_pattern =
    "^Norm_siRNAxCtrl_(Cyto|Nuc)_",
  
  treat_pattern =
    "^Norm_siRNAxNMD3xKD_(Cyto|Nuc)_"
)
# =============================================================================
# PART 3: NUCLEAR RETENTION FOUR-GROUP STATISTICS (REUSE LIMMA FITS)
# =============================================================================


# =============================================================================
# Nuclear retention: All proteins, Other, RPS, RPL
# XPO1i vs DMSO and siNMD3 vs siCtrl
#
# TWO outputs from one complete analysis:
#   1) Violin/boxplots + signed median difference + P values
#   2) Same violin/boxplots WITHOUT effect size or P annotations
#
# Positive DeltaDelta = (N-C)Treat - (N-C)Ctrl = relative nuclear retention
# Delta median = median(Group2) - median(Group1), ORIGINAL y-axis units
# NOTE: "Delta median" is an unstandardized effect size, not Cohen's d.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(limma)
  library(ggplot2)
})

# In-house DIA protein workbook was read once in the core analysis.
stopifnot(exists("prot_raw"), "Gene" %in% names(prot_raw))

# Shared limma interaction results already obtained in core Figure 4D;
# no additional fitted model and no second workbook read.
stopifnot(exists("res_dd_kpt"), exists("res_dd_si"))
res_kpt_4groups <- res_dd_kpt
res_si_4groups <- res_dd_si

# -----------------------------------------------------------------------------
# Step 4. RP families: final definition used in Figure 4D
#          RPL = ^RPL; RPS = ^RPS excluding ^RPS6K (kinases)
# -----------------------------------------------------------------------------

classify_family <- function(gene) {
  g <- toupper(trimws(as.character(gene)))
  dplyr::case_when(
    grepl("^RPL", g) ~ "RPL",
    grepl("^RPS", g) & !grepl("^RPS6K", g) ~ "RPS",
    TRUE ~ "Other"
  )
}

group_levels <- c("All proteins", "Other", "RPS", "RPL")
comparison_levels <- c("XPO1i vs DMSO", "siNMD3 vs siCtrl")

prepare_groups <- function(results, comparison) {
  one_each <- results %>%
    dplyr::filter(!is.na(Gene), is.finite(logFC)) %>%
    dplyr::mutate(Family = classify_family(Gene))
  
  dplyr::bind_rows(
    one_each %>% dplyr::mutate(Group = "All proteins"),
    one_each %>% dplyr::mutate(Group = Family)
  ) %>%
    dplyr::mutate(
      Group = factor(Group, levels = group_levels),
      Comparison = factor(comparison, levels = comparison_levels)
    )
}

violin_df <- dplyr::bind_rows(
  prepare_groups(res_kpt_4groups, comparison_levels[1]),
  prepare_groups(res_si_4groups, comparison_levels[2])
)

# -----------------------------------------------------------------------------
# Step 5. Descriptive statistics
# -----------------------------------------------------------------------------

group_summary <- violin_df %>%
  dplyr::group_by(Comparison, Group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_DeltaDelta = median(logFC),
    mean_DeltaDelta = mean(logFC),
    SD_DeltaDelta = sd(logFC),
    .groups = "drop"
  )

cat("\n===== Group summary =====\n")
print(as.data.frame(group_summary), row.names = FALSE)

# -----------------------------------------------------------------------------
# Step 6. Five comparisons per treatment
#
# Other/RPS/RPL groups are disjoint: two-sided Wilcoxon rank-sum P.
# "All proteins" overlaps with each RP group. For comparisons vs All,
# use the same subset-permutation test as the previous four-group figure.
# Delta median is always median(Group2) - median(Group1).
# -----------------------------------------------------------------------------

set.seed(20261023)
B_permutation <- 10000L

subset_perm_p <- function(all_values, subset_values, B = B_permutation) {
  med_all <- median(all_values)
  observed <- abs(median(subset_values) - med_all)
  null_dist <- replicate(B, {
    s <- sample(all_values, length(subset_values), replace = FALSE)
    abs(median(s) - med_all)
  })
  (1 + sum(null_dist >= observed - 1e-12)) / (B + 1)
}

# The first line is above the top of the RPL violin in the supplied PDF.
# Five separate horizontal lines, fully inside the -4 to +4 plotting range.
comparison_pairs <- data.frame(
  Group1 = c("RPS", "Other", "Other", "All proteins", "All proteins"),
  Group2 = c("RPL", "RPS", "RPL", "RPS", "RPL"),
  y_line = c(2.45, 2.76, 3.07, 3.38, 3.69),
  stringsAsFactors = FALSE
)

calculate_group_stats <- function(comp) {
  dat <- violin_df %>%
    dplyr::filter(as.character(Comparison) == comp)
  
  out <- lapply(seq_len(nrow(comparison_pairs)), function(i) {
    g1 <- comparison_pairs$Group1[i]
    g2 <- comparison_pairs$Group2[i]
    x <- dat$logFC[as.character(dat$Group) == g1]
    y <- dat$logFC[as.character(dat$Group) == g2]
    stopifnot(length(x) > 1L, length(y) > 1L)
    
    if (g1 == "All proteins") {
      p <- subset_perm_p(x, y)
      method <- "Subset permutation (overlapping All group)"
    } else {
      p <- stats::wilcox.test(x, y, exact = FALSE)$p.value
      method <- "Wilcoxon rank-sum"
    }
    
    data.frame(
      Comparison = comp,
      Group1 = g1,
      Group2 = g2,
      n1 = length(x),
      n2 = length(y),
      median1 = median(x),
      median2 = median(y),
      Delta_median_Group2_minus_Group1 = median(y) - median(x),
      Test = method,
      P_value = p,
      y_line = comparison_pairs$y_line[i],
      stringsAsFactors = FALSE
    )
  })
  
  out <- dplyr::bind_rows(out)
  out$FDR_BH <- stats::p.adjust(out$P_value, method = "BH")
  out
}

stats_4groups <- dplyr::bind_rows(
  lapply(comparison_levels, calculate_group_stats)
)

format_p3 <- function(p) {
  ifelse(p < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", p)))
}

stats_4groups <- stats_4groups %>%
  dplyr::mutate(
    Comparison = factor(Comparison, levels = comparison_levels),
    x_start = match(Group1, group_levels),
    x_end = match(Group2, group_levels),
    x_middle = (x_start + x_end) / 2,
    y_label = y_line + 0.12,
    Stat_label = paste0(
      "\u0394median = ",
      sprintf("%+.3f", Delta_median_Group2_minus_Group1),
      "; ", format_p3(P_value)
    )
  )

cat("\n===== Effect size and P values =====\n")
print(as.data.frame(stats_4groups[, c(
  "Comparison", "Group1", "Group2", "n1", "n2",
  "Delta_median_Group2_minus_Group1", "Test", "P_value", "FDR_BH"
)]), row.names = FALSE)

# -----------------------------------------------------------------------------
# Step 7. Shared plotting layers: same appearance for the two PDFs
# -----------------------------------------------------------------------------

group_colors <- c(
  "All proteins" = "#858585",
  "Other" = "#65A392",
  "RPS" = "#4575B4",
  "RPL" = "#D73027"
)

make_retention_plot <- function(show_stats = TRUE, comp_filter = NULL) {
  plot_violin_df <- if (is.null(comp_filter)) violin_df else
    violin_df[as.character(violin_df$Comparison) == comp_filter, , drop = FALSE]
  plot_group_stats <- if (is.null(comp_filter)) stats_4groups else
    stats_4groups[as.character(stats_4groups$Comparison) == comp_filter, , drop = FALSE]
  p <- ggplot2::ggplot(
    plot_violin_df,
    ggplot2::aes(x = Group, y = logFC, fill = Group)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey55",
      linewidth = 0.5
    ) +
    ggplot2::geom_violin(
      width = 0.85,
      alpha = 0.30,
      color = NA,
      scale = "width",
      trim = TRUE
    ) +
    # Points first so the black boxplot is visible on top.
    ggplot2::geom_point(
      ggplot2::aes(color = Group),
      position = ggplot2::position_jitter(
        width = 0.10, height = 0, seed = 42
      ),
      size = 0.45,
      alpha = 0.12
    ) +
    ggplot2::geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.8,
      alpha = 0.85
    ) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_y_continuous(breaks = seq(-4, 4, by = 2)) +
    ggplot2::coord_cartesian(ylim = c(-4, 4), expand = FALSE) +
    ggplot2::facet_grid(Comparison ~ ., scales = "fixed") +
    ggplot2::labs(
      x = NULL,
      y = expression(Delta*Delta~"nuclear-cytoplasmic protein signal"),
      title = "Nuclear retention: All proteins, Other, RPS, RPL"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
      panel.spacing.y = grid::unit(0.6, "cm"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
  
  if (show_stats) {
    p <- p +
      # Horizontal lines only: no brackets or vertical tips.
      ggplot2::geom_segment(
        data = plot_group_stats,
        mapping = ggplot2::aes(
          x = x_start, xend = x_end,
          y = y_line, yend = y_line
        ),
        inherit.aes = FALSE,
        linewidth = 0.45, color = "black"
      ) +
      ggplot2::geom_text(
        data = plot_group_stats,
        mapping = ggplot2::aes(
          x = x_middle, y = y_label, label = Stat_label
        ),
        inherit.aes = FALSE,
        size = 2.8,
        color = "black",
        vjust = 0
      )
  }
  
  p
}

# -----------------------------------------------------------------------------
# Step 8. WITH effect size and P values
# -----------------------------------------------------------------------------

# Same underlying 10,000-permutation tests; split the two PPT panels only.
p_with_effect_P <- make_retention_plot(show_stats = TRUE)

p_fig3e <- make_retention_plot(TRUE, comparison_levels[1L])
p_fig3f <- make_retention_plot(TRUE, comparison_levels[2L])
nmd3_save_pdf("Fig3E_4Groups.pdf", p_fig3e,
                width = 10.5, height = 5.7, device = grDevices::cairo_pdf)
nmd3_save_pdf("Fig3F_siNMD3_4Groups.pdf", p_fig3f,
                width = 10.5, height = 5.7, device = grDevices::cairo_pdf)

# -----------------------------------------------------------------------------
# Step 9. WITHOUT effect size or P values (same data and scales)
# -----------------------------------------------------------------------------

p_without_stats <- make_retention_plot(show_stats = FALSE)

# -----------------------------------------------------------------------------
# Step 10. Export all numerical summaries / tests
# -----------------------------------------------------------------------------

utils::write.csv(
  group_summary,
  "Nuclear_Retention_4Groups_GroupSummary.csv",
  row.names = FALSE
)
utils::write.csv(
  stats_4groups,
  "Nuclear_Retention_4Groups_EffectSize_Pvalues.csv",
  row.names = FALSE
)

}))

# =============================================================================
# Fig4A
# =============================================================================
nmd3_step("Fig4A", c("Fig3A"), quote({
# Fig4A — analysis below
# Shared preparation
#
# IMPORTANT:
# Reuse the finalized Figure 4 siNMD3 RNA-protein analysis.
#
# Required Figure 4 objects:
#   rna_si
#   prot_si
#   multi_si
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(forcats)
  library(stringr)
  library(msigdbr)
  library(clusterProfiler)
  library(ggplot2)
  library(patchwork)
  library(readr)
})


# ============================================================
# 1. Reuse Figure 4 objects
# ============================================================

stopifnot(
  exists("rna_si"),
  exists("prot_si"),
  exists("multi_si")
)


stopifnot(
  all(
    c(
      "gene",
      "log2FC_RNA",
      "pval_RNA",
      "padj_RNA",
      "log2FC_Prot",
      "pval_Prot",
      "padj_Prot"
    ) %in%
      colnames(multi_si)
  )
)


stopifnot(
  all(
    c(
      "gene",
      "log2FC_Prot",
      "pval_Prot",
      "padj_Prot"
    ) %in%
      colnames(prot_si)
  )
)


# ============================================================
# 2. Matched RNA-protein universe for S6A-B
#
# Same siNMD3 vs siCTRL cytoplasmic dataset
# already used in Figure 4.
# ============================================================

df_S6 <- multi_si %>%
  
  dplyr::transmute(
    
    gene = toupper(
      trimws(
        as.character(gene)
      )
    ),
    
    log2FC_RNA =
      suppressWarnings(
        as.numeric(log2FC_RNA)
      ),
    
    pvalue_RNA =
      suppressWarnings(
        as.numeric(pval_RNA)
      ),
    
    padj_RNA =
      suppressWarnings(
        as.numeric(padj_RNA)
      ),
    
    log2FC_Prot =
      suppressWarnings(
        as.numeric(log2FC_Prot)
      ),
    
    pvalue_Prot =
      suppressWarnings(
        as.numeric(pval_Prot)
      ),
    
    padj_Prot =
      suppressWarnings(
        as.numeric(padj_Prot)
      )
  ) %>%
  
  dplyr::filter(
    gene != "",
    is.finite(log2FC_RNA),
    is.finite(log2FC_Prot)
  ) %>%
  
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )


stopifnot(
  nrow(df_S6) > 100,
  !anyDuplicated(df_S6$gene)
)


cat(
  "\nS6 matched RNA-protein genes:",
  nrow(df_S6),
  "\n"
)


# ============================================================
# 3. Protein universe for S6C
#
# Reuse Figure 4 prot_si directly.
# This retains the complete-case proteomics analysis.
# ============================================================

protein_result_S6 <- prot_si %>%
  
  dplyr::transmute(
    
    gene = toupper(
      trimws(
        as.character(gene)
      )
    ),
    
    log2FC_Prot =
      suppressWarnings(
        as.numeric(log2FC_Prot)
      ),
    
    pvalue_Prot =
      suppressWarnings(
        as.numeric(pval_Prot)
      ),
    
    padj_Prot =
      suppressWarnings(
        as.numeric(padj_Prot)
      )
  ) %>%
  
  dplyr::filter(
    gene != "",
    is.finite(log2FC_Prot)
  ) %>%
  
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )


cat(
  "S6C complete-case proteins:",
  nrow(protein_result_S6),
  "\n"
)


# ============================================================
# 4. Hallmark gene sets
# ============================================================

Hallmark_sets_S6 <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "H"
) %>%
  
  dplyr::transmute(
    
    pathway = gs_name,
    
    gene = toupper(
      gene_symbol
    )
  ) %>%
  
  dplyr::distinct()


Hallmark_names_S6 <- sort(
  unique(
    Hallmark_sets_S6$pathway
  )
)


universe_S6 <- unique(
  df_S6$gene
)


cat(
  "Hallmark pathways:",
  length(Hallmark_names_S6),
  "\n"
)


# ============================================================
# 5. Clean pathway labels
# ============================================================

clean_Hallmark_label_S6 <- function(x) {
  
  x %>%
    
    stringr::str_remove(
      "^HALLMARK_"
    ) %>%
    
    stringr::str_replace_all(
      "_",
      " "
    ) %>%
    
    stringr::str_to_title() %>%
    
    stringr::str_replace_all(
      "\\bMyc\\b",
      "MYC"
    ) %>%
    
    stringr::str_replace_all(
      "\\bDna\\b",
      "DNA"
    ) %>%
    
    stringr::str_replace_all(
      "\\bE2f\\b",
      "E2F"
    ) %>%
    
    stringr::str_replace_all(
      "\\bUv\\b",
      "UV"
    ) %>%
    
    stringr::str_replace_all(
      "\\bTnf\\b",
      "TNF"
    ) %>%
    
    stringr::str_replace_all(
      "\\bNfkb\\b",
      "NFκB"
    ) %>%
    
    stringr::str_replace_all(
      "\\bIl6\\b",
      "IL6"
    ) %>%
    
    stringr::str_replace_all(
      "\\bJak\\b",
      "JAK"
    ) %>%
    
    stringr::str_replace_all(
      "\\bStat3\\b",
      "STAT3"
    ) %>%
    
    stringr::str_replace_all(
      "\\bPi3k\\b",
      "PI3K"
    ) %>%
    
    stringr::str_replace_all(
      "\\bAkt\\b",
      "AKT"
    ) %>%
    
    stringr::str_replace_all(
      "\\bMtorc1\\b",
      "MTORC1"
    ) %>%
    
    stringr::str_replace_all(
      "\\bKras\\b",
      "KRAS"
    )
}


# ============================================================
# 6. Hallmark pathway statistics
#
# pathway genes vs all other matched genes
# ============================================================

calculate_pathway_statistics_S6 <- function(
    pathway_name
) {
  
  pathway_genes <- Hallmark_sets_S6 %>%
    
    dplyr::filter(
      pathway == pathway_name
    ) %>%
    
    dplyr::pull(
      gene
    ) %>%
    
    unique()
  
  
  pathway_genes <- intersect(
    pathway_genes,
    universe_S6
  )
  
  
  other_genes <- setdiff(
    universe_S6,
    pathway_genes
  )
  
  
  if (
    length(pathway_genes) < 10 ||
    length(other_genes) < 10
  ) {
    
    return(
      tibble::tibble(
        pathway = pathway_name,
        n_genes = length(pathway_genes),
        delta_rna = NA_real_,
        p_rna = NA_real_,
        delta_protein = NA_real_,
        p_protein = NA_real_
      )
    )
  }
  
  
  RNA_in <- df_S6$log2FC_RNA[
    df_S6$gene %in%
      pathway_genes
  ]
  
  
  RNA_out <- df_S6$log2FC_RNA[
    df_S6$gene %in%
      other_genes
  ]
  
  
  protein_in <- df_S6$log2FC_Prot[
    df_S6$gene %in%
      pathway_genes
  ]
  
  
  protein_out <- df_S6$log2FC_Prot[
    df_S6$gene %in%
      other_genes
  ]
  
  
  tibble::tibble(
    
    pathway = pathway_name,
    
    n_genes =
      length(RNA_in),
    
    
    delta_rna =
      stats::median(
        RNA_in,
        na.rm = TRUE
      ) -
      stats::median(
        RNA_out,
        na.rm = TRUE
      ),
    
    
    p_rna =
      suppressWarnings(
        stats::wilcox.test(
          RNA_in,
          RNA_out,
          exact = FALSE
        )$p.value
      ),
    
    
    delta_protein =
      stats::median(
        protein_in,
        na.rm = TRUE
      ) -
      stats::median(
        protein_out,
        na.rm = TRUE
      ),
    
    
    p_protein =
      suppressWarnings(
        stats::wilcox.test(
          protein_in,
          protein_out,
          exact = FALSE
        )$p.value
      )
  )
}


# ============================================================
# 7. Calculate all Hallmark pathways
# ============================================================

pathway_statistics_S6 <- purrr::map_dfr(
  Hallmark_names_S6,
  calculate_pathway_statistics_S6
) %>%
  
  dplyr::mutate(
    
    fdr_rna = stats::p.adjust(
      p_rna,
      method = "BH"
    ),
    
    fdr_protein = stats::p.adjust(
      p_protein,
      method = "BH"
    ),
    
    pathway_label =
      clean_Hallmark_label_S6(
        pathway
      )
  )


cat(
  "\nHallmark pathway statistics calculated:",
  nrow(pathway_statistics_S6),
  "\n"
)


readr::write_csv(
  pathway_statistics_S6,
  "FigureS6_Hallmark_pathway_statistics.csv"
)


cat(
  "\nS6A-C shared preparation complete.\n"
)

# ============================================================
# Fig4A — analysis below
# Hallmark pathways with RNA FDR >= 0.05
# ============================================================

plot_data_S6A <- pathway_statistics_S6 %>%
  filter(
    !is.na(delta_rna),
    !is.na(fdr_rna),
    fdr_rna >= 0.05
  ) %>%
  mutate(
    pathway_label = fct_reorder(
      pathway_label,
      delta_rna
    )
  )

cat(
  "RNA-stable Hallmark pathways:",
  nrow(plot_data_S6A),
  "\n"
)

p_S6A <- ggplot(
  plot_data_S6A,
  aes(
    x = delta_rna,
    y = pathway_label
  )
) +
  annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey92"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    stroke = 0.45,
    fill = "grey75",
    color = "black"
  ) +
  labs(
    title = expression(
      bolditalic(siNMD3)~"vs siCtrl"
    ),
    subtitle = "RNA pathway-level FDR ≥ 0.05",
    x = expression(
      Delta~"median RNA log"[2]*"FC"
    ),
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    axis.text.y = element_text(
      size = 8,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    panel.grid = element_blank()
  )


# ============================================================

nmd3_save_pdf("Fig4A_RNA_Stable_Hallmark.pdf", p_S6A, width=6.5, height=8, device=grDevices::cairo_pdf)

}))

# =============================================================================
# Fig4B
# =============================================================================
nmd3_step("Fig4B", c("Fig4A"), quote({
# Fig4B — analysis below
# Residual changes among RNA-stable Hallmark pathways
# ============================================================

# ------------------------------------------------------------
# 5. Global RNA-protein regression
# ------------------------------------------------------------

RNA_protein_model_S6 <- lm(
  log2FC_Prot ~ log2FC_RNA,
  data = df_S6
)

print(
  summary(RNA_protein_model_S6)
)

df_residual_S6 <- df_S6 %>%
  mutate(
    expected_log2FC_Prot = predict(
      RNA_protein_model_S6,
      newdata = df_S6
    ),
    protein_RNA_residual =
      log2FC_Prot -
      expected_log2FC_Prot
  )

RNA_stable_pathways_S6 <-
  pathway_statistics_S6 %>%
  filter(
    !is.na(fdr_rna),
    fdr_rna >= 0.05
  ) %>%
  pull(pathway) %>%
  unique()

calculate_residual_statistics_S6 <- function(
    pathway_name
) {
  
  pathway_genes <- Hallmark_sets_S6 %>%
    filter(
      pathway == pathway_name
    ) %>%
    pull(gene) %>%
    unique()
  
  pathway_genes <- intersect(
    pathway_genes,
    df_residual_S6$gene
  )
  
  other_genes <- setdiff(
    df_residual_S6$gene,
    pathway_genes
  )
  
  if (
    length(pathway_genes) < 10 ||
    length(other_genes) < 10
  ) {
    
    return(
      tibble(
        pathway = pathway_name,
        n_genes = length(pathway_genes),
        delta_residual = NA_real_,
        p_residual = NA_real_
      )
    )
  }
  
  residual_in <-
    df_residual_S6$protein_RNA_residual[
      df_residual_S6$gene %in%
        pathway_genes
    ]
  
  residual_out <-
    df_residual_S6$protein_RNA_residual[
      df_residual_S6$gene %in%
        other_genes
    ]
  
  tibble(
    pathway = pathway_name,
    n_genes = length(residual_in),
    
    delta_residual =
      median(
        residual_in,
        na.rm = TRUE
      ) -
      median(
        residual_out,
        na.rm = TRUE
      ),
    
    p_residual = suppressWarnings(
      wilcox.test(
        residual_in,
        residual_out
      )$p.value
    )
  )
}

residual_statistics_S6 <- map_dfr(
  RNA_stable_pathways_S6,
  calculate_residual_statistics_S6
) %>%
  mutate(
    fdr_residual = p.adjust(
      p_residual,
      method = "BH"
    )
  ) %>%
  left_join(
    pathway_statistics_S6 %>%
      dplyr::select(
        pathway,
        pathway_label,
        delta_rna,
        fdr_rna
      ),
    by = "pathway"
  )

write_csv(
  residual_statistics_S6,
  "FigureS6B_residual_statistics.csv"
)

plot_data_S6B <- residual_statistics_S6 %>%
  filter(
    !is.na(delta_residual),
    !is.na(fdr_residual),
    fdr_residual < 0.05
  ) %>%
  mutate(
    residual_direction = ifelse(
      delta_residual > 0,
      "Residual higher",
      "Residual lower"
    ),
    
    FDR_label = paste0(
      "FDR = ",
      signif(
        fdr_residual,
        2
      )
    ),
    
    label_y =
      delta_residual +
      ifelse(
        delta_residual > 0,
        0.045,
        -0.045
      ),
    
    pathway_label = fct_reorder(
      pathway_label,
      delta_residual
    )
  )

cat(
  "Significant residual pathways:",
  nrow(plot_data_S6B),
  "\n"
)

print(
  plot_data_S6B,
  n = Inf
)

if (nrow(plot_data_S6B) == 0) {
  stop(
    "No residual-significant RNA-stable pathways."
  )
}

p_S6B <- ggplot(
  plot_data_S6B,
  aes(x = pathway_label)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      xend = pathway_label,
      y = delta_rna,
      yend = delta_residual
    ),
    color = "grey65",
    linewidth = 0.7
  ) +
  geom_point(
    aes(y = delta_rna),
    shape = 21,
    size = 4.4,
    stroke = 0.9,
    fill = "white",
    color = "black"
  ) +
  geom_point(
    aes(
      y = delta_residual,
      fill = residual_direction
    ),
    shape = 21,
    size = 4.6,
    stroke = 0.4,
    color = "black"
  ) +
  geom_text(
    aes(
      y = label_y,
      label = FDR_label
    ),
    size = 3.2,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Residual lower" = "#4DBBD5",
      "Residual higher" = "#E64B35"
    ),
    name = NULL
  ) +
  labs(
    title =
      "RNA-stable Hallmark pathways",
    subtitle =
      "Open: RNA change | Filled: protein–RNA residual",
    x = NULL,
    y = expression(
      Delta~"median change"
    )
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    panel.grid = element_blank(),
    plot.margin = margin(
      18, 12, 8, 12
    )
  )


# ============================================================

nmd3_save_pdf("Fig4B_Protein_RNA_Residual_Pathways.pdf", p_S6B, width=7.5, height=5.5, device=grDevices::cairo_pdf)

}))

# =============================================================================
# Fig4C
# =============================================================================
nmd3_step("Fig4C", c("Fig4B"), quote({
# Fig4C — analysis below
# Original complete-case proteomics GSEA
# ============================================================

# ------------------------------------------------------------
# 6. Pure protein-log2FC ranking
# ------------------------------------------------------------

rank_df_S6C <- protein_result_S6 %>%
  transmute(
    gene = toupper(
      trimws(gene)
    ),
    rank_value = log2FC_Prot
  ) %>%
  filter(
    gene != "",
    is.finite(rank_value)
  ) %>%
  group_by(gene) %>%
  summarise(
    rank_value = median(
      rank_value,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

ranking_S6C <- rank_df_S6C$rank_value
names(ranking_S6C) <- rank_df_S6C$gene

ranking_S6C <- sort(
  ranking_S6C,
  decreasing = TRUE
)

cat(
  "Proteins in S6C ranking:",
  length(ranking_S6C),
  "\n"
)

# ------------------------------------------------------------
# 7. Original safe GSEA function
# ------------------------------------------------------------

run_GSEA_S6C <- function(
    gene_ranking,
    term2gene,
    minSize = 10,
    maxSize = 3000,
    seed_number = 12345,
    jitter_sd = 1e-6
) {
  
  term2gene_clean <- term2gene %>%
    transmute(
      term = as.character(term),
      gene = toupper(
        as.character(gene)
      )
    ) %>%
    distinct()
  
  pathway_list <- split(
    term2gene_clean$gene,
    term2gene_clean$term
  )
  
  pathway_list <- lapply(
    pathway_list,
    function(x) {
      intersect(
        names(gene_ranking),
        x
      )
    }
  )
  
  pathway_sizes <- vapply(
    pathway_list,
    length,
    integer(1)
  )
  
  pathway_list <- pathway_list[
    pathway_sizes >= minSize &
      pathway_sizes <= maxSize
  ]
  
  has_variation <- vapply(
    pathway_list,
    function(x) {
      
      values <- gene_ranking[x]
      
      values <- values[
        is.finite(values)
      ]
      
      length(values) >= minSize &&
        length(unique(values)) > 1
    },
    logical(1)
  )
  
  pathway_list <-
    pathway_list[has_variation]
  
  if (length(pathway_list) == 0) {
    stop("No usable pathways.")
  }
  
  set.seed(seed_number)
  
  jitter_values <- rnorm(
    length(gene_ranking),
    sd = jitter_sd
  )
  
  names(jitter_values) <-
    names(gene_ranking)
  
  ranking_jittered <-
    gene_ranking +
    jitter_values
  
  ranking_jittered <- sort(
    ranking_jittered,
    decreasing = TRUE
  )
  
  matched_TERM2GENE <- utils::stack(
    pathway_list
  )
  
  colnames(matched_TERM2GENE) <- c(
    "gene",
    "term"
  )
  
  matched_TERM2GENE <-
    matched_TERM2GENE %>%
    distinct()
  
  clusterProfiler::GSEA(
    geneList = ranking_jittered,
    TERM2GENE = matched_TERM2GENE[
      ,
      c("term", "gene")
    ],
    minGSSize = minSize,
    maxGSSize = maxSize,
    pvalueCutoff = 1,
    seed = TRUE,
    verbose = FALSE
  )
}

# ------------------------------------------------------------
# 8. Gene-set collections
# ------------------------------------------------------------

GO_BP_S6C <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
) %>%
  transmute(
    term = gs_name,
    gene = gene_symbol
  ) %>%
  distinct()

GO_CC_S6C <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:CC"
) %>%
  transmute(
    term = gs_name,
    gene = gene_symbol
  ) %>%
  distinct()

Reactome_S6C <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
) %>%
  transmute(
    term = gs_name,
    gene = gene_symbol
  ) %>%
  distinct()

Hallmark_S6C <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "H"
) %>%
  transmute(
    term = gs_name,
    gene = gene_symbol
  ) %>%
  distinct()

# ------------------------------------------------------------
# 9. Run GSEA
# ------------------------------------------------------------

gsea_GO_BP_S6C <- run_GSEA_S6C(
  ranking_S6C,
  GO_BP_S6C
)

gsea_GO_CC_S6C <- run_GSEA_S6C(
  ranking_S6C,
  GO_CC_S6C
)

gsea_Reactome_S6C <- run_GSEA_S6C(
  ranking_S6C,
  Reactome_S6C
)

gsea_Hallmark_S6C <- run_GSEA_S6C(
  ranking_S6C,
  Hallmark_S6C
)

# ------------------------------------------------------------
# 10. Prespecified degradation/proteostasis pathways
# ------------------------------------------------------------

selected_terms_S6C <- tribble(
  ~Source, ~Term, ~Label,
  
  "GO:BP",
  "GOBP_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS",
  "Proteasomal protein\ncatabolic process",
  
  "GO:CC",
  "GOCC_PROTEASOME_COMPLEX",
  "Proteasome\ncomplex",
  
  "Reactome",
  "REACTOME_PROTEASOME_ASSEMBLY",
  "Proteasome\nassembly",
  
  "GO:BP",
  "GOBP_PROTEIN_CATABOLIC_PROCESS",
  "Protein catabolic\nprocess",
  
  "GO:CC",
  "GOCC_AUTOPHAGOSOME",
  "Autophagosome",
  
  "Hallmark",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "Unfolded protein\nresponse"
)

extract_GSEA_table_S6C <- function(
    gsea_object,
    source
) {
  
  result <- as.data.frame(
    gsea_object@result
  )
  
  if (nrow(result) == 0) {
    return(tibble())
  }
  
  as_tibble(result) %>%
    mutate(
      Source = source,
      Term = as.character(ID),
      core_enrichment =
        as.character(core_enrichment)
    )
}

all_GSEA_results_S6C <- bind_rows(
  extract_GSEA_table_S6C(
    gsea_GO_BP_S6C,
    "GO:BP"
  ),
  extract_GSEA_table_S6C(
    gsea_GO_CC_S6C,
    "GO:CC"
  ),
  extract_GSEA_table_S6C(
    gsea_Reactome_S6C,
    "Reactome"
  ),
  extract_GSEA_table_S6C(
    gsea_Hallmark_S6C,
    "Hallmark"
  )
)

plot_data_S6C <- selected_terms_S6C %>%
  left_join(
    all_GSEA_results_S6C,
    by = c(
      "Source",
      "Term"
    )
  )

missing_terms_S6C <- plot_data_S6C %>%
  filter(is.na(NES)) %>%
  dplyr::select(
    Source,
    Term
  )

if (nrow(missing_terms_S6C) > 0) {
  cat("\nMissing S6C pathways:\n")
  print(
    missing_terms_S6C,
    n = Inf
  )
}

plot_data_S6C <- plot_data_S6C %>%
  filter(
    is.finite(NES),
    !is.na(p.adjust)
  ) %>%
  mutate(
    Direction = ifelse(
      NES > 0,
      "Higher in siNMD3",
      "Lower in siNMD3"
    ),
    
    FDR_label = paste0(
      "FDR = ",
      signif(
        p.adjust,
        2
      )
    ),
    
    Label = factor(
      Label,
      levels = rev(
        selected_terms_S6C$Label
      )
    )
  )

print(
  plot_data_S6C %>%
    dplyr::select(
      Source,
      Term,
      Label,
      setSize,
      NES,
      pvalue,
      p.adjust,
      Direction
    ),
  n = Inf
)

p_S6C <- ggplot(
  plot_data_S6C,
  aes(
    x = Label,
    y = NES
  )
) +
  geom_col(
    aes(fill = Direction),
    width = 0.7,
    color = "grey30",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = FDR_label,
      hjust = ifelse(
        NES < 0,
        1.05,
        -0.05
      )
    ),
    size = 3.2,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Lower in siNMD3" = "#4DBBD5",
      "Higher in siNMD3" = "#E64B35"
    ),
    name = NULL
  ) +
  coord_flip(
    clip = "off"
  ) +
  labs(
    title = expression(
      bolditalic(siNMD3)~
        "vs siCtrl cytoplasmic proteome"
    ),
    x = NULL,
    y = "NES"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(
      size = 9,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    panel.grid = element_blank(),
    plot.margin = margin(
      8, 60, 8, 8
    )
  )


# ============================================================
# Part 7. Save individual panels
# ============================================================

nmd3_save_pdf("Fig4C_Proteostasis_GSEA.pdf", p_S6C, width=7, height=4.5, device=grDevices::cairo_pdf)
write_csv(
  pathway_statistics_S6,
  "FigureS6A_all_Hallmark_statistics.csv"
)

write_csv(
  plot_data_S6A %>%
    mutate(
      pathway_label =
        as.character(pathway_label)
    ),
  "FigureS6A_RNA_stable_Hallmark.csv"
)

write_csv(
  residual_statistics_S6,
  "FigureS6B_all_residual_statistics.csv"
)

write_csv(
  plot_data_S6B %>%
    mutate(
      pathway_label =
        as.character(pathway_label)
    ),
  "FigureS6B_significant_residual_pathways.csv"
)

write_csv(
  plot_data_S6C %>%
    transmute(
      Source,
      Term,
      Label = as.character(Label),
      setSize,
      enrichmentScore,
      NES,
      pvalue,
      p.adjust,
      qvalue,
      rank,
      leading_edge =
        as.character(leading_edge),
      core_enrichment =
        as.character(core_enrichment),
      Direction
    ),
  "FigureS6C_selected_GSEA_results.csv"
)

cat("\nFigure S6A-C complete.\n")


}))

# =============================================================================
# Fig4D
# =============================================================================
nmd3_step("Fig4D", c("Fig3A"), quote({
# Fig4D — analysis below
# HALLMARK MYC TARGETS V1 vs all other genes
# RNA log2FC: siNMD3 vs siCTRL
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(msigdbr)
})

# multi_si was generated previously

# ------------------------------------------------------------
# 1. Prepare RNA-protein matched genes
# ------------------------------------------------------------

df_F <- multi_si %>%
  dplyr::select(
    gene,
    log2FC_RNA,
    log2FC_Prot
  ) %>%
  dplyr::filter(
    is.finite(log2FC_RNA),
    is.finite(log2FC_Prot)
  ) %>%
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )

# ------------------------------------------------------------
# 2. HALLMARK MYC TARGETS V1
# ------------------------------------------------------------

myc_genes <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "H"
) %>%
  dplyr::filter(
    gs_name == "HALLMARK_MYC_TARGETS_V1"
  ) %>%
  dplyr::pull(gene_symbol) %>%
  toupper() %>%
  unique()

df_F <- df_F %>%
  dplyr::mutate(
    group = ifelse(
      gene %in% myc_genes,
      "MYC",
      "Other"
    ),
    group = factor(
      group,
      levels = c(
        "Other",
        "MYC"
      )
    )
  )

print(table(df_F$group))

# ------------------------------------------------------------
# 3. Wilcoxon test
# ------------------------------------------------------------

wt_F <- stats::wilcox.test(
  log2FC_RNA ~ group,
  data = df_F,
  exact = FALSE
)

median_other <- median(
  df_F$log2FC_RNA[
    df_F$group == "Other"
  ],
  na.rm = TRUE
)

median_myc <- median(
  df_F$log2FC_RNA[
    df_F$group == "MYC"
  ],
  na.rm = TRUE
)

delta_F <- median_myc - median_other

cat(
  "Wilcoxon p =",
  wt_F$p.value,
  "\nMedian Other =",
  median_other,
  "\nMedian MYC =",
  median_myc,
  "\nDelta median (MYC - Other) =",
  delta_F,
  "\n"
)

# ------------------------------------------------------------
# 4. Statistical label
# ------------------------------------------------------------

p_text_F <- if (
  wt_F$p.value < 0.001
) {
  "p < 0.001"
} else {
  paste0(
    "p = ",
    signif(wt_F$p.value, 2)
  )
}

stat_label_F <- paste0(
  p_text_F,
  "\n\u0394 = ",
  sprintf("%.2f", abs(delta_F))
)

# ------------------------------------------------------------
# 5. Bracket position
# ------------------------------------------------------------


bracket_F <- data.frame(
  group1 = "Other",
  group2 = "MYC",
  y.position = 1.05,
  label = stat_label_F
)

# ------------------------------------------------------------
# 6. Figure F
# ------------------------------------------------------------

p_F <- ggplot(
  df_F,
  aes(
    x = group,
    y = log2FC_RNA,
    fill = group
  )
) +
  geom_violin(
    color = NA,
    alpha = 0.85,
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "Other" = "#4DBBD5",
      "MYC" = "#F39B7F"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Other" = "All other\nproteins",
      "MYC" = "Hallmark MYC-\ntarget proteins"
    )
  ) +
  ggpubr::stat_pvalue_manual(
    bracket_F,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.03,
    bracket.size = 0.8,
    size = 5
  ) +
  coord_cartesian(
    ylim = c(-2.2, 2.2),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = expression(
      "RNA log"[2]*
        "FC (siNMD3 vs siCtrl)"
    )
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    )
  )


# ------------------------------------------------------------
# 7. Save PDF
# ------------------------------------------------------------

nmd3_save_pdf(
  "Fig4D_MYC_Targets_RNA_log2FC.pdf",
  p_F,
  width = 5,
  height = 5,
  device = cairo_pdf
)



}))

# =============================================================================
# Fig4E
# =============================================================================
nmd3_step("Fig4E", c("Fig4D"), quote({
# Fig4E — analysis below
# MYC TARGETS V1 protein-RNA residual
# ============================================================


stopifnot(exists("multi_si"))

df0 <- multi_si %>%
  dplyr::transmute(
    gene = toupper(as.character(gene)),
    log2FC_RNA = as.numeric(log2FC_RNA),
    log2FC_Prot = as.numeric(log2FC_Prot)
  ) %>%
  dplyr::filter(
    is.finite(log2FC_RNA),
    is.finite(log2FC_Prot)
  ) %>%
  dplyr::mutate(
    family = dplyr::case_when(
      grepl("^RPL", gene) ~ "RPL",
      grepl("^RPS", gene) ~ "RPS",
      TRUE ~ "Other"
    )
  )

cat("Genes used for Figure G:", nrow(df0), "\n")
winsor_vec <- function(x, p = 0.01) {
  lower <- stats::quantile(x, p, na.rm = TRUE)
  upper <- stats::quantile(x, 1 - p, na.rm = TRUE)
  
  pmin(pmax(x, lower), upper)
}

# Use non-ribosomal proteins as the regression background
background_pool <- df0 %>%
  dplyr::filter(
    family == "Other",
    is.finite(log2FC_RNA),
    is.finite(log2FC_Prot)
  )

# Retain the central 20-80% of RNA changes for model fitting
rna_limits <- stats::quantile(
  background_pool$log2FC_RNA,
  probs = c(0.20, 0.80),
  na.rm = TRUE
)

background_fit <- background_pool %>%
  dplyr::filter(
    log2FC_RNA >= rna_limits[1],
    log2FC_RNA <= rna_limits[2]
  ) %>%
  dplyr::transmute(
    x = winsor_vec(log2FC_RNA, p = 0.01),
    y = winsor_vec(log2FC_Prot, p = 0.01)
  ) %>%
  tidyr::drop_na()

# Robust protein-RNA regression
fit_residual_old <- MASS::rlm(
  y ~ x,
  data = background_fit,
  psi = MASS::psi.huber
)

print(coef(fit_residual_old))

# Predict protein change from RNA change
x_all <- winsor_vec(df0$log2FC_RNA, p = 0.01)

predicted_protein <- predict(
  fit_residual_old,
  newdata = data.frame(x = x_all)
)

raw_residual <- df0$log2FC_Prot - predicted_protein

# Center residuals using the median of non-ribosomal proteins
other_median <- median(
  raw_residual[df0$family == "Other"],
  na.rm = TRUE
)

df0$residual_old <- raw_residual - other_median

df_G <- df0 %>%
  dplyr::mutate(
    group = ifelse(
      gene %in% myc_genes,
      "MYC",
      "Other"
    ),
    group = factor(
      group,
      levels = c("Other", "MYC")
    )
  )

# ------------------------------------------------------------
# Statistics
# ------------------------------------------------------------

wt_G <- stats::wilcox.test(
  residual_old ~ group,
  data = df_G,
  exact = FALSE
)

median_other <- median(
  df_G$residual_old[df_G$group == "Other"],
  na.rm = TRUE
)

median_myc <- median(
  df_G$residual_old[df_G$group == "MYC"],
  na.rm = TRUE
)

delta_G <- median_myc - median_other

cat(
  "Wilcoxon p =", wt_G$p.value,
  "\nSigned delta =", delta_G,
  "\nAbsolute delta =", abs(delta_G),
  "\n"
)

stat_label_G <- paste0(
  ifelse(
    wt_G$p.value < 0.001,
    "p < 0.001",
    paste0("p = ", signif(wt_G$p.value, 2))
  ),
  "\n\u0394 = ",
  sprintf("%.2f", abs(delta_G))
)

# ------------------------------------------------------------
# Five strongest highlighted genes
# ------------------------------------------------------------

genes_to_label <- c(
  "GSPT1",
  "EIF3D",
  "EEF1B2",
  "CCT5",
  "KPNB1"
)

df_label_G <- df_G %>%
  dplyr::filter(
    group == "MYC",
    gene %in% genes_to_label
  ) %>%
  dplyr::arrange(
    dplyr::desc(residual_old)
  )

print(
  df_label_G %>%
    dplyr::select(
      gene,
      residual_old
    )
)

# ------------------------------------------------------------
# Statistical bracket
# ------------------------------------------------------------

y_range_G <- range(
  df_G$residual_old,
  na.rm = TRUE
)

y_span_G <- diff(y_range_G)

bracket_G <- data.frame(
  group1 = "Other",
  group2 = "MYC",
  y.position =
    y_range_G[2] +
    0.08 * y_span_G,
  label = stat_label_G
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_G <- ggplot(
  df_G,
  aes(
    x = group,
    y = residual_old,
    fill = group
  )
) +
  geom_violin(
    color = NA,
    alpha = 0.85,
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  geom_point(
    data = df_label_G,
    shape = 21,
    size = 2.2,
    fill = "gold",
    color = "black",
    stroke = 0.4
  ) +
  
  ggrepel::geom_text_repel(
    data = df_label_G,
    aes(label = gene),
    nudge_x = 0.40,
    direction = "y",
    hjust = 0,
    size = 4,
    box.padding = 0.10,
    point.padding = 0.08,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  
  scale_fill_manual(
    values = c(
      "Other" = "#4DBBD5",
      "MYC" = "#F39B7F"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "Other" = "All other\nproteins",
      "MYC" = "Hallmark MYC-\ntarget proteins"
    )
  ) +
  
  ggpubr::stat_pvalue_manual(
    bracket_G,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.03,
    bracket.size = 0.8,
    size = 5
  ) +
  
  coord_cartesian(
    ylim = c(-2.2, 2.4),
    clip = "off"
  ) +
  
  labs(
    x = NULL,
    y = "Protein-RNA residual\n(siNMD3 vs siCtrl)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      face = "bold",
      color = "black"
    ),
    plot.margin = margin(
      10, 70, 10, 10
    )
  )


nmd3_save_pdf(
  filename = "Fig4E_MYC_Targets_Protein_RNA_Residual.pdf",
  plot = p_G,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)





}))

# =============================================================================
# Fig6A
# =============================================================================
nmd3_step("Fig6A", character(), quote({
# Fig6A — analysis below
#
# Fig6A — analysis below
# Fig6A — analysis below
#
# NMD3 protein abundance is extracted in the same public proteomics dataset
# and retained as `nmd3_data` for reuse in the RIGHT PANEL of Figure S4A,
# where ONLY Primary and Metastatic samples are compared.
#
# Statistical results are saved as CSV files but are not displayed in figures.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

# -----------------------------------------------------------------------------
# 1. Input and output paths
# -----------------------------------------------------------------------------

input_file <- file.path(
  "data",
  "public_proteomics",
  "10780432ccr181229-sup-200506_2_supp_4875649_pbgdn5.xlsx"
)

sheet_name <- "Table S1"

output_dir <- file.path(
  "results",
  "public_proteomics"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_file), dir.exists(output_dir))

# -----------------------------------------------------------------------------
# 2. Read the original supplementary protein table
# -----------------------------------------------------------------------------

raw_table <- readxl::read_excel(
  input_file,
  sheet = sheet_name,
  col_names = FALSE,
  .name_repair = "minimal"
)

stopifnot(nrow(raw_table) > 3L, ncol(raw_table) >= 63L)

# Row 2 contains sample/feature names; row 3 contains cohort labels.
true_titles <- as.character(unlist(raw_table[2, ], use.names = FALSE))
true_titles <- make.unique(
  ifelse(
    is.na(true_titles) | trimws(true_titles) == "",
    paste0("unnamed_", seq_along(true_titles)),
    trimws(true_titles)
  )
)

group_labels <- trimws(
  as.character(unlist(raw_table[3, ], use.names = FALSE))
)

protein_table <- raw_table[-seq_len(3), ] %>%
  as.data.frame(stringsAsFactors = FALSE)

colnames(protein_table) <- true_titles

stopifnot("Gene name" %in% colnames(protein_table))

# -----------------------------------------------------------------------------
# 3. Construct and validate sample metadata
# -----------------------------------------------------------------------------

sample_columns <- true_titles[6:63]

stopifnot(
  length(sample_columns) == 58L,
  all(sample_columns %in% colnames(protein_table)),
  !anyDuplicated(sample_columns)
)

sample_metadata <- tibble(
  Sample = sample_columns,
  Source_group = group_labels[6:63]
) %>%
  mutate(
    Group = case_when(
      Source_group == "C" ~ "Normal",
      Source_group == "Loc" ~ "Primary",
      str_detect(Source_group, "^Met") ~ "Metastatic",
      TRUE ~ NA_character_
    ),
    Group = factor(
      Group,
      levels = c("Normal", "Primary", "Metastatic")
    )
  )

sample_counts <- table(sample_metadata$Group, useNA = "ifany")
cat("\nAll sample counts:\n")
print(sample_counts)

stopifnot(
  !anyNA(sample_metadata$Group),
  unname(sample_counts["Normal"]) == 8L,
  unname(sample_counts["Primary"]) == 28L,
  unname(sample_counts["Metastatic"]) == 22L,
  !anyDuplicated(sample_metadata$Sample)
)

write.csv(
  sample_metadata,
  file.path(output_dir, "01_sample_metadata.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Construct one gene-level protein-expression matrix
# -----------------------------------------------------------------------------

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x)
}

expression_wide <- protein_table %>%
  transmute(
    Gene = toupper(trimws(as.character(.data[["Gene name"]]))),
    across(
      all_of(sample_columns),
      ~ suppressWarnings(as.numeric(as.character(.x)))
    )
  ) %>%
  filter(!is.na(Gene), Gene != "", Gene != "NAN") %>%
  group_by(Gene) %>%
  summarise(
    across(all_of(sample_columns), safe_median),
    .groups = "drop"
  )

stopifnot(!anyDuplicated(expression_wide$Gene))

all_protein_matrix <- as.matrix(
  expression_wide[, sample_columns, drop = FALSE]
)

storage.mode(all_protein_matrix) <- "double"
rownames(all_protein_matrix) <- expression_wide$Gene

stopifnot(
  ncol(all_protein_matrix) == 58L,
  identical(colnames(all_protein_matrix), sample_metadata$Sample),
  "NMD3" %in% rownames(all_protein_matrix)
)

cat("\nUnique quantified proteins:", nrow(all_protein_matrix), "\n")

# -----------------------------------------------------------------------------
# 5. Define GO:0000054 and select a fixed complete protein panel
# -----------------------------------------------------------------------------

go_bp_sets <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

ribosome_export_genes <- go_bp_sets %>%
  filter(gs_name == "GOBP_RIBOSOMAL_SUBUNIT_EXPORT_FROM_NUCLEUS") %>%
  transmute(Gene = toupper(gene_symbol)) %>%
  distinct() %>%
  pull(Gene)

stopifnot(length(ribosome_export_genes) > 0L)

represented_export_genes <- intersect(
  ribosome_export_genes,
  rownames(all_protein_matrix)
)

export_matrix <- all_protein_matrix[
  represented_export_genes,
  sample_metadata$Sample,
  drop = FALSE
]

# The same proteins must be available in every one of the 58 samples.
complete_export_genes <- rownames(export_matrix)[
  rowSums(is.finite(export_matrix)) == ncol(export_matrix)
]

complete_export_matrix <- export_matrix[
  complete_export_genes,
  ,
  drop = FALSE
]

protein_sd <- apply(complete_export_matrix, 1, stats::sd)

score_genes <- names(
  protein_sd[is.finite(protein_sd) & protein_sd > 0]
)

score_matrix <- complete_export_matrix[score_genes, , drop = FALSE]

stopifnot(
  length(represented_export_genes) >= 5L,
  length(score_genes) >= 5L,
  all(is.finite(score_matrix))
)

cat(
  "\nGO:0000054 genes:", length(ribosome_export_genes),
  "\nRepresented in proteomics:", length(represented_export_genes),
  "\nUsed in the 58-sample score:", length(score_genes),
  "\n"
)
print(score_genes)

gene_set_qc <- tibble(
  Gene = ribosome_export_genes,
  Represented_in_proteomics = Gene %in% represented_export_genes,
  Used_in_score = Gene %in% score_genes
)

write.csv(
  gene_set_qc,
  file.path(output_dir, "02_ribosome_export_gene_set_QC.csv"),
  row.names = FALSE
)

write.csv(
  tibble(Gene = score_genes),
  file.path(output_dir, "03_proteins_used_for_score.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Calculate ribosome-export protein score across all 58 samples
# -----------------------------------------------------------------------------

# Z-standardize each protein across the same 58 samples, then take the
# unweighted mean of the protein z-scores for each sample.
protein_z_matrix <- t(
  scale(t(score_matrix), center = TRUE, scale = TRUE)
)

stopifnot(all(is.finite(protein_z_matrix)))

ribosome_export_score <- colMeans(protein_z_matrix)

ribosome_data <- tibble(
  Sample = names(ribosome_export_score),
  Value = as.numeric(ribosome_export_score)
) %>%
  left_join(sample_metadata, by = "Sample") %>%
  transmute(
    Sample,
    Source_group,
    Group,
    Measure = "Ribosome-export protein score",
    Value
  )

stopifnot(
  nrow(ribosome_data) == 58L,
  all(is.finite(ribosome_data$Value))
)

# -----------------------------------------------------------------------------
# 7. Extract NMD3 protein abundance without imputation
# -----------------------------------------------------------------------------

nmd3_data <- tibble(
  Sample = sample_metadata$Sample,
  Value = as.numeric(
    all_protein_matrix["NMD3", sample_metadata$Sample]
  )
) %>%
  left_join(sample_metadata, by = "Sample") %>%
  filter(is.finite(Value)) %>%
  transmute(
    Sample,
    Source_group,
    Group,
    Measure = "NMD3 protein",
    Value
  )

nmd3_missing_samples <- sample_metadata %>%
  filter(!Sample %in% nmd3_data$Sample)

cat("\nNMD3 valid sample counts:\n")
print(table(nmd3_data$Group))
cat("\nSamples missing NMD3:\n")
print(nmd3_missing_samples, n = Inf)

stopifnot(
  nrow(nmd3_data) == 56L,
  sum(nmd3_data$Group == "Normal") == 8L,
  sum(nmd3_data$Group == "Primary") == 28L,
  sum(nmd3_data$Group == "Metastatic") == 20L,
  setequal(
    nmd3_missing_samples$Sample,
    c("Ratio PC32", "Ratio PC34")
  )
)

write.csv(
  nmd3_missing_samples,
  file.path(output_dir, "04_NMD3_missing_samples.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8. Ribosome-export source data for Figures 2K-2L
#
# IMPORTANT: NMD3 is NOT included in the Figure 2K-2L three-group analysis.
# `nmd3_data` is retained only for the separate two-group Figure S4A analysis.
# -----------------------------------------------------------------------------

measure_data <- ribosome_data %>%
  mutate(
    Group = factor(Group, levels = c("Normal", "Primary", "Metastatic")),
    Measure = factor(
      Measure,
      levels = "Ribosome-export protein score"
    )
  )

stopifnot(
  nrow(measure_data) == 58L,
  all(as.character(measure_data$Measure) == "Ribosome-export protein score")
)

write.csv(
  measure_data,
  file.path(output_dir, "05_Figure2K_2L_ribosome_export_source_data.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 9. Official seven-pair mapping from Supplementary Table S1
# -----------------------------------------------------------------------------

pair_key <- tribble(
  ~Pair_ID, ~Normal,      ~Primary,
  "a",      "Ratio PC2",  "Ratio PC3",
  "b",      "Ratio PC11", "Ratio PC12",
  "c",      "Ratio PC18", "Ratio PC19",
  "d",      "Ratio PC25", "Ratio PC26",
  "e",      "Ratio PC57", "Ratio PC58",
  "f",      "Ratio PC61", "Ratio PC62",
  "g",      "Ratio PC63", "Ratio PC64"
)

pair_design <- pair_key %>%
  pivot_longer(
    cols = c(Normal, Primary),
    names_to = "Group",
    values_to = "Sample"
  )

paired_data <- pair_design %>%
  left_join(
    measure_data %>% mutate(Group = as.character(Group)),
    by = c("Sample", "Group")
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Normal", "Primary")),
    Measure = factor(
      Measure,
      levels = c("Ribosome-export protein score", "NMD3 protein")
    )
  ) %>%
  arrange(Measure, Pair_ID, Group)

pair_check <- paired_data %>%
  dplyr::count(Measure, Pair_ID, name = "n_samples")

stopifnot(
  nrow(paired_data) == 14L,
  !anyNA(paired_data$Value),
  nrow(pair_check) == 7L,
  all(pair_check$n_samples == 2L)
)

paired_ribosome_data <- paired_data %>%
  dplyr::filter(Measure == "Ribosome-export protein score")

stopifnot(
  nrow(paired_ribosome_data) == 14L,
  all(is.finite(paired_ribosome_data$Value))
)

write.csv(
  paired_ribosome_data,
  file.path(output_dir, "06_ribosome_export_seven_paired_source_data.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 10. Descriptive statistics for Figures 2K-2L ribosome-export score
# -----------------------------------------------------------------------------

three_group_descriptive <- measure_data %>%
  group_by(Measure, Group) %>%
  summarise(
    n = n(),
    median = median(Value),
    mean = mean(Value),
    SD = sd(Value),
    .groups = "drop"
  )

paired_descriptive <- paired_ribosome_data %>%
  group_by(Measure, Group) %>%
  summarise(
    n = n(),
    median = median(Value),
    mean = mean(Value),
    SD = sd(Value),
    .groups = "drop"
  )

write.csv(
  three_group_descriptive,
  file.path(output_dir, "07_three_group_descriptive_statistics.csv"),
  row.names = FALSE
)

write.csv(
  paired_descriptive,
  file.path(output_dir, "08_paired_descriptive_statistics.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 11. Seven-pair two-sided Wilcoxon signed-rank tests
# -----------------------------------------------------------------------------

run_paired_wilcoxon <- function(measure_name) {
  wide_data <- paired_data %>%
    filter(Measure == measure_name) %>%
    dplyr::select(Pair_ID, Group, Value) %>%
    pivot_wider(names_from = Group, values_from = Value) %>%
    arrange(Pair_ID)

  stopifnot(
    nrow(wide_data) == 7L,
    all(complete.cases(wide_data))
  )

  paired_difference <- wide_data$Primary - wide_data$Normal

  # Base R can calculate an exact signed-rank P value only without zero
  # differences or tied absolute differences.
  exact_allowed <-
    !any(paired_difference == 0) &&
    !anyDuplicated(abs(paired_difference))

  test_result <- stats::wilcox.test(
    x = wide_data$Primary,
    y = wide_data$Normal,
    paired = TRUE,
    alternative = "two.sided",
    exact = exact_allowed
  )

  tibble(
    Measure = measure_name,
    Comparison = "Primary vs Normal",
    Design = "Two-sided paired Wilcoxon signed-rank test",
    n_pairs = nrow(wide_data),
    Median_Normal = median(wide_data$Normal),
    Median_Primary = median(wide_data$Primary),
    Median_paired_change = median(paired_difference),
    Wilcoxon_V = unname(test_result$statistic),
    Exact_P_used = exact_allowed,
    Raw_P = test_result$p.value
  )
}

paired_statistics <- run_paired_wilcoxon(
  "Ribosome-export protein score"
)

write.csv(
  paired_statistics,
  file.path(
    output_dir,
    "09_ribosome_export_seven_paired_Wilcoxon_statistics.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 12. Prespecified statistics for Figure 2L ribosome-export score
#
# Normal–Primary uses only the seven confirmed pairs.
# Primary–Metastatic and Normal–Metastatic use unpaired rank-sum tests.
# No ordinary Kruskal–Wallis test is used because the design is partially paired.
# NMD3 is intentionally excluded here and is analyzed separately in Figure S4A.
# -----------------------------------------------------------------------------

run_unpaired_wilcoxon <- function(measure_name, group_1, group_2) {
  current_data <- measure_data %>%
    filter(Measure == measure_name)

  x <- current_data %>%
    filter(Group == group_1) %>%
    pull(Value)

  y <- current_data %>%
    filter(Group == group_2) %>%
    pull(Value)

  stopifnot(
    length(x) > 0L,
    length(y) > 0L,
    all(is.finite(x)),
    all(is.finite(y))
  )

  test_result <- stats::wilcox.test(
    x = x,
    y = y,
    paired = FALSE,
    alternative = "two.sided",
    exact = FALSE
  )

  tibble(
    Measure = measure_name,
    Comparison = paste(group_1, "vs", group_2),
    Design = "Two-sided unpaired Wilcoxon rank-sum test",
    n_group_1 = length(x),
    n_group_2 = length(y),
    Median_group_1 = median(x),
    Median_group_2 = median(y),
    Median_difference = median(x) - median(y),
    Wilcoxon_W = unname(test_result$statistic),
    Raw_P = test_result$p.value
  )
}

three_group_statistics <- bind_rows(
  lapply(
    levels(measure_data$Measure),
    function(measure_name) {
      paired_row <- run_paired_wilcoxon(measure_name) %>%
        transmute(
          Measure,
          Comparison,
          Design,
          n_group_1 = n_pairs,
          n_group_2 = n_pairs,
          Median_group_1 = Median_Primary,
          Median_group_2 = Median_Normal,
          Median_difference = Median_paired_change,
          Wilcoxon_W = Wilcoxon_V,
          Raw_P
        )

      bind_rows(
        paired_row,
        run_unpaired_wilcoxon(
          measure_name,
          "Primary",
          "Metastatic"
        ),
        run_unpaired_wilcoxon(
          measure_name,
          "Normal",
          "Metastatic"
        )
      )
    }
  )
) %>%
  group_by(Measure) %>%
  mutate(BH_adjusted_P = p.adjust(Raw_P, method = "BH")) %>%
  ungroup()

write.csv(
  three_group_statistics,
  file.path(output_dir, "10_three_group_prespecified_statistics.csv"),
  row.names = FALSE
)

cat("\nSeven-pair statistics:\n")
print(paired_statistics, n = Inf)

cat("\nThree-group prespecified statistics:\n")
print(three_group_statistics, n = Inf)

# -----------------------------------------------------------------------------
# 13. Plot settings
# -----------------------------------------------------------------------------

group_colors <- c(
  "Normal" = "#91D1C2FF",
  "Primary" = "#4DBBD5",
  "Metastatic" = "#E64B35FF"
)

publication_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    plot.margin = margin(8, 10, 8, 8)
  )

# -----------------------------------------------------------------------------
# 14. Seven-pair plotting function
# -----------------------------------------------------------------------------

plot_paired_measure <- function(measure_name, y_axis_title) {
  plot_data <- paired_data %>%
    filter(Measure == measure_name) %>%
    group_by(Pair_ID) %>%
    mutate(
      normal_value = Value[Group == "Normal"][1],
      point_fill = case_when(
        Group == "Normal" ~ "Normal",
        Group == "Primary" & Value > normal_value ~ "Primary_higher",
        Group == "Primary" & Value < normal_value ~ "Primary_lower",
        Group == "Primary" & Value == normal_value ~ "Primary_equal"
      )
    ) %>%
    ungroup()
  
  stopifnot(nrow(plot_data) == 14L)
  
  ggplot(
    plot_data,
    aes(x = Group, y = Value, group = Pair_ID)
  ) +
    geom_line(
      color = "grey55",
      linewidth = 0.65,
      alpha = 0.80
    ) +
    geom_point(
      aes(fill = point_fill),
      shape = 21,
      size = 3.5,
      stroke = 0.75,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "Normal" = "#91D1C2FF",   # normal 保持原样
        "Primary_higher" = "red",            # tumor > normal
        "Primary_lower" = "blue",            # tumor < normal
        "Primary_equal" = group_colors["Primary"]  # 如果刚好相等，保持原来的 primary 色
      ),
      drop = FALSE
    ) +
    scale_x_discrete(limits = c("Normal", "Primary")) +
    labs(x = NULL, y = y_axis_title) +
    publication_theme
}

# -----------------------------------------------------------------------------
# 15. Three-group plotting function
# -----------------------------------------------------------------------------

plot_three_groups <- function(measure_name, y_axis_title) {
  plot_data <- measure_data %>%
    filter(Measure == measure_name)

  ggplot(
    plot_data,
    aes(x = Group, y = Value, fill = Group)
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      width = 0.80,
      color = NA,
      alpha = 0.65
    ) +
    geom_boxplot(
      width = 0.28,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.65
    ) +
    geom_jitter(
      width = 0.13,
      height = 0,
      size = 1.5,
      alpha = 0.75,
      color = "black",
      show.legend = FALSE
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE)  +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_x_discrete(limits = c("Normal", "Primary", "Metastatic")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    labs(x = NULL, y = y_axis_title) +
    publication_theme
}

# -----------------------------------------------------------------------------
# 16. Figure 2K: Seven paired Normal–Primary ribosome-export scores
# -----------------------------------------------------------------------------

p_2K <- plot_paired_measure(
  "Ribosome-export protein score",
  "Ribosome-export protein score"
)


# -----------------------------------------------------------------------------
save_pdf <- function(plot, filename, width, height) {
  nmd3_save_pdf(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}
save_pdf(
  p_2K,
  "Fig6A_Proteomics_Ribosome_Export_Score_7_Pairs.pdf",
  3.4,
  4.3
)

}))

# =============================================================================
# Fig6B
# =============================================================================
nmd3_step("Fig6B", c("Fig6A"), quote({
# 17. Figure 2L: Three-group ribosome-export score
# -----------------------------------------------------------------------------

p_2L <- plot_three_groups(
  "Ribosome-export protein score",
  "Ribosome-export protein score"
)


save_pdf(
  p_2L,
  "Fig6B_Proteomics_Ribosome_Export_Score_Three_Groups.pdf",
  4.5,
  4.5
)

}))

# =============================================================================
# Prepare_TCGA
# =============================================================================
nmd3_step("Prepare_TCGA", character(), quote({
# Prepare_TCGA — analysis below
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(janitor)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Input paths
# -----------------------------------------------------------------------------

sample_sheet_path <- "data/TCGA_PRAD/gdc_sample_sheet.2024-08-29.tsv"
main_folder <- "data/TCGA_PRAD/counts"

stopifnot(
  file.exists(sample_sheet_path),
  dir.exists(main_folder)
)

# -----------------------------------------------------------------------------
# Shared functions
# -----------------------------------------------------------------------------

# Read one GDC count file and collapse IDs duplicated by _PAR_Y
read_gdc_counts <- function(file) {
  read_tsv(
    file,
    skip = 1,
    show_col_types = FALSE
  ) %>%
    filter(!startsWith(gene_id, "N_")) %>%
    transmute(
      gene_id = sub("\\..*$", "", gene_id),
      count = as.numeric(unstranded)
    ) %>%
    group_by(gene_id) %>%
    summarise(
      count = sum(count),
      .groups = "drop"
    )
}

# Build a count matrix without rereading the first sample
build_count_matrix <- function(sample_files, sample_ids) {
  stopifnot(
    length(sample_files) == length(sample_ids),
    length(sample_files) > 0,
    !anyDuplicated(sample_ids)
  )
  
  first_count <- read_gdc_counts(sample_files[1])
  genes_master <- first_count$gene_id
  
  count_matrix <- vapply(
    seq_along(sample_files),
    function(j) {
      count_table <- if (j == 1) {
        first_count
      } else {
        read_gdc_counts(sample_files[j])
      }
      
      aligned <- setNames(
        count_table$count,
        count_table$gene_id
      )[genes_master]
      
      if (anyNA(aligned)) {
        stop("Gene IDs are inconsistent among GDC count files.")
      }
      
      aligned
    },
    numeric(length(genes_master))
  )
  
  dimnames(count_matrix) <- list(
    genes_master,
    sample_ids
  )
  
  storage.mode(count_matrix) <- "integer"
  count_matrix
}

# Map Ensembl IDs to gene symbols and combine duplicate symbols
map_counts_to_symbols <- function(count_matrix) {
  ensembl_to_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = rownames(count_matrix),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  mapped <- !is.na(ensembl_to_symbol)
  
  count_matrix <- rowsum(
    count_matrix[mapped, , drop = FALSE],
    group = ensembl_to_symbol[mapped]
  )
  
  storage.mode(count_matrix) <- "integer"
  count_matrix
}

# Shared GSEA parameters
run_gsea <- function(gene_ranking, term_to_gene) {
  GSEA(
    geneList = gene_ranking,
    TERM2GENE = term_to_gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE,
    seed = TRUE
  )
}

# -----------------------------------------------------------------------------
# Read TCGA sample information
# -----------------------------------------------------------------------------

sample_info <- read_tsv(
  sample_sheet_path,
  show_col_types = FALSE
) %>%
  clean_names()

tsv_files <- list.files(
  main_folder,
  pattern = "\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

stopifnot(
  length(tsv_files) > 0,
  !anyDuplicated(basename(tsv_files))
)

# Prepare_TCGA — analysis below
# =============================================================================

primary_samples <- sample_info %>%
  filter(
    sample_type == "Primary Tumor",
    file_name %in% basename(tsv_files)
  ) %>%
  mutate(
    file = tsv_files[
      match(file_name, basename(tsv_files))
    ]
  ) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  ) %>%
  transmute(
    sample_id,
    file
  )

stopifnot(
  nrow(primary_samples) > 0,
  !anyNA(primary_samples$file)
)

primary_counts <- build_count_matrix(
  sample_files = primary_samples$file,
  sample_ids = primary_samples$sample_id
)

primary_counts <- map_counts_to_symbols(
  primary_counts
)

}))

# =============================================================================
# Fig6C
# =============================================================================
nmd3_step("Fig6C", c("Prepare_TCGA"), quote({
# Fig6C — analysis below
# TCGA-PRAD overall survival by ribosome export ssGSEA score
# Reuses primary_counts generated in Figure 1D
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(GSVA)
  library(msigdbr)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Check the TCGA primary-tumor count matrix
# ------------------------------------------------------------

stopifnot(
  exists("primary_counts"),
  ncol(primary_counts) > 100
)

rownames(primary_counts) <- toupper(
  rownames(primary_counts)
)

storage.mode(primary_counts) <- "integer"

cat(
  "Primary-tumor samples before patient deduplication:",
  ncol(primary_counts),
  "\n"
)

# ------------------------------------------------------------
# 2. Retain one primary-tumor sample per patient
# ------------------------------------------------------------

sample_information_S4E <- tibble::tibble(
  sample_id = colnames(primary_counts),
  patient_id = substr(
    toupper(colnames(primary_counts)),
    1,
    12
  )
) %>%
  dplyr::arrange(
    patient_id,
    sample_id
  ) %>%
  dplyr::distinct(
    patient_id,
    .keep_all = TRUE
  )

counts_S4E <- primary_counts[
  ,
  sample_information_S4E$sample_id,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts_S4E),
    sample_information_S4E$sample_id
  )
)

cat(
  "Primary-tumor patients retained:",
  ncol(counts_S4E),
  "\n"
)

# ------------------------------------------------------------
# 3. Filter low-expression genes
# ------------------------------------------------------------

minimum_count_S4E <- 10
minimum_sample_proportion_S4E <- 0.10

minimum_sample_number_S4E <- ceiling(
  minimum_sample_proportion_S4E *
    ncol(counts_S4E)
)

keep_genes_S4E <- rowSums(
  counts_S4E >= minimum_count_S4E
) >= minimum_sample_number_S4E

counts_filtered_S4E <- counts_S4E[
  keep_genes_S4E,
  ,
  drop = FALSE
]

cat(
  "Genes retained after expression filtering:",
  nrow(counts_filtered_S4E),
  "\n"
)

# ------------------------------------------------------------
# 4. Variance-stabilizing transformation
# ------------------------------------------------------------

coldata_S4E <- data.frame(
  patient_id = factor(
    sample_information_S4E$patient_id
  ),
  row.names = sample_information_S4E$sample_id
)

stopifnot(
  identical(
    colnames(counts_filtered_S4E),
    rownames(coldata_S4E)
  )
)

dds_S4E <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_filtered_S4E,
  colData = coldata_S4E,
  design = ~ 1
)

vsd_S4E <- DESeq2::vst(
  dds_S4E,
  blind = TRUE
)

expression_vst_S4E <- SummarizedExperiment::assay(
  vsd_S4E
)

cat(
  "VST expression matrix:",
  nrow(expression_vst_S4E),
  "genes x",
  ncol(expression_vst_S4E),
  "patients\n"
)

# ------------------------------------------------------------
# 5. Retrieve the ribosomal-subunit-export gene set
# ------------------------------------------------------------

gobp_sets_S4E <- nmd3_cached_msig(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

ribosome_export_genes_S4E <- gobp_sets_S4E %>%
  dplyr::filter(
    gs_name ==
      "GOBP_RIBOSOMAL_SUBUNIT_EXPORT_FROM_NUCLEUS"
  ) %>%
  dplyr::pull(gene_symbol) %>%
  toupper() %>%
  unique()

genes_used_S4E <- intersect(
  ribosome_export_genes_S4E,
  rownames(expression_vst_S4E)
)

cat(
  "Original ribosome-export genes:",
  length(ribosome_export_genes_S4E),
  "\n"
)

cat(
  "Ribosome-export genes found in TCGA:",
  length(genes_used_S4E),
  "\n"
)

cat(
  "Genes used:",
  paste(
    sort(genes_used_S4E),
    collapse = ", "
  ),
  "\n"
)

stopifnot(
  length(genes_used_S4E) >= 5
)

gene_set_S4E <- list(
  Ribosome_Export = genes_used_S4E
)

# ------------------------------------------------------------
# 6. Calculate single-sample GSEA scores
# ------------------------------------------------------------

ssgsea_parameters_S4E <- GSVA::ssgseaParam(
  exprData = expression_vst_S4E,
  geneSets = gene_set_S4E
)

ssgsea_matrix_S4E <- GSVA::gsva(
  ssgsea_parameters_S4E,
  verbose = FALSE
)

ssgsea_matrix_S4E <- as.matrix(
  ssgsea_matrix_S4E
)

stopifnot(
  "Ribosome_Export" %in%
    rownames(ssgsea_matrix_S4E)
)

ssgsea_data_S4E <- tibble::tibble(
  sample_id = colnames(ssgsea_matrix_S4E),
  patient_id = substr(
    toupper(colnames(ssgsea_matrix_S4E)),
    1,
    12
  ),
  RiboExport_Score = as.numeric(
    ssgsea_matrix_S4E[
      "Ribosome_Export",
      colnames(ssgsea_matrix_S4E)
    ]
  )
) %>%
  dplyr::filter(
    is.finite(RiboExport_Score)
  ) %>%
  dplyr::distinct(
    patient_id,
    .keep_all = TRUE
  )

cat(
  "Patients with ssGSEA scores:",
  nrow(ssgsea_data_S4E),
  "\n"
)

# ------------------------------------------------------------
# 7. Read TCGA clinical survival information
# ------------------------------------------------------------

clinical_path_S4E <- paste0(
  "data/",
  "clinical.project-tcga-prad.2024-09-06/clinical.tsv"
)

clinical_raw_S4E <- readr::read_tsv(
  clinical_path_S4E,
  show_col_types = FALSE
)

required_clinical_columns_S4E <- c(
  "case_submitter_id",
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up"
)

stopifnot(
  all(
    required_clinical_columns_S4E %in%
      colnames(clinical_raw_S4E)
  )
)

survival_information_S4E <- clinical_raw_S4E %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ ifelse(
        .x %in% c(
          "'--",
          "--",
          "Not Reported",
          "not reported"
        ),
        NA,
        .x
      )
    )
  ) %>%
  dplyr::transmute(
    patient_id = toupper(
      trimws(
        as.character(case_submitter_id)
      )
    ),
    vital_status = toupper(
      trimws(
        as.character(vital_status)
      )
    ),
    days_to_death = suppressWarnings(
      as.numeric(days_to_death)
    ),
    days_to_last_follow_up = suppressWarnings(
      as.numeric(days_to_last_follow_up)
    ),
    OS_time = dplyr::if_else(
      !is.na(days_to_death),
      days_to_death,
      days_to_last_follow_up
    ),
    OS_status = dplyr::case_when(
      vital_status == "DEAD" ~ 1,
      vital_status == "ALIVE" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(
    !is.na(patient_id),
    patient_id != "",
    is.finite(OS_time),
    OS_time > 0,
    OS_status %in% c(0, 1)
  ) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::slice_max(
    order_by = OS_time,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

cat(
  "Patients with valid survival information:",
  nrow(survival_information_S4E),
  "\n"
)

cat("Survival status:\n")

print(
  table(
    survival_information_S4E$OS_status
  )
)

# ------------------------------------------------------------
# 8. Merge ssGSEA scores with survival information
# ------------------------------------------------------------

survival_data_S4E <- ssgsea_data_S4E %>%
  dplyr::inner_join(
    survival_information_S4E,
    by = "patient_id"
  ) %>%
  dplyr::filter(
    is.finite(RiboExport_Score),
    is.finite(OS_time),
    OS_time > 0,
    OS_status %in% c(0, 1)
  )

cat(
  "Patients used for Figure S4E:",
  nrow(survival_data_S4E),
  "\n"
)

stopifnot(
  nrow(survival_data_S4E) > 10,
  length(
    unique(
      survival_data_S4E$RiboExport_Score
    )
  ) > 1
)

# ------------------------------------------------------------
# 9. Median-based High/Low grouping
# ------------------------------------------------------------

median_score_S4E <- median(
  survival_data_S4E$RiboExport_Score,
  na.rm = TRUE
)

cat(
  "Median ribosome-export score:",
  median_score_S4E,
  "\n"
)

survival_data_S4E <- survival_data_S4E %>%
  dplyr::mutate(
    Risk_Group = ifelse(
      RiboExport_Score >= median_score_S4E,
      "High",
      "Low"
    ),
    Risk_Group = factor(
      Risk_Group,
      levels = c(
        "Low",
        "High"
      )
    )
  )

group_counts_S4E <- table(
  survival_data_S4E$Risk_Group
)

cat("Median-based groups:\n")
print(group_counts_S4E)

# ------------------------------------------------------------
# 10. Kaplan-Meier and log-rank analyses
# ------------------------------------------------------------

fit_S4E <- survival::survfit(
  survival::Surv(
    OS_time,
    OS_status
  ) ~ Risk_Group,
  data = survival_data_S4E
)

logrank_S4E <- survival::survdiff(
  survival::Surv(
    OS_time,
    OS_status
  ) ~ Risk_Group,
  data = survival_data_S4E
)

logrank_p_S4E <- 1 - stats::pchisq(
  logrank_S4E$chisq,
  df = length(logrank_S4E$n) - 1
)

cat(
  "Log-rank p-value:",
  format(
    logrank_p_S4E,
    scientific = TRUE,
    digits = 8
  ),
  "\n"
)

legend_labels_S4E <- c(
  paste0(
    "Low n = ",
    unname(group_counts_S4E["Low"])
  ),
  paste0(
    "High n = ",
    unname(group_counts_S4E["High"])
  )
)

# ------------------------------------------------------------
# 11. Final Figure S4E
# ------------------------------------------------------------

p_S4E <- survminer::ggsurvplot(
  fit_S4E,
  data = survival_data_S4E,
  
  # Statistical results are calculated but not displayed
  pval = FALSE,
  pval.method = FALSE,
  
  conf.int = TRUE,
  risk.table = FALSE,
  censor = TRUE,
  
  title = "Overall Survival by Ribosome Export Score",
  subtitle = paste0(
    "TCGA-PRAD primary tumors (n = ",
    nrow(survival_data_S4E),
    ")"
  ),
  
  legend.title = "Ribosome Export",
  legend.labs = legend_labels_S4E,
  
  xlab = "Time (days)",
  ylab = "Survival Probability",
  
  palette = c(
    "#1F78B4",
    "#E31A1C"
  ),
  
  font.x = c(
    20,
    "bold",
    "black"
  ),
  font.y = c(
    20,
    "bold",
    "black"
  ),
  font.tickslab = c(
    18,
    "plain",
    "black"
  ),
  font.legend = c(
    18,
    "plain",
    "black"
  ),
  font.main = c(
    20,
    "bold",
    "black"
  ),
  
  ggtheme = theme_minimal(base_size = 18) +
    theme(
      panel.grid.major = element_line(
        color = "grey85",
        linewidth = 0.6
      ),
      panel.grid.minor = element_line(
        color = "grey92",
        linewidth = 0.4
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      axis.text = element_text(
        color = "black"
      ),
      axis.title = element_text(
        color = "black"
      )
    )
)


# ------------------------------------------------------------
# 12. Save outputs
# ------------------------------------------------------------

nmd3_save_pdf(
  filename = "Fig6C_Ribosome_Export_Survival.pdf",
  plot = p_S4E$plot,
  width = 7.5,
  height = 6,
  units = "in",
  device = cairo_pdf
)

write.csv(
  survival_data_S4E,
  "FigureS4E_TCGA_Ribosome_Export_patient_data.csv",
  row.names = FALSE
)

statistics_S4E <- data.frame(
  cohort = "TCGA-PRAD primary tumors",
  total_n = nrow(survival_data_S4E),
  events = sum(
    survival_data_S4E$OS_status == 1
  ),
  gene_set = paste(
    sort(genes_used_S4E),
    collapse = ";"
  ),
  gene_set_size = length(genes_used_S4E),
  cutoff_method = "Median",
  cutoff_value = median_score_S4E,
  low_n = unname(
    group_counts_S4E["Low"]
  ),
  high_n = unname(
    group_counts_S4E["High"]
  ),
  logrank_p = logrank_p_S4E,
  stringsAsFactors = FALSE
)

print(statistics_S4E)

write.csv(
  statistics_S4E,
  "FigureS4E_TCGA_Ribosome_Export_statistics.csv",
  row.names = FALSE
)










}))

# =============================================================================
# Fig6D
# =============================================================================
nmd3_step("Fig6D", character(), quote({
# Fig6D — analysis below
# PEG10 and ENO2 expression in two independent CRPC cohorts
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(Seurat)
  library(pheatmap)
  library(grid)
  library(gridExtra)
})

genes_6A <- c("PEG10", "ENO2")


# ============================================================
# Cohort 1: WCM NEPC 2016
# CRPC-Adeno vs CRPC-NE
# ============================================================

expr_file <- "data/nepc_wcm_2016/data_mrna_seq_fpkm.txt"

clinical_file <- "data/nepc_wcm_2016/data_clinical_sample.txt"


# ------------------------------------------------------------
# 1. Read expression and clinical data
# ------------------------------------------------------------

expr_wcm <- read.delim(
  expr_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

clinical_wcm <- read.delim(
  clinical_file,
  sep = "\t",
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  "Hugo_Symbol" %in% colnames(expr_wcm),
  "SAMPLE_ID" %in% colnames(clinical_wcm),
  "DISEASE_CODE" %in% colnames(clinical_wcm)
)


# ------------------------------------------------------------
# 2. Extract PEG10 and ENO2
# ------------------------------------------------------------

expr_cols <- colnames(expr_wcm)[
  3:ncol(expr_wcm)
]

mat_wcm_df <- expr_wcm %>%
  
  dplyr::filter(
    Hugo_Symbol %in% genes_6A
  ) %>%
  
  dplyr::select(
    Hugo_Symbol,
    dplyr::all_of(expr_cols)
  ) %>%
  
  tibble::column_to_rownames(
    "Hugo_Symbol"
  ) %>%
  
  as.data.frame(
    check.names = FALSE
  )

mat_wcm_df <- mat_wcm_df[
  !duplicated(rownames(mat_wcm_df)),
  ,
  drop = FALSE
]

stopifnot(
  all(
    genes_6A %in%
      rownames(mat_wcm_df)
  )
)

rn <- rownames(mat_wcm_df)
cn <- colnames(mat_wcm_df)

mat_wcm <- as.matrix(
  data.frame(
    lapply(
      mat_wcm_df,
      function(x) {
        suppressWarnings(
          as.numeric(
            as.character(x)
          )
        )
      }
    ),
    check.names = FALSE
  )
)

rownames(mat_wcm) <- rn
colnames(mat_wcm) <- cn

mat_wcm <- mat_wcm[
  genes_6A,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 3. WCM clinical groups
# ------------------------------------------------------------

ann_wcm <- clinical_wcm %>%
  
  dplyr::filter(
    DISEASE_CODE %in%
      c(
        "CRPC-Adeno",
        "CRPC-NE"
      )
  ) %>%
  
  dplyr::transmute(
    
    SAMPLE_ID = SAMPLE_ID,
    
    Group = factor(
      DISEASE_CODE,
      levels = c(
        "CRPC-Adeno",
        "CRPC-NE"
      )
    )
  ) %>%
  
  as.data.frame()

rownames(ann_wcm) <-
  ann_wcm$SAMPLE_ID

ann_wcm$SAMPLE_ID <- NULL


# ------------------------------------------------------------
# 4. Align WCM samples
# ------------------------------------------------------------

common_wcm <- intersect(
  colnames(mat_wcm),
  rownames(ann_wcm)
)

mat_wcm <- mat_wcm[
  ,
  common_wcm,
  drop = FALSE
]

ann_wcm <- ann_wcm[
  common_wcm,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 5. Check WCM sample numbers
# ------------------------------------------------------------

cat("\nWCM NEPC 2016:\n")

print(
  table(
    ann_wcm$Group
  )
)

stopifnot(
  sum(ann_wcm$Group == "CRPC-Adeno") == 34,
  sum(ann_wcm$Group == "CRPC-NE") == 15
)


# ------------------------------------------------------------
# 6. WCM row z-score
# ------------------------------------------------------------

mat_wcm_z <- t(
  scale(
    t(mat_wcm)
  )
)

mat_wcm_z[
  is.na(mat_wcm_z)
] <- 0


# ------------------------------------------------------------
# 7. Order WCM samples
# Group first, then PEG10 + ENO2 score
# ------------------------------------------------------------

score_wcm <- colSums(
  mat_wcm_z
)

ord_wcm <- order(
  ann_wcm$Group,
  score_wcm
)

mat_wcm_z <- mat_wcm_z[
  ,
  ord_wcm,
  drop = FALSE
]

ann_wcm <- ann_wcm[
  ord_wcm,
  ,
  drop = FALSE
]

gap_wcm <- which(
  diff(
    as.integer(
      ann_wcm$Group
    )
  ) != 0
)


# ============================================================
# Cohort 2: SU2C/PCF
# ARPC vs NEPC
# ============================================================


# ------------------------------------------------------------
# 8. Load original processed SU2C object
# ------------------------------------------------------------

loaded_objects <- load(
  "data/prad_su2c_2019/ARPCvNEPCvDNPC_revision.RData"
)

stopifnot(
  "ARNEDNPC" %in%
    loaded_objects
)

obj_su2c <- ARNEDNPC


# ------------------------------------------------------------
# 9. Use original AR/NE classification
# ------------------------------------------------------------

stopifnot(
  "ARNE" %in%
    colnames(
      obj_su2c@meta.data
    )
)

Seurat::Idents(
  obj_su2c
) <- "ARNE"

obj_su2c <- subset(
  obj_su2c,
  idents = c(
    "ARPC",
    "NEPC"
  )
)


# ------------------------------------------------------------
# 10. Check SU2C sample numbers
# ------------------------------------------------------------

cat("\nSU2C/PCF:\n")

print(
  table(
    Seurat::Idents(
      obj_su2c
    )
  )
)

stopifnot(
  sum(
    Seurat::Idents(obj_su2c) ==
      "ARPC"
  ) == 177,
  
  sum(
    Seurat::Idents(obj_su2c) ==
      "NEPC"
  ) == 13
)


# ------------------------------------------------------------
# 11. Extract normalized PEG10 and ENO2 expression
# ------------------------------------------------------------

assay_use <- Seurat::DefaultAssay(
  obj_su2c
)

mat_su2c <- Seurat::GetAssayData(
  obj_su2c,
  assay = assay_use,
  layer = "data"
)

stopifnot(
  all(
    genes_6A %in%
      rownames(mat_su2c)
  )
)

mat_su2c <- mat_su2c[
  genes_6A,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. SU2C annotation
# ------------------------------------------------------------

ann_su2c <- data.frame(
  
  Group = factor(
    as.character(
      Seurat::Idents(
        obj_su2c
      )
    ),
    
    levels = c(
      "ARPC",
      "NEPC"
    )
  )
)

rownames(ann_su2c) <-
  colnames(mat_su2c)


# ------------------------------------------------------------
# 13. SU2C row z-score
# ------------------------------------------------------------

mat_su2c_z <- t(
  scale(
    t(
      as.matrix(
        mat_su2c
      )
    )
  )
)

mat_su2c_z[
  is.na(mat_su2c_z)
] <- 0


# ------------------------------------------------------------
# 14. Order SU2C samples
# Group first, then PEG10 + ENO2 score
# ------------------------------------------------------------

score_su2c <- colSums(
  mat_su2c_z
)

ord_su2c <- order(
  ann_su2c$Group,
  score_su2c
)

mat_su2c_z <- mat_su2c_z[
  ,
  ord_su2c,
  drop = FALSE
]

ann_su2c <- ann_su2c[
  ord_su2c,
  ,
  drop = FALSE
]

gap_su2c <- which(
  diff(
    as.integer(
      ann_su2c$Group
    )
  ) != 0
)


# ============================================================
# Common heatmap settings
# ============================================================

breaks_6A <- seq(
  -2.5,
  2.5,
  length.out = 101
)

colors_6A <- colorRampPalette(
  c(
    "#2C7BB6",
    "white",
    "#D7191C"
  )
)(
  length(breaks_6A) - 1
)


# ============================================================
# 15. WCM heatmap
# ============================================================

ann_colors_wcm <- list(
  
  Group = c(
    "CRPC-Adeno" = "#4DBBD5",
    "CRPC-NE" = "#E64B35"
  )
)

p_wcm <- pheatmap::pheatmap(
  
  mat_wcm_z,
  
  annotation_col =
    ann_wcm,
  
  annotation_colors =
    ann_colors_wcm,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_colnames = FALSE,
  show_rownames = TRUE,
  
  fontsize_row = 12,
  
  color = colors_6A,
  breaks = breaks_6A,
  
  border_color = NA,
  
  gaps_col = gap_wcm,
  
  main = "WCM NEPC 2016",
  
  silent = TRUE
)


# ============================================================
# 16. SU2C/PCF heatmap
# ============================================================

ann_colors_su2c <- list(
  
  Group = c(
    "ARPC" = "#4DBBD5",
    "NEPC" = "#E64B35"
  )
)

p_su2c <- pheatmap::pheatmap(
  
  mat_su2c_z,
  
  annotation_col =
    ann_su2c,
  
  annotation_colors =
    ann_colors_su2c,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_colnames = FALSE,
  show_rownames = TRUE,
  
  fontsize_row = 12,
  
  color = colors_6A,
  breaks = breaks_6A,
  
  border_color = NA,
  
  gaps_col = gap_su2c,
  
  main = "SU2C/PCF",
  
  silent = TRUE
)


# ============================================================
# 17. Show BOTH cohorts
# ============================================================

grid::grid.newpage()

gridExtra::grid.arrange(
  p_wcm$gtable,
  p_su2c$gtable,
  ncol = 1,
  heights = c(
    1,
    1
  )
)


# ============================================================
# 18. Save individual PDFs
# ============================================================

grid::grid.newpage()
grid::grid.draw(
  p_wcm$gtable
)

nmd3_dev_off()

grid::grid.newpage()
grid::grid.draw(
  p_su2c$gtable
)

nmd3_dev_off()


# ============================================================
# 19. Save combined Figure 6A
# ============================================================

nmd3_pdf(
  "Fig6D_WCM_SU2C_PEG10_ENO2_Heatmaps.pdf",
  width = 10,
  height = 5,
  useDingbats = FALSE
)

gridExtra::grid.arrange(
  p_wcm$gtable,
  p_su2c$gtable,
  ncol = 1,
  heights = c(
    1,
    1
  )
)

nmd3_dev_off()



}))

# =============================================================================
# Prepare_GSE35988
# =============================================================================
nmd3_step("Prepare_GSE35988", character(), quote({
# Prepare_GSE35988 — analysis below
# Nuclear export ssGSEA and XPO1 expression
# Complete reproducible version
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(limma)
  library(nlme)
  library(msigdbr)
  library(GSVA)
  library(ggplot2)
})

set.seed(123)

# -----------------------------------------------------------------------------
# 1. Input and output paths
# -----------------------------------------------------------------------------

matrix_file <- paste0(
  "data/GSE35988/",
  "GSE35988-GPL6480_series_matrix.txt.gz"
)

output_dir <- paste0(
  "results/",
  "GSE35988_GPL6480_Figure1EF"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(matrix_file),
  dir.exists(output_dir)
)

# -----------------------------------------------------------------------------
# 2. Load the GPL6480 series matrix
# -----------------------------------------------------------------------------

gpl6480_eset <- GEOquery::getGEO(
  filename = matrix_file,
  getGPL = FALSE
)

if (is.list(gpl6480_eset)) {
  stopifnot(length(gpl6480_eset) == 1L)
  gpl6480_eset <- gpl6480_eset[[1]]
}

expression_probe <- Biobase::exprs(gpl6480_eset)
sample_metadata  <- Biobase::pData(gpl6480_eset)
storage.mode(expression_probe) <- "double"

stopifnot(
  ncol(expression_probe) == 88L,
  nrow(sample_metadata) == 88L,
  identical(colnames(expression_probe), rownames(sample_metadata)),
  "title" %in% colnames(sample_metadata)
)

# The series matrix contains normalized two-color Agilent
# log2(test/reference) values. Do not apply another log2 transformation.

# Keep all 88 samples and one common complete probe set.
probe_qc <- data.frame(
  Probe_ID = rownames(expression_probe),
  n_missing = rowSums(!is.finite(expression_probe)),
  missing_fraction = rowMeans(!is.finite(expression_probe)),
  stringsAsFactors = FALSE
)

sample_qc <- data.frame(
  Sample = colnames(expression_probe),
  n_missing = colSums(!is.finite(expression_probe)),
  missing_fraction = colMeans(!is.finite(expression_probe)),
  stringsAsFactors = FALSE
)

utils::write.csv(
  probe_qc,
  file.path(output_dir, "GSE35988_GPL6480_probe_missingness_QC.csv"),
  row.names = FALSE
)

utils::write.csv(
  sample_qc,
  file.path(output_dir, "GSE35988_GPL6480_sample_missingness_QC.csv"),
  row.names = FALSE
)

complete_probe <- probe_qc$n_missing == 0L
expression_probe <- expression_probe[complete_probe, , drop = FALSE]

stopifnot(
  nrow(expression_probe) == 28796L,
  ncol(expression_probe) == 88L,
  all(is.finite(expression_probe)),
  identical(colnames(expression_probe), rownames(sample_metadata))
)

cat(
  "Complete expression matrix:",
  nrow(expression_probe), "probes x",
  ncol(expression_probe), "samples\n"
)

# -----------------------------------------------------------------------------
# 3. Define groups and patient IDs from the official sample titles
#
# Eight patients have both a Normal and a Primary sample. Those two samples
# receive the same Patient_ID. Every unmatched Normal or Primary sample, and
# every Metastatic CRPC sample, receives its own group-specific Patient_ID.
# -----------------------------------------------------------------------------

sample_info_base <- sample_metadata %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::transmute(
    Sample = as.character(Sample),
    Sample_title = trimws(as.character(title)),
    Group = dplyr::case_when(
      grepl("^N[0-9]", Sample_title, ignore.case = TRUE)  ~ "Normal",
      grepl("^T[0-9]", Sample_title, ignore.case = TRUE)  ~ "Primary",
      grepl("^WA[0-9]", Sample_title, ignore.case = TRUE) ~ "Metastatic CRPC",
      TRUE ~ "Unclassified"
    ),
    Patient_number = sub(
      "^[A-Za-z]+([0-9]+).*$",
      "\\1",
      Sample_title
    )
  )

stopifnot(
  !any(sample_info_base$Group == "Unclassified"),
  !anyNA(sample_info_base$Patient_number),
  !anyDuplicated(sample_info_base$Sample),
  identical(sample_info_base$Sample, colnames(expression_probe))
)

normal_patient_numbers <- sample_info_base %>%
  dplyr::filter(Group == "Normal") %>%
  dplyr::pull(Patient_number) %>%
  unique()

primary_patient_numbers <- sample_info_base %>%
  dplyr::filter(Group == "Primary") %>%
  dplyr::pull(Patient_number) %>%
  unique()

matched_patient_numbers <- intersect(
  normal_patient_numbers,
  primary_patient_numbers
)

stopifnot(length(matched_patient_numbers) == 8L)

sample_info <- sample_info_base %>%
  dplyr::mutate(
    Patient_ID = dplyr::case_when(
      Group %in% c("Normal", "Primary") &
        Patient_number %in% matched_patient_numbers ~
        paste0("NP_PAIR_", Patient_number),
      Group == "Normal" ~ paste0("NORMAL_ONLY_", Patient_number),
      Group == "Primary" ~ paste0("PRIMARY_ONLY_", Patient_number),
      Group == "Metastatic CRPC" ~ paste0("MCRPC_ONLY_", Patient_number),
      TRUE ~ NA_character_
    ),
    Group = factor(
      Group,
      levels = c("Normal", "Primary", "Metastatic CRPC")
    )
  )

stopifnot(
  !anyNA(sample_info$Patient_ID),
  !anyNA(sample_info$Group),
  nrow(sample_info) == 88L,
  dplyr::n_distinct(sample_info$Patient_ID) == 80L
)

group_check <- sample_info %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_patients = dplyr::n_distinct(Patient_ID),
    .groups = "drop"
  )

print(as.data.frame(group_check), row.names = FALSE)

stopifnot(
  group_check$n_samples[group_check$Group == "Normal"] == 12L,
  group_check$n_samples[group_check$Group == "Primary"] == 49L,
  group_check$n_samples[group_check$Group == "Metastatic CRPC"] == 27L,
  group_check$n_patients[group_check$Group == "Normal"] == 12L,
  group_check$n_patients[group_check$Group == "Primary"] == 49L,
  group_check$n_patients[group_check$Group == "Metastatic CRPC"] == 27L
)

matched_normal_primary <- intersect(
  sample_info$Patient_ID[sample_info$Group == "Normal"],
  sample_info$Patient_ID[sample_info$Group == "Primary"]
)

matched_pair_check <- sample_info %>%
  dplyr::filter(
    Patient_ID %in% matched_normal_primary,
    Group %in% c("Normal", "Primary")
  ) %>%
  dplyr::count(Patient_ID, Group, name = "n") %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::arrange(Patient_ID)

stopifnot(
  length(matched_normal_primary) == 8L,
  nrow(matched_pair_check) == 8L,
  all(matched_pair_check$Normal == 1L),
  all(matched_pair_check$Primary == 1L)
)

normal_primary_id_check <- sample_info %>%
  dplyr::filter(Group %in% c("Normal", "Primary")) %>%
  dplyr::count(Patient_ID, name = "n_samples")

stopifnot(
  nrow(normal_primary_id_check) == 53L,
  sum(normal_primary_id_check$n_samples == 2L) == 8L,
  sum(normal_primary_id_check$n_samples == 1L) == 45L,
  all(normal_primary_id_check$n_samples %in% c(1L, 2L))
)

cat("Matched Normal-Primary patients:", length(matched_normal_primary), "\n")
print(as.data.frame(matched_pair_check), row.names = FALSE)
cat(
  "Normal-Primary analysis: 61 samples from 53 patients;",
  "8 paired patients and 45 unpaired patients\n"
)

utils::write.csv(
  sample_info,
  file.path(output_dir, "GSE35988_GPL6480_sample_metadata.csv"),
  row.names = FALSE
)

utils::write.csv(
  normal_primary_id_check,
  file.path(output_dir, "GSE35988_GPL6480_Normal_Primary_Patient_ID_QC.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Map GPL6480 probes to gene symbols
# -----------------------------------------------------------------------------

gpl6480 <- GEOquery::getGEO(
  "GPL6480",
  destdir = dirname(matrix_file)
)

gpl_annotation <- GEOquery::Table(gpl6480)

stopifnot(all(c("ID", "GENE_SYMBOL") %in% colnames(gpl_annotation)))

probe_map <- gpl_annotation %>%
  dplyr::transmute(
    Probe_ID = as.character(ID),
    Gene_symbol = toupper(trimws(as.character(GENE_SYMBOL)))
  ) %>%
  dplyr::filter(
    !is.na(Probe_ID),
    !is.na(Gene_symbol),
    Probe_ID != "",
    Gene_symbol != "",
    Gene_symbol != "---"
  ) %>%
  dplyr::distinct(Probe_ID, .keep_all = TRUE)

stopifnot(!anyDuplicated(probe_map$Probe_ID))

probe_symbol <- probe_map$Gene_symbol[
  match(rownames(expression_probe), probe_map$Probe_ID)
]

valid_symbol <- (
  !is.na(probe_symbol) &
    probe_symbol != "" &
    probe_symbol != "---" &
    !grepl("///|;|,", probe_symbol)
)

stopifnot(sum(valid_symbol) > 10000L)

expression_gene <- limma::avereps(
  expression_probe[valid_symbol, , drop = FALSE],
  ID = probe_symbol[valid_symbol]
)

storage.mode(expression_gene) <- "double"

stopifnot(
  !anyDuplicated(rownames(expression_gene)),
  identical(colnames(expression_gene), sample_info$Sample),
  all(is.finite(expression_gene)),
  "XPO1" %in% rownames(expression_gene)
)

cat(
  "Gene-level expression matrix:",
  nrow(expression_gene), "genes x",
  ncol(expression_gene), "samples\n"
)

xpo1_probe_annotation <- probe_map %>%
  dplyr::filter(
    Gene_symbol == "XPO1",
    Probe_ID %in% rownames(expression_probe)
  )

stopifnot(nrow(xpo1_probe_annotation) > 0L)

utils::write.csv(
  xpo1_probe_annotation,
  file.path(output_dir, "GSE35988_GPL6480_XPO1_probes_used.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Statistical testing
#
# Normal vs Primary: linear mixed-effects model using all 61 samples.
#   Fixed effect: Group
#   Random intercept: Patient_ID
# The 8 matched Normal-Primary pairs share Patient_ID; the other 45 samples
# have unique Patient_ID values. This uses all observations while accounting
# for within-patient dependence in the matched pairs.
#
# Primary vs Metastatic CRPC: unpaired Wilcoxon, 49 vs 27 patients.
# Normal vs Metastatic CRPC: unpaired Wilcoxon, 12 vs 27 patients.
# -----------------------------------------------------------------------------

run_endpoint_tests <- function(endpoint_data) {
  stopifnot(
    all(c("Sample", "Patient_ID", "Group", "Value") %in% colnames(endpoint_data)),
    nrow(endpoint_data) == 88L,
    !anyDuplicated(endpoint_data$Sample),
    all(is.finite(endpoint_data$Value))
  )
  
  normal_primary_data <- endpoint_data %>%
    dplyr::filter(Group %in% c("Normal", "Primary")) %>%
    dplyr::mutate(
      Group = factor(as.character(Group), levels = c("Normal", "Primary")),
      Patient_ID = factor(Patient_ID)
    ) %>%
    dplyr::arrange(Patient_ID, Group)
  
  normal_primary_cluster_check <- normal_primary_data %>%
    dplyr::count(Patient_ID, name = "n_samples")
  
  stopifnot(
    nrow(normal_primary_data) == 61L,
    dplyr::n_distinct(normal_primary_data$Patient_ID) == 53L,
    sum(normal_primary_cluster_check$n_samples == 2L) == 8L,
    sum(normal_primary_cluster_check$n_samples == 1L) == 45L,
    all(normal_primary_cluster_check$n_samples %in% c(1L, 2L))
  )
  
  matched_np <- normal_primary_data %>%
    dplyr::filter(
      Patient_ID %in% matched_normal_primary,
      Group %in% c("Normal", "Primary")
    ) %>%
    dplyr::select(Patient_ID, Group, Value) %>%
    tidyr::pivot_wider(names_from = Group, values_from = Value) %>%
    dplyr::arrange(Patient_ID)
  
  stopifnot(
    nrow(matched_np) == 8L,
    !anyNA(matched_np$Normal),
    !anyNA(matched_np$Primary)
  )
  
  normal_primary_model <- nlme::lme(
    fixed = Value ~ Group,
    random = ~1 | Patient_ID,
    data = normal_primary_data,
    method = "ML",
    na.action = stats::na.fail,
    control = nlme::lmeControl(
      opt = "optim",
      returnObject = FALSE
    )
  )
  
  normal_primary_null_model <- nlme::lme(
    fixed = Value ~ 1,
    random = ~1 | Patient_ID,
    data = normal_primary_data,
    method = "ML",
    na.action = stats::na.fail,
    control = nlme::lmeControl(
      opt = "optim",
      returnObject = FALSE
    )
  )
  
  normal_primary_lrt <- stats::anova(
    normal_primary_null_model,
    normal_primary_model
  )
  
  stopifnot(
    nrow(normal_primary_lrt) == 2L,
    all(c("L.Ratio", "p-value") %in% colnames(normal_primary_lrt))
  )
  
  normal_primary_lrt_statistic <- unname(
    normal_primary_lrt[2L, "L.Ratio"]
  )
  normal_primary_lrt_p <- unname(
    normal_primary_lrt[2L, "p-value"]
  )
  
  normal_primary_coefficients <- summary(normal_primary_model)$tTable
  
  stopifnot(
    "GroupPrimary" %in% rownames(normal_primary_coefficients),
    all(c("Value", "Std.Error", "DF", "t-value", "p-value") %in%
          colnames(normal_primary_coefficients))
  )
  
  normal_primary_result <- normal_primary_coefficients[
    "GroupPrimary",
    ,
    drop = FALSE
  ]
  
  normal_primary_variance <- nlme::VarCorr(normal_primary_model)
  patient_intercept_sd <- as.numeric(
    normal_primary_variance["(Intercept)", "StdDev"]
  )
  residual_sd <- as.numeric(
    normal_primary_variance["Residual", "StdDev"]
  )
  within_patient_correlation <- patient_intercept_sd^2 /
    (patient_intercept_sd^2 + residual_sd^2)
  
  primary_values <- endpoint_data$Value[endpoint_data$Group == "Primary"]
  normal_values <- endpoint_data$Value[endpoint_data$Group == "Normal"]
  metastatic_values <- endpoint_data$Value[
    endpoint_data$Group == "Metastatic CRPC"
  ]
  
  stopifnot(
    length(normal_values) == 12L,
    length(primary_values) == 49L,
    length(metastatic_values) == 27L
  )
  
  primary_metastatic_test <- stats::wilcox.test(
    x = metastatic_values,
    y = primary_values,
    paired = FALSE,
    alternative = "two.sided",
    exact = FALSE
  )
  
  normal_metastatic_test <- stats::wilcox.test(
    x = metastatic_values,
    y = normal_values,
    paired = FALSE,
    alternative = "two.sided",
    exact = FALSE
  )
  
  statistics <- data.frame(
    Comparison = c(
      "Primary vs Normal",
      "Metastatic CRPC vs Primary",
      "Metastatic CRPC vs Normal"
    ),
    Design = c(
      paste0(
        "Linear mixed-effects model (ML): Value ~ Group + ",
        "(1 | Patient_ID); likelihood-ratio test"
      ),
      "Two-sided unpaired Wilcoxon rank-sum test",
      "Two-sided unpaired Wilcoxon rank-sum test"
    ),
    Group_1 = c("Primary", "Metastatic CRPC", "Metastatic CRPC"),
    Group_1_n = c(
      length(primary_values),
      length(metastatic_values),
      length(metastatic_values)
    ),
    Group_1_n_patients = c(
      dplyr::n_distinct(
        endpoint_data$Patient_ID[endpoint_data$Group == "Primary"]
      ),
      dplyr::n_distinct(
        endpoint_data$Patient_ID[endpoint_data$Group == "Metastatic CRPC"]
      ),
      dplyr::n_distinct(
        endpoint_data$Patient_ID[endpoint_data$Group == "Metastatic CRPC"]
      )
    ),
    Group_2 = c("Normal", "Primary", "Normal"),
    Group_2_n = c(
      length(normal_values),
      length(primary_values),
      length(normal_values)
    ),
    Group_2_n_patients = c(
      dplyr::n_distinct(
        endpoint_data$Patient_ID[endpoint_data$Group == "Normal"]
      ),
      dplyr::n_distinct(
        endpoint_data$Patient_ID[endpoint_data$Group == "Primary"]
      ),
      dplyr::n_distinct(
        endpoint_data$Patient_ID[endpoint_data$Group == "Normal"]
      )
    ),
    Shared_patient_n = c(8L, 0L, 0L),
    Effect_estimate = c(
      unname(normal_primary_result[1, "Value"]),
      median(metastatic_values) - median(primary_values),
      median(metastatic_values) - median(normal_values)
    ),
    Effect_estimate_type = c(
      "Model-estimated mean difference",
      "Difference in group medians",
      "Difference in group medians"
    ),
    Median_difference = c(
      median(primary_values) - median(normal_values),
      median(metastatic_values) - median(primary_values),
      median(metastatic_values) - median(normal_values)
    ),
    Standard_error = c(
      unname(normal_primary_result[1, "Std.Error"]),
      NA_real_,
      NA_real_
    ),
    Degrees_of_freedom = c(
      1,
      NA_real_,
      NA_real_
    ),
    Statistic_name = c("Likelihood-ratio chi-square", "W", "W"),
    Statistic = c(
      normal_primary_lrt_statistic,
      unname(primary_metastatic_test$statistic),
      unname(normal_metastatic_test$statistic)
    ),
    P_value = c(
      normal_primary_lrt_p,
      primary_metastatic_test$p.value,
      normal_metastatic_test$p.value
    ),
    stringsAsFactors = FALSE
  )
  
  statistics$P_adjust_BH <- stats::p.adjust(
    statistics$P_value,
    method = "BH"
  )
  
  normal_primary_model_details <- data.frame(
    Comparison = "Primary vs Normal",
    n_samples = nrow(normal_primary_data),
    n_patients = dplyr::n_distinct(normal_primary_data$Patient_ID),
    n_paired_patients = sum(
      normal_primary_cluster_check$n_samples == 2L
    ),
    n_unpaired_samples = sum(
      normal_primary_cluster_check$n_samples == 1L
    ),
    Likelihood_ratio_statistic = normal_primary_lrt_statistic,
    Likelihood_ratio_df = 1L,
    P_value = normal_primary_lrt_p,
    Patient_random_intercept_SD = patient_intercept_sd,
    Residual_SD = residual_sd,
    Estimated_within_patient_correlation = within_patient_correlation,
    stringsAsFactors = FALSE
  )
  
  list(
    statistics = statistics,
    normal_primary_model_data = normal_primary_data %>%
      dplyr::mutate(
        Group = as.character(Group),
        Patient_ID = as.character(Patient_ID)
      ),
    matched_normal_primary = matched_np,
    normal_primary_model = normal_primary_model,
    normal_primary_null_model = normal_primary_null_model,
    normal_primary_model_details = normal_primary_model_details,
    patient_intercept_sd = patient_intercept_sd,
    residual_sd = residual_sd,
    within_patient_correlation = within_patient_correlation
  )
}

# -----------------------------------------------------------------------------
# 6. Plotting function; statistics are not displayed on the figures
# -----------------------------------------------------------------------------

plot_endpoint <- function(endpoint_data, y_label) {
  ggplot2::ggplot(
    endpoint_data,
    ggplot2::aes(x = Group, y = Value, fill = Group)
  ) +
    ggplot2::geom_violin(
      trim = FALSE,
      scale = "width",
      width = 0.8,
      adjust = 1.3,
      color = NA,
      alpha = 0.55
    ) +
    ggplot2::geom_boxplot(
      width = 0.30,
      outlier.shape = NA,
      color = "black",
      alpha = 0.75,
      linewidth = 0.6
    ) +
    ggplot2::geom_jitter(
      width = 0.14,
      size = 1.5,
      alpha = 0.70,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Normal" = "#91D1C2FF",
        "Primary" = "#4DBBD5FF",
        "Metastatic CRPC" = "#E64B35FF"
      )
    )  +
    ggplot2::labs(x = NULL, y = y_label) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 13),
      legend.position = "none"
    )
}


}))

# =============================================================================
# Fig6E_F
# =============================================================================
nmd3_step("Fig6E_F", c("Prepare_GSE35988"), quote({
# Fig6E_F — analysis below
# PEG10 and ENO2 expression
# Primary vs Metastatic CRPC
#
# Required upstream objects:
#   expression_gene
#   sample_info
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

set.seed(123)


# -----------------------------------------------------------------------------
# 1. Output directory
# -----------------------------------------------------------------------------

output_dir_6EF <- paste0(
  "results/",
  "GSE35988_GPL6480_Figure6EF"
)

dir.create(
  output_dir_6EF,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  dir.exists(output_dir_6EF),
  exists("expression_gene"),
  exists("sample_info")
)


# ============================================================
# Fig6E_F — analysis below
# PEG10 and ENO2 expression
#
# Reuse the validated GSE35988 metadata generated in Figure 1
# ============================================================

genes_6EF <- c(
  "PEG10",
  "ENO2"
)

stopifnot(
  exists("expression_gene"),
  exists("sample_info"),
  all(
    genes_6EF %in%
      rownames(expression_gene)
  )
)

# ------------------------------------------------------------
# 1. Reuse Figure 1 GSE35988 metadata
# ------------------------------------------------------------

sample_info_6EF <- sample_info

stopifnot(
  all(
    c(
      "Sample",
      "Group",
      "Patient_ID"
    ) %in%
      colnames(sample_info_6EF)
  )
)

# Sample_title is only descriptive; create it if Figure 1
# metadata does not contain this column
if (!"Sample_title" %in% colnames(sample_info_6EF)) {
  sample_info_6EF$Sample_title <-
    as.character(
      sample_info_6EF$Sample
    )
}

sample_info_6EF <- sample_info_6EF %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Patient_ID = as.character(Patient_ID),
    Group = as.character(Group)
  )

# ------------------------------------------------------------
# 2. Confirm that Figure 1 metadata matches expression matrix
# ------------------------------------------------------------

expression_samples_6EF <- as.character(
  colnames(expression_gene)
)

stopifnot(
  setequal(
    expression_samples_6EF,
    sample_info_6EF$Sample
  ),
  !anyDuplicated(
    sample_info_6EF$Sample
  ),
  !anyNA(
    sample_info_6EF$Patient_ID
  ),
  !anyNA(
    sample_info_6EF$Group
  )
)

# Align metadata exactly to expression matrix
sample_info_6EF <- sample_info_6EF[
  match(
    expression_samples_6EF,
    sample_info_6EF$Sample
  ),
  ,
  drop = FALSE
]

rownames(sample_info_6EF) <- NULL

stopifnot(
  identical(
    expression_samples_6EF,
    sample_info_6EF$Sample
  )
)

# ------------------------------------------------------------
# 3. Confirm GSE35988 group sizes
# ------------------------------------------------------------

group_check_6EF_upstream <- sample_info_6EF %>%
  dplyr::filter(
    Group %in%
      c(
        "Normal",
        "Primary",
        "Metastatic CRPC"
      )
  ) %>%
  dplyr::group_by(
    Group
  ) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_patients = dplyr::n_distinct(
      Patient_ID
    ),
    .groups = "drop"
  )

print(
  as.data.frame(
    group_check_6EF_upstream
  ),
  row.names = FALSE
)

stopifnot(
  group_check_6EF_upstream$n_samples[
    group_check_6EF_upstream$Group ==
      "Normal"
  ] == 12L,
  
  group_check_6EF_upstream$n_samples[
    group_check_6EF_upstream$Group ==
      "Primary"
  ] == 49L,
  
  group_check_6EF_upstream$n_samples[
    group_check_6EF_upstream$Group ==
      "Metastatic CRPC"
  ] == 27L,
  
  group_check_6EF_upstream$n_patients[
    group_check_6EF_upstream$Group ==
      "Normal"
  ] == 12L,
  
  group_check_6EF_upstream$n_patients[
    group_check_6EF_upstream$Group ==
      "Primary"
  ] == 49L,
  
  group_check_6EF_upstream$n_patients[
    group_check_6EF_upstream$Group ==
      "Metastatic CRPC"
  ] == 27L
)

# ------------------------------------------------------------
# 4. Confirm the eight matched Normal-Primary patients
# ------------------------------------------------------------

normal_patient_ids_6EF <- sample_info_6EF %>%
  dplyr::filter(
    Group == "Normal"
  ) %>%
  dplyr::pull(
    Patient_ID
  ) %>%
  unique()

primary_patient_ids_all_6EF <- sample_info_6EF %>%
  dplyr::filter(
    Group == "Primary"
  ) %>%
  dplyr::pull(
    Patient_ID
  ) %>%
  unique()

shared_patient_ids_6EF <- intersect(
  normal_patient_ids_6EF,
  primary_patient_ids_all_6EF
)

cat(
  "\nMatched Normal-Primary patients:",
  length(shared_patient_ids_6EF),
  "\n"
)

stopifnot(
  length(
    shared_patient_ids_6EF
  ) == 8L
)

# ------------------------------------------------------------
# 5. Extract PEG10 and ENO2 expression
# ------------------------------------------------------------

expression_6EF <- expression_gene[
  genes_6EF,
  ,
  drop = FALSE
]

storage.mode(
  expression_6EF
) <- "double"

stopifnot(
  all(
    is.finite(
      expression_6EF
    )
  )
)

df_6EF <- data.frame(
  Sample = colnames(
    expression_6EF
  ),
  
  PEG10 = as.numeric(
    expression_6EF[
      "PEG10",
      ,
      drop = TRUE
    ]
  ),
  
  ENO2 = as.numeric(
    expression_6EF[
      "ENO2",
      ,
      drop = TRUE
    ]
  ),
  
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    sample_info_6EF %>%
      dplyr::select(
        Sample,
        Sample_title,
        Patient_ID,
        Group
      ),
    by = "Sample"
  ) %>%
  dplyr::filter(
    Group %in%
      c(
        "Primary",
        "Metastatic CRPC"
      )
  ) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c(
        "Primary",
        "Metastatic CRPC"
      )
    )
  )

# ------------------------------------------------------------
# 6. Final checks for Figure 6E-F
# ------------------------------------------------------------

group_check_6EF <- df_6EF %>%
  dplyr::group_by(
    Group
  ) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_patients = dplyr::n_distinct(
      Patient_ID
    ),
    .groups = "drop"
  )

print(
  as.data.frame(
    group_check_6EF
  ),
  row.names = FALSE
)

stopifnot(
  nrow(df_6EF) == 76L,
  
  group_check_6EF$n_samples[
    group_check_6EF$Group ==
      "Primary"
  ] == 49L,
  
  group_check_6EF$n_samples[
    group_check_6EF$Group ==
      "Metastatic CRPC"
  ] == 27L,
  
  group_check_6EF$n_patients[
    group_check_6EF$Group ==
      "Primary"
  ] == 49L,
  
  group_check_6EF$n_patients[
    group_check_6EF$Group ==
      "Metastatic CRPC"
  ] == 27L,
  
  all(
    is.finite(
      df_6EF$PEG10
    )
  ),
  
  all(
    is.finite(
      df_6EF$ENO2
    )
  )
)

cat(
  "\nPrimary vs Metastatic CRPC:",
  "76 samples from 76 independent patients\n"
)



# -----------------------------------------------------------------------------
# 5. Run Primary vs Metastatic CRPC test
# -----------------------------------------------------------------------------

run_gene_test_6EF <- function(
    data,
    gene_name
) {
  
  primary_values <- data[
    data$Group == "Primary",
    gene_name,
    drop = TRUE
  ]
  
  metastatic_values <- data[
    data$Group == "Metastatic CRPC",
    gene_name,
    drop = TRUE
  ]
  
  stopifnot(
    length(primary_values) == 49L,
    length(metastatic_values) == 27L,
    all(is.finite(primary_values)),
    all(is.finite(metastatic_values))
  )
  
  test_result <- stats::wilcox.test(
    x = metastatic_values,
    y = primary_values,
    paired = FALSE,
    alternative = "two.sided",
    exact = FALSE
  )
  
  data.frame(
    Gene = gene_name,
    
    Comparison =
      "Metastatic CRPC vs Primary",
    
    Design =
      "Two-sided unpaired Wilcoxon rank-sum test",
    
    Group_1 =
      "Metastatic CRPC",
    
    Group_1_n =
      length(metastatic_values),
    
    Group_1_n_patients =
      dplyr::n_distinct(
        data$Patient_ID[
          data$Group == "Metastatic CRPC"
        ]
      ),
    
    Group_1_median =
      stats::median(
        metastatic_values
      ),
    
    Group_2 =
      "Primary",
    
    Group_2_n =
      length(primary_values),
    
    Group_2_n_patients =
      dplyr::n_distinct(
        data$Patient_ID[
          data$Group == "Primary"
        ]
      ),
    
    Group_2_median =
      stats::median(
        primary_values
      ),
    
    Median_difference =
      stats::median(
        metastatic_values
      ) -
      stats::median(
        primary_values
      ),
    
    Statistic_name =
      "Wilcoxon W",
    
    Statistic =
      unname(
        test_result$statistic
      ),
    
    P_value =
      test_result$p.value,
    
    stringsAsFactors = FALSE
  )
}


statistics_6EF <- dplyr::bind_rows(
  run_gene_test_6EF(
    df_6EF,
    "PEG10"
  ),
  run_gene_test_6EF(
    df_6EF,
    "ENO2"
  )
)

statistics_6EF$P_adjust_BH <- stats::p.adjust(
  statistics_6EF$P_value,
  method = "BH"
)

cat(
  "\nFigure 6E-F statistics:\n"
)

print(
  statistics_6EF,
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 6. Save statistics and source data
# -----------------------------------------------------------------------------

utils::write.csv(
  statistics_6EF,
  file.path(
    output_dir_6EF,
    "Figure6EF_GSE35988_PEG10_ENO2_statistics.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  df_6EF,
  file.path(
    output_dir_6EF,
    "Figure6EF_GSE35988_PEG10_ENO2_source_data.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 7. Plotting function
# -----------------------------------------------------------------------------

plot_gene_6EF <- function(
    data,
    gene_name,
    y_label
) {
  
  plot_data <- data %>%
    dplyr::transmute(
      Sample = Sample,
      Sample_title = Sample_title,
      Patient_ID = Patient_ID,
      Group = Group,
      Value = .data[[gene_name]]
    )
  
  stopifnot(
    nrow(plot_data) == 76L,
    all(is.finite(plot_data$Value))
  )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Group,
      y = Value,
      fill = Group
    )
  ) +
    
    ggplot2::geom_violin(
      trim = FALSE,
      scale = "width",
      width = 0.8,
      adjust = 1.3,
      color = NA,
      alpha = 0.55
    ) +
    
    ggplot2::geom_boxplot(
      width = 0.30,
      outlier.shape = NA,
      color = "black",
      alpha = 0.75,
      linewidth = 0.6
    ) +
    
    ggplot2::geom_jitter(
      width = 0.14,
      size = 1.5,
      alpha = 0.70,
      color = "black"
    ) +
    
    ggplot2::scale_fill_manual(
      values = c(
        "Primary" = "#4DBBD5FF",
        "Metastatic CRPC" = "#E64B35FF"
      )
    ) +
    
    
    ggplot2::scale_x_discrete(
      labels = c(
        "Primary" = "Primary",
        "Metastatic CRPC" = "Metastatic"
      )
    ) +
    
    ggplot2::labs(
      x = NULL,
      y = y_label
    ) +
    
    ggplot2::theme_classic(
      base_size = 14
    ) +
    
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold"
      ),
      
      axis.text.y = ggplot2::element_text(
        size = 12
      ),
      
      axis.title.y = ggplot2::element_text(
        size = 13
      ),
      
      legend.position = "none"
    )
}


# =============================================================================
# Fig6E_F — analysis below
# PEG10 expression: Primary vs Metastatic CRPC
# =============================================================================

p_6E <- plot_gene_6EF(
  data = df_6EF,
  gene_name = "PEG10",
  y_label = "PEG10 mRNA expression"
)


nmd3_save_pdf(
  filename = file.path(
    output_dir_6EF,
    "Fig6E_GSE35988_PEG10_Primary_vs_Metastatic.pdf"
  ),
  plot = p_6E,
  width = 5,
  height = 5,
  device = cairo_pdf
)


# =============================================================================
# Fig6E_F — analysis below
# ENO2 expression: Primary vs Metastatic CRPC
# =============================================================================

p_6F <- plot_gene_6EF(
  data = df_6EF,
  gene_name = "ENO2",
  y_label = "ENO2 mRNA expression"
)


nmd3_save_pdf(
  filename = file.path(
    output_dir_6EF,
    "Fig6F_GSE35988_ENO2_Primary_vs_Metastatic.pdf"
  ),
  plot = p_6F,
  width = 5,
  height = 5,
  device = cairo_pdf
)

cat(
  "\nFigure 6E-F completed.\n",
  "Output directory: ",
  output_dir_6EF,
  "\n",
  sep = ""
)




}))

# =============================================================================
# Fig6G_H
# =============================================================================
nmd3_step("Fig6G_H", c("Fig6D"), quote({
# Fig6G_H — analysis below
# Patient-level Spearman correlations only
#
# Fig6G_H — analysis below
#   Ribosome-export score vs neuroendocrine score
#
# Fig6G_H — analysis below
#   Ribosome-export score vs PEG10 expression
#
# Main statistical unit:
#   One value per patient. Multiple tumors from the same patient
#   are averaged before correlation analysis (15 tumors, 10 patients).
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# 1. Reuse the WCM data loaded for Figure 6A
# ============================================================

stopifnot(
  exists("expr_wcm"),
  exists("clinical_wcm")
)

# Use the complete tables, not mat_wcm from Figure 6A.
# mat_wcm contains only PEG10 and ENO2 and therefore cannot be used
# to calculate the ribosome-export and neuroendocrine signatures.

expression_raw_6CD <- as.data.frame(
  expr_wcm,
  check.names = FALSE
)

clinical_raw_6CD <- as.data.frame(
  clinical_wcm,
  check.names = FALSE
)

output_dir_6CD <- file.path(
  getwd(),
  "Figure6CD_WCM_patient_correlation"
)

dir.create(
  output_dir_6CD,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Prespecified signatures
# ============================================================

ribosome_export_genes_6CD <- c(
  "EIF6", "LSG1", "LTV1", "MDN1", "NMD3",
  "NOP9", "NPM1", "NUP88", "RAN", "RBM10",
  "RIOK2", "RPS15", "SDAD1", "XPO1"
)

neuroendocrine_genes_6CD <- c(
  "SOGA3", "RGS7", "JAKMIP2", "PROX1", "SYT11",
  "CAND2", "TRIM9", "SEC11C", "EZH2", "ETV5",
  "AURKA", "DNMT1", "MAP10", "SVOP", "ST8SIA3",
  "MYCN", "NRSN1", "KCNB2", "BRINP1", "SCG3",
  "SEZ6", "KIAA0408", "LRRC16B", "GPX2", "PCSK1",
  "KCND2", "ASXL3", "GNAO1", "C7ORF76"
)

stopifnot(
  length(ribosome_export_genes_6CD) == 14L,
  length(neuroendocrine_genes_6CD) == 29L,
  !anyDuplicated(ribosome_export_genes_6CD),
  !anyDuplicated(neuroendocrine_genes_6CD)
)

required_genes_6CD <- unique(c(
  ribosome_export_genes_6CD,
  neuroendocrine_genes_6CD,
  "PEG10"
))

# ============================================================
# 3. Check the reused expression and clinical annotation
# ============================================================

stopifnot(
  "Hugo_Symbol" %in% colnames(expression_raw_6CD),
  all(
    c("SAMPLE_ID", "PATIENT_ID", "DISEASE_CODE") %in%
      colnames(clinical_raw_6CD)
  )
)

clinical_map_6CD <- clinical_raw_6CD %>%
  dplyr::transmute(
    Sample_ID = trimws(as.character(SAMPLE_ID)),
    Patient_ID = trimws(as.character(PATIENT_ID)),
    Group = trimws(as.character(DISEASE_CODE))
  ) %>%
  dplyr::filter(
    Group %in% c("CRPC-Adeno", "CRPC-NE"),
    !is.na(Sample_ID),
    Sample_ID != "",
    !is.na(Patient_ID),
    Patient_ID != ""
  ) %>%
  dplyr::distinct(
    Sample_ID,
    .keep_all = TRUE
  )

common_samples_6CD <- clinical_map_6CD$Sample_ID[
  clinical_map_6CD$Sample_ID %in%
    colnames(expression_raw_6CD)
]

clinical_map_6CD <- clinical_map_6CD %>%
  dplyr::filter(
    Sample_ID %in% common_samples_6CD
  ) %>%
  dplyr::arrange(
    match(Sample_ID, common_samples_6CD)
  )

stopifnot(
  nrow(clinical_map_6CD) == 49L,
  sum(clinical_map_6CD$Group == "CRPC-Adeno") == 34L,
  sum(clinical_map_6CD$Group == "CRPC-NE") == 15L,
  !anyDuplicated(clinical_map_6CD$Sample_ID)
)

# ============================================================
# 4. Prepare log2(FPKM + 1) expression
# ============================================================

expression_selected_6CD <- expression_raw_6CD %>%
  dplyr::transmute(
    Gene = toupper(
      trimws(
        as.character(Hugo_Symbol)
      )
    ),
    dplyr::across(
      dplyr::all_of(common_samples_6CD),
      ~ suppressWarnings(
        as.numeric(as.character(.x))
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != "",
    Gene %in% required_genes_6CD
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(common_samples_6CD),
      ~ if (all(is.na(.x))) {
        NA_real_
      } else {
        mean(.x, na.rm = TRUE)
      }
    ),
    .groups = "drop"
  )

expression_matrix_6CD <- expression_selected_6CD %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

storage.mode(expression_matrix_6CD) <- "double"

expression_matrix_6CD <- expression_matrix_6CD[
  ,
  common_samples_6CD,
  drop = FALSE
]

# FPKM cannot be negative. Record negative numerical values for QC
# and set them to zero before log2(FPKM + 1).

negative_index_6CD <- which(
  expression_matrix_6CD < 0,
  arr.ind = TRUE
)

if (nrow(negative_index_6CD) > 0L) {
  negative_fpkm_qc_6CD <- tibble::tibble(
    Gene = rownames(expression_matrix_6CD)[
      negative_index_6CD[, 1]
    ],
    Sample_ID = colnames(expression_matrix_6CD)[
      negative_index_6CD[, 2]
    ],
    Original_value = expression_matrix_6CD[
      negative_index_6CD
    ],
    Corrected_value = 0
  ) %>%
    dplyr::arrange(Original_value)
} else {
  negative_fpkm_qc_6CD <- tibble::tibble(
    Gene = character(),
    Sample_ID = character(),
    Original_value = numeric(),
    Corrected_value = numeric()
  )
}

cat(
  "Negative FPKM values corrected to zero:",
  nrow(negative_fpkm_qc_6CD),
  "of",
  length(expression_matrix_6CD),
  "values\n"
)

expression_fpkm_6CD <- expression_matrix_6CD

expression_fpkm_6CD[
  expression_fpkm_6CD < 0
] <- 0

log_expression_6CD <- log2(
  expression_fpkm_6CD + 1
)

stopifnot(
  all(is.finite(log_expression_6CD))
)

# ============================================================
# 5. Gene-wise z-scores across all 49 tumors
# ============================================================

gene_means_6CD <- rowMeans(
  log_expression_6CD,
  na.rm = TRUE
)

gene_sds_6CD <- apply(
  log_expression_6CD,
  1,
  stats::sd,
  na.rm = TRUE
)

variable_z_genes_6CD <- names(gene_sds_6CD)[
  is.finite(gene_means_6CD) &
    is.finite(gene_sds_6CD) &
    gene_sds_6CD > 0
]

zero_variance_genes_6CD <- setdiff(
  rownames(log_expression_6CD),
  variable_z_genes_6CD
)

z_expression_6CD <- sweep(
  log_expression_6CD[
    variable_z_genes_6CD,
    ,
    drop = FALSE
  ],
  1,
  gene_means_6CD[variable_z_genes_6CD],
  FUN = "-"
)

z_expression_6CD <- sweep(
  z_expression_6CD,
  1,
  gene_sds_6CD[variable_z_genes_6CD],
  FUN = "/"
)

present_ribosome_genes_6CD <- intersect(
  ribosome_export_genes_6CD,
  rownames(z_expression_6CD)
)

present_ne_genes_6CD <- intersect(
  neuroendocrine_genes_6CD,
  rownames(z_expression_6CD)
)

cat(
  "Ribosome-export genes used:",
  length(present_ribosome_genes_6CD),
  "/ 14\n"
)

cat(
  "Neuroendocrine genes used:",
  length(present_ne_genes_6CD),
  "/ 29\n"
)

cat(
  "Missing ribosome-export genes:",
  paste(
    setdiff(
      ribosome_export_genes_6CD,
      rownames(log_expression_6CD)
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Missing neuroendocrine genes:",
  paste(
    setdiff(
      neuroendocrine_genes_6CD,
      rownames(log_expression_6CD)
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Zero-variance genes excluded:",
  paste(zero_variance_genes_6CD, collapse = ", "),
  "\n"
)

stopifnot(
  length(present_ribosome_genes_6CD) > 0L,
  length(present_ne_genes_6CD) > 0L,
  "PEG10" %in% rownames(log_expression_6CD),
  all(
    is.finite(
      z_expression_6CD[
        present_ribosome_genes_6CD,
        ,
        drop = FALSE
      ]
    )
  ),
  all(
    is.finite(
      z_expression_6CD[
        present_ne_genes_6CD,
        ,
        drop = FALSE
      ]
    )
  )
)

# ============================================================
# 6. Calculate sample-level scores
# ============================================================

ribosome_score_6CD <- colMeans(
  z_expression_6CD[
    present_ribosome_genes_6CD,
    ,
    drop = FALSE
  ]
)

ne_score_6CD <- colMeans(
  z_expression_6CD[
    present_ne_genes_6CD,
    ,
    drop = FALSE
  ]
)

sample_scores_all_6CD <- clinical_map_6CD %>%
  dplyr::mutate(
    Ribosome_Export_Score = as.numeric(
      ribosome_score_6CD[Sample_ID]
    ),
    Neuroendocrine_Score = as.numeric(
      ne_score_6CD[Sample_ID]
    ),
    PEG10_log2_FPKM_plus1 = as.numeric(
      log_expression_6CD[
        "PEG10",
        Sample_ID
      ]
    )
  )

sample_scores_ne_6CD <- sample_scores_all_6CD %>%
  dplyr::filter(
    Group == "CRPC-NE"
  ) %>%
  dplyr::arrange(
    Patient_ID,
    Sample_ID
  )

patient_sample_counts_6CD <- sample_scores_ne_6CD %>%
  dplyr::count(
    Patient_ID,
    name = "Number_of_samples"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Number_of_samples),
    Patient_ID
  )

cat("\nCRPC-NE patient/sample structure:\n")
print(patient_sample_counts_6CD)

stopifnot(
  nrow(sample_scores_ne_6CD) == 15L,
  dplyr::n_distinct(
    sample_scores_ne_6CD$Patient_ID
  ) == 10L,
  sum(patient_sample_counts_6CD$Number_of_samples) -
    nrow(patient_sample_counts_6CD) == 5L,
  all(
    is.finite(
      sample_scores_ne_6CD$Ribosome_Export_Score
    )
  ),
  all(
    is.finite(
      sample_scores_ne_6CD$Neuroendocrine_Score
    )
  ),
  all(
    is.finite(
      sample_scores_ne_6CD$PEG10_log2_FPKM_plus1
    )
  )
)

# ============================================================
# 7. Aggregate to one value per patient
# ============================================================

patient_mean_data_6CD <- sample_scores_ne_6CD %>%
  dplyr::group_by(Patient_ID) %>%
  dplyr::summarise(
    Ribosome_Export_Score = mean(
      Ribosome_Export_Score,
      na.rm = TRUE
    ),
    Neuroendocrine_Score = mean(
      Neuroendocrine_Score,
      na.rm = TRUE
    ),
    PEG10_log2_FPKM_plus1 = mean(
      PEG10_log2_FPKM_plus1,
      na.rm = TRUE
    ),
    Number_of_samples = dplyr::n(),
    Samples = paste(
      Sample_ID,
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Patient_ID)

stopifnot(
  nrow(patient_mean_data_6CD) == 10L,
  !anyDuplicated(patient_mean_data_6CD$Patient_ID)
)

cat("\nPatient-level mean data used for Figure 6C-D:\n")
print(patient_mean_data_6CD)

# ============================================================
# 8. Patient-level Spearman correlations
# ============================================================

cor_6C <- stats::cor.test(
  patient_mean_data_6CD$Ribosome_Export_Score,
  patient_mean_data_6CD$Neuroendocrine_Score,
  method = "spearman",
  exact = FALSE
)

cor_6D <- stats::cor.test(
  patient_mean_data_6CD$Ribosome_Export_Score,
  patient_mean_data_6CD$PEG10_log2_FPKM_plus1,
  method = "spearman",
  exact = FALSE
)

cat("\nFigure 6C Spearman correlation:\n")
print(cor_6C)

cat("\nFigure 6D Spearman correlation:\n")
print(cor_6D)

correlation_statistics_6CD <- tibble::tibble(
  Figure = c("Figure 6C", "Figure 6D"),
  X_variable = c(
    "Ribosome Export Score",
    "Ribosome Export Score"
  ),
  Y_variable = c(
    "Neuroendocrine Score",
    "PEG10 log2(FPKM + 1)"
  ),
  Statistical_unit = "Patient",
  n_patients = nrow(patient_mean_data_6CD),
  Spearman_rho = c(
    unname(cor_6C$estimate),
    unname(cor_6D$estimate)
  ),
  P_value = c(
    cor_6C$p.value,
    cor_6D$p.value
  )
)

print(correlation_statistics_6CD)

# ============================================================
# 9. Figure 6C-D
# ============================================================

make_correlation_plot_6CD <- function(
    data,
    endpoint,
    y_axis_title,
    panel_color,
    correlation_text_size
) {
  x_values <- data$Ribosome_Export_Score
  y_values <- data[[endpoint]]
  
  x_range <- range(x_values, finite = TRUE)
  y_range <- range(y_values, finite = TRUE)
  
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  
  if (!is.finite(x_span) || x_span == 0) {
    x_span <- 1
  }
  
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- 1
  }
  
  correlation_label_x <- x_range[1] +
    0.03 * x_span
  
  correlation_label_y <- y_range[2] +
    0.18 * y_span
  
  ggplot(
    data,
    aes(
      x = Ribosome_Export_Score,
      y = .data[[endpoint]]
    )
  ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      fill = panel_color,
      color = "grey30",
      alpha = 0.15,
      linewidth = 0.8
    ) +
    geom_point(
      fill = panel_color,
      shape = 21,
      size = 3.5,
      color = "white",
      stroke = 0.6,
      alpha = 0.9
    ) +
    ggpubr::stat_cor(
      method = "spearman",
      label.x = correlation_label_x,
      label.y = correlation_label_y,
      size = correlation_text_size,
      fontface = "bold",
      color = panel_color,
      hjust = 0
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.05, 0.10)
      )
    ) +
    labs(
      x = "Ribosome Export Score",
      y = y_axis_title
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.2
      ),
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      axis.ticks.length = grid::unit(
        0.2,
        "cm"
      ),
      axis.text = element_text(
        color = "black",
        size = 12,
        face = "bold"
      ),
      axis.title = element_text(
        color = "black",
        size = 14,
        face = "bold"
      ),
      legend.position = "none"
    )
}

p_6C <- make_correlation_plot_6CD(
  data = patient_mean_data_6CD,
  endpoint = "Neuroendocrine_Score",
  y_axis_title = "Neuroendocrine Score",
  panel_color = "#4DBBD5FF",
  correlation_text_size = 6
)

p_6D <- make_correlation_plot_6CD(
  data = patient_mean_data_6CD,
  endpoint = "PEG10_log2_FPKM_plus1",
  y_axis_title = "PEG10 mRNA expression",
  panel_color = "#E64B35FF",
  correlation_text_size = 7
)

p_6CD <- p_6C + p_6D +
  patchwork::plot_layout(ncol = 2)




# ============================================================
# 10. Source Data export
# ============================================================

gene_set_qc_6CD <- dplyr::bind_rows(
  tibble::tibble(
    Signature = "Ribosome export",
    Gene = ribosome_export_genes_6CD,
    Measured = ribosome_export_genes_6CD %in%
      rownames(log_expression_6CD),
    Used = ribosome_export_genes_6CD %in%
      present_ribosome_genes_6CD
  ),
  tibble::tibble(
    Signature = "Neuroendocrine",
    Gene = neuroendocrine_genes_6CD,
    Measured = neuroendocrine_genes_6CD %in%
      rownames(log_expression_6CD),
    Used = neuroendocrine_genes_6CD %in%
      present_ne_genes_6CD
  )
)

write.csv(
  negative_fpkm_qc_6CD,
  file.path(
    output_dir_6CD,
    "00_negative_FPKM_values_corrected_to_zero.csv"
  ),
  row.names = FALSE
)

write.csv(
  gene_set_qc_6CD,
  file.path(
    output_dir_6CD,
    "01_gene_sets_used.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_scores_all_6CD,
  file.path(
    output_dir_6CD,
    "02_all_49_tumor_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_scores_ne_6CD,
  file.path(
    output_dir_6CD,
    "03_CRPC_NE_15_tumor_source_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  patient_sample_counts_6CD,
  file.path(
    output_dir_6CD,
    "04_CRPC_NE_patient_sample_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  patient_mean_data_6CD,
  file.path(
    output_dir_6CD,
    "05_patient_mean_source_data_n10.csv"
  ),
  row.names = FALSE
)

write.csv(
  correlation_statistics_6CD,
  file.path(
    output_dir_6CD,
    "06_patient_level_Spearman_statistics.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 11. Save figures
# ============================================================

nmd3_save_pdf(
  filename = file.path(
    output_dir_6CD,
    "Fig6G_level_NE_score.pdf"
  ),
  plot = p_6C,
  width = 5,
  height = 4.8,
  device = cairo_pdf
)

nmd3_save_pdf(
  filename = file.path(
    output_dir_6CD,
    "Fig6H_level_PEG10.pdf"
  ),
  plot = p_6D,
  width = 5,
  height = 4.8,
  device = cairo_pdf
)

cat(
  "\nFigure 6C-D patient-level correlations completed.\n",
  "Main statistical unit: 10 independent patients.\n",
  "All outputs were saved to: ",
  output_dir_6CD,
  "\n",
  sep = ""
)



}))

# =============================================================================
# Fig6I
# =============================================================================
nmd3_step("Fig6I", c("Fig6A"), quote({
# Fig6I — analysis below
# Three-group KPNB1 and CCT5 protein abundance
# Normal / Primary / Metastatic
#
# Reuses the public-proteomics objects generated in Figures 2K-2L:
#   all_protein_matrix
#   sample_metadata
#   pair_design
#   group_colors
#   publication_theme
#   save_pdf
#   output_dir
# =============================================================================


# -----------------------------------------------------------------------------
# 18. Check upstream objects
# -----------------------------------------------------------------------------

# Fig6A was completed earlier in this same main-figure run.
# Fetch its objects explicitly; vapply(..., exists) in a nested call can
# check the wrong calling frame and falsely report missing objects.
fig6a_env <- .nmd3_state$modules[["Fig6A"]]
if (!is.environment(fig6a_env)) {
  stop("Fig6I needs a completed Fig6A module in this run")
}
required_upstream_objects <- c(
  "all_protein_matrix", "sample_metadata", "pair_design",
  "group_colors", "publication_theme", "save_pdf", "output_dir"
)
missing_upstream_objects <- required_upstream_objects[
  !vapply(required_upstream_objects, function(x) {
    base::exists(x, envir = fig6a_env, inherits = FALSE)
  }, logical(1))
]
if (length(missing_upstream_objects)) {
  stop("Fig6A did not create: ", paste(missing_upstream_objects, collapse = ", "))
}
# Resolve each object against the actual Fig6A module, not Global Environment.
all_protein_matrix <- base::get("all_protein_matrix", envir = fig6a_env, inherits = FALSE)
sample_metadata   <- base::get("sample_metadata",   envir = fig6a_env, inherits = FALSE)
pair_design       <- base::get("pair_design",       envir = fig6a_env, inherits = FALSE)
group_colors      <- base::get("group_colors",      envir = fig6a_env, inherits = FALSE)
publication_theme <- base::get("publication_theme", envir = fig6a_env, inherits = FALSE)
save_pdf          <- base::get("save_pdf",          envir = fig6a_env, inherits = FALSE)
output_dir        <- base::get("output_dir",        envir = fig6a_env, inherits = FALSE)

genes_KPNB1_CCT5 <- c(
  "KPNB1",
  "CCT5"
)

stopifnot(
  all(
    genes_KPNB1_CCT5 %in%
      rownames(all_protein_matrix)
  ),
  identical(
    colnames(all_protein_matrix),
    sample_metadata$Sample
  )
)


# -----------------------------------------------------------------------------
# 19. Extract KPNB1 and CCT5 protein abundance
# -----------------------------------------------------------------------------

KPNB1_CCT5_matrix <- all_protein_matrix[
  genes_KPNB1_CCT5,
  sample_metadata$Sample,
  drop = FALSE
]

storage.mode(KPNB1_CCT5_matrix) <- "double"

stopifnot(
  identical(
    rownames(KPNB1_CCT5_matrix),
    genes_KPNB1_CCT5
  ),
  identical(
    colnames(KPNB1_CCT5_matrix),
    sample_metadata$Sample
  )
)


KPNB1_CCT5_all_data <- as.data.frame(
  t(KPNB1_CCT5_matrix),
  stringsAsFactors = FALSE
) %>%
  tibble::rownames_to_column(
    "Sample"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      genes_KPNB1_CCT5
    ),
    names_to = "Gene",
    values_to = "Value"
  ) %>%
  dplyr::left_join(
    sample_metadata,
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c(
        "Normal",
        "Primary",
        "Metastatic"
      )
    ),
    
    Measure = dplyr::recode(
      Gene,
      "KPNB1" = "KPNB1 protein",
      "CCT5" = "CCT5 protein"
    ),
    
    Measure = factor(
      Measure,
      levels = c(
        "KPNB1 protein",
        "CCT5 protein"
      )
    )
  )

stopifnot(
  nrow(KPNB1_CCT5_all_data) == 116L,
  !anyNA(KPNB1_CCT5_all_data$Group),
  !anyNA(KPNB1_CCT5_all_data$Measure)
)


# -----------------------------------------------------------------------------
# 20. Audit missing protein values
# -----------------------------------------------------------------------------

KPNB1_CCT5_missing <- KPNB1_CCT5_all_data %>%
  dplyr::filter(
    !is.finite(Value)
  ) %>%
  dplyr::select(
    Gene,
    Measure,
    Sample,
    Source_group,
    Group,
    Value
  ) %>%
  dplyr::arrange(
    Gene,
    Group,
    Sample
  )

cat(
  "\nKPNB1 and CCT5 missing-value audit:\n"
)

if (nrow(KPNB1_CCT5_missing) == 0L) {
  
  cat(
    "No missing KPNB1 or CCT5 values.\n"
  )
  
} else {
  
  print(
    KPNB1_CCT5_missing,
    n = Inf
  )
}

write.csv(
  KPNB1_CCT5_missing,
  file.path(
    output_dir,
    "Figure6I_KPNB1_CCT5_missing_value_audit.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 21. Retain measured values without imputation
# -----------------------------------------------------------------------------

KPNB1_CCT5_data <- KPNB1_CCT5_all_data %>%
  dplyr::filter(
    is.finite(Value)
  ) %>%
  dplyr::select(
    Sample,
    Source_group,
    Group,
    Gene,
    Measure,
    Value
  )

KPNB1_CCT5_group_counts <- KPNB1_CCT5_data %>%
  dplyr::count(
    Measure,
    Group,
    name = "n_samples"
  ) %>%
  dplyr::arrange(
    Measure,
    Group
  )

cat(
  "\nKPNB1 and CCT5 valid sample counts:\n"
)

print(
  KPNB1_CCT5_group_counts,
  n = Inf
)

write.csv(
  KPNB1_CCT5_data,
  file.path(
    output_dir,
    "SourceData_Figure6I_KPNB1_CCT5_Three_Group.csv"
  ),
  row.names = FALSE
)

write.csv(
  KPNB1_CCT5_group_counts,
  file.path(
    output_dir,
    "SourceData_Figure6I_KPNB1_CCT5_valid_sample_counts.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 22. Descriptive statistics
# -----------------------------------------------------------------------------

KPNB1_CCT5_descriptive <- KPNB1_CCT5_data %>%
  dplyr::group_by(
    Measure,
    Group
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median = stats::median(Value),
    mean = base::mean(Value),
    SD = stats::sd(Value),
    .groups = "drop"
  )

write.csv(
  KPNB1_CCT5_descriptive,
  file.path(
    output_dir,
    "SourceData_Figure6I_KPNB1_CCT5_descriptive_statistics.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 23. Paired Normal-Primary Wilcoxon function
#
# Only the seven confirmed Normal-Primary pairs are considered.
# If one member of a pair is missing for a protein, that pair is excluded
# from that protein's paired test.
# -----------------------------------------------------------------------------

run_KPNB1_CCT5_paired_test <- function(
    measure_name
) {
  
  current_data <- KPNB1_CCT5_data %>%
    dplyr::filter(
      Measure == measure_name
    ) %>%
    dplyr::mutate(
      Group = as.character(Group)
    )
  
  paired_wide <- pair_design %>%
    dplyr::left_join(
      current_data %>%
        dplyr::select(
          Sample,
          Group,
          Value
        ),
      by = c(
        "Sample",
        "Group"
      )
    ) %>%
    dplyr::select(
      Pair_ID,
      Group,
      Value
    ) %>%
    tidyr::pivot_wider(
      names_from = Group,
      values_from = Value
    ) %>%
    dplyr::filter(
      is.finite(Normal),
      is.finite(Primary)
    ) %>%
    dplyr::arrange(
      Pair_ID
    )
  
  stopifnot(
    nrow(paired_wide) > 0L,
    nrow(paired_wide) <= 7L,
    all(
      complete.cases(
        paired_wide
      )
    )
  )
  
  paired_difference <- (
    paired_wide$Primary -
      paired_wide$Normal
  )
  
  exact_allowed <- (
    !any(paired_difference == 0) &&
      !anyDuplicated(
        abs(paired_difference)
      )
  )
  
  test_result <- stats::wilcox.test(
    x = paired_wide$Primary,
    y = paired_wide$Normal,
    paired = TRUE,
    alternative = "two.sided",
    exact = exact_allowed
  )
  
  tibble::tibble(
    Measure = measure_name,
    
    Comparison =
      "Primary vs Normal",
    
    Design =
      "Two-sided paired Wilcoxon signed-rank test",
    
    n_group_1 =
      nrow(paired_wide),
    
    n_group_2 =
      nrow(paired_wide),
    
    Median_group_1 =
      stats::median(
        paired_wide$Primary
      ),
    
    Median_group_2 =
      stats::median(
        paired_wide$Normal
      ),
    
    Median_difference =
      stats::median(
        paired_difference
      ),
    
    Wilcoxon_W =
      unname(
        test_result$statistic
      ),
    
    Exact_P_used =
      exact_allowed,
    
    Raw_P =
      test_result$p.value
  )
}


# -----------------------------------------------------------------------------
# 24. Unpaired Wilcoxon function
# -----------------------------------------------------------------------------

run_KPNB1_CCT5_unpaired_test <- function(
    measure_name,
    group_1,
    group_2
) {
  
  current_data <- KPNB1_CCT5_data %>%
    dplyr::filter(
      Measure == measure_name
    )
  
  x <- current_data %>%
    dplyr::filter(
      Group == group_1
    ) %>%
    dplyr::pull(
      Value
    )
  
  y <- current_data %>%
    dplyr::filter(
      Group == group_2
    ) %>%
    dplyr::pull(
      Value
    )
  
  stopifnot(
    length(x) > 0L,
    length(y) > 0L,
    all(is.finite(x)),
    all(is.finite(y))
  )
  
  test_result <- stats::wilcox.test(
    x = x,
    y = y,
    paired = FALSE,
    alternative = "two.sided",
    exact = FALSE
  )
  
  tibble::tibble(
    Measure = measure_name,
    
    Comparison =
      paste(
        group_1,
        "vs",
        group_2
      ),
    
    Design =
      "Two-sided unpaired Wilcoxon rank-sum test",
    
    n_group_1 =
      length(x),
    
    n_group_2 =
      length(y),
    
    Median_group_1 =
      stats::median(x),
    
    Median_group_2 =
      stats::median(y),
    
    Median_difference =
      stats::median(x) -
      stats::median(y),
    
    Wilcoxon_W =
      unname(
        test_result$statistic
      ),
    
    Exact_P_used =
      FALSE,
    
    Raw_P =
      test_result$p.value
  )
}


# -----------------------------------------------------------------------------
# 25. Run all prespecified comparisons
# -----------------------------------------------------------------------------

KPNB1_CCT5_statistics <- dplyr::bind_rows(
  lapply(
    levels(
      KPNB1_CCT5_data$Measure
    ),
    function(measure_name) {
      
      dplyr::bind_rows(
        run_KPNB1_CCT5_paired_test(
          measure_name
        ),
        
        run_KPNB1_CCT5_unpaired_test(
          measure_name,
          "Primary",
          "Metastatic"
        ),
        
        run_KPNB1_CCT5_unpaired_test(
          measure_name,
          "Normal",
          "Metastatic"
        )
      )
    }
  )
) %>%
  dplyr::group_by(
    Measure
  ) %>%
  dplyr::mutate(
    BH_adjusted_P = stats::p.adjust(
      Raw_P,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

cat(
  "\nKPNB1 and CCT5 statistics:\n"
)

print(
  KPNB1_CCT5_statistics,
  n = Inf
)

write.csv(
  KPNB1_CCT5_statistics,
  file.path(
    output_dir,
    "SourceData_Figure6I_KPNB1_CCT5_statistics.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 26. Audit the confirmed pairs used for each protein
# -----------------------------------------------------------------------------

KPNB1_CCT5_pair_audit <- dplyr::bind_rows(
  lapply(
    levels(
      KPNB1_CCT5_data$Measure
    ),
    function(measure_name) {
      
      current_data <- KPNB1_CCT5_data %>%
        dplyr::filter(
          Measure == measure_name
        ) %>%
        dplyr::mutate(
          Group = as.character(Group)
        )
      
      pair_design %>%
        dplyr::left_join(
          current_data %>%
            dplyr::select(
              Sample,
              Group,
              Value
            ),
          by = c(
            "Sample",
            "Group"
          )
        ) %>%
        dplyr::select(
          Pair_ID,
          Group,
          Sample,
          Value
        ) %>%
        tidyr::pivot_wider(
          names_from = Group,
          values_from = c(
            Sample,
            Value
          )
        ) %>%
        dplyr::mutate(
          Measure = measure_name,
          
          Complete_pair = (
            is.finite(Value_Normal) &
              is.finite(Value_Primary)
          )
        ) %>%
        dplyr::select(
          Measure,
          Pair_ID,
          Sample_Normal,
          Sample_Primary,
          Value_Normal,
          Value_Primary,
          Complete_pair
        )
    }
  )
)

write.csv(
  KPNB1_CCT5_pair_audit,
  file.path(
    output_dir,
    "SourceData_Figure6I_KPNB1_CCT5_confirmed_pair_audit.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 27. Three-group plotting function
# -----------------------------------------------------------------------------

plot_KPNB1_CCT5_three_groups <- function(
    measure_name,
    y_axis_title
) {
  
  plot_data <- KPNB1_CCT5_data %>%
    dplyr::filter(
      Measure == measure_name
    )
  
  stopifnot(
    nrow(plot_data) > 0L,
    all(is.finite(plot_data$Value))
  )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Group,
      y = Value,
      fill = Group
    )
  ) +
    ggplot2::geom_violin(
      trim = FALSE,
      scale = "width",
      width = 0.80,
      color = NA,
      alpha = 0.65
    ) +
    ggplot2::geom_boxplot(
      width = 0.28,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.65
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(
        color = Group
      ),
      width = 0.13,
      height = 0,
      size = 1.5,
      alpha = 0.75,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = group_colors,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = c(
        "Normal",
        "Primary",
        "Metastatic"
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.05,
          0.10
        )
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_axis_title
    ) +
    publication_theme
}


# -----------------------------------------------------------------------------
# 28. Generate KPNB1 and CCT5 plots
# -----------------------------------------------------------------------------

three_group_KPNB1_plot <- plot_KPNB1_CCT5_three_groups(
  measure_name =
    "KPNB1 protein",
  
  y_axis_title =
    "KPNB1 protein abundance"
)

three_group_CCT5_plot <- plot_KPNB1_CCT5_three_groups(
  measure_name =
    "CCT5 protein",
  
  y_axis_title =
    "CCT5 protein abundance"
)

three_group_KPNB1_CCT5_combined_plot <- (
  three_group_KPNB1_plot +
    three_group_CCT5_plot +
    patchwork::plot_layout(
      ncol = 2
    )
)

print(three_group_KPNB1_plot)
print(three_group_CCT5_plot)
print(three_group_KPNB1_CCT5_combined_plot)


# -----------------------------------------------------------------------------
# 29. Save KPNB1 and CCT5 figures
# -----------------------------------------------------------------------------

save_pdf(
  three_group_KPNB1_CCT5_combined_plot,
  "Fig6I_KPNB1_CCT5_Three_Group_Protein_Combined.pdf",
  9.0,
  4.5
)


# -----------------------------------------------------------------------------
# 30. Completion message
# -----------------------------------------------------------------------------

cat(
  "\nFigure 6I KPNB1 and CCT5 analysis completed.\n",
  "Normal-Primary comparison used confirmed paired samples only.\n",
  "Primary-Metastatic and Normal-Metastatic comparisons used ",
  "two-sided unpaired Wilcoxon rank-sum tests.\n",
  "Missing KPNB1 or CCT5 values were excluded without imputation.\n",
  "P values were saved to CSV files and were not displayed in figures.\n",
  "Output directory: ",
  output_dir,
  "\n",
  sep = ""
)

}))


# =============================================================================
# SUPPLEMENTARY FIGURES S1–S8 (PPT pages 7–14)
# Continues in the same analysis run immediately after Fig6I.
# FigS4H and experimental-only panels are not regenerated.
# FigS7D uses all QC-eligible RPL/RPS genes, as in Fig2F.
# FigS2A and FigS6B plotting are updated as documented in README.md.
# =============================================================================

# 1) Extend the Main figure index with Supplementary panels. This must happen
#    before saving any FigS PDF, so they go to FIGURES_NEW_ORDER, not EXTRA_ANALYSES.
if (!exists(".nmd3_panels", inherits = TRUE) ||
    !exists(".nmd3_state", inherits = TRUE) ||
    !"Fig6I" %in% names(.nmd3_state$modules)) {
  stop("Supplementary modules require completed Fig6I in the same analysis run.")
}
if (any(startsWith(.nmd3_panels$Panel, "FigS"))) {
  stop("Supplementary index is already installed; refusing duplicate panels.")
}
.nmd3_supp_panels <- data.frame(do.call(rbind, list(
  c("FigS1A","","Experimental design schematic","non-R schematic"),
  c("FigS1B","FigS1B_NMD3_Promoted_GSEA.pdf","siNMD3 cytoplasmic protein GSEA","computational"),
  c("FigS1C","","siNMD3 viability and knockdown WB","experimental"),
  c("FigS1D","","shNMD3 colony assay","experimental"),
  c("FigS1E","","Portimine A viability","experimental"),
  c("FigS1F","","BPH1 puromycin pulse","experimental"),
  c("FigS1G","FigS1G_NMD3_Cell_Lines_GSE4016_GSE21032.pdf","NMD3 expression across cell lines (computational upper panels; WB excluded)","partial"),
  c("FigS2A","FigS2A_NMD3_High_vs_Low_Three_GSEA.pdf","TCGA NMD3-high KEGG/Reactome/export GSEA, 3 panels","computational"),
  c("FigS2B","FigS2B_XPO1_NMD3_Spearman_TCGA_MSKCC.pdf","TCGA/MSKCC NMD3-XPO1 correlation","computational"),
  c("FigS2C","FigS2C_NMD3_Protein_Unadjusted.pdf","Unadjusted NMD3 abundance, 28 primary / 20 metastatic","computational"),
  c("FigS2D","FigS2D_RB_Adjusted_NMD3_Protein_BlackPoints.pdf","Public proteome RB-adjusted NMD3 protein","computational"),
  c("FigS2E","FigS2E_RPL_Adjusted_NMD3_Protein_BlackPoints.pdf","Public proteome RPL-adjusted NMD3 protein","computational"),
  c("FigS2F","FigS2F_XPO1_High_Low_NMD3_expression.pdf","TCGA XPO1-low/high NMD3 expression","computational"),
  c("FigS2G","FigS2G_SU2C_NMD3_Survival.pdf","SU2C NMD3 high/low survival (duplicate use of one plot)","computational"),
  c("FigS3A","","XPO1/NMD3 nuclear/cytoplasmic fractionation WB","experimental"),
  c("FigS3B","FigS3B_NMD3_Nuclear_Cytoplasmic_Protein_Ratio.pdf","NMD3 nucleus/cytoplasm protein ratio after KPT-330","computational"),
  c("FigS3C","","KPT-330 cell viability plus WB","experimental"),
  c("FigS3D","","siXPO1 cell viability plus WB","experimental"),
  c("FigS3E","","siXPO1 colony formation","experimental"),
  c("FigS4A","FigS4A_transport_correlation.pdf","TCGA nuclear transport Spearman correlations","computational"),
  c("FigS4B","FigS4B_TRANSLATION_GSEA.pdf","TCGA REACTOME_TRANSLATION tumor-normal GSEA","computational"),
  c("FigS4C","FigS4C_NUCLEAR_TRANSPORT_GSEA.pdf","TCGA GO NUCLEAR_TRANSPORT tumor-normal GSEA","computational"),
  c("FigS4D","FigS4D_XPO1_Nuclear_Transport_GSEA.pdf","TCGA XPO1 high/low nuclear-transport enrichment","computational"),
  c("FigS4E","","RNA/protein compartment masses from fractionation (raw experimental data missing)","experimental"),
  c("FigS4F","FigS4F_KPT330_RNAseq_PCA.pdf","KPT-330 nucleus/cytoplasm RNA PCA","computational"),
  c("FigS4G","FigS4G_KPT330_Proteomics_PCA.pdf","KPT-330 nucleus/cytoplasm DIA PCA","computational"),
  c("FigS4H","","Manually prepared retained-protein family summary (do not regenerate)","user-provided figure"),
  c("FigS4I","FigS4I_KPT330_vs_DMSO_Two_GO_GSEA.pdf","GO cytoplasmic translation initiation and ribosome assembly GSEA","computational"),
  c("FigS5A","FigS5A_KPT330_Nucleus_Cytoplasm_Volcano.pdf","KPT-330 nuclear/cytoplasmic RNA volcanoes","computational"),
  c("FigS5B","FigS5B_KPT330_Suppressed_GOBP_KEGG_Pathways.pdf","KPT-330 suppressed GO/KEGG RNA pathways","computational"),
  c("FigS5C","FigS5C_KPT330_Suppressed_Hallmark_Pathways.pdf","KPT-330 suppressed Hallmark RNA pathways","computational"),
  c("FigS5D","","Nuclear/cytoplasmic KPT-330 and siXPO1 WB","experimental"),
  c("FigS6A","FigS6A_RNA_Proteomics_Eight_Group_PCA.pdf","siNMD3 RNA and proteomics PCA","computational"),
  c("FigS6B","FigS6B_siNMD3_RNA_Volcano_Nucleus_Cytoplasm.pdf","siNMD3 nuclear/cytoplasmic RNA volcano; FDR<0.05 and absolute log2FC >=1; no PEG10 label","computational"),
  c("FigS6C","FigS6C_siNMD3_Proteomics_Volcano_Nucleus_Cytoplasm.pdf","siNMD3 nuclear/cytoplasmic protein volcano","computational"),
  c("FigS6D","FigS6D_final.pdf","Venn overlaps, four sets of nuclear/cytoplasm RNA and protein perturbations","computational"),
  c("FigS6E-top","FigS6E_top_NoRB_Spearman.pdf","Spearman correlations excluding GO:0042254","computational"),
  c("FigS6E-bottom","FigS6E_bottom_NoRPL_Spearman.pdf","Spearman correlations excluding RPL, including statistics","computational"),
  c("FigS7A","FigS7A_RNA_Protein.pdf","KPT-330 cytoplasmic RNA/protein RRHO","computational"),
  c("FigS7B","FigS7B_H1_H2_GO_BP_Top10.pdf","KPT-330 RRHO H1/H2 GO enrichment","computational"),
  c("FigS7C","FigS7C_siNMD3_RRHO_RPL_RPS.pdf","siNMD3 cytoplasmic RNA/protein RRHO","computational"),
  c("FigS7D","FigS7D_siNMD3_RPL_RPS_Heatmaps.pdf","ALL eligible RPL/RPS, siNMD3 nuclear/cytoplasmic protein heatmaps","computational"),
  c("FigS7E","FigS7E_MYC_Residual_NoRB_NoRPLRPS.pdf","MYC residual excluding RB and RPL/RPS","computational"),
  c("FigS7F","FigS7F_MYC_HalfLife_NoRPLRPS.pdf","MYC V1 B-cell reference half-life without RPL/RPS","computational"),
  c("FigS8A","FigS8A_siNMD3_Cytoplasmic_Downregulated_Protein_ORA.pdf","siNMD3 cytoplasmic downregulated protein ORA","computational"),
  c("FigS8B","","NMD3/XPO1 clinical IHC H-scores","experimental"),
  c("FigS8C","","NMD3/XPO1 IHC expression correlation","experimental"),
  c("FigS8D","FigS8D_signature_GSEA.pdf","NEPC signature GSEA after siNMD3","computational"),
  c("FigS8E","FigS8E_Adeno_signature_GSEA.pdf","CRPC Adeno signature GSEA after siNMD3","computational"),
  c("FigS8F","FigS8F_SU2C_NMD3_PEG10_ENO2_Signature_Survival.pdf","SU2C NMD3+PEG10+ENO2 signature survival","computational")
)), stringsAsFactors = FALSE)
names(.nmd3_supp_panels) <- c("Panel", "PDF", "Description", "Coverage")
.nmd3_supp_panels$Status <- ifelse(nzchar(.nmd3_supp_panels$PDF), "PENDING", "NOT_R_INPUT")
.nmd3_supp_panels$Path <- ""
.nmd3_panels <- rbind(.nmd3_panels, .nmd3_supp_panels)
nmd3_write_status()
nmd3_log("[SUPPLEMENTARY] Index added: ", nrow(.nmd3_supp_panels),
         " panels; saved PDFs will follow Main in the SAME folder.")

# 2) Locate ONLY the five additional public input files needed for S1–S8.
#    Main's preflight intentionally checked only Main-figure inputs.
#    We use the same input paths expected by the original full analysis;
#    never substitute a different cohort or change any data contents.
.nmd3_supp_inputs <- c(
  "data/prad_mskcc/data_mrna_agilent_microarray.txt",
  "data/prad_mskcc/data_clinical_sample.txt",
  "data/GSE4016_series_matrix.txt",
  "data/prad_su2c_2019/expression_with_clinical.xlsx",
  "data/prad_su2c_2019/data_clinical_sample.txt"
)
# Supplementary files are also expected under the project data root.
# This keeps the same original input datasets and prevents cohort-name collisions.
.nmd3_supp_missing <- character(0)
for (.nmd3_supp_path in .nmd3_supp_inputs) {
  .nmd3_supp_real <- file.path(.nmd3_data_root,
    sub("^data/", "", .nmd3_supp_path))
  if (!.nmd3_is_present(.nmd3_supp_real)) {
    .nmd3_supp_missing <- c(.nmd3_supp_missing, .nmd3_supp_path)
  } else {
    .nmd3_link(.nmd3_supp_path, .nmd3_supp_real)
  }
}
if (length(.nmd3_supp_missing)) {
  stop("Required supplementary input files missing: ",
       paste(.nmd3_supp_missing, collapse=", "),
       ". See README.md and place them below NMD3_DATA_ROOT.")
}

# 3) Keep completion checks and describe the combined run.
.nmd3_finish_core <- nmd3_finish
nmd3_finish <- function() {
  .nmd3_finish_core()
  writeLines(c(
    "MAIN + SUPPLEMENTARY: Figures 1–6 followed by S1–S8, in PPT sequence.",
    "All statistical procedures and thresholds come from original analyses.",
    "Fig2F and FigS7D include all data-eligible RPL/RPS, applying original QC.",
    "FigS4H and experimental/schematic-only panels were not regenerated.",
    "Each PDF was saved and closed immediately after its panel analysis."
  ), file.path(.nmd3_run_root, "SOURCE_NOTES.txt"))
  nmd3_log("[COMPLETE] MAIN + SUPPLEMENTARY generated in ", .nmd3_order_folder)
  invisible(NULL)
}

# 4) Original supplementary modules in exact page/panel order follow.

# =============================================================================
# FigS1B
# =============================================================================
nmd3_step("FigS1B", character(), quote({
  # FigS1B — analysis below
  # Cytoplasmic protein GSEA: siNMD3 vs siCTRL
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(msigdbr)
    library(clusterProfiler)
    library(ggplot2)
  })
  
  # prot_raw was loaded in the previous figures
  
  # ------------------------------------------------------------
  # 1. Extract siCTRL / siNMD3 cytoplasmic protein results
  # ------------------------------------------------------------
  
  ctrl_cols <- grep(
    "^Norm_siRNAxCtrl_Cyto_",
    colnames(prot_raw),
    value = TRUE
  )
  
  nmd3_cols <- grep(
    "^Norm_siRNAxNMD3xKD_Cyto_",
    colnames(prot_raw),
    value = TRUE
  )
  
  stopifnot(
    length(ctrl_cols) == 3,
    length(nmd3_cols) == 3
  )
  
  expr_cols <- c(ctrl_cols, nmd3_cols)
  
  fc_col <- grep(
    "Log2FC.*siRNAxNMD3xKD_Cyto.*siRNAxCtrl_Cyto",
    colnames(prot_raw),
    value = TRUE
  )
  
  flip_sign <- FALSE
  
  if (length(fc_col) != 1) {
    
    fc_col <- grep(
      "Log2FC.*siRNAxCtrl_Cyto.*siRNAxNMD3xKD_Cyto",
      colnames(prot_raw),
      value = TRUE
    )
    
    flip_sign <- TRUE
  }
  
  stopifnot(length(fc_col) == 1)
  
  prot_5A <- prot_raw %>%
    dplyr::filter(!is.na(Gene)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(expr_cols),
        ~ suppressWarnings(as.numeric(.x))
      )
    ) %>%
    dplyr::filter(
      rowSums(is.na(dplyr::across(dplyr::all_of(expr_cols)))) == 0
    ) %>%
    dplyr::transmute(
      gene = toupper(as.character(Gene)),
      log2FC = if (flip_sign) {
        -as.numeric(.data[[fc_col]])
      } else {
        as.numeric(.data[[fc_col]])
      }
    ) %>%
    dplyr::filter(is.finite(log2FC)) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      log2FC = median(log2FC),
      .groups = "drop"
    )
  
  
  # ------------------------------------------------------------
  # 2. Protein ranking for GSEA
  # ------------------------------------------------------------
  
  gene_ranking <- prot_5A$log2FC
  names(gene_ranking) <- prot_5A$gene
  
  gene_ranking <- sort(
    gene_ranking,
    decreasing = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 3. GSEA function
  # ------------------------------------------------------------
  
  run_gsea <- function(gene_ranking, term2gene,
                       min_size = 10,
                       max_size = 3000) {
    
    term2gene <- term2gene %>%
      dplyr::transmute(
        term = as.character(term),
        gene = toupper(as.character(gene))
      ) %>%
      dplyr::distinct()
    
    pathways <- split(
      term2gene$gene,
      term2gene$term
    )
    
    pathways <- lapply(
      pathways,
      function(x) intersect(x, names(gene_ranking))
    )
    
    sizes <- lengths(pathways)
    
    pathways <- pathways[
      sizes >= min_size &
        sizes <= max_size
    ]
    
    # Preserve the small tie-breaking jitter used in the original analysis
    set.seed(12345)
    
    ranking2 <- gene_ranking +
      rnorm(length(gene_ranking), sd = 1e-6)
    
    ranking2 <- sort(ranking2, decreasing = TRUE)
    
    t2g <- stack(pathways)
    colnames(t2g) <- c("gene", "term")
    
    clusterProfiler::GSEA(
      geneList = ranking2,
      TERM2GENE = t2g[, c("term", "gene")],
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = 1,
      seed = TRUE,
      verbose = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # 4. GO gene sets
  # ------------------------------------------------------------
  
  go_bp <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  ) %>%
    dplyr::transmute(
      term = gs_name,
      gene = gene_symbol
    ) %>%
    dplyr::distinct()
  
  go_cc <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:CC"
  ) %>%
    dplyr::transmute(
      term = gs_name,
      gene = gene_symbol
    ) %>%
    dplyr::distinct()
  
  go_mf <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:MF"
  ) %>%
    dplyr::transmute(
      term = gs_name,
      gene = gene_symbol
    ) %>%
    dplyr::distinct()
  
  
  # ------------------------------------------------------------
  # 5. Run GSEA
  # ------------------------------------------------------------
  
  gsea_bp <- run_gsea(gene_ranking, go_bp)
  gsea_cc <- run_gsea(gene_ranking, go_cc)
  gsea_mf <- run_gsea(gene_ranking, go_mf)
  
  
  # ------------------------------------------------------------
  # 6. Pathways shown in Figure 5A
  # ------------------------------------------------------------
  
  wanted_bp <- c(
    "GOBP_CYTOPLASMIC_TRANSLATION",
    "GOBP_TRANSLATIONAL_INITIATION",
    "GOBP_REGULATION_OF_TRANSLATIONAL_INITIATION",
    "GOBP_REGULATION_OF_TRANSLATION",
    "GOBP_TRANSLATION"
  )
  
  wanted_cc <- c(
    "GOCC_CYTOSOLIC_RIBOSOME",
    "GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT",
    "GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT",
    "GOCC_PRERIBOSOME"
  )
  
  wanted_mf <- c(
    "GOMF_TRANSLATION_INITIATION_FACTOR_ACTIVITY",
    "GOMF_TRANSLATION_FACTOR_ACTIVITY",
    "GOMF_TRANSLATION_REGULATOR_ACTIVITY",
    "GOMF_NUCLEAR_EXPORT_SIGNAL_RECEPTOR_ACTIVITY",
    "GOMF_RIBOSOME_BINDING"
  )
  
  
  # ------------------------------------------------------------
  # 7. Extract selected pathways
  # ------------------------------------------------------------
  
  get_selected <- function(gsea_object, pathways, category) {
    
    as.data.frame(gsea_object@result) %>%
      tibble::as_tibble() %>%
      dplyr::filter(ID %in% pathways) %>%
      dplyr::transmute(
        ID = ID,
        NES = -NES,
        FDR = p.adjust,
        Category = category
      )
  }
  
  bp_plot <- get_selected(
    gsea_bp,
    wanted_bp,
    "GO:BP"
  )
  
  cc_plot <- get_selected(
    gsea_cc,
    wanted_cc,
    "GO:CC"
  )
  
  mf_plot <- get_selected(
    gsea_mf,
    wanted_mf,
    "GO:MF"
  )
  
  
  # ------------------------------------------------------------
  # 8. Labels
  # ------------------------------------------------------------
  
  label_map <- c(
    "GOBP_CYTOPLASMIC_TRANSLATION" =
      "Cytoplasmic translation",
    
    "GOBP_TRANSLATIONAL_INITIATION" =
      "Translation initiation",
    
    "GOBP_REGULATION_OF_TRANSLATIONAL_INITIATION" =
      "Regulation of translation initiation",
    
    "GOBP_REGULATION_OF_TRANSLATION" =
      "Regulation of translation",
    
    "GOBP_TRANSLATION" =
      "Translation",
    
    "GOCC_CYTOSOLIC_RIBOSOME" =
      "Cytoplasmic ribosome",
    
    "GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT" =
      "Cytosolic large ribosomal subunit",
    
    "GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT" =
      "Cytosolic small ribosomal subunit",
    
    "GOCC_PRERIBOSOME" =
      "Preribosome",
    
    "GOMF_TRANSLATION_INITIATION_FACTOR_ACTIVITY" =
      "Translation initiation factor activity",
    
    "GOMF_TRANSLATION_FACTOR_ACTIVITY" =
      "Translation factor activity",
    
    "GOMF_TRANSLATION_REGULATOR_ACTIVITY" =
      "Translation regulator activity",
    
    "GOMF_NUCLEAR_EXPORT_SIGNAL_RECEPTOR_ACTIVITY" =
      "Nuclear export signal receptor activity",
    
    "GOMF_RIBOSOME_BINDING" =
      "Ribosome binding"
  )
  
  
  plot_df <- dplyr::bind_rows(
    bp_plot,
    cc_plot,
    mf_plot
  ) %>%
    dplyr::mutate(
      Pathway = unname(label_map[ID]),
      neglog10FDR = -log10(FDR),
      Category = factor(
        Category,
        levels = c("GO:BP", "GO:CC", "GO:MF")
      )
    ) %>%
    dplyr::group_by(Category) %>%
    dplyr::arrange(Category, dplyr::desc(NES)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Pathway = factor(
        Pathway,
        levels = rev(unique(Pathway))
      )
    )
  
  
  # ------------------------------------------------------------
  # 9. Figure 5A
  # ------------------------------------------------------------
  
  p_5A <- ggplot(
    plot_df,
    aes(
      x = NES,
      y = Pathway,
      fill = neglog10FDR
    )
  ) +
    geom_col(width = 0.75) +
    
    scale_fill_gradientn(
      colours = c(
        "#4575B4",
        "#74ADD1",
        "#ABD9E9",
        "#FEE090",
        "#FDAE61",
        "#F46D43",
        "#D73027",
        "#A50026"
      ),
      name = expression(-log[10]~FDR)
    ) +
    
    facet_grid(
      Category ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    
    labs(
      title = "NMD3 Promoted",
      x = "NES",
      y = NULL
    ) +
    
    theme_classic(base_size = 13) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      
      axis.text.y = element_text(
        size = 11,
        color = "black",
        face = "bold"
      ),
      
      axis.text.x = element_text(
        size = 11,
        color = "black"
      ),
      
      axis.title.x = element_text(
        size = 12,
        face = "bold"
      ),
      
      strip.text = element_text(
        size = 10,
        face = "bold"
      ),
      
      strip.background = element_rect(
        fill = "grey90",
        color = "black"
      ),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      
      panel.grid.major.x = element_line(
        color = "grey90",
        linewidth = 0.3
      )
    )
  
  
  # ------------------------------------------------------------
  # 10. Save
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    "FigS1B_NMD3_Promoted_GSEA.pdf",
    p_5A,
    width = 8,
    height = 10,
    device = cairo_pdf
  )
  
  
}))

# =============================================================================
# FigS1G
# =============================================================================
nmd3_step("FigS1G", character(), quote({
  msk_raw_S4D <- data.table::fread("data/prad_mskcc/data_mrna_agilent_microarray.txt", sep="\t",data.table=FALSE,check.names=FALSE)
  # FigS1G — analysis below
  # NMD3 expression across prostate cancer cell lines
  #
  # GSE4016:
  #   NMD3 probes identified from GPL3341 annotation
  #
  # GSE21032/MSKCC:
  #   NMD3 identified using Entrez Gene ID 51068
  # ============================================================
  
  suppressPackageStartupMessages({
    library(GEOquery)
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
  })
  
  # ============================================================
  # 1. Read GSE4016 series matrix
  # ============================================================
  
  gse4016_file_S4H <- paste0(
    "data/",
    "GSE4016_series_matrix.txt"
  )
  
  raw_lines_S4H <- readLines(
    gse4016_file_S4H,
    warn = FALSE
  )
  
  table_start_S4H <- grep(
    "^!series_matrix_table_begin",
    raw_lines_S4H
  )
  
  table_end_S4H <- grep(
    "^!series_matrix_table_end",
    raw_lines_S4H
  )
  
  stopifnot(
    length(table_start_S4H) == 1,
    length(table_end_S4H) == 1,
    table_end_S4H >
      table_start_S4H + 1
  )
  
  gse4016_table_lines_S4H <- raw_lines_S4H[
    (table_start_S4H + 1):
      (table_end_S4H - 1)
  ]
  
  gse4016_data_S4H <- read.table(
    text = gse4016_table_lines_S4H,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  
  cat(
    "GSE4016 expression matrix:",
    nrow(gse4016_data_S4H),
    "probes x",
    ncol(gse4016_data_S4H),
    "samples\n"
  )
  
  # ============================================================
  # 2. Extract GEO platform ID
  # ============================================================
  
  platform_lines_S4H <- grep(
    "^!Sample_platform_id",
    raw_lines_S4H,
    value = TRUE
  )
  
  stopifnot(
    length(platform_lines_S4H) >= 1
  )
  
  platform_matches_S4H <- regmatches(
    platform_lines_S4H,
    gregexpr(
      "GPL[0-9]+",
      platform_lines_S4H
    )
  )
  
  platform_ids_S4H <- unique(
    unlist(
      platform_matches_S4H
    )
  )
  
  stopifnot(
    length(platform_ids_S4H) == 1
  )
  
  gpl_id_S4H <- platform_ids_S4H[1]
  
  cat(
    "GSE4016 platform:",
    gpl_id_S4H,
    "\n"
  )
  
  # ============================================================
  # 3. Retrieve platform annotation
  # ============================================================
  
  gpl_object_S4H <- GEOquery::getGEO(
    GEO = gpl_id_S4H,
    AnnotGPL = TRUE
  )
  
  gpl_table_S4H <- GEOquery::Table(
    gpl_object_S4H
  )
  
  stopifnot(
    "ID" %in%
      colnames(gpl_table_S4H)
  )
  
  # ------------------------------------------------------------
  # Identify gene-symbol column
  # ------------------------------------------------------------
  
  symbol_column_priority_S4H <- c(
    "Gene Symbol",
    "GENE_SYMBOL",
    "Gene symbol",
    "Gene_Symbol",
    "SYMBOL",
    "Symbol"
  )
  
  symbol_columns_S4H <- intersect(
    symbol_column_priority_S4H,
    colnames(gpl_table_S4H)
  )
  
  if (length(symbol_columns_S4H) == 0) {
    
    symbol_columns_S4H <- grep(
      "gene.?symbol|symbol",
      colnames(gpl_table_S4H),
      value = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(symbol_columns_S4H) == 0) {
    stop(
      paste0(
        "No gene-symbol column was found in ",
        gpl_id_S4H,
        "."
      )
    )
  }
  
  symbol_column_S4H <-
    symbol_columns_S4H[1]
  
  cat(
    "Gene-symbol column:",
    symbol_column_S4H,
    "\n"
  )
  
  # ============================================================
  # 4. Identify probes annotated to NMD3
  # ============================================================
  
  gpl_symbols_S4H <- as.character(
    gpl_table_S4H[[symbol_column_S4H]]
  )
  
  contains_NMD3_symbol_S4H <- function(x) {
    
    tokens <- strsplit(
      toupper(
        as.character(x)
      ),
      "[^A-Z0-9-]+"
    )
    
    vapply(
      tokens,
      function(z) {
        "NMD3" %in% z
      },
      logical(1)
    )
  }
  
  is_NMD3_probe_S4H <-
    contains_NMD3_symbol_S4H(
      gpl_symbols_S4H
    )
  
  nmd3_probe_ids_S4H <- unique(
    as.character(
      gpl_table_S4H$ID[
        is_NMD3_probe_S4H
      ]
    )
  )
  
  # Retain probes that are actually present in the series matrix
  
  nmd3_probe_ids_S4H <- intersect(
    nmd3_probe_ids_S4H,
    rownames(gse4016_data_S4H)
  )
  
  if (length(nmd3_probe_ids_S4H) == 0) {
    stop(
      paste0(
        "No GPL-annotated NMD3 probe was found ",
        "in the GSE4016 expression matrix."
      )
    )
  }
  
  cat(
    "GSE4016 NMD3 probes:",
    paste(
      nmd3_probe_ids_S4H,
      collapse = ", "
    ),
    "\n"
  )
  
  # ============================================================
  # 5. Extract GSE4016 cell-line expression
  # ============================================================
  
  gse4016_sample_map_S4H <- c(
    "GSM91921" = "LNCaP",
    "GSM91918" = "LAPC4",
    "GSM91922" = "22Rv1",
    "GSM91920" = "DU145",
    "GSM91917" = "PC3"
  )
  
  gse4016_cell_order_S4H <- c(
    "LNCaP",
    "LAPC4",
    "22Rv1",
    "DU145",
    "PC3"
  )
  
  gse4016_samples_S4H <- names(
    gse4016_sample_map_S4H
  )
  
  stopifnot(
    all(
      gse4016_samples_S4H %in%
        colnames(gse4016_data_S4H)
    )
  )
  
  nmd3_probe_df_GSE4016_S4H <-
    gse4016_data_S4H[
      nmd3_probe_ids_S4H,
      gse4016_samples_S4H,
      drop = FALSE
    ]
  
  nmd3_probe_df_GSE4016_S4H[] <- lapply(
    nmd3_probe_df_GSE4016_S4H,
    function(x) {
      suppressWarnings(
        as.numeric(
          as.character(x)
        )
      )
    }
  )
  
  nmd3_probe_matrix_GSE4016_S4H <- as.matrix(
    nmd3_probe_df_GSE4016_S4H
  )
  
  stopifnot(
    nrow(
      nmd3_probe_matrix_GSE4016_S4H
    ) >= 1,
    all(
      is.finite(
        nmd3_probe_matrix_GSE4016_S4H
      )
    )
  )
  
  # Average multiple probes if more than one probe maps to NMD3
  
  nmd3_raw_GSE4016_S4H <- colMeans(
    nmd3_probe_matrix_GSE4016_S4H,
    na.rm = TRUE
  )
  
  # Preserve the original relative-expression transformation
  
  nmd3_shifted_GSE4016_S4H <-
    nmd3_raw_GSE4016_S4H -
    min(
      nmd3_raw_GSE4016_S4H,
      na.rm = TRUE
    ) +
    0.01
  
  df_GSE4016_S4H <- data.frame(
    Dataset = "GSE4016",
    Sample = gse4016_samples_S4H,
    CellLine = unname(
      gse4016_sample_map_S4H
    ),
    RawExpression = as.numeric(
      nmd3_raw_GSE4016_S4H[
        gse4016_samples_S4H
      ]
    ),
    Expression = as.numeric(
      nmd3_shifted_GSE4016_S4H[
        gse4016_samples_S4H
      ]
    ),
    NMD3_Probe_Count =
      length(nmd3_probe_ids_S4H),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      CellLine = factor(
        CellLine,
        levels =
          gse4016_cell_order_S4H
      )
    ) %>%
    dplyr::arrange(
      CellLine
    )
  
  stopifnot(
    nrow(df_GSE4016_S4H) == 5,
    all(
      is.finite(
        df_GSE4016_S4H$Expression
      )
    )
  )
  
  cat("\nGSE4016 NMD3 expression:\n")
  
  print(df_GSE4016_S4H)
  
  # ============================================================
  # 6. Reuse the MSKCC expression matrix loaded in Figure S4D
  # ============================================================
  
  stopifnot(
    exists("msk_raw_S4D", inherits = TRUE),
    "Entrez_Gene_Id" %in% colnames(msk_raw_S4D)
  )
  
  msk_expression_S4H <- msk_raw_S4D
  
  stopifnot(
    "Entrez_Gene_Id" %in%
      colnames(msk_expression_S4H)
  )
  
  # Human NMD3 Entrez Gene ID
  
  nmd3_entrez_S4H <- 51068
  
  nmd3_rows_MSK_S4H <- which(
    as.character(
      msk_expression_S4H$Entrez_Gene_Id
    ) ==
      as.character(
        nmd3_entrez_S4H
      )
  )
  
  stopifnot(
    length(nmd3_rows_MSK_S4H) > 0
  )
  
  # ============================================================
  # 7. Extract MSKCC cell-line expression
  # ============================================================
  
  msk_cell_order_S4H <- c(
    "LNCaP104S",
    "LNCaP104R",
    "LNCaP",
    "DU145",
    "PC3"
  )
  
  stopifnot(
    all(
      msk_cell_order_S4H %in%
        colnames(msk_expression_S4H)
    )
  )
  
  nmd3_matrix_MSK_S4H <-
    msk_expression_S4H[
      nmd3_rows_MSK_S4H,
      msk_cell_order_S4H,
      drop = FALSE
    ]
  
  nmd3_matrix_MSK_S4H[] <- lapply(
    nmd3_matrix_MSK_S4H,
    function(x) {
      suppressWarnings(
        as.numeric(
          as.character(x)
        )
      )
    }
  )
  
  nmd3_matrix_MSK_S4H <- as.matrix(
    nmd3_matrix_MSK_S4H
  )
  
  if (nrow(nmd3_matrix_MSK_S4H) == 1) {
    
    nmd3_expression_MSK_S4H <- as.numeric(
      nmd3_matrix_MSK_S4H[1, ]
    )
    
  } else {
    
    nmd3_expression_MSK_S4H <- colMeans(
      nmd3_matrix_MSK_S4H,
      na.rm = TRUE
    )
  }
  
  names(nmd3_expression_MSK_S4H) <-
    msk_cell_order_S4H
  
  stopifnot(
    all(
      is.finite(
        nmd3_expression_MSK_S4H
      )
    )
  )
  
  df_MSK_S4H <- data.frame(
    Dataset = "GSE21032",
    Sample = msk_cell_order_S4H,
    CellLine = factor(
      msk_cell_order_S4H,
      levels = msk_cell_order_S4H
    ),
    RawExpression = as.numeric(
      nmd3_expression_MSK_S4H[
        msk_cell_order_S4H
      ]
    ),
    Expression = as.numeric(
      nmd3_expression_MSK_S4H[
        msk_cell_order_S4H
      ]
    ),
    NMD3_Probe_Count =
      length(nmd3_rows_MSK_S4H),
    stringsAsFactors = FALSE
  )
  
  stopifnot(
    nrow(df_MSK_S4H) == 5,
    all(
      is.finite(
        df_MSK_S4H$Expression
      )
    )
  )
  
  cat("\nGSE21032 NMD3 expression:\n")
  
  print(df_MSK_S4H)
  
  # ============================================================
  # 8. Colors
  # ============================================================
  
  colors_GSE4016_S4H <- c(
    "LNCaP" = "#F8766D",
    "LAPC4" = "#A3A500",
    "22Rv1" = "#00BF7D",
    "DU145" = "#00B0F6",
    "PC3" = "#E76BF3"
  )
  
  colors_MSK_S4H <- c(
    "LNCaP104S" = "#F8766D",
    "LNCaP104R" = "#A3A500",
    "LNCaP" = "#00BF7D",
    "DU145" = "#00B0F6",
    "PC3" = "#E76BF3"
  )
  
  # ============================================================
  # 9. GSE4016 plot
  # ============================================================
  
  p_GSE4016_S4H <- ggplot2::ggplot(
    df_GSE4016_S4H,
    ggplot2::aes(
      x = CellLine,
      y = Expression,
      fill = CellLine
    )
  ) +
    ggplot2::geom_col(
      width = 0.58,
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::scale_fill_manual(
      values = colors_GSE4016_S4H
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(0, 0.08)
      )
    ) +
    ggplot2::labs(
      title = "GSE4016",
      x = NULL,
      y = expression(
        italic(NMD3)~Expression
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title =
        ggplot2::element_text(
          face = "bold",
          hjust = 0.5,
          size = 15
        ),
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          color = "black",
          size = 11,
          face = "bold"
        ),
      axis.text.y =
        ggplot2::element_text(
          color = "black",
          size = 11
        ),
      axis.title.y =
        ggplot2::element_text(
          color = "black",
          size = 13,
          face = "bold"
        ),
      axis.line =
        ggplot2::element_line(
          color = "black",
          linewidth = 0.7
        ),
      axis.ticks =
        ggplot2::element_line(
          color = "black",
          linewidth = 0.6
        )
    )
  
  # ============================================================
  # 10. GSE21032 plot
  # ============================================================
  
  p_MSK_S4H <- ggplot2::ggplot(
    df_MSK_S4H,
    ggplot2::aes(
      x = CellLine,
      y = Expression,
      fill = CellLine
    )
  ) +
    ggplot2::geom_col(
      width = 0.58,
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::scale_fill_manual(
      values = colors_MSK_S4H
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(0, 0.08)
      )
    ) +
    ggplot2::labs(
      title = "GSE21032",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title =
        ggplot2::element_text(
          face = "bold",
          hjust = 0.5,
          size = 15
        ),
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          color = "black",
          size = 11,
          face = "bold"
        ),
      axis.text.y =
        ggplot2::element_text(
          color = "black",
          size = 11
        ),
      axis.line =
        ggplot2::element_line(
          color = "black",
          linewidth = 0.7
        ),
      axis.ticks =
        ggplot2::element_line(
          color = "black",
          linewidth = 0.6
        )
    )
  
  # ============================================================
  # 11. Combine Figure S4H
  # ============================================================
  
  p_S4H <- p_GSE4016_S4H +
    p_MSK_S4H +
    patchwork::plot_layout(
      nrow = 1,
      widths = c(1, 1)
    )
  
  
  
  
  # ============================================================
  # 12. Save outputs
  # ============================================================
  
  nmd3_save_pdf(
    filename =
      "FigS1G_NMD3_Cell_Lines_GSE4016_GSE21032.pdf",
    plot = p_S4H,
    width = 8,
    height = 4.8,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  data_S4H <- dplyr::bind_rows(
    df_GSE4016_S4H,
    df_MSK_S4H
  )
  
  utils::write.csv(
    data_S4H,
    file =
      "FigureS4H_NMD3_Cell_Line_expression_data.csv",
    row.names = FALSE
  )
  
  probe_information_S4H <- data.frame(
    Dataset = c(
      rep(
        "GSE4016",
        length(nmd3_probe_ids_S4H)
      ),
      rep(
        "GSE21032",
        length(nmd3_rows_MSK_S4H)
      )
    ),
    Probe_or_ID = c(
      nmd3_probe_ids_S4H,
      as.character(
        msk_expression_S4H$Entrez_Gene_Id[
          nmd3_rows_MSK_S4H
        ]
      )
    ),
    Mapping_method = c(
      rep(
        paste0(
          gpl_id_S4H,
          " gene-symbol annotation"
        ),
        length(nmd3_probe_ids_S4H)
      ),
      rep(
        "Entrez Gene ID 51068",
        length(nmd3_rows_MSK_S4H)
      )
    ),
    stringsAsFactors = FALSE
  )
  
  utils::write.csv(
    probe_information_S4H,
    file =
      "FigureS4H_NMD3_probe_mapping.csv",
    row.names = FALSE
  )
  
  cat(
    "\nFigure S4H completed and saved.\n"
  )
  
  
  # ============================================================
  
  
  
}))

# =============================================================================
# FigS2A
# =============================================================================
nmd3_step("FigS2A", c("Prepare_TCGA"), quote({
  # FigS2A — three horizontal GSEA plots
  # TCGA-PRAD primary tumors
  # NMD3-high vs NMD3-low GSEA
  #
  # Left:   KEGG prostate cancer
  # Middle: Reactome translation
  # Right:  Ribosomal large subunit export from nucleus
  #
  # Reuses primary_counts generated in Figure 1D
  # ============================================================
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tibble)
    library(msigdbr)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(patchwork)
  })
  
  stopifnot(
    exists("primary_counts"),
    "NMD3" %in% rownames(primary_counts),
    ncol(primary_counts) > 100
  )
  
  cat(
    "TCGA primary tumor samples:",
    ncol(primary_counts),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 1. Normalize NMD3 expression
  # ------------------------------------------------------------
  
  coldata_norm_S4B <- data.frame(
    row.names = colnames(primary_counts)
  )
  
  dds_norm_S4B <- DESeq2::DESeqDataSetFromMatrix(
    countData = primary_counts,
    colData = coldata_norm_S4B,
    design = ~ 1
  )
  
  dds_norm_S4B <- DESeq2::estimateSizeFactors(
    dds_norm_S4B
  )
  
  nmd3_expression_S4B <- DESeq2::counts(
    dds_norm_S4B,
    normalized = TRUE
  )["NMD3", ]
  
  nmd3_median_S4B <- median(
    nmd3_expression_S4B,
    na.rm = TRUE
  )
  
  nmd3_group_S4B <- factor(
    ifelse(
      nmd3_expression_S4B >
        nmd3_median_S4B,
      "NMD3_high",
      "NMD3_low"
    ),
    levels = c(
      "NMD3_low",
      "NMD3_high"
    )
  )
  
  cat(
    "NMD3 median normalized expression:",
    nmd3_median_S4B,
    "\n"
  )
  
  cat("NMD3 groups:\n")
  print(table(nmd3_group_S4B))
  
  # ------------------------------------------------------------
  # 2. Differential expression
  # ------------------------------------------------------------
  
  coldata_S4B <- data.frame(
    group = nmd3_group_S4B,
    row.names = names(nmd3_group_S4B)
  )
  
  dds_S4B <- DESeq2::DESeqDataSetFromMatrix(
    countData = primary_counts,
    colData = coldata_S4B,
    design = ~ group
  )
  
  dds_S4B <- DESeq2::DESeq(
    dds_S4B,
    quiet = TRUE
  )
  
  res_nmd3_S4B <- DESeq2::results(
    dds_S4B,
    contrast = c(
      "group",
      "NMD3_high",
      "NMD3_low"
    )
  )
  
  res_nmd3_df_S4B <- as.data.frame(
    res_nmd3_S4B
  ) %>%
    tibble::rownames_to_column(
      "gene"
    ) %>%
    dplyr::arrange(
      padj
    )
  
  # Positive Wald statistic means NMD3-high > NMD3-low
  
  gene_ranking_nmd3_S4B <- res_nmd3_S4B$stat
  
  names(gene_ranking_nmd3_S4B) <- toupper(
    rownames(res_nmd3_S4B)
  )
  
  gene_ranking_nmd3_S4B <- gene_ranking_nmd3_S4B[
    is.finite(gene_ranking_nmd3_S4B)
  ]
  
  gene_ranking_nmd3_S4B <- sort(
    gene_ranking_nmd3_S4B,
    decreasing = TRUE
  )
  
  cat(
    "Genes used for Figure S4B GSEA:",
    length(gene_ranking_nmd3_S4B),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 3. Retrieve pathway gene sets
  # ------------------------------------------------------------
  
  kegg_sets_S4B <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY"
  )
  
  reactome_sets_S4B <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:REACTOME"
  )
  
  gobp_sets_S4B <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  )
  
  term_prostate_S4B <- kegg_sets_S4B %>%
    dplyr::filter(
      gs_exact_source == "hsa05215"
    ) %>%
    dplyr::transmute(
      term = "KEGG_PROSTATE_CANCER",
      gene = toupper(gene_symbol)
    ) %>%
    dplyr::distinct()
  
  term_translation_S4B <- reactome_sets_S4B %>%
    dplyr::filter(
      gs_name == "REACTOME_TRANSLATION"
    ) %>%
    dplyr::transmute(
      term = "REACTOME_TRANSLATION",
      gene = toupper(gene_symbol)
    ) %>%
    dplyr::distinct()
  
  term_large_export_S4B <- gobp_sets_S4B %>%
    dplyr::filter(
      gs_exact_source == "GO:0000055"
    ) %>%
    dplyr::transmute(
      term =
        "GOBP_RIBOSOMAL_LARGE_SUBUNIT_EXPORT_FROM_NUCLEUS",
      gene = toupper(gene_symbol)
    ) %>%
    dplyr::distinct()
  
  stopifnot(
    nrow(term_prostate_S4B) > 0,
    nrow(term_translation_S4B) > 0,
    nrow(term_large_export_S4B) > 0
  )
  
  cat(
    "KEGG prostate cancer genes:",
    nrow(term_prostate_S4B),
    "\n"
  )
  
  cat(
    "Reactome translation genes:",
    nrow(term_translation_S4B),
    "\n"
  )
  
  cat(
    "Large subunit export genes:",
    nrow(term_large_export_S4B),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 4. KEGG prostate cancer GSEA
  # ------------------------------------------------------------
  
  set.seed(123)
  
  gsea_prostate_S4B <- clusterProfiler::GSEA(
    geneList = gene_ranking_nmd3_S4B,
    TERM2GENE = term_prostate_S4B,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE,
    seed = TRUE
  )
  
  # ------------------------------------------------------------
  # 5. Reactome translation GSEA
  # ------------------------------------------------------------
  
  set.seed(123)
  
  gsea_translation_S4B <- clusterProfiler::GSEA(
    geneList = gene_ranking_nmd3_S4B,
    TERM2GENE = term_translation_S4B,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE,
    seed = TRUE
  )
  
  # ------------------------------------------------------------
  # 6. Large subunit export GSEA
  # ------------------------------------------------------------
  
  set.seed(123)
  
  gsea_large_export_S4B <- clusterProfiler::GSEA(
    geneList = gene_ranking_nmd3_S4B,
    TERM2GENE = term_large_export_S4B,
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE,
    seed = TRUE
  )
  
  # ------------------------------------------------------------
  # 7. Print results
  # ------------------------------------------------------------
  
  cat("\nKEGG prostate cancer:\n")
  
  print(
    gsea_prostate_S4B@result[
      ,
      c(
        "ID",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
  cat("\nReactome translation:\n")
  
  print(
    gsea_translation_S4B@result[
      ,
      c(
        "ID",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
  cat("\nRibosomal large subunit export:\n")
  
  print(
    gsea_large_export_S4B@result[
      ,
      c(
        "ID",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
  # ------------------------------------------------------------
  # 8. Individual GSEA plots
  # ------------------------------------------------------------
  
  p_S4B_prostate <- enrichplot::gseaplot2(
    gsea_prostate_S4B,
    geneSetID = "KEGG_PROSTATE_CANCER",
    subplots = 1:3,
    pvalue_table = FALSE,
    title = "PROSTATE CANCER",
    base_size = 11
  )
  
  p_S4B_translation <- enrichplot::gseaplot2(
    gsea_translation_S4B,
    geneSetID = "REACTOME_TRANSLATION",
    subplots = 1:3,
    pvalue_table = FALSE,
    title = "TRANSLATION",
    base_size = 11
  )
  
  p_S4B_large_export <- enrichplot::gseaplot2(
    gsea_large_export_S4B,
    geneSetID =
      "GOBP_RIBOSOMAL_LARGE_SUBUNIT_EXPORT_FROM_NUCLEUS",
    subplots = 1:3,
    pvalue_table = FALSE,
    title = "LARGE SUBUNIT EXPORT FROM NUCLEUS",
    base_size = 11
  )
  
  
  
  
  # ------------------------------------------------------------
  # 9. Combine Figure S4B in the requested order
  # ------------------------------------------------------------
  
  # Capture the FULL three-subplot gseaplot2 object as a grid grob.
  # gseaplot2() may return an aplot gglist, which patchwork::wrap_plots()
  # cannot reliably nest across package versions. Capture each complete
  # graph, then place the three full graphs side-by-side on ONE PDF page.
  .s2a_capture <- function(p) {
    g <- grid::grid.grabExpr(print(p), wrap=TRUE)
    if (is.null(g) || !inherits(g, "grob")) {
      stop("FigS2A could not capture a complete GSEA plot")
    }
    g
  }
  p_S4B <- gridExtra::arrangeGrob(
    grobs = list(
      .s2a_capture(p_S4B_prostate),
      .s2a_capture(p_S4B_translation),
      .s2a_capture(p_S4B_large_export)
    ),
    ncol = 3
  )

  # ------------------------------------------------------------
  # 10. Save combined Figure S4B
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS2A_NMD3_High_vs_Low_Three_GSEA.pdf",
    plot = p_S4B,
    width = 16,
    height = 5.7,
    units = "in",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  
  # ------------------------------------------------------------
  # 11. Additional individual PDFs intentionally not generated
  # ------------------------------------------------------------
  
  # ------------------------------------------------------------
  # 12. Save differential-expression and GSEA results
  # ------------------------------------------------------------
  
  utils::write.csv(
    res_nmd3_df_S4B,
    file = "FigS2A_TCGA_NMD3_High_vs_Low_DESeq2.csv",
    row.names = FALSE
  )
  
  gsea_results_S4B <- dplyr::bind_rows(
    gsea_prostate_S4B@result %>%
      dplyr::mutate(
        Panel = "1_KEGG_prostate_cancer"
      ),
    gsea_translation_S4B@result %>%
      dplyr::mutate(
        Panel = "2_Reactome_translation"
      ),
    gsea_large_export_S4B@result %>%
      dplyr::mutate(
        Panel = "3_Large_subunit_export"
      )
  )
  
  utils::write.csv(
    gsea_results_S4B,
    file = "FigS2A_TCGA_NMD3_High_vs_Low_GSEA_results.csv",
    row.names = FALSE
  )
  
  
  
}))

# =============================================================================
# FigS2B
# =============================================================================
nmd3_step("FigS2B", character(), quote({
  # FigS2B — analysis below
  # Spearman correlation between XPO1 and NMD3
  #
  # TCGA-PRAD:
  #   TCGA RNASeqV2 RSEM-normalized expression
  #   obtained directly from cBioPortal API (study: prad_tcga)
  #   Primary tumors only, one value per patient
  #
  # MSKCC:
  #   Clinically annotated tumor samples only
  #   One independent observation per patient
  #   (unchanged from original script)
  # ============================================================
  
  
  suppressPackageStartupMessages({
    library(cBioPortalData)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
  })
  
  # ============================================================
  # 1. TCGA-PRAD primary tumors
  # RSEM-normalized expression via cBioPortal API,
  # patient-level aggregation
  # ============================================================
  
  study_id_S4D <- "prad_tcga"
  genes_S4D <- c("XPO1", "NMD3")
  
  cbio_S4D <- cBioPortal()
  
  cat("\nConnected to cBioPortal API.\n")
  
  # ------------------------------------------------------------
  # Confirm study
  # ------------------------------------------------------------
  
  studies_S4D <- getStudies(cbio_S4D)
  
  study_check_S4D <- studies_S4D %>%
    dplyr::filter(studyId == study_id_S4D)
  
  if (nrow(study_check_S4D) == 0) {
    stop("Study prad_tcga was not found in cBioPortal.")
  }
  
  cat("\nStudy:\n")
  print(study_check_S4D %>% dplyr::select(studyId, name, allSampleCount))
  
  # ------------------------------------------------------------
  # Identify RNASeqV2 RSEM molecular profile (non-Z-score)
  # ------------------------------------------------------------
  
  profiles_S4D <- molecularProfiles(cbio_S4D, study_id_S4D)
  
  rsem_profile_candidates_S4D <- profiles_S4D$molecularProfileId[
    grepl("rna_seq_v2_mrna$", profiles_S4D$molecularProfileId, ignore.case = TRUE) &
      !grepl("zscore", profiles_S4D$molecularProfileId, ignore.case = TRUE)
  ]
  
  if (length(rsem_profile_candidates_S4D) == 0) {
    stop("Could not identify the non-Z-score RNASeqV2 RSEM expression profile.")
  }
  
  rsem_profile_S4D <- rsem_profile_candidates_S4D[1]
  
  cat("\nRNA expression profile selected:\n", rsem_profile_S4D, "\n", sep = "")
  
  # ------------------------------------------------------------
  # Download XPO1 and NMD3 RSEM values
  # ------------------------------------------------------------
  
  expression_list_S4D <- getDataByGenes(
    api = cbio_S4D,
    studyId = study_id_S4D,
    genes = genes_S4D,
    by = "hugoGeneSymbol",
    molecularProfileIds = rsem_profile_S4D
  )
  
  stopifnot(rsem_profile_S4D %in% names(expression_list_S4D))
  
  expression_long_S4D <- expression_list_S4D[[rsem_profile_S4D]]
  
  stopifnot(
    all(c("sampleId", "value", "hugoGeneSymbol") %in% colnames(expression_long_S4D))
  )
  
  expression_long_S4D <- expression_long_S4D %>%
    dplyr::mutate(
      hugoGeneSymbol = toupper(as.character(hugoGeneSymbol)),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    dplyr::filter(
      hugoGeneSymbol %in% genes_S4D,
      is.finite(value)
    )
  
  cat("\nGenes returned:\n")
  print(table(expression_long_S4D$hugoGeneSymbol))
  
  # ------------------------------------------------------------
  # One row per sample
  # ------------------------------------------------------------
  
  expression_wide_S4D <- expression_long_S4D %>%
    dplyr::select(sampleId, patientId, hugoGeneSymbol, value) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = hugoGeneSymbol, values_from = value)
  
  stopifnot(all(c("XPO1", "NMD3") %in% colnames(expression_wide_S4D)))
  
  expression_wide_S4D <- expression_wide_S4D %>%
    dplyr::filter(is.finite(XPO1), is.finite(NMD3))
  
  cat("\nSamples with both XPO1 and NMD3:", nrow(expression_wide_S4D), "\n")
  
  # ------------------------------------------------------------
  # Keep primary tumor samples only
  # TCGA sample type code 01 = Primary Solid Tumor
  # ------------------------------------------------------------
  
  expression_wide_S4D <- expression_wide_S4D %>%
    dplyr::mutate(Sample_Type_Code = substr(sampleId, 14, 15))
  
  cat("\nTCGA sample-type codes:\n")
  print(table(expression_wide_S4D$Sample_Type_Code, useNA = "ifany"))
  
  primary_S4D <- expression_wide_S4D %>%
    dplyr::filter(Sample_Type_Code == "01")
  
  cat("\nPrimary tumor samples:", nrow(primary_S4D), "\n")
  cat("Unique patients:", dplyr::n_distinct(primary_S4D$patientId), "\n")
  
  # ------------------------------------------------------------
  # Check for patients with >1 primary tumor sample
  # ------------------------------------------------------------
  
  duplicated_tcga_patients_S4D <- primary_S4D %>%
    dplyr::count(patientId, name = "Number_of_samples") %>%
    dplyr::filter(Number_of_samples > 1)
  
  cat(
    "TCGA patients represented by multiple samples:",
    nrow(duplicated_tcga_patients_S4D),
    "\n"
  )
  
  if (nrow(duplicated_tcga_patients_S4D) > 0) {
    print(duplicated_tcga_patients_S4D)
  }
  
  # ------------------------------------------------------------
  # log2(RSEM + 1) and aggregate to patient level
  # (mean of log2 values, same aggregation logic as before;
  #  patients with one sample remain unchanged)
  # ------------------------------------------------------------
  
  tcga_sample_level_S4D <- primary_S4D %>%
    dplyr::transmute(
      Sample = sampleId,
      Patient_ID = patientId,
      XPO1 = log2(XPO1 + 1),
      NMD3 = log2(NMD3 + 1)
    ) %>%
    dplyr::filter(is.finite(XPO1), is.finite(NMD3))
  
  tcga_cor_S4D <- tcga_sample_level_S4D %>%
    dplyr::group_by(Patient_ID) %>%
    dplyr::summarise(
      XPO1 = mean(XPO1, na.rm = TRUE),
      NMD3 = mean(NMD3, na.rm = TRUE),
      Number_of_samples = dplyr::n(),
      Samples = paste(Sample, collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(XPO1), is.finite(NMD3))
  
  stopifnot(
    nrow(tcga_cor_S4D) > 10,
    !anyDuplicated(tcga_cor_S4D$Patient_ID)
  )
  
  cat("Independent TCGA patients after aggregation:", nrow(tcga_cor_S4D), "\n")
  
  # ------------------------------------------------------------
  # Spearman correlation
  # ------------------------------------------------------------
  
  cor_tcga_S4D <- stats::cor.test(
    tcga_cor_S4D$XPO1,
    tcga_cor_S4D$NMD3,
    method = "spearman",
    exact = FALSE
  )
  
  cat("\nTCGA-PRAD patient-level result (RSEM normalized):\n")
  cat("Spearman rho =", round(unname(cor_tcga_S4D$estimate), 4), "\n")
  cat("P value =", format(cor_tcga_S4D$p.value, scientific = TRUE, digits = 6), "\n")
  cat("Independent patients =", nrow(tcga_cor_S4D), "\n")
  
  # ============================================================
  # 2. MSKCC prostate cancer
  # (unchanged from original script)
  # ============================================================
  
  msk_expression_file_S4D <- paste0(
    "data/",
    "prad_mskcc/",
    "data_mrna_agilent_microarray.txt"
  )
  
  msk_clinical_file_S4D <- paste0(
    "data/",
    "prad_mskcc/",
    "data_clinical_sample.txt"
  )
  
  msk_raw_S4D <- data.table::fread(
    msk_expression_file_S4D,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
  )
  
  stopifnot("Entrez_Gene_Id" %in% colnames(msk_raw_S4D))
  
  msk_clinical_S4D <- read.delim(
    msk_clinical_file_S4D,
    sep = "\t",
    header = TRUE,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  stopifnot(all(c("SAMPLE_ID", "PATIENT_ID") %in% colnames(msk_clinical_S4D)))
  
  sample_type_candidates_S4D <- intersect(
    c("SAMPLE_TYPE", "SAMPLE_CLASS", "TUMOR_TYPE", "DISEASE_CODE"),
    colnames(msk_clinical_S4D)
  )
  
  if (length(sample_type_candidates_S4D) == 0) {
    stop("No sample-type column was found in the MSKCC clinical sample table.")
  }
  
  sample_type_column_S4D <- sample_type_candidates_S4D[1]
  
  cat("\nMSKCC sample-type column:", sample_type_column_S4D, "\n")
  cat("Original MSKCC sample types:\n")
  print(table(msk_clinical_S4D[[sample_type_column_S4D]], useNA = "ifany"))
  
  msk_sample_map_S4D <- msk_clinical_S4D %>%
    dplyr::transmute(
      Sample = trimws(as.character(SAMPLE_ID)),
      Patient_ID = trimws(as.character(PATIENT_ID)),
      Sample_Type = trimws(as.character(.data[[sample_type_column_S4D]]))
    ) %>%
    dplyr::filter(
      !is.na(Sample),
      Sample != "",
      !is.na(Patient_ID),
      Patient_ID != "",
      Sample %in% colnames(msk_raw_S4D),
      !is.na(Sample_Type),
      Sample_Type != "",
      !grepl("normal|benign|cell.?line", Sample_Type, ignore.case = TRUE)
    ) %>%
    dplyr::distinct(Sample, .keep_all = TRUE)
  
  known_cell_lines_S4D <- c("LNCaP104S", "LNCaP104R", "LNCaP", "DU145", "PC3")
  
  msk_sample_map_S4D <- msk_sample_map_S4D %>%
    dplyr::filter(!Sample %in% known_cell_lines_S4D)
  
  stopifnot(
    nrow(msk_sample_map_S4D) > 10,
    !anyDuplicated(msk_sample_map_S4D$Sample)
  )
  
  cat("MSKCC tumor samples retained:", nrow(msk_sample_map_S4D), "\n")
  cat("Retained sample types:\n")
  print(table(msk_sample_map_S4D$Sample_Type))
  
  # ============================================================
  # 3. Extract MSKCC XPO1 and NMD3 expression
  # ============================================================
  
  xpo1_rows_S4D <- which(as.character(msk_raw_S4D$Entrez_Gene_Id) == "7514")
  nmd3_rows_S4D <- which(as.character(msk_raw_S4D$Entrez_Gene_Id) == "51068")
  
  stopifnot(length(xpo1_rows_S4D) > 0, length(nmd3_rows_S4D) > 0)
  
  extract_msk_gene_S4D <- function(data, row_index, sample_columns) {
    gene_matrix <- data[row_index, sample_columns, drop = FALSE]
    gene_matrix[] <- lapply(gene_matrix, function(x) {
      suppressWarnings(as.numeric(as.character(x)))
    })
    gene_matrix <- as.matrix(gene_matrix)
    if (nrow(gene_matrix) == 1) {
      values <- as.numeric(gene_matrix[1, ])
    } else {
      values <- colMeans(gene_matrix, na.rm = TRUE)
    }
    names(values) <- sample_columns
    values
  }
  
  msk_samples_S4D <- msk_sample_map_S4D$Sample
  
  xpo1_msk_S4D <- extract_msk_gene_S4D(msk_raw_S4D, xpo1_rows_S4D, msk_samples_S4D)
  nmd3_msk_S4D <- extract_msk_gene_S4D(msk_raw_S4D, nmd3_rows_S4D, msk_samples_S4D)
  
  msk_sample_expression_S4D <- msk_sample_map_S4D %>%
    dplyr::mutate(
      XPO1 = as.numeric(xpo1_msk_S4D[Sample]),
      NMD3 = as.numeric(nmd3_msk_S4D[Sample])
    ) %>%
    dplyr::filter(is.finite(XPO1), is.finite(NMD3))
  
  # ============================================================
  # 4. Aggregate MSKCC samples to patient level
  # ============================================================
  
  msk_cor_S4D <- msk_sample_expression_S4D %>%
    dplyr::group_by(Patient_ID) %>%
    dplyr::summarise(
      XPO1 = mean(XPO1, na.rm = TRUE),
      NMD3 = mean(NMD3, na.rm = TRUE),
      Number_of_samples = dplyr::n(),
      Sample_types = paste(sort(unique(Sample_Type)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(XPO1), is.finite(NMD3))
  
  stopifnot(
    nrow(msk_cor_S4D) > 10,
    !anyDuplicated(msk_cor_S4D$Patient_ID)
  )
  
  cat("\nMSKCC independent patients:", nrow(msk_cor_S4D), "\n")
  cat(
    "Patients represented by multiple tumor samples:",
    sum(msk_cor_S4D$Number_of_samples > 1),
    "\n"
  )
  
  cor_msk_S4D <- stats::cor.test(
    msk_cor_S4D$XPO1,
    msk_cor_S4D$NMD3,
    method = "spearman",
    exact = FALSE
  )
  
  cat("\nMSKCC patient-level correlation:\n")
  cat("Spearman rho =", round(unname(cor_msk_S4D$estimate), 4), "\n")
  cat("P value =", format(cor_msk_S4D$p.value, scientific = TRUE, digits = 6), "\n")
  
  # ============================================================
  # 5. Correlation labels
  # ============================================================
  
  format_p_S4D <- function(p_value) {
    if (p_value < .Machine$double.eps) {
      return("p < 2.2e-16")
    }
    paste0("p = ", formatC(p_value, format = "e", digits = 2))
  }
  
  label_tcga_S4D <- paste0(
    "rho = ", sprintf("%.2f", unname(cor_tcga_S4D$estimate)),
    "\n", format_p_S4D(cor_tcga_S4D$p.value),
    "\nn = ", nrow(tcga_cor_S4D)
  )
  
  label_msk_S4D <- paste0(
    "rho = ", sprintf("%.2f", unname(cor_msk_S4D$estimate)),
    "\n", format_p_S4D(cor_msk_S4D$p.value),
    "\nn = ", nrow(msk_cor_S4D)
  )
  
  # ============================================================
  # 6. Plotting function (unchanged)
  # ============================================================
  
  make_correlation_plot_S4D <- function(
    plot_data,
    label_text,
    panel_title,
    x_label,
    y_label
  ) {
    x_position <- min(plot_data$XPO1, na.rm = TRUE)
    y_position <- max(plot_data$NMD3, na.rm = TRUE)
    
    ggplot2::ggplot(plot_data, ggplot2::aes(x = XPO1, y = NMD3)) +
      ggplot2::geom_point(
        shape = 1, color = "#1F78B4", size = 2.3, alpha = 0.70, stroke = 0.8
      ) +
      ggplot2::geom_smooth(
        method = "lm", formula = y ~ x, color = "black", fill = "grey80",
        se = TRUE, linewidth = 0.8, alpha = 0.50
      ) +
      ggplot2::annotate(
        "text", x = x_position, y = y_position, label = label_text,
        hjust = 0, vjust = 1, size = 4.5, color = "black"
      ) +
      ggplot2::labs(title = panel_title, x = x_label, y = y_label) +
      ggplot2::theme_classic(base_size = 13, base_family = "Arial") +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
        axis.text = ggplot2::element_text(color = "black", size = 11),
        axis.title = ggplot2::element_text(color = "black", size = 13, face = "bold"),
        plot.title = ggplot2::element_text(color = "black", size = 14, face = "bold", hjust = 0.5)
      )
  }
  
  # ============================================================
  # 7. Draw Figure S4D
  # ============================================================
  
  p_S4D_TCGA <- make_correlation_plot_S4D(
    plot_data = tcga_cor_S4D,
    label_text = label_tcga_S4D,
    panel_title = "TCGA-PRAD primary tumors",
    x_label = expression(
      italic(XPO1) ~ log[2] * "(RSEM normalized expression + 1)"
    ),
    y_label = expression(
      italic(NMD3) ~ log[2] * "(RSEM normalized expression + 1)"
    )
  )
  
  p_S4D_MSK <- make_correlation_plot_S4D(
    plot_data = msk_cor_S4D,
    label_text = label_msk_S4D,
    panel_title = "MSKCC prostate cancer",
    x_label = expression(italic(XPO1) ~ expression),
    y_label = expression(italic(NMD3) ~ expression)
  )
  
  p_S4D <- p_S4D_TCGA + p_S4D_MSK + patchwork::plot_layout(nrow = 1)
  
  
  
  
  # ============================================================
  # 8. Save PDFs
  # ============================================================
  
  nmd3_save_pdf(
    filename = "FigS2B_XPO1_NMD3_Spearman_TCGA_MSKCC.pdf",
    plot = p_S4D,
    width = 10, height = 4.8, units = "in",
    device = grDevices::cairo_pdf
  )
  
  # ============================================================
  # 9. Save source data
  # ============================================================
  
  utils::write.csv(
    tcga_cor_S4D,
    file = "FigureS4D_TCGA_XPO1_NMD3_patient_data.csv",
    row.names = FALSE
  )
  
  utils::write.csv(
    msk_cor_S4D,
    file = "FigureS4D_MSKCC_XPO1_NMD3_patient_data.csv",
    row.names = FALSE
  )
  
}))

# =============================================================================
# FigS2C_E
# =============================================================================
nmd3_step("FigS2C_E", c("Fig6A"), quote({
  # =============================================================================
  # PART 5: PUBLIC PROTEOMICS RB/RPL ADJUSTED NMD3; REUSE PUBLIC MATRIX
  # =============================================================================
  
  
  # NEW PUBLIC-PROTEOMICS SENSITIVITY (ALL / RB-adjusted / RPL-adjusted)
  # Shares the already-loaded 58-specimen public matrix, metadata, and NMD3 values.
  stopifnot(exists("all_protein_matrix"), exists("sample_metadata"),
            exists("nmd3_data"), exists("go_bp_sets"))
  output_dir <- file.path(.nmd3_run_root, "extra_sensitivity", "Public_NMD3_RB_RPL_Adjusted")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  # 3. GENE SETS: RB is MSigDB C5 GO:BP GO:0042254 (without NMD3).
  # RPL is a user-defined full measured protein family: symbol starts with RPL.
  msig_args <- if ("collection" %in% names(formals(msigdbr::msigdbr))) {
    list(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  } else {
    list(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  }
  gobp <- go_bp_sets
  rb_mask <- !is.na(gobp$gs_name) & gobp$gs_name == "GOBP_RIBOSOME_BIOGENESIS"
  if ("gs_exact_source" %in% colnames(gobp)) {
    rb_mask <- rb_mask | (!is.na(gobp$gs_exact_source) &
                            gobp$gs_exact_source == "GO:0042254")
  }
  rb_gene_set <- setdiff(sort(unique(toupper(trimws(
    as.character(gobp$gene_symbol[rb_mask])
  )))), "NMD3")
  rb_gene_set <- rb_gene_set[!is.na(rb_gene_set) & nzchar(rb_gene_set)]
  rpl_gene_set <- setdiff(sort(unique(
    rownames(all_protein_matrix)[grepl("^RPL", rownames(all_protein_matrix))]
  )), "NMD3")
  stopifnot(length(rb_gene_set) >= 30L, length(rpl_gene_set) >= 5L)
  
  # 4. SCORE — match original exactly: each protein is row-Z-scaled across
  # ALL 58 specimens (8 Normal, 28 Primary, 22 Metastatic); ONLY proteins with
  # nonmissing data in all 58 and nonzero variance enter each score.
  make_protein_score <- function(mat, gene_set, name) {
    represented <- intersect(gene_set, rownames(mat))
    complete <- represented[
      rowSums(is.finite(mat[represented, , drop = FALSE])) == ncol(mat)
    ]
    genes <- complete[vapply(complete, function(g) {
      ss <- stats::sd(as.numeric(mat[g, ]))
      is.finite(ss) && ss > 0
    }, logical(1))]
    if (length(genes) < 5L) stop(name, " score: fewer than five usable proteins")
    zmat <- t(scale(t(mat[genes, , drop = FALSE])))
    stopifnot(is.matrix(zmat), all(is.finite(zmat)))
    score <- colMeans(zmat)
    names(score) <- colnames(mat)
    utils::write.csv(
      data.frame(Gene = gene_set,
                 Represented = gene_set %in% represented,
                 Complete_58 = gene_set %in% complete,
                 Used = gene_set %in% genes),
      file.path(output_dir, paste0(name, "_gene_coverage.csv")),
      row.names = FALSE
    )
    message(name, ": ", length(genes), " proteins used")
    score
  }
  rb_score <- make_protein_score(all_protein_matrix, rb_gene_set, "RB")
  rpl_score <- make_protein_score(all_protein_matrix, rpl_gene_set, "RPL")
  
  # 5. FIXED COMPARISON GROUP: 28 Primary + 20 Metastatic (48 total).
  plot_data <- nmd3_data %>%
    dplyr::mutate(Sample = as.character(Sample), Group = as.character(Group)) %>%
    dplyr::filter(Group %in% c("Primary", "Metastatic"), is.finite(Value)) %>%
    dplyr::transmute(
      Sample,
      Group = factor(Group, levels = c("Primary", "Metastatic")),
      NMD3 = as.numeric(Value),
      RB = as.numeric(rb_score[Sample]),
      RPL = as.numeric(rpl_score[Sample])
    )
  stopifnot(nrow(plot_data) == 48L,
            sum(plot_data$Group == "Primary") == 28L,
            sum(plot_data$Group == "Metastatic") == 20L,
            all(is.finite(as.matrix(plot_data[, c("NMD3", "RB", "RPL")]))))
  
  z <- function(x) {
    stopifnot(all(is.finite(x)), stats::sd(x) > 0)
    as.numeric(scale(x))
  }
  plot_data$Metastatic <- as.integer(plot_data$Group == "Metastatic")
  plot_data$RB_z  <- z(plot_data$RB)
  plot_data$RPL_z <- z(plot_data$RPL)
  
  # Original plot: two-sided Wilcoxon (same as manuscript figure S4A).
  original_test <- stats::wilcox.test(
    NMD3 ~ Group, data = plot_data,
    alternative = "two.sided", exact = FALSE
  )
  
  # Adjusted P values: group coefficient from the covariate-adjusted model.
  # These are NOT P values from a t-test/Wilcoxon on the residuals.
  model_rb <- stats::lm(NMD3 ~ Metastatic + RB_z, data = plot_data)
  model_rpl <- stats::lm(NMD3 ~ Metastatic + RPL_z, data = plot_data)
  p_original <- original_test$p.value
  p_rb <- summary(model_rb)$coefficients["Metastatic", "Pr(>|t|)"]
  p_rpl <- summary(model_rpl)$coefficients["Metastatic", "Pr(>|t|)"]
  
  # Data plotted in adjusted panels: score-only residual, shifted back by the
  # overall original mean so all three plots share an interpretable center.
  plot_data$Y1  <- plot_data$NMD3
  plot_data$Y11 <- as.numeric(residuals(lm(NMD3 ~ RB_z, data = plot_data))) +
    mean(plot_data$NMD3)
  plot_data$Y12 <- as.numeric(residuals(lm(NMD3 ~ RPL_z, data = plot_data))) +
    mean(plot_data$NMD3)
  
  stats_out <- data.frame(
    Figure = c("1", "1.1", "1.2"),
    Comparison = "Primary vs Metastatic",
    Method = c("Wilcoxon, two-sided",
               "LM NMD3 ~ Metastatic + RB_z",
               "LM NMD3 ~ Metastatic + RPL_z"),
    Primary_n = 28L,
    Metastatic_n = 20L,
    P_value = c(p_original, p_rb, p_rpl),
    stringsAsFactors = FALSE
  )
  print(stats_out, row.names = FALSE)
  
  # 6. DRAW VIOLINS — original colors, same boxplot/jitter/bracket/theme.
  # *** ONLY CHANGE from previous script: geom_point(color="black") ***
  violin_colors <- c(Primary = "#4DBBD5", Metastatic = "#E64B35FF")
  p_label <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", p))
  }
  
  plot_violin <- function(ycol, ylabel, pv, filename) {
    dd <- plot_data
    dd$Value <- as.numeric(dd[[ycol]])
    rr <- range(dd$Value, finite = TRUE)
    span <- diff(rr)
    if (!is.finite(span) || span <= 0) span <- 1
    by <- rr[2] + 0.58 * span
    ty <- by + 0.04 * span
    
    p <- ggplot2::ggplot(dd, ggplot2::aes(x = Group, y = Value, fill = Group)) +
      ggplot2::geom_violin(trim = FALSE, scale = "width", width = 0.80,
                           color = NA, alpha = 0.65) +
      ggplot2::geom_boxplot(width = 0.28, outlier.shape = NA,
                            fill = "white", color = "black", linewidth = 0.65) +
      ggplot2::geom_point(
        color = "black",  # <-- requested change: every individual dot is black
        position = ggplot2::position_jitter(width = 0.13, height = 0, seed = 123),
        size = 1.7, alpha = 0.80, show.legend = FALSE
      ) +
      ggplot2::annotate("segment", x = 1, xend = 2, y = by, yend = by,
                        linewidth = 0.6, color = "black") +
      ggplot2::annotate("segment", x = 1, xend = 1, y = by,
                        yend = by - 0.025 * span, linewidth = 0.6,
                        color = "black") +
      ggplot2::annotate("segment", x = 2, xend = 2, y = by,
                        yend = by - 0.025 * span, linewidth = 0.6,
                        color = "black") +
      ggplot2::annotate("text", x = 1.5, y = ty,
                        label = p_label(pv), size = 4.2,
                        fontface = "bold", color = "black") +
      ggplot2::scale_fill_manual(values = violin_colors, drop = FALSE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.20))) +
      ggplot2::labs(x = NULL, y = ylabel) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA,
                                             linewidth = 0.8),
        axis.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = ggplot2::element_text(size = 11, color = "black"),
        axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "black"),
        legend.position = "none",
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    nmd3_save_pdf(file.path(output_dir, filename), p,
                  width = 4.6, height = 5.4, units = "in",
                  device = grDevices::cairo_pdf)
    invisible(p)
  }
  
  # Unadjusted NMD3 violin exists in the main public-proteomics workflow;
  # it is NOT plotted again here. Adjusted panels below use the same raw input.
  p_1 <- plot_violin("Y1", "NMD3 protein abundance", p_original,
                     "FigS2C_NMD3_Protein_Unadjusted.pdf")
  p_11 <- plot_violin("Y11", "RB-adjusted NMD3 protein abundance", p_rb,
                      "FigS2D_RB_Adjusted_NMD3_Protein_BlackPoints.pdf")
  p_12 <- plot_violin("Y12", "RPL-adjusted NMD3 protein abundance", p_rpl,
                      "FigS2E_RPL_Adjusted_NMD3_Protein_BlackPoints.pdf")
  
  utils::write.csv(plot_data,
                   file.path(output_dir, "3_Plots_48_Samples_Source_Data.csv"),
                   row.names = FALSE)
  utils::write.csv(stats_out,
                   file.path(output_dir, "3_Plots_Statistics.csv"),
                   row.names = FALSE)
  utils::write.csv(data.frame(Gene = rb_gene_set),
                   file.path(output_dir, "RB_GO0042254_Gene_Set.csv"),
                   row.names = FALSE)
  
  message("\nDONE. Three PDF figures and source/statistics CSVs saved in:\n",
          normalizePath(output_dir))
  
  
  
  
  
  
}))

# =============================================================================
# Prepare_XPO1_Groups
# =============================================================================
nmd3_step("Prepare_XPO1_Groups", c("Prepare_TCGA"), quote({
  primary_coldata <- data.frame(
    row.names = colnames(primary_counts)
  )
  
  dds_primary_norm <- DESeqDataSetFromMatrix(
    countData = primary_counts,
    colData = primary_coldata,
    design = ~ 1
  )
  
  dds_primary_norm <- estimateSizeFactors(
    dds_primary_norm
  )
  
  stopifnot(
    "XPO1" %in% rownames(dds_primary_norm)
  )
  
  xpo1_expression <- counts(
    dds_primary_norm,
    normalized = TRUE
  )["XPO1", ]
  
  xpo1_median <- median(
    xpo1_expression,
    na.rm = TRUE
  )
  
  xpo1_group <- factor(
    ifelse(
      xpo1_expression > xpo1_median,
      "XPO1_high",
      "XPO1_low"
    ),
    levels = c(
      "XPO1_low",
      "XPO1_high"
    )
  )
  
  primary_coldata <- data.frame(
    group = xpo1_group,
    row.names = colnames(primary_counts)
  )
  
}))

# =============================================================================
# FigS2F
# =============================================================================
nmd3_step("FigS2F", c("Prepare_XPO1_Groups"), quote({
  # FigS2F — analysis below
  # NMD3 expression in TCGA-PRAD XPO1-low vs XPO1-high tumors
  # Reuses objects generated in Figure 1D
  # ============================================================
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
  })
  
  stopifnot(
    exists("dds_primary_norm"),
    exists("xpo1_group"),
    "NMD3" %in% rownames(dds_primary_norm)
  )
  
  # ------------------------------------------------------------
  # 1. Extract normalized NMD3 expression
  # ------------------------------------------------------------
  
  normalized_counts_S4C <- DESeq2::counts(
    dds_primary_norm,
    normalized = TRUE
  )
  
  sample_ids_S4C <- colnames(normalized_counts_S4C)
  
  stopifnot(
    all(sample_ids_S4C %in% names(xpo1_group))
  )
  
  df_S4C <- data.frame(
    sample_id = sample_ids_S4C,
    XPO1_group = as.character(xpo1_group[sample_ids_S4C]),
    NMD3_expression = log2(
      normalized_counts_S4C["NMD3", sample_ids_S4C] + 1
    ),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      XPO1_group %in% c("XPO1_low", "XPO1_high"),
      is.finite(NMD3_expression)
    ) %>%
    dplyr::mutate(
      XPO1_group = factor(
        XPO1_group,
        levels = c("XPO1_low", "XPO1_high"),
        labels = c("XPO1-low", "XPO1-high")
      )
    )
  
  cat("Figure S4C sample numbers:\n")
  print(table(df_S4C$XPO1_group))
  
  # ------------------------------------------------------------
  # 2. Statistical comparison
  # ------------------------------------------------------------
  
  test_S4C <- stats::wilcox.test(
    NMD3_expression ~ XPO1_group,
    data = df_S4C,
    exact = FALSE,
    alternative = "two.sided"
  )
  
  cat(
    "Wilcoxon p-value:",
    format(test_S4C$p.value, scientific = TRUE, digits = 4),
    "\n"
  )
  
  p_label_S4C <- if (test_S4C$p.value < 0.001) {
    expression(italic(p) < 0.001)
  } else {
    bquote(
      italic(p) == .(
        format.pval(
          test_S4C$p.value,
          digits = 2,
          eps = 0.001
        )
      )
    )
  }
  
  # ------------------------------------------------------------
  # 3. Plot
  # ------------------------------------------------------------
  
  y_range_S4C <- range(df_S4C$NMD3_expression, na.rm = TRUE)
  y_span_S4C  <- diff(y_range_S4C)
  
  bracket_y_S4C <- y_range_S4C[2] + 0.08 * y_span_S4C
  tick_y_S4C    <- bracket_y_S4C - 0.05 * y_span_S4C
  label_y_S4C   <- bracket_y_S4C + 0.07 * y_span_S4C
  
  group_colors_S4C <- c(
    "XPO1-low"  = "#4DBBD5",
    "XPO1-high" = "#F8766D"
  )
  
  set.seed(123)
  
  p_S4C <- ggplot(
    df_S4C,
    aes(
      x = XPO1_group,
      y = NMD3_expression,
      fill = XPO1_group,
      color = XPO1_group
    )
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      width = 0.82,
      alpha = 0.60,
      color = NA
    ) +
    geom_boxplot(
      width = 0.32,
      outlier.shape = NA,
      fill = NA,
      color = "black",
      linewidth = 0.7
    ) +
    geom_jitter(
      width = 0.14,
      size = 1.25,
      alpha = 0.65
    ) +
    geom_segment(
      aes(
        x = 1,
        xend = 2,
        y = bracket_y_S4C,
        yend = bracket_y_S4C
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_segment(
      aes(
        x = 1,
        xend = 1,
        y = tick_y_S4C,
        yend = bracket_y_S4C
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_segment(
      aes(
        x = 2,
        xend = 2,
        y = tick_y_S4C,
        yend = bracket_y_S4C
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = 1.5,
      y = label_y_S4C,
      label = p_label_S4C,
      parse = FALSE,
      size = 5
    ) +
    scale_fill_manual(values = group_colors_S4C) +
    scale_color_manual(values = group_colors_S4C) +
    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.18))
    ) +
    labs(
      x = NULL,
      y = expression(italic(NMD3)~expression)
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.7
      ),
      axis.text.x = element_text(
        color = "black",
        size = 13,
        face = "bold.italic"
      ),
      axis.text.y = element_text(
        color = "black",
        size = 12
      ),
      axis.title.y = element_text(
        color = "black",
        size = 14,
        face = "bold"
      ),
      plot.margin = margin(12, 12, 12, 12)
    )
  
  
  # ------------------------------------------------------------
  # 4. Save
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS2F_XPO1_High_Low_NMD3_expression.pdf",
    plot = p_S4C,
    width = 4.2,
    height = 5.2,
    units = "in",
    device = cairo_pdf
  )
  
  
}))

# =============================================================================
# Prepare_SU2C
# =============================================================================
nmd3_step("Prepare_SU2C", character(), quote({
  # SHARED SU2C/PCF 2019 SURVIVAL DATA
  # Used by Figures 1G, 6K, S2G, and S4F
  #
  # The official polyA cohort contains 81 eligible samples from
  # 80 unique patients. Repeated samples from the same patient are
  # averaged before survival analysis so that each patient contributes
  # one survival record. The same 80-patient cohort is used for all
  # four SU2C survival panels.
  # ============================================================
  
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(survival)
    library(survminer)
    library(ggplot2)
  })
  
  # ------------------------------------------------------------
  # 1. Input files
  # ------------------------------------------------------------
  
  su2c_expression_file <- paste0(
    "data/prad_su2c_2019/",
    "expression_with_clinical.xlsx"
  )
  
  su2c_clinical_sample_file <- paste0(
    "data/prad_su2c_2019/",
    "data_clinical_sample.txt"
  )
  
  stopifnot(
    file.exists(su2c_expression_file),
    file.exists(su2c_clinical_sample_file)
  )
  
  # ------------------------------------------------------------
  # 2. Genes required across all four SU2C survival panels
  # ------------------------------------------------------------
  
  signature_genes_S2G <- c(
    "BRCA1",
    "BRCA2",
    "EZH2",
    "AURKA",
    "AURKB",
    "FOXM1",
    "UBE2C",
    "BIRC5"
  )
  
  signature_genes_6K <- c(
    "NMD3",
    "PEG10",
    "ENO2"
  )
  
  su2c_expression_genes <- unique(
    c(
      "XPO1",
      "NMD3",
      signature_genes_S2G,
      signature_genes_6K
    )
  )
  
  required_features_SU2C <- c(
    "OS_STATUS",
    "OS_MONTHS",
    su2c_expression_genes
  )
  
  # ------------------------------------------------------------
  # 3. Read expression and clinical Excel file
  # ------------------------------------------------------------
  
  su2c_raw <- readxl::read_excel(
    su2c_expression_file,
    .name_repair = "minimal"
  )
  
  stopifnot(
    "Hugo_Symbol" %in%
      colnames(su2c_raw),
    nrow(su2c_raw) > 1,
    ncol(su2c_raw) > 1
  )
  
  su2c_cleaned <- su2c_raw[
    -1,
    -1
  ]
  
  su2c_feature_names <- as.character(
    su2c_raw$Hugo_Symbol[-1]
  )
  
  required_feature_counts_SU2C <- vapply(
    required_features_SU2C,
    function(x) {
      sum(
        !is.na(su2c_feature_names) &
          su2c_feature_names == x
      )
    },
    integer(1)
  )
  
  if (
    any(
      required_feature_counts_SU2C != 1
    )
  ) {
    stop(
      "Required feature counts are incorrect: ",
      paste(
        paste0(
          names(required_feature_counts_SU2C),
          "=",
          required_feature_counts_SU2C
        ),
        collapse = ", "
      )
    )
  }
  
  transposed_data <- as.data.frame(
    t(su2c_cleaned),
    stringsAsFactors = FALSE
  )
  
  colnames(transposed_data) <-
    su2c_feature_names
  
  # ------------------------------------------------------------
  # 4. Numeric and survival-status parsers
  # ------------------------------------------------------------
  
  to_numeric_SU2C <- function(x) {
    suppressWarnings(
      as.numeric(
        as.character(x)
      )
    )
  }
  
  parse_status_SU2C <- function(x) {
    
    x <- toupper(
      trimws(
        as.character(x)
      )
    )
    
    dplyr::case_when(
      grepl(
        "DECEASED|DEAD|^1$|^1:",
        x
      ) ~ 1,
      
      grepl(
        "LIVING|ALIVE|^0$|^0:",
        x
      ) ~ 0,
      
      TRUE ~ suppressWarnings(
        as.numeric(x)
      )
    )
  }
  
  # ------------------------------------------------------------
  # 5. Build the sample-level expression table
  # ------------------------------------------------------------
  
  su2c_expression_data <- as.data.frame(
    lapply(
      su2c_expression_genes,
      function(gene_name) {
        to_numeric_SU2C(
          transposed_data[[gene_name]]
        )
      }
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(su2c_expression_data) <-
    su2c_expression_genes
  
  su2c_sample_data <- data.frame(
    SAMPLE_ID = rownames(
      transposed_data
    ),
    
    OS_STATUS = parse_status_SU2C(
      transposed_data[["OS_STATUS"]]
    ),
    
    OS_MONTHS = to_numeric_SU2C(
      transposed_data[["OS_MONTHS"]]
    ),
    
    su2c_expression_data,
    
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  stopifnot(
    !anyDuplicated(
      su2c_sample_data$SAMPLE_ID
    )
  )
  
  # ------------------------------------------------------------
  # 6. Read official SAMPLE_ID-to-PATIENT_ID mapping
  # ------------------------------------------------------------
  
  su2c_clinical_sample <- read.delim(
    su2c_clinical_sample_file,
    sep = "\t",
    header = TRUE,
    comment.char = "#",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  stopifnot(
    all(
      c(
        "SAMPLE_ID",
        "PATIENT_ID"
      ) %in%
        colnames(su2c_clinical_sample)
    )
  )
  
  su2c_sample_patient_map <-
    su2c_clinical_sample %>%
    dplyr::transmute(
      SAMPLE_ID = trimws(
        as.character(SAMPLE_ID)
      ),
      PATIENT_ID = trimws(
        as.character(PATIENT_ID)
      )
    ) %>%
    dplyr::filter(
      !is.na(SAMPLE_ID),
      SAMPLE_ID != "",
      !is.na(PATIENT_ID),
      PATIENT_ID != ""
    ) %>%
    dplyr::distinct()
  
  mapping_conflicts_SU2C <-
    su2c_sample_patient_map %>%
    dplyr::count(
      SAMPLE_ID,
      name = "Number_of_patients"
    ) %>%
    dplyr::filter(
      Number_of_patients > 1
    )
  
  if (nrow(mapping_conflicts_SU2C) > 0) {
    print(mapping_conflicts_SU2C)
    
    stop(
      "At least one SU2C sample maps ",
      "to multiple patients."
    )
  }
  
  su2c_sample_data <-
    su2c_sample_data %>%
    dplyr::left_join(
      su2c_sample_patient_map,
      by = "SAMPLE_ID"
    )
  
  unmapped_samples_SU2C <-
    su2c_sample_data %>%
    dplyr::filter(
      is.na(PATIENT_ID)
    ) %>%
    dplyr::select(SAMPLE_ID)
  
  if (nrow(unmapped_samples_SU2C) > 0) {
    print(unmapped_samples_SU2C)
    
    stop(
      "At least one expression sample could ",
      "not be mapped to a patient."
    )
  }
  
  # ------------------------------------------------------------
  # 7. Select samples complete for all four panels
  # ------------------------------------------------------------
  
  su2c_complete_samples <-
    su2c_sample_data %>%
    dplyr::filter(
      is.finite(OS_MONTHS),
      OS_MONTHS > 0,
      OS_STATUS %in% c(0, 1),
      
      dplyr::if_all(
        dplyr::all_of(
          su2c_expression_genes
        ),
        is.finite
      )
    )
  
  # Check repeated patients for consistent survival records.
  
  survival_conflicts_SU2C <-
    su2c_complete_samples %>%
    dplyr::group_by(PATIENT_ID) %>%
    dplyr::summarise(
      Number_of_OS_times =
        dplyr::n_distinct(OS_MONTHS),
      
      Number_of_OS_statuses =
        dplyr::n_distinct(OS_STATUS),
      
      .groups = "drop"
    ) %>%
    dplyr::filter(
      Number_of_OS_times > 1 |
        Number_of_OS_statuses > 1
    )
  
  if (nrow(survival_conflicts_SU2C) > 0) {
    print(survival_conflicts_SU2C)
    
    stop(
      "Repeated samples have conflicting ",
      "survival information."
    )
  }
  
  # ------------------------------------------------------------
  # 8. Collapse repeated samples to one patient
  # ------------------------------------------------------------
  
  su2c_patient_data <-
    su2c_complete_samples %>%
    dplyr::group_by(PATIENT_ID) %>%
    dplyr::summarise(
      SAMPLE_ID = paste(
        sort(
          unique(SAMPLE_ID)
        ),
        collapse = ";"
      ),
      
      Number_of_samples =
        dplyr::n_distinct(SAMPLE_ID),
      
      OS_MONTHS =
        dplyr::first(OS_MONTHS),
      
      OS_STATUS =
        dplyr::first(OS_STATUS),
      
      dplyr::across(
        dplyr::all_of(
          su2c_expression_genes
        ),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
    )
  
  cat(
    "Eligible SU2C samples:",
    nrow(su2c_complete_samples),
    "\n"
  )
  
  cat(
    "Unique SU2C patients:",
    nrow(su2c_patient_data),
    "\n"
  )
  
  cat(
    "Observed events:",
    sum(
      su2c_patient_data$OS_STATUS == 1
    ),
    "\n"
  )
  
  cat("Patients contributing multiple samples:\n")
  
  print(
    su2c_patient_data %>%
      dplyr::filter(
        Number_of_samples > 1
      ) %>%
      dplyr::select(
        PATIENT_ID,
        SAMPLE_ID,
        Number_of_samples
      )
  )
  
  stopifnot(
    nrow(su2c_complete_samples) == 81,
    nrow(su2c_patient_data) == 80,
    sum(
      su2c_patient_data$OS_STATUS == 1
    ) == 46
  )
  
  utils::write.csv(
    su2c_patient_data,
    file = paste0(
      "SourceData_SU2C_Common_",
      "Patient_Level_Data.csv"
    ),
    row.names = FALSE
  )
  
  # Ordinary log-rank P-value function.
  
  logrank_p_SU2C <- function(x) {
    stats::pchisq(
      x$chisq,
      df = length(x$n) - 1,
      lower.tail = FALSE
    )
  }
  
}))

# =============================================================================
# FigS2G
# =============================================================================
nmd3_step("FigS2G", c("Prepare_SU2C"), quote({
  # FigS2G — analysis below
  # NMD3 expression and overall survival
  #
  # Algorithm preserved:
  # outcome-optimized cutoff + ordinary log-rank P
  # ============================================================
  
  stopifnot(
    exists("su2c_patient_data"),
    exists("logrank_p_SU2C")
  )
  
  survival_data_S4F <- su2c_patient_data %>%
    dplyr::select(
      PATIENT_ID,
      SAMPLE_ID,
      Number_of_samples,
      OS_MONTHS,
      OS_STATUS,
      NMD3
    )
  
  cutpoint_S4F <- survminer::surv_cutpoint(
    survival_data_S4F,
    time = "OS_MONTHS",
    event = "OS_STATUS",
    variables = "NMD3",
    minprop = 0.10
  )
  
  cutoff_value_S4F <- as.numeric(
    cutpoint_S4F$cutpoint[
      "NMD3",
      "cutpoint"
    ]
  )
  
  survival_data_S4F$NMD3_Group <-
    factor(
      ifelse(
        survival_data_S4F$NMD3 >
          cutoff_value_S4F,
        "NMD3-high",
        "NMD3-low"
      ),
      levels = c(
        "NMD3-high",
        "NMD3-low"
      )
    )
  
  fit_S4F <- survival::survfit(
    survival::Surv(
      OS_MONTHS,
      OS_STATUS
    ) ~ NMD3_Group,
    data = survival_data_S4F
  )
  
  logrank_S4F <- survival::survdiff(
    survival::Surv(
      OS_MONTHS,
      OS_STATUS
    ) ~ NMD3_Group,
    data = survival_data_S4F
  )
  
  p_value_S4F <- logrank_p_SU2C(
    logrank_S4F
  )
  
  group_counts_S4F <- table(
    survival_data_S4F$NMD3_Group
  )
  
  cat(
    "\nFigure S4F (SU2C) NMD3 cutoff:",
    cutoff_value_S4F,
    "\n"
  )
  
  cat(
    "Figure S4F (SU2C) ordinary log-rank P:",
    p_value_S4F,
    "\n"
  )
  
  p_label_S4F <- if (
    p_value_S4F < 0.001
  ) {
    "p < 0.001"
  } else {
    paste0(
      "p = ",
      sprintf(
        "%.3f",
        p_value_S4F
      )
    )
  }
  
  legend_labels_S4F <- c(
    paste0(
      "NMD3-high n = ",
      group_counts_S4F["NMD3-high"]
    ),
    paste0(
      "NMD3-low n = ",
      group_counts_S4F["NMD3-low"]
    )
  )
  
  km_S4F <- survminer::ggsurvplot(
    fit_S4F,
    data = survival_data_S4F,
    
    pval = FALSE,
    conf.int = TRUE,
    risk.table = FALSE,
    censor = TRUE,
    
    xlim = c(
      0,
      60
    ),
    break.time.by = 10,
    
    legend.title = "",
    legend.labs = legend_labels_S4F,
    
    palette = c(
      "#F04444",
      "#1F78B4"
    ),
    
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    
    ggtheme =
      ggplot2::theme_minimal(
        base_size = 15
      )
  )
  
  p_S4F <- km_S4F$plot +
    ggplot2::annotate(
      "text",
      x = 2,
      y = 0.12,
      label = p_label_S4F,
      hjust = 0,
      size = 5,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = 2,
      y = 0.04,
      label = "SU2C/PCF",
      hjust = 0,
      size = 5,
      fontface = "bold"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        0,
        60
      ),
      ylim = c(
        0,
        1.03
      ),
      clip = "off"
    ) +
    ggplot2::theme(
      panel.border =
        ggplot2::element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.8
        ),
      
      axis.text =
        ggplot2::element_text(
          color = "black"
        ),
      
      legend.position =
        c(
          0.72,
          0.90
        ),
      
      legend.background =
        ggplot2::element_blank()
    )
  
  
  source_data_S4F <- survival_data_S4F %>%
    dplyr::mutate(
      NMD3_Cutoff =
        cutoff_value_S4F,
      
      Cutoff_Method =
        "Outcome-optimized surv_cutpoint, minprop=0.10",
      
      P_Value_Method =
        "Ordinary log-rank",
      
      Logrank_P =
        p_value_S4F
    )
  
  utils::write.csv(
    source_data_S4F,
    file = paste0(
      "SourceData_FigureS4F_SU2C_NMD3_",
      "patient_level.csv"
    ),
    row.names = FALSE
  )
  
  nmd3_save_pdf(
    filename = paste0(
      "FigS2G_SU2C_NMD3_Survival.pdf"
    ),
    plot = p_S4F,
    width = 6,
    height = 5.5,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  # ------------------------------------------------------------
  # S4F summary statistics
  # ------------------------------------------------------------
  
  statistics_S4F <- data.frame(
    Cohort = "SU2C/PCF 2019",
    Cutoff_Method = "Outcome-optimized surv_cutpoint, minprop=0.10",
    NMD3_Cutoff = cutoff_value_S4F,
    Total_Patients = nrow(survival_data_S4F),
    NMD3_High_n = unname(
      group_counts_S4F["NMD3-high"]
    ),
    NMD3_Low_n = unname(
      group_counts_S4F["NMD3-low"]
    ),
    Events = sum(
      survival_data_S4F$OS_STATUS == 1
    ),
    Logrank_P = p_value_S4F,
    stringsAsFactors = FALSE
  )
  
  print(statistics_S4F)
  
  utils::write.csv(
    statistics_S4F,
    "SourceData_FigureS4F_SU2C_NMD3_statistics.csv",
    row.names = FALSE
  )
  
}))

# =============================================================================
# FigS3B
# =============================================================================
nmd3_step("FigS3B", character(), quote({
  # FigS3B — analysis below
  # NMD3 nuclear-to-cytoplasmic protein ratio
  # DMSO vs KPT-330
  #
  # Each point represents one biological replicate.
  # Reuses prot_raw generated in Figure 2.
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
  })
  
  stopifnot(
    exists("prot_raw"),
    "Gene" %in% colnames(prot_raw)
  )
  
  # ------------------------------------------------------------
  # 1. Select normalized DMSO and KPT-330 protein columns
  # ------------------------------------------------------------
  
  sample_columns_S4A <- grep(
    "^Norm_.*(DMSO|XPO1inhibitor|KPT).*(Cyto|Nuc)",
    colnames(prot_raw),
    value = TRUE,
    ignore.case = TRUE
  )
  
  stopifnot(length(sample_columns_S4A) > 0)
  
  cat("Protein samples considered for Figure S4A:\n")
  print(sample_columns_S4A)
  
  # ------------------------------------------------------------
  # 2. Extract NMD3 protein abundance
  # ------------------------------------------------------------
  
  nmd3_long_S4A <- prot_raw %>%
    dplyr::filter(
      !is.na(Gene),
      toupper(as.character(Gene)) == "NMD3"
    ) %>%
    dplyr::select(
      Gene,
      dplyr::all_of(sample_columns_S4A)
    ) %>%
    tidyr::pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Expression"
    ) %>%
    dplyr::mutate(
      Expression = suppressWarnings(
        as.numeric(Expression)
      ),
      Treatment = dplyr::case_when(
        stringr::str_detect(
          Sample,
          stringr::regex(
            "DMSO",
            ignore_case = TRUE
          )
        ) ~ "Ctrl",
        
        stringr::str_detect(
          Sample,
          stringr::regex(
            "XPO1inhibitor|KPT",
            ignore_case = TRUE
          )
        ) ~ "KPT-330",
        
        TRUE ~ NA_character_
      ),
      Compartment = dplyr::case_when(
        stringr::str_detect(
          Sample,
          stringr::regex(
            "Nuc",
            ignore_case = TRUE
          )
        ) ~ "Nucleus",
        
        stringr::str_detect(
          Sample,
          stringr::regex(
            "Cyto|Cyt",
            ignore_case = TRUE
          )
        ) ~ "Cytoplasm",
        
        TRUE ~ NA_character_
      ),
      Replicate = stringr::str_match(
        Sample,
        "_[DX]([123])[CN]$"
      )[, 2]
    ) %>%
    dplyr::filter(
      !is.na(Treatment),
      !is.na(Compartment),
      !is.na(Replicate),
      is.finite(Expression)
    ) %>%
    dplyr::group_by(
      Treatment,
      Replicate,
      Compartment
    ) %>%
    dplyr::summarise(
      Expression = mean(
        Expression,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  cat("\nNMD3 normalized protein abundance:\n")
  print(nmd3_long_S4A)
  
  # ------------------------------------------------------------
  # 3. Calculate NMD3 log2(Nucleus/Cytoplasm)
  # ------------------------------------------------------------
  
  nmd3_ratio_S4A <- nmd3_long_S4A %>%
    tidyr::pivot_wider(
      names_from = Compartment,
      values_from = Expression
    )
  
  stopifnot(
    all(
      c("Nucleus", "Cytoplasm") %in%
        colnames(nmd3_ratio_S4A)
    )
  )
  
  nmd3_ratio_S4A <- nmd3_ratio_S4A %>%
    dplyr::filter(
      is.finite(Nucleus),
      is.finite(Cytoplasm),
      Nucleus >= 0,
      Cytoplasm >= 0
    ) %>%
    dplyr::mutate(
      log2_NC = log2(
        (Nucleus + 1e-6) /
          (Cytoplasm + 1e-6)
      ),
      Treatment = factor(
        Treatment,
        levels = c(
          "Ctrl",
          "KPT-330"
        )
      )
    ) %>%
    dplyr::arrange(
      Treatment,
      Replicate
    )
  
  cat("\nNMD3 nuclear/cytoplasmic ratios:\n")
  print(nmd3_ratio_S4A)
  
  cat("\nSamples per treatment:\n")
  print(table(nmd3_ratio_S4A$Treatment))
  
  test_S4A <- stats::t.test(
    log2_NC ~ Treatment,
    data = nmd3_ratio_S4A,
    alternative = "two.sided"
  )
  
  print(test_S4A)
  
  pvalue_S4A <- test_S4A$p.value
  
  pvalue_label_S4A <- paste0(
    "p = ",
    format.pval(
      pvalue_S4A,
      digits = 2,
      eps = 0.001
    )
  )
  
  cat(
    "Welch t-test P value:",
    pvalue_S4A,
    "\n"
  )
  # ------------------------------------------------------------
  # 5. Figure S4A
  # ------------------------------------------------------------
  
  colors_S4A <- c(
    "Ctrl" = "#E68600",
    "KPT-330" = "#00A6E8"
  )
  
  y_range_S4A <- range(
    nmd3_ratio_S4A$log2_NC,
    finite = TRUE
  )
  
  p_S4A <- ggplot2::ggplot(
    nmd3_ratio_S4A,
    ggplot2::aes(
      x = Treatment,
      y = log2_NC,
      color = Treatment
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.45,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.8
    ) +
    ggplot2::geom_jitter(
      width = 0.06,
      size = 2.7,
      alpha = 1
    ) +
    ggplot2::annotate(
      geom = "text",
      x = 2.25,
      y = y_range_S4A[1],
      label = pvalue_label_S4A,
      hjust = 1,
      vjust = -0.5,
      size = 4.2,
      fontface = "bold",
      color = "black"
    ) +
    ggplot2::scale_color_manual(
      values = colors_S4A
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.12,
          0.10
        )
      )
    ) +
    ggplot2::labs(
      title = NULL,
      x = NULL,
      y = expression(
        Log[2]~"NMD3 (Nuclear / Cytoplasm)"
      )
    ) +
    ggplot2::theme_classic(
      base_size = 13
    ) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      axis.ticks = ggplot2::element_line(
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 11,
        color = "black"
      ),
      axis.title.y = ggplot2::element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(
        10,
        10,
        10,
        10
      )
    )
  
  
  # ------------------------------------------------------------
  # 6. Save Figure S4A
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS3B_NMD3_Nuclear_Cytoplasmic_Protein_Ratio.pdf",
    plot = p_S4A,
    width = 4,
    height = 5,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  utils::write.csv(
    nmd3_ratio_S4A,
    file = "FigureS4A_NMD3_Nuclear_Cytoplasmic_Protein_Ratio.csv",
    row.names = FALSE
  )
  
}))

# =============================================================================
# Prepare_TCGA_Paired
# =============================================================================
nmd3_step("Prepare_TCGA_Paired", c("Prepare_TCGA"), quote({
  # Prepare_TCGA_Paired — analysis below
  # =============================================================================
  
  all_samples <- sample_info %>%
    filter(
      file_name %in% basename(tsv_files),
      sample_type %in% c(
        "Primary Tumor",
        "Solid Tissue Normal"
      )
    ) %>%
    mutate(
      file = tsv_files[
        match(file_name, basename(tsv_files))
      ],
      patient = substr(sample_id, 1, 12)
    ) %>%
    group_by(patient) %>%
    filter(n_distinct(sample_type) == 2) %>%
    arrange(patient, sample_type, sample_id) %>%
    group_by(patient, sample_type) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      sample_id,
      file,
      condition = if_else(
        sample_type == "Primary Tumor",
        "Tumor",
        "Normal"
      ),
      patient
    )
  
  stopifnot(
    nrow(all_samples) == 104,
    n_distinct(all_samples$patient) == 52,
    all(table(all_samples$patient) == 2),
    !anyNA(all_samples$file)
  )
  
  counts_mat <- build_count_matrix(
    sample_files = all_samples$file,
    sample_ids = all_samples$sample_id
  )
  
  counts_mat <- map_counts_to_symbols(
    counts_mat
  )
  
  col_data <- all_samples %>%
    dplyr::select(
      sample_id,
      condition,
      patient
    ) %>%
    mutate(
      condition = factor(
        condition,
        levels = c("Normal", "Tumor")
      ),
      patient = factor(patient)
    ) %>%
    column_to_rownames("sample_id")
  
  counts_mat <- counts_mat[
    ,
    rownames(col_data),
    drop = FALSE
  ]
  
  keep_genes <- rowSums(
    counts_mat >= 10
  ) >= ceiling(
    0.1 * ncol(counts_mat)
  )
  
  counts_mat <- counts_mat[
    keep_genes,
    ,
    drop = FALSE
  ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = col_data,
    design = ~ patient + condition
  )
  
  dds <- DESeq(
    dds,
    quiet = TRUE
  )
  
  res <- results(
    dds,
    contrast = c(
      "condition",
      "Tumor",
      "Normal"
    )
  )
  
  gene_list <- res$stat
  names(gene_list) <- toupper(rownames(res))
  
  gene_list <- sort(
    gene_list[is.finite(gene_list)],
    decreasing = TRUE
  )
  
  stopifnot(
    length(gene_list) > 0,
    !anyDuplicated(names(gene_list))
  )
  
  # -----------------------------------------------------------------------------
  
}))

# =============================================================================
# FigS4A
# =============================================================================
nmd3_step("FigS4A", c("Prepare_TCGA_Paired"), quote({
  # FigS4A — analysis below
  # =============================================================================
  
  signature_genes <- c(
    "XPO1", "RAN", "RANBP2",
    "NUP214", "NUP88", "NUP62", "NUP153", "NUP50",
    "NUP98", "NUP133", "NUP35", "NUP155", "NUP160", "NUP205"
  )
  
  anchor_set <- c(
    "XPO1",
    "NUP214"
  )
  
  # Use DESeq2 size-factor-normalized counts from the 52 paired patients
  normalized_counts <- counts(
    dds,
    normalized = TRUE
  )
  
  stopifnot(
    all(
      signature_genes %in%
        rownames(normalized_counts)
    )
  )
  
  expr_mat <- normalized_counts[
    signature_genes,
    ,
    drop = FALSE
  ]
  
  tumor_sids <- rownames(col_data)[
    col_data$condition == "Tumor"
  ]
  
  normal_sids <- rownames(col_data)[
    col_data$condition == "Normal"
  ]
  
  mat_tumor <- expr_mat[
    ,
    tumor_sids,
    drop = FALSE
  ]
  
  mat_normal <- expr_mat[
    ,
    normal_sids,
    drop = FALSE
  ]
  
  stopifnot(
    ncol(mat_tumor) == 52,
    ncol(mat_normal) == 52
  )
  
  get_gene_corr <- function(mat_expression) {
    suppressWarnings(
      cor(
        t(mat_expression),
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )
  }
  
  corr_tumor <- get_gene_corr(
    mat_tumor
  )
  
  corr_normal <- get_gene_corr(
    mat_normal
  )
  
  # Common gene order for both groups
  corr_avg <- (
    corr_tumor + corr_normal
  ) / 2
  
  hc <- hclust(
    as.dist(1 - corr_avg),
    method = "average"
  )
  
  gene_order <- hc$labels[
    hc$order
  ]
  
  corr_to_long <- function(corr_mat, group_label) {
    corr_df <- as.data.frame(
      as.table(corr_mat)
    )
    
    colnames(corr_df) <- c(
      "GeneX",
      "GeneY",
      "Correlation"
    )
    
    corr_df$Group <- group_label
    corr_df
  }
  
  df_long <- bind_rows(
    corr_to_long(
      corr_tumor,
      "Tumor"
    ),
    corr_to_long(
      corr_normal,
      "Normal"
    )
  ) %>%
    mutate(
      GeneX = factor(
        GeneX,
        levels = gene_order,
        ordered = TRUE
      ),
      GeneY = factor(
        GeneY,
        levels = gene_order,
        ordered = TRUE
      ),
      Group = factor(
        Group,
        levels = c("Tumor", "Normal")
      )
    )
  
  df_plot <- df_long %>%
    filter(
      as.integer(GeneX) <
        as.integer(GeneY)
    ) %>%
    mutate(
      mark =
        as.character(GeneX) %in% anchor_set |
        as.character(GeneY) %in% anchor_set
    )
  
  # -----------------------------------------------------------------------------
  # Quantile-enhanced shared color scale
  # -----------------------------------------------------------------------------
  
  corr_vals <- df_plot$Correlation
  
  # Always include zero in the shared color range
  corr_min <- min(
    c(0, corr_vals),
    na.rm = TRUE
  )
  
  corr_max <- max(
    c(0, corr_vals),
    na.rm = TRUE
  )
  
  pad <- (corr_max - corr_min) * 0.05
  
  corr_min <- corr_min - pad
  corr_max <- corr_max + pad
  
  q_breaks <- quantile(
    corr_vals,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )
  
  color_stops <- sort(
    unique(
      c(
        corr_min,
        q_breaks,
        0,
        corr_max
      )
    )
  )
  
  n_below <- sum(color_stops < 0)
  n_above <- sum(color_stops > 0)
  
  pal_below <- if (n_below > 0) {
    colorRampPalette(
      c("#08306B", "#4575B4", "white")
    )(n_below + 1)[seq_len(n_below)]
  } else {
    character(0)
  }
  
  pal_above <- if (n_above > 0) {
    colorRampPalette(
      c("white", "#D73027", "#7F0000")
    )(n_above + 1)[-1]
  } else {
    character(0)
  }
  
  color_map <- rep(
    NA_character_,
    length(color_stops)
  )
  
  color_map[color_stops < 0] <- pal_below
  color_map[color_stops == 0] <- "white"
  color_map[color_stops > 0] <- pal_above
  
  # -----------------------------------------------------------------------------
  # Plot
  # -----------------------------------------------------------------------------
  
  fig1c <- ggplot(
    df_plot,
    aes(
      x = GeneX,
      y = GeneY,
      fill = Correlation
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.2
    ) +
    geom_tile(
      data = filter(df_plot, mark),
      fill = NA,
      color = "#5CB85C",
      linewidth = 1.1,
      alpha = 0.65
    ) +
    coord_fixed() +
    facet_wrap(
      ~ Group,
      nrow = 1
    ) +
    scale_fill_gradientn(
      colours = color_map,
      values = scales::rescale(
        color_stops,
        from = c(corr_min, corr_max)
      ),
      limits = c(
        corr_min,
        corr_max
      ),
      breaks = scales::pretty_breaks(n = 6)(
        c(corr_min, corr_max)
      ),
      oob = scales::squish,
      name = "Spearman r"
    ) +
    labs(
      title = "Nuclear transport signature correlation",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        color = "black"
      ),
      axis.text.y = element_text(
        color = "black"
      ),
      strip.text = element_text(
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      )
    )
  
  
  nmd3_save_pdf(
    "FigS4A_transport_correlation.pdf",
    plot = fig1c,
    width = 9,
    height = 5,
    units = "in"
  )
  
}))

# =============================================================================
# FigS4B
# =============================================================================
nmd3_step("FigS4B", c("Prepare_TCGA_Paired"), quote({
  # FigS4B — analysis below
  # -----------------------------------------------------------------------------
  
  translation_set <- msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:REACTOME"
  ) %>%
    filter(
      gs_name == "REACTOME_TRANSLATION"
    ) %>%
    transmute(
      term = "REACTOME_TRANSLATION",
      gene = toupper(gene_symbol)
    ) %>%
    distinct()
  
  stopifnot(nrow(translation_set) > 0)
  
  gsea_translation <- run_gsea(
    gene_ranking = gene_list,
    term_to_gene = translation_set
  )
  
  fig1a <- gseaplot2(
    gsea_translation,
    geneSetID = 1,
    subplots = 1:3,
    pvalue_table = FALSE
  )
  
  nmd3_save_pdf(
    "FigS4B_TRANSLATION_GSEA.pdf",
    plot = fig1a,
    width = 8,
    height = 6,
    units = "in"
  )
  
  print(
    gsea_translation@result[
      ,
      c(
        "ID",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
}))

# =============================================================================
# FigS4C
# =============================================================================
nmd3_step("FigS4C", c("Prepare_TCGA_Paired"), quote({
  # FigS4C — analysis below
  # -----------------------------------------------------------------------------
  
  nuclear_transport_set <- msigdbr(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  ) %>%
    filter(
      gs_exact_source == "GO:0051169"
    ) %>%
    transmute(
      term = "GOBP_NUCLEAR_TRANSPORT",
      gene = toupper(gene_symbol)
    ) %>%
    distinct()
  
  stopifnot(nrow(nuclear_transport_set) > 0)
  
  gsea_nuclear_transport <- run_gsea(
    gene_ranking = gene_list,
    term_to_gene = nuclear_transport_set
  )
  
  fig1b <- gseaplot2(
    gsea_nuclear_transport,
    geneSetID = 1,
    subplots = 1:3,
    pvalue_table = FALSE
  )
  
  nmd3_save_pdf(
    "FigS4C_NUCLEAR_TRANSPORT_GSEA.pdf",
    plot = fig1b,
    width = 8,
    height = 6,
    units = "in"
  )
  
  print(
    gsea_nuclear_transport@result[
      ,
      c(
        "ID",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
  
}))

# =============================================================================
# FigS4D
# =============================================================================
nmd3_step("FigS4D", c("Prepare_XPO1_Groups","FigS4C"), quote({
  dds_xpo1 <- DESeqDataSetFromMatrix(
    countData = primary_counts,
    colData = primary_coldata,
    design = ~ group
  )
  
  dds_xpo1 <- DESeq(
    dds_xpo1,
    quiet = TRUE
  )
  
  res_xpo1 <- results(
    dds_xpo1,
    contrast = c(
      "group",
      "XPO1_high",
      "XPO1_low"
    )
  )
  
  gene_ranking_xpo1 <- res_xpo1$stat
  
  names(gene_ranking_xpo1) <- toupper(
    rownames(res_xpo1)
  )
  
  gene_ranking_xpo1 <- sort(
    gene_ranking_xpo1[
      is.finite(gene_ranking_xpo1)
    ],
    decreasing = TRUE
  )
  
  stopifnot(
    length(gene_ranking_xpo1) > 0,
    !anyDuplicated(names(gene_ranking_xpo1))
  )
  
  gsea_xpo1_nuclear_transport <- run_gsea(
    gene_ranking = gene_ranking_xpo1,
    term_to_gene = nuclear_transport_set
  )
  
  fig1d <- gseaplot2(
    gsea_xpo1_nuclear_transport,
    geneSetID = 1,
    subplots = 1:3,
    pvalue_table = FALSE
  )
  
  
  print(
    gsea_xpo1_nuclear_transport@result[
      ,
      c(
        "ID",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
  nmd3_save_pdf(
    "FigS4D_XPO1_Nuclear_Transport_GSEA.pdf",
    plot = fig1d,
    width = 8,
    height = 6,
    units = "in"
  )
  
  
  
}))

# =============================================================================
# FigS4F
# =============================================================================
nmd3_step("FigS4F", character(), quote({
  # FigS4F — analysis below
  # PCA of DMSO and KPT-330 cytoplasmic/nuclear RNA-seq samples
  #
  # Reuses expr_matrix generated in Figure 2
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggforce)
  })
  
  stopifnot(
    exists("expr_matrix"),
    "Length" %in% colnames(expr_matrix)
  )
  
  # ------------------------------------------------------------
  # 1. Convert raw counts to TPM
  # ------------------------------------------------------------
  
  gene_length_kb_S2A <- as.numeric(
    expr_matrix$Length
  ) / 1000
  
  count_matrix_S2A <- as.matrix(
    expr_matrix[
      ,
      setdiff(colnames(expr_matrix), "Length"),
      drop = FALSE
    ]
  )
  
  storage.mode(count_matrix_S2A) <- "numeric"
  
  valid_genes_S2A <- (
    is.finite(gene_length_kb_S2A) &
      gene_length_kb_S2A > 0
  )
  
  gene_length_kb_S2A <- gene_length_kb_S2A[
    valid_genes_S2A
  ]
  
  count_matrix_S2A <- count_matrix_S2A[
    valid_genes_S2A,
    ,
    drop = FALSE
  ]
  
  rpk_S2A <- sweep(
    count_matrix_S2A,
    1,
    gene_length_kb_S2A,
    FUN = "/"
  )
  
  scaling_factors_S2A <- colSums(
    rpk_S2A,
    na.rm = TRUE
  )
  
  stopifnot(all(scaling_factors_S2A > 0))
  
  tpm_S2A <- sweep(
    rpk_S2A,
    2,
    scaling_factors_S2A,
    FUN = "/"
  ) * 1e6
  
  log2_tpm_S2A <- log2(
    tpm_S2A + 1
  )
  
  # ------------------------------------------------------------
  # 2. Select DMSO and KPT-330 samples
  # ------------------------------------------------------------
  
  samples_S2A <- grep(
    "^[DX][0-9]+[CN]$",
    colnames(log2_tpm_S2A),
    value = TRUE
  )
  
  stopifnot(length(samples_S2A) > 0)
  
  log2_tpm_S2A <- log2_tpm_S2A[
    ,
    samples_S2A,
    drop = FALSE
  ]
  
  cat("Samples used in Figure S2A:\n")
  print(samples_S2A)
  
  # ------------------------------------------------------------
  # 3. Remove genes with zero variance
  # ------------------------------------------------------------
  
  gene_sd_S2A <- apply(
    log2_tpm_S2A,
    1,
    stats::sd,
    na.rm = TRUE
  )
  
  keep_genes_S2A <- (
    is.finite(gene_sd_S2A) &
      gene_sd_S2A > 1e-6
  )
  
  log2_tpm_S2A <- log2_tpm_S2A[
    keep_genes_S2A,
    ,
    drop = FALSE
  ]
  
  cat(
    "Genes used for PCA:",
    nrow(log2_tpm_S2A),
    "\n"
  )
  
  cat(
    "Samples used for PCA:",
    ncol(log2_tpm_S2A),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 4. PCA
  # ------------------------------------------------------------
  
  pca_result_S2A <- stats::prcomp(
    t(log2_tpm_S2A),
    center = TRUE,
    scale. = TRUE
  )
  
  variance_S2A <- summary(
    pca_result_S2A
  )$importance[2, 1:2] * 100
  
  # ------------------------------------------------------------
  # 5. Sample annotation
  # ------------------------------------------------------------
  
  pca_df_S2A <- as.data.frame(
    pca_result_S2A$x[, 1:2, drop = FALSE]
  ) %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(
      Prefix = substr(Sample, 1, 1),
      Suffix = substr(Sample, nchar(Sample), nchar(Sample)),
      Compartment = dplyr::if_else(
        Suffix == "C",
        "Cytoplasm",
        "Nucleus"
      ),
      Treatment = dplyr::case_when(
        Prefix == "D" ~ "DMSO",
        Prefix == "X" ~ "KPT-330",
        TRUE ~ NA_character_
      ),
      GroupLabel = paste0(
        Treatment,
        "-",
        Compartment
      ),
      GroupLabel = factor(
        GroupLabel,
        levels = c(
          "DMSO-Cytoplasm",
          "DMSO-Nucleus",
          "KPT-330-Cytoplasm",
          "KPT-330-Nucleus"
        )
      )
    )
  
  cat("Samples per group:\n")
  print(table(pca_df_S2A$GroupLabel))
  
  # ------------------------------------------------------------
  # 6. Colors and shapes
  # ------------------------------------------------------------
  
  colors_S2A <- c(
    "DMSO-Cytoplasm"    = "#E31A1C",
    "DMSO-Nucleus"      = "#FF7F00",
    "KPT-330-Cytoplasm" = "#1F78B4",
    "KPT-330-Nucleus"   = "#33A02C"
  )
  
  shapes_S2A <- c(
    "DMSO-Cytoplasm"    = 15,
    "DMSO-Nucleus"      = 3,
    "KPT-330-Cytoplasm" = 18,
    "KPT-330-Nucleus"   = 7
  )
  
  # ------------------------------------------------------------
  # 7. Figure S2A
  # ------------------------------------------------------------
  
  p_S2A <- ggplot2::ggplot(
    pca_df_S2A,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = GroupLabel,
      shape = GroupLabel
    )
  ) +
    ggforce::geom_mark_ellipse(
      ggplot2::aes(
        group = GroupLabel,
        fill = GroupLabel
      ),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 4,
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = colors_S2A
    ) +
    ggplot2::scale_fill_manual(
      values = colors_S2A
    ) +
    ggplot2::scale_shape_manual(
      values = shapes_S2A
    ) +
    ggplot2::coord_fixed(
      ratio = 1.5
    ) +
    ggplot2::labs(
      title = NULL,
      x = paste0(
        "PC1 (",
        round(variance_S2A[1], 1),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(variance_S2A[2], 1),
        "%)"
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        color = "black",
        face = "bold"
      ),
      axis.line = ggplot2::element_line(
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        color = "black"
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = "right",
      panel.grid = ggplot2::element_blank()
    )
  
  
  # ------------------------------------------------------------
  # 8. Save Figure S2A
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS4F_KPT330_RNAseq_PCA.pdf",
    plot = p_S2A,
    width = 10,
    height = 7,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  
  
}))

# =============================================================================
# FigS4G
# =============================================================================
nmd3_step("FigS4G", character(), quote({
  # FigS4G — analysis below
  # PCA of DMSO and KPT-330 cytoplasmic/nuclear proteomics
  #
  # Reuses prot_raw generated in Figure 2
  #
  # CHANGE FROM ORIGINAL:
  #   Proteins with any missing value across the selected samples
  #   are now excluded before PCA, instead of being imputed as 0.
  #   Filling NA with 0 + scale.=TRUE can create artificial outliers
  #   that dominate PCA loadings, especially when missingness is
  #   non-random (e.g. low-abundance proteins, or missingness that
  #   differs systematically by compartment/treatment).
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggforce)
  })
  
  stopifnot(
    exists("prot_raw"),
    "Gene" %in% colnames(prot_raw)
  )
  
  # ------------------------------------------------------------
  # 1. Select DMSO and KPT-330 normalized protein columns
  # ------------------------------------------------------------
  
  protein_data_S3A <- as.data.frame(
    prot_raw,
    check.names = FALSE
  )
  
  protein_data_S3A <- protein_data_S3A[
    !is.na(protein_data_S3A$Gene) &
      protein_data_S3A$Gene != "",
    ,
    drop = FALSE
  ]
  
  normalized_columns_S3A <- grep(
    "^Norm_",
    colnames(protein_data_S3A),
    value = TRUE
  )
  
  selected_samples_S3A <- normalized_columns_S3A[
    grepl(
      "DMSO|XPO1inhibitor|KPT",
      normalized_columns_S3A,
      ignore.case = TRUE
    )
  ]
  
  stopifnot(length(selected_samples_S3A) > 0)
  
  cat("Proteomics samples used in Figure S3A:\n")
  print(selected_samples_S3A)
  
  # ------------------------------------------------------------
  # 2. Build the protein-expression matrix
  # ------------------------------------------------------------
  
  protein_matrix_df_S3A <- protein_data_S3A[
    ,
    selected_samples_S3A,
    drop = FALSE
  ]
  
  protein_matrix_df_S3A[] <- lapply(
    protein_matrix_df_S3A,
    function(x) {
      suppressWarnings(
        as.numeric(as.character(x))
      )
    }
  )
  
  protein_matrix_S3A <- as.matrix(
    protein_matrix_df_S3A
  )
  
  rownames(protein_matrix_S3A) <- make.unique(
    toupper(
      as.character(
        protein_data_S3A$Gene
      )
    )
  )
  
  storage.mode(protein_matrix_S3A) <- "numeric"
  
  # ------------------------------------------------------------
  # 3. Remove proteins with ANY missing value across selected
  #    samples (complete-case filtering, no imputation), then
  #    remove proteins/samples with zero variance
  # ------------------------------------------------------------
  
  n_proteins_before_S3A <- nrow(protein_matrix_S3A)
  
  missing_per_protein_S3A <- rowSums(
    is.na(protein_matrix_S3A) | !is.finite(protein_matrix_S3A)
  )
  
  complete_proteins_S3A <- missing_per_protein_S3A == 0
  
  protein_matrix_S3A <- protein_matrix_S3A[
    complete_proteins_S3A,
    ,
    drop = FALSE
  ]
  
  cat(
    "Proteins removed due to missing values in >=1 sample:",
    n_proteins_before_S3A - nrow(protein_matrix_S3A),
    "\n"
  )
  
  cat(
    "Proteins retained with complete data across all samples:",
    nrow(protein_matrix_S3A),
    "\n"
  )
  
  protein_sd_S3A <- apply(
    protein_matrix_S3A,
    1,
    stats::sd,
    na.rm = TRUE
  )
  
  keep_proteins_S3A <- (
    is.finite(protein_sd_S3A) &
      protein_sd_S3A > 0
  )
  
  protein_matrix_S3A <- protein_matrix_S3A[
    keep_proteins_S3A,
    ,
    drop = FALSE
  ]
  
  sample_sd_S3A <- apply(
    protein_matrix_S3A,
    2,
    stats::sd,
    na.rm = TRUE
  )
  
  keep_samples_S3A <- (
    is.finite(sample_sd_S3A) &
      sample_sd_S3A > 0
  )
  
  protein_matrix_S3A <- protein_matrix_S3A[
    ,
    keep_samples_S3A,
    drop = FALSE
  ]
  
  if (
    nrow(protein_matrix_S3A) < 2 ||
    ncol(protein_matrix_S3A) < 2
  ) {
    stop(
      "The proteomics matrix is too small for PCA."
    )
  }
  
  cat(
    "Proteins used for PCA:",
    nrow(protein_matrix_S3A),
    "\n"
  )
  
  cat(
    "Samples used for PCA:",
    ncol(protein_matrix_S3A),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 4. PCA
  # ------------------------------------------------------------
  
  pca_result_S3A <- stats::prcomp(
    t(protein_matrix_S3A),
    center = TRUE,
    scale. = TRUE
  )
  
  variance_explained_S3A <- (
    pca_result_S3A$sdev^2 /
      sum(pca_result_S3A$sdev^2)
  )
  
  pc1_label_S3A <- sprintf(
    "PC1 (%.1f%%)",
    100 * variance_explained_S3A[1]
  )
  
  pc2_label_S3A <- sprintf(
    "PC2 (%.1f%%)",
    100 * variance_explained_S3A[2]
  )
  
  # ------------------------------------------------------------
  # 5. Sample annotation
  # ------------------------------------------------------------
  
  pca_df_S3A <- as.data.frame(
    pca_result_S3A$x[
      ,
      1:2,
      drop = FALSE
    ]
  ) %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(
      Treatment = dplyr::case_when(
        grepl(
          "XPO1inhibitor|KPT",
          Sample,
          ignore.case = TRUE
        ) ~ "KPT-330",
        grepl(
          "DMSO",
          Sample,
          ignore.case = TRUE
        ) ~ "DMSO",
        TRUE ~ NA_character_
      ),
      Compartment = dplyr::case_when(
        grepl(
          "Cyto|Cyt",
          Sample,
          ignore.case = TRUE
        ) ~ "Cytoplasm",
        grepl(
          "Nuc|Nucleus",
          Sample,
          ignore.case = TRUE
        ) ~ "Nucleus",
        TRUE ~ NA_character_
      ),
      GroupLabel = paste0(
        Treatment,
        "-",
        Compartment
      ),
      GroupLabel = factor(
        GroupLabel,
        levels = c(
          "DMSO-Cytoplasm",
          "DMSO-Nucleus",
          "KPT-330-Cytoplasm",
          "KPT-330-Nucleus"
        )
      )
    )
  
  if (
    any(
      is.na(pca_df_S3A$Treatment) |
      is.na(pca_df_S3A$Compartment)
    )
  ) {
    stop(
      "At least one proteomics sample could not be assigned ",
      "to a treatment or compartment."
    )
  }
  
  cat("Samples per group:\n")
  
  print(
    table(
      pca_df_S3A$GroupLabel,
      useNA = "ifany"
    )
  )
  
  # ------------------------------------------------------------
  # 6. Colors and shapes
  # ------------------------------------------------------------
  
  colors_S3A <- c(
    "DMSO-Cytoplasm" = "#E31A1C",
    "DMSO-Nucleus" = "#FF7F00",
    "KPT-330-Cytoplasm" = "#1F78B4",
    "KPT-330-Nucleus" = "#33A02C"
  )
  
  shapes_S3A <- c(
    "DMSO-Cytoplasm" = 15,
    "DMSO-Nucleus" = 3,
    "KPT-330-Cytoplasm" = 18,
    "KPT-330-Nucleus" = 7
  )
  
  # ------------------------------------------------------------
  # 7. Figure S3A
  # ------------------------------------------------------------
  
  p_S3A <- ggplot2::ggplot(
    pca_df_S3A,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = GroupLabel,
      shape = GroupLabel
    )
  ) +
    ggforce::geom_mark_ellipse(
      ggplot2::aes(
        group = GroupLabel,
        fill = GroupLabel
      ),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 4,
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = colors_S3A
    ) +
    ggplot2::scale_fill_manual(
      values = colors_S3A
    ) +
    ggplot2::scale_shape_manual(
      values = shapes_S3A
    ) +
    ggplot2::coord_fixed(
      ratio = 1.6
    ) +
    ggplot2::labs(
      title = NULL,
      x = pc1_label_S3A,
      y = pc2_label_S3A
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        color = "black",
        face = "bold"
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = "right"
    )
  
  
  # ------------------------------------------------------------
  # 8. Save Figure S3A
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS4G_KPT330_Proteomics_PCA.pdf",
    plot = p_S3A,
    width = 10,
    height = 7,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  
}))

# FigS4H is already manually drawn. No computational reconstruction is run.

# =============================================================================
# Prepare_KPT_RNA_Protein
# =============================================================================
nmd3_step("Prepare_KPT_RNA_Protein", c("Fig2A"), quote({
  # 1. Cytoplasmic RNA and protein differential-expression results
  # -----------------------------------------------------------------------------
  
  res_rna_df <- rna_cyto_df %>%
    dplyr::transmute(
      gene = toupper(gene),
      log2FC_RNA = log2FC_RNA_cyto,
      padj_RNA = padj_RNA_cyto
    )
  
  expr_cols <- c(
    "Norm_DMSO_Cyto_D1C",
    "Norm_DMSO_Cyto_D2C",
    "Norm_DMSO_Cyto_D3C",
    "Norm_XPO1inhibitor_Cyto_X1C",
    "Norm_XPO1inhibitor_Cyto_X2C",
    "Norm_XPO1inhibitor_Cyto_X3C"
  )
  
  stopifnot(
    length(expr_cols) == 6L,
    !anyDuplicated(expr_cols),
    all(expr_cols %in% colnames(prot_raw))
  )
  
  prot_mat <- prot_raw %>%
    dplyr::filter(!is.na(Gene)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(expr_cols),
        ~ suppressWarnings(as.numeric(.x))
      )
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(expr_cols),
        ~ if (all(is.na(.x))) NA_real_ else median(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(expr_cols), ~ !is.na(.x))
    )
  
  prot_expr <- as.matrix(prot_mat[, expr_cols])
  rownames(prot_expr) <- toupper(prot_mat$Gene)
  storage.mode(prot_expr) <- "double"
  
  stopifnot(
    nrow(prot_expr) > 0L,
    ncol(prot_expr) == 6L,
    !anyDuplicated(rownames(prot_expr)),
    all(is.finite(prot_expr))
  )
  
  prot_group <- factor(
    c(rep("DMSO", 3L), rep("KPT", 3L)),
    levels = c("DMSO", "KPT")
  )
  prot_design <- model.matrix(~ prot_group)
  fit <- limma::lmFit(prot_expr, prot_design)
  contrast <- limma::makeContrasts(
    KPTvsDMSO = prot_groupKPT,
    levels = prot_design
  )
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2, robust = TRUE, trend = TRUE)
  
  prot_result <- limma::topTable(
    fit2,
    coef = "KPTvsDMSO",
    number = Inf,
    sort.by = "none"
  )
  
  res_prot_df <- prot_result %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::transmute(
      gene = toupper(gene),
      log2FC_Prot = as.numeric(logFC),
      pval_Prot = as.numeric(P.Value),
      padj_Prot = as.numeric(adj.P.Val)
    )
  
}))

# =============================================================================
# FigS4I
# =============================================================================
nmd3_step("FigS4I", c("Prepare_KPT_RNA_Protein"), quote({
  # FigS4I — analysis below
  # KPT-330 vs DMSO cytoplasmic RNA-seq GSEA
  #
  # GO:0002183 - Cytoplasmic translational initiation
  # GO:0042255 - Ribosomal subunit assembly
  #
  # Reuses res_rna_df generated in Figure 2
  # Ranking metric: apeglm-shrunken log2FC
  # Positive: increased in KPT-330
  # Negative: decreased in KPT-330
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(msigdbr)
    library(clusterProfiler)
    library(ggplot2)
  })
  
  stopifnot(
    exists("res_rna_df"),
    all(
      c("gene", "log2FC_RNA") %in%
        colnames(res_rna_df)
    )
  )
  
  # ------------------------------------------------------------
  # 1. Build the ranked gene list
  # ------------------------------------------------------------
  
  ranking_df_S3C <- res_rna_df %>%
    dplyr::transmute(
      gene = toupper(as.character(gene)),
      log2FC = as.numeric(log2FC_RNA)
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      is.finite(log2FC)
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(
      order_by = abs(log2FC),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
  
  gene_ranking_S3C <- ranking_df_S3C$log2FC
  names(gene_ranking_S3C) <- ranking_df_S3C$gene
  
  gene_ranking_S3C <- sort(
    gene_ranking_S3C,
    decreasing = TRUE
  )
  
  cat(
    "Genes used in Figure S3C:",
    length(gene_ranking_S3C),
    "\n"
  )
  
  cat(
    "Ranking range:",
    round(min(gene_ranking_S3C), 3),
    "to",
    round(max(gene_ranking_S3C), 3),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 2. Retrieve the two GO:BP gene sets
  # ------------------------------------------------------------
  
  go_ids_S3C <- c(
    "GO:0002183",
    "GO:0042255"
  )
  
  
  gobp_sets_S3C <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  )
  
  term2gene_S3C <- gobp_sets_S3C %>%
    dplyr::filter(
      gs_exact_source %in% go_ids_S3C
    ) %>%
    dplyr::transmute(
      term = gs_exact_source,
      gene = toupper(gene_symbol)
    ) %>%
    dplyr::distinct()
  
  stopifnot(
    all(
      go_ids_S3C %in%
        term2gene_S3C$term
    )
  )
  
  cat("Genes in each GO pathway:\n")
  
  print(
    table(
      term2gene_S3C$term
    )
  )
  
  # ------------------------------------------------------------
  # 3. Standard GSEA analysis
  # ------------------------------------------------------------
  
  set.seed(123)
  
  gsea_S3C <- clusterProfiler::GSEA(
    geneList = gene_ranking_S3C,
    TERM2GENE = term2gene_S3C,
    exponent = 1,
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    verbose = FALSE,
    seed = TRUE
  )
  
  stopifnot(
    nrow(gsea_S3C@result) == 2
  )
  
  cat("\nFigure S3C GSEA results:\n")
  
  print(
    gsea_S3C@result[
      ,
      c(
        "ID",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
  # ------------------------------------------------------------
  # 4. Standard weighted running enrichment-score function
  # ------------------------------------------------------------
  
  compute_weighted_ES_S3C <- function(
    ranked_stats,
    pathway_genes,
    exponent = 1
  ) {
    
    ranked_stats <- sort(
      ranked_stats,
      decreasing = TRUE
    )
    
    ranked_genes <- names(ranked_stats)
    
    pathway_genes <- intersect(
      toupper(pathway_genes),
      ranked_genes
    )
    
    hits <- ranked_genes %in%
      pathway_genes
    
    number_hits <- sum(hits)
    number_genes <- length(ranked_stats)
    
    if (
      number_hits == 0 ||
      number_hits == number_genes
    ) {
      stop(
        "The pathway contains zero genes or all ranked genes."
      )
    }
    
    hit_weights <- abs(
      ranked_stats
    )^exponent
    
    hit_weights[!hits] <- 0
    
    hit_weights <- hit_weights /
      sum(hit_weights)
    
    miss_weights <- ifelse(
      hits,
      0,
      1 / (number_genes - number_hits)
    )
    
    running_ES <- cumsum(
      hit_weights - miss_weights
    )
    
    data.frame(
      Rank = seq_len(number_genes),
      ES = running_ES,
      Hit = hits,
      stringsAsFactors = FALSE
    )
  }
  
  # ------------------------------------------------------------
  # 5. Calculate running ES curves
  # ------------------------------------------------------------
  
  pathway_genes_S3C <- split(
    term2gene_S3C$gene,
    term2gene_S3C$term
  )
  
  curve_initiation_S3C <- compute_weighted_ES_S3C(
    ranked_stats = gene_ranking_S3C,
    pathway_genes =
      pathway_genes_S3C[["GO:0002183"]],
    exponent = 1
  ) %>%
    dplyr::mutate(
      Pathway =
        "Cytoplasmic translational initiation"
    )
  
  curve_assembly_S3C <- compute_weighted_ES_S3C(
    ranked_stats = gene_ranking_S3C,
    pathway_genes =
      pathway_genes_S3C[["GO:0042255"]],
    exponent = 1
  ) %>%
    dplyr::mutate(
      Pathway =
        "Ribosomal subunit assembly"
    )
  
  curve_data_S3C <- dplyr::bind_rows(
    curve_initiation_S3C,
    curve_assembly_S3C
  )
  
  # Verify curve ES values
  
  get_curve_ES_S3C <- function(x) {
    
    maximum_ES <- max(
      x,
      na.rm = TRUE
    )
    
    minimum_ES <- min(
      x,
      na.rm = TRUE
    )
    
    if (
      abs(maximum_ES) >= abs(minimum_ES)
    ) {
      maximum_ES
    } else {
      minimum_ES
    }
  }
  
  cat(
    "\nRunning-curve ES for GO:0002183:",
    get_curve_ES_S3C(
      curve_initiation_S3C$ES
    ),
    "\n"
  )
  
  cat(
    "Running-curve ES for GO:0042255:",
    get_curve_ES_S3C(
      curve_assembly_S3C$ES
    ),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 6. Figure S3C
  # ------------------------------------------------------------
  
  p_S3C <- ggplot2::ggplot(
    curve_data_S3C,
    ggplot2::aes(
      x = Rank,
      y = ES
    )
  ) +
    ggplot2::geom_line(
      data = curve_initiation_S3C,
      color = "#D73027",
      linewidth = 1.1
    ) +
    ggplot2::geom_line(
      data = curve_assembly_S3C,
      color = "#2166AC",
      linewidth = 1.1
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::labs(
      title = NULL,
      x = "Rank in ordered gene list",
      y = "Running enrichment score"
    ) +
    ggplot2::theme_classic(
      base_size = 13
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        color = "black",
        face = "bold"
      ),
      legend.position = "none",
      panel.grid = ggplot2::element_blank()
    )
  
  
  # ------------------------------------------------------------
  # 7. Save Figure S3C
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS4I_KPT330_vs_DMSO_Two_GO_GSEA.pdf",
    plot = p_S3C,
    width = 10,
    height = 7,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  
  
  
  
}))

# =============================================================================
# FigS5A
# =============================================================================
nmd3_step("FigS5A", character(), quote({
  # FigS5A — analysis below
  # KPT-330 vs DMSO volcano plots in nucleus and cytoplasm
  #
  # Red labels: eight predefined KPT-330-decreased genes
  # Blue labels: top 10 KPT-330-increased genes in each compartment
  # Purple labels: genes increased in nucleus and decreased in cytoplasm
  #
  # Reuses expr_matrix generated in Figure 2
  # ============================================================
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
  })
  
  stopifnot(
    exists("expr_matrix"),
    "Length" %in% colnames(expr_matrix)
  )
  
  # ------------------------------------------------------------
  # 1. Prepare the raw count matrix
  # ------------------------------------------------------------
  
  count_matrix_S2C <- expr_matrix[
    ,
    setdiff(colnames(expr_matrix), "Length"),
    drop = FALSE
  ]
  
  count_matrix_S2C <- as.matrix(count_matrix_S2C)
  storage.mode(count_matrix_S2C) <- "integer"
  
  # ------------------------------------------------------------
  # 2. Sample annotation
  # ------------------------------------------------------------
  
  sample_info_S2C <- data.frame(
    Sample = colnames(count_matrix_S2C),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      grepl("^[DX][0-9]+[CN]$", Sample)
    ) %>%
    dplyr::mutate(
      Prefix = substr(Sample, 1, 1),
      Compartment = dplyr::case_when(
        endsWith(Sample, "N") ~ "Nucleus",
        endsWith(Sample, "C") ~ "Cytoplasm",
        TRUE ~ NA_character_
      ),
      Treatment = dplyr::recode(
        Prefix,
        "D" = "DMSO",
        "X" = "KPT-330"
      )
    )
  
  cat("Samples used in Figure S2C:\n")
  print(
    table(
      sample_info_S2C$Treatment,
      sample_info_S2C$Compartment
    )
  )
  
  stopifnot(
    all(
      table(
        sample_info_S2C$Treatment,
        sample_info_S2C$Compartment
      ) >= 2
    )
  )
  
  # ------------------------------------------------------------
  # 3. Run DESeq2 separately for each compartment
  # ------------------------------------------------------------
  
  run_deseq2_S2C <- function(compartment_name) {
    
    sample_table <- sample_info_S2C %>%
      dplyr::filter(
        Compartment == compartment_name
      )
    
    samples <- sample_table$Sample
    
    counts <- count_matrix_S2C[
      ,
      samples,
      drop = FALSE
    ]
    
    storage.mode(counts) <- "integer"
    
    coldata <- data.frame(
      row.names = samples,
      condition = factor(
        sample_table$Treatment,
        levels = c(
          "DMSO",
          "KPT-330"
        )
      )
    )
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = coldata,
      design = ~ condition
    )
    
    dds <- DESeq2::DESeq(
      dds,
      quiet = TRUE
    )
    
    result <- DESeq2::results(
      dds,
      contrast = c(
        "condition",
        "KPT-330",
        "DMSO"
      )
    )
    
    result_df <- as.data.frame(result) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::transmute(
        gene = toupper(as.character(gene)),
        baseMean = as.numeric(baseMean),
        log2FoldChange = dplyr::if_else(
          is.na(log2FoldChange),
          0,
          as.numeric(log2FoldChange)
        ),
        pvalue = dplyr::if_else(
          is.na(pvalue),
          1,
          as.numeric(pvalue)
        ),
        padj = dplyr::if_else(
          is.na(padj),
          1,
          as.numeric(padj)
        )
      )
    
    result_df
  }
  
  res_nucleus_S2C <- run_deseq2_S2C(
    "Nucleus"
  )
  
  res_cytoplasm_S2C <- run_deseq2_S2C(
    "Cytoplasm"
  )
  
  # ------------------------------------------------------------
  # 4. Check the five nucleus-increased/cytoplasm-decreased genes
  # ------------------------------------------------------------
  
  cross_genes_S2C <- c(
    "SHROOM2",
    "HIC2",
    "RMRP",
    "TLN1",
    "COL1A1"
  )
  
  cross_check_S2C <- res_nucleus_S2C %>%
    dplyr::select(
      gene,
      log2FC_Nucleus = log2FoldChange,
      padj_Nucleus = padj
    ) %>%
    dplyr::inner_join(
      res_cytoplasm_S2C %>%
        dplyr::select(
          gene,
          log2FC_Cytoplasm = log2FoldChange,
          padj_Cytoplasm = padj
        ),
      by = "gene"
    ) %>%
    dplyr::filter(
      gene %in% cross_genes_S2C
    ) %>%
    dplyr::mutate(
      Nucleus_up_Cytoplasm_down = (
        log2FC_Nucleus > 0 &
          log2FC_Cytoplasm < 0
      )
    ) %>%
    dplyr::arrange(
      match(gene, cross_genes_S2C)
    )
  
  cat(
    "\nFive nucleus-increased/cytoplasm-decreased genes:\n"
  )
  
  print(cross_check_S2C)
  
  if (
    !all(cross_check_S2C$Nucleus_up_Cytoplasm_down)
  ) {
    warning(
      "At least one of the five genes does not show ",
      "positive nuclear log2FC and negative cytoplasmic log2FC."
    )
  }
  
  # ------------------------------------------------------------
  # 5. Label definitions
  # ------------------------------------------------------------
  
  red_genes_S2C <- c(
    "BRCA2",
    "EZH2",
    "AURKA",
    "UBE2C",
    "AURKB",
    "FOXM1",
    "BIRC5",
    "BRCA1"
  )
  
  lfc_threshold_S2C <- 1
  y_cap_S2C <- 300
  
  colors_S2C <- c(
    "KPT-330 increased" = "#4575B4",
    "KPT-330 decreased" = "#D73027",
    "NS" = "gray80"
  )
  
  # ------------------------------------------------------------
  # 6. Function for spreading fixed labels vertically
  # ------------------------------------------------------------
  
  spread_labels_S2C <- function(
    data,
    minimum_gap = 14,
    lower_limit = 5,
    upper_limit = 300
  ) {
    
    if (nrow(data) == 0) {
      return(data)
    }
    
    data <- data %>%
      dplyr::arrange(y_point) %>%
      dplyr::mutate(
        y_label = y_point
      )
    
    if (nrow(data) > 1) {
      for (i in 2:nrow(data)) {
        data$y_label[i] <- max(
          data$y_label[i],
          data$y_label[i - 1] + minimum_gap
        )
      }
    }
    
    if (max(data$y_label) > upper_limit) {
      data$y_label <- (
        data$y_label -
          (max(data$y_label) - upper_limit)
      )
    }
    
    if (min(data$y_label) < lower_limit) {
      data$y_label <- (
        data$y_label +
          (lower_limit - min(data$y_label))
      )
    }
    
    data
  }
  
  # ------------------------------------------------------------
  # 7. Volcano plot function
  # ------------------------------------------------------------
  
  make_volcano_S2C <- function(
    result_df,
    compartment_name,
    show_y_axis = TRUE
  ) {
    
    plot_df <- result_df %>%
      dplyr::mutate(
        minus_log10_padj = pmin(
          -log10(
            pmax(
              padj,
              1e-300
            )
          ),
          y_cap_S2C
        ),
        Direction = dplyr::case_when(
          log2FoldChange > lfc_threshold_S2C ~
            "KPT-330 increased",
          log2FoldChange < -lfc_threshold_S2C ~
            "KPT-330 decreased",
          TRUE ~ "NS"
        ),
        Direction = factor(
          Direction,
          levels = c(
            "KPT-330 increased",
            "KPT-330 decreased",
            "NS"
          )
        )
      )
    
    # Fixed eight red genes
    
    red_label_df <- plot_df %>%
      dplyr::filter(
        gene %in% red_genes_S2C
      ) %>%
      dplyr::transmute(
        gene,
        x_point = log2FoldChange,
        y_point = minus_log10_padj
      )
    
    missing_red <- setdiff(
      red_genes_S2C,
      red_label_df$gene
    )
    
    if (length(missing_red) > 0) {
      warning(
        compartment_name,
        ": missing red-label genes: ",
        paste(missing_red, collapse = ", ")
      )
    }
    
    red_label_df <- spread_labels_S2C(
      red_label_df,
      minimum_gap = 14,
      lower_limit = 20,
      upper_limit = y_cap_S2C
    ) %>%
      dplyr::mutate(
        x_label = -4.1
      )
    
    # Top 10 KPT-330-increased genes
    
    blue_top10_df <- plot_df %>%
      dplyr::filter(
        log2FoldChange > lfc_threshold_S2C
      ) %>%
      dplyr::arrange(
        padj,
        dplyr::desc(abs(log2FoldChange))
      ) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::anti_join(
        red_label_df %>%
          dplyr::select(gene),
        by = "gene"
      ) %>%
      dplyr::transmute(
        gene,
        x_point = log2FoldChange,
        y_point = minus_log10_padj
      )
    
    blue_top10_df <- spread_labels_S2C(
      blue_top10_df,
      minimum_gap = 14,
      lower_limit = 20,
      upper_limit = y_cap_S2C
    ) %>%
      dplyr::mutate(
        x_label = 5.4
      )
    
    cat(
      "\n",
      compartment_name,
      " blue top 10 genes:\n",
      sep = ""
    )
    
    print(blue_top10_df$gene)
    
    # Five nucleus-up/cytoplasm-down genes
    
    cross_label_df <- plot_df %>%
      dplyr::filter(
        gene %in% cross_genes_S2C
      ) %>%
      dplyr::mutate(
        label_order = match(
          gene,
          cross_genes_S2C
        )
      ) %>%
      dplyr::arrange(label_order) %>%
      dplyr::transmute(
        gene,
        x_point = log2FoldChange,
        y_point = minus_log10_padj,
        label_order
      ) %>%
      dplyr::mutate(
        x_label = -0.90,
        y_label = c(
          82,
          68,
          54,
          40,
          24
        )[label_order]
      )
    
    # Plot
    
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = log2FoldChange,
        y = minus_log10_padj
      )
    ) +
      ggplot2::geom_hline(
        yintercept = -log10(0.05),
        linetype = "dotted",
        color = "gray35",
        linewidth = 0.4
      ) +
      ggplot2::geom_vline(
        xintercept = c(
          -lfc_threshold_S2C,
          lfc_threshold_S2C
        ),
        linetype = "dotted",
        color = "gray35",
        linewidth = 0.4
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          color = Direction
        ),
        size = 2.2,
        alpha = 0.80
      ) +
      ggplot2::geom_rect(
        data = data.frame(
          xmin = -1.45,
          xmax = 1.25,
          ymin = 0,
          ymax = 90
        ),
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = NA,
        color = "#8E24AA",
        linewidth = 0.8
      ) +
      
      # Eight red genes
      
      ggplot2::geom_point(
        data = red_label_df,
        ggplot2::aes(
          x = x_point,
          y = y_point
        ),
        inherit.aes = FALSE,
        color = "#D73027",
        size = 2.8
      ) +
      ggplot2::geom_segment(
        data = red_label_df,
        ggplot2::aes(
          x = x_point,
          y = y_point,
          xend = x_label,
          yend = y_label
        ),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.3
      ) +
      ggplot2::geom_text(
        data = red_label_df,
        ggplot2::aes(
          x = x_label,
          y = y_label,
          label = gene
        ),
        inherit.aes = FALSE,
        hjust = 1,
        size = 3.7,
        family = "Arial",
        fontface = "italic"
      ) +
      
      # Blue top 10 genes
      
      ggplot2::geom_point(
        data = blue_top10_df,
        ggplot2::aes(
          x = x_point,
          y = y_point
        ),
        inherit.aes = FALSE,
        color = "#4575B4",
        size = 2.8
      ) +
      ggplot2::geom_segment(
        data = blue_top10_df,
        ggplot2::aes(
          x = x_point,
          y = y_point,
          xend = x_label,
          yend = y_label
        ),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.3
      ) +
      ggplot2::geom_text(
        data = blue_top10_df,
        ggplot2::aes(
          x = x_label,
          y = y_label,
          label = gene
        ),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3.7,
        family = "Arial",
        fontface = "italic"
      ) +
      
      # Five compartment-shift genes
      
      ggplot2::geom_point(
        data = cross_label_df,
        ggplot2::aes(
          x = x_point,
          y = y_point
        ),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.7,
        size = 3
      ) +
      ggplot2::geom_segment(
        data = cross_label_df,
        ggplot2::aes(
          x = x_point,
          y = y_point,
          xend = x_label,
          yend = y_label
        ),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.3
      ) +
      ggplot2::geom_text(
        data = cross_label_df,
        ggplot2::aes(
          x = x_label,
          y = y_label,
          label = gene
        ),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3.4,
        family = "Arial",
        fontface = "italic"
      ) +
      ggplot2::scale_color_manual(
        values = colors_S2C,
        drop = FALSE
      ) +
      ggplot2::scale_x_continuous(
        limits = c(-6, 8),
        breaks = c(-4, 0, 4, 8),
        expand = ggplot2::expansion(
          mult = 0.02
        )
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, y_cap_S2C + 10),
        breaks = seq(
          0,
          y_cap_S2C,
          by = 50
        ),
        labels = function(x) {
          ifelse(
            x == y_cap_S2C,
            paste0(">", y_cap_S2C),
            x
          )
        },
        expand = ggplot2::expansion(
          mult = 0.01
        )
      ) +
      ggplot2::labs(
        title = compartment_name,
        x = expression(log[2]~"Fold change"),
        y = if (show_y_axis) {
          expression(-log[10]~"Adjusted p-value")
        } else {
          NULL
        },
        color = NULL
      ) +
      ggplot2::theme_classic(
        base_size = 12,
        base_family = "Arial"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 17,
          face = "bold",
          hjust = 0.5
        ),
        axis.title = ggplot2::element_text(
          size = 12,
          face = "bold"
        ),
        axis.text = ggplot2::element_text(
          size = 10,
          color = "black"
        ),
        legend.position = "top",
        legend.text = ggplot2::element_text(
          size = 10,
          face = "bold"
        ),
        panel.grid = ggplot2::element_blank()
      )
    
    p
  }
  
  # ------------------------------------------------------------
  # 8. Create nucleus and cytoplasm panels
  # ------------------------------------------------------------
  
  p_S2C_nucleus <- make_volcano_S2C(
    result_df = res_nucleus_S2C,
    compartment_name = "Nucleus",
    show_y_axis = TRUE
  )
  
  p_S2C_cytoplasm <- make_volcano_S2C(
    result_df = res_cytoplasm_S2C,
    compartment_name = "Cytoplasm",
    show_y_axis = FALSE
  )
  
  # ------------------------------------------------------------
  # 9. Combine Figure S2C
  # ------------------------------------------------------------
  
  p_S2C <- (
    p_S2C_nucleus +
      p_S2C_cytoplasm
  ) +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "top"
    )
  
  
  # ------------------------------------------------------------
  # 10. Save Figure S2C
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS5A_KPT330_Nucleus_Cytoplasm_Volcano.pdf",
    plot = p_S2C,
    width = 14,
    height = 7,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  # ------------------------------------------------------------
  # 11. Save differential-expression results
  # ------------------------------------------------------------
  
  utils::write.csv(
    res_nucleus_S2C,
    file = "FigureS2C_Nucleus_KPT330_vs_DMSO_DESeq2.csv",
    row.names = FALSE
  )
  
  utils::write.csv(
    res_cytoplasm_S2C,
    file = "FigureS2C_Cytoplasm_KPT330_vs_DMSO_DESeq2.csv",
    row.names = FALSE
  )
  
  utils::write.csv(
    cross_check_S2C,
    file = "FigureS2C_Nucleus_Up_Cytoplasm_Down_Genes.csv",
    row.names = FALSE
  )
  
  
}))

# =============================================================================
# Prepare_KPT_GSEA
# =============================================================================
nmd3_step("Prepare_KPT_GSEA", c("Prepare_KPT_RNA_Protein"), quote({
  # Prepare_KPT_GSEA — analysis below
  # KEGG prostate cancer GSEA
  # KPT-330 vs DMSO cytoplasmic RNA-seq
  #
  # Reuses res_rna_df generated in Figure 2
  # Ranking metric: apeglm-shrunken DESeq2 log2FC
  # Positive values: KPT-330 > DMSO
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(msigdbr)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
  })
  
  stopifnot(
    exists("res_rna_df"),
    all(
      c("gene", "log2FC_RNA") %in%
        colnames(res_rna_df)
    )
  )
  
  # ------------------------------------------------------------
  # 1. Build the ranked gene list
  # ------------------------------------------------------------
  
  ranking_df_S2B <- res_rna_df %>%
    dplyr::transmute(
      gene = toupper(as.character(gene)),
      log2FC = as.numeric(log2FC_RNA)
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      is.finite(log2FC)
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(
      order_by = abs(log2FC),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
  
  gene_ranking_S2B <- ranking_df_S2B$log2FC
  
  names(gene_ranking_S2B) <-
    ranking_df_S2B$gene
  
  gene_ranking_S2B <- sort(
    gene_ranking_S2B,
    decreasing = TRUE
  )
  
  cat(
    "Genes used for Figure S2B GSEA:",
    length(gene_ranking_S2B),
    "\n"
  )
  
  cat(
    "Ranking range:",
    round(min(gene_ranking_S2B), 3),
    "to",
    round(max(gene_ranking_S2B), 3),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 2. Retrieve the complete KEGG legacy collection
  # ------------------------------------------------------------
  
  kegg_sets_S2B <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY"
  )
  
  term2gene_kegg_S2B <- kegg_sets_S2B %>%
    dplyr::transmute(
      term = gs_name,
      gene = toupper(gene_symbol)
    ) %>%
    dplyr::distinct()
  
  stopifnot(
    "KEGG_PROSTATE_CANCER" %in%
      term2gene_kegg_S2B$term
  )
  
  cat(
    "KEGG legacy pathways:",
    dplyr::n_distinct(
      term2gene_kegg_S2B$term
    ),
    "\n"
  )
  
  cat(
    "Genes in KEGG_PROSTATE_CANCER:",
    sum(
      term2gene_kegg_S2B$term ==
        "KEGG_PROSTATE_CANCER"
    ),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 3. Run GSEA across the complete KEGG legacy collection
  # ------------------------------------------------------------
  
  set.seed(123)
  
  gsea_S2B <- clusterProfiler::GSEA(
    geneList = gene_ranking_S2B,
    TERM2GENE = term2gene_kegg_S2B,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    nPermSimple = 10000,
    verbose = FALSE,
    seed = TRUE
  )
  
  # ------------------------------------------------------------
  # 4. Extract the KEGG prostate cancer result
  # ------------------------------------------------------------
  
  result_S2B <- gsea_S2B@result %>%
    dplyr::filter(
      ID == "KEGG_PROSTATE_CANCER"
    )
  
  stopifnot(
    nrow(result_S2B) == 1
  )
  
  cat("\nFigure S2B result:\n")
  
  print(
    result_S2B[
      ,
      c(
        "ID",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
  # ------------------------------------------------------------
  # 5. GSEA plot
  # ------------------------------------------------------------
  
  p_S2B <- enrichplot::gseaplot2(
    gsea_S2B,
    geneSetID = "KEGG_PROSTATE_CANCER",
    subplots = 1:3,
    pvalue_table = FALSE,
    title = "",
    base_size = 18
  )
  
  
  # ------------------------------------------------------------
  
}))

# =============================================================================
# FigS5B_C
# =============================================================================
nmd3_step("FigS5B_C", c("Prepare_KPT_GSEA"), quote({
  # FigS5B_C — analysis below
  # KPT-330-suppressed cytoplasmic pathways
  #
  # FigS5B_C — analysis below
  # FigS5B_C — analysis below
  #
  # Reuses gene_ranking_S2B and gsea_S2B generated in Figure S2B
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(msigdbr)
    library(clusterProfiler)
    library(ggplot2)
    library(patchwork)
    library(scales)
    library(stringr)
  })
  
  stopifnot(
    exists("gene_ranking_S2B"),
    length(gene_ranking_S2B) > 0,
    exists("gsea_S2B")
  )
  
  cat(
    "Genes in the reused ranking:",
    length(gene_ranking_S2B),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 1. Retrieve GO:BP and Hallmark gene sets
  # ------------------------------------------------------------
  
  gobp_sets_S2DE <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  )
  
  hallmark_sets_S2DE <- nmd3_cached_msig(
    species = "Homo sapiens",
    collection = "H"
  )
  
  term2gene_gobp_S2DE <- gobp_sets_S2DE %>%
    dplyr::transmute(
      term = gs_name,
      gene = toupper(gene_symbol)
    ) %>%
    dplyr::distinct()
  
  term2gene_hallmark_S2DE <- hallmark_sets_S2DE %>%
    dplyr::transmute(
      term = gs_name,
      gene = toupper(gene_symbol)
    ) %>%
    dplyr::distinct()
  
  cat(
    "GO:BP pathways:",
    dplyr::n_distinct(term2gene_gobp_S2DE$term),
    "\n"
  )
  
  cat(
    "Hallmark pathways:",
    dplyr::n_distinct(term2gene_hallmark_S2DE$term),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 2. Run GO:BP GSEA
  # ------------------------------------------------------------
  
  set.seed(123)
  
  gsea_gobp_S2D <- clusterProfiler::GSEA(
    geneList = gene_ranking_S2B,
    TERM2GENE = term2gene_gobp_S2DE,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    nPermSimple = 10000,
    verbose = FALSE,
    seed = TRUE
  )
  
  # ------------------------------------------------------------
  # 3. Reuse the complete KEGG GSEA from Figure S2B
  # ------------------------------------------------------------
  
  gsea_kegg_S2D <- gsea_S2B
  
  # ------------------------------------------------------------
  # 4. Run Hallmark GSEA
  # ------------------------------------------------------------
  
  set.seed(123)
  
  gsea_hallmark_S2E <- clusterProfiler::GSEA(
    geneList = gene_ranking_S2B,
    TERM2GENE = term2gene_hallmark_S2DE,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    nPermSimple = 10000,
    verbose = FALSE,
    seed = TRUE
  )
  
  # ------------------------------------------------------------
  # 5. Function for pathway labels
  # ------------------------------------------------------------
  
  format_pathway_label_S2DE <- function(x) {
    
    x <- gsub(
      "^(GOBP_|KEGG_|HALLMARK_)",
      "",
      x
    )
    
    x <- gsub(
      "_",
      " ",
      x
    )
    
    x <- stringr::str_to_sentence(
      tolower(x)
    )
    
    replacements <- c(
      "\\bDna\\b" = "DNA",
      "\\bRna\\b" = "RNA",
      "\\bMyc\\b" = "MYC",
      "\\bE2f\\b" = "E2F",
      "\\bG2m\\b" = "G2M",
      "\\bMtorc1\\b" = "Mtorc1",
      "\\bKegg\\b" = "KEGG",
      "\\bGo\\b" = "GO"
    )
    
    for (pattern in names(replacements)) {
      x <- gsub(
        pattern,
        replacements[[pattern]],
        x
      )
    }
    
    x
  }
  
  # ------------------------------------------------------------
  # 6. Extract the top 10 KPT-330-suppressed pathways
  # ------------------------------------------------------------
  
  extract_suppressed_top10_S2DE <- function(
    gsea_object,
    database_name
  ) {
    
    result <- gsea_object@result %>%
      dplyr::filter(
        is.finite(NES),
        is.finite(p.adjust),
        NES < 0
      ) %>%
      dplyr::arrange(
        p.adjust,
        dplyr::desc(abs(NES))
      ) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::transmute(
        ID = as.character(ID),
        Description = as.character(Description),
        Pathway = format_pathway_label_S2DE(
          as.character(ID)
        ),
        NES = as.numeric(NES),
        NES_absolute = abs(as.numeric(NES)),
        pvalue = as.numeric(pvalue),
        FDR = as.numeric(p.adjust),
        neg_log10_FDR = -log10(
          pmax(
            as.numeric(p.adjust),
            1e-300
          )
        ),
        Gene_set_size = as.numeric(setSize),
        Database = database_name
      )
    
    result
  }
  
  top10_gobp_S2D <- extract_suppressed_top10_S2DE(
    gsea_gobp_S2D,
    "GO:BP"
  )
  
  top10_kegg_S2D <- extract_suppressed_top10_S2DE(
    gsea_kegg_S2D,
    "KEGG"
  )
  
  top10_hallmark_S2E <- extract_suppressed_top10_S2DE(
    gsea_hallmark_S2E,
    "Hallmark"
  )
  
  cat("\nFigure S2D — GO:BP top 10:\n")
  
  print(
    top10_gobp_S2D %>%
      dplyr::select(
        Pathway,
        NES,
        FDR
      )
  )
  
  cat("\nFigure S2D — KEGG top 10:\n")
  
  print(
    top10_kegg_S2D %>%
      dplyr::select(
        Pathway,
        NES,
        FDR
      )
  )
  
  cat("\nFigure S2E — Hallmark top 10:\n")
  
  print(
    top10_hallmark_S2E %>%
      dplyr::select(
        Pathway,
        NES,
        FDR
      )
  )
  
  # ============================================================
  # FigS5B_C — analysis below
  # GO:BP and KEGG
  # ============================================================
  
  plot_data_S2D <- dplyr::bind_rows(
    top10_gobp_S2D,
    top10_kegg_S2D
  ) %>%
    dplyr::mutate(
      Database = factor(
        Database,
        levels = c(
          "GO:BP",
          "KEGG"
        )
      )
    ) %>%
    dplyr::arrange(
      Database,
      NES_absolute
    ) %>%
    dplyr::mutate(
      Pathway_unique = paste(
        Pathway,
        Database,
        sep = "___"
      ),
      Pathway_unique = factor(
        Pathway_unique,
        levels = unique(Pathway_unique)
      )
    )
  
  color_limits_S2D <- range(
    plot_data_S2D$neg_log10_FDR,
    finite = TRUE
  )
  
  size_limits_S2D <- range(
    plot_data_S2D$Gene_set_size,
    finite = TRUE
  )
  
  nes_limits_S2D <- c(
    floor(
      min(
        plot_data_S2D$NES_absolute,
        na.rm = TRUE
      ) * 10
    ) / 10,
    ceiling(
      max(
        plot_data_S2D$NES_absolute,
        na.rm = TRUE
      ) * 10
    ) / 10
  )
  
  p_S2D <- ggplot2::ggplot(
    plot_data_S2D,
    ggplot2::aes(
      x = NES_absolute,
      y = Pathway_unique
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = Gene_set_size,
        color = neg_log10_FDR
      ),
      alpha = 0.95
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Database),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) {
        sub(
          "___.*$",
          "",
          x
        )
      }
    ) +
    ggplot2::scale_color_gradient(
      low = "#3B82F6",
      high = "#E63946",
      limits = color_limits_S2D,
      oob = scales::squish,
      name = expression(-log[10]~FDR)
    ) +
    ggplot2::scale_size_continuous(
      range = c(3, 10),
      limits = size_limits_S2D,
      name = "Gene set size"
    ) +
    ggplot2::scale_x_continuous(
      limits = nes_limits_S2D,
      expand = ggplot2::expansion(
        mult = c(0.04, 0.08)
      )
    ) +
    ggplot2::labs(
      x = "Normalized Enrichment Score (|NES|)",
      y = NULL,
      title = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 13,
      base_family = "Arial"
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        color = "black",
        size = 11
      ),
      axis.text.y = ggplot2::element_text(
        color = "black",
        size = 11,
        face = "bold"
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 13,
        face = "bold"
      ),
      strip.background = ggplot2::element_rect(
        fill = "gray90",
        color = "black",
        linewidth = 0.8
      ),
      strip.text.y.right = ggplot2::element_text(
        angle = 270,
        face = "bold",
        color = "black"
      ),
      panel.spacing.y = grid::unit(
        0.10,
        "lines"
      ),
      legend.position = "right",
      legend.title = ggplot2::element_text(
        size = 10
      ),
      legend.text = ggplot2::element_text(
        size = 9
      )
    )
  
  
  nmd3_save_pdf(
    filename = "FigS5B_KPT330_Suppressed_GOBP_KEGG_Pathways.pdf",
    plot = p_S2D,
    width = 10,
    height = 8,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  # ============================================================
  # FigS5B_C — analysis below
  # Hallmark pathways
  # ============================================================
  
  highlight_hallmark_S2E <- c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_SPERMATOGENESIS"
  )
  
  plot_data_S2E <- top10_hallmark_S2E %>%
    dplyr::arrange(
      NES_absolute
    ) %>%
    dplyr::mutate(
      Pathway = factor(
        Pathway,
        levels = unique(Pathway)
      ),
      Label_color = dplyr::if_else(
        ID %in% highlight_hallmark_S2E,
        "red",
        "black"
      )
    )
  
  color_limits_S2E <- range(
    plot_data_S2E$neg_log10_FDR,
    finite = TRUE
  )
  
  size_limits_S2E <- range(
    plot_data_S2E$Gene_set_size,
    finite = TRUE
  )
  
  nes_lower_S2E <- floor(
    min(
      plot_data_S2E$NES_absolute,
      na.rm = TRUE
    ) * 10
  ) / 10
  
  nes_upper_S2E <- ceiling(
    max(
      plot_data_S2E$NES_absolute,
      na.rm = TRUE
    ) * 10
  ) / 10
  
  # Label position outside the left side of the plotting panel
  
  label_x_S2E <- nes_lower_S2E - 0.08
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  p_S2E <- ggplot2::ggplot(
    plot_data_S2E,
    ggplot2::aes(
      x = NES_absolute,
      y = Pathway
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = Gene_set_size,
        color = neg_log10_FDR
      ),
      alpha = 0.95
    ) +
    
    # Black pathway labels
    
    ggplot2::geom_text(
      data = plot_data_S2E %>%
        dplyr::filter(Label_color == "black"),
      ggplot2::aes(
        x = label_x_S2E,
        y = Pathway,
        label = Pathway
      ),
      inherit.aes = FALSE,
      hjust = 1,
      size = 4.1,
      family = "Arial",
      fontface = "bold",
      color = "black"
    ) +
    
    # Red highlighted pathway labels
    
    ggplot2::geom_text(
      data = plot_data_S2E %>%
        dplyr::filter(Label_color == "red"),
      ggplot2::aes(
        x = label_x_S2E,
        y = Pathway,
        label = Pathway
      ),
      inherit.aes = FALSE,
      hjust = 1,
      size = 4.1,
      family = "Arial",
      fontface = "bold",
      color = "#D00000"
    ) +
    ggplot2::scale_color_gradient(
      low = "#0000FF",
      high = "#FF1A1A",
      limits = color_limits_S2E,
      oob = scales::squish,
      name = expression(-log[10]~FDR)
    ) +
    ggplot2::scale_size_continuous(
      range = c(3, 10),
      limits = size_limits_S2E,
      breaks = pretty(
        size_limits_S2E,
        n = 3
      ),
      name = "Gene set size"
    ) +
    ggplot2::scale_x_continuous(
      breaks = pretty(
        c(
          nes_lower_S2E,
          nes_upper_S2E
        ),
        n = 5
      ),
      expand = ggplot2::expansion(
        mult = c(0.02, 0.08)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        nes_lower_S2E,
        nes_upper_S2E
      ),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "KPT-330 suppressed cytoplasmic\npathways",
      x = "Normalized Enrichment Score (|NES|)",
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 13,
      base_family = "Arial"
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        color = "black",
        size = 11
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 13,
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(
        size = 10
      ),
      legend.text = ggplot2::element_text(
        size = 9
      ),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 10,
        b = 10,
        l = 180
      )
    )
  
  
  nmd3_save_pdf(
    filename = "FigS5C_KPT330_Suppressed_Hallmark_Pathways.pdf",
    plot = p_S2E,
    width = 9,
    height = 6.5,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  
}))

# =============================================================================
# FigS6A
# =============================================================================
nmd3_step("FigS6A", character(), quote({
  # FigS6A — analysis below
  # PCA of all eight RNA-seq and proteomics groups
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggforce)
    library(patchwork)
  })
  
  # ------------------------------------------------------------
  # 1. Group definitions
  # ------------------------------------------------------------
  
  group_levels_S5A <- c(
    "Ctrl-Cytoplasm",
    "KPT-330-Cytoplasm",
    "Ctrl-Nucleus",
    "KPT-330-Nucleus",
    "siCtrl-Cytoplasm",
    "siNMD3-Cytoplasm",
    "siCtrl-Nucleus",
    "siNMD3-Nucleus"
  )
  
  group_colors_S5A <- c(
    "Ctrl-Cytoplasm" = "#E31A1C",
    "Ctrl-Nucleus" = "#E31A1C",
    "KPT-330-Cytoplasm" = "#1F78B4",
    "KPT-330-Nucleus" = "#1F78B4",
    "siCtrl-Cytoplasm" = "#33A02C",
    "siCtrl-Nucleus" = "#33A02C",
    "siNMD3-Cytoplasm" = "#6A3D9A",
    "siNMD3-Nucleus" = "#6A3D9A"
  )
  
  group_shapes_S5A <- c(
    "Ctrl-Cytoplasm" = 15,
    "KPT-330-Cytoplasm" = 15,
    "siCtrl-Cytoplasm" = 15,
    "siNMD3-Cytoplasm" = 15,
    "Ctrl-Nucleus" = 17,
    "KPT-330-Nucleus" = 17,
    "siCtrl-Nucleus" = 17,
    "siNMD3-Nucleus" = 17
  )
  
  # ============================================================
  # 2. RNA-seq PCA
  # ============================================================
  
  stopifnot(
    exists("expr_matrix"),
    "Length" %in% colnames(expr_matrix)
  )
  
  rna_sample_columns_S5A <- grep(
    "^[CDNX][0-9]+[CN]$",
    colnames(expr_matrix),
    value = TRUE
  )
  
  stopifnot(
    length(rna_sample_columns_S5A) >= 8
  )
  
  cat(
    "RNA-seq samples used:",
    length(rna_sample_columns_S5A),
    "\n"
  )
  
  print(rna_sample_columns_S5A)
  
  # ------------------------------------------------------------
  # 2.1 Calculate TPM
  # ------------------------------------------------------------
  
  gene_length_kb_S5A <- suppressWarnings(
    as.numeric(expr_matrix$Length) / 1000
  )
  
  rna_counts_S5A <- as.matrix(
    expr_matrix[
      ,
      rna_sample_columns_S5A,
      drop = FALSE
    ]
  )
  
  storage.mode(rna_counts_S5A) <- "numeric"
  
  valid_length_S5A <- is.finite(
    gene_length_kb_S5A
  ) &
    gene_length_kb_S5A > 0
  
  rna_counts_S5A <- rna_counts_S5A[
    valid_length_S5A,
    ,
    drop = FALSE
  ]
  
  gene_length_kb_S5A <- gene_length_kb_S5A[
    valid_length_S5A
  ]
  
  rna_rpk_S5A <- sweep(
    rna_counts_S5A,
    1,
    gene_length_kb_S5A,
    FUN = "/"
  )
  
  rna_scaling_factors_S5A <- colSums(
    rna_rpk_S5A,
    na.rm = TRUE
  )
  
  stopifnot(
    all(
      is.finite(rna_scaling_factors_S5A)
    ),
    all(
      rna_scaling_factors_S5A > 0
    )
  )
  
  rna_tpm_S5A <- sweep(
    rna_rpk_S5A,
    2,
    rna_scaling_factors_S5A,
    FUN = "/"
  ) * 1e6
  
  rna_log2_tpm_S5A <- log2(
    rna_tpm_S5A + 1
  )
  
  # ------------------------------------------------------------
  # 2.2 Filter genes
  # ------------------------------------------------------------
  
  keep_rna_genes_S5A <- apply(
    rna_log2_tpm_S5A,
    1,
    function(x) {
      all(is.finite(x)) &&
        stats::sd(x) > 1e-6
    }
  )
  
  rna_pca_matrix_S5A <- rna_log2_tpm_S5A[
    keep_rna_genes_S5A,
    ,
    drop = FALSE
  ]
  
  cat(
    "RNA-seq genes used for PCA:",
    nrow(rna_pca_matrix_S5A),
    "\n"
  )
  
  stopifnot(
    nrow(rna_pca_matrix_S5A) > 2,
    ncol(rna_pca_matrix_S5A) > 2
  )
  
  # ------------------------------------------------------------
  # 2.3 Run RNA PCA
  # ------------------------------------------------------------
  
  pca_rna_S5A <- stats::prcomp(
    t(rna_pca_matrix_S5A),
    center = TRUE,
    scale. = TRUE
  )
  
  rna_variance_S5A <- (
    pca_rna_S5A$sdev^2 /
      sum(pca_rna_S5A$sdev^2)
  ) * 100
  
  pca_rna_df_S5A <- as.data.frame(
    pca_rna_S5A$x[
      ,
      1:2,
      drop = FALSE
    ]
  ) %>%
    tibble::rownames_to_column(
      "Sample"
    ) %>%
    dplyr::mutate(
      Prefix = substr(
        Sample,
        1,
        1
      ),
      Suffix = substr(
        Sample,
        nchar(Sample),
        nchar(Sample)
      ),
      Treatment = dplyr::case_when(
        Prefix == "D" ~ "Ctrl",
        Prefix == "X" ~ "KPT-330",
        Prefix == "C" ~ "siCtrl",
        Prefix == "N" ~ "siNMD3",
        TRUE ~ NA_character_
      ),
      Compartment = dplyr::case_when(
        Suffix == "C" ~ "Cytoplasm",
        Suffix == "N" ~ "Nucleus",
        TRUE ~ NA_character_
      ),
      GroupLabel = paste0(
        Treatment,
        "-",
        Compartment
      ),
      GroupLabel = factor(
        GroupLabel,
        levels = group_levels_S5A
      )
    ) %>%
    dplyr::filter(
      !is.na(GroupLabel)
    )
  
  cat("RNA-seq groups:\n")
  
  print(
    table(
      pca_rna_df_S5A$GroupLabel
    )
  )
  
  # ============================================================
  # 3. Proteomics PCA
  # ============================================================
  
  if (exists("prot_raw")) {
    
    protein_source_S5A <- prot_raw
    
  } else if (exists("protein_data")) {
    
    protein_source_S5A <- protein_data
    
  } else {
    
    stop(
      "Neither prot_raw nor protein_data exists. Run the proteomics loading section first."
    )
  }
  
  stopifnot(
    "Gene" %in% colnames(protein_source_S5A)
  )
  
  all_normalized_columns_S5A <- grep(
    "^Norm_",
    colnames(protein_source_S5A),
    value = TRUE
  )
  
  # ------------------------------------------------------------
  # 3.1 Identify the 24 samples from the four treatments
  # ------------------------------------------------------------
  
  protein_sample_information_S5A <- data.frame(
    Sample = all_normalized_columns_S5A,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Treatment = dplyr::case_when(
        grepl(
          "siRNAxNMD3xKD",
          Sample,
          ignore.case = TRUE
        ) ~ "siNMD3",
        
        grepl(
          "siRNAxCtrl",
          Sample,
          ignore.case = TRUE
        ) ~ "siCtrl",
        
        grepl(
          "XPO1inhibitor|KPT|Selinexor",
          Sample,
          ignore.case = TRUE
        ) ~ "KPT-330",
        
        grepl(
          "DMSO",
          Sample,
          ignore.case = TRUE
        ) ~ "Ctrl",
        
        TRUE ~ NA_character_
      ),
      
      Compartment = dplyr::case_when(
        grepl(
          "Cyto|Cyt",
          Sample,
          ignore.case = TRUE
        ) ~ "Cytoplasm",
        
        grepl(
          "Nuc|Nucleus",
          Sample,
          ignore.case = TRUE
        ) ~ "Nucleus",
        
        TRUE ~ NA_character_
      ),
      
      GroupLabel = paste0(
        Treatment,
        "-",
        Compartment
      ),
      
      GroupLabel = factor(
        GroupLabel,
        levels = group_levels_S5A
      )
    ) %>%
    dplyr::filter(
      !is.na(Treatment),
      !is.na(Compartment),
      !is.na(GroupLabel)
    )
  
  protein_sample_columns_S5A <-
    protein_sample_information_S5A$Sample
  
  cat(
    "Proteomics samples used:",
    length(protein_sample_columns_S5A),
    "\n"
  )
  
  cat("Proteomics groups:\n")
  
  print(
    table(
      protein_sample_information_S5A$GroupLabel
    )
  )
  
  stopifnot(
    length(protein_sample_columns_S5A) >= 8
  )
  
  # ------------------------------------------------------------
  # 3.2 Construct the proteomics expression matrix
  # ------------------------------------------------------------
  
  protein_table_S5A <- protein_source_S5A %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != ""
    ) %>%
    dplyr::select(
      Gene,
      dplyr::all_of(
        protein_sample_columns_S5A
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(
          protein_sample_columns_S5A
        ),
        ~ suppressWarnings(
          as.numeric(.x)
        )
      )
    )
  
  protein_pca_matrix_S5A <- as.matrix(
    protein_table_S5A[
      ,
      protein_sample_columns_S5A,
      drop = FALSE
    ]
  )
  
  rownames(protein_pca_matrix_S5A) <- make.unique(
    as.character(
      protein_table_S5A$Gene
    )
  )
  
  # Preserve the original missing-value handling
  
  protein_pca_matrix_S5A[
    is.na(protein_pca_matrix_S5A)
  ] <- 0
  
  # ------------------------------------------------------------
  # 3.3 Filter constant genes and samples
  # ------------------------------------------------------------
  
  keep_protein_genes_S5A <- apply(
    protein_pca_matrix_S5A,
    1,
    function(x) {
      all(is.finite(x)) &&
        stats::sd(x) > 0
    }
  )
  
  protein_pca_matrix_S5A <-
    protein_pca_matrix_S5A[
      keep_protein_genes_S5A,
      ,
      drop = FALSE
    ]
  
  keep_protein_samples_S5A <- apply(
    protein_pca_matrix_S5A,
    2,
    function(x) {
      all(is.finite(x)) &&
        stats::sd(x) > 0
    }
  )
  
  protein_pca_matrix_S5A <-
    protein_pca_matrix_S5A[
      ,
      keep_protein_samples_S5A,
      drop = FALSE
    ]
  
  cat(
    "Proteins used for PCA:",
    nrow(protein_pca_matrix_S5A),
    "\n"
  )
  
  stopifnot(
    nrow(protein_pca_matrix_S5A) > 2,
    ncol(protein_pca_matrix_S5A) > 2
  )
  
  # ------------------------------------------------------------
  # 3.4 Run proteomics PCA
  # ------------------------------------------------------------
  
  pca_protein_S5A <- stats::prcomp(
    t(protein_pca_matrix_S5A),
    center = TRUE,
    scale. = TRUE
  )
  
  protein_variance_S5A <- (
    pca_protein_S5A$sdev^2 /
      sum(pca_protein_S5A$sdev^2)
  ) * 100
  
  pca_protein_df_S5A <- as.data.frame(
    pca_protein_S5A$x[
      ,
      1:2,
      drop = FALSE
    ]
  ) %>%
    tibble::rownames_to_column(
      "Sample"
    ) %>%
    dplyr::left_join(
      protein_sample_information_S5A %>%
        dplyr::select(
          Sample,
          Treatment,
          Compartment,
          GroupLabel
        ),
      by = "Sample"
    ) %>%
    dplyr::mutate(
      GroupLabel = factor(
        GroupLabel,
        levels = group_levels_S5A
      )
    ) %>%
    dplyr::filter(
      !is.na(GroupLabel)
    )
  
  # ============================================================
  # 4. PCA plotting function
  # ============================================================
  
  make_pca_plot_S5A <- function(
    pca_data,
    variance_values,
    panel_title
  ) {
    
    ggplot(
      pca_data,
      aes(
        x = PC1,
        y = PC2,
        color = GroupLabel,
        shape = GroupLabel
      )
    ) +
      ggforce::geom_mark_ellipse(
        aes(
          group = GroupLabel,
          fill = GroupLabel
        ),
        alpha = 0.12,
        color = NA,
        expand = grid::unit(
          2,
          "mm"
        ),
        show.legend = FALSE
      ) +
      geom_point(
        size = 3.5,
        alpha = 0.95
      ) +
      scale_color_manual(
        values = group_colors_S5A,
        breaks = group_levels_S5A,
        name = NULL
      ) +
      scale_fill_manual(
        values = group_colors_S5A,
        breaks = group_levels_S5A,
        guide = "none"
      ) +
      scale_shape_manual(
        values = group_shapes_S5A,
        breaks = group_levels_S5A,
        name = NULL
      ) +
      guides(
        color = guide_legend(
          nrow = 4,
          byrow = TRUE,
          override.aes = list(
            size = 3.5,
            alpha = 1
          )
        ),
        shape = guide_legend(
          nrow = 4,
          byrow = TRUE
        )
      ) +
      labs(
        title = panel_title,
        x = sprintf(
          "PC1 (%.1f%%)",
          variance_values[1]
        ),
        y = sprintf(
          "PC2 (%.1f%%)",
          variance_values[2]
        )
      ) +
      theme_classic(
        base_size = 13
      ) +
      theme(
        aspect.ratio = 1,
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.8
        ),
        axis.line = element_blank(),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.6
        ),
        axis.text = element_text(
          color = "black",
          size = 11
        ),
        axis.title = element_text(
          color = "black",
          size = 13,
          face = "bold"
        ),
        plot.title = element_text(
          color = "black",
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        legend.position = "bottom",
        legend.text = element_text(
          size = 10
        ),
        plot.margin = margin(
          8,
          8,
          8,
          8
        )
      )
  }
  
  # ============================================================
  # 5. Generate RNA-seq and proteomics panels
  # ============================================================
  
  p_RNA_S5A <- make_pca_plot_S5A(
    pca_data = pca_rna_df_S5A,
    variance_values = rna_variance_S5A,
    panel_title = "RNA-seq"
  )
  
  p_Protein_S5A <- make_pca_plot_S5A(
    pca_data = pca_protein_df_S5A,
    variance_values = protein_variance_S5A,
    panel_title = "Proteomics"
  )
  
  
  
  # ============================================================
  # 6. Combine Figure S5A
  # ============================================================
  
  p_S5A <- (
    p_RNA_S5A /
      p_Protein_S5A
  ) +
    patchwork::plot_layout(
      guides = "collect",
      heights = c(1, 1)
    ) &
    theme(
      legend.position = "bottom"
    )
  
  
  # ============================================================
  # 7. Save
  # ============================================================
  
  nmd3_save_pdf(
    filename = "FigS6A_RNA_Proteomics_Eight_Group_PCA.pdf",
    plot = p_S5A,
    width = 7.5,
    height = 11,
    units = "in",
    device = cairo_pdf
  )
  
  write.csv(
    pca_rna_df_S5A,
    "FigureS5A_RNA_PCA_coordinates.csv",
    row.names = FALSE
  )
  
  write.csv(
    pca_protein_df_S5A,
    "FigureS5A_Proteomics_PCA_coordinates.csv",
    row.names = FALSE
  )
  
  variance_S5A <- data.frame(
    Dataset = c(
      "RNA-seq",
      "Proteomics"
    ),
    PC1_percent = c(
      rna_variance_S5A[1],
      protein_variance_S5A[1]
    ),
    PC2_percent = c(
      rna_variance_S5A[2],
      protein_variance_S5A[2]
    )
  )
  
  print(variance_S5A)
  
  write.csv(
    variance_S5A,
    "FigureS5A_PCA_variance_explained.csv",
    row.names = FALSE
  )
  
  
}))

# =============================================================================
# FigS6B
# =============================================================================
nmd3_step("FigS6B", c("Prepare_Fig3B_S6"), quote({
  # 3. Select RNA labels dynamically
  # ------------------------------------------------------------
  
  select_rna_labels_S5B <- function(
    result_table
  ) {
    
    top_increased <- result_table %>%
      dplyr::filter(
        Direction == "siNMD3 increased"
      ) %>%
      dplyr::arrange(
        padj,
        dplyr::desc(abs(log2FC))
      ) %>%
      dplyr::slice_head(
        n = rna_label_number_S5BC
      )
    
    top_decreased <- result_table %>%
      dplyr::filter(
        Direction == "siNMD3 decreased"
      ) %>%
      dplyr::arrange(
        padj,
        dplyr::desc(abs(log2FC))
      ) %>%
      dplyr::slice_head(
        n = rna_label_number_S5BC
      )
    
    # Labels are selected by the same significance definition as colors.
    dplyr::bind_rows(
      top_increased,
      top_decreased
    ) %>%
      dplyr::distinct(
        gene,
        .keep_all = TRUE
      ) %>%
      dplyr::mutate(
        LabelColor = "black"  # labels derive only from significance-ranked genes
      )
  }
  
  rna_labels_nucleus_S5B <-
    select_rna_labels_S5B(
      rna_nucleus_S5B
    )
  
  rna_labels_cytoplasm_S5B <-
    select_rna_labels_S5B(
      rna_cytoplasm_S5B
    )
  
  cat("RNA nucleus labels:\n")
  print(rna_labels_nucleus_S5B$gene)
  
  cat("RNA cytoplasm labels:\n")
  print(rna_labels_cytoplasm_S5B$gene)
  
  # ------------------------------------------------------------
  # 4. RNA volcano plotting function
  # ------------------------------------------------------------
  
  make_rna_volcano_S5B <- function(
    result_table,
    label_table,
    panel_title,
    show_x_title = TRUE,
    show_legend = TRUE
  ) {
    
    black_labels <- label_table %>%
      dplyr::filter(
        LabelColor == "black"
      )
    
    p <- ggplot(
      result_table,
      aes(
        x = log2FC,
        y = neg_log10_padj
      )
    ) +
      geom_hline(
        yintercept = -log10(
          rna_fdr_threshold_S5BC
        ),
        linetype = "dotted",
        color = "gray40",
        linewidth = 0.45
      ) +
      geom_vline(
        xintercept = c(
          -rna_lfc_threshold_S5BC,
          rna_lfc_threshold_S5BC
        ),
        linetype = "dotted",
        color = "gray40",
        linewidth = 0.45
      ) +
      geom_point(
        data = result_table %>%
          dplyr::filter(
            Direction == "NS"
          ),
        color = "gray80",
        size = 1.5,
        alpha = 0.75
      ) +
      geom_point(
        data = result_table %>%
          dplyr::filter(
            Direction != "NS"
          ),
        aes(
          color = Direction
        ),
        size = 1.7,
        alpha = 0.90
      ) +
      ggrepel::geom_text_repel(
        data = black_labels,
        aes(
          label = gene
        ),
        color = "black",
        size = 3.2,
        fontface = "italic",
        box.padding = 0.35,
        point.padding = 0.20,
        segment.color = "gray40",
        segment.linewidth = 0.25,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = volcano_colors_S5BC,
        breaks = c(
          "siNMD3 increased",
          "siNMD3 decreased",
          "NS"
        ),
        drop = FALSE,
        name = NULL
      ) +
      coord_cartesian(
        xlim = c(-6, 6),
        clip = "off"
      ) +
      labs(
        title = panel_title,
        x = if (show_x_title) {
          expression(log[2]~"Fold change")
        } else {
          NULL
        },
        y = expression(
          -log[10]~"adjusted p-value"
        )
      ) +
      theme_classic(
        base_size = 12
      ) +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(
          color = "black",
          linewidth = 0.7
        ),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.6
        ),
        axis.text = element_text(
          color = "black",
          size = 10
        ),
        axis.title = element_text(
          color = "black",
          size = 12,
          face = "bold"
        ),
        plot.title = element_text(
          color = "black",
          size = 14,
          face = "bold",
          hjust = 0.05
        ),
        legend.position = if (
          show_legend
        ) {
          "top"
        } else {
          "none"
        },
        legend.text = element_text(
          size = 9
        ),
        plot.margin = margin(
          8,
          8,
          8,
          8
        )
      )
    
    p
  }
  
  p_RNA_nucleus_S5B <- make_rna_volcano_S5B(
    result_table = rna_nucleus_S5B,
    label_table = rna_labels_nucleus_S5B,
    panel_title = "Nucleus",
    show_x_title = FALSE,
    show_legend = TRUE
  )
  
  p_RNA_cytoplasm_S5B <- make_rna_volcano_S5B(
    result_table = rna_cytoplasm_S5B,
    label_table = rna_labels_cytoplasm_S5B,
    panel_title = "Cytoplasm",
    show_x_title = TRUE,
    show_legend = FALSE
  )
  
  p_S5B <- (
    p_RNA_nucleus_S5B /
      p_RNA_cytoplasm_S5B
  ) +
    patchwork::plot_annotation(
      title = "RNA-seq",
      theme = theme(
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        )
      )
    )
  
  
  nmd3_save_pdf(
    filename = "FigS6B_siNMD3_RNA_Volcano_Nucleus_Cytoplasm.pdf",
    plot = p_S5B,
    width = 6.5,
    height = 9,
    units = "in",
    device = cairo_pdf
  )
  
}))

# =============================================================================
# FigS6C
# =============================================================================
nmd3_step("FigS6C", c("FigS6B"), quote({
  # FigS6C — analysis below
  # ============================================================
  
  # ------------------------------------------------------------
  # 5. Reuse the proteomics report
  # ------------------------------------------------------------
  
  if (exists("prot_raw")) {
    
    protein_source_S5C <- prot_raw
    
  } else if (exists("protein_data")) {
    
    protein_source_S5C <- protein_data
    
  } else {
    
    stop(
      "Neither prot_raw nor protein_data exists. Run the proteomics loading section first."
    )
  }
  
  stopifnot(
    "Gene" %in% colnames(protein_source_S5C)
  )
  
  # ------------------------------------------------------------
  # 6. Build proteomics result tables
  # Positive log2FC = higher in siNMD3
  # ------------------------------------------------------------
  
  make_protein_result_S5C <- function(
    compartment_code,
    compartment_label
  ) {
    
    log2fc_column <- paste0(
      "Log2FC_siRNAxNMD3xKD_",
      compartment_code,
      ".siRNAxCtrl_",
      compartment_code
    )
    
    pvalue_column <- paste0(
      "Tukey_siRNAxNMD3xKD_",
      compartment_code,
      ".siRNAxCtrl_",
      compartment_code
    )
    
    stopifnot(
      all(
        c(
          log2fc_column,
          pvalue_column
        ) %in% colnames(protein_source_S5C)
      )
    )
    
    result_table <- protein_source_S5C %>%
      dplyr::filter(
        !is.na(Gene),
        Gene != ""
      ) %>%
      dplyr::transmute(
        gene = as.character(Gene),
        log2FC = suppressWarnings(
          as.numeric(
            .data[[log2fc_column]]
          )
        ),
        pvalue = suppressWarnings(
          as.numeric(
            .data[[pvalue_column]]
          )
        ),
        Compartment = compartment_label
      ) %>%
      dplyr::filter(
        is.finite(log2FC),
        is.finite(pvalue)
      )
    
    minimum_positive_p <- suppressWarnings(
      min(
        result_table$pvalue[
          result_table$pvalue > 0
        ],
        na.rm = TRUE
      )
    )
    
    if (!is.finite(minimum_positive_p)) {
      minimum_positive_p <-
        .Machine$double.xmin
    }
    
    result_table %>%
      dplyr::mutate(
        pvalue = ifelse(
          pvalue <= 0,
          minimum_positive_p / 10,
          pvalue
        ),
        pvalue_plot = pmax(
          pvalue,
          .Machine$double.xmin
        ),
        neg_log10_pvalue = -log10(
          pvalue_plot
        ),
        Direction = dplyr::case_when(
          pvalue < protein_p_threshold_S5BC &
            log2FC > 0 ~
            "siNMD3 increased",
          
          pvalue < protein_p_threshold_S5BC &
            log2FC < 0 ~
            "siNMD3 decreased",
          
          TRUE ~ "NS"
        ),
        Direction = factor(
          Direction,
          levels = c(
            "siNMD3 increased",
            "siNMD3 decreased",
            "NS"
          )
        )
      )
  }
  
  protein_nucleus_S5C <- make_protein_result_S5C(
    compartment_code = "Nuc",
    compartment_label = "Nucleus"
  )
  
  protein_cytoplasm_S5C <- make_protein_result_S5C(
    compartment_code = "Cyto",
    compartment_label = "Cytoplasm"
  )
  
  # ------------------------------------------------------------
  # 7. Select protein labels dynamically
  # ------------------------------------------------------------
  
  select_protein_labels_S5C <- function(
    result_table
  ) {
    
    top_increased <- result_table %>%
      dplyr::filter(
        Direction == "siNMD3 increased"
      ) %>%
      dplyr::arrange(
        pvalue,
        dplyr::desc(abs(log2FC))
      ) %>%
      dplyr::slice_head(
        n = protein_label_number_S5BC
      )
    
    top_decreased <- result_table %>%
      dplyr::filter(
        Direction == "siNMD3 decreased"
      ) %>%
      dplyr::arrange(
        pvalue,
        dplyr::desc(abs(log2FC))
      ) %>%
      dplyr::slice_head(
        n = protein_label_number_S5BC
      )
    
    nmd3_target <- result_table %>%
      dplyr::filter(
        toupper(gene) == "NMD3"
      )
    
    dplyr::bind_rows(
      top_increased,
      top_decreased,
      nmd3_target
    ) %>%
      dplyr::distinct(
        gene,
        .keep_all = TRUE
      )
  }
  
  protein_labels_nucleus_S5C <-
    select_protein_labels_S5C(
      protein_nucleus_S5C
    )
  
  protein_labels_cytoplasm_S5C <-
    select_protein_labels_S5C(
      protein_cytoplasm_S5C
    )
  
  cat("Proteomics nucleus labels:\n")
  print(protein_labels_nucleus_S5C$gene)
  
  cat("Proteomics cytoplasm labels:\n")
  print(protein_labels_cytoplasm_S5C$gene)
  
  # ------------------------------------------------------------
  # 8. Protein volcano plotting function
  # ------------------------------------------------------------
  
  make_protein_volcano_S5C <- function(
    result_table,
    label_table,
    panel_title,
    show_x_title = TRUE,
    show_legend = TRUE
  ) {
    
    ggplot(
      result_table,
      aes(
        x = log2FC,
        y = neg_log10_pvalue
      )
    ) +
      geom_hline(
        yintercept = -log10(
          protein_p_threshold_S5BC
        ),
        linetype = "dotted",
        color = "gray50",
        linewidth = 0.4
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dotted",
        color = "gray50",
        linewidth = 0.4
      ) +
      geom_vline(
        xintercept = c(-1, 1),
        linetype = "dotted",
        color = "gray60",
        linewidth = 0.4
      ) +
      geom_point(
        data = result_table %>%
          dplyr::filter(
            Direction == "NS"
          ),
        color = "gray80",
        size = 1.2,
        alpha = 0.65
      ) +
      geom_point(
        data = result_table %>%
          dplyr::filter(
            Direction != "NS"
          ),
        aes(
          color = Direction
        ),
        size = 1.4,
        alpha = 0.85
      ) +
      ggrepel::geom_text_repel(
        data = label_table,
        aes(
          label = gene
        ),
        color = "black",
        size = 2.8,
        box.padding = 0.30,
        point.padding = 0.15,
        segment.color = "gray40",
        segment.linewidth = 0.20,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = volcano_colors_S5BC,
        breaks = c(
          "siNMD3 increased",
          "siNMD3 decreased",
          "NS"
        ),
        drop = FALSE,
        name = NULL
      ) +
      coord_cartesian(
        ylim = c(0, 10),
        clip = "off"
      ) +
      labs(
        title = panel_title,
        x = if (show_x_title) {
          expression(log[2]~"Fold change")
        } else {
          NULL
        },
        y = expression(
          -log[10]~italic(p)*"-value"
        )
      ) +
      theme_classic(
        base_size = 12
      ) +
      theme(
        axis.line = element_line(
          color = "black",
          linewidth = 0.7
        ),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.6
        ),
        axis.text = element_text(
          color = "black",
          size = 10
        ),
        axis.title = element_text(
          color = "black",
          size = 12,
          face = "bold"
        ),
        plot.title = element_text(
          color = "black",
          size = 14,
          face = "bold",
          hjust = 0.05
        ),
        legend.position = if (
          show_legend
        ) {
          "top"
        } else {
          "none"
        },
        legend.text = element_text(
          size = 9
        ),
        plot.margin = margin(
          8,
          8,
          8,
          8
        )
      )
  }
  
  p_Protein_nucleus_S5C <- make_protein_volcano_S5C(
    result_table = protein_nucleus_S5C,
    label_table = protein_labels_nucleus_S5C,
    panel_title = "Nucleus",
    show_x_title = FALSE,
    show_legend = TRUE
  )
  
  p_Protein_cytoplasm_S5C <- make_protein_volcano_S5C(
    result_table = protein_cytoplasm_S5C,
    label_table = protein_labels_cytoplasm_S5C,
    panel_title = "Cytoplasm",
    show_x_title = TRUE,
    show_legend = FALSE
  )
  
  p_S5C <- (
    p_Protein_nucleus_S5C /
      p_Protein_cytoplasm_S5C
  ) +
    patchwork::plot_annotation(
      title = "Proteomics",
      theme = theme(
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        )
      )
    )
  
  
  p_S5BC <- p_S5B | p_S5C
  
  
  # ============================================================
  # 10. Save figures
  # ============================================================
  
  
  
  nmd3_save_pdf(
    filename = "FigS6C_siNMD3_Proteomics_Volcano_Nucleus_Cytoplasm.pdf",
    plot = p_S5C,
    width = 6.5,
    height = 9,
    units = "in",
    device = cairo_pdf
  )
  
  # ============================================================
  # 11. Save source data
  # ============================================================
  
  rna_results_S5B <- dplyr::bind_rows(
    rna_nucleus_S5B,
    rna_cytoplasm_S5B
  )
  
  protein_results_S5C <- dplyr::bind_rows(
    protein_nucleus_S5C,
    protein_cytoplasm_S5C
  )
  
  rna_labels_S5B <- dplyr::bind_rows(
    rna_labels_nucleus_S5B,
    rna_labels_cytoplasm_S5B
  )
  
  protein_labels_S5C <- dplyr::bind_rows(
    protein_labels_nucleus_S5C,
    protein_labels_cytoplasm_S5C
  )
  
  write.csv(
    rna_results_S5B,
    "FigureS5B_RNA_volcano_source_data.csv",
    row.names = FALSE
  )
  
  write.csv(
    protein_results_S5C,
    "FigureS5C_Proteomics_volcano_source_data.csv",
    row.names = FALSE
  )
  
  write.csv(
    rna_labels_S5B,
    "FigureS5B_RNA_labeled_genes.csv",
    row.names = FALSE
  )
  
  write.csv(
    protein_labels_S5C,
    "FigureS5C_Proteomics_labeled_genes.csv",
    row.names = FALSE
  )
  
}))

# =============================================================================
# FigS6D
# =============================================================================
nmd3_step("FigS6D", c("Prepare_Fig3B_S6"), quote({
  # 6. Build all eight Figure S5D panels
  # ============================================================
  
  p_rna_nmd3_promoted_S5D <- make_venn_S5D(
    sets_S5D$RNA_NMD3_promoted_Nucleus,
    sets_S5D$RNA_NMD3_promoted_Cytoplasm,
    "RNA-seq: NMD3-promoted",
    expression(
      abs(log[2] * FC) > 1 ~
        "and DESeq2 BH-FDR < 0.05"
    )
  )
  
  p_rna_nmd3_suppressed_S5D <- make_venn_S5D(
    sets_S5D$RNA_NMD3_suppressed_Nucleus,
    sets_S5D$RNA_NMD3_suppressed_Cytoplasm,
    "RNA-seq: NMD3-suppressed",
    expression(
      abs(log[2] * FC) > 1 ~
        "and DESeq2 BH-FDR < 0.05"
    )
  )
  
  p_rna_xpo1_promoted_S5D <- make_venn_S5D(
    sets_S5D$RNA_XPO1_promoted_Nucleus,
    sets_S5D$RNA_XPO1_promoted_Cytoplasm,
    "RNA-seq: XPO1-promoted",
    expression(
      abs(log[2] * FC) > 1 ~
        "and DESeq2 BH-FDR < 0.05"
    )
  )
  
  p_rna_xpo1_suppressed_S5D <- make_venn_S5D(
    sets_S5D$RNA_XPO1_suppressed_Nucleus,
    sets_S5D$RNA_XPO1_suppressed_Cytoplasm,
    "RNA-seq: XPO1-suppressed",
    expression(
      abs(log[2] * FC) > 1 ~
        "and DESeq2 BH-FDR < 0.05"
    )
  )
  
  p_prot_nmd3_promoted_S5D <- make_venn_S5D(
    sets_S5D$Protein_NMD3_promoted_Nucleus,
    sets_S5D$Protein_NMD3_promoted_Cytoplasm,
    "Proteomics: NMD3-promoted",
    expression(
      log[2] * FC < 0 ~
        "and Tukey p < 0.05"
    )
  )
  
  p_prot_nmd3_suppressed_S5D <- make_venn_S5D(
    sets_S5D$Protein_NMD3_suppressed_Nucleus,
    sets_S5D$Protein_NMD3_suppressed_Cytoplasm,
    "Proteomics: NMD3-suppressed",
    expression(
      log[2] * FC > 0 ~
        "and Tukey p < 0.05"
    )
  )
  
  p_prot_xpo1_promoted_S5D <- make_venn_S5D(
    sets_S5D$Protein_XPO1_promoted_Nucleus,
    sets_S5D$Protein_XPO1_promoted_Cytoplasm,
    "Proteomics: XPO1-promoted",
    expression(
      log[2] * FC < 0 ~
        "and Tukey p < 0.05"
    )
  )
  
  p_prot_xpo1_suppressed_S5D <- make_venn_S5D(
    sets_S5D$Protein_XPO1_suppressed_Nucleus,
    sets_S5D$Protein_XPO1_suppressed_Cytoplasm,
    "Proteomics: XPO1-suppressed",
    expression(
      log[2] * FC > 0 ~
        "and Tukey p < 0.05"
    )
  )
  
  
  # ============================================================
  # 7. Assemble Figure S5D
  # ============================================================
  
  p_S5D_final <- (
    p_rna_nmd3_promoted_S5D |
      p_rna_nmd3_suppressed_S5D
  ) /
    (
      p_rna_xpo1_promoted_S5D |
        p_rna_xpo1_suppressed_S5D
    ) /
    (
      p_prot_nmd3_promoted_S5D |
        p_prot_nmd3_suppressed_S5D
    ) /
    (
      p_prot_xpo1_promoted_S5D |
        p_prot_xpo1_suppressed_S5D
    ) +
    patchwork::plot_annotation(
      title = "PC3",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 18
        )
      )
    )
  
  # ============================================================
  # 8. Save Figure S5D
  # ============================================================
  
  nmd3_save_pdf(
    filename = "FigS6D_final.pdf",
    plot = p_S5D_final,
    width = 10,
    height = 15,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  # Optional PNG export intentionally omitted: the PPT workflow requires PDFs only.
  
  
  # ============================================================
  # 9. Print and save overlap counts
  # ============================================================
  
  summarise_pair_S5D <- function(
    analysis_name,
    nucleus_set,
    cytoplasm_set
  ) {
    
    data.frame(
      Analysis = analysis_name,
      Nucleus_only = length(
        setdiff(
          nucleus_set,
          cytoplasm_set
        )
      ),
      Shared = length(
        intersect(
          nucleus_set,
          cytoplasm_set
        )
      ),
      Cytoplasm_only = length(
        setdiff(
          cytoplasm_set,
          nucleus_set
        )
      ),
      stringsAsFactors = FALSE
    )
  }
  
  overlap_summary_S5D <- dplyr::bind_rows(
    
    summarise_pair_S5D(
      "RNA_NMD3_promoted",
      sets_S5D$RNA_NMD3_promoted_Nucleus,
      sets_S5D$RNA_NMD3_promoted_Cytoplasm
    ),
    
    summarise_pair_S5D(
      "RNA_NMD3_suppressed",
      sets_S5D$RNA_NMD3_suppressed_Nucleus,
      sets_S5D$RNA_NMD3_suppressed_Cytoplasm
    ),
    
    summarise_pair_S5D(
      "RNA_XPO1_promoted",
      sets_S5D$RNA_XPO1_promoted_Nucleus,
      sets_S5D$RNA_XPO1_promoted_Cytoplasm
    ),
    
    summarise_pair_S5D(
      "RNA_XPO1_suppressed",
      sets_S5D$RNA_XPO1_suppressed_Nucleus,
      sets_S5D$RNA_XPO1_suppressed_Cytoplasm
    ),
    
    summarise_pair_S5D(
      "Protein_NMD3_promoted",
      sets_S5D$Protein_NMD3_promoted_Nucleus,
      sets_S5D$Protein_NMD3_promoted_Cytoplasm
    ),
    
    summarise_pair_S5D(
      "Protein_NMD3_suppressed",
      sets_S5D$Protein_NMD3_suppressed_Nucleus,
      sets_S5D$Protein_NMD3_suppressed_Cytoplasm
    ),
    
    summarise_pair_S5D(
      "Protein_XPO1_promoted",
      sets_S5D$Protein_XPO1_promoted_Nucleus,
      sets_S5D$Protein_XPO1_promoted_Cytoplasm
    ),
    
    summarise_pair_S5D(
      "Protein_XPO1_suppressed",
      sets_S5D$Protein_XPO1_suppressed_Nucleus,
      sets_S5D$Protein_XPO1_suppressed_Cytoplasm
    )
  )
  
  print(
    overlap_summary_S5D
  )
  
  utils::write.csv(
    overlap_summary_S5D,
    "FigureS5D_overlap_counts.csv",
    row.names = FALSE
  )
  
  set_membership_S5D <- utils::stack(
    sets_S5D
  )
  
  colnames(
    set_membership_S5D
  ) <- c(
    "Gene",
    "Set"
  )
  
  utils::write.csv(
    set_membership_S5D,
    "FigureS5D_gene_sets.csv",
    row.names = FALSE
  )
  
}))

# =============================================================================
# FigS6E
# =============================================================================
nmd3_step("FigS6E", c("Fig3A","Fig3C"), quote({
  
  go_bp_sets <- nmd3_cached_msig(species="Homo sapiens",collection="C5",subcollection="GO:BP")
  multi_kpt_noRPL <- multi_kpt[!grepl("^RPL",toupper(multi_kpt$gene)),,drop=FALSE]
  multi_si_noRPL <- multi_si[!grepl("^RPL",toupper(multi_si$gene)),,drop=FALSE]
  # =============================================================================
  # PART 6: RNA–PROTEIN NoRB VS NoRPL CENTRAL 95% + SPEARMAN
  # =============================================================================
  
  
  
  # =============================================================================
  # NO-RB GENE-SET SENSITIVITY: current Supplementary Fig. S6
  # =============================================================================
  # Important: NoRB = GO:0042254; NOT the same as NoRPL or NoRPL/RPS.
  # Reuses MSigDB and core gene-level fold changes, WITHOUT rerunning DESeq2/limma.
  # DO NOT substitute the public-proteomics RB score set: it intentionally excludes NMD3.
  .nmd3_rb_rows <- go_bp_sets$gs_name == "GOBP_RIBOSOME_BIOGENESIS"
  if ("gs_exact_source" %in% names(go_bp_sets))
    .nmd3_rb_rows <- .nmd3_rb_rows | go_bp_sets$gs_exact_source == "GO:0042254"
  .nmd3_rb_rows[is.na(.nmd3_rb_rows)] <- FALSE
  .nmd3_full_rb_genes <- unique(toupper(as.character(
    go_bp_sets$gene_symbol[.nmd3_rb_rows])))
  .nmd3_full_rb_genes <- .nmd3_full_rb_genes[!is.na(.nmd3_full_rb_genes) &
                                               nzchar(.nmd3_full_rb_genes)]
  stopifnot(length(.nmd3_full_rb_genes) > 20)
  stopifnot(exists("multi_kpt_noRPL"),
            exists("multi_si_noRPL"), exists("multi_kpt"), exists("multi_si"))
  .nmd3_rb_kpt <- multi_kpt[!(toupper(multi_kpt$gene) %in% .nmd3_full_rb_genes), , drop=FALSE]
  .nmd3_rb_si  <- multi_si[!(toupper(multi_si$gene) %in% .nmd3_full_rb_genes), , drop=FALSE]
  stopifnot(nrow(.nmd3_rb_kpt) > 100, nrow(.nmd3_rb_si) > 100)
  .nmd3_rb_dir <- file.path(.nmd3_run_root, "extra_sensitivity", "Central95_NoRB")
  dir.create(.nmd3_rb_dir, showWarnings=FALSE, recursive=TRUE)
  .nmd3_scatter4 <- function(xpo, sirna, with_text=TRUE) {
    top1 <- plot_center95(xpo, get_95_bounds(xpo), "red",
                          "RNA log2FC (KPT-330 vs DMSO)","Protein log2FC (KPT-330 vs DMSO)",
                          x_lim=c(-3,3), y_lim=c(-3,3))
    top2 <- plot_center95(sirna, get_95_bounds(sirna), "#2448E8",
                          "RNA log2FC (siNMD3 vs siCTRL)","Protein log2FC (siNMD3 vs siCTRL)",
                          x_lim=c(-1.8,1.8), y_lim=c(-3,3))
    plot_bottom <- function(df, color, xlab, ylab, top_label=FALSE) {
      p <- ggplot2::ggplot(df, ggplot2::aes(log2FC_RNA, log2FC_Prot)) +
        ggplot2::geom_point(size=1.5, alpha=0.5, color="grey30") +
        ggplot2::geom_smooth(method="lm", color=color, fill="grey70",
                             linewidth=0.8, se=TRUE, alpha=0.3) +
        ggplot2::labs(x=xlab,y=ylab) + ggplot2::theme_classic(base_size=14) +
        ggplot2::theme(axis.line=ggplot2::element_blank(),
                       panel.border=ggplot2::element_rect(color="black",fill=NA,linewidth=0.8),
                       axis.text=ggplot2::element_text(color="black"))
      if (with_text) {
        ct <- suppressWarnings(stats::cor.test(df$log2FC_RNA,df$log2FC_Prot,
                                               method="spearman",exact=FALSE))
        x <- min(df$log2FC_RNA) + .05*diff(range(df$log2FC_RNA))
        y <- if(top_label) max(df$log2FC_Prot)-.05*diff(range(df$log2FC_Prot))
        else min(df$log2FC_Prot)+.08*diff(range(df$log2FC_Prot))
        p <- p + ggplot2::annotate("text",x=x,y=y,
                                   label=sprintf("rho = %.3f, p = %.2e\nn = %d",
                                                 unname(ct$estimate),ct$p.value,nrow(df)),
                                   hjust=0,vjust=if(top_label) 1 else 0,size=5)
      }
      p
    }
    bot1 <- plot_bottom(xpo,"red","RNA log2FC (KPT-330 vs DMSO)",
                        "Protein log2FC (KPT-330 vs DMSO)")
    bot2 <- plot_bottom(sirna,"#2456FF","RNA log2FC (siNMD3 vs siCTRL)",
                        "Protein log2FC (siNMD3 vs siCTRL)",TRUE)
    bot1 | bot2
  }
  .nmd3_rb_tests <- function(df,label) {
    ct <- suppressWarnings(stats::cor.test(df$log2FC_RNA,df$log2FC_Prot,
                                           method="spearman",exact=FALSE))
    out <- data.frame(Subset=label, N=nrow(df), Spearman_rho=unname(ct$estimate),
                      P_value=ct$p.value, t(get_95_bounds(df)), check.names=FALSE)
    return(out)
  }
  .nmd3_rb_summary <- dplyr::bind_rows(
    .nmd3_rb_tests(.nmd3_rb_kpt,"KPT-330_NoRB"),
    .nmd3_rb_tests(.nmd3_rb_si,"siNMD3_NoRB"),
    .nmd3_rb_tests(multi_kpt_noRPL,"KPT-330_NoRPL"),
    .nmd3_rb_tests(multi_si_noRPL,"siNMD3_NoRPL"))
  utils::write.csv(.nmd3_rb_summary,file.path(.nmd3_rb_dir,
                                              "NoRB_NoRPL_Central95_Spearman_Statistics.csv"),row.names=FALSE)
  for (.nmd3_option in TRUE) {
    .nmd3_suffix <- if(.nmd3_option) "WithStatistics" else "NoText"
    .nmd3_p <- .nmd3_scatter4(.nmd3_rb_kpt,.nmd3_rb_si,.nmd3_option)
    nmd3_save_pdf(file.path(.nmd3_rb_dir,
                            "FigS6E_top_NoRB_Spearman.pdf"),
                  .nmd3_p,width=12,height=5,device=grDevices::cairo_pdf)
  }
  # NoRPL (already computed): no-label visualization only, NO recomputations.
  .nmd3_p <- .nmd3_scatter4(multi_kpt_noRPL,multi_si_noRPL,TRUE)
  nmd3_save_pdf(file.path(.nmd3_rb_dir,"FigS6E_bottom_NoRPL_Spearman.pdf"),
                .nmd3_p,width=12,height=5,device=grDevices::cairo_pdf)
  
  
  
  
}))

# =============================================================================
# FigS7A
# =============================================================================
nmd3_step("FigS7A", c("Prepare_KPT_RNA_Protein"), quote({
  # 2. Matched RNA-protein universe and original signed ranking
  # -----------------------------------------------------------------------------
  
  tiny_p <- 1e-300
  
  prepare_rrho <- function(data) {
    data %>%
      dplyr::mutate(
        gene = toupper(gene),
        log2FC = as.numeric(log2FC),
        padj = as.numeric(padj),
        padj = ifelse(
          is.na(padj),
          1,
          pmin(pmax(padj, tiny_p), 1)
        )
      ) %>%
      dplyr::filter(!is.na(gene), is.finite(log2FC)) %>%
      dplyr::group_by(gene) %>%
      dplyr::slice_max(abs(log2FC), n = 1L, with_ties = FALSE) %>%
      dplyr::ungroup()
  }
  
  rna_rrho <- res_rna_df %>%
    dplyr::transmute(gene, log2FC = log2FC_RNA, padj = padj_RNA) %>%
    prepare_rrho()
  prot_rrho <- res_prot_df %>%
    dplyr::transmute(gene, log2FC = log2FC_Prot, padj = padj_Prot) %>%
    prepare_rrho()
  
  universe <- base::intersect(rna_rrho$gene, prot_rrho$gene)
  rna_rrho <- rna_rrho %>% dplyr::filter(gene %in% universe)
  prot_rrho <- prot_rrho %>% dplyr::filter(gene %in% universe)
  
  rank_genes <- function(data) {
    data %>%
      dplyr::mutate(stat = log2FC * -log10(padj)) %>%
      dplyr::arrange(dplyr::desc(stat))
  }
  
  rna_ranked <- rank_genes(rna_rrho)
  prot_ranked <- rank_genes(prot_rrho)
  
  rna_up <- rna_ranked %>%
    dplyr::filter(stat > 0) %>%
    dplyr::arrange(dplyr::desc(stat)) %>%
    dplyr::pull(gene)
  rna_down <- rna_ranked %>%
    dplyr::filter(stat < 0) %>%
    dplyr::arrange(stat) %>%
    dplyr::pull(gene)
  prot_up <- prot_ranked %>%
    dplyr::filter(stat > 0) %>%
    dplyr::arrange(dplyr::desc(stat)) %>%
    dplyr::pull(gene)
  prot_down <- prot_ranked %>%
    dplyr::filter(stat < 0) %>%
    dplyr::arrange(stat) %>%
    dplyr::pull(gene)
  
  N <- length(universe)
  grid_points <- 60L
  max_rank_frac <- 0.50
  cap_log10p <- 50
  
  stopifnot(
    N > 100L,
    all(c(length(rna_up), length(rna_down), length(prot_up), length(prot_down)) > 0L)
  )
  
  # Original RRHO grid retained exactly.
  rrho_pair <- function(listA, listB, panel_label) {
    maxK <- max(50L, floor(N * max_rank_frac))
    Ks <- unique(as.integer(floor(
      seq(50, maxK, length.out = grid_points)
    )))
    output <- vector("list", length(Ks)^2)
    counter <- 1L
    
    for (i in seq_along(Ks)) {
      A <- head(listA, Ks[i])
      for (j in seq_along(Ks)) {
        B <- head(listB, Ks[j])
        overlap <- length(base::intersect(A, B))
        kA <- length(A)
        kB <- length(B)
        p_value <- stats::phyper(
          overlap - 1L,
          kA,
          N - kA,
          kB,
          lower.tail = FALSE
        )
        output[[counter]] <- data.frame(
          i_rank = Ks[i],
          j_rank = Ks[j],
          log10p = -log10(p_value + tiny_p),
          overlap = overlap,
          panel = panel_label
        )
        counter <- counter + 1L
      }
    }
    dplyr::bind_rows(output)
  }
  
  panel_levels <- c("Down-Down", "Down-Up", "Up-Down", "Up-Up")
  rrho_all <- dplyr::bind_rows(
    rrho_pair(rna_down, prot_down, "Down-Down"),
    rrho_pair(rna_down, prot_up, "Down-Up"),
    rrho_pair(rna_up, prot_down, "Up-Down"),
    rrho_pair(rna_up, prot_up, "Up-Up")
  ) %>%
    dplyr::mutate(
      panel = factor(panel, levels = panel_levels),
      log10p_cap = pmin(log10p, cap_log10p)
    )
  
  H1 <- rrho_all %>%
    dplyr::filter(panel == "Up-Up") %>%
    dplyr::slice_max(log10p, n = 1L, with_ties = FALSE) %>%
    dplyr::mutate(H_label = "H1")
  H2 <- rrho_all %>%
    dplyr::filter(panel == "Down-Down") %>%
    dplyr::slice_max(log10p, n = 1L, with_ties = FALSE) %>%
    dplyr::mutate(H_label = "H2")
  hotspots <- dplyr::bind_rows(H1, H2)
  print(hotspots)
  
  # Diagnostic only: no result is modified. A failure indicates that H1/H2 was
  # selected from a rank label beyond the corresponding real directional list.
  stopifnot(
    H1$i_rank <= length(rna_up),
    H1$j_rank <= length(prot_up),
    H2$i_rank <= length(rna_down),
    H2$j_rank <= length(prot_down)
  )
  
  # -----------------------------------------------------------------------------
  # 3. RPL/RPS positions and Figure 2G
  # -----------------------------------------------------------------------------
  
  is_RPL_strict <- function(x) grepl("^RPL(\\d+[A-Z]?|P[0-9]+)$", x)
  is_RPS_strict <- function(x) grepl("^RPS(\\d+[A-Z]?|A|4Y1)$", x)
  
  matched_data <- dplyr::inner_join(
    rna_rrho %>% dplyr::rename(log2FC_RNA = log2FC),
    prot_rrho %>% dplyr::rename(log2FC_Prot = log2FC),
    by = "gene"
  )
  
  family_data <- matched_data %>%
    dplyr::mutate(
      Family = dplyr::case_when(
        is_RPL_strict(gene) ~ "RPL",
        is_RPS_strict(gene) ~ "RPS",
        TRUE ~ NA_character_
      ),
      panel = dplyr::case_when(
        log2FC_RNA < 0 & log2FC_Prot < 0 ~ "Down-Down",
        log2FC_RNA < 0 & log2FC_Prot > 0 ~ "Down-Up",
        log2FC_RNA > 0 & log2FC_Prot < 0 ~ "Up-Down",
        log2FC_RNA > 0 & log2FC_Prot > 0 ~ "Up-Up"
      )
    ) %>%
    dplyr::filter(!is.na(Family), !is.na(panel))
  
  get_family_rank <- function(gene, panel) {
    if (panel == "Down-Down") {
      c(match(gene, rna_down), match(gene, prot_down))
    } else if (panel == "Down-Up") {
      c(match(gene, rna_down), match(gene, prot_up))
    } else if (panel == "Up-Down") {
      c(match(gene, rna_up), match(gene, prot_down))
    } else {
      c(match(gene, rna_up), match(gene, prot_up))
    }
  }
  
  family_coords <- lapply(seq_len(nrow(family_data)), function(i) {
    xy <- get_family_rank(family_data$gene[i], family_data$panel[i])
    data.frame(
      gene = family_data$gene[i],
      Family = family_data$Family[i],
      panel = family_data$panel[i],
      i_rank = xy[1],
      j_rank = xy[2]
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(
      is.finite(i_rank),
      is.finite(j_rank),
      i_rank <= floor(N * max_rank_frac),
      j_rank <= floor(N * max_rank_frac)
    ) %>%
    dplyr::mutate(panel = factor(panel, levels = panel_levels))
  
  contour_breaks <- stats::quantile(
    rrho_all$overlap,
    c(0.50, 0.75, 0.90, 0.95, 0.99),
    na.rm = TRUE
  )
  
  fig2g <- ggplot2::ggplot(
    rrho_all,
    ggplot2::aes(x = i_rank, y = j_rank, fill = log10p_cap)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_contour(
      ggplot2::aes(z = overlap),
      breaks = contour_breaks,
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      data = family_coords,
      ggplot2::aes(x = i_rank, y = j_rank, shape = Family),
      inherit.aes = FALSE,
      color = "black",
      fill = "white",
      size = 1.8,
      stroke = 0.5
    ) +
    ggplot2::geom_label(
      data = hotspots,
      ggplot2::aes(x = i_rank, y = j_rank, label = H_label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.5,
      fill = "white"
    ) +
    ggplot2::scale_shape_manual(
      values = c("RPL" = 24, "RPS" = 22),
      name = "Family"
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = c(0, cap_log10p),
      name = expression(-log[10] ~ P)
    ) +
    ggplot2::facet_wrap(~ panel, ncol = 2) +
    ggplot2::labs(x = "RNA rank", y = "Protein rank") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 1
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  
  nmd3_save_pdf(
    filename = file.path(output_dir, "FigS7A_RNA_Protein.pdf"),
    plot = fig2g,
    width = 7.2,
    height = 5.2,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  utils::write.csv(
    matched_data,
    file.path(output_dir, "Fig2G_RNA_Protein_matched_source_data.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    rrho_all,
    file.path(output_dir, "Fig2G_RRHO_grid_statistics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    hotspots,
    file.path(output_dir, "Fig2G_RRHO_hotspots.csv"),
    row.names = FALSE
  )
  
}))

# =============================================================================
# FigS7B
# =============================================================================
nmd3_step("FigS7B", c("FigS7A"), quote({
  # FigS7B — analysis below
  # =============================================================================
  
  go_palette <- c(
    "#4575B4", "#74ADD1", "#ABD9E9", "#FEE090",
    "#FDAE61", "#F46D43", "#D73027", "#A50026"
  )
  
  topk <- function(x, k) head(x, as.integer(k))
  genes_H1 <- unique(base::intersect(
    topk(rna_up, H1$i_rank),
    topk(prot_up, H1$j_rank)
  ))
  genes_H2 <- unique(base::intersect(
    topk(rna_down, H2$i_rank),
    topk(prot_down, H2$j_rank)
  ))
  
  stopifnot(
    length(genes_H1) == H1$overlap,
    length(genes_H2) == H2$overlap,
    length(genes_H1) > 0L,
    length(genes_H2) > 0L
  )
  
  cat("H1 genes:", length(genes_H1), "\n")
  cat("H2 genes:", length(genes_H2), "\n")
  
  # The original enrichGO settings, including its default background, are kept.
  run_hotspot_GO <- function(genes, hotspot_name) {
    ego <- clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    )
    ego_df <- as.data.frame(ego)
    stopifnot(nrow(ego_df) > 0L)
    
    ego_df <- ego_df %>%
      dplyr::mutate(
        GeneRatio_numeric =
          as.numeric(sub("/.*$", "", GeneRatio)) /
          as.numeric(sub("^.*/", "", GeneRatio)),
        neglog10FDR = -log10(pmax(p.adjust, .Machine$double.xmin))
      )
    
    top10 <- ego_df %>%
      dplyr::arrange(
        p.adjust,
        pvalue,
        dplyr::desc(GeneRatio_numeric),
        Description
      ) %>%
      dplyr::slice_head(n = 10L) %>%
      dplyr::arrange(GeneRatio_numeric) %>%
      dplyr::mutate(
        Description = factor(Description, levels = Description)
      )
    
    plot <- ggplot2::ggplot(
      top10,
      ggplot2::aes(
        x = GeneRatio_numeric,
        y = Description,
        fill = neglog10FDR
      )
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::scale_fill_gradientn(
        colours = go_palette,
        name = expression(-log[10] ~ FDR)
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::percent_format(accuracy = 0.1),
        expand = ggplot2::expansion(mult = c(0.01, 0.05))
      ) +
      ggplot2::labs(
        title = paste0(hotspot_name, " GO Biological Process"),
        x = "Gene Ratio",
        y = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 13
        ),
        axis.text.y = ggplot2::element_text(
          size = 11, color = "black", face = "bold"
        ),
        axis.text.x = ggplot2::element_text(size = 11, color = "black"),
        axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        panel.border = ggplot2::element_rect(
          color = "black", fill = NA, linewidth = 1.1
        )
      )
    
    list(enrichment = ego, all_pathways = ego_df, top10 = top10, plot = plot)
  }
  
  H1_GO <- run_hotspot_GO(genes_H1, "H1")
  H2_GO <- run_hotspot_GO(genes_H2, "H2")
  fig2h <- H1_GO$plot + H2_GO$plot + patchwork::plot_layout(ncol = 2)
  
  nmd3_save_pdf(
    filename = file.path(output_dir, "FigS7B_H1_H2_GO_BP_Top10.pdf"),
    plot = fig2h,
    width = 12,
    height = 5.5,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  utils::write.csv(
    dplyr::bind_rows(
      H1_GO$all_pathways %>% dplyr::mutate(Hotspot = "H1", .before = 1L),
      H2_GO$all_pathways %>% dplyr::mutate(Hotspot = "H2", .before = 1L)
    ),
    file.path(output_dir, "Fig2H_GO_BP_all_results.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    dplyr::bind_rows(
      H1_GO$top10 %>%
        dplyr::mutate(Hotspot = "H1", Description = as.character(Description)),
      H2_GO$top10 %>%
        dplyr::mutate(Hotspot = "H2", Description = as.character(Description))
    ),
    file.path(output_dir, "Fig2H_GO_BP_plotted_results.csv"),
    row.names = FALSE
  )
  
}))

# =============================================================================
# FigS7C
# =============================================================================
nmd3_step("FigS7C", c("Fig3A"), quote({
  # FigS7C — analysis below
  # siNMD3 vs siCTRL cytoplasmic RNA–protein RRHO
  #
  # Same statistical logic as the XPO1 RRHO analysis:
  # signed statistic + right-tailed hypergeometric test
  #
  # Input: multi_si generated in Figure 4
  # RPL = triangle; RPS = square
  # No H1/H2
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
    library(readr)
    library(scales)
  })
  
  # ------------------------------------------------------------
  # 1. Parameters
  # ------------------------------------------------------------
  
  grid_points_S6D   <- 60
  max_rank_frac_S6D <- 0.50
  cap_log10p_S6D    <- 50
  tiny_p_S6D        <- 1e-300
  fc_eps_S6D        <- 0.01
  
  # ------------------------------------------------------------
  # 2. Use the frozen Figure 4 siNMD3 RNA–protein result
  # ------------------------------------------------------------
  
  stopifnot(
    exists("multi_si"),
    all(
      c(
        "gene",
        "log2FC_RNA",
        "padj_RNA",
        "log2FC_Prot",
        "padj_Prot"
      ) %in% colnames(multi_si)
    )
  )
  
  input_S6D <- multi_si %>%
    dplyr::transmute(
      gene = toupper(trimws(as.character(gene))),
      log2FC_RNA = suppressWarnings(
        as.numeric(log2FC_RNA)
      ),
      padj_RNA = suppressWarnings(
        as.numeric(padj_RNA)
      ),
      log2FC_Prot = suppressWarnings(
        as.numeric(log2FC_Prot)
      ),
      padj_Prot = suppressWarnings(
        as.numeric(padj_Prot)
      )
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      is.finite(log2FC_RNA),
      is.finite(log2FC_Prot),
      is.finite(padj_RNA),
      is.finite(padj_Prot)
    ) %>%
    dplyr::mutate(
      padj_RNA = ifelse(
        padj_RNA > 0,
        pmin(padj_RNA, 1),
        1
      ),
      padj_Prot = ifelse(
        padj_Prot > 0,
        pmin(padj_Prot, 1),
        1
      )
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(
      order_by = abs(log2FC_RNA) + abs(log2FC_Prot),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      RNA_stat =
        log2FC_RNA *
        (-log10(padj_RNA + tiny_p_S6D)),
      
      Protein_stat =
        log2FC_Prot *
        (-log10(padj_Prot + tiny_p_S6D))
    )
  
  N_S6D <- nrow(input_S6D)
  
  if (N_S6D < 200) {
    stop(
      "RNA–protein matched universe is too small: ",
      N_S6D
    )
  }
  
  cat(
    "Figure S6D RNA–protein universe:",
    N_S6D,
    "genes\n"
  )
  
  print(
    input_S6D %>%
      dplyr::summarise(
        RNA_up = sum(RNA_stat > 0),
        RNA_down = sum(RNA_stat < 0),
        RNA_zero = sum(RNA_stat == 0),
        Protein_up = sum(Protein_stat > 0),
        Protein_down = sum(Protein_stat < 0),
        Protein_zero = sum(Protein_stat == 0)
      )
  )
  
  # ------------------------------------------------------------
  # 3. Construct directional ranked lists
  # ------------------------------------------------------------
  
  # Up = higher in siNMD3
  # Down = lower in siNMD3
  
  rna_up_S6D <- input_S6D %>%
    dplyr::filter(
      RNA_stat > 0
    ) %>%
    dplyr::arrange(
      dplyr::desc(RNA_stat)
    ) %>%
    dplyr::pull(gene)
  
  rna_down_S6D <- input_S6D %>%
    dplyr::filter(
      RNA_stat < 0
    ) %>%
    dplyr::arrange(
      RNA_stat
    ) %>%
    dplyr::pull(gene)
  
  protein_up_S6D <- input_S6D %>%
    dplyr::filter(
      Protein_stat > 0
    ) %>%
    dplyr::arrange(
      dplyr::desc(Protein_stat)
    ) %>%
    dplyr::pull(gene)
  
  protein_down_S6D <- input_S6D %>%
    dplyr::filter(
      Protein_stat < 0
    ) %>%
    dplyr::arrange(
      Protein_stat
    ) %>%
    dplyr::pull(gene)
  
  cat("\nDirectional list sizes:\n")
  cat("RNA up:       ", length(rna_up_S6D), "\n")
  cat("RNA down:     ", length(rna_down_S6D), "\n")
  cat("Protein up:   ", length(protein_up_S6D), "\n")
  cat("Protein down: ", length(protein_down_S6D), "\n")
  
  # ------------------------------------------------------------
  # 4. RRHO function
  # Right-tailed hypergeometric test
  # ------------------------------------------------------------
  
  run_RRHO_S6D <- function(
    RNA_list,
    protein_list,
    panel_name
  ) {
    
    maxK <- max(
      50L,
      floor(N_S6D * max_rank_frac_S6D)
    )
    
    rank_grid <- unique(
      as.integer(
        floor(
          seq(
            from = 50,
            to = maxK,
            length.out = grid_points_S6D
          )
        )
      )
    )
    
    RNA_prefixes <- lapply(
      rank_grid,
      function(k) {
        RNA_list[
          seq_len(
            min(k, length(RNA_list))
          )
        ]
      }
    )
    
    protein_prefixes <- lapply(
      rank_grid,
      function(k) {
        protein_list[
          seq_len(
            min(k, length(protein_list))
          )
        ]
      }
    )
    
    P_matrix <- matrix(
      NA_real_,
      nrow = length(rank_grid),
      ncol = length(rank_grid)
    )
    
    overlap_matrix <- matrix(
      NA_real_,
      nrow = length(rank_grid),
      ncol = length(rank_grid)
    )
    
    for (i in seq_along(rank_grid)) {
      
      RNA_top <- RNA_prefixes[[i]]
      k_RNA <- length(RNA_top)
      
      for (j in seq_along(rank_grid)) {
        
        protein_top <- protein_prefixes[[j]]
        k_protein <- length(protein_top)
        
        overlap_n <- length(
          intersect(
            RNA_top,
            protein_top
          )
        )
        
        rrho_p <- stats::phyper(
          q = overlap_n - 1,
          m = k_RNA,
          n = N_S6D - k_RNA,
          k = k_protein,
          lower.tail = FALSE
        )
        
        P_matrix[i, j] <- -log10(
          rrho_p + tiny_p_S6D
        )
        
        overlap_matrix[i, j] <- overlap_n
      }
    }
    
    expand.grid(
      i = seq_along(rank_grid),
      j = seq_along(rank_grid)
    ) %>%
      dplyr::transmute(
        RNA_rank = rank_grid[i],
        Protein_rank = rank_grid[j],
        neg_log10_P = as.numeric(
          P_matrix[
            cbind(i, j)
          ]
        ),
        Overlap = as.numeric(
          overlap_matrix[
            cbind(i, j)
          ]
        ),
        Panel = panel_name
      ) %>%
      dplyr::mutate(
        neg_log10_P = pmin(
          neg_log10_P,
          cap_log10p_S6D
        )
      )
  }
  
  # ------------------------------------------------------------
  # 5. Calculate four RRHO quadrants
  # ------------------------------------------------------------
  
  rrho_DownDown_S6D <- run_RRHO_S6D(
    RNA_list = rna_down_S6D,
    protein_list = protein_down_S6D,
    panel_name = "Down-Down"
  )
  
  rrho_DownUp_S6D <- run_RRHO_S6D(
    RNA_list = rna_down_S6D,
    protein_list = protein_up_S6D,
    panel_name = "Down-Up"
  )
  
  rrho_UpDown_S6D <- run_RRHO_S6D(
    RNA_list = rna_up_S6D,
    protein_list = protein_down_S6D,
    panel_name = "Up-Down"
  )
  
  rrho_UpUp_S6D <- run_RRHO_S6D(
    RNA_list = rna_up_S6D,
    protein_list = protein_up_S6D,
    panel_name = "Up-Up"
  )
  
  rrho_grid_S6D <- dplyr::bind_rows(
    rrho_DownDown_S6D,
    rrho_DownUp_S6D,
    rrho_UpDown_S6D,
    rrho_UpUp_S6D
  ) %>%
    dplyr::mutate(
      Panel = factor(
        Panel,
        levels = c(
          "Down-Down",
          "Down-Up",
          "Up-Down",
          "Up-Up"
        )
      )
    )
  
  # ------------------------------------------------------------
  # 6. Automatic RPL/RPS identification
  # Same strict definition as the XPO1 RRHO
  # ------------------------------------------------------------
  
  is_RPL_S6D <- function(x) {
    grepl(
      "^RPL[0-9]+[A-Z]?$",
      toupper(x),
      perl = TRUE
    )
  }
  
  is_RPS_S6D <- function(x) {
    grepl(
      "^RPS[0-9]+[A-Z]?$",
      toupper(x),
      perl = TRUE
    )
  }
  
  # ------------------------------------------------------------
  # 7. Obtain exact RPL/RPS coordinates
  # ------------------------------------------------------------
  
  maxK_S6D <- floor(
    N_S6D * max_rank_frac_S6D
  )
  
  family_coordinates_S6D <- input_S6D %>%
    dplyr::mutate(
      Family = dplyr::case_when(
        is_RPL_S6D(gene) ~ "RPL",
        is_RPS_S6D(gene) ~ "RPS",
        TRUE ~ NA_character_
      ),
      
      Panel = dplyr::case_when(
        log2FC_RNA < -fc_eps_S6D &
          log2FC_Prot < -fc_eps_S6D ~
          "Down-Down",
        
        log2FC_RNA < -fc_eps_S6D &
          log2FC_Prot > fc_eps_S6D ~
          "Down-Up",
        
        log2FC_RNA > fc_eps_S6D &
          log2FC_Prot < -fc_eps_S6D ~
          "Up-Down",
        
        log2FC_RNA > fc_eps_S6D &
          log2FC_Prot > fc_eps_S6D ~
          "Up-Up",
        
        TRUE ~ NA_character_
      ),
      
      RNA_rank = dplyr::case_when(
        Panel %in% c(
          "Down-Down",
          "Down-Up"
        ) ~ match(
          gene,
          rna_down_S6D
        ),
        
        Panel %in% c(
          "Up-Down",
          "Up-Up"
        ) ~ match(
          gene,
          rna_up_S6D
        ),
        
        TRUE ~ NA_integer_
      ),
      
      Protein_rank = dplyr::case_when(
        Panel %in% c(
          "Down-Down",
          "Up-Down"
        ) ~ match(
          gene,
          protein_down_S6D
        ),
        
        Panel %in% c(
          "Down-Up",
          "Up-Up"
        ) ~ match(
          gene,
          protein_up_S6D
        ),
        
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::filter(
      !is.na(Family),
      !is.na(Panel),
      !is.na(RNA_rank),
      !is.na(Protein_rank),
      RNA_rank <= maxK_S6D,
      Protein_rank <= maxK_S6D
    ) %>%
    dplyr::mutate(
      Panel = factor(
        Panel,
        levels = c(
          "Down-Down",
          "Down-Up",
          "Up-Down",
          "Up-Up"
        )
      ),
      Family = factor(
        Family,
        levels = c(
          "RPL",
          "RPS"
        )
      )
    )
  
  family_marker_summary_S6D <- family_coordinates_S6D %>%
    dplyr::group_by(
      Panel,
      Family,
      .drop = FALSE
    ) %>%
    dplyr::summarise(
      n = dplyr::n(),
      .groups = "drop"
    )
  
  cat("\nRPL/RPS markers displayed:\n")
  print(family_marker_summary_S6D)
  
  # ------------------------------------------------------------
  # 8. Contour levels
  # Use overlap-count quantiles from concordant quadrants
  # ------------------------------------------------------------
  
  contour_values_S6D <- rrho_grid_S6D %>%
    dplyr::filter(
      Panel %in% c(
        "Down-Down",
        "Up-Up"
      ),
      Overlap > 0
    ) %>%
    dplyr::pull(Overlap)
  
  if (length(contour_values_S6D) >= 4) {
    
    contour_breaks_S6D <- stats::quantile(
      contour_values_S6D,
      probs = c(
        0.50,
        0.75,
        0.90,
        0.95,
        0.99
      ),
      na.rm = TRUE,
      names = FALSE
    )
    
    contour_breaks_S6D <- sort(
      unique(
        as.numeric(contour_breaks_S6D)
      )
    )
    
    contour_breaks_S6D <- contour_breaks_S6D[
      is.finite(contour_breaks_S6D) &
        contour_breaks_S6D > 0
    ]
    
  } else {
    
    contour_breaks_S6D <- numeric(0)
  }
  
  # ------------------------------------------------------------
  # 9. Plotting function
  # ------------------------------------------------------------
  
  make_panel_S6D <- function(
    panel_name,
    panel_title,
    x_title,
    y_title
  ) {
    
    panel_grid <- rrho_grid_S6D %>%
      dplyr::filter(
        as.character(Panel) == panel_name
      )
    
    panel_markers <- family_coordinates_S6D %>%
      dplyr::filter(
        as.character(Panel) == panel_name
      )
    
    p <- ggplot2::ggplot(
      panel_grid,
      ggplot2::aes(
        x = RNA_rank,
        y = Protein_rank
      )
    ) +
      ggplot2::geom_raster(
        ggplot2::aes(
          fill = neg_log10_P
        ),
        interpolate = FALSE
      )
    
    if (length(contour_breaks_S6D) > 0) {
      
      p <- p +
        ggplot2::geom_contour(
          ggplot2::aes(
            z = Overlap
          ),
          breaks = contour_breaks_S6D,
          color = "black",
          linewidth = 0.28
        )
    }
    
    p +
      ggplot2::geom_point(
        data = panel_markers,
        ggplot2::aes(
          x = RNA_rank,
          y = Protein_rank,
          shape = Family
        ),
        inherit.aes = FALSE,
        fill = "white",
        color = "black",
        size = 2.6,
        stroke = 0.6
      ) +
      ggplot2::scale_shape_manual(
        values = c(
          "RPL" = 24,
          "RPS" = 22
        ),
        limits = c(
          "RPL",
          "RPS"
        ),
        drop = FALSE,
        name = "Family"
      ) +
      ggplot2::scale_fill_viridis_c(
        option = "C",
        direction = 1,
        limits = c(
          0,
          cap_log10p_S6D
        ),
        oob = scales::squish,
        name = expression(-log[10]~P)
      ) +
      ggplot2::scale_x_continuous(
        limits = c(
          0,
          maxK_S6D
        ),
        labels = scales::comma,
        expand = ggplot2::expansion(
          mult = c(0, 0.01)
        )
      ) +
      ggplot2::scale_y_continuous(
        limits = c(
          0,
          maxK_S6D
        ),
        labels = scales::comma,
        expand = ggplot2::expansion(
          mult = c(0, 0.01)
        )
      ) +
      ggplot2::labs(
        title = panel_title,
        x = x_title,
        y = y_title
      ) +
      ggplot2::theme_classic(
        base_size = 12
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold.italic",
          size = 13,
          hjust = 0.5
        ),
        axis.title = ggplot2::element_text(
          face = "bold",
          size = 11,
          color = "black"
        ),
        axis.text = ggplot2::element_text(
          face = "bold",
          size = 10,
          color = "black"
        ),
        axis.line = ggplot2::element_line(
          color = "black",
          linewidth = 0.8
        ),
        axis.ticks = ggplot2::element_line(
          color = "black",
          linewidth = 0.7
        ),
        panel.grid = ggplot2::element_blank()
      )
  }
  
  # ------------------------------------------------------------
  # 10. Four panels
  # ------------------------------------------------------------
  
  p_DownDown_S6D <- make_panel_S6D(
    panel_name = "Down-Down",
    panel_title = "RNA Down - Protein Down",
    x_title = "RNA rank (most downregulated)",
    y_title = "Protein rank (most downregulated)"
  )
  
  p_DownUp_S6D <- make_panel_S6D(
    panel_name = "Down-Up",
    panel_title = "RNA Down - Protein Up",
    x_title = "RNA rank (most downregulated)",
    y_title = "Protein rank (most upregulated)"
  )
  
  p_UpDown_S6D <- make_panel_S6D(
    panel_name = "Up-Down",
    panel_title = "RNA Up - Protein Down",
    x_title = "RNA rank (most upregulated)",
    y_title = "Protein rank (most downregulated)"
  )
  
  p_UpUp_S6D <- make_panel_S6D(
    panel_name = "Up-Up",
    panel_title = "RNA Up - Protein Up",
    x_title = "RNA rank (most upregulated)",
    y_title = "Protein rank (most upregulated)"
  )
  
  p_S6D <- (
    p_DownDown_S6D |
      p_DownUp_S6D
  ) / (
    p_UpDown_S6D |
      p_UpUp_S6D
  ) +
    patchwork::plot_layout(
      guides = "collect"
    )
  
  p_S6D <- p_S6D &
    ggplot2::theme(
      legend.position = "right"
    )
  
  
  # ------------------------------------------------------------
  # 11. Save results
  # ------------------------------------------------------------
  
  nmd3_save_pdf(
    filename = "FigS7C_siNMD3_RRHO_RPL_RPS.pdf",
    plot = p_S6D,
    width = 9,
    height = 7.5,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  # Optional PNG export intentionally omitted: the PPT workflow requires PDFs only.
  
  readr::write_csv(
    rrho_grid_S6D,
    "FigureS6D_RRHO_grid.csv"
  )
  
  readr::write_csv(
    family_coordinates_S6D,
    "FigureS6D_RPL_RPS_coordinates.csv"
  )
  
  readr::write_csv(
    family_marker_summary_S6D,
    "FigureS6D_RPL_RPS_summary.csv"
  )
  
  cat(
    "\nFigure S6D completed.\n",
    "Input: multi_si from Figure 4\n",
    "Ranking: signed log2FC × -log10 adjusted P\n",
    "Test: right-tailed hypergeometric RRHO\n",
    "No H1/H2 calculated or displayed.\n",
    sep = ""
  )
  
  
  
  
  
}))

# =============================================================================
# FigS7D
# =============================================================================
nmd3_step("FigS7D", c("Fig2F"), quote({
  # ============================================================================
  # 3. siNMD3 vs siCtrl — nuclear and cytoplasmic protein heatmaps
  # ============================================================================
  
  rp_sort_sample_cols <- function(cols) {
    # Original suffix: C1N, C2N, C3N / N1N, N2N, N3N etc.
    idx <- suppressWarnings(as.numeric(sub(".*[CN]([123])[CN]$", "\\1", cols)))
    if (anyNA(idx)) {
      stop("Cannot parse replicate IDs from: ", paste(cols[is.na(idx)], collapse = ", "))
    }
    cols[order(idx, cols)]
  }
  
  rp_si_nuc_kd <- rp_sort_sample_cols(grep(
    "^Norm_siRNAxNMD3xKD_Nuc_", colnames(prot_raw), value = TRUE
  ))
  rp_si_nuc_ctrl <- rp_sort_sample_cols(grep(
    "^Norm_siRNAxCtrl_Nuc_", colnames(prot_raw), value = TRUE
  ))
  rp_si_cyto_kd <- rp_sort_sample_cols(grep(
    "^Norm_siRNAxNMD3xKD_Cyto_", colnames(prot_raw), value = TRUE
  ))
  rp_si_cyto_ctrl <- rp_sort_sample_cols(grep(
    "^Norm_siRNAxCtrl_Cyto_", colnames(prot_raw), value = TRUE
  ))
  
  stopifnot(
    length(rp_si_nuc_kd) == 3, length(rp_si_nuc_ctrl) == 3,
    length(rp_si_cyto_kd) == 3, length(rp_si_cyto_ctrl) == 3
  )
  
  rp_si_nuc_cols <- c(rp_si_nuc_kd, rp_si_nuc_ctrl)
  rp_si_cyto_cols <- c(rp_si_cyto_kd, rp_si_cyto_ctrl)
  rp_si_all_cols <- c(rp_si_nuc_cols, rp_si_cyto_cols)
  
  rp_si_all_protein <- rp_protein_matrix(prot_raw, rp_si_all_cols)
  rp_si_nuc <- rp_si_all_protein[, rp_si_nuc_cols, drop = FALSE]
  rp_si_cyto <- rp_si_all_protein[, rp_si_cyto_cols, drop = FALSE]
  
  # Original Figure S6E: transform each compartment separately, if median > 50.
  if (stats::median(rp_si_nuc, na.rm = TRUE) > 50) {
    rp_si_nuc <- log2(rp_si_nuc + 1)
  }
  if (stats::median(rp_si_cyto, na.rm = TRUE) > 50) {
    rp_si_cyto <- log2(rp_si_cyto + 1)
  }
  
  rp_si_candidates <- rownames(rp_si_all_protein)
  rp_si_candidates <- rp_si_candidates[
    rp_is_rpl(rp_si_candidates) | rp_is_rps(rp_si_candidates)
  ]
  
  # Same gene universe in Nucleus and Cytoplasm; all 6 values must be finite
  # with positive SD in EACH compartment.
  rp_si_keep <- rp_si_candidates[vapply(rp_si_candidates, function(g) {
    rp_valid(rp_si_nuc[g, ]) && rp_valid(rp_si_cyto[g, ])
  }, logical(1))]
  rp_si_genes <- rp_sort_numeric(rp_si_keep)
  
  stopifnot(length(rp_si_genes) > 1L,
            any(rp_is_rpl(rp_si_genes)), any(rp_is_rps(rp_si_genes)))
  
  rp_si_nuc_z <- rp_zscore(rp_si_nuc[rp_si_genes, , drop = FALSE])
  rp_si_cyto_z <- rp_zscore(rp_si_cyto[rp_si_genes, , drop = FALSE])
  
  # Original Figure S6E arrangement: 6 samples in rows, genes in columns.
  rp_si_h_nuc <- t(rp_si_nuc_z)
  rp_si_h_cyto <- t(rp_si_cyto_z)
  
  rp_si_sample_labels <- c(
    "siNMD3-1", "siNMD3-2", "siNMD3-3",
    "siCtrl-1", "siCtrl-2", "siCtrl-3"
  )
  rownames(rp_si_h_nuc) <- rp_si_sample_labels
  rownames(rp_si_h_cyto) <- rp_si_sample_labels
  stopifnot(identical(colnames(rp_si_h_nuc), colnames(rp_si_h_cyto)))
  
  rp_si_gap <- sum(rp_is_rpl(rp_si_genes))
  
  rp_draw_si <- function(mat, title) {
    pheatmap::pheatmap(
      mat = mat,
      color = rp_colors,
      breaks = rp_breaks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 10,
      fontsize_col = 7.5,
      angle_col = 90,
      border_color = NA,
      gaps_row = 3,
      gaps_col = rp_si_gap,
      legend = FALSE,
      annotation_legend = FALSE,
      main = title,
      silent = TRUE
    )
  }
  
  rp_si_p_nuc <- rp_draw_si(rp_si_h_nuc, "Nucleus")
  rp_si_p_cyto <- rp_draw_si(rp_si_h_cyto, "Cytoplasm")
  
  rp_si_shared_widths <- grid::unit.pmax(
    rp_si_p_nuc$gtable$widths,
    rp_si_p_cyto$gtable$widths
  )
  rp_si_p_nuc$gtable$widths <- rp_si_shared_widths
  rp_si_p_cyto$gtable$widths <- rp_si_shared_widths
  
  rp_si_plot <- gridExtra::arrangeGrob(
    rp_si_p_nuc$gtable, rp_si_p_cyto$gtable,
    ncol = 1, heights = c(1, 1),
    top = grid::textGrob(
      "RPL/RPS proteomics: siNMD3 vs siCtrl",
      gp = grid::gpar(fontsize = 15, fontface = "bold")
    )
  )
  
  grid::grid.newpage()
  grid::grid.draw(rp_si_plot)
  
  rp_si_pdf <- file.path(rp_out, "FigS7D_siNMD3_RPL_RPS_Heatmaps.pdf")
  nmd3_cairo_pdf(
    rp_si_pdf,
    width = max(12, length(rp_si_genes) * 0.18),
    height = 9
  )
  grid::grid.newpage()
  grid::grid.draw(rp_si_plot)
  nmd3_dev_off()
  
  rp_save_genelist(rp_si_genes, "siNMD3_RPL_RPS_GeneOrder.csv")
  rp_save_matrix(rp_si_h_nuc, "siNMD3_Nucleus_RowZ.csv")
  rp_save_matrix(rp_si_h_cyto, "siNMD3_Cytoplasm_RowZ.csv")
  
  cat("\nsiNMD3: ", sum(rp_is_rpl(rp_si_genes)), " RPL; ",
      sum(rp_is_rps(rp_si_genes)), " RPS\n", sep = "")
  cat("siNMD3 excluded after QC: ",
      paste(setdiff(rp_si_candidates, rp_si_keep), collapse = ", "), "\n", sep = "")
  
  # ============================================================================
  # 4. Compare both datasets without forcing identical measured-gene sets
  # ============================================================================
  
  rp_shared <- base::intersect(rp_xpo_genes, rp_si_genes)
  stopifnot(identical(
    rp_xpo_genes[rp_xpo_genes %in% rp_shared],
    rp_si_genes[rp_si_genes %in% rp_shared]
  ))
  
  cat("\nShared eligible genes: ", length(rp_shared), "\n", sep = "")
  cat("Only in XPO1i: ", paste(setdiff(rp_xpo_genes, rp_si_genes), collapse = ", "), "\n", sep = "")
  cat("Only in siNMD3: ", paste(setdiff(rp_si_genes, rp_xpo_genes), collapse = ", "), "\n", sep = "")
  
  rp_qc_report <- data.frame(
    Gene = rp_sort_numeric(unique(c(rp_xpo_candidates, rp_si_candidates))),
    stringsAsFactors = FALSE
  )
  rp_qc_report$XPO1i_PresentAndQC <- rp_qc_report$Gene %in% rp_xpo_genes
  rp_qc_report$siNMD3_PresentAndQC <- rp_qc_report$Gene %in% rp_si_genes
  utils::write.csv(
    rp_qc_report, file.path(rp_out, "RPL_RPS_QC_Membership_Comparison.csv"),
    row.names = FALSE
  )
  
  cat("\nSaved PDFs:\n", rp_xpo_pdf, "\n", rp_si_pdf, "\n", sep = "")
  cat("Analysis complete: numeric RP order, no legends or annotations.\n")
  
  
  
  
  
}))

# =============================================================================
# FigS7E_F
# =============================================================================
nmd3_step("FigS7E_F", c("Fig4D","Fig4E"), quote({
  # =============================================================================
  # PART 7: MYC TARGETS V1 RESIDUAL AND B-CELL HALF-LIFE SENSITIVITY
  # =============================================================================
  
  
  # NEW MYC V1 SENSITIVITY: uses the core original df_F and df_G residuals.
  # Primary Fig.4D/E full-set plots were already calculated by core workflow.
  stopifnot(exists("df_F"), exists("df_G"))
  myc4_saved <- list(df_F = df_F, df_G = df_G)
  myc4_out <- file.path(.nmd3_run_root, "extra_sensitivity", "MYC_V1_NoRPLRPS_HalfLife")
  dir.create(myc4_out, recursive=TRUE, showWarnings=FALSE)
  myc4_half_life_override <- Sys.getenv("NMD3_HALF_LIFE_FILE", unset = "")
  myc4_half_life_url <- paste0("https://zenodo.org/records/14827610/files/",
                               "Savitski_halflives.csv?download=1")
  
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(ggrepel)
    library(tibble)
  })
  
  myc4_out <- file.path(.nmd3_run_root, "extra_sensitivity", "MYC_V1_NoRPLRPS_HalfLife")
  dir.create(myc4_out, recursive = TRUE, showWarnings = FALSE)
  cat("Output folder: ", myc4_out, "\n", sep = "")
  
  # Half-life file: override ONLY if your local copy is elsewhere.
  # Previous workflow used this CSV filename, not a fixed absolute Mac path.
  myc4_half_life_override <- Sys.getenv("NMD3_HALF_LIFE_FILE", unset = "")
  myc4_half_life_url <- paste0(
    "https://zenodo.org/records/14827610/files/",
    "Savitski_halflives.csv?download=1"
  )
  
  # 1. REUSE ORIGINAL FIG 4F/4G GENE UNIVERSES ------------------------------
  # df_F stores the original RNA log2FC and MYC V1 classification.
  # df_G stores the actual original robust-fit residual, NOT to be recalculated.
  myc4_F <- as.data.frame(myc4_saved[["df_F"]]) %>%
    dplyr::transmute(
      gene = toupper(trimws(as.character(gene))),
      group = as.character(group),
      log2FC_RNA = as.numeric(log2FC_RNA)
    ) %>%
    dplyr::filter(!is.na(gene), nzchar(gene),
                  group %in% c("Other", "MYC"), is.finite(log2FC_RNA)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(group = factor(group, levels = c("Other", "MYC")))
  
  myc4_G <- as.data.frame(myc4_saved[["df_G"]]) %>%
    dplyr::transmute(
      gene = toupper(trimws(as.character(gene))),
      group = as.character(group),
      residual_old = as.numeric(residual_old)
    ) %>%
    dplyr::filter(!is.na(gene), nzchar(gene),
                  group %in% c("Other", "MYC"), is.finite(residual_old)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(group = factor(group, levels = c("Other", "MYC")))
  
  stopifnot(!anyDuplicated(myc4_F$gene), !anyDuplicated(myc4_G$gene),
            all(table(myc4_F$group) >= 3), all(table(myc4_G$group) >= 3))
  myc4_overlap <- inner_join(
    dplyr::select(myc4_F, gene, group_F = group),
    dplyr::select(myc4_G, gene, group_G = group),
    by = "gene"
  )
  if (any(as.character(myc4_overlap$group_F) !=
          as.character(myc4_overlap$group_G))) {
    stop("Original df_F/df_G MYC gene memberships disagree!")
  }
  cat("Original Figure 4F matched genes: ", nrow(myc4_F), "\n", sep = "")
  cat("Original Figure 4G matched genes: ", nrow(myc4_G), "\n", sep = "")
  if (!setequal(myc4_F$gene, myc4_G$gene)) {
    warning("Original Fig4F/4G gene universes differ slightly; preserving each ",
            "figure's saved universe instead of forcing an intersection.")
  }
  
  myc4_is_RP <- function(g) {
    g <- toupper(trimws(as.character(g)))
    grepl("^RPL", g) | (grepl("^RPS", g) & !grepl("^RPS6K", g))
  }
  myc4_G$excluded_RPL_RPS <- myc4_is_RP(myc4_G$gene)
  myc4_G_noRP <- myc4_G %>% dplyr::filter(!excluded_RPL_RPS)
  cat("RPL/RPS excluded for analyses 3 and 4: ",
      sum(myc4_G$excluded_RPL_RPS), "\n", sep = "")
  print(table(myc4_G_noRP$group))
  stopifnot(all(table(myc4_G_noRP$group) >= 3))
  
  # 2. COMMON STATISTICS AND FIGURE 4F/G VISUAL STYLE -----------------------
  myc4_colors <- c("Other" = "#4DBBD5", "MYC" = "#F39B7F")
  myc4_labels <- c(
    "Other" = "All other\nproteins",
    "MYC" = "Hallmark MYC-\ntarget proteins"
  )
  myc4_p_label <- function(p) {
    if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
  }
  
  myc4_test <- function(df, value_col, analysis, label) {
    stopifnot(value_col %in% names(df), all(table(df$group) >= 3))
    a <- df[[value_col]][df$group == "Other"]
    b <- df[[value_col]][df$group == "MYC"]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    wt <- stats::wilcox.test(a, b, alternative = "two.sided", exact = FALSE)
    res <- data.frame(
      Analysis = analysis,
      Description = label,
      N_Other = length(a), N_MYC = length(b),
      Median_Other = stats::median(a), Median_MYC = stats::median(b),
      Delta_MYC_minus_Other = stats::median(b) - stats::median(a),
      Wilcoxon_P = unname(wt$p.value),
      stringsAsFactors = FALSE
    )
    cat("\n========== ", label, " ==========\n", sep = "")
    print(res, row.names = FALSE)
    utils::write.csv(res, file.path(myc4_out, paste0(analysis, "_Statistics.csv")),
                     row.names = FALSE)
    res
  }
  
  myc4_draw_violin <- function(df, value_col, ylabel, test, y_limits,
                               bracket_y, highlight = FALSE) {
    statlabel <- paste0(
      myc4_p_label(test$Wilcoxon_P[1]), "\n\u0394 = ",
      sprintf("%.2f", abs(test$Delta_MYC_minus_Other[1]))
    )
    bracket <- data.frame(group1 = "Other", group2 = "MYC",
                          y.position = bracket_y, label = statlabel)
    g <- ggplot2::ggplot(
      df, ggplot2::aes(x = group, y = .data[[value_col]], fill = group)
    ) +
      ggplot2::geom_violin(color = NA, alpha = 0.85, trim = TRUE,
                           scale = "width") +
      ggplot2::geom_boxplot(width = 0.15, fill = "white", color = "black",
                            linewidth = 0.4, outlier.shape = NA) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey50", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = myc4_colors) +
      ggplot2::scale_x_discrete(labels = myc4_labels) +
      ggpubr::stat_pvalue_manual(
        bracket, label = "label", xmin = "group1", xmax = "group2",
        y.position = "y.position", tip.length = 0.03,
        bracket.size = 0.8, size = 5
      ) +
      ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
      ggplot2::labs(x = NULL, y = ylabel) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(face = "bold", color = "black"),
        axis.text.y = ggplot2::element_text(color = "black"),
        plot.margin = ggplot2::margin(10, 75, 10, 10)
      )
    if (highlight) {
      label_genes <- c("GSPT1", "EIF3D", "EEF1B2", "CCT5", "KPNB1")
      labels <- df %>%
        dplyr::filter(group == "MYC", gene %in% label_genes) %>%
        dplyr::arrange(dplyr::desc(.data[[value_col]]))
      if (nrow(labels) > 0L) {
        g <- g +
          ggplot2::geom_point(
            data = labels, ggplot2::aes(x = group, y = .data[[value_col]]),
            inherit.aes = FALSE, shape = 21, size = 2.2,
            fill = "gold", color = "black", stroke = 0.4
          ) +
          ggrepel::geom_text_repel(
            data = labels, ggplot2::aes(label = gene),
            nudge_x = 0.40, direction = "y", hjust = 0,
            size = 4, box.padding = 0.10, point.padding = 0.08,
            min.segment.length = 0, max.overlaps = Inf, seed = 123
          )
      }
    }
    g
  }
  myc4_save <- function(plot, filename, width = 5, height = 5) {
    print(plot)
    nmd3_save_pdf(file.path(myc4_out, filename), plot = plot,
                  width = width, height = height, units = "in",
                  device = grDevices::cairo_pdf)
  }
  
  
  
  myc4_res_range <- range(myc4_G$residual_old)
  myc4_res_bracket <- max(myc4_res_range) + 0.08 * diff(myc4_res_range)
  myc4_res_ylim <- c(-2.2, max(2.4, myc4_res_bracket + 0.2))
  # Baseline full-universe Wilcoxon tests were ALREADY run in original Fig.4F/G.
  # Reuse wt_F and wt_G from the core, with descriptive counts/medians only.
  stopifnot(exists("wt_F"), exists("wt_G"))
  myc4_from_existing <- function(dd, value, analysis, label, stored_p) {
    a <- dd[[value]][dd$group == "Other"]
    b <- dd[[value]][dd$group == "MYC"]
    data.frame(Analysis=analysis,Description=label,
               N_Other=length(a),N_MYC=length(b),
               Median_Other=median(a),Median_MYC=median(b),
               Delta_MYC_minus_Other=median(b)-median(a),
               Wilcoxon_P=unname(stored_p),stringsAsFactors=FALSE)
  }
  myc4_stats1 <- myc4_from_existing(myc4_F,"log2FC_RNA","1",
                                    "RNA log2FC — all matched genes",wt_F$p.value)
  myc4_stats2 <- myc4_from_existing(myc4_G,"residual_old","2",
                                    "Protein-RNA residual — all matched genes",wt_G$p.value)
  
  
  # 5. ANALYSIS 3: SAME SAVED RESIDUAL, EXCLUDING RPL/RPS -------------------
  # DO NOT refit! Only subset the ORIGINAL df_G.
  myc4_stats3 <- myc4_test(myc4_G_noRP, "residual_old", "3",
                           "Protein-RNA residual — No RPL/RPS")
  myc4_plot3 <- myc4_draw_violin(
    myc4_G_noRP, "residual_old",
    "Protein-RNA residual\n(siNMD3 vs siCtrl)",
    myc4_stats3, y_limits = myc4_res_ylim,
    bracket_y = myc4_res_bracket, highlight = TRUE
  )
  
  
  # FigS7E in the revised manuscript excludes GO:0042254 together with RPL/RPS.
  # Preserve the original residual definition; only change the membership filter.
  myc_rb_sets <- nmd3_cached_msig(species="Homo sapiens",collection="C5",subcollection="GO:BP")
  myc_rb_genes <- unique(toupper(myc_rb_sets$gene_symbol[myc_rb_sets$gs_name=="GOBP_RIBOSOME_BIOGENESIS"]))
  stopifnot(length(myc_rb_genes)>20)
  myc4_G_noRB <- myc4_G_noRP[!myc4_G_noRP$gene %in% myc_rb_genes,,drop=FALSE]
  myc4_stats_RB <- myc4_test(myc4_G_noRB,"residual_old","S7E_NoRB_NoRPLRPS","MYC residual: exclude RB and RPL/RPS")
  myc4_plot_RB <- myc4_draw_violin(myc4_G_noRB,"residual_old",
                                   "Non-RB protein-RNA residual (siNMD3 vs siCtrl)",myc4_stats_RB,
                                   y_limits=myc4_res_ylim,bracket_y=myc4_res_bracket,highlight=FALSE)
  myc4_save(myc4_plot_RB,"FigS7E_MYC_Residual_NoRB_NoRPLRPS.pdf")
  utils::write.csv(myc4_G_noRB,"FigS7E_MYC_Residual_Source.csv",row.names=FALSE)
  
  # 6. FIND THE REFERENCE HALF-LIFE CSV -------------------------------------
  # The original scripts used Savitski_halflives.csv without a fixed local path.
  # Check the project data directory or current run cache; otherwise download from Zenodo.
  myc4_candidates <- unique(c(
    myc4_half_life_override,
    file.path(.nmd3_data_root, "references", "Savitski_halflives.csv"),
    file.path(myc4_out, "Savitski_halflives.csv")
  ))
  myc4_candidates <- myc4_candidates[nzchar(myc4_candidates)]
  myc4_found <- myc4_candidates[file.exists(myc4_candidates)]
  if (length(myc4_found)) {
    myc4_hl_file <- myc4_found[1]
    message("Using half-life CSV: ", normalizePath(myc4_hl_file))
  } else {
    myc4_hl_file <- file.path(myc4_out, "Savitski_halflives.csv")
    message("Reference half-life CSV not found locally. Downloading from Zenodo:")
    message(myc4_half_life_url)
    status <- tryCatch(
      utils::download.file(myc4_half_life_url, myc4_hl_file,
                           mode = "wb", quiet = FALSE),
      error = function(e) e
    )
    if (inherits(status, "error") ||
        !file.exists(myc4_hl_file) ||
        file.info(myc4_hl_file)$size < 1000L) {
      stop("Automatic download failed. Download Savitski_halflives.csv manually ",
           "from https://zenodo.org/records/14827610 and set ",
           "myc4_half_life_override to its local path.\n",
           "Original error: ", if (inherits(status, "error"))
             conditionMessage(status) else "missing or very small file")
    }
  }
  
  myc4_hl_raw <- utils::read.csv(myc4_hl_file,
                                 check.names = TRUE, stringsAsFactors = FALSE)
  myc4_req_cols <- c("gene_name", "Bcells.replicate.1.half_life",
                     "Bcells.replicate.2.half_life")
  myc4_missing_cols <- setdiff(myc4_req_cols, colnames(myc4_hl_raw))
  if (length(myc4_missing_cols)) {
    stop("Unexpected reference half-life file format. Missing columns: ",
         paste(myc4_missing_cols, collapse = ", "))
  }
  
  myc4_rep1 <- suppressWarnings(as.numeric(
    myc4_hl_raw[["Bcells.replicate.1.half_life"]]))
  myc4_rep2 <- suppressWarnings(as.numeric(
    myc4_hl_raw[["Bcells.replicate.2.half_life"]]))
  for (i in 1:2) {
    qc_col <- paste0("Bcells.replicate.", i, ".dataQual")
    if (qc_col %in% names(myc4_hl_raw)) {
      good <- !is.na(myc4_hl_raw[[qc_col]]) &
        tolower(trimws(as.character(myc4_hl_raw[[qc_col]]))) == "good"
      if (i == 1) myc4_rep1[!good] <- NA_real_
      if (i == 2) myc4_rep2[!good] <- NA_real_
    } else {
      warning("QC column ", qc_col,
              " not in CSV. Using both valid positive replicates as in prior code.")
    }
  }
  
  # Replicate mean for each original record; median if duplicate gene symbols.
  myc4_hl_ref <- tibble::tibble(
    gene = toupper(trimws(as.character(myc4_hl_raw$gene_name))),
    rep1 = myc4_rep1, rep2 = myc4_rep2
  ) %>%
    dplyr::filter(!is.na(gene), nzchar(gene), !grepl(";|,", gene),
                  is.finite(rep1), is.finite(rep2), rep1 > 0, rep2 > 0) %>%
    dplyr::mutate(half_life_h = (rep1 + rep2) / 2) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(half_life_h = median(half_life_h), .groups = "drop")
  stopifnot(nrow(myc4_hl_ref) > 0, !anyDuplicated(myc4_hl_ref$gene))
  
  # 7. ANALYSIS 4: B-CELL HALF-LIFE, EXCLUDING RPL/RPS ----------------------
  # Start from SAME original Fig4G universe, then remove RPL/RPS; keep NAs for QC.
  myc4_hl_merged <- myc4_G_noRP %>%
    dplyr::left_join(myc4_hl_ref, by = "gene")
  myc4_hl_coverage <- myc4_hl_merged %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_with_half_life = sum(is.finite(half_life_h) & half_life_h > 0),
      coverage_pct = 100 * n_with_half_life / n_total,
      .groups = "drop"
    )
  cat("\n========== Half-life reference coverage ==========\n")
  print(myc4_hl_coverage)
  myc4_hl_plotdata <- myc4_hl_merged %>%
    dplyr::filter(is.finite(half_life_h), half_life_h > 0)
  stopifnot(all(table(myc4_hl_plotdata$group) >= 3))
  myc4_stats4 <- myc4_test(myc4_hl_plotdata, "half_life_h", "4",
                           "B-cell protein half-life — No RPL/RPS")
  
  myc4_n_other <- sum(myc4_hl_plotdata$group == "Other")
  myc4_n_myc <- sum(myc4_hl_plotdata$group == "MYC")
  myc4_plot4 <- ggplot2::ggplot(
    myc4_hl_plotdata,
    ggplot2::aes(x = group, y = half_life_h, fill = group)
  ) +
    ggplot2::geom_violin(color = NA, alpha = 0.75,
                         trim = TRUE, scale = "width") +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA,
                          fill = "white", color = "black", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = myc4_colors) +
    ggplot2::scale_x_discrete(labels = c(
      "Other" = paste0("Other\n(n = ", myc4_n_other, ")"),
      "MYC" = paste0("MYC targets V1\n(n = ", myc4_n_myc, ")")
    )) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = NULL,
      y = "Protein half-life (hours; log10 scale)",
      title = "MYC targets V1 vs Other (without RPL/RPS)",
      subtitle = paste0("Human B-cell reference | Wilcoxon ",
                        myc4_p_label(myc4_stats4$Wilcoxon_P[1]))
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  myc4_save(myc4_plot4, "FigS7F_MYC_HalfLife_NoRPLRPS.pdf", 5.7, 5.2)
  
  # 8. EXPORT REPRODUCIBLE SOURCE DATA AND SUMMARY --------------------------
  utils::write.csv(myc4_F,
                   file.path(myc4_out, "1_RNA_All_Source.csv"), row.names = FALSE)
  utils::write.csv(myc4_G,
                   file.path(myc4_out, "2_Residual_All_Source.csv"), row.names = FALSE)
  utils::write.csv(myc4_G_noRP,
                   file.path(myc4_out, "3_Residual_NoRPLRPS_Source.csv"),
                   row.names = FALSE)
  utils::write.csv(myc4_hl_merged,
                   file.path(myc4_out, "4_HalfLife_NoRPLRPS_Source_With_NA.csv"),
                   row.names = FALSE)
  utils::write.csv(myc4_hl_coverage,
                   file.path(myc4_out, "4_HalfLife_GroupCoverage.csv"),
                   row.names = FALSE)
  
  myc4_results <- dplyr::bind_rows(
    myc4_stats1, myc4_stats2, myc4_stats3, myc4_stats4
  )
  utils::write.csv(myc4_results,
                   file.path(myc4_out, "All_Four_Analyses_Summary.csv"),
                   row.names = FALSE)
  cat("\n========== ALL FOUR COMPARISONS ==========\n")
  print(myc4_results, row.names = FALSE)
  cat("\nCOMPLETED: four PDFs + per-gene CSVs + summary in:\n",
      myc4_out, "\n", sep = "")
  
  
  # ============================================================
  
  
  
}))

# =============================================================================
# FigS8A
# =============================================================================
nmd3_step("FigS8A", character(), quote({
  suppressPackageStartupMessages({
    library(readxl)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(ggrepel)
  })
  
  
  # ============================================================
  # 1. Derive the 60 proteins exactly as in Figure S5D Venn
  #    siNMD3 vs siCTRL, cytoplasmic fraction
  #    Tukey P < 0.05 and log2FC < 0
  # ============================================================
  
  # Reuse the in-house proteomics table loaded once for Figures 2G-2I.
  stopifnot(
    exists("prot_raw", inherits = TRUE),
    "Gene" %in% colnames(prot_raw)
  )
  
  prot_raw_5B <- as.data.frame(
    prot_raw,
    check.names = FALSE
  )
  
  stopifnot("Gene" %in% colnames(prot_raw_5B))
  
  log2fc_column_5B <-
    "Log2FC_siRNAxNMD3xKD_Cyto.siRNAxCtrl_Cyto"
  
  pvalue_column_5B <-
    "Tukey_siRNAxNMD3xKD_Cyto.siRNAxCtrl_Cyto"
  
  stopifnot(
    all(
      c(log2fc_column_5B, pvalue_column_5B) %in%
        colnames(prot_raw_5B)
    )
  )
  
  protein_statistics_5B <- prot_raw_5B %>%
    dplyr::transmute(
      Gene = toupper(trimws(as.character(Gene))),
      log2FC = suppressWarnings(
        as.numeric(.data[[log2fc_column_5B]])
      ),
      P.Value = suppressWarnings(
        as.numeric(.data[[pvalue_column_5B]])
      )
    ) %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      is.finite(log2FC),
      is.finite(P.Value),
      P.Value >= 0,
      P.Value <= 1
    )
  
  downregulated_proteins_5B <- protein_statistics_5B %>%
    dplyr::filter(
      P.Value < 0.05,
      log2FC < 0
    ) %>%
    dplyr::arrange(P.Value, log2FC) %>%
    dplyr::distinct(Gene, .keep_all = TRUE)
  
  genes_5B <- sort(
    unique(downregulated_proteins_5B$Gene)
  )
  
  stopifnot(
    length(genes_5B) == 60L,
    !anyNA(genes_5B),
    !anyDuplicated(genes_5B),
    all(nzchar(genes_5B))
  )
  
  cat(
    "Figure 5B cytoplasmic proteins with Tukey P < 0.05 and log2FC < 0:",
    length(genes_5B),
    "\n"
  )
  
  utils::write.csv(
    protein_statistics_5B,
    "Figure5B_siNMD3_vs_siCTRL_cytoplasmic_protein_statistics.csv",
    row.names = FALSE
  )
  
  utils::write.csv(
    downregulated_proteins_5B,
    "Figure5B_60_downregulated_cytoplasmic_proteins.csv",
    row.names = FALSE
  )
  
  # ============================================================
  # 2. Map gene symbols to Entrez IDs
  # ============================================================
  gene_entrez_5B <- clusterProfiler::bitr(
    genes_5B,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  stopifnot(
    nrow(gene_entrez_5B) > 0L,
    "ENTREZID" %in% colnames(gene_entrez_5B)
  )
  
  cat(
    "Successfully mapped proteins:",
    dplyr::n_distinct(gene_entrez_5B$SYMBOL),
    "of",
    length(genes_5B),
    "\n"
  )
  # ------------------------------------------------------------
  # 3. GO Biological Process ORA
  # ------------------------------------------------------------
  
  ego_5B <- clusterProfiler::enrichGO(
    gene = unique(gene_entrez_5B$ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  bp_results_5B <- as.data.frame(
    ego_5B
  )
  
  stopifnot(
    nrow(bp_results_5B) > 0L
  )
  
  # ------------------------------------------------------------
  # 4. Retain pathways with FDR < 0.25
  # ------------------------------------------------------------
  
  bp_filtered_5B <- bp_results_5B %>%
    dplyr::filter(
      is.finite(qvalue),
      qvalue < 0.25
    ) %>%
    dplyr::arrange(qvalue) %>%
    dplyr::mutate(
      desc_lower = tolower(Description)
    )
  
  stopifnot(
    nrow(bp_filtered_5B) > 0L
  )
  
  # ------------------------------------------------------------
  # 5. Pathways displayed in Figure 5B
  # ------------------------------------------------------------
  
  selected_pathways_5B <- c(
    "cytoplasmic translation",
    "cytoplasmic translational initiation",
    "translational initiation",
    "positive regulation of translation",
    "regulation of translational initiation",
    "regulation of mRNA stability",
    "regulation of translation",
    "regulation of RNA stability"
  )
  
  pick_pathway_5B <- function(
    target,
    result_table
  ) {
    
    target_lower <- tolower(target)
    
    exact_index <- which(
      result_table$desc_lower == target_lower
    )
    
    if (length(exact_index) > 0L) {
      
      exact_index <- exact_index[
        which.min(
          result_table$qvalue[
            exact_index
          ]
        )
      ]
      
      result <- result_table[
        exact_index,
        ,
        drop = FALSE
      ]
      
      result$SelectedLabel <- target
      
      return(result)
    }
    
    partial_index <- grep(
      target_lower,
      result_table$desc_lower,
      fixed = TRUE
    )
    
    if (length(partial_index) > 0L) {
      
      partial_index <- partial_index[
        which.min(
          result_table$qvalue[
            partial_index
          ]
        )
      ]
      
      result <- result_table[
        partial_index,
        ,
        drop = FALSE
      ]
      
      result$SelectedLabel <- target
      
      return(result)
    }
    
    NULL
  }
  
  selected_list_5B <- lapply(
    selected_pathways_5B,
    pick_pathway_5B,
    result_table = bp_filtered_5B
  )
  
  sel_df_5B <- dplyr::bind_rows(
    Filter(
      Negate(is.null),
      selected_list_5B
    )
  )
  
  stopifnot(
    nrow(sel_df_5B) == length(selected_pathways_5B),
    !anyDuplicated(sel_df_5B$SelectedLabel)
  )
  
  # ------------------------------------------------------------
  # 6. Pathway–gene connections
  # ------------------------------------------------------------
  
  edge_tbl_5B <- sel_df_5B %>%
    dplyr::select(
      SelectedLabel,
      geneID,
      qvalue
    ) %>%
    dplyr::mutate(
      geneID = as.character(geneID)
    ) %>%
    dplyr::filter(
      !is.na(geneID),
      nchar(geneID) > 0
    ) %>%
    tidyr::separate_rows(
      geneID,
      sep = "/"
    ) %>%
    dplyr::transmute(
      term = SelectedLabel,
      gene = toupper(
        trimws(geneID)
      ),
      qvalue = qvalue
    ) %>%
    dplyr::distinct(
      term,
      gene,
      .keep_all = TRUE
    )
  
  stopifnot(
    nrow(edge_tbl_5B) > 0L
  )
  
  # ------------------------------------------------------------
  # 7. Gene degree
  # ------------------------------------------------------------
  
  degree_df_5B <- edge_tbl_5B %>%
    dplyr::count(
      gene,
      name = "degree"
    ) %>%
    dplyr::arrange(
      dplyr::desc(degree),
      gene
    )
  
  print(degree_df_5B)
  
  # ------------------------------------------------------------
  # 8. Display labels
  # ------------------------------------------------------------
  
  display_labels_5B <- c(
    "cytoplasmic translation" =
      "Cytoplasmic translation",
    
    "cytoplasmic translational initiation" =
      "Cytoplasmic translation initiation",
    
    "translational initiation" =
      "Translation initiation",
    
    "positive regulation of translation" =
      "Positive regulation of translation",
    
    "regulation of translational initiation" =
      "Regulation of translation initiation",
    
    "regulation of mRNA stability" =
      "Regulation of mRNA stability",
    
    "regulation of translation" =
      "Regulation of translation",
    
    "regulation of RNA stability" =
      "Regulation of RNA stability"
  )
  
  # ------------------------------------------------------------
  # 9. Nodes
  # ------------------------------------------------------------
  
  term_nodes_5B <- sel_df_5B %>%
    dplyr::transmute(
      name = SelectedLabel,
      type = "term",
      FDR = pmax(
        qvalue,
        .Machine$double.xmin
      )
    ) %>%
    dplyr::distinct(
      name,
      .keep_all = TRUE
    )
  
  gene_nodes_5B <- degree_df_5B %>%
    dplyr::transmute(
      name = gene,
      type = "gene",
      FDR = NA_real_
    )
  
  nodes_5B <- dplyr::bind_rows(
    term_nodes_5B,
    gene_nodes_5B
  )
  
  # ------------------------------------------------------------
  # 10. Two-column node positions
  # ------------------------------------------------------------
  
  term_order_5B <- sel_df_5B %>%
    dplyr::arrange(qvalue) %>%
    dplyr::distinct(SelectedLabel) %>%
    dplyr::pull(SelectedLabel)
  
  gene_order_5B <- degree_df_5B$gene
  
  term_y_center_5B <- length(
    gene_order_5B
  ) / 2
  
  term_pos_5B <- data.frame(
    name = term_order_5B,
    x = 0,
    y = seq(
      from =
        term_y_center_5B +
        length(term_order_5B) / 2,
      to =
        term_y_center_5B -
        length(term_order_5B) / 2 + 1,
      length.out =
        length(term_order_5B)
    )
  )
  
  gene_pos_5B <- data.frame(
    name = gene_order_5B,
    x = 1,
    y = rev(
      seq_along(gene_order_5B)
    )
  )
  
  position_df_5B <- dplyr::bind_rows(
    term_pos_5B,
    gene_pos_5B
  )
  
  node_df_5B <- nodes_5B %>%
    dplyr::left_join(
      position_df_5B,
      by = "name"
    ) %>%
    dplyr::mutate(
      neglog10FDR = ifelse(
        type == "term",
        -log10(FDR),
        NA_real_
      ),
      label = ifelse(
        type == "term",
        unname(
          display_labels_5B[name]
        ),
        name
      )
    )
  
  stopifnot(
    !anyNA(node_df_5B$x),
    !anyNA(node_df_5B$y),
    !anyNA(node_df_5B$label)
  )
  
  # ------------------------------------------------------------
  # 11. Edge plotting table
  # ------------------------------------------------------------
  
  edge_plot_5B <- edge_tbl_5B %>%
    dplyr::rename(
      from = term,
      to = gene
    ) %>%
    dplyr::left_join(
      node_df_5B %>%
        dplyr::select(
          name,
          x,
          y
        ),
      by = c(
        "from" = "name"
      )
    ) %>%
    dplyr::rename(
      x_from = x,
      y_from = y
    ) %>%
    dplyr::left_join(
      node_df_5B %>%
        dplyr::select(
          name,
          x,
          y
        ),
      by = c(
        "to" = "name"
      )
    ) %>%
    dplyr::rename(
      x_to = x,
      y_to = y
    )
  
  stopifnot(
    !anyNA(edge_plot_5B$x_from),
    !anyNA(edge_plot_5B$y_from),
    !anyNA(edge_plot_5B$x_to),
    !anyNA(edge_plot_5B$y_to)
  )
  
  # ------------------------------------------------------------
  # 12. Figure 5B
  # ------------------------------------------------------------
  
  palette_5B <- c(
    "#4575B4",
    "#74ADD1",
    "#ABD9E9",
    "#FEE090",
    "#FDAE61",
    "#F46D43",
    "#D73027",
    "#A50026"
  )
  
  p_5B <- ggplot() +
    geom_segment(
      data = edge_plot_5B,
      aes(
        x = x_from,
        y = y_from,
        xend = x_to,
        yend = y_to
      ),
      linewidth = 0.35,
      alpha = 0.45,
      color = "grey30"
    ) +
    geom_point(
      data = subset(
        node_df_5B,
        type == "term"
      ),
      aes(
        x = x,
        y = y,
        fill = neglog10FDR
      ),
      shape = 21,
      color = "black",
      stroke = 0.6,
      size = 4.8
    ) +
    geom_point(
      data = subset(
        node_df_5B,
        type == "gene"
      ),
      aes(
        x = x,
        y = y
      ),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.6,
      size = 3
    ) +
    ggrepel::geom_text_repel(
      data = subset(
        node_df_5B,
        type == "term"
      ),
      aes(
        x = x,
        y = y,
        label = label
      ),
      direction = "y",
      nudge_x = -0.32,
      hjust = 1,
      size = 3.8,
      fontface = "bold",
      segment.color = "grey60",
      segment.size = 0.28,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.66,
      max.overlaps = Inf
    ) +
    ggrepel::geom_text_repel(
      data = subset(
        node_df_5B,
        type == "gene"
      ),
      aes(
        x = x,
        y = y,
        label = label
      ),
      direction = "y",
      nudge_x = 0.28,
      hjust = 0,
      size = 3.1,
      segment.color = "grey70",
      segment.size = 0.25,
      min.segment.length = 0,
      box.padding = 0.12,
      point.padding = 0.4,
      max.overlaps = Inf
    ) +
    scale_fill_gradientn(
      colours = palette_5B,
      name = expression(
        -log[10]~FDR
      ),
      na.value = "white"
    ) +
    scale_x_continuous(
      limits = c(-0.9, 1.6),
      breaks = NULL,
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      breaks = NULL
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    labs(
      title = expression(
        italic("siNMD3")~
          "cytoplasmic downregulated proteins"
      ),
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 12.5
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.1
      ),
      plot.margin = margin(
        10,
        28,
        10,
        40
      )
    )
  
  
  nmd3_save_pdf(
    filename =
      "FigS8A_siNMD3_Cytoplasmic_Downregulated_Protein_ORA.pdf",
    plot = p_5B,
    width = 12,
    height = 5,
    units = "in",
    device = cairo_pdf
  )
  
}))

# =============================================================================
# FigS8D_E
# =============================================================================
nmd3_step("FigS8D_E", character(), quote({
  # FigS8D_E — analysis below
  # Cytoplasmic RNA GSEA: siNMD3 vs siCTRL
  #
  # S6F: NEPC signature
  # S6G: CRPC-Adeno signature
  # ============================================================
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tibble)
    library(fgsea)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
  })
  
  # ============================================================
  # 1. Cytoplasmic RNA samples
  # ============================================================
  
  # expr_matrix was loaded previously from Count_Matrix.txt.
  
  stopifnot(exists("expr_matrix"))
  
  cyto_samples <- c(
    "C1C", "C2C", "C3C",
    "N1C", "N2C", "N3C"
  )
  
  stopifnot(
    all(cyto_samples %in% colnames(expr_matrix))
  )
  
  cyto_counts <- as.matrix(
    expr_matrix[
      ,
      cyto_samples,
      drop = FALSE
    ]
  )
  
  storage.mode(cyto_counts) <- "integer"
  
  # ============================================================
  # 2. Original expression filter
  # Mean raw count >= 8 across the six cytoplasmic samples
  # ============================================================
  
  keep <- rowMeans(cyto_counts) >= 8
  
  cyto_counts <- cyto_counts[
    keep,
    ,
    drop = FALSE
  ]
  
  cat(
    "Genes retained after filtering:",
    nrow(cyto_counts),
    "\n"
  )
  
  # ============================================================
  # 3. DESeq2: siNMD3 vs siCTRL
  # ============================================================
  
  condition <- factor(
    c(
      "siCTRL", "siCTRL", "siCTRL",
      "siNMD3", "siNMD3", "siNMD3"
    ),
    levels = c("siCTRL", "siNMD3")
  )
  
  coldata <- data.frame(
    row.names = cyto_samples,
    condition = condition
  )
  
  dds_lineage <- DESeq2::DESeqDataSetFromMatrix(
    countData = cyto_counts,
    colData = coldata,
    design = ~ condition
  )
  
  dds_lineage <- DESeq2::DESeq(
    dds_lineage,
    quiet = TRUE
  )
  
  res_lineage <- DESeq2::results(
    dds_lineage,
    contrast = c(
      "condition",
      "siNMD3",
      "siCTRL"
    )
  )
  
  # ============================================================
  # 4. Original GSEA ranking
  # DESeq2 Wald statistic
  # Positive = higher in siNMD3
  # Negative = lower in siNMD3
  # ============================================================
  
  rank_df <- as.data.frame(res_lineage) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::transmute(
      gene = toupper(as.character(gene)),
      stat = as.numeric(stat)
    ) %>%
    dplyr::filter(is.finite(stat))
  
  ranks <- rank_df$stat
  names(ranks) <- rank_df$gene
  ranks <- sort(ranks, decreasing = TRUE)
  
  cat(
    "Genes used for ranking:",
    length(ranks),
    "\n"
  )
  
  # ============================================================
  # 5. Original NEPC 29-gene signature
  # ============================================================
  
  NEPC_29 <- c(
    "SOGA3", "RGS7", "JAKMIP2", "PROX1", "SYT11",
    "CAND2", "TRIM9", "SEC11C", "EZH2", "ETV5",
    "AURKA", "DNMT1", "MAP10", "SVOP", "ST8SIA3",
    "MYCN", "NRSN1", "KCNB2", "BRINP1", "SCG3",
    "SEZ6", "KIAA0408", "LRRC16B", "GPX2", "PCSK1",
    "KCND2", "ASXL3", "GNAO1", "C7ORF76"
  )
  
  # ============================================================
  # 6. Original CRPC-Adeno 41-gene signature
  # ============================================================
  
  CRPC_ADENO_41 <- c(
    "AR", "CCND1", "CREBBP", "CSDE1", "CYLD", "DICER1",
    "FHIT", "FOXP1", "HERPUD1", "MYH9", "NUP93", "PAX8",
    "RBBP6", "TRIM33", "GATA2", "MAPKAPK3", "SPDEF",
    "EFNA4", "EPN3", "EVPL", "HOXB13", "NKX3-1", "OPHN1",
    "PIEZO1", "PSCA", "RAB27B", "RGS10", "RIPK2",
    "SLC25A37", "TC2N", "RB1",
    
    # Removed by the original expression filter:
    "CATSPERB", "MMP2", "KLK3", "PRR5-ARHGAP8", "ARHGAP8",
    "CIITA", "LMAN1L", "KLK4", "UPK2", "SLC44A4"
  )
  
  stopifnot(
    length(NEPC_29) == 29L,
    length(CRPC_ADENO_41) == 41L
  )
  
  # ============================================================
  # 7. GSEA pathways
  # ============================================================
  
  pathways <- list(
    NEPC_signature_29 = NEPC_29,
    CRPC_Adeno_signature_41 = CRPC_ADENO_41
  )
  
  # ============================================================
  # 8. Check genes entering GSEA
  # ============================================================
  
  gene_check_summary <- dplyr::bind_rows(
    data.frame(
      signature = "NEPC_signature_29",
      gene = NEPC_29
    ),
    data.frame(
      signature = "CRPC_Adeno_signature_41",
      gene = CRPC_ADENO_41
    )
  ) %>%
    dplyr::mutate(
      used_in_GSEA = gene %in% names(ranks)
    ) %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(
      total_genes = dplyr::n(),
      used_in_GSEA = sum(used_in_GSEA),
      removed = total_genes - used_in_GSEA,
      used_genes = paste(
        gene[used_in_GSEA],
        collapse = ", "
      ),
      removed_genes = paste(
        gene[!used_in_GSEA],
        collapse = ", "
      ),
      .groups = "drop"
    )
  
  print(gene_check_summary, width = Inf)
  
  # Historical sanity check.
  
  stopifnot(
    sum(NEPC_29 %in% names(ranks)) == 12L,
    sum(CRPC_ADENO_41 %in% names(ranks)) == 31L
  )
  
  # ============================================================
  # 9. fgsea
  # Original analysis: 100,000 permutations
  # ============================================================
  
  set.seed(123)
  
  fgsea_res <- fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    nperm = 100000
  ) %>%
    as.data.frame() %>%
    dplyr::arrange(padj)
  
  fgsea_res_simple <- fgsea_res %>%
    dplyr::select(
      pathway,
      size,
      ES,
      NES,
      pval,
      padj,
      leadingEdge
    )
  
  print(fgsea_res_simple)
  
  # ============================================================
  # 10. Save GSEA statistics
  # ============================================================
  
  fgsea_res_export <- fgsea_res_simple %>%
    dplyr::mutate(
      leadingEdge = vapply(
        leadingEdge,
        function(x) paste(x, collapse = ";"),
        character(1)
      )
    )
  
  write.csv(
    fgsea_res_export,
    "FigureS6FG_lineage_GSEA_fgsea_results.csv",
    row.names = FALSE
  )
  
  write.csv(
    gene_check_summary,
    "FigureS6FG_lineage_GSEA_genes_used.csv",
    row.names = FALSE
  )
  
  # ============================================================
  # 11. clusterProfiler object for gseaplot2
  # fgsea remains the source of NES/FDR statistics.
  # This object is used only to draw the enrichment curves.
  # ============================================================
  
  term2gene_lineage <- dplyr::bind_rows(
    data.frame(
      term = "NEPC_signature_29",
      gene = NEPC_29
    ),
    data.frame(
      term = "CRPC_Adeno_signature_41",
      gene = CRPC_ADENO_41
    )
  )
  
  set.seed(123)
  
  gsea_plot_object <- clusterProfiler::GSEA(
    geneList = ranks,
    TERM2GENE = term2gene_lineage,
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 1,
    seed = TRUE,
    verbose = FALSE
  )
  
  # ============================================================
  # 12. Extract fgsea statistics for labels
  # ============================================================
  
  nepc_result <- fgsea_res_simple %>%
    dplyr::filter(
      pathway == "NEPC_signature_29"
    )
  
  adeno_result <- fgsea_res_simple %>%
    dplyr::filter(
      pathway == "CRPC_Adeno_signature_41"
    )
  
  stopifnot(
    nrow(nepc_result) == 1L,
    nrow(adeno_result) == 1L
  )
  
  nepc_title <- sprintf(
    "NEPC signature\nNES = %.2f, FDR = %.3f",
    nepc_result$NES,
    nepc_result$padj
  )
  
  adeno_title <- sprintf(
    "CRPC-Adeno signature\nNES = %.2f, FDR = %.3f",
    adeno_result$NES,
    adeno_result$padj
  )
  
  # ============================================================
  # 13. Figure S6F: NEPC GSEA plot
  # ============================================================
  
  p_S6F <- enrichplot::gseaplot2(
    gsea_plot_object,
    geneSetID = "NEPC_signature_29",
    title = nepc_title,
    color = "#4DBBD5FF",
    base_size = 14,
    subplots = 1:3,
    pvalue_table = FALSE
  )
  
  
  # ============================================================
  # 14. Figure S6G: CRPC-Adeno GSEA plot
  # ============================================================
  
  p_S6G <- enrichplot::gseaplot2(
    gsea_plot_object,
    geneSetID = "CRPC_Adeno_signature_41",
    title = adeno_title,
    color = "#E64B35FF",
    base_size = 14,
    subplots = 1:3,
    pvalue_table = FALSE
  )
  
  
  # ============================================================
  # 15. Save PDFs
  # ============================================================
  
  nmd3_save_pdf(
    "FigS8D_signature_GSEA.pdf",
    p_S6F,
    width = 6,
    height = 6,
    device = cairo_pdf
  )
  
  nmd3_save_pdf(
    "FigS8E_Adeno_signature_GSEA.pdf",
    p_S6G,
    width = 6,
    height = 6,
    device = cairo_pdf
  )
  
  
  
  
  
  
  
  
  
  # =============================================================================
  # END OF SCRIPT
  # =============================================================================
  
  
  
  
  
}))

# =============================================================================
# FigS8F
# =============================================================================
nmd3_step("FigS8F", c("Prepare_SU2C"), quote({
  # FigS8F — analysis below
  # NMD3/PEG10/ENO2 signature and overall survival
  #
  # Algorithm preserved:
  # gene-wise z-score across 80 patients
  # + unweighted mean + median cutoff
  # + ordinary log-rank P
  # ============================================================
  
  stopifnot(
    exists("su2c_patient_data"),
    exists("logrank_p_SU2C"),
    exists("signature_genes_6K")
  )
  
  survival_data_6K <- su2c_patient_data %>%
    dplyr::select(
      PATIENT_ID,
      SAMPLE_ID,
      Number_of_samples,
      OS_MONTHS,
      OS_STATUS,
      dplyr::all_of(
        signature_genes_6K
      )
    )
  
  signature_z_matrix_6K <- scale(
    survival_data_6K[
      ,
      signature_genes_6K,
      drop = FALSE
    ]
  )
  
  stopifnot(
    all(
      is.finite(
        signature_z_matrix_6K
      )
    )
  )
  
  colnames(signature_z_matrix_6K) <-
    paste0(
      signature_genes_6K,
      "_z"
    )
  
  survival_data_6K <- dplyr::bind_cols(
    survival_data_6K,
    as.data.frame(
      signature_z_matrix_6K
    )
  )
  
  survival_data_6K$Signature_Score <-
    rowMeans(
      signature_z_matrix_6K
    )
  
  cutoff_value_6K <- median(
    survival_data_6K$Signature_Score
  )
  
  survival_data_6K$Signature_Group <-
    factor(
      ifelse(
        survival_data_6K$Signature_Score >
          cutoff_value_6K,
        "High",
        "Low"
      ),
      levels = c(
        "High",
        "Low"
      )
    )
  
  fit_6K <- survival::survfit(
    survival::Surv(
      OS_MONTHS,
      OS_STATUS
    ) ~ Signature_Group,
    data = survival_data_6K
  )
  
  logrank_6K <- survival::survdiff(
    survival::Surv(
      OS_MONTHS,
      OS_STATUS
    ) ~ Signature_Group,
    data = survival_data_6K
  )
  
  p_value_6K <- logrank_p_SU2C(
    logrank_6K
  )
  
  cat(
    "\nFigure 6K (SU2C) median cutoff:",
    cutoff_value_6K,
    "\n"
  )
  
  cat(
    "Figure 6K (SU2C) ordinary log-rank P:",
    p_value_6K,
    "\n"
  )
  
  fig6k <- survminer::ggsurvplot(
    fit_6K,
    data = survival_data_6K,
    
    pval = TRUE,
    conf.int = TRUE,
    risk.table = FALSE,
    censor = TRUE,
    
    legend.title =
      "NMD3/PEG10/ENO2",
    
    legend.labs = c(
      "High",
      "Low"
    ),
    
    palette = c(
      "High" = "#E31A1C",
      "Low" = "#1F78B4"
    ),
    
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    
    ggtheme =
      ggplot2::theme_minimal(),
    
    pval.size = 5,
    pval.coord = c(
      100,
      0.9
    )
  )
  
  print(fig6k$plot)
  
  source_data_6K <- survival_data_6K %>%
    dplyr::mutate(
      Signature_Cutoff =
        cutoff_value_6K,
      
      Cutoff_Method =
        "Median",
      
      P_Value_Method =
        "Ordinary log-rank",
      
      Logrank_P =
        p_value_6K
    )
  
  utils::write.csv(
    source_data_6K,
    file = paste0(
      "SourceData_Figure6K_SU2C_NMD3_PEG10_ENO2_",
      "signature_patient_level.csv"
    ),
    row.names = FALSE
  )
  
  nmd3_save_pdf(
    filename =
      "FigS8F_SU2C_NMD3_PEG10_ENO2_Signature_Survival.pdf",
    plot = fig6k$plot,
    width = 8,
    height = 6,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  
  
  
  
  # ============================================================
  
  
}))

# FINAL: only after all main and supplementary modules have executed.
nmd3_finish()
invisible(.nmd3_run_root)
}

# Run once, from Figure 1 through Supplementary Figure S8.
nmd3_run_main()
