export_precise_results <- function(stage, output, d) {
  figure_dir <- file.path(output, "figures")
  table_dir <- file.path(output, "tables")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  copy_result <- function(source, destination) {
    if (!file.exists(source)) stop("Required result was not generated: ", source)
    if (!file.copy(source, destination, overwrite = TRUE)) stop("Could not write result: ", destination)
  }
  for (i in 2:3) {
    name <- sprintf("Fig_%02d.pdf", i)
    copy_result(file.path(stage, "safety", "figures", name), file.path(figure_dir, name))
  }
  for (i in 3:15) {
    module <- if (i <= 11L) "core" else if (i == 12L) "recent" else "timing"
    copy_result(file.path(stage, module, "figures", sprintf("Supplemental_Fig_%02d.pdf", i)),
                file.path(figure_dir, sprintf("Supplemental_Fig_%02d.pdf", i)))
  }
  draw_precise_ecg(d, file.path(figure_dir, "Supplemental_Fig_01.pdf"))
  draw_precise_utilization(d, file.path(figure_dir, "Supplemental_Fig_02.pdf"))
  for (i in 1:3)
    copy_result(file.path(stage, "core", "tables", sprintf("Main_Table_%d.csv", i)),
                file.path(table_dir, sprintf("Table_%d.csv", i)))
  supplementary <- c("core/tables/Supplemental_Table_01_Propensity_Model_Diagnostics.csv",
                     "core/tables/Supplemental_Table_02_Propensity_Distribution.csv",
                     sprintf("core/tables/Supplementary_Table_%d.csv", 3:7),
                     "core/tables/Supplemental_Table_08_Risk_History_Baseline.csv",
                     "recent/analysis/Supplemental_Table_09_Recent_ECG_Baseline.csv")
  for (i in seq_along(supplementary))
    copy_result(file.path(stage, supplementary[i]), file.path(table_dir, sprintf("Supplemental_Table_%d.csv", i)))
  invisible(output)
}

run_precise_ai <- function(input, output = "output", project = ".") {
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  d <- read_precise_data(input)
  stage <- tempfile("precise_ai_")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  run_analysis_engines(d, stage)
  export_precise_results(stage, output, d)
  message("Saved 17 PDF figures and 12 CSV tables in ", output)
  invisible(output)
}
