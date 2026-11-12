# Run from the repository root: Rscript scripts/render-tutorial.R
install.packages(".", repos = NULL, type = "source")

html_env <- new.env(parent = globalenv())
pdf_env <- new.env(parent = globalenv())

rmarkdown::render(
  "vignettes/JointShocks.Rmd",
  output_format = "rmarkdown::html_vignette",
  output_file = "index.html",
  output_dir = normalizePath("docs"),
  envir = html_env,
  clean = TRUE
)

rmarkdown::render(
  "vignettes/JointShocks.Rmd",
  output_format = "pdf_document",
  output_file = "JointShocks.pdf",
  output_dir = normalizePath("vignettes"),
  envir = pdf_env,
  clean = TRUE
)

objects <- c(
  "jIRF", "gIRF", "jfevd", "jIRF_CI", "gIRF_sum_CI", "gIRF_CI",
  "jointspill_results", "genspill_results",
  "joint_spill_results_roll", "gen_spill_result_roll"
)
checks <- vapply(objects, function(name) {
  identical(
    get(name, envir = html_env, inherits = FALSE),
    get(name, envir = pdf_env, inherits = FALSE)
  )
}, logical(1))
report <- paste(names(checks), ifelse(checks, "MATCH", "DIFFERENT"))
writeLines(c("HTML and PDF results from the same R Markdown source:", report),
           "consistency-check.txt")
if (!all(checks)) stop("HTML and PDF numerical results differ.")
message("All HTML and PDF numerical results match exactly.")
