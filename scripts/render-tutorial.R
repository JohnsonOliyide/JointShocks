# Run from the repository root: Rscript scripts/render-tutorial.R
# Install the current checkout before rendering its examples.
install.packages(".", repos = NULL, type = "source")
rmarkdown::render(
  "vignettes/JointShocks.Rmd",
  output_format = "rmarkdown::html_vignette",
  output_file = "index.html",
  output_dir = normalizePath("docs"),
  envir = new.env(parent = globalenv()),
  clean = TRUE
)
