if (requireNamespace("pkgload", quietly = TRUE) && !"joinme" %in% loadedNamespaces()) {
  # Use source tree during ad-hoc test runs (non-CRAN).
  repo_root <- normalizePath(file.path("..", ".."), mustWork = FALSE)
  if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
    repo_root <- normalizePath(".", mustWork = FALSE)
  }
  if (file.exists(file.path(repo_root, "DESCRIPTION"))) {
    try(pkgload::load_all(repo_root, quiet = TRUE, compile = FALSE), silent = TRUE)
  }
}
