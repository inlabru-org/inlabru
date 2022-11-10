install.packages("DiagrammeR")
#'  Update the mermaid library in the package
#' @param version A [character] of the desired version (the latest by default)
#' @return See value returned by [download.file]
updateMermaid <- function(version = "") {
  url <- "https://cdn.jsdelivr.net/npm/mermaid@version/dist/mermaid.min.js"
  if (version != "") {
    stopifnot(grepl("^[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+$", version))
    version <- paste0("@", version)
  }
  url <- gsub("@version", version, url)
  try(
    download.file(
      url,
      system.file("htmlwidgets/lib/mermaid/dist/mermaid.slim.min.js",
        package = "DiagrammeR"
      )
    )
  )
}
updateMermaid()
