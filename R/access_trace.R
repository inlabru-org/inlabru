access_trace <- function(name) {
  method <- fmesher::fm_caller_name(-1)
  access_method <- if (grepl("\\$", method)) {
    "$"
  } else if (grepl("\\[\\[", method)) {
    "[["
  } else {
    NA_character_
  }
  method_text <- if (access_method == "$") {
    "\\$"
  } else if (access_method == "[[") {
    "\\[\\["
  } else {
    ""
  }
  open_text <- if (access_method == "$") {
    "$"
  } else if (access_method == "[[") {
    "[['"
  } else {
    ""
  }
  close_text <- if (access_method == "$") {
    ""
  } else if (access_method == "[[") {
    "]]'"
  } else {
    ""
  }
  class_name <- sub(paste0("^", method_text, "\\."), "", method)
  caller <- fmesher::fm_caller_name(-3)
  print(glue::glue("{caller}: {class_name}{open_text}{name}{close_text}"))
}
