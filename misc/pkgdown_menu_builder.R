library(glue)

get_info_original <- function(lines, pattern) {
  found <- grepl(
    pattern = pattern,
    x = lines,
    fixed = TRUE
  )
  if (!any(found)) {
    return(NULL)
  }
  begin <- min(which(found))
  version <- sub(
    pattern = "</small>",
    replacement = "",
    x = sub(
      pattern = pattern,
      replacement = "",
      x = trimws(lines[begin]),
      fixed = TRUE
    ),
    fixed = TRUE
  )
  list(original = TRUE, begin = begin, end = begin, version = version)
}
get_info_updated <- function(
  lines,
  pattern_begin,
  pattern_end,
  pattern_endtag
) {
  found <- grepl(
    pattern = pattern_begin,
    x = lines,
    fixed = TRUE
  )
  if (!any(found)) {
    return(NULL)
  }
  begin <- min(which(found))
  found <- grepl(
    pattern = pattern_end,
    x = lines,
    fixed = TRUE
  )
  if (!any(found)) {
    return(NULL)
  }
  end <- which(found)
  end <- end[end > begin]
  end <- min(end)
  if (length(end) == 0) {
    return(NULL)
  }
  version <- trimws(sub(
    pattern = pattern_endtag,
    replacement = "",
    x = sub(
      pattern = pattern_begin,
      replacement = "",
      x = trimws(lines[begin]),
      fixed = TRUE
    ),
    fixed = TRUE
  ))
  list(original = FALSE, begin = begin, end = end, version = version)
}

get_info_lines <- function(lines) {
  pattern <- '<small class="nav-text text-muted me-auto" data-bs-toggle="tooltip" data-bs-placement="bottom" title="">'
  info <- get_info_original(lines, pattern)
  if (!is.null(info)) {
    return(info)
  }
  pattern <- '<small class="nav-text text-muted me-auto" data-bs-toggle="tooltip" data-bs-placement="bottom" title="Released version">'
  info <- get_info_original(lines, pattern)
  if (!is.null(info)) {
    return(info)
  }
  pattern <- '<small class="nav-text text-danger me-auto" data-bs-toggle="tooltip" data-bs-placement="bottom" title="In-development version">'
  info <- get_info_original(lines, pattern)
  if (!is.null(info)) {
    return(info)
  }
  info <- get_info_updated(
    lines,
    "<!-- pkgdown version menu begin.",
    "<!-- pkgdown version menu end.",
    "-->"
  )
}

get_lines <- function(path) {
  lines <- tryCatch(
    readLines(path),
    error = function(e) {
      NULL
    },
    warning = function(w) {
      NULL
    }
  )
  lines
}
get_info_path <- function(path) {
  lines <- get_lines(path)
  if (is.null(lines)) {
    return(NULL)
  }
  get_info_lines(lines)
}

get_versions <- function(base_path) {
  current <- get_info_path(file.path(base_path, "index.html"))
  devel <- get_info_path(file.path(base_path, "dev", "index.html"))

  paths <- list.dirs(path = base_path, recursive = FALSE, full.names = TRUE) |>
    basename() |>
    grep(pattern = "^v[0-9]+\\.[0-9]+\\.[0-9]+", value = TRUE) |>
    sub(pattern = "^(v[0-9]+\\.[0-9]+\\.[0-9]+)", replacement = "\\1")
  info <- lapply(paths, function(x) {
    get_info_path(file.path(base_path, x, "index.html"))
  })

  old <- vapply(info, function(x) x$version, "")
  old <- old[order(package_version(old), decreasing = TRUE)]
  list(
    current = current$version,
    devel = devel$version,
    old = old
  )
}

get_version_paths <- function(versions) {
  c(
    get_version_path(versions$current, versions),
    get_version_path(versions$devel, versions),
    vapply(versions$old, function(x) get_version_path(x, versions), "")
  ) |>
    setNames(c(versions$current, versions$devel, versions$old))
}
get_version_path <- function(version, versions) {
  if (is.null(version)) {
    return(NULL)
  }
  if (identical(version, versions$current)) {
    ""
  } else if (identical(version, versions$devel)) {
    "dev/"
  } else if (version %in% versions$old) {
    glue("v{version}/")
  } else {
    NULL
  }
}

get_version_names <- function(versions) {
  c(
    get_version_name(versions$current, versions),
    get_version_name(versions$devel, versions),
    vapply(versions$old, function(x) get_version_name(x, versions), "")
  ) |>
    setNames(c(versions$current, versions$devel, versions$old))
}
get_version_name <- function(version, versions) {
  if (is.null(version)) {
    return(NULL)
  }
  if (identical(version, versions$current)) {
    glue("{version} (Latest release)")
  } else if (identical(version, versions$devel)) {
    glue("{version} (Development version)")
  } else if (version %in% versions$old) {
    glue("{version}")
  } else {
    NULL
  }
}


as_menu_line <- function(x, this_version, this_subpath, versions) {
  if (is.null(x)) {
    return(NULL)
  }
  this_dirname <- dirname(this_subpath)
  this_basename <- basename(this_subpath)

  if (identical(this_dirname, ".")) {
    this_dirname <- ""
    levels <- 0L
  } else {
    this_dirname <- paste0(this_dirname, "/")
    levels <- length(strsplit(this_dirname, split = "/")[[1]])
  }
  path_to_root <- glue_collapse(
    rep("../", levels + !identical(this_version, versions$current)),
    sep = ""
  )
  if (x == this_version) {
    path <- this_basename
  } else {
    version_path <- get_version_path(x, versions)
    path <- glue("{path_to_root}{version_path}{this_dirname}{this_basename}")
  }
  version_name <- get_version_name(x, versions)
  glue(
    "<li><a class=\"dropdown-item\" href=\"{path}\" title=\"{version_name}\">{version_name}</a></li>"
  )
}
get_menu_lines <- function(this_version, this_subpath, versions) {
  list(
    current = as_menu_line(
      versions$current,
      this_version,
      this_subpath,
      versions
    ),
    devel = as_menu_line(versions$devel, this_version, this_subpath, versions),
    old = vapply(
      versions$old,
      function(x) {
        as_menu_line(x, this_version, this_subpath, versions)
      },
      ""
    )
  )
}
get_menu <- function(this_version, this_subpath, versions) {
  menu_lines <- get_menu_lines(this_version, this_subpath, versions)
  if (length(menu_lines$old) == 0) {
    menu_old_versions <- NULL
  } else {
    menu_old_versions <- glue(
      '    <li><hr class="dropdown-divider"></li>
    <li><h6 class="dropdown-header" data-toc-skip>Older releases</h6></li>
    {glue_collapse(menu_lines$old, sep="\n    ")}
'
    )
  }
  version_name <- get_version_name(this_version, versions)
  menu <- glue(
    '
<!-- pkgdown version menu begin. {this_version} -->
<ul class="navbar-nav me-auto">
<li class="nav-item dropdown">
  <button class="nav-link dropdown-toggle" type="button" id="dropdown-versions"
    data-bs-toggle="dropdown" aria-expanded="false"
    aria-haspopup="true" title="{version_name}">{this_version}</button>
  <ul class="dropdown-menu" aria-labelledby="dropdown-versions">
    {menu_lines$current %||% ""}
    {menu_lines$devel %||% ""}
    {menu_old_versions %||% ""}
  </ul>
</li>
</ul>
<!-- pkgdown version menu end. -->
'
  )
  strsplit(menu, split = "\n")[[1]]
}

get_subpaths <- function(base_path, this_version, versions) {
  version_paths <- get_version_paths(versions)
  files <- dir(base_path, recursive = TRUE, include.dirs = FALSE)
  files <- files[grepl(pattern = "\\.html$", x = files)]
  if (identical(this_version, versions$current)) {
    for (v_path in version_paths[!(names(version_paths) %in% this_version)]) {
      files <- files[!startsWith(prefix = v_path, files)]
    }
  } else {
    files <- files[startsWith(prefix = version_paths[this_version], files)]
    files <- sub(version_paths[this_version], "", x = files, fixed = TRUE)
  }
  files
}

update_menus_for_subpath <- function(
  base_path,
  this_version,
  this_subpath,
  versions
) {
  version_path <- get_version_path(this_version, versions)
  path <- file.path(base_path, glue("{version_path}{this_subpath}"))
  lines <- get_lines(path)
  info <- get_info_lines(lines)
  if (!is.null(info)) {
    menu <- get_menu(this_version, this_subpath, versions)
    lines <- c(
      lines[seq_len(info$begin - 1L)],
      menu,
      lines[-seq_len(info$end)]
    )
    writeLines(lines, con = path)
  }
}

update_menus_for_version <- function(base_path, this_version, versions) {
  subpaths <- get_subpaths(base_path, this_version, versions)
  for (this_subpath in subpaths) {
    update_menus_for_subpath(base_path, this_version, this_subpath, versions)
  }
}

update_menus <- function(base_path, versions = get_versions(base_path)) {
  for (this_version in c(versions$current, versions$devel, versions$old)) {
    update_menus_for_version(base_path, this_version, versions)
  }
}

#config <- yaml::read_yaml("_pkgdown.yml")
#base_url <- config$url
#update_menus("docs")
