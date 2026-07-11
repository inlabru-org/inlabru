# Load the package
library(plyr)
library(dplyr)
library(stringr)
devtools::load_all()

run_example <- function(file, pkg, run_donttest) {
  rd_file <- tools::parse_Rd(file)

  # 2. Extract the tags list
  tags <- sapply(rd_file, attr, "Rd_tag")

  # 3. Find the \examples section and convert it to text
  examples <- rd_file[tags == "\\examples"]
  if (length(examples) == 0) {
    return(NULL)
  }
  stopifnot(length(examples) == 1L)
  examples <- examples[[1]]

  dontrun <- sapply(examples, attr, "Rd_tag") == "\\dontrun"
  donttest <- sapply(examples, attr, "Rd_tag") == "\\donttest"
  examples[!donttest & !dontrun] <- lapply(
    examples[!donttest & !dontrun],
    function(x) {
      if (is.character(x)) {
        return(x)
      } else if (inherits(x, "Rd")) {
        return(paste0(unlist(x), collapse = ""))
      } else {
        return(as.character(x))
      }
    }
  )
  examples[donttest] <- lapply(examples[donttest], function(x) {
    if (!run_donttest) {
      return("\n")
    } else if (is.character(x)) {
      return(x)
    } else if (inherits(x, "Rd")) {
      return(paste0(unlist(x), collapse = ""))
    } else {
      return(as.character(x))
    }
  })
  examples[dontrun] <- lapply(examples[dontrun], function(x) {
    "\n"
  })

  code <- trimws(paste0(unlist(examples), collapse = ""))

  # Run example and collect runtime
  start_time <- Sys.time()
  tryCatch(
    expr = eval(parse(text = code, keep.source = FALSE)),
    error = function(e) print(paste("Error running example:", e))
  )
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Add result to data frame
  results <- data.frame(
    Package = pkg,
    File = basename(file),
    Example = code,
    Runtime = runtime,
    RunDontTest = run_donttest,
    stringsAsFactors = FALSE
  )
  return(results)
}

# Define a function to run examples and collect runtime
run_examples <- function(pkg = NULL, run_donttest = NULL) {
  # Get all Rd files
  if (is.null(pkg)) {
    pkg <- desc::desc(file = "DESCRIPTION")$get("Package")
    path <- file.path(getwd(), "man")
  } else {
    path <- file.path(find.package(pkg), "man")
  }
  rd_files <- dir(path, pattern = "*.Rd", full.names = TRUE)

  # Initialize data frame to store results
  results <- list()

  # Iterate over Rd files
  k <- 0
  for (file in rd_files) {
    if (is.null(run_donttest) || isTRUE(run_donttest)) {
      k <- k + 1
      results[[k]] <- run_example(file, pkg, TRUE)
    }
    if (is.null(run_donttest) || isFALSE(run_donttest)) {
      k <- k + 1
      results[[k]] <- run_example(file, pkg, FALSE)
    }
  }
  results <- do.call(rbind, results)

  rownames(results) <- NULL
  return(results)
}

# Example usage
results <- run_examples()

results <- results |>
  group_by(RunDontTest) |>
  arrange(desc(Runtime)) |>
  mutate(Rank = seq_len(n()), RelRank = Rank / n()) |>
  mutate(
    RunFrac = Runtime / sum(Runtime),
    CumRuntime = cumsum(Runtime),
    RelCumRuntime = CumRuntime / sum(Runtime)
  ) |>
  ungroup()
results |> dev
dplyr::select(File, RunDontTest, Runtime, Rank, RunFrac) |>
  arrange(desc(Runtime)) |>
  head(10)

results |>
  group_by(RunDontTest) |>
  summarise(TotalRuntime = sum(Runtime), .groups = "drop") |>
  arrange(desc(TotalRuntime))

library(ggplot2)
ggplot(results) +
  geom_point(aes(RelRank * 100, CumRuntime)) +
  coord_cartesian(ylim = c(0, max(results$CumRuntime))) +
  labs(x = "Percent of examples", y = "Cumulative runtime (seconds)") +
  geom_vline(xintercept = 10) +
  facet_wrap(~RunDontTest, ncol = 1)

ggplot(results) +
  geom_point(aes(RelRank * 100, Runtime, color = as.factor(RunDontTest))) +
  coord_cartesian(ylim = c(1e-5, max(results$Runtime))) +
  labs(x = "Percent of examples", y = "Runtime (seconds)") +
  scale_y_log10()
