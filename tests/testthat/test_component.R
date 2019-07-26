context("Component construction (test_component.R)")
test_that("Component construction: linear model", {

  # Only check the following if the new backend is present
  if (!exists("component")) {
    skip("Test non-functional using old backend.")
  } else {
    cmp <- component(~ beta(main = x, model = "linear", values = 1))[[1]]

    expect_equal(cmp$label, "beta")
    expect_equal(cmp$main$model, "linear")
    expect_equal(as.character(cmp$main$input$input), "x")

    # Covariate mapping
    df <- data.frame(x = 1:10)
    expect_equal(input_eval(cmp, part = "main", data = df), 1:10)

    # TODO: The following must be preceeded by an auto-model initialisation based on
    
    # Index generator
    # likelihood data.
    idx <- index_eval(cmp)
    expect_is(idx, "data.frame")
    expect_equal(colnames(idx)[1], "beta")
    expect_equal(colnames(idx)[2], "beta.group")
    expect_equal(colnames(idx)[3], "beta.repl")
    expect_equal(idx$beta, 1)

    # A-matrix
    A <- amatrix_eval(cmp, data = df)
    expect_is(A, "ddiMatrix")
    expect_equal(nrow(A), 10)
    expect_equal(ncol(A), 1)
    expect_equal(as.vector(A), 1:10)

    # Value
    v <- value(cmp, data = df, state = 2)
    expect_equal(v, 2 * df$x)

    v <- value(cmp, data = df, state = list(beta = 2))
    expect_equal(v, 2 * df$x)

    v <- value(cmp, data = df, state = list(beta = 2), A = A)
    expect_equal(v, 2 * df$x)
  }
})



test_that("Component construction: default mesh/mapping construction", {

  # Only check the following if the new backend is present
  if (!exists("make.default.mesh")) {
    skip("Test non-functional using old backend.")
  } else {
    cmp <- list(label = "testlabel")
    class(cmp) <- c("component", "list")

    # Check for failure/success on valid/invalid inputs
    expect_error(make.default.mesh(cmp,
      model = NULL, model.type = "iid",
      fvals = NULL
    ))
    expect_error(make.default.mesh(cmp,
      model = NULL, model.type = "iid",
      fvals = list()
    ))
    expect_error(
      make.default.mesh(cmp,
        model = NULL, model.type = "iid",
        fvals = list(n = 2)
      ),
      NA
    )

    expect_error(make.default.mesh(cmp,
      model = NULL, model.type = "seasonal",
      fvals = NULL
    ))
    expect_error(make.default.mesh(cmp,
      model = NULL, model.type = "seasonal",
      fvals = list()
    ))
    expect_error(
      make.default.mesh(cmp,
        model = NULL, model.type = "seasonal",
        fvals = list(season.length = 2)
      ),
      NA
    )

    for (model.type in c("rw1", "rw2", "ar", "ar1", "ou")) {
      expect_error(make.default.mesh(cmp,
        model = NULL, model.type = model.type,
        fvals = NULL
      ))
      expect_error(make.default.mesh(cmp,
        model = NULL, model.type = model.type,
        fvals = list()
      ))
      expect_error(
        make.default.mesh(cmp,
          model = NULL, model.type = model.type,
          fvals = list(values = 1:2)
        ),
        NA
      )
    }
  }
})
