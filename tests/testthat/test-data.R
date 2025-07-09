# Test for data input with mismatching input/output sizes (e.g. for regional
# integration models, etc)

test_that("Component construction: default mesh/mapping, data is list", {
  skip_on_cran()
  local_bru_safe_inla()

  lik <- bru_obs("gaussian",
    formula = y ~ .,
    data = list(x = c(1, 1.5, 2, NA, 3, 4), y = 11:15),
    used = bru_used(effect = "effect"),
    allow_combine = TRUE
  )

  cmp1 <- bru_comp_list(~ effect(c(1, 1.5, 2, 3, 4), model = "iid") - 1)
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(unique(lik$data$x), na.last = NA)
  )

  cmp1 <- bru_comp_list(~ effect(x, model = "rw2") - 1)
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(unique(lik$data$x), na.last = NA)
  )

  mesh1 <- fm_mesh_1d(
    sort(unique(lik$data$x), na.last = NA)
  )
  expect_error(
    bru_comp_list(
      ~ effect(x, model = "rw2", mapper = mesh1) - 1
    ),
    regexp = "Unknown mapper"
  )

  cmp1 <- bru_comp_list(
    ~ effect(x,
      model = "rw2",
      mapper = bru_mapper(mesh1, indexed = FALSE)
    ) - 1
  )
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    sort(unique(lik$data$x), na.last = NA)
  )

  cmp1 <- bru_comp_list(
    ~ effect(x,
      model = "rw2",
      mapper = bru_mapper(mesh1, indexed = TRUE)
    ) - 1
  )
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    seq_along(sort(unique(lik$data$x), na.last = NA))
  )
})



test_that("Component construction: unsafe intercepts, data is list", {
  local_bru_safe_inla()

  cmp <- bru_comp_list(~ something_unknown - 1)
  lik <- bru_obs(
    formula = response ~ ., data = list(response = 1:5),
    allow_combine = TRUE
  )
  lik <- bru_used_update(lik, labels = names(cmp))
  expect_error(
    object = {
      model <- bru_model(cmp, bru_obs_list(list(lik)))
    },
    paste0(
      "The input evaluation 'something_unknown' for 'something_unknown' ",
      "failed. Perhaps the data object doesn't contain the needed variables?"
    )
  )
})



test_that("Component construction: separate response_data input", {
  skip_on_cran()
  local_bru_safe_inla()

  lik1 <- bru_obs("gaussian",
    formula = y ~ c(sum(effect), sum(effect^2)),
    data = list(x = c(1, 1.5, 2, 3, 4)),
    response_data = data.frame(y = c(11, 12)),
    allow_combine = TRUE
  )

  lik2 <- bru_obs("gaussian",
    formula = y ~ c(sum(effect), sum(effect^2)),
    data = data.frame(x = c(1, 1.5, 2, 3, 4)),
    response_data = data.frame(y = c(11, 12)),
    allow_combine = TRUE
  )


  cmp1 <- bru_comp_list(~ effect(x, model = "iid") - 1)

  fit1 <- bru(components = cmp1, lik1)
  fit2 <- bru(components = cmp1, lik2)


  cmp1 <- bru_comp_list(~ effect(x, model = "rw2") - 1)
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik2)))
  expect_equal(ibm_values(cmp2$effect$mapper, multi = 1)$main, lik2$data$x)

  mesh1 <- fm_mesh_1d(lik1$data$x)
  expect_error(
    bru_comp_list(
      ~ effect(x, model = "rw2", mapper = mesh1) - 1
    ),
    regexp = "Unknown mapper"
  )

  cmp1 <- bru_comp_list(
    ~ effect(x,
      model = "rw2",
      mapper = bru_mapper(mesh1, indexed = FALSE)
    ) - 1
  )
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik1)))
  expect_equal(ibm_values(cmp2$effect$mapper, multi = 1)$main, lik1$data$x)

  cmp1 <- bru_comp_list(
    ~ effect(x,
      model = "rw2",
      mapper = bru_mapper(mesh1, indexed = TRUE)
    ) - 1
  )
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik1)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    seq_along(lik1$data$x)
  )
  cmp2 <- add_mappers(cmp1, lhoods = bru_obs_list(list(lik2)))
  expect_equal(
    ibm_values(cmp2$effect$mapper, multi = 1)$main,
    seq_along(lik2$data$x)
  )
})
