local_bru_testthat_setup()

test_that("Discrete integration", {
  domain <- 2:5
  samplers <- 3:7
  ips_ <- data.frame(x = 3:5, weight = rep(1, 3), .block = 1L:3L)

  ips <- fm_int(domain, samplers = samplers)
  expect_identical(ips, ips_)

  domain <- as.character(domain)
  samplers <- as.character(samplers)
  ips <- fm_int(domain, samplers)
  ips_$x <- as.character(ips_$x)
  expect_identical(ips, ips_)

  domain <- factor(domain, levels = domain)
  samplers <- factor(samplers, levels = samplers)
  ips <- fm_int(domain, samplers)
  ips_$x <- factor(ips_$x, levels = domain)
  expect_identical(ips, ips_)
})


test_that("Continuous integration", {
  local_bru_safe_inla()

  domain <- INLA::inla.mesh.1d(2:5)

  samplers <- c(3, 7)
  ips_ <- data.frame(
    x = c(3:5, 3.5, 4.5),
    weight = c(1 / 6, 1 / 3, 1 / 6, 2 / 3, 2 / 3),
    .block = 1L
  )

  ips <- fm_int(domain, samplers = samplers)
  expect_identical(ips, ips_)

  # Check blockwise integration
  samplers <- rbind(c(3, 7), c(0, 10))
  ips <- fm_int(domain, samplers = samplers)
  sums <- as.vector(by(ips$weight, ips$.block, sum))
  expect_equal(sums, c(2, 3))
})


test_that("Tensor space integration", {
  local_bru_safe_inla()

  mesh_time <- INLA::inla.mesh.1d(1:5)
  mesh_space <- INLA::inla.mesh.1d(c(0, 5, 10))
  domain <- list(space = mesh_space, time = mesh_time)
  samplers1 <- tibble::tibble(
    time = rbind(c(1, 3), c(2, 4), c(3, 5)),
    space = rbind(c(0, 10), c(0, 5), c(5, 10)),
    weight = c(1, 10, 100)
  )

  ips1 <- fm_int(domain, samplers1)

  samplers2 <- tibble::tibble(
    space = samplers1$space,
    time = samplers1$time,
    weight = samplers1$weight
  )

  ips2 <- fm_int(domain, samplers2)

  expect_equal(sort(names(ips1)), sort(names(ips2)))
  expect_equal(
    dplyr::arrange(ips1, .block, time, space),
    dplyr::arrange(ips2[names(ips1)], .block, time, space)
  )
})
