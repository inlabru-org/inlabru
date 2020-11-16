if (requireNamespace("INLA", quietly = TRUE)) {
  INLA::inla.setOption(num.threads = testthat_inla_num_threads)
}
