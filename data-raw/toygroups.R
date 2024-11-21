make.toygroups <- function() {
  set.seed(123)

  xx <- seq(-10, 10, length.out = 300)

  # Sample group locations
  mesh1D <- fm_mesh_1d(xx)
  g.lambda <- function(x) 60 * dnorm(0.2 * (x - 10)) + 20
  groups <- sample.lgcp(mesh1D, log(g.lambda(mesh1D$loc)))

  # Sample group size
  rate <- function(x) 1.5 * dnorm(0.1 * (x - 10))
  gsize <- function(x) {
    rexp(length(x), rate = rate(x))
  }
  groups$size <- gsize(groups$x)

  # Some data for easier plotting
  df.intensity <- data.frame(x = xx, g.lambda = g.lambda(xx))
  df.rate <- data.frame(x = xx, rate = rate(xx))
  sz <- seq(0, 50, length.out = 100)
  df.size <- data.frame(size = sz, dexp = dexp(sz, rate = 0.29))

  ##### Binned data
  # breaks = seq(-10,10,length.out = 16)
  # groups$bin = findInterval(groups$x, breaks)
  # hst = hist(groups$x, breaks = breaks, plot = FALSE)
  # hst = data.frame(
  #   bin = 1:(length(breaks)-1),
  #   nanimals = 0,
  #   ngroups = hst$counts
  # )
  # agg = aggregate(groups$size, by = list(bin = groups$bin), sum)
  # hst$nanimals[agg$bin] = agg$x
  # hst$x = breaks[1:(length(breaks)-1)]

  toygroups <- list(
    groups = groups,
    df.rate = df.rate,
    df.size = df.size,
    df.intensity = df.intensity
  )
  toygroups
}

# toygroups <- make.toygroups()
# usethis::use_data(toygroups, overwrite = TRUE)
