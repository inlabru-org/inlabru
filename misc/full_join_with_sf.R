library(tidyverse)
library(sf)

X <-
  data.frame(
    x = 1:3,
    y = 2:4,
    time = 101:103,
    w = 1:3,
    grp = c(1, 1, 2)
  )
Y <-
  data.frame(
    x = 11:14,
    y = 22:25,
    z = c("A", "B", "C", "A"),
    w = 1:4,
    grp = c(1, 2, 1, 2)
  )

full_join(X, Y, by = "grp")

X <- st_as_sf(X, coords = c("x", "y"))
Y <- st_as_sf(Y, coords = c("x", "y"))
st_geometry(Y) <- "something"
Y

Z <- full_join(as_tibble(X), as_tibble(Y), by = "grp")
st_as_sf(Z)

# This doesn't do what we want:
st_join(X, Y, by = "grp")
