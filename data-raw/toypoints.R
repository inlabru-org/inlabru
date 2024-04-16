# Generates toypoints.rda package data file
library(inlabru)
library(fmesher)
library(INLA)
library(sf)
library(ggplot2)

# Define domain and mesh

inner_boundary <-
  st_sf(
    geometry = list(st_sfc(
      st_polygon(
        list(
          matrix(
            c(
              pi, -pi,
              pi, pi,
              -pi, pi,
              -pi, -pi,
              pi, -pi
            ),
            byrow = TRUE,
            nrow = 5
          )
        )
      )
    ))
  )

outer_boundary <- st_buffer(inner_boundary, pi)

ggplot() +
  geom_sf(
    data = outer_boundary,
    fill = "blue",
    alpha = 0.1
  ) +
  geom_sf(
    data = inner_boundary,
    fill = "red",
    alpha = 0.1
  )

boundary <- fm_as_segm_list(list(inner_boundary, outer_boundary))
mesh <- fm_mesh_2d_inla(
  boundary = boundary,
  max.edge = c(0.5, 2),
  offset = c(pi, pi)
)

ggplot() +
  geom_sf(
    data = outer_boundary,
    fill = "blue",
    alpha = 0.1
  ) +
  geom_sf(
    data = inner_boundary,
    fill = "red",
    alpha = 0.1
  ) +
  geom_fm(data = mesh, alpha = 0)

# Simulate point data
set.seed(13927)
locs <- st_sf(geometry = st_sample(inner_boundary, 100))

ggplot() +
  geom_sf(
    data = outer_boundary,
    fill = "blue",
    alpha = 0.1
  ) +
  geom_sf(
    data = inner_boundary,
    fill = "red",
    alpha = 0.1
  ) +
  geom_fm(data = mesh, alpha = 0) +
  geom_sf(data = locs, colour = "green")

rho <- 3
sigma_matern <- 2
beta_0 <- 3
sigma_noise <- 0.2
set.seed(1320)
samp <- fm_matern_sample(mesh, rho = rho, sigma = sigma_matern)[, 1]
locs$z <- beta_0 + fm_evaluate(mesh, loc = locs, field = samp) +
  rnorm(n = nrow(locs),
        sd = sigma_noise)

ggplot() +
  geom_sf(
    data = inner_boundary,
    alpha = 0.1
  ) +
  geom_sf(
    data = locs,
    aes(colour = z)
  ) +
  scale_colour_viridis_c()

pred_locs <- st_sf(
  geometry = st_sample(inner_boundary,
    type = "regular",
    size = 100 * 100
  )
)

toypoints <- list(
  points = locs,
  mesh = mesh,
  boundary = inner_boundary,
  pred_locs = pred_locs
)

# usethis::use_data(toypoints, overwrite=TRUE)
