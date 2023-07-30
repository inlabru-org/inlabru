# convert the gorillas dataset to sf format
# to use for testing

library(inlabru)
library(sf)

data(gorillas, package = "inlabru")
gorillas_sf <- list()

gorillas_sf$nests <- st_as_sf(gorillas$nests)
gorillas_sf$mesh <- gorillas$mesh
gorillas_sf$boundary <- st_as_sf(gorillas$boundary)
gorillas_sf$gcov <- lapply(gorillas$gcov, st_as_sf)
gorillas_sf$plotsample <- lapply(gorillas$plotsample, st_as_sf)

saveRDS(
  gorillas_sf,
  here::here("sf", "Data", "gorillas_sf.RDS")
)
