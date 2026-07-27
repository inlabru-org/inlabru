library(terra)
library(INLA)
library(inlabru)

# Set up template raster
testRast <- rast(
  res = 10000, # 10km raster
  xmin = 0,
  xmax = 700000,
  ymin = -100000,
  ymax = 1300000,
  crs = "EPSG:27700"
)

# Set up spatial covariate data spatRast with 2 'time points'
spCovar <- c(testRast, testRast)
spCovar[[1]] <- runif(ncell(testRast), 0, 10)
spCovar[[2]] <- runif(ncell(testRast), 0, 10)
names(spCovar) <- c("i1", "i2")

# Create 'where'
n <- 100
where <- spatSample(spCovar[[1]], n, as.points = TRUE)
names(where) <- "time"
timeVector <- sample(seq_len(2), n, replace = TRUE)
where[["time"]][] <- timeVector

### Testing

eval_spatial(spCovar, where = sf::st_as_sf(where), selector = "time")
eval_spatial(spCovar, where = sf::st_as_sf(where), layer = timeVector)

eval_spatial(spCovar, where = where, selector = "time")
eval_spatial(spCovar, where = where, layer = timeVector)

# From 2.13.0.9024: non-earth radius crs:es handled
eval_spatial(
  spCovar,
  where = fm_transform(sf::st_as_sf(where), fm_crs("sphere")),
  selector = "time"
)
eval_spatial(
  spCovar,
  where = fm_transform(sf::st_as_sf(where), fm_crs("sphere")),
  layer = timeVector
)

# Would need fmesher support (not yet in 0.5.0.9014), but terra 1.8-86 doesn't
# store points in +proj=geocent format, and silently drops one of the
# coordinates, so it would be limited to non-geocent crs specifications.
# But as _input_ from terra, it would be ok to transform between projections.
# Not OK for terra (unless it detects geocent and returns sf!):
eval_spatial(
  spCovar,
  where = fm_transform(where, fm_crs("sphere")),
  selector = "time"
)
eval_spatial(
  spCovar,
  where = fm_transform(where, fm_crs("sphere")),
  layer = timeVector
)
# OK for terra:
eval_spatial(
  spCovar,
  where = fm_transform(where, fm_crs("longlat_globe")),
  selector = "time"
)
eval_spatial(
  spCovar,
  where = fm_transform(where, fm_crs("longlat_globe")),
  layer = timeVector
)
