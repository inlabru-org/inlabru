\donttest{
# Load rgl library (needed due to a bug in sphereplot library)

library(rgl)

# Load pantropoical dolphin data

data("mexdolphin")

# Show the globe

globe()

# Add mesh, ship transects and dolphin sightings stored
# as inla.mesh, SpatialLines and SpatialPoints objects, respectively

glplot(mexdolphin$mesh)
glplot(mexdolphin$samplers)
glplot(mexdolphin$points)
}
