\donttest{
# Load Gorilla data

data("gorillas", package = "inlabru")

# Plot Gorilla elevation covariate provided as SpatialPixelsDataFrame. 
# The same syntax applies to SpatialGridDataFrame objects.

ggplot() + gg(gorillas$gcov$elevation)

# Add Gorilla survey boundary and nest sightings

ggplot() + 
  gg(gorillas$gcov$elevation) + 
  gg(gorillas$boundary) +
  gg(gorillas$nests)

# Load pantropical dolphin data

data("mexdolphin")

# Plot the pantropiical survey boundary, ship transects and dolphin sightings

ggplot() + 
  gg(mexdolphin$ppoly) + # survey boundary as SpatialPolygon
  gg(mexdolphin$samplers) + # ship transects as SpatialLines
  gg(mexdolphin$points)  # dolphin sightings as SpatialPoints

# Change color

ggplot() + 
  gg(mexdolphin$ppoly, color = "green") + # survey boundary as SpatialPolygon
  gg(mexdolphin$samplers, color = "red") + # ship transects as SpatialLines
  gg(mexdolphin$points, color = "blue")  # dolphin sightings as SpatialPoints


# Visualize data annotations: line width by segment number

names(mexdolphin$samplers) # 'seg' holds the segment number
ggplot() + gg(mexdolphin$samplers, aes(color = seg))

# Visualize data annotations: point size by dolphin group size

names(mexdolphin$points) # 'size' holds the group size
ggplot() + gg(mexdolphin$points, aes(size = size))
}
