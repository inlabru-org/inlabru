\donttest{
  if (require("ggplot2", quietly = TRUE) &&
      requireNamespace("terra", quietly = TRUE) &&
      require("tidyterra", quietly = TRUE)) {
    # Load Gorilla data

    gorillas <- inlabru::gorillas_sf
    gorillas$gcov <- gorillas_sf_gcov()

    # Plot Gorilla elevation covariate provided as terra::rast.

    ggplot() +
      gg(gorillas$gcov$elevation)

    # Add Gorilla survey boundary and nest sightings

    ggplot() +
      gg(gorillas$gcov$elevation) +
      gg(gorillas$boundary, alpha = 0) +
      gg(gorillas$nests)

    # Load pantropical dolphin data

    mexdolphin <- inlabru::mexdolphin_sf

    # Plot the pantropical survey boundary, ship transects and dolphin sightings

    ggplot() +
      gg(mexdolphin$ppoly, alpha = 0.5) + # survey boundary
      gg(mexdolphin$samplers) + # ship transects
      gg(mexdolphin$points) # dolphin sightings

    # Change color

    ggplot() +
      gg(mexdolphin$ppoly, color = "green", alpha = 0.5) + # survey boundary
      gg(mexdolphin$samplers, color = "red") + # ship transects
      gg(mexdolphin$points, color = "blue") # dolphin sightings


    # Visualize data annotations: line width by segment number

    names(mexdolphin$samplers) # 'seg' holds the segment number
    ggplot() +
      gg(mexdolphin$samplers, aes(color = seg))

    # Visualize data annotations: point size by dolphin group size

    names(mexdolphin$points) # 'size' holds the group size
    ggplot() +
      gg(mexdolphin$points, aes(size = size))
  }
}
