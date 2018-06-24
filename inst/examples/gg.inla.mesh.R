\donttest{
# Load Gorilla data
data("gorillas", package = "inlabru")

# Plot mesh using default edge colors

ggplot() + gg(gorillas$mesh)

# Don't show interior and exterior boundaries

ggplot() + gg(gorillas$mesh, interior = FALSE, exterior = FALSE)

# Change the edge colors

ggplot() + gg(gorillas$mesh, 
              edge.color = "green",
              int.color = "black",
              ext.color = "blue"
              )

# Use the x-coordinate of the vertices to colorize the triangles and
# mask the plotted area by the survey boundary, i.e. only plot the inside

xcoord = gorillas$mesh$loc[,1]
ggplot() + 
  gg(gorillas$mesh, color = (xcoord-580), mask = gorillas$boundary) +
  gg(gorillas$boundary)
}
