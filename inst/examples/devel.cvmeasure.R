\donttest{
if (require("INLA", quietly = TRUE)) {
  
# Load Gorilla data

data("gorillas", package = "inlabru")

# Use RColorBrewer 

library(RColorBrewer)

# Fit a model with two components:
# 1) A spatial smoothe SPDE 
# 2) A spatial covariate effect (vegetation)

pcmatern <- inla.spde2.pcmatern(gorillas$mesh, 
                                prior.sigma = c(0.1, 0.01), 
                                prior.range = c(5, 0.01))

cmp <- coordinates ~ vegetation(map = gorillas$gcov$vegetation, model = "factor") +
  spde(map = coordinates, model = pcmatern, mesh = gorillas$mesh) -
  Intercept

fit <- lgcp(cmp, gorillas$nests, samplers = gorillas$boundary)

# Predict SPDE and vegetation at the mesh vertex locations

vrt = vertices(gorillas$mesh)
joint <- predict(fit, vrt,  ~ spde + vegetation)
field <- predict(fit, vrt, ~ spde)
veg <- predict(fit, vrt, ~ vegetation)

# Plot component mean 

multiplot(ggplot() + gg(gorillas$mesh, color = joint$mean) + 
            coord_equal() + theme(legend.position= "bottom"),
          ggplot() + gg(gorillas$mesh, color = field$mean) + 
            coord_equal() + theme(legend.position= "bottom"),
          ggplot() + gg(gorillas$mesh, color = veg$mean) + 
            coord_equal() + theme(legend.position = "bottom"),
          cols = 3)

# Plot component variance 

multiplot(ggplot() + gg(gorillas$mesh, color = joint$var) + 
            coord_equal() + theme(legend.position = "bottom"),
          ggplot() + gg(gorillas$mesh, color = field$var) + 
            coord_equal() + theme(legend.position = "bottom"),
          ggplot() + gg(gorillas$mesh, color = veg$var) + 
            coord_equal() + theme(legend.position = "bottom"),
          cols = 3)

# Calculate variance and correlation measure

vm <- devel.cvmeasure(joint, field, veg)
lprange <- range(vm$var.joint, vm$var1,vm$var2)

# Variance contribution of the components

csc <- scale_fill_gradientn(colours = brewer.pal(9,"YlOrRd"), limits = lprange)
boundary <- gorillas$boundary

plot.1 <- ggplot() + gg(gorillas$mesh, color = vm$var.joint, mask = boundary) + 
  csc + coord_equal() + ggtitle("joint") + theme(legend.position = "bottom")
plot.2 <- ggplot() + gg(gorillas$mesh, color = vm$var1, mask = boundary) + 
  csc + coord_equal() + ggtitle("SPDE") + theme(legend.position = "bottom")
plot.3 <- ggplot() + gg(gorillas$mesh, color = vm$var2, mask = boundary) + 
  csc + coord_equal() + ggtitle("vegetation") + theme(legend.position = "bottom")

multiplot(plot.1, plot.2, plot.3, cols = 3)

# Covariance of SPDE field and vegetation

ggplot() + gg(gorillas$mesh, color = vm$cov) 

# Correlation between field and vegetation

ggplot() + gg(gorillas$mesh, color = vm$cor) 

# Variance and correlation integrated over space

vm.int <- devel.cvmeasure(joint, field, veg, 
                          samplers = ipoints(gorillas$boundary, gorillas$mesh),
                          mesh = gorillas$mesh)
vm.int

}
}
