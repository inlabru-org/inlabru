sf_obj1 <- sf::st_as_sf(data.frame(x=1:3, y=3:5), coords=c("x","y"))
sf_obj2 <- sf::st_as_sf(data.frame(x=4:6, y=6:8), coords=c("x","y"))

library(sp)
sp_obj1 <- data.frame(x=1:3, y=3:5)
sp::coordinates(sp_obj1) <- ~x+y
# sp_obj1 <- as(sf_obj1, 'Spatial')

sp_obj1 <- as(sf_obj1, 'Spatial')
sp_obj2 <- as(sf_obj2, 'Spatial')

class(sp_obj1)
sp_obj1

cprod(sp_obj1, sp_obj2)
# what does it originally do? check cprod 20230108
# Error message from original cprod in non-overlapping dimension
# Error in .local(x, y, ...) : non-unique matches detected

ips <- ipoints(c(0, 10), 50, name = "myDim")
ips <- ipoints(matrix(c(0, 3, 5, 10), nrow = 2, byrow = TRUE), 50)
mesh <- inla.mesh.1d(seq(0, 10, by = 1))
ips <- ipoints(mesh, name = "time")

ips <- cprod(
       ipoints(c(0, 10), 10, name = "x"),
       ipoints(c(1, 5), 10, name = "y")
       )

ips1 <- ipoints(rbind(c(0, 3), c(3, 8)), 17, name = "myDim")
ips2 <- ipoints(domain = c(1, 2, 4), name = "myDiscreteDim")
cprod(ips1, ips2)


data(gorillas, package = "inlabru")
ips1 <- ipoints(gorillas$boundary)
cprod(ips1,ips1)
ips2 <- ipoints(gorillas$boundary, domain = gorillas$mesh)
