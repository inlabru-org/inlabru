

rm(list=ls())
region.coords = matrix(c(0, 0,
                         0, 5,
                         5, 5,
                         5, 0,
                         0, 0), ncol = 2,byrow=T)
region.points = sp::SpatialPoints(coords=region.coords)
region.polygon = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords=region.coords)),'0')))
mesh = INLA::inla.mesh.2d(loc.domain = region.points,max.edge = 1,offset = 1)

ggplot()+gg(region.polygon)+gg(mesh)


ips_default <- ipoints(region = region.polygon,domain = mesh)

ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_default, aes(size = weight))

ips_nsub20 <- ipoints(region = region.polygon,domain = mesh,nsub = 20)
ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_nsub20, aes(size = weight))

all.equal(coordinates(ips_default),coordinates(ips_nsub20))

ips_diff = ips_default
ips_diff$weight = NULL
ips_diff$absweightdiff = abs(ips_default$weight-ips_nsub20$weight)
ips_diff$sign = as.factor(sign(ips_default$weight-ips_nsub20$weight))

ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_diff, aes(size = absweightdiff,col=sign))

ips <- ipoints(domain, 50)
plot(ips)

data(gorillas, package = "inlabru")
ips <- ipoints(region = gorillas$boundary,domain = gorillas$mesh)
ips_nsub20 <- ipoints(region = gorillas$boundary,domain = gorillas$mesh,nsub=20)

all.equal(ips,ips_nsub20)

ggplot() + gg(gorillas$boundary) + gg(gorillas$mesh)
  
  + gg(ips, aes(size = weight))


