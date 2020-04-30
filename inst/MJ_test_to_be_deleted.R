
library(inlabru)

rm(list=ls())
# Defining test variables
region.coords = matrix(c(0, 0,
                         0, 5,
                         5, 5,
                         5, 0,
                         0, 0), ncol = 2,byrow=T)
region.points = sp::SpatialPoints(coords=region.coords)
region.polygon = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords=region.coords)),'0')))
mesh = INLA::inla.mesh.2d(loc.domain = region.points,max.edge = 1,offset = 1)

nsub_400_options = list(int.args=list(method="stable",nsub=19),
                        num.threads=1, blas.num.threads=1)
reproduce_options = list(num.threads=1, blas.num.threads=1)

# Default integration scheme
ips_default <- ipoints(region = region.polygon,domain = mesh)
ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_default, aes(size = weight))

# More points spread in each triangle as part of the stable integration method than default (nsub=9)
ips_nsub_400 <- ipoints(region = region.polygon,domain = mesh,int.args = nsub_400_options$int.args)
#ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_nsub20, aes(size = weight))

# Checking that integration points are identical
all.equal(coordinates(ips_default),coordinates(ips_nsub_400))

# Plotting the difference between the methods
ips_diff = ips_default
ips_diff$weight = NULL
ips_diff$absweightdiff = abs(ips_default$weight-ips_nsub_400$weight)
ips_diff$sign = as.factor(sign(ips_default$weight-ips_nsub_400$weight))

ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_diff, aes(size = absweightdiff,col=sign))

#### Now trying the lgcp function with the two methods ####

# Sample point pattern
set.seed(123)
spatstat::spatstat.options(npixel=100)
win <- spatstat::owin(range(region.coords[,1]),range(region.coords[,2]))
lg.s = spatstat::rpoispp(lambda=1,win=win)
obs = sp::SpatialPoints(cbind(x=lg.s$x,y=lg.s$y))

# Just use default matern prior
matern = INLA::inla.spde2.matern(mesh)

# Define model to fit
cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept

# Fit default model
set.seed(123)
start_default = proc.time()
mod_default = lgcp(cmp,data = obs,samplers = region.polygon,options=reproduce_options)
proc.time()-start_default
#user  system elapsed 
#1.953   0.353   2.306 


# Checking reproducability
set.seed(123)
mod_default_2 = lgcp(cmp,data = obs,samplers = region.polygon,options=reproduce_options)
all.equal(mod_default$summary.fixed,mod_default_2$summary.fixed)
#TRUE

# Fit model with nsub_400
set.seed(123)
start_400 = proc.time()
mod_nsub_400 = lgcp(cmp,data= obs, samplers = region.polygon,options=nsub_400_options)
proc.time()-start_400
#user  system elapsed 
#1.998   0.377   2.375 

# Fit with nsub_900
nsub_900_options = list(int.args=list(method="stable",nsub=29),
                        num.threads=1, blas.num.threads=1)
set.seed(123)
start_900 = proc.time()
mod_nsub_900 = lgcp(cmp,data= obs, samplers = region.polygon,options=nsub_900_options)
proc.time()-start_900
#user  system elapsed 
#2.385   0.473   2.856 

mod_default$summary.fixed
mod_nsub_400$summary.fixed
mod_nsub_900$summary.fixed

#> mod_default$summary.fixed
#mean        sd 0.025quant   0.5quant 0.975quant       mode         kld
#Intercept -0.1356498 0.2237671  -0.598403 -0.1270052  0.2777677 -0.1095957 2.06608e-05
#> mod_nsub_400$summary.fixed
#mean        sd 0.025quant   0.5quant 0.975quant       mode         kld
#Intercept -0.1357184 0.2237489 -0.5984426 -0.1270741  0.2776714 -0.1096651 2.05794e-05
#> mod_nsub_900$summary.fixed
#mean        sd 0.025quant   0.5quant 0.975quant       mode          kld
#Intercept -0.1357376 0.2237433  -0.598453 -0.1270934  0.2776439 -0.1096847 2.055263e-05


range(mod_default$summary.random$mySmooth$mean)
range(mod_nsub_400$summary.random$mySmooth$mean)
range(mod_nsub_900$summary.random$mySmooth$mean)

#> range(mod_default$summary.random$mySmooth$mean)
#[1] -0.01471897  0.03742867
#> range(mod_nsub_400$summary.random$mySmooth$mean)
#[1] -0.01471391  0.03736488
#> range(mod_nsub_900$summary.random$mySmooth$mean)
#[1] -0.01471252  0.03734367


range(mod_default$summary.random$mySmooth$sd)
range(mod_nsub_400$summary.random$mySmooth$sd)
range(mod_nsub_900$summary.random$mySmooth$sd)

#> range(mod_default$summary.random$mySmooth$sd)
#[1] 0.06702159 0.27098332
#> range(mod_nsub_400$summary.random$mySmooth$sd)
#[1] 0.0681257 0.2708596
#> range(mod_nsub_900$summary.random$mySmooth$sd)
#[1] 0.06157466 0.27083319

