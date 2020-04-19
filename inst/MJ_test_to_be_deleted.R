
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

nsub_400_options = list(int.args=list(method="stable",nsub=19))


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
mod_default = lgcp(cmp,data = obs,samplers = region.polygon)
proc.time()-start_default
#user  system elapsed 
#5.052   5.866   1.665 


# Checking reproducability
set.seed(123)
mod_default_2 = lgcp(cmp,data = obs,samplers = region.polygon)
all.equal(mod_default$summary.fixed,mod_default_2$summary.fixed)

# Unfortunatly, I am not able to get reproducable results running lgcp, so the 
# below comparison does not make sense, but I include it anyway.
# IS there a way to get reprocucable results with lgcp? I can't find an appropriate
# seed argument which can be passed to lgcp() through bru.options or similar

# Fit model with nsub_20
set.seed(123)
start_400 = proc.time()
mod_nsub_400 = lgcp(cmp,data= obs, samplers = region.polygon,options=nsub_400_options)
proc.time()-start_400
#user  system elapsed 
#4.868   6.659   1.772 

nsub_900_options = list(int.args=list(method="stable",nsub=29))
set.seed(123)
start_900 = proc.time()
mod_nsub_900 = lgcp(cmp,data= obs, samplers = region.polygon,options=nsub_900_options)
proc.time()-start_900
#user  system elapsed 
#5.258   5.621   2.129 

mod_default$summary.fixed
mod_nsub_400$summary.fixed
#> mod_default$summary.fixed
#mean        sd 0.025quant   0.5quant 0.975quant       mode        kld
#Intercept -0.141376 0.4902942 -0.8431297 -0.1280812  0.4278482 -0.1093226 0.07765077
#> mod_nsub_400$summary.fixed
#mean        sd 0.025quant   0.5quant 0.975quant       mode        kld
#Intercept -0.1414562 0.4934652 -0.8422086 -0.1281776  0.4319441 -0.1094175 0.09025917


range(mod_default$summary.random$mySmooth$mean)
range(mod_nsub_400$summary.random$mySmooth$mean)
#> range(mod_default$summary.random$mySmooth$mean)
#[1] -0.05809279  0.08567435
#> range(mod_nsub_400$summary.random$mySmooth$mean)
#[1] -0.05841278  0.08915449


range(mod_default$summary.random$mySmooth$sd)
range(mod_nsub_400$summary.random$mySmooth$sd)
#> range(mod_default$summary.random$mySmooth$sd)
#[1] 0.2039991 0.7429675
#> range(mod_nsub_400$summary.random$mySmooth$sd)
#[1] 0.2011438 0.7526814

