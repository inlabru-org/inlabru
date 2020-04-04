
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

nsub_20_options = list(int.args=list(method="stable",nsub=19))


# Default integration scheme
ips_default <- ipoints(region = region.polygon,domain = mesh)
ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_default, aes(size = weight))

# More points spread in each triangle as part of the stable integration method than default (nsub=9)
ips_nsub20 <- ipoints(region = region.polygon,domain = mesh,int.args = nsub_20_options$int.args)
#ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_nsub20, aes(size = weight))

# Checking that integration poitns are identical
all.equal(coordinates(ips_default),coordinates(ips_nsub20))

# Plotting the difference between the methods
ips_diff = ips_default
ips_diff$weight = NULL
ips_diff$absweightdiff = abs(ips_default$weight-ips_nsub20$weight)
ips_diff$sign = as.factor(sign(ips_default$weight-ips_nsub20$weight))

ggplot()+gg(region.polygon)+gg(mesh)  + gg(ips_diff, aes(size = absweightdiff,col=sign))

#### Now trying the lgcp function with the two methods ####

# Sample point pattern
set.seed(123)
spatstat::spatstat.options(npixel=100)
win <- spatstat::owin(range(region.coords[,1]),range(region.coords[,2]))
lg.s = spatstat::rpoispp(lambda=1,win=win)
obs = sp::SpatialPoints(cbind(x=lg.s$x,y=lg.s$y))

#lg.s <- spatstat::rLGCP('matern', mu = 0,nu=1,win=win)
#Lam = attr(lg.s, 'Lambda')
#rf.s = log(Lam$v)
#fields::image.plot(rf.s)

# Just use default matern prior
matern = INLA::inla.spde2.matern(mesh)

# Define model to fit
cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept

# Fit default model
mod_default = lgcp(cmp,data = obs,samplers = region.polygon)

# Fit model with nsub_20
mod_nsub_20 = lgcp(cmp,data= obs, samplers = region.polygon,options=nsub_20_options)

mod_default$summary.fixed
mod_nsub_20$summary.fixed
#> mod_default$summary.fixed
#mean        sd 0.025quant   0.5quant 0.975quant       mode        kld
#Intercept -0.141376 0.4902942 -0.8431297 -0.1280812  0.4278482 -0.1093226 0.07765077
#> mod_nsub_20$summary.fixed
#mean        sd 0.025quant   0.5quant 0.975quant       mode        kld
#Intercept -0.1414562 0.4934652 -0.8422086 -0.1281776  0.4319441 -0.1094175 0.09025917


range(mod_default$summary.random$mySmooth$mean)
range(mod_nsub_20$summary.random$mySmooth$mean)
#> range(mod_default$summary.random$mySmooth$mean)
#[1] -0.05809279  0.08567435
#> range(mod_nsub_20$summary.random$mySmooth$mean)
#[1] -0.05841278  0.08915449


range(mod_default$summary.random$mySmooth$sd)
range(mod_nsub_20$summary.random$mySmooth$sd)
#> range(mod_default$summary.random$mySmooth$sd)
#[1] 0.2039991 0.7429675
#> range(mod_nsub_20$summary.random$mySmooth$sd)
#[1] 0.2011438 0.7526814


# Predict the spatial intensity surfaces
# increasing the number of samples to make sure differences are not due to randomness in sampling
lambda_default <- predict(mod_default, pixels(mesh), ~ exp(mySmooth + Intercept),n.samples = 10000)
lambda_nsub_20 <- predict(mod_nsub_20, pixels(mesh), ~ exp(mySmooth + Intercept),n.samples = 10000)

# Summaries of the mean fields. Just showing that these are different
summary(lambda_default$mean)
summary(lambda_nsub_20$mean)
#> summary(lambda_default$mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8642  0.8779  0.8868  0.8990  0.9134  0.9863 
#> summary(lambda_nsub_20$mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8701  0.8863  0.8937  0.9081  0.9182  1.0390 


# Visualizing how different these mean fields may be 
lambda_mean_diff = lambda_default
lambda_mean_diff$diffmean = lambda_default$mean - lambda_nsub_20$mean

lim = max(abs(range(lambda_mean_diff$diffmean)))

gg1 = ggplot()+  gg(lambda_default)+gg(region.polygon)+gg(mesh)+gg(obs) + ggtitle("Mean field default")
gg2 = ggplot()+  gg(lambda_nsub_20)+gg(region.polygon)+gg(mesh)+gg(obs) + ggtitle("Mean field nsub_20")
gg3 = ggplot()+  gg(lambda_mean_diff,aes(fill=diffmean))+gg(region.polygon)+gg(mesh)+gg(obs) +
  scale_fill_gradientn(colors = c("red","white","green"),limits = c(-lim,lim)) +
  ggtitle("Mean field difference (default - nsub_20)")

gridExtra::grid.arrange(gg1,gg2,gg3,ncol=2)
  