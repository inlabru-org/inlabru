  
  # Simulate observation distances from truncated half normal
  
    ns = 200
    tr = 3 # truncation
    d = abs(rnorm(n=ns,mean=0,sd=1))
    distances = d[d<tr]
  
    
  # run half normal detection function estimation via INLA   
    
    source("dfun_halfnormal.R")
    dTrafo = function(x){x^2}
    iMod=dfun_halfnormal(data.frame(distance=distances),truncation=tr,dTrafo=dTrafo,nIntegrate=500)
    
    source("dfun_plot.R")
    plot(iMod)
  

  # Compare with half normal estimated by the "Distance" package
    
    # run Distance
    library(Distance)
    dMod = ds(data.frame(distance=d),truncation=tr,adjustment=NULL)
    
    # plot versus INLA result
    plot(dMod)

    
  # Sample values of the detection function for our data
    source("inla_essmod.R")
    nsmp = 20
    svalues = dfun_halfnormal.sample.value(iMod,nsmp)
    plot(iMod)
    par(new=TRUE)
    plot(distances,svalues[[1]],ylab="", xlab="",xaxt='n', yaxt='n',ylim=c(0,1),xlim=c(0,iMod$truncation))
    par(new=TRUE)
    for (i in 2:nsmp) {
      par(new=TRUE)
      plot(distances,svalues[[i]],ylab="", xlab="",xaxt='n', yaxt='n',ylim=c(0,1),xlim=c(0,iMod$truncation))
    }
    

  # Apply INLA half normal DF to whales data
    source("io_whales.R")
    whales = io_whales.pkgdata()
    hist(whales$distance)
    iwMod=dfun_halfnormal(whales)
    plot(iwMod)
