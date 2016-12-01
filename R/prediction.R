print.prediction = function(prd) {
  print(data.frame(mean = prd$mean,
                   sd = prd$sd, 
                   cv = prd$cv, 
                   q0.025 = prd$quantiles[1],
                   q0.5 = prd$quantiles[2],
                   q0.975 = prd$quantiles[3],
                   mce = prd$mce))
}
