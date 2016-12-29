#############################################################################################
### Shiny app for 1D Poisson process fitting
#############################################################################################

shiny.practical.fit1d = function(){

library(shiny)  
  
shinyApp(server = shinyServer(function(input, output) {
  
  output$distPlot <- renderPlot({
    
    # input = list()
    # input$bins = 50
    # input$imult = 1
    
    xmin = 1
    xmax = 55
    breaks = seq(xmin, xmax, length.out = input$bins+1)
    
    do.gam = 1 %in% input$method
    do.bru = 2 %in% input$method
    do.lgcp = 3 %in% input$method
    do.brupcp = 4 %in% input$method
    do.lgcppcp = 5 %in% input$method
    
    show.qtls = 1 %in% input$show.qtls
    
   
    range0 = input$range0
    sigma0 = input$sigma0
    
    # Seed
    set.seed(input$seed)
    
    #' Prediction data
    preddata = data.frame(x = seq(xmin, xmax, length.out = 200))
    
    #' Load the intensity function shown in slides:
    data("lgcp1D")
    lambda = function(x) lambda_1D(x) * input$imult
    lambda = function(x) eval(parse(text = input$ifun)) * input$imult
    dflambda = data.frame(x = preddata$x,
                          uq = NA, lq = NA, median = NA, sd = NA,
                          mean = lambda(preddata$x), method = "lambda")
    
    #' Sample from intensity
    x = seq(xmin, xmax, length.out=100)
    y = lambda(x)
    smesh <- inla.mesh.1d(x, degree = 1, boundary = input$boundary)
    pts = sample.lgcp(smesh, log(lambda(x)))
    
    #' Histogram
    binwidth = breaks[2] - breaks[1]
    hst = hist(pts$x, breaks = breaks, plot = FALSE)
    hst = data.frame(x = hst$mid, count = hst$count, exposure = binwidth)
    
    #' Mesh for INLA methods
    mesh <- inla.mesh.1d(seq(xmin, xmax, length.out=input$n.mesh), degree = 1, boundary = input$boundary)
    
    #' GAM inference
    if (do.gam) {
      library(mgcv)
      fit.gam = gam(count ~ s(x), offset=log(exposure), data = hst, family = poisson())
      rgam = exp(predict(fit.gam, newdata = preddata))
      dfgam = data.frame(x = preddata$x, mean = rgam, method = "gam", uq = NA, lq = NA, median = NA, sd = NA)
    } else {dfgam = NULL}
    
    
    # bru inference
    if (do.bru) {
      mdl = ~ mySPDE(map=x, model = inla.spde2.matern(mesh), mesh = mesh) + Intercept -1
      prd = count ~ mySPDE + Intercept
      rbru = bru(hst, model = mdl, predictor = prd, E = hst$exposure, n = 1, family = "poisson")
      dbru = predict(rbru, ~ exp(mySPDE + Intercept), points = preddata)
      dfbru = data.frame(dbru, method = "bru")
    } else {dfbru = NULL}
    
    
    
    # LGCP inference
    if (do.lgcp) {
      mdl = ~ mySPDE(map=x, model = inla.spde2.matern(mesh), mesh = mesh) + Intercept -1
      prd = x ~ mySPDE + Intercept
      rlgcp = lgcp(points = pts, model = mdl, predictor = prd, n = 1)
      dlgcp = predict(rlgcp, ~ exp(mySPDE + Intercept), points = preddata)
      dflgcp = data.frame(dlgcp, method = "lgcp")
    } else {dflgcp = NULL}
    
    # bru with pc prior inference
    if (do.brupcp) {
      mdl = ~ mySPDE(map=x, model = inla.spde2.pcmatern(mesh, prior.range = c(range0, NA), prior.sigma = c(1,NA)), mesh = mesh) + Intercept -1
      prd = count ~ mySPDE + Intercept
      rbru = bru(hst, model = mdl, predictor = prd, E = hst$exposure, n = 1, family = "poisson")
      dbrupcp = predict(rbru, ~ exp(mySPDE + Intercept), points = preddata)
      dfbrupcp = data.frame(dbrupcp, method = "bru pcprior")
    } else {dfbrupcp = NULL}
    
    # LGCP with pc prior inference
    if (do.lgcppcp) {
      mdl = ~ mySPDE(map=x, model = inla.spde2.pcmatern(mesh, prior.range = c(range0, NA), prior.sigma = c(1,NA)), mesh = mesh) + Intercept -1
      prd = x ~ mySPDE + Intercept
      rlgcp = lgcp(points = pts, model = mdl, predictor = prd, n = 1)
      dlgcppcp = predict(rlgcp, ~ exp(mySPDE + Intercept), points = preddata)
      dflgcppcp = data.frame(dlgcppcp, method = "lgcp pcprior")
    } else {dflgcppcp = NULL}
    
    df = rbind(dflambda, dfgam, dfbru, dflgcp, dfbrupcp, dflgcppcp)
    min.y = min(0.9*c(df$mean, df$lq, hst$count/hst$exposure), na.rm = TRUE)
    max.y = max(1.1*c(df$mean, df$uq, hst$count/hst$exposure), na.rm = TRUE)
    
    ggp = ggplot() + 
      geom_line(data = df,  mapping = aes(x,y=mean,color=method)) +
      geom_point(data=pts,aes(x=x), y=min.y,shape="|",cex=4,alpha=0.4) +
      geom_step(data=hst, aes(x=x-0.5*exposure,y=count/exposure), lwd=0.5, col="black", alpha = 0.5) +
      ylim(min.y, max.y)
    
    if (show.qtls & (do.bru | do.lgcp | do.brupcp | do.lgcppcp)) {
      ggp = ggp + geom_ribbon(data = df, 
                              mapping = aes(x = x, ymin = lq, ymax = uq, color = method), 
                              alpha = 0.1)
    }
    
    ggp
    
  })
  
}),



## UI

ui = shinyUI(fluidPage(
  
  # Application title
  titlePanel("Poisson process fitting"),
  
  # Sidebar with a slider input for number of binsa
  sidebarLayout(
    sidebarPanel(
      sliderInput("bins",
                  "Number of bins:",
                  min = 5,
                  max = 100,
                  value = 26),
      sliderInput("imult",
                  "Intensity multiplier:",
                  min = 0,
                  max = 10,
                  value = 0.4,
                  step = 0.1),
      
      numericInput("seed","Seed", 1234, min = 1, max = 5000, step = 1),
      
      textInput("ifun","Lambda", value = "lambda_1D(x)"),
      
      checkboxGroupInput("method", label = "Method", 
                         choices = list("gam" = 1, "bru" = 2,"lgcp" = 3, "bru with fixed pc prior" = 4, "lgcp with fixed pc prior" = 5),
                         selected = c()),
      
      checkboxGroupInput("show.qtls", label = "", 
                         choices = list("Show quantiles" = 1),
                         selected = c()),
      
      
      sliderInput("range0",
                  "PC prior range_0",
                  min = 1,
                  max = 50,
                  value = 5),
      
      sliderInput("sigma0",
                  "PC prior sigma_0",
                  min = 0.1,
                  max = 10,
                  value = 1,
                  step = 0.1),
      
      selectInput("boundary", "bru/lgcp boundary",
                  c("free" = "free",
                    "neumann" = "neumann")),
      
      sliderInput("n.mesh",
                  "Mesh nodes:",
                  min = 10,
                  max = 100,
                  value = 50,
                  step = 1)
      
      
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
)} # Close shinyApp and top level function



