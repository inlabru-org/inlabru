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
    
    #' Prediction data
    preddata = data.frame(x = seq(xmin, xmax, length.out = 200))
    
    #' Load the intensity function shown in slides:
    data("lgcp1D")
    lambda = function(x) lambda_1D(x) * input$imult
    lambda = function(x) eval(parse(text = input$ifun)) * input$imult
    dflambda = data.frame(x = preddata$x, intensity = lambda(preddata$x), method = "lambda")
    
    #' Sample from intensity
    #set.seed(123)
    x = seq(xmin, xmax, length.out=100)
    y = lambda(x)
    mesh <- inla.mesh.1d(x, degree = 1)
    pts = sample.lgcp(mesh, log(lambda(x)))
    
    #' Histogram
    binwidth = breaks[2] - breaks[1]
    hst = hist(pts$x, breaks = breaks, plot = FALSE)
    hst = data.frame(x = hst$mid, count = hst$count, exposure = binwidth)
    
    
    #' GAM inference
    if (do.gam) {
      library(mgcv)
      fit.gam = gam(count ~ s(x), offset=log(exposure), data = hst, family = poisson())
      rgam = exp(predict(fit.gam, newdata = preddata))
      dfgam = data.frame(x = preddata$x, intensity = rgam, method = "gam")
    } else {dfgam = NULL}
    
    
    # bru inference
    if (do.bru) {
      
      mdl = ~ mySPDE(map=x, model = inla.spde2.matern(mesh), mesh = mesh) + Intercept -1
      prd = count ~ mySPDE + Intercept
      rbru = bru(hst, model = mdl, predictor = prd, E = hst$exposure, n = 1, family = "poisson")
      
      dbru = predict(rbru, ~ exp(mySPDE + Intercept), points = preddata)
      dfbru = data.frame(x = dbru$x, intensity = dbru$mean, method = "bru")
    } else {dfbru = NULL}
    
    # LGCP inference
    if (do.lgcp) {
      mdl = ~ mySPDE(map=x, model = inla.spde2.matern(mesh), mesh = mesh) + Intercept -1
      prd = x ~ mySPDE + Intercept
      rlgcp = lgcp(points = pts, model = mdl, predictor = prd, n = 1)
      dlgcp = predict(rlgcp, ~ exp(mySPDE + Intercept), points = preddata)
      dflgcp = data.frame(x = dlgcp$x, intensity = dlgcp$mean, method = "lgcp")
    } else {dflgcp = NULL}
    
    df = rbind(dflambda, dfgam, dfbru, dflgcp)
    ggplot() + 
      geom_line(data = df,  mapping = aes(x,y=intensity,color=method)) +
      geom_point(data=pts,aes(x=x), y=0.2,shape="|",cex=4,alpha=0.4) +
      geom_step(data=hst, aes(x=x-0.5*exposure,y=count/exposure), lwd=0.5, col="black", alpha = 0.5)
    
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
                  min = 20,
                  max = 100,
                  value = 26),
      sliderInput("imult",
                  "Intensity multiplier:",
                  min = 0,
                  max = 10,
                  value = 0.4,
                  step = 0.1),
      checkboxGroupInput("method", label = h3("Method"), 
                         choices = list("gam" = 1, "bru" = 2,"lgcp" = 3),
                         selected = c()),
      textInput("ifun","Intensity function", value = "lambda_1D(x)")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
)} # Close shinyApp and top level function



