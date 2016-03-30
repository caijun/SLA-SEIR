shinyUI(fluidPage(
  titlePanel("Forecasts by state space SEIR model"),
  
  sidebarLayout(
    sidebarPanel(h3("Parameters"), 
                 selectInput("sel.season", label = "Season", 
                             choices = list("2003/2004"), 
                             selected = "2003-04"),
                 selectInput("sel.step", label = "Steps ahead", 
                             choices = list("1" = 1), 
                             selected = 1),
                 sliderInput("slider", label = "Week number", 
                             min = 1, max = 36, value = 1, 
                             animate = TRUE, step = 1),
                 width = 3),
    mainPanel(plotOutput("plot"), width = 9)
  )
))
