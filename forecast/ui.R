shinyUI(fluidPage(
  titlePanel("1-week-ahead forecasts by state space model"),
  
  sidebarLayout(
    sidebarPanel(h3("Parameters"), 
                 uiOutput("selectInput"),
                 uiOutput("sliderInput"),
                 width = 3),
    
    mainPanel(plotOutput("plot"), width = 9)
  )
))
