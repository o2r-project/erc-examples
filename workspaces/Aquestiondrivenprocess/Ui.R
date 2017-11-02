# ui.R

shinyUI(fluidPage(
  titlePanel("ReservoirOperations"),
  
  sidebarLayout(
    sidebarPanel(
      
      sliderInput("Kp", 
                  label = "Select hedging slope, Kp:",
                  min = 1.0, max = 3.0, value = 2.0, step = 0.1),
      
      selectInput("StreamSeq", 
                  label = "Select streamflow sequence:",
                  choices = list("Q1", "Q2", "Q3"),
                  selected = "Q1"),
      br(),
      p("This model explores a hypothesized feedback between water supply 
        reliability and per capita demand change on a surface  water 
        system with a reservoir. Users can explore the impact 
        of operating policy, either Standard Operating Policy (SOP) or Hedging 
        Policy (HP), by selecting the degree of hedging. Users can also test 
        the impact of multiple synthetic streamflow sequences. Please see the 
        accompanying paper, 'A Question Driven Socio-hydrological Modeling 
        Process', for model development details.")
      ),
    
    mainPanel(
      plotOutput(outputId = "Plot", width = "100%")
    )
      )
))