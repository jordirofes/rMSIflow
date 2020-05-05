#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(plotly)
library(rMSIproc)
library(shiny)
library(shinyFiles)

ui <- fluidPage(
    tabsetPanel(id = "apptab", selected = "Initial Parameters",type = "tabs", 
                tabPanel("Initial Parameters", mainPanel(
                    textInput(inputId= "peakMatrix", label= "Peak Matrix Path"))),
                tabPanel("M/Z Images", sidebarLayout(sidebarPanel(sliderInput(inputId= "image", label = "Select an image", min = 1, 
                                                                              max = length(peakM$names), value = 1, step = 1),
                                                                  sliderInput(inputId= "mz", label = "Select a m/z", min = 1, 
                                                                              max = length(peakM$mass), value = 1, step = 1),
                                                                  radioButtons(inputId = "norm", label = "Normalizations", choiceNames = c("TIC", "RMS", "AcqTic", "No Normalization"), 
                                                                               choiceValues = c("tic", "rms", "acqtic", "nono") ,selected = "nono")),
                                                                  mainPanel(plotlyOutput(outputId = "mzPlot"))))
                # tabPanel("PCA", sidebarLayout(sidebarPanel( )
                                              
                                              
                                              # ))
                )
    )
    
    


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$mzPlot <- renderPlotly({
        normal <- switch(input$norm,
            tic = peakM$normalizations$TIC,
            rms = peakM$normalizations$RMS,
            acqtic = peakM$normalizations$AcqTic,
            nono = NA, nono
        )
        
        initialPlotter(peakM,matrixList, peakM$mass[input$mz], input$image, normal, F, NA)

    })
# 
#     output$distPlot <- renderPlot({
#         # generate bins based on input$bins from ui.R
#         x    <- faithful[, 2]
#         bins <- seq(min(x), max(x), length.out = input$bins + 1)
# 
#         # draw the histogram with the specified number of bins
#         hist(x, breaks = bins, col = 'darkgray', border = 'white')
#     })
}

# Run the application 
shinyApp(ui = ui, server = server)
