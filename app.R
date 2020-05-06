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
library(shinythemes)

ui <- navbarPage("rMSIflow",
                 theme = shinytheme("flatly"),
        tabPanel("Introduction", icon = icon("info-circle"), mainPanel(
          h1("Introduction", style = "font-size:24pt"),
          p(tags$b("rMSI"), "and", tags$b("rMSIproc"), "are R packages that allow the visualization of mass spectometry imaging data (rMSI)
            and the following data processing from all the spectra for each pixel to a peak matrix (rMSIproc).
            This tutorial wants to provide a flexible and customizable workflow for scientists analyzing different
            rMSIproc peak matrixes obtained from different images."),
          p("This tutorial will explain  how to:"),
          tags$ul(tags$li("Explore the peak matrix to visualize the different m/z intensities on the images"),
                  tags$li("What normalization steps should be utilized when making comparisons"),
                  tags$li("How to visualize the medium spectrum of the raw data"),
                  tags$li("The application of multidimensional reduction methods (PCA)"),
                  tags$li("The utilization of clustering methods like Kmeans and SOM to find clusters of
                           similar spectrums that can be related to differences between different tissue zones"),
                  tags$li("Visualization and comparison of those clusters")
          ),
        style = "font-size: 16pt")),
        tabPanel("Initial Parameters", icon = icon("sliders-h"), mainPanel(
            p(),
            sidebarLayout(sidebarPanel =sidebarPanel(selectInput("drive", "Select the root directory", choices = getVolumes()(), selected = getVolumes()()[1]),
                                                     shinyFilesButton(id = "peakMatrix", title = "Choose your peak matrix file: ", label= "Select Peak Matrix File",
                                                         multiple = F),
                                                     p(),
                                                     width = 6),
                          mainPanel = mainPanel(p("This is your selected file, please click the button to load it: "),
                                                verbatimTextOutput(outputId = "selectedfile", placeholder = T),
                                                actionButton("loadfile", label = "Load file!", icon = icon("bolt")),
                                                width = 6)
            ) ,
            hr(), width = 12
            )
        )
        # tabPanel("M/Z Images", sidebarLayout(sidebarPanel(sliderInput(inputId= "image", label = "Select an image", min = 1,
        #                                                               max = length(peakM$names), value = 1, step = 1),
        #                                                   sliderInput(inputId= "mz", label = "Select a m/z", min = 1,
        #                                                               max = length(peakM$mass), value = 1, step = 1),
        #                                                   radioButtons(inputId = "norm", label = "Normalizations", choiceNames = c("TIC", "RMS", "AcqTic", "No Normalization"),
        #                                                                choiceValues = c("tic", "rms", "acqtic", "nono") ,selected = "nono")),
        #                                                   mainPanel(plotlyOutput(outputId = "mzPlot"))))
        # tabPanel("PCA", sidebarLayout(sidebarPanel( )


                                      # ))
        )





# Define server logic required to draw a histogram
server <- function(input, output) {
    shinyFileChoose(input, "peakMatrix", roots= getVolumes(), filetypes = "RData")

    global <- reactiveValues(datapath = "No data selected", savedatapath = getwd())
    output$selectedfile <- renderText({
      paste(unlist(global$datapath), collapse = "\n")
    })
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                   input$peakMatrix
                   input$drive
                 },
                 handlerExpr = {
                   # browser()
                   global$datapath <-
                     strsplit(paste(unlist(input$peakMatrix[-length(input$peakMatrix)])[-1], collapse = "/"), split = "//")
                   global$datapath <- lapply(global$datapath, function(x){return(paste0(input$drive, x))})
                 })

    observeEvent(ignoreNULL = T,
                 eventExpr = {input$loadfile},
                 handlerExpr = {
                   tryCatch({
                     load(paste(unlist(global$datapath), collapse = "\n"))
                     peakm <- LoadPeakMatrix(paste(unlist(global$datapath), collapse = "\n"))
                     showNotification(paste("The file", paste(unlist(global$datapath), collapse = "\n"), "has been loaded!"), type = "message")
                     }, error = function(cond){
                     showNotification("Oops! There has been an error, please check the file address and type are valid", type = "error")}
                   )
                   }
                 )

    # output$mzPlot <- renderPlotly({
    #     normal <- switch(input$norm,
    #         tic = peakM$normalizations$TIC,
    #         rms = peakM$normalizations$RMS,
    #         acqtic = peakM$normalizations$AcqTic,
    #         nono = NA, nono
    #     )
    #
    #     initialPlotter(peakM,matrixList, peakM$mass[input$mz], input$image, normal, F, NA)
    #
    # })
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
