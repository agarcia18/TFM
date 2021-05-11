# Load required packages 
library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)

# SIZE FILTERING
circ_size <- 5000
range_small <- c(250,2500)
range_big <- c(10000,1000000)

# User interface
ui <- dashboardPage(
    
    dashboardHeader(title = "Visual EccDNA"),
    
    dashboardSidebar(
        sidebarMenu(
            id="tabs",
            menuItem("Get started",tabName="about",icon=icon("cloud-upload")),
            menuItem("View results", tabName ="results",icon=icon("bar-chart"), startExpanded = TRUE,
                menuSubItem(paste0("Circles under ", circ_size,"bp"), tabName = "smallcirc"),
                menuSubItem(paste0("Circles over ", circ_size,"bp"), tabName = "bigcirc")),
            menuItem("Circle Info",tabName="circle",icon=icon("circle-o")),
            menuItem("GitHub",tabName = "github",icon=icon("github"))
        )
        ),
    
    dashboardBody(
        tabItems(
            tabItem(
                tabName = "about",
                    fluidRow(HTML('<center><img src="eccdna.jpg"></center>'),
                             br(),br(),
                        box(
                        title = "Wellcome to Visual eccDNA2",width = 12, solidHeader = TRUE,
                        "This is an app for dynamic data visualization of eccDNA output files. Start selecting a BED file and this page will show an overview of the results. For better understanding of the output, circles will be separated by size. You can press the buttons View Results, and you will be able to see further information about the DNA circles contained in your file."
                    )),
                    fluidRow(
                            box(title= "Upload a Circle BED file with output:",fileInput("bedfile","Choose file:"), status="primary"),
                            uiOutput("total")                        
                            ),
                    fluidRow(uiOutput("quality_gg"),
                             uiOutput("click_small"),
                             uiOutput("click_big")
                             ),
                    fluidRow(uiOutput("table"))
            ),
            tabItem(
                tabName = "smallcirc",
                    fluidRow(uiOutput("plot_results_small"),
                             uiOutput("filters_small")),
                    ),
            tabItem(
                tabName = "bigcirc",
                fluidRow(uiOutput("plot_results_big"),
                         uiOutput("filters_big")),
            ),
            tabItem(
                tabName="circle",
                    fluidRow(uiOutput("selected_circle"))
                ),
            
            tabItem(
                tabName="github",
                tags$a(href="https://github.com/agarcia18/eccDNA","https://github.com/agarcia18/eccDNA")
            )
        )
    )
)

server <- function(input, output,session){
    # Change options to process bigger files (30MB)
    options(shiny.maxRequestSize=30*1024^2)
    
    
    # Import data from file and add labels
    dataframe<-reactive({
        if (is.null(input$bedfile))
            return(NULL)                
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        
        # Add circle id
        circ <- circ %>% mutate(circle_id = 1:n()) %>% select(circle_id, everything())
        
        # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        
        # Add size of circle (number of bp)
        circ$size_bp <- circ$end - circ$start
        
        circ
    })

    
    ##### "Get started" tab server functions ####
    
    # TOTAL CIRCLES VALUEBOX
    output$total <- renderUI({
        
        # Read data and count()
        req(input$bedfile)
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        
        valueBox(value=circ %>% count(),tags$b("TOTAL eccDNA CIRCLES"), color="navy",width=6,
                 icon = tags$i(class = "far fa-circle", style="color:white"))
    })
    # RESULTS BUTTON
    output$click_small<-renderUI({
        req(input$bedfile)
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        perc<-as.data.frame(circ %>% group_by(size_bp < circ_size) %>% count() %>% mutate(pct_tot = n/nrow(circ)*100))
        perc_small <- perc[2,3]
        
        box(background = "maroon",height=150,width=3,
            tags$b("Small eccDNA circles:"),
            tags$h2(paste0(round(perc_small,2),"%")),
            actionButton("click_smallcirc","View Results"))
    })
    
    output$click_big<-renderUI({
        req(input$bedfile)
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        perc<-as.data.frame(circ %>% group_by(size_bp < circ_size) %>% count() %>% mutate(pct_tot = n/nrow(circ)*100))
        perc_big <- perc[1,3]
        
        box(background = "purple",height=150,width=3,
            tags$b("Big eccDNA circles:"),
            tags$h2(paste0(round(perc_big,2),"%")),
            actionButton("click_bigcirc","View Results"))
    })
    

    # Enable going to each tab of results after click in Value Box
    observeEvent(input$click_smallcirc, {
        newtab <- switch(input$tabs,
                         "about" = "smallcirc")
        
        updateTabItems(session, "tabs", newtab)
    })
    observeEvent(input$click_bigcirc, {
        newtab2 <- switch(input$tabs,
                         "about" = "bigcirc")
        
        updateTabItems(session, "tabs", newtab2)
    })
    
    # QUALITY PLOT
    output$quality_gg<-renderUI({
        
        # Read data and add labels
        req(input$bedfile)
        if (is.null(input$bedfile))
            return(NULL)                
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        
        # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        
        # Make the box and the plot 
        box(title="Quality in eccDNA circles",height=300,width=6,status="primary",
            renderPlot(
            ggplot(circ, aes(x = quality)) +  
                geom_bar(aes(y = (..count..)/sum(..count..),fill=quality),width=0.4,alpha=0.8,position = position_stack(reverse = TRUE))+
                scale_fill_manual(values=c("#FF5733","#FFC300","#DAF7A6","#6EBC63"))+
                scale_y_continuous(labels=scales::percent)+
                xlab("")+ylab("")+coord_flip()+
                theme_minimal()+
                guides(fill="none")
            ,height=200),
            p(em("Quality is calculated from score: Bad (< 10), Low (10-50), Medium (50-200), Good (> 200).",style = "font-size:12px;")))
    })
    
             # DATA TABLE

        # Make the data table
        output$data <- DT::renderDataTable({
            DT::datatable(dataframe(),options = list(
                searching = FALSE,
                pageLength = 5,
                lengthMenu = c(5, 10, 15, 20),
                scrollX=TRUE,
                columnDefs = list(list(className = 'dt-center', targets = 0:4))
            ))
        })
    
    # Render UI - Appear after file is uploaded
    output$table <- renderUI({
        req(input$bedfile)
        box(width = NULL, solidHeader = TRUE,
            dataTableOutput("data"))
    })
    
    
    
    ##### "View Results" tab server functions - Smaller circles ####  
    
    # PLOT OUTPUT 
    
    # Import data and make it reactive
    data_circ_small <- reactive({
        if (is.null(input$bedfile))
            return(NULL) 
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        
        names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        
        # Set chromosome as factor and fix random chr
        circ$chr <- substr(circ$chr, start = 1, stop = 5)
        circ$chr<-str_remove(circ$chr,"_")
        circ$chr<-str_remove(circ$chr,"chr")
        circ$chr <- factor(circ$chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))
        
        
        # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        
        # Add circle id
        circ <- circ %>% mutate(circle_id = 1:n()) %>% select(circle_id, everything())
        
        # Add size of circle (number of bp)
        circ$size_bp <- circ$end - circ$start
        
        # Remove bad coverage circles and wrong discordant reads outputs
        circ <- filter(circ,coverage_cont < 0.5 & discordant_reads < size_bp)
        
        # Incorporate filters input (size and quality)
        data_circ_small<- filter(circ, between(size_bp, input$size_small[1], input$size_small[2]), quality == input$quality_small)
    })
    
    # Plot for small circles
    output$plot_results_small <- renderUI({
        req(input$bedfile)
        
        box(title=span(icon("fal fa-dna"), paste0("Circles under ", circ_size,"bp")),width=8, height = 480, background = "maroon",solidHeader = TRUE,
            
            renderPlotly({
                fig <- data_circ_small() %>%
                    plot_ly(type = 'scatter',
                            source="smallcircleSource",
                            mode = 'markers',
                            marker = list(
                                size = 8),
                            color = ~chr,
                            x = ~chr,
                            y = ~size_bp,
                            text = ~discordant_reads,
                            customdata=~split_reads,
                            hovertemplate = paste("<b>Discordant Reads: %{text}<br>",
                                                  "Split Reads: %{customdata}<br>",
                                                  "Size: %{y:.0} bp <br>"),
                            showlegend = FALSE
                    )
                fig <- fig  %>% layout(xaxis = list(
                    title = "Chromosome of origin"), yaxis = list(title="Size (number of base pairs)"))
                
                fig
            })
        )
    })
    
    
    # Box of widgets for extra filtering the plot
    output$filters_small<-renderUI({
        req(input$bedfile)
        box(title=span(icon("filter"),"Filters"),width=3,
            sliderInput(inputId = "size_small",
                        "By size (number of base pairs):",
                        min = 10,
                        max = circ_size,
                        value = range_small),br(),
            checkboxGroupInput(inputId = "quality_small",
                               "By quality:",
                               choices = c("Bad","Low","Medium","Good"),
                               selected= c("Medium","Good")),
        )
    })
    
    
    ##### "View Results" tab server functions - Bigger circles ####  
    
    # PLOT OUTPUT 
    
    # Import data and make it reactive
    data_circ_big <- reactive({
        if (is.null(input$bedfile))
            return(NULL) 
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        
        names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        
        # Set chromosome as factor and fix random chr
        circ$chr <- substr(circ$chr, start = 1, stop = 5)
        circ$chr<-str_remove(circ$chr,"_")
        circ$chr<-str_remove(circ$chr,"chr")
        circ$chr <- factor(circ$chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))
        
        
        # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        
        # Add circle id
        circ <- circ %>% mutate(circle_id = 1:n()) %>% select(circle_id, everything())
        
        # Add size of circle (number of bp)
        circ$size_bp <- circ$end - circ$start
        
        # Remove bad coverage circles and wrong discordant reads outputs
        circ <- filter(circ,coverage_cont < 0.5 & discordant_reads < size_bp)
        
        # Incorporate filters input (size and quality)
        data_circ_big<- filter(circ, between(size_bp, input$size_big[1], input$size_big[2]), quality == input$quality_big)
    })
    
    # Plot for big circles
    output$plot_results_big <- renderUI({
        req(input$bedfile)
        
        box(title=span(icon("fal fa-dna"), paste0( "Circles over ", circ_size,"bp")),width=8, height = 480, background = "purple",solidHeader = TRUE,
            
            renderPlotly({
                fig <- data_circ_big() %>%
                    plot_ly(type = 'scatter',
                            source="bigcircleSource",
                            mode = 'markers',
                            marker = list(
                                size = 22),
                            color = ~chr,
                            x = ~chr,
                            y = ~size_bp,
                            text = ~discordant_reads,
                            customdata=~split_reads,
                            hovertemplate = paste("<b>Discordant Reads: %{text}<br>",
                                                  "Split Reads: %{customdata}<br>",
                                                  "Size: %{y:.0} bp <br>"),
                            showlegend = FALSE
                    )
                fig <- fig  %>% layout(xaxis = list(
                    title = "Chromosome of origin"), yaxis = list(title="Size (number of base pairs)"))
                
                fig
            })
        )
    })
    
    
    # Box of widgets for extra filtering the plot
    output$filters_big<-renderUI({
        req(input$bedfile)
        box(title=span(icon("filter"), "Filters"),width=3,
            sliderInput(inputId = "size_big",
                        "By size (number of base pairs):",
                        min = circ_size,
                        max = 2000000,
                        value = range_big),
            checkboxGroupInput(inputId = "quality_big",
                               "By quality:",
                               choices = c("Bad","Low","Medium","Good"),
                               selected= c("Medium","Good")),
        )
    })
    
    
    
    
    ##### "Circle info" tab server functions ####  
    
    observeEvent(event_data("plotly_click",source="smallcircleSource"), {
        circletab <- switch(input$tabs,
                            "smallcirc" = "circle")
        
        updateTabItems(session, "tabs", circletab)
    })
    
    clickData <- reactive({
        currentEventData <- unlist(event_data(event = "plotly_click", source = "smallcircleSource", priority = "event"))
    })
    
    output$clickDataOut <- renderText({
        paste("Click data:", paste(names(clickData()), unlist(clickData()), sep = ": ", collapse = " | "))
    })
    
    observeEvent(event_data("plotly_click",source="smallcircleSource"), {
        output$selected_circle<- renderUI({
            box(title = "Circle", status="primary",width = 12, solidHeader = TRUE,
                htmlOutput("clickDataOut"))
            
        })
    })
    
    observeEvent(event_data("plotly_click",source="bigcircleSource"), {
        circletab2 <- switch(input$tabs,
                            "bigcirc" = "circle")
        
        updateTabItems(session, "tabs", circletab2)
    })
    
    clickDataBig <- reactive({
        currentEventData <- unlist(event_data(event = "plotly_click", source = "bigcircleSource", priority = "event"))
    })
    
    output$clickDataOutBig <- renderText({
        paste("Click data:", paste(names(clickDataBig()), unlist(clickDataBig()), sep = ": ", collapse = " | "))
    })
    
    observeEvent(event_data("plotly_click",source="bigcircleSource"), {
        output$selected_circle<- renderUI({
            box(title = "Circle", status="primary",width = 12, solidHeader = TRUE,
                htmlOutput("clickDataOutBig"))
            
        })
    })
    
}

shinyApp(ui=ui, server = server)
