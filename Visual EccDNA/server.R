server <- function(input, output,session){
  # Change options to process bigger files (30MB)
  options(shiny.maxRequestSize=30*1024^2)
  
##### "Get started" tab server functions ####
    
  # TOTAL CIRCLES VALUEBOX
    output$total <- renderUI({
      
      # Read data and count()
      req(input$bedfile)
      circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
      
      valueBox(value=circ %>% count(),tags$b("TOTAL DNA CIRCLES"), color="navy",width=6,
               icon = tags$i(class = "fas fa-dna", style="color:white"))
      })
    
  # RESULTS BUTTON
    output$start<-renderUI({
      req(input$bedfile)
      box(background = "light-blue",height=120,width=4,
          tags$b("Click this button to explore your file:"),br(),br(),
          actionButton("click","View Results")
          )
      })
    
    # Enable going to View Results tab after click
    observeEvent(input$click, {
      newtab <- switch(input$tabs,
                       "about" = "results")
      
      updateTabItems(session, "tabs", newtab)
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
      box(title="Quality in DNA circles",height=250,renderPlot(
        ggplot(circ, aes(x = quality)) +  
          geom_bar(aes(y = (..count..)/sum(..count..),fill=quality),width=0.4,alpha=0.8,position = position_stack(reverse = TRUE))+
          scale_fill_manual(values=c("#FF5733","#FFC300","#DAF7A6","#6EBC63"))+
          scale_y_continuous(labels=scales::percent)+
          xlab("")+ylab("")+coord_flip()+
          theme_minimal()+
          guides(fill="none")
        ,height=200))
    })
  
  # DATA TABLE
    
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
    

    
##### "View Results" tab server functions ####  
  
  # PLOT OUTPUT 
    
    # Import data and make it reactive
    plot_circ <- reactive({
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
      circ <- filter(circ,coverage_cont < 0.5 & discordant_reads < 5000)
      
      # Incorporate filters input (size and quality)
      plot_circ<- filter(circ, between(size_bp, input$size[1], input$size[2]), quality == input$quality)
    })
    
   # Plot for circles under 5000bp
    output$plot_results <- renderUI({
      req(input$click)
      
      box(title="Circles under 5000bp",status="primary", width=8,solidHeader = TRUE,
          renderGirafe({
               gg_results <-ggplot(plot_circ(), aes( x = chr, y = size_bp,color=chr))+
                geom_point_interactive (aes(tooltip=discordant_reads),size=4)+
                 xlab("Original Chromosome")+
                 ylab("Size")+
                 theme_minimal()+
                 guides(color="none")
                
              girafe(ggobj = gg_results)
    })
      )
    })
    
    # Box of widgets for extra filtering the plot
    output$filters<-renderUI({
      req(input$click)
      box(title="Filters (for plot only)", icon="filter",width=4,status="warning",
        sliderInput(inputId = "size",
                    "By size (number of base pairs):",
                    min = 50,
                    max= 5000,
                    value = c(250,2500)),
        checkboxGroupInput(inputId = "quality",
                           "By quality:",
                           choices = c("Bad","Low","Medium","Good"),
                           selected= c("Medium","Good")),
    )
    })
    
  # TABLE OUTPUT
    table_results<-reactive({
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
      
      # Remove bad coverage circles and outliers with wrong discordant reads outputs
      circ <- circ %>% 
        filter(coverage_cont < 0.5 & size_bp >5000)  %>%
        select (chr,start,end,discordant_reads,split_reads,score,quality,size_bp)
      
      # Display only some columns
      names(circ) <- c("Original Chromosome","Start","End","Discordant Reads","Split Reads","Score","Quality","Size")
      circ
      })
    
    # Make the data table for circles over 5000bp
    output$bigcircles <- DT::renderDataTable({
      DT::datatable(table_results(),options = list(
        searching = FALSE,
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        scrollX=TRUE,
        columnDefs = list(list(className = 'dt-center', targets = 0:4))
      ))
    })
    
    # Render UI - Appear after click on the first page
    output$table_bigcircles <- renderUI({
      req(input$click)
      box(title="Circles over 5000bp",status="danger",width = NULL, solidHeader = TRUE,
          dataTableOutput("bigcircles"))
    })
    
}
##### "Circle info" tab server functions ####  

