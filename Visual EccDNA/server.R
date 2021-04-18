server <- function(input, output,session){
  # Change options to process bigger files (30MB)
  options(shiny.maxRequestSize=30*1024^2)
  
##### "Get started" tab server functions ####
    
  # Total number of circles ValueBox
    output$total <- renderUI({
      req(input$bedfile)
      circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
      
      valueBox(value=circ %>% count(),tags$b("TOTAL DNA CIRCLES"), color="navy",width=6,
               icon = tags$i(class = "fas fa-dna", style="color:white"))
      })
  # Get started button
    output$start<-renderUI({
      req(input$bedfile)
      box(status="primary", background = "light-blue",height=100,width=4,
          "Click this button to go to View Results:",br(),
          actionButton("click","Get Started!")
          )
      })
  # Quality plot
    output$quality_gg<-renderUI({
      req(input$bedfile)
      if (is.null(input$bedfile))
        return(NULL)                
      circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
      
      # Add column names
      names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
      
      # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
      circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
      
      box(collapsible=TRUE,title="Quality in DNA circles",height=250,renderPlot(
        ggplot(circ, aes(x = quality)) +  
          geom_bar(aes(y = (..count..)/sum(..count..),fill=quality),width=0.4,alpha=0.8,position = position_stack(reverse = TRUE))+
          scale_fill_manual(values=c("#FF5733","#FFC300","#DAF7A6","#6EBC63"))+
          scale_y_continuous(labels=scales::percent)+
          xlab("")+ylab("")+coord_flip()+
          theme_minimal()+
          guides(fill="none")
        ,height=200))
    })
  
  # Data table
    
    # Import data from file and add labels
    dataframe<-reactive({
      if (is.null(input$bedfile))
        return(NULL)                
      circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
      
      # Add column names
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
        scrollX=TRUE
      ))
    })
    
    # Render UI - Appear after file is uploaded
    output$table <- renderUI({
      req(input$bedfile)
      box(width = NULL, solidHeader = TRUE,
          dataTableOutput("data"))
    })
    
    
##### "Results" tab server functions ####  
    
  output$results <- renderPlot({
   
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
    
    # Add size of circle (number of bp)
    circ$size_bp <- circ$end - circ$start
    
    circ
    
      ggplot(circ, aes(x = chr, y = size_bp))+
        geom_point (aes(color=quality))+
        theme_minimal()
    })
}
