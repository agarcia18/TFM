server <- function(input, output){
  
    dataframe<-reactive({
      if (is.null(input$bedfile))
        return(NULL)                
     circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
     names(circ) <- c("chr","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
     
     circ
     
    })
    
    output$ui_total <- renderUI({
      req(input$bedfile)
      circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
      
      valueBox(value=circ %>% count(),tags$b("TOTAL DNA CIRCLES"), color="maroon",
                   icon = tags$i(class = "fas fa-dna", style="font-size: 50px")
          )
      })
    
    output$table <- renderUI({
      req(input$bedfile)
      circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
      
      box(title = "Summary", width = NULL, solidHeader = TRUE,
          tableOutput("head"))
    })
    
    
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
    
  
  output$head <- renderTable({
    head(dataframe())
  })
}
