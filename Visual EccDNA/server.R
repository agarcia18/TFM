server <- function(input, output){
  
  total_circles<-circ %>% count()
  p_bad <-((circ %>% filter (quality=="Bad") %>% count())/total_circles)*100
  p_low <- ((circ %>% filter(quality=="Low") %>% count())/total_circles)*100
  p_medium <- ((circ %>% filter(quality=="Medium") %>% count())/total_circles)*100
  p_good <- ((circ %>% filter(quality=="Good") %>% count())/total_circles)*100
  
  filter_circ <- reactive(circ %>% filter(quality %in% input$quality))
  filter_circ <- reactive(circ[circ$size_bp>=input$size[1]&circ$size_bp<=input$size[2],])
  
  output$results <- renderPlot({
    ggplot(filter_circ(), aes( x = size_bp, y = coverage_cont))+
      geom_point(colour="#E84E84")+
      geom_point(colour="white",size=0.6)+
      theme_bw() + theme(legend.position =c(0.15,0.85),legend.title = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  })
  
}
