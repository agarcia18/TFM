# Load required packages
library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(stringr)
library(plotly)
library(DT)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(org.Hs.eg.db)
library(RITANdata)
library(RITAN)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

launch_VisualEccDNA<-function(){
# SIZE FILTERING
circ_size <- 5000
range_small <- c(250,2500)
range_big <- c(10000,1000000)

######################################################### UI #########################################################

ui <- dashboardPage(

    dashboardHeader(title = "Visual eccDNA"),

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
                fluidRow(HTML('<center><img src="https://i.ibb.co/7tj8wsG/logo.png"></center>'),
                             br(),br(),
                        box(width = 12, solidHeader = TRUE,
                        "Wellcome to Visual eccDNA. This is an app for dynamic data visualization of Extrachromosomal Circular DNA (eccDNA) output files. Start selecting a BED file and this page will show an overview of the results. For better understanding of the output, circles will be separated by size. You can press the buttons View Results, and you will be able to see further information about the DNA circles contained in your file."
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
                    fluidRow(
                        uiOutput("plot_results_small"),
                             uiOutput("filters_small"),
                             uiOutput("backbuttonsmall"),

                        ),
                    fluidRow(uiOutput("table_small"),
                             uiOutput("pathway_small"))
                    ),
            tabItem(
                tabName = "bigcirc",
                fluidRow(uiOutput("plot_results_big"),
                         uiOutput("filters_big"),
                         uiOutput("backbuttonbig")
                         ),
                fluidRow(uiOutput("table_big"),
                         uiOutput("pathway_big"))
            ),
            tabItem(
                tabName="circle",
                    fluidRow(uiOutput("selected_circle"))
                ),

            tabItem(
                tabName="github",
                "You can follow the development and version control of this app in:",br(),
                tags$a(href="https://github.com/agarcia18/eccDNA","https://github.com/agarcia18/eccDNA")
            )
        )
    )
)

######################################################### SERVER #########################################################

server <- function(input, output,session){
    # Change options to process bigger files (30MB)
    options(shiny.maxRequestSize=30*1024^2)


    # Import data from file and add labels
    dataframe<-reactive({
        if (is.null(input$bedfile))
            return(NULL)
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")

        # Add circle id
        circ <- circ %>% mutate(circle_id = 1:n()) %>% dplyr::select(circle_id, everything())

        # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

        # Add size of circle (number of bp)
        circ$size_bp <- circ$end - circ$start

        circ })


#################################### "Get Started" tab server functions ####################################

    # TOTAL CIRCLES VALUEBOX
    output$total <- renderUI({
        req(input$bedfile)
        valueBox(value=dataframe() %>% count(),tags$b("TOTAL eccDNA CIRCLES"), color="navy",width=6,
                 icon = tags$i(class = "far fa-circle", style="color:white"))
    })


    # VIEW RESULTS BUTTONS - SMALL & BIG CIRCLES

    # Small circles % and click button
    output$click_small<-renderUI({
        req(input$bedfile)
        perc<-as.data.frame(dataframe() %>% group_by(size_bp < circ_size) %>% count() %>% mutate(pct_tot = n/nrow(dataframe())*100))
        perc_small <- perc[2,3]

        box(background = "maroon",height=150,width=3,
            h4("Small eccDNA circles:",align="center"),
            h1(paste0(round(perc_small,2),"%"),align="center"),
            actionButton("click_smallcirc","View Results"),align="center")
    })

    # Click event - Change to the small circles plot tab
    observeEvent(input$click_smallcirc, {
        newtab <- switch(input$tabs,
                         "about" = "smallcirc")

        updateTabItems(session, "tabs", newtab)
    })


   # Big circles % and click button
    output$click_big<-renderUI({
        req(input$bedfile)
        perc<-as.data.frame(dataframe() %>% group_by(size_bp < circ_size) %>% count() %>% mutate(pct_tot = n/nrow(dataframe())*100))
        perc_big <- perc[1,3]

        box(background = "purple",height=150,width=3,
            h4("Big eccDNA circles:",align="center"),
            h1(paste0(round(perc_big,2),"%"),align="center"),
            actionButton("click_bigcirc","View Results"),align="center")
    })
    # Click event - Change to the big circles plot tab
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
        box(title="Quality in eccDNA circles",height=300,width=6,status="primary",
            renderPlot(
            ggplot(dataframe(), aes(x = quality)) +
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
    output$table <- renderUI({
        req(input$bedfile)
        box(width = NULL, solidHeader = TRUE,
            DT::renderDataTable({
            DT::datatable(dataframe(),options = list(
                searching = FALSE,
                pageLength =8 ,
                lengthMenu = c(5, 10, 15, 20),
                scrollX=TRUE,
                columnDefs = list(list(className = 'dt-center', targets = 0:4))
            ))
        })
        )
    })

#################################### "View Results" tab server functions - Small circles ####################################

    # PLOT OUTPUT - SMALL CIRCLES

    # Import data for plot and make it reactive
    data_circ_small <- reactive({
        if (is.null(input$bedfile))
            return(NULL)
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)

        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")

        # Set chromosome as factor and fix random outputs
        circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
        circ$chrom<-str_remove(circ$chrom,"_")
        circ$chrom<-str_remove(circ$chrom,"chr")
        circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))


        # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

        # Add circle id
        circ <- circ %>% mutate(circle_id = 1:n()) %>% dplyr::select(circle_id, everything())

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
        box(title=span(icon("circle-o"), paste0("Circles under ", circ_size,"bp")),width=9, height = 460,solidHeader = TRUE,background="maroon",
            renderPlotly({
                fig <- data_circ_small() %>%
                    plot_ly(type = 'scatter',
                            source="smallcircleSource",
                            mode = 'markers',
                            marker = list(
                                size = 8),
                            color = ~chrom,
                            x = ~chrom,
                            y = ~size_bp,
                            text = ~discordant_reads,
                            customdata=~split_reads,
                            hovertemplate = paste("<b>Discordant Reads: %{text}<br>",
                                                  "Split Reads: %{customdata}<br>",
                                                  "Size: %{y:.0} bp <br>"),
                            showlegend = FALSE)
                fig <- fig  %>% layout(xaxis = list(
                    title = "Chromosome of origin"), yaxis = list(title="Size (number of base pairs)"))
                fig
            })
        )
    })


    # BOX OF FILTERS FOR THE PLOT - SMALL CIRCLES
    output$filters_small<-renderUI({
        req(input$bedfile)
        box(title=span(icon("filter"),"Filters (for plot)"),width=3,
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

    # "BACK TO GET STARTED" BUTTON

    # Render the button
    output$backbuttonsmall<-renderUI({
        req(input$bedfile)
        box(width=2,
        actionButton("click_backsmallcirc","Back to Get Started"),align="center")
    })

    # Click event - Change tab
    observeEvent(input$click_backsmallcirc, {
        back1 <- switch(input$tabs,
                         "smallcirc" = "about")

        updateTabItems(session, "tabs", back1)
    })

    # TABLE OF GENES IN SMALL CIRCLES
    output$table_small <- renderUI({
        req(input$bedfile)

        # Import the data
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

        # Process it to get gene list
        coords<- circ %>% filter (size_bp < circ_size & quality !="Bad") %>% dplyr::select(chrom,start,end) %>% makeGRangesFromDataFrame
        genes <-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
        genes_df <- as.data.frame(subsetByOverlaps(genes,coords))
        entrezid <- genes_df$gene_id
        hs <- org.Hs.eg.db
        genesymbol_df<- AnnotationDbi::select(hs,
                                              keys = entrezid,
                                              columns = c("ENTREZID", "SYMBOL"),
                                              keytype = "ENTREZID")
        names(genesymbol_df)[names(genesymbol_df) == "ENTREZID"] <- "gene_id"
        gene_list_small<-left_join(genes_df,genesymbol_df)
        names(gene_list_small)<-c("Chromosome","Start","End","Width","Strand","Gene ID", "Gene Symbol")

        # Render the data table
        box(title=span(icon("fal fa-dna"),"Genes in circles"), width = 7, solidHeader = TRUE,
            DT::renderDataTable({
                DT::datatable(gene_list_small[, names(gene_list_small) != "Gene ID"],options = list(
                    autoWidth = TRUE,
                    searching = FALSE,
                    pageLength =6 ,
                    lengthMenu = c(5, 10, 15, 20),
                    scrollX=TRUE,
                    columnDefs = list(list(className = 'dt-center', targets = 0:5))
                )) %>% DT::formatStyle(columns = colnames(.), fontSize = '8pt')
            })
        )
        })

    # PLOT - PATHWAY ENRICHMENT - SMALL CIRCLES
    output$pathway_small<-renderUI({
        req(input$bedfile)

        # Import the data
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

        # Get genes list
        coords<- circ %>% filter (size_bp < circ_size & quality!="Bad") %>% dplyr::select(chrom,start,end) %>% makeGRangesFromDataFrame
        genes <-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
        genes_df <- as.data.frame(subsetByOverlaps(genes,coords))
        entrezid <- genes_df$gene_id
        hs <- org.Hs.eg.db
        genesymbol_df<- AnnotationDbi::select(hs,
                                              keys = entrezid,
                                              columns = c("ENTREZID", "SYMBOL"),
                                              keytype = "ENTREZID")
        names(genesymbol_df)[names(genesymbol_df) == "ENTREZID"] <- "gene_id"
        gene_list_small<-left_join(genes_df,genesymbol_df)


        # Enrichment to find pathways
        enrich_small <- term_enrichment(gene_list_small$SYMBOL, resources = "ReactomePathways", all_symbols = cached_coding_genes)

        # Fix name of pathways for the plot
        enrich_small$name <-sapply(strsplit(enrich_small$name, split='.', fixed=TRUE), function(x) (x[2]))

        box(title=span(icon("microscope"),"Enriched pathways (Reactome)"),width=5,
            renderPlot(ggplot(enrich_small[1:5,],aes(x=q,y=name,fill=p))+
                           geom_col()+
                           scale_fill_gradient(low="#F7CFDA",high="#DA1853")+
                           xlab("-log10 (q-value)")+ylab("")+labs(fill = "p-value")+
                           theme_minimal()

                )
            )
    })

#################################### "View Results" tab server functions - Big circles #####################################

    # PLOT OUTPUT - BIG CIRCLES

    # Import data and make it reactive
    data_circ_big <- reactive({
        if (is.null(input$bedfile))
            return(NULL)
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)

        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")

        # Set chromosome as factor and fix random outputs
        circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
        circ$chrom<-str_remove(circ$chrom,"_")
        circ$chrom<-str_remove(circ$chrom,"chr")
        circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))


        # Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

        # Add circle id
        circ <- circ %>% mutate(circle_id = 1:n()) %>% dplyr::select(circle_id, everything())

        # Add size of circle (number of bp)
        circ$size_bp <- circ$end - circ$start

        # Remove bad coverage circles and wrong discordant reads outputs
        circ <- filter(circ,coverage_cont < 0.5 & discordant_reads < size_bp)

        # Incorporate filters input (size and quality)
        data_circ_big<- filter(circ, between(size_bp, input$size_big[1], input$size_big[2]), quality == input$quality_big)
    })

    # Plot of big circles
    output$plot_results_big <- renderUI({
        req(input$bedfile)

        box(title=span(icon("circle-o"), paste0( "Circles over ", circ_size,"bp")),width=9, height = 480, background = "purple",solidHeader = TRUE,

            renderPlotly({
                fig <- data_circ_big() %>%
                    plot_ly(type = 'scatter',
                            source="bigcircleSource",
                            mode = 'markers',
                            marker = list(
                                size = 22),
                            color = ~chrom,
                            x = ~chrom,
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


    # BOX OF FILTERS FOR THE PLOT - BIG CIRCLES
    output$filters_big<-renderUI({
        req(input$bedfile)
        box(title=span(icon("filter"), "Filters (for plot)"),width=3,
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

    # "BACK TO GET STARTED" BUTTON

    #Render the button
    output$backbuttonbig<-renderUI({
        req(input$bedfile)
        box(width=2,
            actionButton("click_backbigcirc","Back to Get Started"),align="center")
    })

    # Click event - Back to "Get Started"
    observeEvent(input$click_backbigcirc, {
        back2 <- switch(input$tabs,
                        "bigcirc" = "about")

        updateTabItems(session, "tabs", back2)
    })

    # TABLE OF GENES IN BIG CIRCLES
    output$table_big <- renderUI({
        req(input$bedfile)

        # Import the data
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

        # Process it to get gene list
        coords<- circ %>% filter (size_bp > circ_size & quality!="Bad") %>% dplyr::select(chrom,start,end) %>% makeGRangesFromDataFrame
        genes <-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
        genes_df <- as.data.frame(subsetByOverlaps(genes,coords))
        entrezid <- genes_df$gene_id
        hs <- org.Hs.eg.db
        genesymbol_df<- AnnotationDbi::select(hs,
                                              keys = entrezid,
                                              columns = c("ENTREZID", "SYMBOL"),
                                              keytype = "ENTREZID")
        names(genesymbol_df)[names(genesymbol_df) == "ENTREZID"] <- "gene_id"
        gene_list_big<-left_join(genes_df,genesymbol_df)
        names(gene_list_big)<-c("Chromosome","Start","End","Width","Strand","Gene ID", "Gene Symbol")

        # Render the data table
        box(title=span(icon("fal fa-dna"),"Genes in circles"), width = 7, solidHeader = TRUE,
            DT::renderDataTable({
                DT::datatable(gene_list_big[, names(gene_list_big) != "Gene ID"],options = list(
                    autoWidth = TRUE,
                    searching = FALSE,
                    pageLength =6 ,
                    lengthMenu = c(5, 10, 15, 20),
                    scrollX=TRUE,
                    columnDefs = list(list(className = 'dt-center', targets = 0:4))
                ))
            })
        )
    })

    # PLOT - PATHWAY ENRICHMENT
    output$pathway_big<-renderUI({
        req(input$bedfile)

        # Import the data
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

        # Get genes list
        coords<- circ %>% filter (size_bp > circ_size & quality!="Bad") %>% dplyr::select(chrom,start,end) %>% makeGRangesFromDataFrame
        genes <-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
        genes_df <- as.data.frame(subsetByOverlaps(genes,coords))
        entrezid <- genes_df$gene_id
        hs <- org.Hs.eg.db
        genesymbol_df<- AnnotationDbi::select(hs,
                                              keys = entrezid,
                                              columns = c("ENTREZID", "SYMBOL"),
                                              keytype = "ENTREZID")
        names(genesymbol_df)[names(genesymbol_df) == "ENTREZID"] <- "gene_id"
        gene_list_big<-left_join(genes_df,genesymbol_df)


        # Enrichment to find pathways
        enrich_big <- term_enrichment(gene_list_big$SYMBOL, resources = "ReactomePathways", all_symbols = cached_coding_genes)

        # Fix name of pathways for the plot
        enrich_big$name <-sapply(strsplit(enrich_big$name, split='.', fixed=TRUE), function(x) (x[2]))

        box(title=span(icon("microscope"),"Enriched pathways (Reactome)"),width=5,
            renderPlot(ggplot(enrich_big[1:5,],aes(x=q,y=name,fill=p))+
                           geom_col()+
                           scale_fill_gradient(low="#CEB9DF",high="#51119E")+
                           xlab("-log10 (q-value)")+ylab("")+labs(fill = "p-value")+
                           theme_minimal()

            )
        )
    })


####################################  "Circle info" tab server functions ####################################

# As we have two plots of output (small & big circles), this page will be reactive to both of them

    ## SOURCE: SMALL CIRCLES PLOT
    # Change tab and render info when clicking in the plot of small circles
    observeEvent(event_data("plotly_click",source="smallcircleSource"), {
        circletab1 <- switch(input$tabs,
                             "smallcirc" = "circle")

        updateTabItems(session, "tabs", circletab1)
    })

    observeEvent(event_data("plotly_click",source="smallcircleSource"), {
        output$selected_circle<- renderUI({
            fluidRow(
                box(title = "CIRCLE", status="primary",width = 12, solidHeader = TRUE,
                    tableOutput("TableDataOutSmall"),
                    plotOutput("CirclePlotSmall")),
                box(title=span(icon("fal fa-dna"),"GENOMIC DATA"), status="primary",width = 12, solidHeader = TRUE,
                    dataTableOutput("GenesSmall"))
            )

        })
    })


    # Get the data from the plot click
    clickDataSmall <- reactive({unlist(event_data(event = "plotly_click", source = "smallcircleSource", priority = "event"))
    })

    # SELECTED CIRCLE INFORMATION
    output$TableDataOutSmall <- renderTable({
        # Import the data file and process it
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
        circ$chrom<-str_remove(circ$chrom,"_")
        circ$chrom<-str_remove(circ$chrom,"chr")
        circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))

        # Use data from the plot to filter
        chrom_filter_small <- unlist(clickDataSmall()[3])
        size_filter_small <- unlist(clickDataSmall()[4])
        sr_filter_small <-unlist(clickDataSmall()[5])

        circ %>%
            dplyr::filter (chrom==chrom_filter_small & size_bp==size_filter_small & split_reads==sr_filter_small) %>%
            dplyr::select(chrom,start,end,discordant_reads,split_reads,score,coverage_mean,size_bp) %>%
            dplyr::rename(Chromosome=chrom,Start=start,End=end,DiscordantReads=discordant_reads,SplitReads=split_reads,Score=score,Coverage=coverage_mean,Size=size_bp)

    })


    # CIRCULAR PLOT
    output$CirclePlotSmall <-renderPlot({
        # Import the data file and process it
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
        circ$chrom<-str_remove(circ$chrom,"_")
        circ$chrom<-str_remove(circ$chrom,"chr")
        circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))

        # Use data from the plot to filter
        chrom_filter_small <- unlist(clickDataSmall()[3])
        size_filter_small <- unlist(clickDataSmall()[4])
        sr_filter_small <-unlist(clickDataSmall()[5])

        selected_circ<-circ %>% dplyr::filter (chrom==chrom_filter_small & size_bp==size_filter_small & split_reads==sr_filter_small)

        # Use Biostrings to get sequence
        my.dnastring <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0("chr",selected_circ[,1]), selected_circ[,2], selected_circ[,3])
        n<-as.vector(alphabetFrequency(my.dnastring)[1:4])

        # Make a data frame
        base <-c("A","C","G","T")
        data<-data.frame(base,n)
        data$fraction <- data$n / sum(data$n) # Percentages
        data$ymax <- cumsum(data$fraction) # Cumulative percentages (top of each rectangle)
        data$ymin <- c(0, head(data$ymax, n=-1)) # Bottom of each rectangle
        data$labelPosition <- (data$ymax + data$ymin) / 2 # Label position
        data$label <- paste0(data$base,":\n ",(round((data$fraction)*100,2)),"%") # Label

        # Use data frame to make circular plot
        ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=base)) +
            geom_rect() +
            geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
            scale_fill_brewer(palette=4) +
            coord_polar(theta="y") +
            xlim(c(1, 4)) +
            theme_void() +
            theme(legend.position = "none")
    })

    # TABLE WITH GENES IN CIRCLE
    output$GenesSmall<-  DT::renderDataTable({
        # Import an process the data file
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
        circ$chrom<-str_remove(circ$chrom,"_")
        circ$chrom<-str_remove(circ$chrom,"chr")
        circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))

        # Use the data from the previuos plot to filter
        chrom_filter_small <- unlist(clickDataSmall()[3])
        size_filter_small <- unlist(clickDataSmall()[4])
        sr_filter_small <-unlist(clickDataSmall()[5])

        selected_circ <- circ %>% dplyr::filter (chrom==chrom_filter_small & size_bp==size_filter_small & split_reads==sr_filter_small)

        # Get coordenates
        coords<-GRanges(paste0("chr",selected_circ[,1],":",selected_circ[,2],"-",selected_circ[,3]))

        # Find genes and annotate with gene symbols
        genes <-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
        genes_df <- as.data.frame(subsetByOverlaps(genes,coords))
        entrezid <- genes_df$gene_id
        hs <- org.Hs.eg.db
        genesymbol_df<- AnnotationDbi::select(hs,
                                              keys = entrezid,
                                              columns = c("ENTREZID", "SYMBOL"),
                                              keytype = "ENTREZID")
        names(genesymbol_df)[names(genesymbol_df) == "ENTREZID"] <- "gene_id"
        genes_smallcircle<-left_join(genes_df,genesymbol_df)

        # Fix names of columns for the table
        names(genes_smallcircle)<-c("Chromosome","Start","End","Width","Strand","Gene ID","Gene Symbol")

        # Make a null data frame
        null_gene<- as.data.frame("No genes (or parts of genes) found in this circle.")

        # Show results
        if(nrow(genes_smallcircle) == 0){
            DT::datatable(null_gene,
                          options = list(
                              searching = FALSE))
        }else{
            DT::datatable(genes_smallcircle,options = list(
                searching = FALSE))
        }
    })


    ## SOURCE: BIG CIRCLES PLOT
    # Change tab and render info when clicking in the plot of big circles
    observeEvent(event_data("plotly_click",source="bigcircleSource"), {
        circletab2 <- switch(input$tabs,
                            "bigcirc" = "circle")

        updateTabItems(session, "tabs", circletab2)
    })

    observeEvent(event_data("plotly_click",source="bigcircleSource"), {
        output$selected_circle<- renderUI({
            fluidRow(
                box(title=span(icon("circle-o"),"CIRCLE"), status="primary",width = 12, solidHeader = TRUE,
                    tableOutput("TableDataOutBig"),
                    plotOutput("CirclePlotBig")),
                box(title=span(icon("fal fa-dna"),"GENOMIC DATA"), status="primary",width = 12, solidHeader = TRUE,
                    dataTableOutput("GenesBig")))

        })
    })

    # Get the data from the plot click
    clickDataBig <- reactive({unlist(event_data(event = "plotly_click", source = "bigcircleSource", priority = "event"))
    })

    # SELECTED CIRCLE INFORMATION
    output$TableDataOutBig <- renderTable({
        # Import an process the data file
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
        circ$chrom<-str_remove(circ$chrom,"_")
        circ$chrom<-str_remove(circ$chrom,"chr")
        circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))

        # Use the data from the previuos plot to filter
        chrom_filter_big <- unlist(clickDataBig()[3])
        size_filter_big <- unlist(clickDataBig()[4])
        sr_filter_big <-unlist(clickDataBig()[5])
        circ %>%
            dplyr::filter (chrom==chrom_filter_big & size_bp==size_filter_big & split_reads==sr_filter_big) %>%
            dplyr::select(chrom,start,end,discordant_reads,split_reads,score,coverage_mean,size_bp) %>%
            dplyr::rename(Chromosome=chrom,Start=start,End=end,DiscordantReads=discordant_reads,SplitReads=split_reads,Score=score,Coverage=coverage_mean,Size=size_bp)
    })

    # CIRCULAR PLOT
    output$CirclePlotBig <-renderPlot({
        # Import an process the data file
        circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
        names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
        circ$size_bp <- circ$end - circ$start
        circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
        circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
        circ$chrom<-str_remove(circ$chrom,"_")
        circ$chrom<-str_remove(circ$chrom,"chr")
        circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))

        # Use the data from the previuos plot to filter
        chrom_filter_big <- unlist(clickDataBig()[3])
        size_filter_big <- unlist(clickDataBig()[4])
        sr_filter_big <-unlist(clickDataBig()[5])

        selected_circ<-circ %>% dplyr::filter (chrom==chrom_filter_big & size_bp==size_filter_big & split_reads==sr_filter_big)

        # Use Biostrings to get sequence
        my.dnastring <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0("chr",selected_circ[,1]), selected_circ[,2], selected_circ[,3])
        n<-as.vector(alphabetFrequency(my.dnastring)[1:4])

        # Make a data frame
        base <-c("A","C","G","T")
        data<-data.frame(base,n)
        data$fraction <- data$n / sum(data$n) # Percentage
        data$ymax <- cumsum(data$fraction) # Cumulative percentages (top of each rectangle)
        data$ymin <- c(0, head(data$ymax, n=-1)) # Bottom of each rectangle
        data$labelPosition <- (data$ymax + data$ymin) / 2 # Label position
        data$label <- paste0(data$base,":\n ",(round((data$fraction)*100,2)),"%") # Label

        # Use data frame to make circular plot
        ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=base)) +
            geom_rect() +
            geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
            scale_fill_brewer(palette=4) +
            coord_polar(theta="y") +
            xlim(c(1, 4)) +
            theme_void() +
            theme(legend.position = "none")

    })


    #TABLE WITH GENES IN CIRCLE
        output$GenesBig<-  DT::renderDataTable({
            # Import an process the data file
            circ <-read.table(input$bedfile$datapath,header = FALSE, sep="\t",stringsAsFactors=FALSE)
            names(circ) <- c("chrom","start","end","discordant_reads","split_reads","score","coverage_mean","coverage_sd","coverage_start", "coverage_end","coverage_cont")
            circ$size_bp <- circ$end - circ$start
            circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)
            circ$chrom <- substr(circ$chrom, start = 1, stop = 5)
            circ$chrom<-str_remove(circ$chrom,"_")
            circ$chrom<-str_remove(circ$chrom,"chr")
            circ$chrom <- factor(circ$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))

            # Use the data from the previuos plot to filter
            chrom_filter_big <- unlist(clickDataBig()[3])
            size_filter_big <- unlist(clickDataBig()[4])
            sr_filter_big <-unlist(clickDataBig()[5])

            selected_circ <- circ %>% dplyr::filter (chrom==chrom_filter_big & size_bp==size_filter_big & split_reads==sr_filter_big)

            # Get coordenates
            coords<-GRanges(paste0("chr",selected_circ[,1],":",selected_circ[,2],"-",selected_circ[,3]))

            # Find genes and annotate with gene symbols
            genes <-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
            genes_df <- as.data.frame(subsetByOverlaps(genes,coords))
            entrezid <- genes_df$gene_id
            hs <- org.Hs.eg.db
            genesymbol_df<- AnnotationDbi::select(hs,
                                                  keys = entrezid,
                                                  columns = c("ENTREZID", "SYMBOL"),
                                                  keytype = "ENTREZID")
            names(genesymbol_df)[names(genesymbol_df) == "ENTREZID"] <- "gene_id"
            genes_bigcircle<-left_join(genes_df,genesymbol_df)

            # Fix names of columns for the table
            names(genes_bigcircle)<-c("Chromosome","Start","End","Width","Strand","Gene ID","Gene Symbol")

            # Make a null data frame
            null_gene<- as.data.frame("No genes (or parts of genes) found in this circle.")

            # Show results
           if(nrow(genes_bigcircle) == 0){
               DT::datatable(null_gene,
                             options = list(
                   searching = FALSE))
            }else{
                DT::datatable(genes_bigcircle,options = list(
                    searching = FALSE))
                }

        })

}

shinyApp(ui=ui, server = server)
}
