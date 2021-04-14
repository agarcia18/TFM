# Load required packages
library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggplot2)
library(plotly)
library(randDNA)

circ <- as.data.frame(read.table("C:/Users/ainho/eccDNA/BED/SRR6315412_unknown_circle.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))

# Change names
names(circ)[1] <- "chr"
names(circ)[2] <- "start"
names(circ)[3] <- "end"
names(circ)[4] <- "discordant_reads"
names(circ)[5] <- "split_reads"
names(circ)[6] <- "score"
names(circ)[7] <- "coverage_mean"
names(circ)[8] <- "coverage_sd"
names(circ)[9] <- "coverage_start"
names(circ)[10] <- "coverage_end"
names(circ)[11] <- "coverage_cont"

# Set chromosome as factor
circ$chr <- substr(circ$chr, start = 1, stop = 5)
circ$chr<-str_remove(circ$chr,"_")
circ$chr<-str_remove(circ$chr,"chr")
circ$chr <- factor(circ$chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Un","M","X","Y"))

# Add quality levels (score < 10 = Bad, score < 50 = Low, score < 200 = Medium, score > 200 = Good)
circ$quality <- cut(circ$score,breaks=c(-Inf,10,50,200,Inf),labels= c("Bad","Low", "Medium", "Good"),right = FALSE)

# Add size of circle (number of bp)
size_bp <- circ$end - circ$start
circ<- add_column(circ,size_bp)
head(circ)



# Get Chromosome Input choices
all <- "All"
chr_choices <- append(all,levels(circ$chr))

# Random sequence

# User interface
ui <- dashboardPage(
    dashboardHeader(title = "Visual EccDNA"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("About",tabName="about",icon=icon("tasks")),
            menuItem("Results", tabName ="results", icon=icon("bar-chart")),
            menuItem("Circle Info",tabName="circle",icon=icon("circle-o"))
        )
    ),
    
    dashboardBody(
        tabItems(
            tabItem(
                tabName = "about",
                    fluidRow(
                        h2("Wellcome!")),
                    fluidRow(
                        box(title= "Upload a Circle BED file with output:",fileInput("bedfile","Choose file:"))
                    ),
            ),
            tabItem(
                tabName = "results",
                    fluidRow(
                        valueBox(total_circles,tags$b("TOTAL DNA CIRCLES"), color="maroon",
                                 icon = tags$i(class = "fas fa-dna", style="font-size: 50px")
                                 ),
                        box(width = 3,background = "olive",tags$h4("GOOD:",paste0(round(p_good,2),"%"))),
                        box(width = 3, background = "yellow",tags$h4("LOW:", paste0(round(p_low,2),"%"))),
                        box(width = 3, background = "light-blue",tags$h4("MEDIUM:",paste0(round(p_medium,2),"%"))),
                        box(width = 3, background = "red",tags$h4("BAD:", paste0(round(p_bad,2),"%"))),
                        ),
                    fluidRow(
                        box(plotOutput("results"),width=7,
                            status = "primary",
                            ),
                        box(title= span(icon("filter"), "Filters"), status = "warning", width = 5,
                            checkboxGroupInput(inputId = "quality",
                                               "Quality of your circles:",
                                               choices = c("Bad","Low","Medium","Good"),
                                               selected= c("Medium","Good")
                            ),
                            sliderInput(inputId = "size",
                                        "By size (number of base pairs):",
                                        min = min(circ$size_bp),
                                        max= max(circ$size_bp),
                                        value = c(min(circ$size_bp),max(circ$size_bp))),
                            selectInput(inputId = "chr",
                                        "By source chromosome:",
                                        choices = chr_choices),
                           
    )
                ),
            ),
            tabItem(
                tabName="circle",
                        h2("Here goes the info of each circle")
                    )
        )
    )
)
