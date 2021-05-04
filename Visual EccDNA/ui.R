# Load required packages 
library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)

# User interface
ui <- dashboardPage(
    dashboardHeader(title = "Visual EccDNA"),
    
    dashboardSidebar(
        sidebarMenu(
            id="tabs",
            menuItem("Get started",tabName="about",icon=icon("cloud-upload")),
            menuItem("View results", tabName ="results", icon=icon("bar-chart")),
            menuItem("Circle Info",tabName="circle",icon=icon("circle-o")),
            menuItem("GitHub",tabName = "github",icon=icon("github"))
        )
    ),
    
    dashboardBody(
        tabItems(
            tabItem(
                tabName = "about",
                    fluidRow(
                        box(
                        title = "Wellcome to Visual eccDNA", status="primary",width = 12, solidHeader = TRUE,
                        "This is an app for dynamic data visualization of eccDNA output files. Start selecting a BED file and this page will show an overview of the results. Then you can press the button View Results, and you will be able to see further information about the DNA circles contained in your file."
                    )),
                    fluidRow(
                            box(title= "Upload a Circle BED file with output:",fileInput("bedfile","Choose file:"),status="primary"),
                            uiOutput("start")
                        ),
                    fluidRow(uiOutput("total"),
                             uiOutput("quality_gg")
                             ),
                    fluidRow(uiOutput("table"))
            ),
            tabItem(
                tabName = "results",
                    fluidRow(uiOutput("plot_results"),
                             uiOutput("filters")),
                    fluidRow(uiOutput("table_bigcircles")),
                        
                    ),

            tabItem(
                tabName="circle",
                        h2("Here goes the info of each circle"),
                tags$img(src="under_construction.png",align="center")
                ),
            
            tabItem(
                tabName="github",
                tags$a(href="https://github.com/agarcia18/eccDNA","https://github.com/agarcia18/eccDNA")
            )
        )
    )
)


