# Load required packages
library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggplot2)
library(ggiraph)
library(plotly)
library(DT)

# User interface
ui <- dashboardPage(
    dashboardHeader(title = "Visual EccDNA"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Get started",tabName="about",icon=icon("tasks")),
            menuItem("View results", tabName ="results", icon=icon("bar-chart")),
            menuItem("Circle Info",tabName="circle",icon=icon("circle-o"))
        )
    ),
    
    dashboardBody(
        tabItems(
            tabItem(
                tabName = "about",
                    fluidRow(
                        box(
                        title = "Wellcome to Visual eccDNA", status="primary",width = 12, solidHeader = TRUE,
                        "This is an app for dynamic data visualization of eccDNA output files. Start selecting a BED file and this page will show an overview of the results. Then press Get Started! button, and you will be able to see further information about the DNA circles contained in your file, in View Results page."
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
                    fluidRow(
                        box(plotOutput("results"),width=7,
                            status = "primary",
                            ),
                    ),
            ),
            tabItem(
                tabName="circle",
                        h2("Here goes the info of each circle")
                )
        )
    )
)


