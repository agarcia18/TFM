# Load required packages
library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggplot2)
library(ggiraph)
library(plotly)


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
                    fluidRow(h2("Wellcome!")),
                    fluidRow(
                        box(title= "Upload a Circle BED file with output:",fileInput("bedfile","Choose file:"),status="primary"),
                        uiOutput("ui_total")
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


