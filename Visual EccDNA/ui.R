# Load required packages
library(shiny)
library(shinydashboard)


ui <- dashboardPage(
    dashboardHeader(title = "Visual EccDNA"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("About", tabName ="about", icon=icon("clipboard"))
        )
    ),
    dashboardBody(
        # Boxes need to be put in a row (or column)
        fluidRow(
            box(plotOutput("plot1", height = 250)),
            
            box(
                title = "Controls",
                sliderInput("slider", "Number of observations:", 1, 100, 50)
            )
        )
    )
)
