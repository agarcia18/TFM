if("shiny" %in% rownames(installed.packages()) == FALSE) {install.packages("shiny")}
if("shinydashboard" %in% rownames(installed.packages()) == FALSE) {install.packages("shinydashboard")}

shiny::runApp("./visualeccdna_app")
