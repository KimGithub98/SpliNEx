App_info_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Info", 
           titlePanel(
             title = span(
               style = "display: flex; align-items: flex-end; gap: 10px;",
               img(src = "SpliNEx_full_logo.jpg", height = "50px"),
               " : Splicing, NMD-sensitive splicing and gene expression in iPSC-derived cells of the Brain"
             )
           ),
           uiOutput(ns("citation_request")),
           uiOutput(ns("splinex_features")),
           uiOutput(ns("data_generation")),
           uiOutput(ns("explenation_concepts")),
           hr(),
           HTML("
           <b>References & data sources:</b><br>
           1. Wijnant et al. 2025<br>
           <b>Please cite: Wijnant et al. 2025</b><br>
           Build with Shiny & Shinythemes in R version 4.4.1")
  )
}