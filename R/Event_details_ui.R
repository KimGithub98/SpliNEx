event_details_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Event details", 
           fluidPage(
             # Input + submit on one line
             fluidRow(
               column(width = 12,
                      div(style = "display: flex; align-items: center;",
                          textInput(ns("event_input"), "AS event:", "Event id", width = "300px"),
                          actionButton(ns("submit_event"), "Submit", style = "margin-left: 10px; height: 50px;")
                      )
               )
             ),
             
             # PSI plot
             fluidRow(
               column(width = 5, plotOutput(ns("event_PSI_plot"), height = "500px")),
               column(width = 6,  # full width
                      DT::dataTableOutput(ns("event_psi_DT"))
               )
               
             ),
             # Sashimi plot and UCSC button
             fluidRow(
               column(width = 6, 
                      actionButton(ns("go_button"), "View transcripts in UCSC genome browser", class = "btn btn-primary")
               ),
               column(width = 6,
                      checkboxInput(ns("show_related_junctions"), "Show related junctions in region", FALSE)
               ),
               fluidRow(
                 column(width = 11, 
                        plotOutput(ns("Faux_sashimi_plot"), click = ns("faux_click"), height = "800px"),
                        fluidRow(
                          column(width = 8,
                                 verbatimTextOutput(ns("faux_info"))
                          ),
                          column(width = 4,
                                 uiOutput(ns("goto_event_ui"))
                          )
                        )
                 ),
                 verbatimTextOutput(ns("nmd_gene_output"))
               )
             )
           ),
           hr(),
           HTML("
           <b>References & data sources:</b><br>
           1. Wijnant et al. 2025<br><br>
           <b>Please city: Wijnant et al. 2025</b><br>
           Build with Shiny & Shinythemes in R version 4.4.1")
  )
  
}