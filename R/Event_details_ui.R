event_details_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Event details", 
           fluidPage(
             # Input + submit on one line
             fluidRow(
               column(width = 4,
                      div(style = "display: flex; align-items: center;",
                          textInput(ns("event_input"), "Alternative splice event (Event ID):", "Event id", width = "300px"),
                          actionButton(ns("submit_event"), "Submit", style = "margin-left: 10px; height: 50px;")
                      )
               ),
               column(width = 4,
                      uiOutput(ns("message"))),
               column(width = 4, 
                      uiOutput(ns("ucsc_button_ui"))
               )
             ),
             
             #NMD_prediction/event details 
             fluidRow(
               column(width = 11,
                      uiOutput(ns("event_details_print"))
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
                 column(width = 6,
                        uiOutput(ns("checkbox_related_junctions"))
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