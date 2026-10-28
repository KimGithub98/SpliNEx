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
               column(width = 7,
                      uiOutput(ns("event_details_print"))
                      ),
               column(width = 4, div(title = "Shows the relative event coverage for the selected event per long-read sequencing sample.", 
                                     plotOutput(ns("perc_region_plot"), height = "280px"))),
               ),
             
             # PSI plot
             fluidRow(
               column(width = 5, div(title = "Shows the percent spliced in (PSI) for the selected event across samples.",
                                     plotOutput(ns("event_PSI_plot"), height = "500px"))),
               column(width = 6,  # full width
                      DT::dataTableOutput(ns("event_psi_DT"))
               )
               
               
             ),
             # Sashimi plot and UCSC button
             fluidRow(
               fluidRow(
                 column(width = 11, div(title = "Shows the junctions and range of junction-spanning read counts for the selected event and any selected interfering events across short-read RNA-seq samples. CLICK ON JUNCTION LABEL to view the corresponding junction and event information", 
                        plotOutput(ns("Faux_sashimi_plot"), click = ns("faux_click"), height = "800px")),
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
             ),
             
             #related junctions table
             fluidRow(
               tags$h4("Table of interfering AS events:", style = "margin-top: 25px;"),
               column(12,  # full width
                      DT::dataTableOutput(ns("related_junction_DT"))
               )
             )
             
           ),
           hr(),
           tags$div(
             tags$b("References & data sources:"), tags$br(),
             "1. ", tags$a("Wijnant et al. Mapping NMD-coupled alternative splicing in iPSC-derived brain cells: a resource for therapeutic discovery. bioRxiv. 2026 Jan 22", href = "https://doi.org/10.64898/2026.01.19.700316", target = "_blank"), tags$br(),tags$br(),
             tags$b("Please cite: Wijnant et al. 2026"), tags$br(),
             "Build with Shiny & Shinythemes in R version 4.4.1"
           )
  )
  
}