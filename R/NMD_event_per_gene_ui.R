nmd_event_per_gene_ui <- function(id) {
  ns <- NS(id)
  tabPanel("NMD events", 
           sidebarLayout(
             sidebarPanel(
               textInput(ns("nmd_gene"), "Gene:", "HGNC symbol"),
              
               uiOutput(ns("message_NMD_event_per_gene")),
               checkboxGroupInput(ns("selected_NMD_event_types"), "NMD AS types to show (Splice-context aware method):",
                                  choices = NMD_event_types, 
                                  selected = NMD_event_types),
               checkboxGroupInput(ns("selected_NMD_only_types"), "NMD AS types to show (NMD-exclusive method):",
                                  choices = NMD_only_types, 
                                  selected = NMD_only_types),
               checkboxGroupInput(ns("selected_event_types"), "AS event types to show:",
                                  choices = event_types, 
                                  selected = event_types),
               actionButton(ns("submit_nmd"), "Submit")
             ),
             
             mainPanel(
               uiOutput(ns("gene_title")),
               fluidRow(
                 column(width = 6, plotOutput(ns("nmd_LR_gene_plot"), height = "400px")),
                 column(width = 6, plotOutput(ns("nmd_SR_gene_plot"), height = "400px"))
               ),
               uiOutput(ns("coverage_legend")),
               conditionalPanel(condition = sprintf("output['%s'] == 'true'", ns("data_ready")),
                                fluidRow(div(
                                  style = "text-align: right; margin-bottom: 10px;",
                                  downloadButton(ns("download_nmd_event_per_gene_table"), "Download table as CSV"),
                                  downloadButton(ns("download_nmd_event_per_gene_full_table"), "Download extended table as CSV"),
                                ))
               ),
               fluidRow(
                 column(12,  # full width
                        DT::dataTableOutput(ns("NMD_table"))
                 )
               ),
               verbatimTextOutput(ns("nmd_gene_output"))
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