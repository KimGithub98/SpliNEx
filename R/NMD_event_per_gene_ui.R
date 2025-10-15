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
           HTML("
           <b>References & data sources:</b><br>
           1. Wijnant et al. 2025<br><br>
           <b>Please city: Wijnant et al. 2025</b><br>
           Build with Shiny & Shinythemes in R version 4.4.1")
  )
}