nmd_event_per_gene_ui <- function(id) {
  ns <- NS(id)
  tabPanel("NMD events", 
           sidebarLayout(
             sidebarPanel(
               textInput(ns("nmd_gene"), "Gene:", "HGNC symbol"),
               checkboxGroupInput(ns("selected_NMD_event_types"), "NMD types to show (sister transcript):",
                                  choices = NMD_event_types, 
                                  selected = NMD_event_types),
               checkboxGroupInput(ns("selected_NMD_only_types"), "NMD types to show (NMD exclusive method):",
                                  choices = NMD_only_types, 
                                  selected = NMD_only_types),
               checkboxGroupInput(ns("selected_event_types"), "Event types to show:",
                                  choices = event_types, 
                                  selected = event_types),
               actionButton(ns("submit_nmd"), "Submit")
             ),
             
             mainPanel(
               fluidRow(
                 column(width = 6, plotOutput(ns("nmd_LR_gene_plot"), height = "400px")),
                 column(width = 6, plotOutput(ns("nmd_SR_gene_plot"), height = "400px"))
               ),
               fluidRow(
                 column(12,  # full width
                        DT::dataTableOutput(ns("NMD_table"))
                 )
               ),
               verbatimTextOutput(ns("nmd_gene_output"))
             )
           )
  )
}