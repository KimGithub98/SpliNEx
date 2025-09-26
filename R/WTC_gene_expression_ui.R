wtc_gene_expression_ui <- function(id) {
  ns <- NS(id)
  tabPanel("WTC gene expression", 
           sidebarLayout(
             sidebarPanel(
               checkboxGroupInput(ns("selected_batches"), "Samples to show:",
                                  choices = names(batch_list), 
                                  selected = names(batch_list)),
               textInput(ns("wtc_gene"), "Gene:", "HGNC symbol"),
               actionButton(ns("submit_wtc"), "Submit")
             ),
             
             mainPanel(
               plotOutput(ns("wtc_gene_plot")),
               verbatimTextOutput(ns("wtc_gene_output"))
             )
           )
  )
}