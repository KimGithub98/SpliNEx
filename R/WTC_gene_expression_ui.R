pretty_names <- c(
  "Astrocytes (CHX, Wijnant et al. 2025)"   = "WTC_Astrocytes_CHX_K",
  "Astrocytes (DMSO, Wijnant et al. 2025)"  = "WTC_Astrocytes_DMSO_K",
  "DIV21 Neurons (Elly et al. 2025)"      = "WTC_DIV21_S",
  "DIV28 Neurons (Emma et al. 2025)"      = "WTC_DIV28_S",
  "DIV49 Neurons (Marina et al. 2025)" = "WTC_DIV49_batch1_S",
  "DIV49 Neurons (Xiuming et al. 2025)" = "WTC_DIV49_batch2_S",
  "DIV49 (CHX, Wijnant et al. 2025)"        = "WTC_DIV49_CHX_K",
  "DIV49 (DMSO, Wijnant et al. 2025)"       = "WTC_DIV49_DMSO_K",
  "iPSC (CHX, Wijnant et al. 2025)"         = "WTC_iPSC_CHX_K",
  "iPSC (DMSO, Wijnant et al. 2025)"        = "WTC_iPSC_DMSO_K",
  "Microglia (CHX, Wijnant et al. 2025)"    = "WTC_Microglia_CHX_K",
  "Microglia (DMSO, Wijnant et al. 2025)"   = "WTC_Microglia_DMSO_K"
)


wtc_gene_expression_ui <- function(id) {
  ns <- NS(id)
  tabPanel("WTC gene expression", 
           sidebarLayout(
             sidebarPanel(
               checkboxGroupInput(ns("selected_batches"), "Samples to show:",
                                  choices = pretty_names, 
                                  selected = pretty_names),
               textInput(ns("wtc_gene"), "Gene:", "HGNC symbol"),
               uiOutput(ns("message")),
               actionButton(ns("submit_wtc"), "Submit")
             ),
             
             mainPanel(
               plotOutput(ns("wtc_gene_plot")),
               verbatimTextOutput(ns("wtc_gene_output"))
             )
           ),
           hr(),
           HTML("
           <b>References & data sources:</b><br>
           1. Wijnant et al. 2025<br>
           2. Elly et al. 2025<br>
           3. Emma et al. 2025<br>
           4. Marina et al. 2025<br>
           5. Xiuming et al. 2025<br><br>
           <b>Please city: Wijnant et al. 2025</b><br>
           Build with Shiny & Shinythemes in R version 4.4.1")
  )
}