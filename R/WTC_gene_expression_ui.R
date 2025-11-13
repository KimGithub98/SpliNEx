pretty_names <- c(
  "Astrocytes (CHX, Wijnant et al. 2025)"   = "WTC_Astrocytes_CHX_K",
  "Astrocytes (DMSO, Wijnant et al. 2025)"  = "WTC_Astrocytes_DMSO_K",
  "Neurons (DIV21, Leonardi et al. 2025)"      = "WTC_DIV21_S",
  "Neurons (DIV28, Emma et al. 2025)"      = "WTC_DIV28_S",
  "Neurons (DIV49, Hommersom et al. 2025)" = "WTC_DIV49_batch1_S",
  "Neurons (DIV49, Yuan et al. 2023)" = "WTC_DIV49_batch2_S",
  "Neurons (DIV49, CHX, Wijnant et al. 2025)"        = "WTC_DIV49_CHX_K",
  "Neurons (DIV49, DMSO, Wijnant et al. 2025)"       = "WTC_DIV49_DMSO_K",
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
               conditionalPanel(condition = sprintf("output['%s'] == 'true'", ns("plot_ready")),
                                fluidRow(downloadButton(ns("download_gene_exp_plot_png"), "Download Plot as PNG"),
                                         downloadButton(ns("download_gene_exp_plot_svg"), "Download Plot as SVG"))
                                ),
               fluidRow(
                 column(8,
                        plotOutput(ns("wtc_gene_plot")),
                        verbatimTextOutput(ns("wtc_gene_output"))
                        )
                 
               )
               
             )
           ),
           hr(),
           tags$div(
             tags$b("References & data sources:"), tags$br(),
             "1. ", tags$a("Wijnant et al. 2025", href = "https://doi.org/YOUR_DOI_1", target = "_blank"), tags$br(),
             "2. ", tags$a("Leonardi O et al. CHD2 Dosage Ties Autolysosomal Pathway to Cortical Maturation in Disease and Evolution. bioRxiv. 2025 Jan 21.", 
                           href = "https://doi.org/10.1101/2025.01.21.634145", target = "_blank"), tags$br(),
             "3. ", tags$a("Emma et al. 2025", href = "https://doi.org/YOUR_DOI_3", target = "_blank"), tags$br(),
             "4. ", tags$a("Hommersom MP et al. CACNA1A haploinsufficiency leads to reduced synaptic function and increased intrinsic excitability. Brain. 2025 Apr 3;148(4):1286-1301.", 
                           href = "https://doi.org/10.1093/brain/awae330", target = "_blank"), tags$br(),
             "5. ", tags$a("Yuan X et al. A human in vitro neuronal model for studying homeostatic plasticity at the network level. Stem Cell Reports. 2023 Nov 14;18(11):2222-2239.", 
                           href = "https://doi.org/10.1016/j.stemcr.2023.09.011", target = "_blank"), tags$br(), tags$br(),
             tags$b("Please cite: Wijnant et al. 2025"), tags$br(),
             "Build with Shiny & Shinythemes in R version 4.4.1"
           )
  )
}