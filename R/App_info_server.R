App_info_server <- function(id, parent_session) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    ##citation information
    output$citation_request <- renderUI({
      fluidRow(
        column(
          width = 12,
          wellPanel(
            style = "background-color: #f9f9f9;",
            tags$strong("Please cite: "), 
            tags$span("Wijnant KN et al. Mapping NMD-coupled alternative splicing in iPSC-derived brain cells: a resource for therapeutic discovery. Journal. Year Month date; Pages. doi: #. PMID: #; PMCID: #."),
          )
        )
      )
    })
    
    ## SpliNEx usage
    output$splinex_features <- renderUI({
      fluidRow(
      column(
        width = 12,
        wellPanel(
          style = "background-color: #ffffff; ",
          tags$div(
            style = "font-size: 20px; font-weight: bold; margin-bottom: 10px;",
            "Use SpliNEx to:"
          ),
          tags$div(style = "margin-top: 10px;"),
            tags$ul(
              tags$li("Explore gene expression in wild-type iPSC-derived iNeurons, astrocytes, and microglia", 
                      actionLink(ns("go_gene_expression_page"), "[view expression for gene]")),
              tags$li("Search for alternative splicing events identified by long-read sequencing",
                      actionLink(ns("go_event_per_gene_page_1"), "[view AS events for gene]")),
              tags$li("Access detailed information on alternative splice events",
                      actionLink(ns("go_event_per_gene_page_2"), "[view AS events for gene]"),
                      "or",
                      actionLink(ns("go_event_detail_page"), "[view AS event details for AS event id]"),
                      ", including:",
                      tags$ul(
                        tags$li("Estimated usage of each splice event across iPSC-derived cells of the brain"),
                        tags$li("Predictions of NMD-sensitive alternative splicing"),
                        tags$li("Related splice junctions and alternative splice events in the same genomic region")
                      )
              )
            )
        )
      )
    )
    })
    
    ##navigate to other tabs
    observeEvent(input$go_gene_expression_page, {
      updateTabsetPanel(parent_session, "main_tabs", selected = "wtc")
    })
    observeEvent(input$go_event_per_gene_page_1, {
      updateTabsetPanel(parent_session, "main_tabs", selected = "nmd_event")
    })
    observeEvent(input$go_event_per_gene_page_2, {
      updateTabsetPanel(parent_session, "main_tabs", selected = "nmd_event")
    })
    observeEvent(input$go_event_detail_page, {
      updateTabsetPanel(parent_session, "main_tabs", selected = "nmd_event_details")
    })
    
    ##Data generation and explanation
    output$data_generation <- renderUI({
      fluidRow(
        column(
          width = 12,
          wellPanel(
            style = "background-color: #ffffff; ",
            tags$div(
              style = "font-size: 20px; font-weight: bold; margin-bottom: 10px;",
              "Alternative splicing data generation:"
            ),
            tags$div(style = "margin-top: 10px;"),
            tags$p(
              "Alternative splicing events were identified using long-read sequencing data from iPSC-derived iNeurons, astrocytes, and microglia, 
   and their quantification was based on short-read sequencing data. 
   For a detailed explanation of alternative splicing data generation: ",
              actionLink("link", "[Wijnant et al. 2025]"), "."
            ),
            img(src = "Data_generation.png", height = "400px"),
          )
        )
      )
    })
    
    ##Explenation of concepts and 
    output$explenation_concepts <- 
      renderUI({
        fluidRow(
          column(
            width = 12,
            wellPanel(
              style = "background-color: #ffffff; ",
              tags$div(
                style = "font-size: 20px; font-weight: bold; margin-bottom: 10px;",
                "Key concepts & terms used in this app:"
              ),
              tags$div(style = "margin-top: 10px;"),
              tags$ul(
                tags$li("Event ID & AS event type:"),
                tags$p(
                  "A unique identifier for each alternative splicing event, specifying its gene, genomic coordinates, AS event type and strand, as defined by the SUPPA tool. A detailed explanation of the event ID and the different AS event types can be found at the",
                  tags$a("SUPPA GitHub page", href = "https://github.com/comprna/SUPPA", target = "_blank")
                ),
                
                tags$li("Percentage spliced in (PSI)"),
                tags$p(
                  "PSI represents the fraction of transcripts in which a given alternative splicing event is included versus excluded, calculated from short-read RNA-seq splice-junction read counts, and ranges from 0 (never included) to 1 (always included). PSI provides an estimate of AS event usage, but alternative splicing events sharing the same splice junction(s) can negatively affect the calculated PSI."
                ),
                
                tags$li("Relative event coverage"),
                tags$p(
                  "Relative event coverage represents the proportion of long-read (LR) transcripts supporting a given alternative splicing (AS) event (either inclusion or exclusion) relative to the total LR coverage in the region, providing an estimate of the relevance of the AS event within the gene (calculated as: # supporting AS event LR reads / # total LR reads spanning the AS event region). Low to medium relative event coverage (0–0.8) may indicate spurious AS event annotation, complex splicing in the region, or inaccurate PSI calculations."
                ),
                
                tags$li("NMD AS prediction methods"),
                tags$p(
                  "Two complementary NMD AS prediction methods: (i) a more inclusive, less stringent method which labels an AS event as NMD-inducing if the AS event occurs exclusively in isoforms annotated as NMD sensitive (NMD-exclusive method), and (ii) a novel, more stringent LR-based splice-context aware method, which compares the protein-coding potential of isoforms with and without the AS event within the same splice context."
                ),
                
                img(src = "NMD_AS_methods.png", height = "300px"),
                
                tags$p(
                  "For all NMD AS event types and a more detailed explanation & evaluation of both NMD AS prediction methods see:",
                  actionLink("link", "Wijnant et al. 2025")
                )
              ),
              tags$p(
                "For a more detailed explanation of terms, methods and calculations see:",
                actionLink("link", "Wijnant et al. 2025")
              )
          )
        )
        )
      })
    
  })
}

