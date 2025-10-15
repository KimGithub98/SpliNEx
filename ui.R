fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),
  
  # Change browser tab title
  tags$head(
    tags$title("SpliNEx"),
    tags$link(rel = "icon", type = "image/png", href = "SpliNEx_small_logo.png")
  ),
  
  navbarPage(title = img(src = "SpliNEx_small_logo.png", height = "30px"),
             id = "main_tabs",
             #info tab
             tabPanel(title = "Info", value = "info",
                      App_info_ui("app_info")  # keep the module UI inside
             ),
             #tab1 wtc gene expression
             tabPanel(title = "Gene Expression", value = "wtc",
                      wtc_gene_expression_ui("wtc")
             ),
             #tab2 nmd event per gene
             tabPanel(title = "Alternative splicing", value = "nmd_event",
                      nmd_event_per_gene_ui("nmd_event")
             ),
             #tab3 nmd event details
             tabPanel(title = "Alternative splicing event details", value = "nmd_event_details",
                      event_details_ui("nmd_event_details")
             )
  )
)