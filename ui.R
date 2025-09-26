fluidPage(
  theme = shinytheme("flatly"),
  
  navbarPage("My App",
             id = "main_tabs",
             #tab1 wtc gene expression
             wtc_gene_expression_ui("wtc"),
             #tab2 nmd event per gene
             nmd_event_per_gene_ui("nmd_event"),
             #tab3 nmd event details
             event_details_ui("nmd_event_details")
  )
)