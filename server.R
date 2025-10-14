server <- function(input, output, session) {
  #tab1 wtc gene expression
  wtc_gene_expression_server("wtc")
  
  #tab2 nmd event per gene
  nmd_event_per_gene_server("nmd_event")
  
  #tab3 nmd event details
  selected_event <- reactive({input$selected_event})  
  event_details_server("nmd_event_details", selected_event = reactive(input$selected_event))
  
  #auto switch from nmd_event_per_gene tab to nmd_event_detail tab when clicked
  observeEvent(input$selected_event, {
    updateTabsetPanel(session, "main_tabs", selected = "Event details")
  })
  
  observe({
    print(paste("Memory used:", pryr::mem_used() / 1024^2, "MB"))
  })
}
