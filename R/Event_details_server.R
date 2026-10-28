event_details_server <- function(id, selected_event = selected_event,parent_session = session) {
  moduleServer(id, function(input, output, session) {
    #--- reactive sates --- 
    event_msg <- reactiveVal(NULL)
    selected_event_id <- reactiveVal(NULL)
    clicked_junction <- reactiveVal(NULL)
    selected_related_events <- reactiveVal(character(0))
    
    #--- colors --- 
    div_colors <- c(
      "iPSCs"              = "#cfcece",
      "Neurons (DIV21)"    = "#cbf4ce",
      "Neurons (DIV28)"    = "#95ba9a",
      "Neurons (DIV49)"    = "#637b64",
      "Astrocytes (DIVx)"  = "#b29ac9",
      "Microglia (DIVx)"   = "#9ebae1"
    )
    fill_colors <- c(
      "iPSC.DMSO" = "#cfcece",
      "iPSC.CHX" = "#727272",
      "Neurons.DMSO" = "#95ba9a",
      "Neurons.CHX" = "#637b64",
      "Astro.DMSO" = "#b29ac9",
      "Astro.CHX" = "#71627f",
      "Micro.DMSO" = "#9ebae1",
      "Micro.CHX" = "#637489"
    )
    
    #--- reactives --- 
    ##get junction event data
    junction_event_data <- reactive({
      req(selected_event_id())
      
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      
      event_type <- row$event_type
      strand <- row$strand
      
      if (event_type == "RI") {
        return(NULL)
      }
      
      get_junction_event_data(
        event,
        event_type,
        row,
        junction_summary_for_app,
        strand
      )
    })
    
    ##gets sashimi plot
    faux_sashimi_plot <- reactive({
      req(selected_event_id())
      
      selected_region_junction_df
      
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      
      event_type <- row$event_type
      chrom <- row$chrom
      
      if (event_type == "RI") {
        return(NULL)
      }
      
      jd <- junction_event_data()
      req(jd)
      
      # Get currently selected event IDs
      selected <- selected_related_events()
      
      region_junction_df <- jd$region_junction_df %>%
        mutate(
          show_interfering_event = event_id %in% selected
        )
      
      exon_df_region <- jd$exon_df_region %>%
        mutate(
          show_interfering_event = event_id %in% selected
        )
      
      x_axis_label <- paste0("Genomic Position (chr", chrom, ")")
      
      get_faux_sashimi_plot_v2(
        event_type,
        junction_df = jd$junction_df,
        region_junction_df = region_junction_df,
        x_axis_label,
        jd$exon_df,
        exon_df_region,
        fill_colors
      )
    })
    
    ##gets selected junction data
    selected_region_junction_df <- reactive({
      jd <- junction_event_data()
      req(jd)
      
      selected <- selected_related_events()
      
      jd$region_junction_df %>%
        mutate(
          selected = event_id %in% selected
        )
    })
    
    #--- observers--- 
    ## Update reactive when user clicks in NMD tab
    observeEvent(selected_event(), {
      selected_related_events(character(0))
      event_msg(NULL)
      selected_event_id(selected_event())
      updateTextInput(session, "event_input", value = selected_event())
      updateTabsetPanel(parent_session, "main_tabs", selected = "nmd_event_details")
    })
    
    ## Update reactive when user submits an event manually
    observeEvent(input$submit_event, {
      selected_related_events(character(0))
      event <- as.character(input$event_input)
      event_msg(NULL)
      if (event %in% NMD_table_for_app$event_id) {
        selected_event_id(event)
      }
      
      if (!(event %in% NMD_table_for_app$event_id)) {
        event_msg(tags$span(style = "color:red;", "Event not found in dataset."))
        return()
      }
    })
    
    ## Update reactive when user clicks a related event
    observeEvent(input$goto_event, {
      req(input$faux_click)
      
      # Get the reactive junction data
      jd <- junction_event_data()
      req(jd)
      
      # Get the current sashimi plot
      p <- faux_sashimi_plot()
      req(p)
      
      selected <- selected_related_events()
      
      region_junction_df <- jd$region_junction_df %>%
        mutate(
          show_interfering_event = event_id %in% selected
        )
      
      # Get plot coordinates
      pb <- ggplot_build(p)
      
      x_range <- pb$layout$panel_params[[1]]$x.range
      y_range <- pb$layout$panel_params[[1]]$y.range
      
      # Determine which junction/event was clicked
      click_selected_row <- return_clicked_junction(
        click = input$faux_click,
        junction_df = jd$junction_df,
        region_junction_df = region_junction_df,
        x_range = x_range,
        y_range = y_range
      )
      
      new_event <- click_selected_row$event_id
      
      selected_related_events(character(0))
      
      clicked_junction(new_event)
      
      # Only update if it actually changed
      if (!identical(selected_event_id(), new_event)) {
        selected_event_id(new_event)
      }
    })
    
    ##links click on button to UCSC genome browser
    observeEvent(input$go_button, {
      req(selected_event_id())
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      url <- get_ucsc_url(event, row$event_type, row$chrom)
      shinyjs::runjs(sprintf("window.open('%s', '_blank')", url))
      shinyjs::disable("go_button") 
    })
    
    ##observe if checkbox is clicked
    observeEvent(input$dtable_related_junction_check, {
      
      info <- input$dtable_related_junction_check
      
      event_id <- info[["event_id"]]
      checked <- as.logical(info[["value"]])
      
      current <- selected_related_events()
      
      if (checked) {
        current <- union(current, event_id)
      } else {
        current <- setdiff(current, event_id)
      }
      
      selected_related_events(current)
    })
    
    #--- Output --- 
    ##get event details
    output$event_details_print <- renderUI({
      req(selected_event_id())
      #get event details
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      event_coverage <- paste0(round(row$perc_region_avg*100, 0), "%")
      event_type <- row$event_type
      event_type <- names(event_choices)[event_choices == event_type]
      
      fluidRow(
        column(                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
          width = 12,
          wellPanel(
            style = "background-color: #f9f9f9; border: 1px solid #ccc;",
            tags$strong("Selected event ID:"),
            tags$span(event),
            tags$br(),
            tags$br(),
            tags$strong("Selected event type:"),
            tags$span(event_type),
            tags$br(),
            tags$br(),
            tags$strong("NMD-prediction (splice-context aware method): "), 
            tags$br(),
            tags$span(row$NMD_sister_transcript_summary),
            tags$br(),
            tags$span("Conclusion: "),
            tags$strong(row$NMD_sister_transcript),
            tags$br(),
            tags$br(),
            tags$strong("Relative event coverage: "), 
            tags$span(event_coverage),
          )
        )
      )
      #print event details
    })
    
    ##get PSI plot 
    output$event_PSI_plot <- renderPlot({
      req(selected_event_id())
      
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event,]
      PSI_df <- get_PSI_event_data(row)
      
      get_PSI_plot(PSI_df, fill_colors)
    })
    
    ##get relative event coverage plot
    output$perc_region_plot <- renderPlot({
      req(selected_event_id())
      
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event,]
      get_perc_region_plot(row)
    })
    
    ##get PSI datatable
    output$event_psi_DT <- DT::renderDT({
      req(selected_event_id())
      
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event,]
      PSI_df <- get_PSI_event_data(row)
      
      DT::datatable(get_PSI_table_for_printing(PSI_df), escape = FALSE, options = list(scrollX = TRUE))
    })
    
    ##output faux sashimi plot
    output$Faux_sashimi_plot <- renderPlot({
      req(selected_event_id())
      
      p <- faux_sashimi_plot()
      
      if (is.null(p)) {
        plot.new()
        text(
          0.5, 0.5,
          "Sashimi plot is not available for retained intron (RI) events.",
          col = "red",
          cex = 1.5
        )
      } else {
        p
      }
      
      print(class(p))
      print(p)
    })
    
    ##output related junction data table
    output$related_junction_DT <- DT::renderDT({
      req(selected_event_id())
      
      jd <- junction_event_data()
      req(jd)
      
      DT::datatable(
        get_junction_table_for_printing_V2(
          jd$junction_df,
          jd$region_junction_df,
          NMD_table_for_app
        ),
        escape = FALSE,
        options = list(
          scrollX = TRUE,
          
          callback = DT::JS(
            paste(
              related_junction_js(
                "dtable",
                session$ns
              ),
              collapse = "\n"
            )
          ),
          
          headerCallback = DT::JS(
            "function(thead, data, start, end, display){",
            "  $(thead).find('th').eq(1).attr('title','Unique SUPPA event ID (Annotated by SUPPA)');",
            "  $(thead).find('th').eq(2).attr('title','% contribution of LR event coverage to LR region coverage. A high score is indicative for a relevant AS-event annotation and accurate PSI estimation; a medium to low score is indicative for an irrelevant AS-event annotation and/or complex splicing in the region, which may lead to less accurate PSI estimations.');",
            "  $(thead).find('th').eq(3).attr('title','Indicates if event causes NMD: predicted using the splice-context aware method');",
            "  $(thead).find('th').eq(4).attr('title','Indicates if event causes NMD: predicted using the NMD-exclusive method');",
            "  $(thead).find('th').eq(5).attr('title','Junctions occuring in the interfering event but not the primary event (range of mean junction counts per condition in brackets).');",
            "  $(thead).find('th').eq(6).attr('title','Junctions occuring in the primary and interfering event (range of mean junction counts per condition in brackets).');",
            "  $(thead).find('th').eq(7).attr('title','Junctions occuring in the primary event but not the interfering event (range of mean junction counts per condition in brackets).');",
            "  $(thead).find('th').eq(8).attr('title','Click the box to show the junctions of the interfering event in the Faux-Sashimi plot');",
            
            "}"
          )
        ),
        
        colnames = c(
          "Event ID",
          "Relative event coverage (%)",
          "NMD prediction (splice-context aware)",
          "NMD prediction (NMD-exclusive)",
          "Non-shared junctions (interfering event)",
          "Shared junctions (interfering event)",
          "Non-shared junctions (event of interest)",
          "Show interfering event"
        )
      )
    })
    
    ##prints text for clicked junction and event
    output$faux_info <- renderText({
      #req(input$faux_click)
      if (is.null(input$faux_click)) {
        return("CLICK ON A JUNCTION LABEL TO SELECT THE JUNCTION & EVENT")
      }
      
      jd <- junction_event_data()
      req(jd)
      
      # Get currently selected event IDs
      selected <- selected_related_events()
      
      region_junction_df <- jd$region_junction_df %>%
        mutate(
          show_interfering_event = event_id %in% selected
        )
      
      p <- faux_sashimi_plot()
      req(p)
      
      pb <- ggplot_build(p)
      
      x_range <- pb$layout$panel_params[[1]]$x.range
      y_range <- pb$layout$panel_params[[1]]$y.range
      
      # Find the junction closest to the click
      click_selected_row <- return_clicked_junction(
        click = input$faux_click,
        junction_df = jd$junction_df,
        region_junction_df = region_junction_df,
        x_range = x_range,
        y_range = y_range
      )
      
      paste0(
        "CLICK ON A JUNCTION LABEL TO SELECT THE JUNCTION & EVENT", 
        "\nselected junction = ", click_selected_row$junction_range,
        "\nselected (interfering) event = ", click_selected_row$event_id,
        "\nsample condition = ", click_selected_row$cell,
        " treated with ", click_selected_row$treatment,
        "\nmean junction count = ", round(click_selected_row$mean_count, 2),
        "\nstandard deviation = ", round(click_selected_row$sd_count, 2)
      )
    })
    
    #create botton to go to clicked event
    output$goto_event_ui <- renderUI({
      req(input$faux_click)  # show only after click
      actionButton(session$ns("goto_event"), "Go to selected event")
    })
    
    #UCSC genome browser button 
    ##create button
    output$ucsc_button_ui <- renderUI({
      req(selected_event_id())   # only show button if an event is selected
      actionButton(session$ns("go_button"), "View transcripts in event region (UCSC genome browser)", class = "btn btn-primary")
    })
    
    ##outputs message if event not found
    output$message <- renderUI({
      event_msg()
    })
  })
}
