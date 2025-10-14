event_details_server <- function(id, selected_event) {
  moduleServer(id, function(input, output, session) {
    event_msg <- reactiveVal(NULL)
    selected_event_id <- reactiveVal(NULL)
    clicked_junction <- reactiveVal(NULL)
    fill_colors <- c(
      "iPSC.DMSO" = "ivory2",
      "iPSC.CHX" = colorspace::darken("ivory2", 0.3),
      "DIV49.DMSO" = "cadetblue",
      "DIV49.CHX" = colorspace::darken("cadetblue", 0.3),
      "Astro.DMSO" = "palevioletred",
      "Astro.CHX" = colorspace::darken("palevioletred", 0.3),
      "Micro.DMSO" = "burlywood4",
      "Micro.CHX" = colorspace::darken("burlywood4", 0.3)
    )
    
    
    # Update reactive when user clicks in NMD tab
    observeEvent(selected_event(), {
      event_msg(NULL)
      selected_event_id(selected_event())
      updateTextInput(session, "event_input", value = selected_event())
      updateTabsetPanel(session, "main_tabs", selected = "Event details")
    })
    
    # Update reactive when user submits an event manually
    observeEvent(input$submit_event, {
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
    
    ## Update reactive when user clicks junction in related junction plot
    observeEvent(input$goto_event, {
      req(clicked_junction())
      selected_event_id(clicked_junction())  # update reactive
      updateTextInput(session, "event_input", value = clicked_junction())
      updateTabsetPanel(session, "main_tabs", selected = "Event details")
    })
    
    # --- react to selected_event_id to generate plots/tables ---
    observeEvent(selected_event_id(), {
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      
      event_type <- row$event_type
      chrom <- row$chrom
      strand <- row$strand
      print(row)
      # event details
      output$event_details_print <- renderUI({
        fluidRow(
          column(
            width = 12,
            wellPanel(
              style = "background-color: #f9f9f9; border: 1px solid #ccc;",
              tags$strong("NMD-prediction (splice-context aware method): "), 
              tags$br(),
              tags$span(row$NMD_sister_transcript_summary),
              tags$br(),
              tags$span("Conclusion: "),
              tags$strong(row$NMD_sister_transcript),
              tags$br(),
              tags$br(),
              tags$strong("NMD-prediction (NMD-exclusive method): "), 
              tags$br(),
              tags$span(row$NMD_exclusive_summary),
              tags$br(),
              tags$span("Conclusion: "),
              tags$strong(row$NMD_exclusive),
            )
          )
        )
      })
      
      # PSI plot
      PSI_df <- get_PSI_event_data(row)
      output$event_PSI_plot <- renderPlot({ get_PSI_plot(PSI_df, fill_colors) })
      
      # PSI table
      output$event_psi_DT <- DT::renderDT({
        DT::datatable(get_PSI_table_for_printing(PSI_df), escape = FALSE, options = list(scrollX = TRUE))
      })
      
      # Junction plots
      junction_event_data <- get_junction_event_data(event, event_type, row, junction_summary_for_app, strand) ##gets junction data
      ##gets sashimi plot
      output$Faux_sashimi_plot <- renderPlot({
        req(selected_event_id())
        show_all_junctions_in_region <- input$show_related_junctions ##checks if to show related junctions
        x_axis_label <- paste("Genomic Position (chr", chrom, ")", sep = "")
        get_faux_sashimi_plot(event_type, junction_event_data$junction_df, 
                              junction_event_data$region_junction_df, x_axis_label, 
                              junction_event_data$exon_df, fill_colors, show_all_junctions_in_region)
      })
      ##gives checkbox to view related junctions
      output$checkbox_related_junctions <- renderUI({
        req(selected_event_id())  # show only after click
        checkboxInput(session$ns("show_related_junctions"), "Show related junctions in region", FALSE)
      })
      
      # Related junctions table
      output$related_junction_DT <- DT::renderDT({
        DT::datatable(get_junction_table_for_printing(junction_event_data$junction_df, junction_event_data$region_junction_df),
                      escape = FALSE, options = list(scrollX = TRUE))
      })
      
      # Clickable plot to select events and junctions
      ##creates text with information on clicked junction
      output$faux_info <- renderText({
        req(input$faux_click)
        click_selected_row <- return_clicked_junction(input$faux_click$x, input$faux_click$y,
                                                      junction_event_data$junction_df,
                                                      junction_event_data$region_junction_df)
        paste0("CLICK ON LABEL TO SELECT JUNCTION", 
               "\nselected junction = ", click_selected_row$junction_range,
               "\nselected event = ", click_selected_row$event_id,
               "\nsample condition = ", click_selected_row$cell, " treated with ", click_selected_row$treatment,
               "\nmean junction count = ", round(click_selected_row$mean_count,2),
               "\nstandard deviation = ", round(click_selected_row$sd_count,2))
      })
      ## creates go-to-clicked event button
      output$goto_event_ui <- renderUI({
        req(input$faux_click)  # show only after click
        actionButton(session$ns("goto_event"), "Go to selected event")
      })
      ## Update reactive when user clicks junction in related junction plot
      observeEvent(input$goto_event, {
        req(input$faux_click)
        click_selected_row <- return_clicked_junction(input$faux_click$x, input$faux_click$y,
                                                      junction_event_data$junction_df,
                                                      junction_event_data$region_junction_df)
        clicked_junction(click_selected_row$event_id)
      })
      
    })
    
    #UCSC genome browser button 
    ##create button
    output$ucsc_button_ui <- renderUI({
      req(selected_event_id())   # only show button if an event is selected
      actionButton(session$ns("go_button"), "View transcripts in event region (UCSC genome browser)", class = "btn btn-primary")
    })
    ##links click on button to UCSC genome browser
    observeEvent(input$go_button, {
      req(selected_event_id())
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      url <- get_ucsc_url(event, row$event_type, row$chrom)
      #browseURL(url)
      shinyjs::disable("go_button") 
    })
    ##outputs message if event not found
    output$message <- renderUI({
      event_msg()
    })
  })
}
