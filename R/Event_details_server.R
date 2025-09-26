event_details_server <- function(id, selected_event) {
  moduleServer(id, function(input, output, session) {
    
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
    
    # --- single reactive value for the currently selected event ---
    selected_event_id <- reactiveVal(NULL)
    
    # Update reactive when user clicks in NMD tab
    observeEvent(selected_event(), {
      selected_event_id(selected_event())
      updateTextInput(session, "event_input", value = selected_event())
      updateTabsetPanel(session, "main_tabs", selected = "Event details")
    })
    
    # Update reactive when user submits an event manually
    observeEvent(input$submit_event, {
      event <- as.character(input$event_input)
      if (event %in% NMD_table_for_app$event_id) {
        selected_event_id(event)
      }
    })
    
    # --- react to selected_event_id to generate plots/tables ---
    observeEvent(selected_event_id(), {
      event <- selected_event_id()
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      event_type <- row$event_type
      chrom <- row$chrom
      strand <- row$strand
      
      # PSI plot
      PSI_df <- get_PSI_event_data(row)
      output$event_PSI_plot <- renderPlot({ get_PSI_plot(PSI_df, fill_colors) })
      
      # PSI table
      output$event_psi_DT <- DT::renderDT({
        DT::datatable(get_PSI_table_for_printing(PSI_df), escape = FALSE, options = list(scrollX = TRUE))
      })
      
      # Junction plots
      junction_event_data <- get_junction_event_data(event, event_type, row, junction_summary_for_app, strand)
      output$Faux_sashimi_plot <- renderPlot({
        show_all_junctions_in_region <- input$show_related_junctions
        x_axis_label <- paste("Genomic Position (chr", chrom, ")", sep = "")
        get_faux_sashimi_plot(event_type, junction_event_data$junction_df, 
                              junction_event_data$region_junction_df, x_axis_label, 
                              junction_event_data$exon_df, fill_colors, show_all_junctions_in_region)
      })
      
      # Related junctions table
      output$related_junction_DT <- DT::renderDT({
        DT::datatable(get_junction_table_for_printing(junction_event_data$junction_df, junction_event_data$region_junction_df),
                      escape = FALSE, options = list(scrollX = TRUE))
      })
      
      # Click info
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
      
      # Go to selected event button
      output$goto_event_ui <- renderUI({
        req(input$faux_click)  # show only after click
        actionButton(session$ns("goto_event"), "Go to selected event")
      })
      
      observeEvent(input$goto_event, {
        req(input$faux_click)
        click_selected_row <- return_clicked_junction(input$faux_click$x, input$faux_click$y,
                                                      junction_event_data$junction_df,
                                                      junction_event_data$region_junction_df)
        selected_event_id(click_selected_row$event_id)  # update reactive
        updateTextInput(session, "event_input", value = click_selected_row$event_id)
        updateTabsetPanel(session, "main_tabs", selected = "Event details")
      })
      
    })
    
    # UCSC button
    observeEvent(input$go_button, {
      event <- selected_event_id()
      req(event)
      row <- NMD_table_for_app[NMD_table_for_app$event_id == event, ]
      url <- get_ucsc_url(event, row$event_type, row$chrom)
      browseURL(url)
      shinyjs::disable("go_button") 
    })
  })
}
