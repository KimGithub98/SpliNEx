wtc_gene_expression_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    gene_msg <- reactiveVal(NULL)
    
    #reactive plot function
    wtc_gene_plot <- reactive({
      req(input$submit_wtc)
      gene <- toupper(input$wtc_gene)
      selected_batches <- input$selected_batches
      selected_samples <- unlist(batch_list[selected_batches], use.names = FALSE)
      
      gene_msg(NULL)
      
      if (!(gene %in% rownames(normalized_counts_df))) {
        gene_msg(tags$span(style = "color:red;", "Gene not found in dataset. Did you use the HGNC gene symbol?"))
        return()
      }
      
      if (is.null(selected_samples) || length(selected_samples) == 0) {
        gene_msg(tags$span(style = "color:red;", "Please select at least one sample."))
        return()
      }
      
      # Subset gene counts for selected samples
      gene_counts <- t(normalized_counts_df[gene, selected_samples, drop = FALSE])
      
      # Subset metadata to selected samples
      metadata <- colData_df[selected_samples, , drop = FALSE]
      
      # Combine into one dataframe
      gene_count_df <- as.data.frame(gene_counts)
      gene_count_df <- cbind(gene_count_df, metadata)
      colnames(gene_count_df)[1] <- "gene_count"
      gene_count_df$DIV <- factor(
        gene_count_df$DIV,
        levels = c("iPSC", "DIV21", "DIV28", "DIV49", "Astrocytes", "Microglia"),
        labels = c("iPSCs","Neurons (DIV21)","Neurons (DIV28)","Neurons (DIV49)","Astrocytes (DIV41)","Microglia (DIV23)"))
      
      # Plot
      ggplot(gene_count_df, aes(x = DIV, y = gene_count)) +
        geom_boxplot(aes(fill = DIV), outlier.shape = NA) +
        geom_point(aes(shape = Data_source_pretty_name, color = treatment),
                   size = 3,
                   position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
        labs(x = "Cell-type", y = "DESeq2 normalized gene counts", fill = "Cell-type", shape = "Data source", color = "Treatment") +
        scale_shape_manual(values = c(15, 16, 17, 18, 8)) +
        scale_fill_manual(values = div_colors) +
        scale_color_manual(values = c("red2", "khaki3", "gray18")) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(gene)
    })
    
    observeEvent(input$submit_wtc, {
      # reset message
      output$wtc_gene_plot <- renderPlot({
        req(wtc_gene_plot())
        wtc_gene_plot()
      })
    })
    
    
    output$download_gene_exp_plot_png <- downloadHandler(
      filename = function(){
        paste0("Gene_expression_", input$wtc_gene, ".png")
      },
      content = function(file){
        ggsave(file, plot = wtc_gene_plot(), device = "png")
      }
    )
    
    output$download_gene_exp_plot_svg <- downloadHandler(
      filename = function() {
        paste0("Gene_expression_", input$wtc_gene, ".svg")
      },
      content = function(file) {
        # Save SVG with explicit window_title to avoid the warning
        svglite::svglite(
          filename = file,
          width = 6,
          height = 4
        )
        
        print(wtc_gene_plot())  # actually draw the plot into the device
        dev.off()               # close the SVG device
      }
    )
    
    output$plot_ready <- renderText({
      if (!is.null(wtc_gene_plot())) {
        "true"
      } else {
        "false"
      }
    })
    
    outputOptions(output, "plot_ready", suspendWhenHidden = FALSE)
    
    
    output$message <- renderUI({
      gene_msg()
    })
  })
}

