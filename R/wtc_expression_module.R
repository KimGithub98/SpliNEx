wtc_module_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("selected_batches", "Samples to show:",
                         choices = names(batch_list), 
                         selected = names(batch_list)),
      textInput("wtc_gene", "Gene:", "HGNC symbol"),
      actionButton("submit_wtc", "Submit")
    ),
    
    mainPanel(
      plotOutput("wtc_gene_plot"),
      verbatimTextOutput("wtc_gene_output")
    )
  )
}

wtc_module_server <- function(id, normalized_counts_df, colData_df) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$submit_wtc, {
      output$wtc_gene_plot <- renderPlot({
        gene <- input$wtc_gene
        selected_batches <- input$selected_batches
        selected_samples <- unlist(batch_list[selected_batches], use.names = FALSE)
        
        if (!(gene %in% rownames(normalized_counts_df))) return(NULL)
        if (length(selected_samples) == 0) return(NULL)
        
        gene_counts <- t(normalized_counts_df[gene, selected_samples, drop = FALSE])
        metadata <- colData_df[selected_samples, , drop = FALSE]
        df <- cbind(as.data.frame(gene_counts), metadata)
        colnames(df)[1] <- "gene_count"
        df$DIV <- factor(df$DIV, levels = names(div_colors))
        
        ggplot(df, aes(DIV, gene_count)) +
          geom_boxplot(aes(fill = DIV), outlier.shape = NA) +
          geom_point(aes(shape = Data_source, color = treatment), position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), size = 3) +
          scale_fill_manual(values = div_colors) +
          scale_color_manual(values = c("red2", "khaki3", "gray18")) +
          scale_shape_manual(values = c(15,16,17,18,8)) +
          theme_classic() +
          theme(legend.position = "bottom") +
          labs(title = gene, y = "Normalized gene counts")
      })
    })
  })
}