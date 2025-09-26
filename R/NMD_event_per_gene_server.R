pretty_LR_coverage_labels <- c(
  "1_D_C" = "iNeurons (CHX)",
  "2_D_D" = "iNeurons (DMSO)",
  "3_i_C" = "iPSC (CHX)",
  "4_i_D" = "iPSC (DMSO)",
  "5_A_C" = "Astrocytes (CHX)",
  "6_A_D" = "Astrocytes (DMSO)",
  "7_M_C" = "Microglia (CHX)",
  "8_M_D" = "Microglia (DMSO)"
)

pretty_SR_coverage_labels <- c(
  "iPSC" = "iPSC", 
  "DIV21" = "iNeurons (DIV21)", 
  "DIV28" = "iNeurons (DIV28)",
  "DIV49" = "iNeurons", 
  "Astrocytes" = "Astrocytes", 
  "Microglia" = "Microglia"
)

nmd_event_per_gene_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$submit_nmd, {
      gene <- input$nmd_gene
      selected_batches <- batch_ids[grepl("DMSO|CHX", batch_ids)]
      selected_samples <- unlist(batch_list[selected_batches], use.names = FALSE)
      NMD_table_for_app_gene <- NMD_table_for_app[grepl(gene, NMD_table_for_app$hgnc_name),]
      
      # Show message if gene not found
      if (!(gene %in% NMD_table_for_app$hgnc_name)) {
        output$message <- renderText({"Gene not found in dataset."})
        return(NULL)
      }
      
      #plot LR coverage
      output$nmd_LR_gene_plot <- renderPlot({
        # LR coverage
        long_df <- NMD_table_for_app_gene %>%
          pivot_longer(
            cols = matches("gene_coverage"),
            names_to = "sample",
            values_to = "gene_count"
          ) %>%
          mutate(sample = str_remove(sample, "gene_coverage_"),
                 sample_pretty = pretty_LR_coverage_labels[sample]) 
        
        long_df <- long_df[!grepl("Ensembl|RefSeq", long_df$sample),]
        print(head(long_df))
        ggplot(long_df, aes(x = sample_pretty, y = log10(gene_count+1))) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = log10(10), fill = "red", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(10), ymax = log10(50), fill = "orange", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(50), ymax = log10(200), fill = "yellow", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(200), ymax = Inf, fill = "green", alpha = 0.2) +
          geom_point(size = 3) +
          labs(x = "", y = "log10(LR gene coverage)") +
          theme_classic() +
          theme(legend.position = "bottom",
                plot.title = element_text(hjust = 0.5, size = 20),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Long-read (LR) gene coverage")
      })
      
      #create NMD table to show
      output$NMD_table <- DT::renderDT({
        selected_NMD_event <- input$selected_NMD_event_types
        selected_NMD_only <- input$selected_NMD_only_types
        selected_event_types <- input$selected_event_types
        dt <- NMD_table_for_app_gene[NMD_table_for_app_gene$NMD_sister_transcript %in% selected_NMD_event &
                                       NMD_table_for_app_gene$NMD_exclusive %in% selected_NMD_only &
                                       NMD_table_for_app_gene$event_type %in% selected_event_types,
                                     c("event_type", "event_id", "NMD_sister_transcript", "NMD_exclusive", "hgnc_name", "hgnc_id", "perc_region_avg", "avg_PSI", "sd_PSI")]
        dt$perc_region_avg <- round(dt$perc_region_avg*100)
        dt$avg_PSI <- round(dt$avg_PSI*100, digits = 1)
        dt$sd_PSI <- round(dt$sd_PSI*100, digits = 1)
        dt$event_id <- paste0(
          '<a href="#" onclick="Shiny.setInputValue(\'selected_event\', \'', dt$event_id, '\', {priority: \'event\'});">',
          dt$event_id,
          '</a>'
        )
        
        dt$perc_region_avg <- ifelse(
          dt$perc_region_avg < 10,
          paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:red;" title="Very low contribution to total gene coverage ( <10% ) — likely biologically irrelevant"></span> (', dt$perc_region_avg, "%)"),
          ifelse(
            dt$perc_region_avg < 60,
            paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:orange;" title="Low contribution total gene coverage (10–49%) — potentially minor biological relevance"></span> (', dt$perc_region_avg, "%)"),
            ifelse(
              dt$perc_region_avg < 80,
              paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:gold;" title="Moderate (50–89%)"></span> (', dt$perc_region_avg, "%)"),
              paste0('<span class="glyphicon glyphicon-ok" style="color:green;" title="High (≥90%)"></span> (', dt$perc_region_avg, "%)")
            )
          )
        )
        
        DT::datatable(dt, 
                      escape = FALSE,  # needed to render HTML
                      options = list(scrollX = TRUE))
        #DT::datatable(dt, options = list(scrollX = TRUE))
      })
      
      #plot short read coverage
      output$nmd_SR_gene_plot <- renderPlot({
        # SR coverage
        gene_counts <- t(normalized_counts_df[gene, unique(selected_samples), drop = FALSE])
        
        # Subset metadata to selected samples
        metadata <- colData_df[unique(selected_samples), , drop = FALSE]
        
        # Combine into one dataframe
        gene_count_df <- as.data.frame(gene_counts)
        gene_count_df <- cbind(gene_count_df, metadata)
        colnames(gene_count_df)[1] <- "gene_count"
        gene_count_df$pretty_sample_name <- pretty_SR_coverage_labels[gene_count_df$DIV]
        gene_count_df$pretty_sample_name <- factor(gene_count_df$pretty_sample_name, levels = c("iPSC", "iNeurons (DIV21)", "iNeurons (DIV28)","iNeurons", "Astrocytes", "Microglia"))
        gene_count_df$DIV <- factor(gene_count_df$DIV, levels = c("iPSC", "DIV21", "DIV28","DIV49", "Astrocytes", "Microglia"))
        
        # Plot
        ggplot(gene_count_df, aes(x = pretty_sample_name, y = log10(gene_count)), fill = "ivory2") +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = log10(10), fill = "red", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(10), ymax = log10(100), fill = "orange", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(100), ymax = log10(1000), fill = "yellow", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(1000), ymax = Inf, fill = "green", alpha = 0.2) +
          geom_boxplot(fill = "ivory2", outlier.shape = NA) +
          geom_jitter(aes(shape = treatment),
                      width = 0.2,   
                      height = 0, 
                      size = 3) + # slight transparency so overlaps are visible
          labs(x = "", y = "log10(Normalized gene counts)", shape = "Treatment") +
          theme_classic() +
          theme(legend.position = "bottom",
                plot.title = element_text(hjust = 0.5, size = 20),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Short-read (SR) gene coverage")
      })
    }
    )
  })
}