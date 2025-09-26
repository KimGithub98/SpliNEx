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
      NMD_table_for_app_gene <- NMD_table_for_app[NMD_table_for_app$hgnc_name == gene,]
      
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
        long_df$sample_pretty <- factor(long_df$sample_pretty, c("iPSC (DMSO)", "iPSC (CHX)","iNeurons (DMSO)","iNeurons (CHX)",
                                                                 "Astrocytes (DMSO)","Astrocytes (CHX)","Microglia (DMSO)","Microglia (CHX)"))
        
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
      ##explanation coverage color zones
      output$coverage_legend <- renderUI({
        req(input$nmd_gene)  # ensures plot is rendered first
        fluidRow(
          column(
            width = 12,
            wellPanel(
              style = "background-color: #f9f9f9; border: 1px solid #ccc;",
              tags$strong("Explanation of color zones for short-read (SR) / long-read (LR) coverage: "), 
              tags$br(),
              tags$span(style = "color:green;", "Green: High SR/LR coverage — robust detection of isoforms and more reliable PSI estimates."),
              tags$br(),
              tags$span(style = "color:goldenrod;", "Yellow: Medium SR/LR coverage — moderate detection of isoforms and moderately reliable PSI values."),
              tags$br(),
              tags$span(style = "color:orange;", "Orange: Low SR/LR coverage — reduced detection of isoforms and potentially less reliable or missing PSI values."),
              tags$br(),
              tags$span(style = "color:red;", "Red: Very low SR/LR coverage — poor detection of isoforms and potentially less reliable or missing PSI values.")
            )
          )
        )
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
          dt$perc_region_avg < 20,
          paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:red;" title="Low contribution of LR event coverage to LR region coverage (<20%): Indicative for an irrelavent AS-event annotation and inacurate PSI estimation"></span> (', dt$perc_region_avg, "%)"),
          ifelse(
            dt$perc_region_avg < 60,
            paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:orange;" title="Low contribution of LR event coverage to LR region coverage (20–49%): Indicative for complex splicing or a less relavent AS-event annotation and less acurate PSI estimation"></span> (', dt$perc_region_avg, "%)"),
            ifelse(
              dt$perc_region_avg < 80,
              paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:gold;" title="Moderate contribution of LR event coverage to LR region coverage (50–89%): Indicative for complex splicing or a less relavent AS-event annotation and less acurate PSI estimation"></span> (', dt$perc_region_avg, "%)"),
              paste0('<span class="glyphicon glyphicon-ok" style="color:green;" title="High contribution of LR event coverage to LR region coverage (≥90%): Indicative for a relavent AS-event annotation and acurate PSI estimation"></span> (', dt$perc_region_avg, "%)")
            )
          )
        )
        
        DT::datatable(dt, 
                      escape = FALSE,  # needed to render HTML
                      options = list(scrollX = TRUE,
                                     headerCallback = DT::JS(
                                       "function(thead, data, start, end, display){",
                                       "  $(thead).find('th').eq(1).attr('title','Type of alternative splicing event (Annotated by SUPPA)');",
                                       "  $(thead).find('th').eq(2).attr('title','Unique SUPPA event ID (Annotated by SUPPA)');",
                                       "  $(thead).find('th').eq(3).attr('title','Indicates if event causes NMD predicted by the splice-context aware method');",
                                       "  $(thead).find('th').eq(4).attr('title','Indicates if event causes NMD predicted by the NMD-exclusive method');",
                                       "  $(thead).find('th').eq(5).attr('title','HGNC gene symbol');",
                                       "  $(thead).find('th').eq(6).attr('title','HGNC gene ID');",
                                       "  $(thead).find('th').eq(7).attr('title','Relative event coverage: (LR reads supporting the inclusion or exclusion of the event)/(LR reads spanning the event region). A high score is indicative for a relavant AS-event annotation and accurate PSI estimation, a medium to low score is indicative for a irrelavent AS-event annotation and/or complex splicing in the region leading which may lead to less-acurate PSI estimations');",
                                       "  $(thead).find('th').eq(8).attr('title','Mean PSI (%) of all samples (iPSCs, iNeurons, Astrocytes and Microglia treated with CHX or DMSO)');",
                                       "  $(thead).find('th').eq(9).attr('title','Standard deviation of the PSI values(%) of all samples (iPSCs, iNeurons, Astrocytes and Microglia treated with CHX or DMSO');",
                                       "}"
                                     )
                                     ),
                      colnames = c(
                        "Event type",
                        "Event ID",
                        "NMD prediction (splice-context aware)",
                        "NMD prediction (NMD-exclusive)",
                        "Gene symbol (HGNC)",
                        "HGNC ID",
                        "Relative event coverage (%)",
                        "Average PSI (%)",
                        "Standard deviation PSI (%)"
                      ))
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