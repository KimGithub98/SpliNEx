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
    gene_msg <- reactiveVal(NULL)
    observeEvent(input$submit_nmd, {
      gene <- toupper(input$nmd_gene)
      selected_batches <- batch_ids[grepl("DMSO|CHX", batch_ids)]
      selected_samples <- unlist(batch_list[selected_batches], use.names = FALSE)
      NMD_table_for_app$hgnc_name <- toupper(NMD_table_for_app$hgnc_name)
      NMD_table_for_app_gene <- NMD_table_for_app[NMD_table_for_app$hgnc_name == gene,]
      
      # reset message
      gene_msg(NULL)
      if (!(gene %in% NMD_table_for_app$hgnc_name) & !(gene %in% rownames(normalized_counts_df))) {
        gene_msg(tags$span(style = "color:red;", "Gene not found in dataset. Did you use the HGNC gene symbol?"))
        return()
      } else if (!(gene %in% NMD_table_for_app$hgnc_name) & (gene %in% rownames(normalized_counts_df))) {
        gene_msg(
          tags$span(
            style = "color:red;",
            tags$strong("No alternative splicing data available for the selected gene in SpliNEx."),
            "This is likely due to low long-read coverage and/or filtering thresholds. ",
            "Additional information may still be present in the unfiltered data available from ",
            tags$a(
              href = "https://doi.org/your-doi-here", 
              "[Wijnant et al., 2025]", 
              target = "_blank",
              style = "color:red; text-decoration: underline;"
            ),
            "."
          )
        )
        return()
      }
      
      output$gene_title <- renderUI({
        req(gene)  # Ensure a gene is selected
        h3(
          paste("Alternative splice events for", gene),
          style = "text-align: center; margin-bottom: 20px;"
        )
      })
      
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
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = log10(10), fill = "darkorange3", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(10), ymax = log10(50), fill = "goldenrod1", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(50), ymax = log10(200), fill = "cadetblue1", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(200), ymax = Inf, fill = "cadetblue", alpha = 0.2) +
          geom_point(size = 3, na.rm = TRUE) +
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
              tags$span(style = "color:#5F9EA0;", "Blue: High SR/LR coverage — Likely robust detection of isoforms and more reliable PSI estimates."),
              tags$br(),
              tags$span(style = "color:#BFDFE8;", "Light blue: Medium SR/LR coverage — Likely moderate detection of isoforms and moderately reliable PSI values."),
              tags$br(),
              tags$span(style = "color:#FFC125;", "Yellow: Low SR/LR coverage — Likely reduced detection of isoforms and potentially less reliable or missing PSI values."),
              tags$br(),
              tags$span(style = "color:#CD6600;", "Dark orange: Very low SR/LR coverage — Likely poor detection of isoforms and potentially less reliable or missing PSI values.")
            )
          )
        )
      })
      
      #function to create table
      nmd_table_filtered <- reactive({
        req(input$nmd_gene)
        
        selected_NMD_event <- input$selected_NMD_event_types
        selected_NMD_only <- input$selected_NMD_only_types
        selected_event_types <- input$selected_event_types
        dt <- NMD_table_for_app_gene[NMD_table_for_app_gene$NMD_sister_transcript %in% selected_NMD_event &
                                       NMD_table_for_app_gene$NMD_exclusive %in% selected_NMD_only &
                                       NMD_table_for_app_gene$event_type %in% selected_event_types,
                                     c("event_type", "event_id", "NMD_sister_transcript", "NMD_exclusive", "hgnc_name", "hgnc_id", "perc_region_avg", "avg_PSI", "sd_PSI")]
        if (nrow(dt) == 0) {
          gene_msg(tags$span(style = "color:red;", "No events found in dataset for your selected options."))
          return()
        }
        dt$perc_region_avg <- round(dt$perc_region_avg*100)
        dt$avg_PSI <- round(dt$avg_PSI*100, digits = 1)
        dt$sd_PSI <- round(dt$sd_PSI*100, digits = 1)
        dt
        #DT::datatable(dt, options = list(scrollX = TRUE))
      })
      
      #function_to_create_extended table
      nmd_table_filtered_extended <- reactive({
        req(input$nmd_gene)
        
        selected_NMD_event <- input$selected_NMD_event_types
        selected_NMD_only <- input$selected_NMD_only_types
        selected_event_types <- input$selected_event_types
        dt <- NMD_table_for_app_gene[NMD_table_for_app_gene$NMD_sister_transcript %in% selected_NMD_event &
                                       NMD_table_for_app_gene$NMD_exclusive %in% selected_NMD_only &
                                       NMD_table_for_app_gene$event_type %in% selected_event_types,]
        cols_to_remove <- paste0("pos", 1:8)
        dt <- dt %>% dplyr::select(-all_of(cols_to_remove))
        
        dt$perc_region_avg <- round(dt$perc_region_avg*100)
        dt$avg_PSI <- round(dt$avg_PSI*100, digits = 1)
        dt$sd_PSI <- round(dt$sd_PSI*100, digits = 1)
        colnames(dt) <- c("event_type","chrom","strand","event_id","NMD_sister_transcript","NMD_exclusive",
                          "hgnc_name","hgnc_id","relative_event_coverage_iNeuron_CHX","relative_event_coverage_iNeuron_DMSO","relative_event_coverage_iPSC_CHX",
                          "relative_event_coverage_4_iPSC_DMSO","relative_event_coverage_Astrocyte_CHX","relative_event_coverage_Astrocyte_DMSO","relative_event_coverage_Microglia_CHX",
                          "relative_event_coverage_Microglia_DMSO","mean_relative_event_coverage","iPSC_CHX_1_PSI","iPSC_CHX_2_PSI","iPSC_CHX_3_PSI",
                          "iPSC_DMSO_1_PSI","iPSC_DMSO_2_PSI","iPSC_DMSO_3_PSI","iNeuron_CHX_1_PSI","iNeuron_CHX_2_PSI",
                          "iNeuron_CHX_3_PSI","iNeuron_DMSO_1_PSI","iNeuron_DMSO_2_PSI","iNeuron_DMSO_3_PSI","Astro_CHX_1_PSI",
                          "Astro_CHX_2_PSI","Astro_CHX_3_PSI","Astro_DMSO_1_PSI","Astro_DMSO_2_PSI","Astro_DMSO_3_PSI",
                          "Micro_CHX_1_PSI","Micro_CHX_2_PSI","Micro_CHX_3_PSI","Micro_DMSO_1_PSI","Micro_DMSO_2_PSI",
                          "Micro_DMSO_3_PSI","iNeuron_PSI_CHX","iNeuron_PSI_DMSO","iPSC_PSI_CHX","iPSC_PSI_DMSO","Astro_PSI_CHX",
                          "Astro_PSI_DMSO","Micro_PSI_CHX","Micro_PSI_DMSO","gene_coverage_Ensembl",
                          "gene_coverage_RefSeq","gene_coverage_Microglia_DMSO","gene_coverage_iPSC_DMSO","gene_coverage_Microglia_CHX",
                          "gene_coverage_Astrocyte_CHX","gene_coverage_iPSC_CHX","gene_coverage_Astrocyte_DMSO","gene_coverage_iNeuron_CHX",
                          "gene_coverage_iNeuron_DMSO","avg_PSI","sd_PSI", "NMD_sister_transcript_summary", "NMD_exclusive_summary")
        
        
        dt
        #DT::datatable(dt, options = list(scrollX = TRUE))
      })
      
      #create NMD table to show
      output$NMD_table <- DT::renderDT({
        req(nmd_table_filtered())
        dt <- nmd_table_filtered()
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
                                       "  $(thead).find('th').eq(3).attr('title','Indicates if event causes NMD: predicted using the splice-context aware method');",
                                       "  $(thead).find('th').eq(4).attr('title','Indicates if event causes NMD: predicted using the NMD-exclusive method');",
                                       "  $(thead).find('th').eq(5).attr('title','HGNC gene symbol');",
                                       "  $(thead).find('th').eq(6).attr('title','HGNC gene ID');",
                                       "  $(thead).find('th').eq(7).attr('title','% contribution of LR event coverage to LR region coverage. A high score is indicative for a relavant AS-event annotation and accurate PSI estimation, a medium to low score is indicative for a irrelavent AS-event annotation and/or complex splicing in the region leading which may lead to less-acurate PSI estimations');",
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
        ggplot(gene_count_df, aes(x = pretty_sample_name, y = log10(gene_count)), fill = "#cfcece") +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = log10(10), fill = "darkorange3", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(10), ymax = log10(100), fill = "goldenrod1", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(100), ymax = log10(1000), fill = "cadetblue1", alpha = 0.2) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(1000), ymax = Inf, fill = "cadetblue", alpha = 0.2) +
          geom_boxplot(fill = "#cfcece", outlier.shape = NA, na.rm = TRUE) +
          geom_jitter(aes(shape = treatment),
                      width = 0.2,   
                      height = 0, 
                      size = 3, na.rm = TRUE) + # slight transparency so overlaps are visible
          labs(x = "", y = "log10(Normalized gene counts)", shape = "Treatment") +
          theme_classic() +
          theme(legend.position = "bottom",
                plot.title = element_text(hjust = 0.5, size = 20),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Short-read (SR) gene coverage")
      })
      
      #download short table
      output$download_nmd_event_per_gene_table <- downloadHandler(
        filename = function() {
          paste0(input$nmd_gene, "_NMD_events.csv")
        },
        content = function(file) {
          df <- nmd_table_filtered()
          req(df)
          if ("event_id" %in% names(df)) {
            df$event_id <- gsub("<.*?>", "", df$event_id)   # remove any tags
          }
          write.csv(df, file, row.names = FALSE)
          }
        )
      
      #download extended table
      output$download_nmd_event_per_gene_full_table <- downloadHandler(
        filename = function() {
          paste0(input$nmd_gene, "_AS_events_extended_table.csv")
        },
        content = function(file) {
          df <- nmd_table_filtered_extended()
          req(df)
          if ("event_id" %in% names(df)) {
            df$event_id <- gsub("<.*?>", "", df$event_id)   # remove any tags
          }
          write.csv(df, file, row.names = FALSE)
        }
      )
      
      
      output$data_ready <- renderText({
        if (!is.null(nmd_table_filtered()) & !is.null(nmd_table_filtered_extended())) {
          "true"
        } else {
          "false"
        }
      })
      
      outputOptions(output, "data_ready", suspendWhenHidden = FALSE)
      
    })
    
    output$message_NMD_event_per_gene <- renderUI({
      gene_msg()
    })
  })
}