get_PSI_event_data <- function(row){
  PSI <- row[, grepl("_\\d_PSI", colnames(row))]
  PSI_df <- data.frame(cbind(as.numeric(PSI), colnames(row)[grepl("_\\d_PSI", colnames(row))]))
  colnames(PSI_df) <- c("PSI", "sample_name")
  PSI_df <- PSI_df %>% mutate(
    PSI = as.numeric(PSI),
    treatment = str_extract(sample_name, "CHX|DMSO"),
    cell_type = str_extract(sample_name, "iPSC|DIV49|Astro|Micro")
  )
  return(PSI_df)
}

get_PSI_table_for_printing <- function(PSI_df){
  PSI_df <- PSI_df %>% mutate(measurement = stringr::str_extract(sample_name, "\\d_PSI"))  # extracts the trailing number
  
  PSI_wide <- PSI_df %>%
    mutate(PSI = round(PSI, 2)) %>%
    dplyr::select(-sample_name) %>%  # remove original sample_name if desired
    pivot_wider(
      names_from = measurement,
      values_from = PSI,
      names_prefix = "replicate_"
    ) %>%
    mutate(mean_PSI = round(rowMeans(dplyr::select(., starts_with("replicate_")), na.rm = TRUE), 2))
  
  return(PSI_wide)
}

get_PSI_plot <- function(PSI_df, fill_colors){
  p <- ggplot(PSI_df, aes(x = cell_type, y = PSI*100, fill = interaction(cell_type, treatment))) +
    geom_bar(stat = "summary", position = position_dodge(width = 0.8), fun = "mean", size = 1, width = 0.7) +
    geom_point(aes(shape = treatment), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15), size = 2, alpha = 0.8) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    scale_shape_manual(values = c(16,17), name = NULL) +
    theme_classic() +
    coord_cartesian(ylim = c(0,100)) +
    theme(
      strip.text = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    ggtitle("PSI in iPSC-derived cells")
  return(p)
}

get_junction_event_data <- function(event, event_type, row, junction_summary_for_app, strand){
  spacing_factor = 3
  row_spacing = 1.5
  
  junction_df <- junction_summary_for_app[junction_summary_for_app$event_id == event,]
  colnames(junction_df)[colnames(junction_df) == "range"] <- "junction_range"
  
  junction_df <- junction_df %>% rowwise() %>% dplyr::mutate(
    start = min(as.numeric(unlist(str_split(junction_range, "-")))),
    end = max(as.numeric(unlist(str_split(junction_range, "-"))))
  )
  
  #labels
  junction_df <- junction_df  %>% mutate(
    ypos = (-((abs(end - start) / 2) * 0.2)),
    xpos = (end-start)/2+start,
    label = paste(round(mean_count), "±", round(sd_count), sep = "")
  )
  junction_df <- junction_df %>%
    group_by(cell, treatment) %>%
    dplyr::mutate(group_id_sample = dplyr::cur_group_id()*2*row_spacing,
           y_pos = group_id_sample + (as.numeric(junction)/5)) %>%  # assigns a unique number to each group
    ungroup()
  
  
  #faux sashimi plot
  region_junction_df <- junction_summary_for_app[junction_summary_for_app$group_id == junction_df$group_id & 
                                                   junction_summary_for_app$event_id != event &
                                                   !junction_summary_for_app$range %in% junction_df$junction_range & 
                                                   !is.na(junction_summary_for_app$range),]
  colnames(region_junction_df)[colnames(region_junction_df) == "range"] <- "junction_range"
  
  region_junction_df <- region_junction_df %>%
    dplyr::distinct(cell, treatment, junction_range, .keep_all = TRUE)
  
  region_junction_df <- region_junction_df %>% rowwise() %>% dplyr::mutate(
    start = min(as.numeric(unlist(str_split(junction_range, "-")))),
    end = max(as.numeric(unlist(str_split(junction_range, "-"))))
  )
  
  #labels
  region_junction_df <- region_junction_df  %>% mutate(
    ypos = (-((abs(end - start) / 2) * 0.2)),
    xpos = (end-start)/2+start,
    label = paste(round(mean_count), "±", round(sd_count), sep = "")
  )
  region_junction_df <- region_junction_df %>%
    group_by(cell, treatment) %>%
    dplyr::mutate(group_id_sample = dplyr::cur_group_id()*2*row_spacing) %>%  # assigns a unique number to each group
    ungroup()
  
  region_junction_df <- region_junction_df %>%
    mutate(
      y_pos = group_id_sample - (as.numeric(junction)/spacing_factor) - 0.2,
      y_pos = jitter(y_pos, amount = 0.2)   # adds small random vertical spread
    )
  
  #scaling line thickness
  group_max <- junction_df %>%
    group_by(cell, treatment) %>%
    summarise(max_count = max(mean_count, na.rm = TRUE), .groups = "drop")
  
  region_junction_df <- region_junction_df %>%
    left_join(group_max, by = c("cell", "treatment")) %>%
    mutate(log_count_group = log(mean_count / max_count * 100 + 1)) %>%
    dplyr::select(-max_count)
  junction_df <- junction_df %>%
    left_join(group_max, by = c("cell", "treatment")) %>%
    mutate(log_count_group = log(mean_count / max_count * 100 + 1)) %>%
    dplyr::select(-max_count)
  
  #create location exon_df
  if (event_type == "MX"){
    xmin = c(row$pos1-20, row$pos2, row$pos6, row$pos8, row$pos1-20)
    xmax = c(row$pos1 ,row$pos3, row$pos7, row$pos8+20, row$pos8+20)
    intron_exon = c("exon", "exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "yes", "no", "no")
    n_group <- max(unique(junction_df$group_id_sample))
    xmin <- rep(xmin, n_group)
    xmax <- rep(xmax, n_group)
    intron_exon <- rep(intron_exon, n_group)
    alternative_exon = rep(alternative_exon, n_group)
    y <- rep(unique(junction_df$group_id_sample), each = 5)
    exon_df <- as.data.frame(cbind(xmin, xmax, intron_exon, alternative_exon, y))
  } else if ((event_type == "A3" & strand == "+")| (event_type == "A5" & strand == "-")){
    xmin = c(row$pos1-20, row$pos2, row$pos4, row$pos1-20)
    xmax = c(row$pos1 ,row$pos4, row$pos4+20, row$pos4+20)
    intron_exon = c("exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "no", "no")
    n_group <- max(unique(junction_df$group_id_sample))
    xmin <- rep(xmin, n_group)
    xmax <- rep(xmax, n_group)
    intron_exon <- rep(intron_exon, n_group)
    alternative_exon = rep(alternative_exon, n_group)
    y <- rep(unique(junction_df$group_id_sample), each = 4)
    exon_df <- as.data.frame(cbind(xmin, xmax, intron_exon, alternative_exon, y))
  } else if ((event_type == "A3" & strand == "-")| (event_type == "A5" & strand == "+")){
    xmin = c(row$pos3-20, row$pos3, row$pos4, row$pos3-20)
    xmax = c(row$pos3 ,row$pos1, row$pos4+20, row$pos4+20)
    intron_exon = c("exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "no", "no")
    n_group <- max(unique(junction_df$group_id_sample))
    xmin <- rep(xmin, n_group)
    xmax <- rep(xmax, n_group)
    intron_exon <- rep(intron_exon, n_group)
    alternative_exon = rep(alternative_exon, n_group)
    y <- rep(unique(junction_df$group_id_sample), each = 4)
    exon_df <- as.data.frame(cbind(xmin, xmax, intron_exon, alternative_exon, y))
  } else if (event_type == "SE"){
    xmin = c(row$pos1-20, row$pos2, row$pos4, row$pos1-20)
    xmax = c(row$pos1 ,row$pos3, row$pos4+20, row$pos4+20)
    intron_exon = c("exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "no", "no")
    n_group <- max(unique(junction_df$group_id_sample))
    xmin <- rep(xmin, n_group)
    xmax <- rep(xmax, n_group)
    intron_exon <- rep(intron_exon, n_group)
    alternative_exon = rep(alternative_exon, n_group)
    y <- rep(unique(junction_df$group_id_sample), each = 4)
    exon_df <- as.data.frame(cbind(xmin, xmax, intron_exon, alternative_exon, y))
  }
  
  exon_df <- exon_df %>% 
    mutate(
      y = as.numeric(y),
      ymax = case_when(
        intron_exon == "exon"   ~ y,
        intron_exon == "intron" ~ y - 0.08
      ),
      ymin = case_when(
        intron_exon == "exon"   ~ y - 0.20,
        intron_exon == "intron" ~ y - 0.12
      )
    )
  exon_df <- left_join(exon_df, junction_df[,c("group_id_sample", "cell", "treatment")], by = c("y" = "group_id_sample"))
  return(list("junction_df" = junction_df, "region_junction_df"= region_junction_df, "exon_df" = exon_df ))

}

get_faux_sashimi_plot <- function(event_type, junction_df, region_junction_df, x_axis_label, exon_df, fill_colors, show_all_junctions_in_region){
  if (event_type == "RI"){
    p <- ggplot(data = junction_df, aes(color = interaction(cell, treatment))) +
      # Coverage plot
      #geom_ribbon(data = all_df_per_base, aes(x = pos, ymax = cov, ymin = ymin, fill = bam_file, group = bam_file, alpha = 0.3), color = NA) +
      scale_size(range = c(0.5, 2)) +
      labs(x = x_axis_label, y = "Coverage") +
      theme_classic() +
      theme(
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18)
      ) +
      #exons_introns
      annotate("text", x = 1, y = 1, label = "Faux-sashimi cannot be created for RI events", size = 8) +
      ggtitle("Faux-sashimi selected event")
  } else {
    if (show_all_junctions_in_region){
      # Automatically generate colors for all unique event_ids
      unique_events <- unique(region_junction_df$event_id)
      
      # Option 1: Using RColorBrewer (good for up to ~12 events)
      n_colors <- length(unique_events)
      if (n_colors <= 12) {
        event_colors <- brewer.pal(n_colors, "Pastel1")
      } else {
        # fallback to hue palette if more than 12
        event_colors <- hue_pal()(n_colors)
      }
      names(event_colors) <- unique_events
      p <- ggplot(data = junction_df) +
        # Coverage plot
        #geom_ribbon(data = all_df_per_base, aes(x = pos, ymax = cov, ymin = ymin, fill = bam_file, group = bam_file, alpha = 0.3), color = NA) +
        scale_size(range = c(0.5, 2)) +
        labs(x = x_axis_label, y = "Coverage") +
        theme_classic() +
        #arches
        geom_segment(
          aes(x = start, xend = end, y = y_pos, yend = y_pos, size = log_count_group), color = "black"
        ) +
        geom_segment(
          data = junction_df,
          aes(x = start, xend = start, y = y_pos, yend = group_id_sample , size = log_count_group), color = "black"
        ) +
        geom_segment(
          data = junction_df,
          aes(x = end, xend = end, y = y_pos, yend = group_id_sample , size = log_count_group), color = "black"
        ) +
        geom_label(data = junction_df,aes(x = xpos, y = y_pos, label = label), color = "black", fill = "white") +
        #scale_color_manual(values = fill_colors, name = NULL) +
        #related arches
        geom_segment(
          data = region_junction_df,
          aes(x = start, xend = end, y = y_pos, yend = y_pos, size = log_count_group, color = event_id),
          inherit.aes = FALSE, show.legend = FALSE
        ) +
        geom_segment(
          data = region_junction_df,
          aes(x = start, xend = start, y = y_pos, yend = group_id_sample- 0.2 , size = log_count_group, color = event_id),
          inherit.aes = FALSE, show.legend = FALSE
        ) +
        geom_segment(
          data = region_junction_df,
          aes(x = end, xend = end, y = y_pos, yend = group_id_sample- 0.2 , size = log_count_group, color = event_id),
          inherit.aes = FALSE, show.legend = FALSE
        ) +
        geom_label(data = region_junction_df,aes(x = xpos, y = y_pos, label = label, color = event_id), fill = "white", show.legend = FALSE) +
        scale_color_manual(values = c(fill_colors, event_colors)) +  # separate palettes for cell/treatment and event_id
        scale_fill_manual(values = fill_colors, name = NULL) +
        theme(
          strip.text = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank()
        ) +
        geom_text(data = region_junction_df %>% filter(start == min(start)) %>% dplyr::select(start, group_id_sample, cell, treatment) %>% distinct(),  aes(y = group_id_sample, x = as.numeric(start)-30, label = paste(cell, treatment)),
                  hjust = 1,  # 0 = left-aligned, 1 = right-aligned
                  vjust = 0.5, # vertically centered
                  size = 3) +
        #exons_introns
        geom_rect(data = exon_df, aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), fill = interaction(cell, treatment))) + 
        ggtitle("Faux-sashimi selected event") +
        coord_cartesian(clip = "off")
        #coord_cartesian(xlim = c(min(as.numeric(exon_df$xmin), na.rm = TRUE), max(as.numeric(exon_df$xmin), na.rm = TRUE)))
    } else {
      p <- ggplot(data = junction_df) +
        # Coverage plot
        #geom_ribbon(data = all_df_per_base, aes(x = pos, ymax = cov, ymin = ymin, fill = bam_file, group = bam_file, alpha = 0.3), color = NA) +
        scale_size(range = c(0.5, 2)) +
        labs(x = x_axis_label, y = "Coverage") +
        theme_classic() +
        #arches
        geom_segment(
          aes(x = start, xend = end, y = y_pos, yend = y_pos, size = log_count_group, color = interaction(cell, treatment))
        ) +
        geom_segment(
          data = junction_df,
          aes(x = start, xend = start, y = y_pos, yend = group_id_sample , size = log_count_group, color = interaction(cell, treatment))
        ) +
        geom_segment(
          data = junction_df,
          aes(x = end, xend = end, y = y_pos, yend = group_id_sample , size = log_count_group, color = interaction(cell, treatment))
        ) +
        geom_label(data = junction_df,aes(x = xpos, y = y_pos, label = label, color = interaction(cell, treatment)), fill = "white") +
        scale_color_manual(values = fill_colors, name = NULL) +
        scale_fill_manual(values = fill_colors, name = NULL) +
        theme(
          strip.text = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank()
        ) +
        geom_text(data = exon_df %>% filter(xmin == min(xmin), intron_exon == "exon") %>% distinct(), aes(y = ymin, x = as.numeric(xmin), label = paste(cell, treatment), color = interaction(cell, treatment)),
                  hjust = 1,  # 0 = left-aligned, 1 = right-aligned
                  vjust = 0, # vertically centered
                  size = 3) +
        #exons_introns
        geom_rect(data = exon_df, aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), fill = interaction(cell, treatment))) + 
        ggtitle("Faux-sashimi selected event") +
        coord_cartesian(xlim = c(min(as.numeric(exon_df$xmin), na.rm = TRUE), max(as.numeric(exon_df$xmin), na.rm = TRUE)), clip = "off")
    }
    
  }
  return(p)
}

assign_interfering_junction <- function(interfering_junction_range, junction_df){
  interfering_positions <- unlist(str_split(interfering_junction_range, "-"))
  interfering_positions <- paste(interfering_positions, collapse = "|")
  junction_df$interfering_junction_range[grepl(interfering_positions,junction_df$junction_range)] <- interfering_junction_range
  junction_df <- junction_df[grepl(interfering_positions,junction_df$junction_range),]
  return(junction_df)
}

get_junction_table_for_printing <- function(junction_df, region_junction_df){
  junction_df <- junction_df[,c("cell", "treatment", "junction_range", "label")]
  region_junction_df <- region_junction_df[,c("event_id", "cell", "treatment", "junction_range", "label")]
  colnames(region_junction_df) <- c("event_id", "cell", "treatment", "interfering_junction_range", "interfering_junction_count")
  
  junction_df_list <- lapply(region_junction_df$interfering_junction_range, function(x) assign_interfering_junction(x, junction_df))
  junction_df_list <- junction_df_list %>% bind_rows() %>% distinct()
  colnames(junction_df_list) <- c("cell", "treatment", "affected_junction_range", "affected_junction_count", "interfering_junction_range")
  
  region_junction_df <- left_join(region_junction_df, junction_df_list) %>% distinct()
  region_junction_df <- region_junction_df %>%
    mutate(sample_condition = paste(cell, treatment, sep = "_")) %>%
    dplyr::select(-cell, -treatment) %>%
    pivot_wider(
      names_from = sample_condition,
      values_from = c(interfering_junction_count, affected_junction_count),
      names_glue = "{sample_condition}_{.value}"
    )
  #reorder columns
  cols <- grep("count", colnames(region_junction_df), value = TRUE)
  non_count_cols <- grep("count", colnames(region_junction_df), value = TRUE, invert = TRUE)
  samples <- unique(gsub("_(interfering_junction_count|affected_junction_count)$", "", cols))
  count_order <- c()
  for (s in samples) {
    count_order <- c(
      count_order, 
      paste0(s, "_interfering_junction_count"), 
      paste0(s, "_affected_junction_count")
    )
  }
  col_order <- c(non_count_cols, count_order)
  region_junction_df <- region_junction_df %>% dplyr::select(all_of(col_order))
  return(region_junction_df)
}

return_clicked_junction <- function(x_coord, y_coord, junction_df, region_junction_df){
  junction_df$related_junction <- "in_event"
  region_junction_df$related_junction <- "in_related_event"
  
  combined_junction_df <- rbind(junction_df, region_junction_df)
  x_range <- max(combined_junction_df$xpos, na.rm = TRUE) - min(combined_junction_df$xpos, na.rm = TRUE)
  y_range <- max(combined_junction_df$y_pos, na.rm = TRUE) - min(combined_junction_df$y_pos, na.rm = TRUE)
  
  combined_junction_df$dist <- sqrt(
    ((combined_junction_df$xpos - x_coord) / x_range)^2 +
      ((combined_junction_df$y_pos - y_coord) / y_range)^2
  )
  
  nearest_row <- combined_junction_df[which.min(combined_junction_df$dist), ]
  return(nearest_row)
}

get_ucsc_url <- function(event, event_type, chrom){
  if (event_type == "SE" | event_type == "RI"){
    region <- paste(unlist(str_split(event, ";|:|-"))[5], unlist(str_split(event, ";|:|-"))[6], sep = "%2D")
  } else if (event_type == "A3"| event_type == "A5"){
    remove_from_region <- unlist(str_split(event, ";|:|-"))[4:7][duplicated(unlist(str_split(event, ";|:|-"))[4:7])]
    region <- paste(unlist(str_split(event, ";|:|-"))[4:7][unlist(str_split(event, ";|:|-"))[4:7] != remove_from_region], collapse = "%2D")
  } else if (event_type == "MX"){
    region <- paste(unlist(str_split(event, ";|:|-"))[5], unlist(str_split(event, ";|:|-"))[10], sep = "%2D")
  }
  url <- paste("https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=3191311462&position=chr", chrom, ":", region, sep = "")
  return(url)
}
