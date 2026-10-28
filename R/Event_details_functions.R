get_PSI_event_data <- function(row){
  PSI <- row[, grepl("_\\d_PSI", colnames(row))]
  PSI_df <- data.frame(cbind(as.numeric(PSI), colnames(row)[grepl("_\\d_PSI", colnames(row))]))
  colnames(PSI_df) <- c("PSI", "sample_name")
  PSI_df <- PSI_df %>% mutate(
    PSI = as.numeric(PSI),
    treatment = str_extract(sample_name, "CHX|DMSO"),
    cell_type = str_extract(sample_name, "iPSC|DIV49|Astro|Micro")
  )
  PSI_df$cell_type <- str_replace(PSI_df$cell_type, "DIV49", "Neurons")
  PSI_df$cell_type <- factor(
    PSI_df$cell_type,
    levels = c("iPSC", "Neurons", "Astro", "Micro")
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
    scale_fill_manual(values = fill_colors, breaks = names(fill_colors), name = NULL) +
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
    labs(x = NULL, y = "PSI * 100")
    ggtitle("PSI in iPSC-derived cells")
  return(p)
}

get_perc_region_plot <- function(row){
  row %>%
    pivot_longer(cols = starts_with("perc_region_") & !ends_with("avg"), names_to = "condition", values_to = "percentage") %>%
    tidyr::separate(condition, into = c("perc", "region", "sample_number", "cell_type", "condition"), sep = "_") %>%
    mutate(cell_type = recode(cell_type, D = "Neurons", i = "iPSCs", A = "Astrocytes", M = "Microglia"),
           percentage_color = cut(percentage * 100, c(-Inf, 20, 60, 80, Inf), labels = c("<20%", "20–59%", "60–79%", "≥80%"))) %>%
    ggplot(aes(cell_type, percentage * 100, colour = percentage_color, shape = condition)) +
    geom_point(size = 4, position = position_jitter(width = 0.2, height = 0)) +
    scale_color_manual(values = c("<20%" = "red", "20–59%" = "orange", "60–79%" = "gold", "≥80%" = "green"), guide = "none") +
    scale_shape_manual(values = c(C = 16, D = 17), guide = "none") +
    theme_classic() +
    theme(strip.text = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 18)) +
    labs(x = NULL, y = "Relative event coverage") +
    ggtitle("Relative event covarage in iPSC-derived cells") + 
    coord_cartesian(ylim = c(0,100))
}

get_unique_positions <- function(input_df) {
  #get position from junction df
  lapply(1:8, function(i) {
    unique(input_df[[paste0("pos", i)]])
  }) |>
    setNames(paste0("pos", 1:8))
}

get_exon_df <- function(strand, input_df, event_id_name){
  input_df <- input_df[input_df$event_id == event_id_name,]
  
  #create location exon_df
  pos1 <- unique(input_df$pos1)
  pos2 <- unique(input_df$pos2)
  pos3 <- unique(input_df$pos3)
  pos4 <- unique(input_df$pos4)
  pos5 <- unique(input_df$pos5)
  pos6 <- unique(input_df$pos6)
  pos7 <- unique(input_df$pos7)
  pos8 <- unique(input_df$pos8)
  
  group_id_samples <- input_df$group_id_sample
  
  event_type <- str_extract(event_id_name, "A3|A5|MX|SE|RI")
  
  if (event_type == "MX"){
    xmin = c(pos1-20, pos2, pos6, pos8, pos1-20)
    xmax = c(pos1 ,pos3, pos7, pos8+20, pos8+20)
    intron_exon = c("exon", "exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "yes", "no", "no")
    i_groups = 5
  } else if ((event_type == "A3" & strand == "+")| (event_type == "A5" & strand == "-")){
    xmin = c(pos1-20, pos2, pos4, pos1-20)
    xmax = c(pos1 ,pos4, pos4+20, pos4+20)
    intron_exon = c("exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "no", "no")
    i_groups = 4
  } else if ((event_type == "A3" & strand == "-")| (event_type == "A5" & strand == "+")){
    xmin = c(pos3-20, pos3, pos4, pos3-20)
    xmax = c(pos3 ,pos1, pos4+20, pos4+20)
    intron_exon = c("exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "no", "no")
    i_groups = 4
  } else if (event_type == "SE"){
    xmin = c(pos1-20, pos2, pos4, pos1-20)
    xmax = c(pos1 ,pos3, pos4+20, pos4+20)
    intron_exon = c("exon", "exon", "exon", "intron")
    alternative_exon = c("no", "yes", "no", "no")
    i_groups = 4
  }
  
  n_group <- max(unique(group_id_samples))
  xmin <- rep(xmin, n_group)
  xmax <- rep(xmax, n_group)
  intron_exon <- rep(intron_exon, n_group)
  alternative_exon = rep(alternative_exon, n_group)
  y <- rep(unique(group_id_samples), each = i_groups)
  exon_df <- as.data.frame(cbind(xmin, xmax, intron_exon, alternative_exon, y))
  
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
  
  return(exon_df)
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
  junction_df$cell <- str_replace(junction_df$cell, "DIV49", "Neurons")
  
  junction_df <- junction_df %>%
    mutate(
      cell = factor(cell,levels = c("iPSC", "Micro", "Neurons", "Astro")),
      treatment = factor(treatment,levels = c("CHX", "DMSO"))) %>%
    group_by(cell, treatment) %>%
    dplyr::mutate(group_id_sample = dplyr::cur_group_id()*2*row_spacing,
           y_pos = group_id_sample + (as.numeric(junction)/5)) %>%  # assigns a unique number to each group
    ungroup()
  
  
  #faux sashimi plot
  region_junction_df <- junction_summary_for_app[junction_summary_for_app$group_id == junction_df$group_id & 
                                                   junction_summary_for_app$event_id != event &
                                                   #!junction_summary_for_app$range %in% junction_df$junction_range & 
                                                   !is.na(junction_summary_for_app$range),]
  #test
  #region_junction_df <- region_junction_df[junction_summary_for_app$event_id == unique(junction_summary_for_app$event_id)[1],]
  colnames(region_junction_df)[colnames(region_junction_df) == "range"] <- "junction_range"
  
  region_junction_df <- region_junction_df %>%
    dplyr::distinct(event_id, cell, treatment, junction_range, .keep_all = TRUE)
  
  region_junction_df <- region_junction_df %>% rowwise() %>% dplyr::mutate(
    start = min(as.numeric(unlist(str_split(junction_range, "-")))),
    end = max(as.numeric(unlist(str_split(junction_range, "-"))))
  )
  
  region_junction_df$cell <- str_replace(region_junction_df$cell, "DIV49", "Neurons")
  
  #labels
  region_junction_df <- region_junction_df  %>% mutate(
    ypos = (-((abs(end - start) / 2) * 0.2)),
    xpos = (end-start)/2+start,
    label = paste(round(mean_count), "±", round(sd_count), sep = "")
  )
  region_junction_df <- region_junction_df %>%
    mutate(
      cell = factor(cell,levels = c("iPSC", "Micro", "Neurons", "Astro")),
      treatment = factor(treatment,levels = c("CHX", "DMSO"))) %>%
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
    mutate(log_count_group = log(mean_count / max_count * 100 + 1),
           log_count_group = ifelse(is.finite(log_count_group), log_count_group, NA_real_)) %>%
    dplyr::select(-max_count)
  
  region_junction_df$show_bool <- FALSE
  
  junction_df <- junction_df %>%
    left_join(group_max, by = c("cell", "treatment")) %>%
    mutate(log_count_group = log(mean_count / max_count * 100 + 1),
           log_count_group = ifelse(is.finite(log_count_group), log_count_group, NA_real_)) %>%
    dplyr::select(-max_count)
  
  exon_df <- get_exon_df(strand, input_df = junction_df, event_id_name = event)
  exon_df <- left_join(exon_df, junction_df[,c("group_id_sample", "cell", "treatment")], by = c("y" = "group_id_sample"), relationship = "many-to-many")
  
  if (nrow(region_junction_df) > 0){
    exon_df_region <- data.frame()
    for (event_id_name in unique(region_junction_df$event_id)){
      exon_df_region_event <- get_exon_df(strand, input_df = region_junction_df, event_id_name = event_id_name)
      exon_df_region_event$event_id <- event_id_name
      exon_df_region <- rbind(exon_df_region, exon_df_region_event)
    }
    exon_df_region <- left_join(exon_df_region, region_junction_df[,c("group_id_sample", "cell", "treatment", "event_id")], by = c("y" = "group_id_sample", "event_id" = "event_id"), relationship = "many-to-many")
  } else {
    exon_df_region <- data.frame(
      xmin = numeric(),
      xmax = numeric(),
      intron_exon = character(),
      alternative_exon = character(),
      y = numeric(),
      ymax = numeric(),
      ymin = numeric(),
      event_id = character(),
      cell = character(),
      treatment = character()
    )
  }

  return(list("junction_df" = junction_df, "region_junction_df"= region_junction_df, "exon_df" = exon_df, "exon_df_region" = exon_df_region))
}

get_faux_sashimi_plot <- function(event_type, junction_df, region_junction_df, x_axis_label, exon_df, fill_colors, show_all_junctions_in_region){
  if (show_all_junctions_in_region){
    # Automatically generate colors for all unique event_ids
    unique_events <- unique(region_junction_df$event_id)
    
    region_junction_df <- region_junction_df %>%
      dplyr::distinct(cell, treatment, junction_range, .keep_all = TRUE)
    
    # Option 1: Using RColorBrewer (good for up to ~12 events)
    n_colors <- length(unique_events)
    if (n_colors <= 12) {
      event_colors <- brewer.pal(n_colors, "Pastel1")
    } else {
      # fallback to hue palette if more than 12
      event_colors <- hue_pal()(n_colors)
    }
    names(event_colors) <- unique_events
    p <- ggplot(data = junction_df, aes(x = start, y = y_pos)) +
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
      geom_label(data = junction_df, aes(x = xpos, y = y_pos, label = label), color = "black", fill = "white") +
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
                size = 3,
                inherit.aes = FALSE) +
      #exons_introns
      geom_rect(data = exon_df, 
                aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), fill = interaction(cell, treatment)),
                inherit.aes = FALSE) + 
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
  return(p)
}

get_faux_sashimi_plot_v2 <- function(event_type, junction_df, region_junction_df, x_axis_label, exon_df, exon_df_region, fill_colors){
  # Automatically generate colors for all unique event_ids
  
  region_junction_df <- region_junction_df[region_junction_df$show_interfering_event,]
  exon_df_region <- exon_df_region[exon_df_region$show_interfering_event,]
  unique_events <- unique(region_junction_df$event_id)
  
  #get color pallete
  n_colors <- length(unique_events)
  if (n_colors <= 8) {
    # Option 1: Using RColorBrewer (for up to ~8 events)
    event_colors <- brewer.pal(n_colors, "Dark2")
  } else {
    # Option 2: fallback to hue palette if more than 8
    event_colors <- hue_pal()(n_colors)
  }
  names(event_colors) <- unique_events
  
  combined_junction_df <- rbind(region_junction_df %>% select(-show_bool, -show_interfering_event), junction_df)
  
  combined_junction_df <- combined_junction_df %>% filter(start == min(start, na.rm = TRUE) & !is.na(start)) %>% dplyr::select(start, group_id_sample, cell, treatment) %>% distinct()
  
  if(nrow(region_junction_df) == 0){
    p <- ggplot(data = junction_df, aes(x = start, y = y_pos)) +
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
        aes(x = start, xend = start, y = y_pos, yend = group_id_sample , size = log_count_group, color =interaction(cell, treatment))
      ) +
      geom_segment(
        data = junction_df,
        aes(x = end, xend = end, y = y_pos, yend = group_id_sample , size = log_count_group, color = interaction(cell, treatment))
      ) +
      geom_label(data = junction_df, aes(x = xpos, y = y_pos, label = label, color = interaction(cell, treatment)), fill = "white") +
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
      scale_fill_manual(values = c(fill_colors, event_colors), name = NULL) +
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
      geom_text(data = combined_junction_df,  aes(y = group_id_sample, x = as.numeric(start)-30, label = paste(cell, treatment)),
                hjust = 1,  # 0 = left-aligned, 1 = right-aligned
                vjust = 0.5, # vertically centered
                size = 3,
                inherit.aes = FALSE) +
      #exons_introns
      geom_rect(data = exon_df %>% filter(intron_exon == "intron"), 
                aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), fill = interaction(cell, treatment)),
                inherit.aes = FALSE) + 
      geom_rect(data = exon_df %>% filter(intron_exon == "exon"), 
                aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), fill = interaction(cell, treatment)),
                inherit.aes = FALSE) + 
      geom_rect(data = exon_df_region %>% filter(intron_exon == "exon"), 
                aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), color= event_id), fill = NA,
                inherit.aes = FALSE) +
      ggtitle("Faux-sashimi selected event") +
      coord_cartesian(clip = "off")
  } else if (nrow(region_junction_df) > 0){
    p <- ggplot(data = junction_df, aes(x = start, y = y_pos)) +
      # Coverage plot
      #geom_ribbon(data = all_df_per_base, aes(x = pos, ymax = cov, ymin = ymin, fill = bam_file, group = bam_file, alpha = 0.3), color = NA) +
      scale_size(range = c(0.5, 2)) +
      labs(x = x_axis_label, y = "Coverage") +
      theme_classic() +
      #arches
      geom_segment(
        aes(x = start, xend = end, y = y_pos, yend = y_pos, size = log_count_group), color = "#cfcece"
      ) +
      geom_segment(
        data = junction_df,
        aes(x = start, xend = start, y = y_pos, yend = group_id_sample , size = log_count_group), color = "#cfcece"
      ) +
      geom_segment(
        data = junction_df,
        aes(x = end, xend = end, y = y_pos, yend = group_id_sample , size = log_count_group), color = "#cfcece"
      ) +
      geom_label(data = junction_df, aes(x = xpos, y = y_pos, label = label), color = "#cfcece", fill = "white") +
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
      scale_fill_manual(values = c(fill_colors, event_colors), name = NULL) +
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
      geom_text(data = combined_junction_df,  aes(y = group_id_sample, x = as.numeric(start)-30, label = paste(cell, treatment)),
                hjust = 1,  # 0 = left-aligned, 1 = right-aligned
                vjust = 0.5, # vertically centered
                size = 3,
                inherit.aes = FALSE) +
      #exons_introns
      geom_rect(data = exon_df %>% filter(intron_exon == "intron"), 
                aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax)), fill = "#cfcece",
                inherit.aes = FALSE) + 
      geom_rect(data = exon_df %>% filter(intron_exon == "exon"), 
                aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax)), fill = "#cfcece",
                inherit.aes = FALSE) + 
      geom_rect(data = exon_df_region %>% filter(intron_exon == "exon"), 
                aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), color= event_id), fill = NA,
                inherit.aes = FALSE) +
      ggtitle("Faux-sashimi selected event") +
      coord_cartesian(clip = "off")
  }
  
  #coord_cartesian(xlim = c(min(as.numeric(exon_df$xmin), na.rm = TRUE), max(as.numeric(exon_df$xmin), na.rm = TRUE)))
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

shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
  }
  inputs
}

related_junction_js <- function(dtid, ns) {
  
  input_id <- ns(sprintf("%s_related_junction_check", dtid))
  
  sprintf(
    "
    $(document).off(
      'change.related_junction',
      'input[id^=\"checkb_related_junctions_\"]'
    );
    
    $(document).on(
      'change.related_junction',
      'input[id^=\"checkb_related_junctions_\"]',
      function() {
        
        var event_id = this.id.replace(
          'checkb_related_junctions_',
          ''
        );
        
        var value = $(this).prop('checked');
        
        console.log('RELATED JUNCTION CLICK');
        console.log('event_id:', event_id);
        console.log('value:', value);
        
        Shiny.setInputValue(
          '%s',
          [{
            event_id: event_id,
            value: value
          }],
          {priority: 'event'}
        );
      }
    );
    ",
    input_id
  )
}

get_junction_table_for_printing_V2 <- function(junction_df, region_junction_df, NMD_table_for_app){
  junction_df <- junction_df[,c("cell", "treatment", "junction_range", "mean_count")]
  junction_df_sum <- junction_df %>% 
    filter(!is.na(junction_range)) %>%
    mutate(mean_count = round(mean_count, 0)) %>%
    group_by(junction_range) %>%
    summarise(count_range = paste(range(mean_count, na.rm = TRUE),collapse = "-"),
              .groups = "drop") %>%
    mutate("non_overlapping_junctions_selected_event" = paste0(junction_range," (",count_range,")"))
  
  region_junction_df <- region_junction_df[,c("event_id", "cell", "treatment", "junction_range", "mean_count")]
  region_junction_df <- left_join(region_junction_df, NMD_table_for_app[,c("event_id", "perc_region_avg", "NMD_sister_transcript", "NMD_exclusive")])
  
  colnames(region_junction_df) <- c("event_id", "cell", "treatment", "interfering_junction_range", "interfering_junction_count", "perc_region_avg", "NMD_sister_transcript", "NMD_exclusive")
  
  region_junction_df$overlapping_junction_range <- region_junction_df$interfering_junction_range %in% junction_df$junction_range
  
  region_junction_df <- region_junction_df %>% 
    mutate(interfering_junction_count = round(interfering_junction_count, 0)) %>%
    group_by(event_id,interfering_junction_range,overlapping_junction_range, perc_region_avg, NMD_sister_transcript, NMD_exclusive) %>%
    summarise(count_range = paste(range(interfering_junction_count, na.rm = TRUE),collapse = "-"),
      .groups = "drop") %>%
    group_by(event_id, perc_region_avg, NMD_sister_transcript, NMD_exclusive) %>%
    summarise(
      ##create event specfic and overlapping junctions
      event_specific_junctions = {
        keep <- !overlapping_junction_range
        if (any(keep)) {
          paste0(interfering_junction_range[keep]," (",count_range[keep],")",collapse = "; ")
        } else {
          NA_character_
        }
      },
      
      overlapping_junctions = {keep <- overlapping_junction_range
        if (any(keep)) {
          paste0(interfering_junction_range[keep]," (",count_range[keep],")",collapse = "; ")
        } else {
          NA_character_
        }
      },
      ##create selected event unique junctions
      non_overlapping_junctions_selected_event = paste0(
        unique(junction_df_sum$non_overlapping_junctions_selected_event[!junction_df_sum$junction_range %in% interfering_junction_range]), collapse = "; "),
      
      .groups = "drop"
    )
  
  region_junction_df$show_interfering_event <- vapply(
    region_junction_df$event_id,
    function(event_id) {
      as.character(
        checkboxInput(
          inputId = paste0("checkb_related_junctions_", event_id),
          label = NULL,
          value = FALSE
        )
      )
    },
    character(1)
  )
  
  #make columns pretty 
  region_junction_df$perc_region_avg <- round(region_junction_df$perc_region_avg*100)
  
  region_junction_df$perc_region_avg <- ifelse(
    region_junction_df$perc_region_avg < 20,
    paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:red;" title="Low contribution of LR event coverage to LR region coverage (<20%): Indicative for an irrelavent AS-event annotation and inacurate PSI estimation"></span> (', region_junction_df$perc_region_avg, "%)"),
    ifelse(
      region_junction_df$perc_region_avg < 60,
      paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:orange;" title="Low contribution of LR event coverage to LR region coverage (20–49%): Indicative for complex splicing or a less relavent AS-event annotation and less acurate PSI estimation"></span> (', region_junction_df$perc_region_avg, "%)"),
      ifelse(
        region_junction_df$perc_region_avg < 80,
        paste0('<span class="glyphicon glyphicon-exclamation-sign" style="color:gold;" title="Moderate contribution of LR event coverage to LR region coverage (50–89%): Indicative for complex splicing or a less relavent AS-event annotation and less acurate PSI estimation"></span> (', region_junction_df$perc_region_avg, "%)"),
        paste0('<span class="glyphicon glyphicon-ok" style="color:green;" title="High contribution of LR event coverage to LR region coverage (≥90%): Indicative for a relavent AS-event annotation and acurate PSI estimation"></span> (', region_junction_df$perc_region_avg, "%)")
      )
    )
  )
  
  return(region_junction_df)
}

shiny_click_to_coordinates <- function(click, x_range, y_range) {
  
  domain_left <- click$domain$left
  domain_right <- click$domain$right
  domain_bottom <- click$domain$bottom
  domain_top <- click$domain$top
  
  x_coord <- x_range[1] +
    ((click$x - domain_left) / (domain_right - domain_left)) *
    diff(x_range)
  
  # Shiny/plot coordinates have y increasing downwards,
  # whereas ggplot data coordinates increase upwards.
  y_coord <- y_range[1] +
    ((click$y - domain_bottom) / (domain_top - domain_bottom)) *
    diff(y_range)
  
  return(list(
    x = x_coord,
    y = y_coord
  ))
}

return_clicked_junction <- function(click, junction_df, region_junction_df, x_range, y_range){
  clicked <- shiny_click_to_coordinates(
    click = click,
    x_range = x_range,
    y_range = y_range
  )
  
  x_coord <- clicked$x
  y_coord <- clicked$y
  
  region_junction_df <- region_junction_df[region_junction_df$show_interfering_event,]
  
  junction_df$related_junction <- "in_event"
  region_junction_df$related_junction <- "in_related_event"
  
  combined_junction_df <- rbind(junction_df, region_junction_df %>% select(colnames(junction_df)))
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
  #url <- paste("https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=3191311462&position=chr", chrom, ":", region, sep = "")
  session_user <- "Kim"
  session_name <- "SpliNEx_LR_data"
  hgsid <- "3339419303"
  
  # Build UCSC Genome Browser URL
  url <- paste0(
    "https://genome.ucsc.edu/cgi-bin/hgTracks?",
    "hgsid=", hgsid,
    "&hgS_doOtherUser=submit",
    "&hgS_otherUserName=", session_user,
    "&hgS_otherUserSessionName=", session_name,
    "&position=chr", chrom, ":", region
  )
  
  # Print or use the URL
  url
  
  return(url)
}
