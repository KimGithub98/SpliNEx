script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(shinyjs)
library(RColorBrewer)
library(bslib)
source(paste(script_dir,"/R/Event_details_functions.R", sep = ""))

#-------------Load data ------------------------------
normalized_counts_df <- read.csv(file.path(script_dir, "/data/SR_sofia_Kim_WTC_normalized_reads.csv"),row.names = 1)
colData_df <- read.csv(file.path(script_dir, "/data/SR_sofia_Kim_WTC_normalized_reads_col_data.csv"),row.names = 1)
NMD_table_for_app <- read.csv(file.path(script_dir, "/data/NMD_table_for_app_Astro_micro.csv"))
junction_summary_for_app <- read.csv(file.path(script_dir, "/data/junction_table_for_app_Astro_micro.csv"))

#-------------Global variables -----------------------
sample_names <- names(normalized_counts_df)
batch_ids <- gsub("_[0-9]+(?=(_|$))", "", sample_names, perl = TRUE)
batch_list <- split(sample_names, batch_ids)

NMD_event_types <- unique(NMD_table_for_app$NMD_sister_transcript)
NMD_only_types <- unique(NMD_table_for_app$NMD_exclusive)
event_types <- unique(NMD_table_for_app$event_type)

#-------------costum palette -----------------------
div_colors <- c(
  "iPSC"       = "ivory2",
  "DIV21"      = "cadetblue1",
  "DIV28"      = "cadetblue2",
  "DIV49"      = "cadetblue",
  "Astrocytes" = "palevioletred",
  "Microglia"  = "burlywood4"
)
