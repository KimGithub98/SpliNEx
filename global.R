#script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(script_dir)
#setwd("C:/Users/Z978196/OneDrive - Radboudumc/z978196/R/Shiny Apps/NMD_event lookup structured")
library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(shinyjs)
library(RColorBrewer)
library(bslib)
library(svglite)
source("R/Event_details_functions.R")

#-------------Load data ------------------------------
normalized_counts_df <- readRDS("data/normalized_counts_df.rds")
rownames(normalized_counts_df) <- toupper(rownames(normalized_counts_df))
colData_df <- readRDS("data/colData_df.rds")
colData_df$Data_source_pretty_name[colData_df$Data_source_pretty_name == " Emma et al. 2025"] <- " Unpublished"
NMD_table_for_app <- readRDS("data/NMD_table_for_app.rds")
junction_summary_for_app <- readRDS("data/junction_summary_for_app.rds")

#-------------Global variables -----------------------
sample_names <- names(normalized_counts_df)
batch_ids <- gsub("_[0-9]+(?=(_|$))", "", sample_names, perl = TRUE)
batch_list <- split(sample_names, batch_ids)

NMD_event_types <- unique(NMD_table_for_app$NMD_sister_transcript)
NMD_only_types <- unique(NMD_table_for_app$NMD_exclusive)
event_types <- unique(NMD_table_for_app$event_type)

#-------------costum palette -----------------------
div_colors <- c(
  "iPSCs"              = "#cfcece",
  "Neurons (DIV21)"    = "#cbf4ce",
  "Neurons (DIV28)"    = "#95ba9a",
  "Neurons (DIV49)"    = "#637b64",
  "Astrocytes (DIVx)"  = "#b29ac9",
  "Microglia (DIVx)"   = "#9ebae1"
)
