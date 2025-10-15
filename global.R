#script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(script_dir)
setwd("C:/Users/Z978196/OneDrive - Radboudumc/z978196/R/Shiny Apps/NMD_event lookup structured")
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
colData_df <- readRDS("data/colData_df.rds")
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
  "iPSC"       = "ivory2",
  "DIV21"      = "cadetblue1",
  "DIV28"      = "cadetblue2",
  "DIV49"      = "cadetblue",
  "Astrocytes" = "palevioletred",
  "Microglia"  = "burlywood4"
)
