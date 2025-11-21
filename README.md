**SpliNEx ShinyApp**

Code used in "Mapping NMD-coupled alternative splicing in iPSC-derived brain cells: a resource for therapeutic discovery" (Wijnant et al. 2026)

This repository contains all code used to generate the SpliNEx Shiny application (link to Shiny app), which provides interactive access to the data and analyses presented in Wijnant et al. (2026) (link to paper).

**Repository Structure**

- Top-level Shiny application files
  - global.R, ui.R, server.R
- R Modules (located in the R/ directory):
  - App_info_server.R & App_info_ui.R: Server and UI logic for the App Info tab.
  - WTC_gene_expression_server.R & WTC_gene_expression_ui.R: Modules for the WTC gene expression tab.
  - NMD_event_per_gene_server.R & NMD_event_per_gene_ui.R: Modules for querying NMD-coupled alternative splicing events per HGNC symbol.
  - Event_details_server.R, Event_details_ui.R, Event_details_functions.R: Server, UI, and helper functions for the Event Details tab.
- www/ directory: Contains images and other static assets used in the SpliNEx application.
- data/ directory: Contains R data files required by the application.

**Software Environment**
All code was generated and tested using:
R version 4.4.1

**Package Dependencies**
- shiny_1.9.1
- shinythemes_1.2.0
- ggplot2_4.0.0
- dplyr_1.1.4
- tidyr_1.3.1
- stringr_1.5.1
- shinyjs_2.1.0
- RColorBrewer_1.1-3
- bslib_0.8.0
- svglite_2.1.3
 


**Please cite:**

Wijnant KN et al. x. x
