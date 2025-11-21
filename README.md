**SpliNEx ShinyApp**

Code used in "Mapping NMD-coupled alternative splicing in iPSC-derived brain cells: a resource for therapeutic discovery" (Wijnant et al. 2026)

This repository contains code used to generate the ShinyApp (link Shiny App) containing the data from (Wijnant et al. 2026) <link> it contains:

- ShinyApp general app files (global.R, ui.R & server.R)
- Modules in the R folder:
  - App_info_server.R & App_info_ui.R: contains server & ui information for the SpliNEx shiny app info tab
  - WTC_gene_expression_server.R & WTC_gene_expression_ui.R : contains server & ui code for the WTC gene expression tab
  - NMD_event_per_gene_server.R & NMD_event_per_gene_ui.R: contains server & ui code for the tab to search for AS events per hgnc_symbol
  - Event_details_server.R & Event_details_ui.R & Event_details_functions.R: contains server & ui & functions for the tab to view details for AS events
-  'www' folder: images used in the SpliNEx shiny app
-  'data' folder: r data files used in the SpliNEx shiny app

This code was generated in: R version 4.4.1. 

Dependencies & packages:
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
