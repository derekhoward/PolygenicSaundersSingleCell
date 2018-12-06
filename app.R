library(shiny)
library(tibble)
library(dplyr)
library(here)
library(tidyr)
library(shinyjs)
library(magrittr)
library(purrr)
source("./AUCFunction.R")
source("./string_processing.R")

#load reference data on session start
saunders_ranks_matrix_types <- readRDS(here('data', 'processed', 'saunders_ranks_matrix.RDS'))
saunders_ranks_matrix_classes <- readRDS(here('data', 'processed', 'saunders_classes_ranks_matrix.RDS'))
meta <- readRDS(here('data' , 'processed', 'meta.RDS'))

  
unique_genes <- saunders_ranks_matrix_types$gene_symbol

apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(includeHTML("google-analytics.html")),
  # App title ----
  titlePanel("Polygenic tester for mouse nervous system (Saunders, et al., dataset)"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      radioButtons(inputId = "data_selection",
                   label = "Restrict analysis to:",
                   choices = c("All regions - cell classes (n=13)" = "classes",
                     "All regions - tissue subclusters (n=565 - slow)" = "subclusters", sort(unique(meta$tissue_name)))
                   ),
      
      textAreaInput(inputId = "genelist",
                    label = "Input your gene list:",
                    #value = 'Mag\nMobp\nMog\nMbp\nOmg',
                    value = 'Opalin\nMag\nMobp\nCarns1\nKlk6\nErmn\nTmem144',
                    rows = 7),
                  
      selectInput(inputId = 'species',
                  label = 'Species of input genes:',
                  choices = c('Mouse', 'Human')),
                
      actionButton(inputId = "submit",
                   label = "Submit")
),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(
        id = "main",
        
        p("This tool is made possible by data from:"),
        tags$p("Saunders A*, Macosko E.Z*, Wysoker A, Goldman M, Krienen, F, de Rivera H, Bien E, Baum M, Wang S, Bortolin L, Goeva A, Nemesh J, Kamitaki N, Brumbaugh S, Kulp D and McCarroll, S.A. 2018. Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain. "),
        tags$a(href="https://doi.org/10.1016/j.cell.2018.07.028", "2018. Cell. 174(4) P1015-1030.E16"),
        br(),
        br()),
        
        verbatimTextOutput("summary"),
        br(),
        dataTableOutput("view")
    )
  )
)


# Define server logic process and output top cortical layers/zones ----
server <- function(input, output) {
  observeEvent(input$submit, {
    start <- Sys.time()
    shinyjs::hide("main")
    shinyjs::disable("submit") 
    
    if (input$data_selection == "classes") {
      saunders_ranks_matrix <- saunders_ranks_matrix_classes
    } else if (input$data_selection == "subclusters") {
      saunders_ranks_matrix <- saunders_ranks_matrix_types
    } else {
      saunders_ranks_matrix <- saunders_ranks_matrix_types
      tissue_subcluster_to_keep <- meta %>% filter(tissue_name == input$data_selection) %>% .$tissue_subcluster
      #tissue_subcluster_to_keep <- meta %>% filter(tissue_name == "Thalamus") %>% .$tissue_subcluster
      saunders_ranks_matrix %<>% select_(.dots=c("gene_symbol", quote(tissue_subcluster_to_keep)))
    }
    
    unique_genes <- saunders_ranks_matrix$gene_symbol

    cleaned_gene_list <- isolate(process_input_genes(input$genelist))
    if (input$species == 'Human') {
      cleaned_gene_list <- convert_genes(cleaned_gene_list)
    }
    
    #for indices - use dplyr for ease
    forIndices <- as_tibble(saunders_ranks_matrix$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    # only columns from cortical zones remain in df
    df <- saunders_ranks_matrix %>%
      select(-gene_symbol)
    #print(cleaned_gene_list)
    #print(head(forIndices))
    #AUROC <- map_df(df, auroc_analytic, as.numeric(targetIndices))
    AUROC <- map_df(df, auroc_analytic, targetIndices)
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <- bind_cols(gather(AUROC, key = tissue_subcluster, value = AUROC), 
                       gather(wilcox_tests, value = pValue)) %>%
      select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    # these are the values for the results table
    table %<>% arrange(-AUROC)
    table %<>% mutate(pValue = signif(pValue, digits = 3), 
                      AUROC = signif(AUROC, digits = 3),
                      adjusted_P = signif(p.adjust(pValue), digits = 3))
    if (input$data_selection != "classes") {
      table <- left_join(table, meta %>% select(tissue_subcluster, full_name, tissue_name, class), by="tissue_subcluster")
    } else {
      table %<>% rename(class = tissue_subcluster)
    }
  
  output$summary <- renderPrint({
    #count of intersection of submitted genes with total gene list
    cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
    cat(paste(
      "\nGenes found in data:",
      sum(cleaned_gene_list %in% unique_genes),
      "of",
      length(cleaned_gene_list)
    ))
  })
  
  
  output$view <- renderDataTable({
    if (isolate(input$data_selection) != "classes") {
      table %<>% select(full_name, class, everything())
    }
    table
  }, escape = FALSE)
  shinyjs::enable("submit")
  
  })
  
}

shinyApp(ui, server)



