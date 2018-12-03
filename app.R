library(shiny)
library(tibble)
library(readr)
library(dplyr)
library(here)
library(tidyr)
library(magrittr)
library(purrr)
source("./AUCFunction.R")
source("./string_processing.R")

#load reference data on session start
saunders_ranks_matrix_types <- readRDS(here('data', 'processed', 'saunders_ranks_matrix.RDS'))
saunders_ranks_matrix_classes <- readRDS(here('data', 'processed', 'saunders_classes_ranks_matrix.RDS'))
  
unique_genes <- saunders_ranks_matrix_types$gene_symbol

apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  # App title ----
  titlePanel("Polygenic tester for mouse nervous system (Saunders, et al., dataset)"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      radioButtons(inputId = "data_selection",
                   label = "Restrict analysis to:",
                   choices = c("Cell classes (n=13)" = "classes",
                     "All tissue subclusters (n=565)" = "subclusters")
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
        verbatimTextOutput("summary"),
        br(),
        dataTableOutput("view")
        #plotlyOutput("dotplot"),
        #verbatimTextOutput("info")
      )
    )
  )
)


# Define server logic process and output top cortical layers/zones ----
server <- function(input, output) {
  observeEvent(input$submit, {
    start <- Sys.time()

    if (input$data_selection == "classes") {
      saunders_ranks_matrix <- saunders_ranks_matrix_classes
    } else {
      saunders_ranks_matrix <- saunders_ranks_matrix_types
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
    if (input$data_selection == "subclusters") {
      meta <- readRDS(here('data' , 'processed', 'meta.RDS'))
      meta %<>% select(tissue_subcluster, full_name, tissue_name, class)
      table <- left_join(table, meta, by="tissue_subcluster")
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
    table
  }, escape = FALSE)
  })
}

shinyApp(ui, server)



