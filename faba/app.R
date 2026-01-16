library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)
library(plotly)
library(Seurat)
library(shinyhelper)
library(data.table)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(gridExtra)
library(fastmatch)
library(DBI)
library(RSQLite)
library(data.table)
library(RColorBrewer)             # load every session
library(cowplot)             # load each session



#set work DIRECTORY
setwd('~/Downloads/faba_spatial_seedD25_shiny/')

# load data
load_app_data <- function() {
  tryCatch({
    data <- list()
    
    data_files <- c("metadata.rds", "reductions.rds", "gene_stats.rds", 
                    "cluster_stats.rds", "expression_matrix_sparse.rds", "basic_info.rds")
    
    missing_files <- data_files[!file.exists(data_files)]
    if(length(missing_files) > 0) {
      stop(paste("lack data:", paste(missing_files, collapse = ", ")))
    }
    
    # load data
    data$metadata <- readRDS("metadata.rds")
    data$reductions <- readRDS("reductions.rds") 
    data$gene_stats <- readRDS("gene_stats.rds")
    data$cluster_stats <- readRDS("cluster_stats.rds")
    data$expression_sparse <- readRDS("expression_matrix_sparse.rds")
    data$basic_info <- readRDS("basic_info.rds")
    
    return(data)
  }, error = function(e) {
    stop(paste("fail to load:", e$message))
  })
}

app_data <- load_app_data()

# get gene expression
get_gene_expression <- function(gene) {
  gene_data <- app_data$expression_sparse[app_data$expression_sparse$gene == gene, ]
  if(nrow(gene_data) == 0) return(NULL)
  
  all_cells <- rownames(app_data$metadata)
  expr_vector <- rep(0, length(all_cells))
  names(expr_vector) <- all_cells
  
  if(nrow(gene_data) > 0) {
    valid_cells <- intersect(gene_data$cell_id, all_cells)
    expr_vector[valid_cells] <- gene_data$expression[gene_data$cell_id %in% valid_cells]
  }
  
  return(expr_vector)
}
anno<-read.table('HedinV2.mercator.functional.txt',sep = '\t',header = T)
cogene <- readRDS("cogene.rds")

ui <- dashboardPage(
  dashboardHeader(title = "D25 seed- Integrated Spatial Transcriptome & Co-expression Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Co-expression Analysis", tabName = "coexp", icon = icon("project-diagram")),
      menuItem("Spatial transcriptome Analysis", tabName = "sc", icon = icon("th"),
               menuSubItem("Cluster Analysis", tabName = "clusters", icon = icon("layer-group")),
               menuSubItem("Gene Expression", tabName = "genes", icon = icon("dna")),
               menuSubItem("Data Info", tabName = "info", icon = icon("info-circle"))
      )
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #f4f4f4; }
        .plot-container { margin: 10px 0; }
      "))
    ),
    
    tabItems(
      # —— coexpression tab ——
      tabItem(tabName = "coexp",
              fluidRow(
                box(
                  width = 12,
                  title = "Gene Co-expression Analysis",
                  status = "primary",
                  solidHeader = TRUE,
                  "Select Target Gene for Co-Expression Analysis Using Pearson Correlation (Based on D25 spatial transcriptome Metacell Dataset)"
                )
              ),
              fluidRow(
                box(width = 4,
                    selectizeInput("targetGene", "Select Target Gene:",
                                   choices = cogene,
                                   selected = "Vfaba.Hedin2.R2.1g010868",
                                   options = list(placeholder = "Type to search genes")
                    ),
                    sliderInput("corThreshold", "Correlation Threshold:",
                                min = -1, max = 1, value = 0.1, step = 0.01
                    ),
                    sliderInput("fdrThreshold", "FDR Threshold:",
                                min = 0, max = 1, value = 0.1, step = 0.001
                    ),
                    numericInput("topGenes", "Top Genes to Display:", value = 50, min = 5, max = 1000),
                    selectInput("plotType", "Plot Type:", choices = c("Bar Plot","Volcano Plot"))
                ),
                box(width = 8,
                    plotlyOutput("correlationPlot", height = "400px")
                )
              ),
              fluidRow(
                box(width = 12,
                    DTOutput("correlationTable")
                )
              )
      ),
      
      # —— cluster tab ——
      tabItem(tabName = "clusters",
              fluidRow(
                box(
                  title = "Cluster Selection", status = "primary", solidHeader = TRUE,
                  width = 3, height = 250,
                  selectInput("clusters", "Select Clusters:",
                              choices = sort(app_data$basic_info$clusters),
                              multiple = TRUE,
                              selected = app_data$basic_info$clusters[1]
                  ),
                  selectInput("reduction_clusters", "Reduction Method:",
                              choices = names(app_data$reductions),
                              selected = names(app_data$reductions)[1]
                  ),
                  actionButton("plot_clusters", "Generate Plot", class = "btn-primary", style = "width: 100%;"),
                  br(), br(),
                  downloadButton("download_cluster_plot", "Download Plot", class = "btn-success", style = "width: 100%;")
                ),
                box(
                  title = "Cluster Summary", status = "info", solidHeader = TRUE,
                  width = 9, height = 250,
                  DT::dataTableOutput("cluster_summary")
                )
              ),
              fluidRow(
                box(
                  title = "Cluster Visualization", status = "success", solidHeader = TRUE,
                  width = 12,
                  div(class = "plot-container",
                      plotOutput("cluster_plot", height = "600px")
                  )
                )
              )
      ),
      
      # —— gene expression tab ——
      tabItem(tabName = "genes",
              fluidRow(
                box(
                  title = "Gene Selection", status = "primary", solidHeader = TRUE,
                  width = 3, height = 250,
                  selectInput("gene", "Select Gene:",
                              choices = sort(unique(app_data$expression_sparse$gene)),
                              selected = app_data$expression_sparse$gene[1]
                  ),
                  selectInput("reduction_genes", "Reduction Method:",
                              choices = names(app_data$reductions),
                              selected = names(app_data$reductions)[1]
                  ),
                  actionButton("plot_gene", "Generate Plot", class = "btn-primary", style = "width: 100%;"),
                  br(), br(),
                  downloadButton("download_gene_plot", "Download Plot", class = "btn-success", style = "width: 100%;")
                ),
                box(
                  title = "Gene Information", status = "info", solidHeader = TRUE,
                  width = 9, height = 250,
                  verbatimTextOutput("gene_info")
                )
              ),
              fluidRow(
                box(
                  title = "Gene Expression Visualization", status = "success", solidHeader = TRUE,
                  width = 12,
                  div(class = "plot-container",
                      plotOutput("gene_plot", height = "600px")
                  )
                )
              )
      ),
      
      # —— data info tab ——
      tabItem(tabName = "info",
              fluidRow(
                box(
                  title = "Dataset Overview", status = "primary", solidHeader = TRUE,
                  width = 6,
                  verbatimTextOutput("data_overview")
                ),
                box(
                  title = "Processing Info", status = "info", solidHeader = TRUE,
                  width = 6,
                  verbatimTextOutput("processing_info")
                )
              ),
              fluidRow(
                box(
                  title = "Available Reductions", status = "success", solidHeader = TRUE,
                  width = 6,
                  DT::dataTableOutput("reductions_table")
                ),
                box(
                  title = "Top Expressed Genes", status = "warning", solidHeader = TRUE,
                  width = 6,
                  DT::dataTableOutput("top_genes_table")
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  
  ## --------------------- ##
  ## co-expression analysis part
  ## --------------------- ##
  correlation_data <- reactive({
    req(input$targetGene)
    conn <- dbConnect(SQLite(), "correlations.db")
    query <- paste0("SELECT * FROM correlations WHERE gene1 = '", input$targetGene, "'")
    result <- dbGetQuery(conn, query)
    dbDisconnect(conn)
    match_idx<-fmatch(result[,2],gsub("\\.?-?[0-9]+$", "", anno$IDENTIFIER))
    result <- cbind(result, anno[match_idx, 4:6])
    rownames(result) <- result$gene2
    result <- result[, -c(1,2)]
    # return(cor_df)
    return(result)
  })
  filtered_data <- reactive({
    correlation_data() %>%
      filter(abs(correlation) >= input$corThreshold,
             fdr <= input$fdrThreshold) %>%
      arrange(desc(abs(correlation))) %>%
      head(input$topGenes)
  })
  
  output$correlationPlot <- renderPlotly({
    data <- filtered_data()
    plot_title <- paste("Gene Co-expression with", input$targetGene)
    if (input$plotType == "Bar Plot") {
      p <- ggplot(data, aes(x = reorder(rownames(data), correlation), y = correlation)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = plot_title)
      ggplotly(p)
    } else {
      p <- ggplot(correlation_data(), aes(x = correlation, y = -log10(fdr))) +
        geom_point(aes(color = correlation)) +
        theme_minimal() +
        labs(title = plot_title)
      ggplotly(p)
    }
  })
  
  output$correlationTable <- renderDT({
    datatable(filtered_data(), options = list(pageLength = 20))
  })
  
  plot_clusters_web <- function(clusters, reduction = "umap") {
    meta_df <- app_data$metadata
    red_coords <- app_data$reductions[[reduction]]
    
    plot_data <- cbind(meta_df, red_coords)
    
    p1 <- ggplot(plot_data, aes(x = imagecol, y = imagerow, color = seurat_clusters)) +
      geom_point(size = 0.5, alpha = 0.7) +
      coord_fixed() +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = "Spatial Distribution (All Clusters)",
           x = "X coordinate", y = "Y coordinate")
    
    p2 <- ggplot(plot_data, aes(x = imagecol, y = imagerow)) +
      geom_point(color = "lightgrey", size = 0.4, alpha = 0.5) +
      geom_point(data = subset(plot_data, seurat_clusters %in% clusters),
                 aes(color = seurat_clusters), size = 0.6, alpha = 0.8) +
      scale_color_manual(values = setNames(brewer.pal(max(3, length(clusters)), "Set1"), clusters)) +
      coord_fixed() +
      theme_minimal() +
      labs(title = paste0("Highlighted Clusters: ", paste(clusters, collapse = ", ")),
           x = "X coordinate", y = "Y coordinate")
    
    coord_names <- colnames(red_coords)
    p3 <- ggplot(plot_data, aes_string(x = coord_names[1], y = coord_names[2], color = "seurat_clusters")) +
      geom_point(size = 0.5, alpha = 0.7) +
      theme_minimal() +
      theme(legend.position = "right") +
      labs(title = paste("DimPlot (All Clusters,", toupper(reduction), ")"),
           x = coord_names[1], y = coord_names[2])+
      guides(color = guide_legend(override.aes = list(size = 3)))
    
    
    highlight_cells <- rownames(plot_data)[plot_data$seurat_clusters %in% clusters]
    plot_data$highlight <- ifelse(rownames(plot_data) %in% highlight_cells, "Highlighted", "Other")
    
    p4 <- ggplot(plot_data, aes_string(x = coord_names[1], y = coord_names[2], color = "highlight")) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(values = c("Other" = "lightgrey", "Highlighted" = "red")) +
      theme_minimal() +
      labs(title = paste("Highlighted:", paste(clusters, collapse = ", ")),
           x = coord_names[1], y = coord_names[2], color = "")
    
    # 2x2
    plot_grid(p1, p2, p3, p4, ncol = 2)
  }
  
  # gene expression visulization
  plot_gene_web <- function(gene, reduction = "umap") {
    meta_df <- app_data$metadata
    red_coords <- app_data$reductions[[reduction]]
    
    # get gene expression
    expr <- get_gene_expression(gene)
    if(is.null(expr)) {
      return(ggplot() + ggtitle("Gene not found"))
    }
    
    # 
    plot_data <- cbind(meta_df, red_coords)
    plot_data$expression <- expr[rownames(plot_data)]
    plot_data$expression[is.na(plot_data$expression)] <- 0
    
    # Plot 1: spatial distribution (all clusters)
    p1 <- ggplot(plot_data, aes(x = imagecol, y = imagerow, color = seurat_clusters)) +
      geom_point(size = 0.5, alpha = 0.7) +
      coord_fixed() +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = "Spatial Distribution (All Clusters)",
           x = "X coordinate", y = "Y coordinate")
    
    # Plot 2: spatial gene expression
    p2 <- ggplot(plot_data, aes(x = imagecol, y = imagerow, color = expression)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_gradient(low = "lightgrey", high = "red", name = "Expression") +
      coord_fixed() +
      theme_minimal() +
      labs(title = paste0("Spatial Expression (", gene, ")"),
           x = "X coordinate", y = "Y coordinate")
    
    # Plot 3: umap (all clusters)
    coord_names <- colnames(red_coords)
    p3 <- ggplot(plot_data, aes_string(x = coord_names[1], y = coord_names[2], color = "seurat_clusters")) +
      geom_point(size = 0.5, alpha = 0.7) +
      theme_minimal() +
      theme(legend.position = "right") +
      labs(title = paste("DimPlot (", toupper(reduction), ")"),
           x = coord_names[1], y = coord_names[2])+
      guides(color = guide_legend(override.aes = list(size = 3)))
    
    
    # Plot 4: gene expression in DimPlot
    p4 <- ggplot(plot_data, aes_string(x = coord_names[1], y = coord_names[2], color = "expression")) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_gradient(low = "lightgrey", high = "red", name = "Expression") +
      theme_minimal() +
      labs(title = paste("Expression (", gene, ")"),
           x = coord_names[1], y = coord_names[2])
    
    # 2x2
    plot_grid(p1, p2, p3, p4, ncol = 2)
  }
  
  # 
  cluster_plot_reactive <- reactive({
    req(input$clusters, input$reduction_clusters)
    plot_clusters_web(input$clusters, input$reduction_clusters)
  })
  
  gene_plot_reactive <- reactive({
    req(input$gene, input$reduction_genes)
    plot_gene_web(input$gene, input$reduction_genes)
  })
  
  # generate cluster plot
  observeEvent(input$plot_clusters, {
    output$cluster_plot <- renderPlot({
      cluster_plot_reactive()
    })
  })
  
  # generate plot
  observeEvent(input$plot_gene, {
    output$gene_plot <- renderPlot({
      gene_plot_reactive()
    })
  })
  
  # output
  output$cluster_summary <- DT::renderDataTable({
    DT::datatable(app_data$cluster_stats, options = list(pageLength = 10))
  })
  
  output$gene_info <- renderText({
    req(input$gene)
    expr <- get_gene_expression(input$gene)
    if(is.null(expr)) return("Gene not found")
    
    paste("Gene:", input$gene, "\n",
          "Mean expression:", round(mean(expr), 4), "\n", 
          "Max expression:", round(max(expr), 4), "\n",
          "Expressing cells:", sum(expr > 0), "/", length(expr),
          "(", round(sum(expr > 0)/length(expr) * 100, 2), "%)")
  })
  
  output$data_overview <- renderText({
    paste("Cells:", app_data$basic_info$n_cells, "\n",
          "Original Features:", app_data$basic_info$n_features, "\n",
          "Available Genes:", length(unique(app_data$expression_sparse$gene)), "\n",
          "Clusters:", app_data$basic_info$n_clusters, "\n",
          "Reductions:", length(app_data$reductions))
  })
  
  output$processing_info <- renderText({
    "This is a lightweight web version using pre-processed data.\n\nData has been filtered to include only expressed genes to reduce size for web deployment.\n\nAll visualizations are generated in real-time from the processed data."
  })
  
  output$reductions_table <- DT::renderDataTable({
    red_info <- data.frame(
      Method = names(app_data$reductions),
      Dimensions = sapply(app_data$reductions, ncol)
    )
    DT::datatable(red_info, options = list(pageLength = 5))
  })
  
  output$top_genes_table <- DT::renderDataTable({
    top_genes <- app_data$gene_stats %>%
      arrange(desc(mean_expr)) %>%
      select(gene, mean_expr, max_expr, pct_cells) %>%
      head(20)
    DT::datatable(top_genes, options = list(pageLength = 10)) %>%
      DT::formatRound(columns = c('mean_expr', 'max_expr', 'pct_cells'), digits = 3)
  })
  
  # download handlers
  output$download_cluster_plot <- downloadHandler(
    filename = function() {
      paste0("cluster_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, cluster_plot_reactive(), width = 12, height = 8, dpi = 300)
    }
  )
  
  output$download_gene_plot <- downloadHandler(
    filename = function() {
      paste0("gene_plot_", input$gene, "_", Sys.Date(), ".png") 
    },
    content = function(file) {
      ggsave(file, gene_plot_reactive(), width = 12, height = 8, dpi = 300)
    }
  )
}

# run app
shinyApp(ui = ui, server = server)
