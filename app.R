
library(shiny)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tools)
library(DT)
library(svglite)
library(readxl)      # for .xlsx/.xls
library(SummarizedExperiment)
# library(DESeq2)

# if (!exists("theme_publish")) theme_publish <- function(...) ggplot2::theme_minimal()
if (!exists("theme_publish")) {
  theme_publish <- function() {
    ggplot2::theme_minimal(base_size = 14, base_family = "sans")
  }
}

ui <- fluidPage(
  titlePanel("Volcano interactive"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "file", "Upload dataset (CSV/TSV/TXT or RDS)",
        accept = c(".csv",".CSV",".tsv",".TSV",".txt",".TXT",".rds",".RDS")
      ),
      # radioButtons(
      #   "filetype", "File type", inline = TRUE,
      #   choices = c("Auto"="auto","CSV"="csv","RDS"="rds"), selected="auto"
      # ),
      helpText("Select columns after upload."),
      
      uiOutput("col_selectors"),
      
      # selectizeInput(
      #   "genes", "Genes to show / label",
      #   choices = NULL, multiple = TRUE,
      #   options = list(placeholder = "Type or pick symbols")
      # ),
      textAreaInput(
        "gene_text", "Paste gene symbols", rows = 3,
        placeholder = "TP53, EGFR, STAT1"
      ),
      
      fluidRow(
        column(4, actionButton("apply_genes", "Use pasted")),
        column(4, actionButton("clear_genes", "Clear")),
        column(4, actionButton("list_sig", "List significant"))
      ),
      
      fluidRow(
        column(6, actionButton("select_sig", "Use all significant")),
        column(6, checkboxInput("only_selected", "Only show selected genes", FALSE))
      ),
      
      checkboxInput("auto_add_brush", "Auto-add brushed genes", TRUE),
      checkboxInput("show_sel_circle", "Highlight selected with circle", value = TRUE),
      
      tags$hr(),
      sliderInput("point_size", "Point size", min = 0.5, max = 5, value = 2.5, step = 0.5),
      sliderInput("label_size", "Label text size", min = 2, max = 8, value = 3, step = 0.5),
      numericInput("y_limit",  "Y cap",          value = 40, min = 5, max = 200, step = 5),
      numericInput("FC_cutoff","FC cutoff",      value = 1,  step = 0.1),
      numericInput("padj_thr", "padj threshold", value = 0.05, step = 0.01),
      textInput("plot_title", "Plot title", value = "Tumor vs Control"),
      textInput("plot_subtitle", "Plot subtitle", value = "NK"),
      numericInput("plot_width",  "Plot width",  value = 8, min = 4, max = 20, step = 0.5),
      numericInput("plot_height", "Plot height", value = 7, min = 4, max = 20, step = 0.5),
      
      tags$hr(),
      fluidRow(
        column(
          6,
          selectInput(
            "download_format", "Download format",
            choices = c("PNG"="png","PDF"="pdf","SVG"="svg"),
            selected = "png"
          )
        ),
        column(
          6,
          downloadButton("download_plot", "Download plot")
        )
      )
    ),
    mainPanel(
      h4("Preview of uploaded data"),
      DTOutput("preview"),
      tags$hr(),
      plotOutput(
        "volcano", height = "720px",
        brush = brushOpts(id = "volcano_brush", resetOnNew = TRUE),
        click = "volcano_click"
      ),
      tags$br(),
      strong("Selected genes:"),
      verbatimTextOutput("selected_preview", placeholder = TRUE)
    )
  )
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

server <- function(input, output, session) {
  raw_ds <- reactive({
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    type <- ext
    
    x <- switch(
      type,
      
      "csv" = {
        hdr <- readLines(input$file$datapath, n = 1)
        delim <- if (grepl("\t", hdr)) "\t" else ","
        read.delim(input$file$datapath, sep = delim, check.names = FALSE)
      },
      
      "rds" = {
        obj <- readRDS(input$file$datapath)
        
        if (is.data.frame(obj)) {
          obj
        } else if (is.matrix(obj)) {
          as.data.frame(obj)
        } else {
          validate(need(
            FALSE,
            paste(
              "This RDS is not a data.frame or matrix.",
              "Please convert it in R before uploading:",
              "",
              "Example:",
              "res_df <- as.data.frame(your_DESeqResults_object)",
              "saveRDS(res_df, 'res_df.rds')"
            )
          ))
        }
      },
      
      
      "xlsx" = {
        readxl::read_excel(input$file$datapath)
      },
      
      "xls" = {
        readxl::read_excel(input$file$datapath)
      },
      
      {
        hdr <- readLines(input$file$datapath, n = 1)
        delim <- if (grepl("\t", hdr)) "\t" else ","
        read.delim(input$file$datapath, sep = delim, check.names = FALSE)
      }
    )
    
    if (is.matrix(x)) x <- as.data.frame(x)
    validate(need(is.data.frame(x), "Uploaded file is not a tabular object."))
    
    if (!is.null(rownames(x)) && any(nzchar(rownames(x))) && !".rownames" %in% names(x)) {
      x$.rownames <- rownames(x)
    }
    
    x
  })
  

  observeEvent(input$file, {
    req(input$file$name)
    fname <- tools::file_path_sans_ext(basename(input$file$name))
    updateTextInput(session, "plot_title", value = fname)
  }, ignoreInit = TRUE)
  
  # output$preview <- renderDT({
  #   req(raw_ds())
  #   datatable(head(raw_ds(), 50), options = list(scrollX = TRUE))
  # })
  # 
  output$preview <- renderDT({
    req(raw_ds())
    datatable(
      head(raw_ds(), 5),
      options = list(scrollX = TRUE, pageLength = 5))
  })
  
  output$col_selectors <- renderUI({
    df <- raw_ds()
    cn <- names(df)
    
    guess_gene <- cn[grep("gene_symbol|gene.symbol|symbol|gene|gene_name", cn, ignore.case = TRUE)][1]
    guess_lfc  <- cn[grep("log2foldchange|log2_fc|log2.fc|logfc|lfc", cn, ignore.case = TRUE)][1]
    guess_padj <- cn[grep("padj|adj.p|fdr|qvalue", cn, ignore.case = TRUE)][1]
    
    tagList(
      selectInput("gene_col", "Gene symbol column", choices = cn, selected = guess_gene),
      selectInput("lfc_col",  "log2FoldChange column", choices = cn, selected = guess_lfc),
      selectInput("padj_col", "padj / FDR column", choices = cn, selected = guess_padj)
    )
  })
  
  ds <- reactive({
    df <- raw_ds()
    req(input$gene_col, input$lfc_col, input$padj_col)
    validate(
      need(input$gene_col  %in% names(df), "Invalid gene column"),
      need(input$lfc_col   %in% names(df), "Invalid log2FC column"),
      need(input$padj_col  %in% names(df), "Invalid padj column")
    )
    
    out <- tibble::tibble(
      Gene_Symbol    = as.character(df[[input$gene_col]]),
      log2FoldChange = suppressWarnings(as.numeric(df[[input$lfc_col]])),
      padj           = suppressWarnings(as.numeric(df[[input$padj_col]]))
    )
    
    validate(
      need(!all(is.na(out$log2FoldChange)), "All log2FoldChange are NA after conversion."),
      need(!all(is.na(out$padj)), "All padj are NA after conversion.")
    )
    
    out
  })
  
  # observeEvent(ds(), {
  #   updateSelectizeInput(
  #     session, "genes",
  #     choices = sort(unique(ds()$Gene_Symbol)),
  #     server = TRUE
  #   )
  # }, ignoreInit = TRUE)
  
  sel_genes <- reactiveVal(character(0))
  
  canonicalize <- function(vec, pool) {
    if (length(vec) == 0) return(character(0))
    upool <- toupper(pool)
    uvec  <- unique(toupper(vec))
    idx   <- match(uvec, upool)
    pool[stats::na.omit(idx)]
  }
  
  # observeEvent(input$genes, {
  #   pool <- req(ds())$Gene_Symbol
  #   sel_genes(canonicalize(input$genes, pool))
  # }, ignoreInit = TRUE)
  
  # observeEvent(input$apply_genes, {
  #   pool <- req(ds())$Gene_Symbol
  #   pasted <- gsub("[,;\\s]+", " ", trimws(input$gene_text))
  #   vec <- unique(unlist(strsplit(pasted, " ")))
  #   vec <- vec[nzchar(vec)]
  #   all <- unique(c(sel_genes(), canonicalize(vec, pool)))
  #   sel_genes(all)
  # })
  
  observeEvent(input$apply_genes, {
    pool <- req(ds())$Gene_Symbol
    
    txt <- input$gene_text %||% ""
    txt <- trimws(txt)
    if (txt == "") return(NULL)
    
    # split on comma, semicolon, or any whitespace (spaces, tabs, newlines)
    tokens <- unlist(strsplit(txt, "[,;[:space:]]+"))
    tokens <- trimws(tokens)
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens) == 0) return(NULL)
    
    # case-insensitive + trim match to pool
    pool_clean   <- trimws(pool)
    upool        <- toupper(pool_clean)
    utokens      <- toupper(tokens)
    idx          <- match(utokens, upool)
    matched_idx  <- idx[!is.na(idx)]
    if (length(matched_idx) == 0) return(NULL)
    
    picked <- pool[matched_idx]
    
    sel_genes(unique(c(sel_genes(), picked)))
  })
  
  
  
  observeEvent(input$clear_genes, {
    sel_genes(character(0))
    updateSelectizeInput(session, "genes", selected = character(0))
  }, ignoreInit = TRUE)
  
  observeEvent(input$select_sig, {
    x <- req(ds())
    sig <- x %>%
      filter(
        !is.na(padj),
        padj <= input$padj_thr,
        !is.na(log2FoldChange),
        abs(log2FoldChange) >= input$FC_cutoff
      )
    g <- unique(sig$Gene_Symbol)
    sel_genes(g)
    updateSelectizeInput(session, "genes", selected = g)
  }, ignoreInit = TRUE)
  
  observeEvent(input$volcano_brush, {
    req(isTRUE(input$auto_add_brush))
    x <- req(ds()) %>% mutate(logP = -log10(padj))
    brushed <- brushedPoints(
      x, input$volcano_brush,
      xvar = "log2FoldChange", yvar = "logP"
    )
    if (nrow(brushed) > 0) {
      all <- unique(c(sel_genes(), brushed$Gene_Symbol))
      sel_genes(all)
      updateSelectizeInput(session, "genes", selected = all)
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$volcano_click, {
    x <- req(ds()) %>% mutate(logP = -log10(padj))
    hit <- nearPoints(
      x, input$volcano_click,
      xvar = "log2FoldChange", yvar = "logP",
      maxpoints = 1, threshold = 10
    )
    if (nrow(hit) == 1) {
      g <- hit$Gene_Symbol[1]
      cur <- sel_genes()
      if (g %in% cur) cur <- setdiff(cur, g) else cur <- c(cur, g)
      sel_genes(cur)
      updateSelectizeInput(session, "genes", selected = cur)
    }
  }, ignoreInit = TRUE)
  
  x_limits <- reactive({
    x <- req(ds())
    r <- range(x$log2FoldChange, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(r))) r <- c(-2, 2)
    # pad <- max(0.06 * diff(r), 0.25)
    pad <- max(0.12 * diff(r), 0.5) 
    c(r[1] - pad, r[2] + pad)
  })
  
  volcano_plot <- reactive({
    x <- ds()
    y_limit   <- input$y_limit
    FC_cutoff <- input$FC_cutoff
    padj_thr  <- input$padj_thr
    selected  <- sel_genes()
    ps        <- input$point_size
    ls        <- input$label_size
    
    main_title <- (input$plot_title %||% "Tumor vs Control")
    sub_title  <- (input$plot_subtitle %||% "NK")
    
    x <- x %>%
      mutate(
        logP = -log10(padj),
        logP_capped = pmin(logP, y_limit + 2),
        is_sel = Gene_Symbol %in% selected
      )
    
    
    plot_data <- if (isTRUE(input$only_selected) && length(selected) > 0) {
      filter(x, is_sel)
    } else x
    
    main <- filter(plot_data, logP <= y_limit)
    high <- filter(plot_data, logP >  y_limit)
    
    n_sig  <- sum(x$padj <= padj_thr, na.rm = TRUE)
    n_up   <- sum(x$log2FoldChange >= FC_cutoff & x$padj <= padj_thr, na.rm = TRUE)
    n_down <- sum(x$log2FoldChange <= -FC_cutoff & x$padj <= padj_thr, na.rm = TRUE)
    
    p <- ggplot() +
      xlab("log2FoldChange") + ylab("-log10(padj)") +
      labs(
        title = main_title,
        subtitle = sub_title,
        caption = paste0(
          "Total significant = ", n_sig,
          "   Up: ", n_up,
          "   Down: ", n_down,
          "   FC = ", FC_cutoff,
          "   padj <= ", padj_thr
        )
      ) +
      geom_vline(xintercept = c(-FC_cutoff, FC_cutoff),
                 colour = "grey", linetype = "longdash") +
      geom_hline(yintercept = -log10(padj_thr),
                 colour = "grey", linetype = "longdash") +
      geom_hline(yintercept = y_limit,
                 colour = "grey", linetype = "dashed") +
      # geom_vline(xintercept = 0, colour = "grey", size = 1) +
      geom_point(
        data = main,
        aes(x = log2FoldChange, y = logP),
        size = ps, color = "#a3a3a3", alpha = 0.3
      ) +
      geom_point(
        data = filter(main, log2FoldChange >= FC_cutoff & padj < padj_thr),
        aes(x = log2FoldChange, y = logP),
        size = ps, color = "#C00000"
      ) +
      geom_point(
        data = filter(main, log2FoldChange <= -FC_cutoff & padj < padj_thr),
        aes(x = log2FoldChange, y = logP),
        size = ps, color = "#4474c4"
      ) +
      geom_point(
        data = high,
        aes(
          x = log2FoldChange,
          y = y_limit,
          color = log2FoldChange > 0
        ),
        size = ps,
        shape = 17
      ) +
      scale_color_manual(values = c("TRUE" = "#C00000", "FALSE" = "#4474c4"), guide = "none")
    
    
    
    if (length(selected) > 0) {
      if (isTRUE(input$show_sel_circle)) {
        p <- p +
          geom_point(
            data = filter(main, is_sel),
            aes(x = log2FoldChange, y = logP),
            shape = 21, fill = NA, stroke = 1.1, size = ps + 0.5
          ) +
          geom_point(
            data = filter(high, is_sel),
            aes(x = log2FoldChange, y = y_limit),
            shape = 17, size = ps + 0.5, color = "black"
          )
      }
      
      p <- p +
        ggrepel::geom_text_repel(
          data = filter(x, is_sel),
          aes(x = log2FoldChange, y = logP_capped, label = Gene_Symbol),
          colour = "black", size = ls,
          max.overlaps = getOption("ggrepel.max.overlaps", 20)
        )
    }
    
  #   p +
  #     coord_cartesian(
  #       xlim = x_limits(),
  #       ylim = c(0, y_limit + 4),
  #       expand = FALSE
  #     ) +
  #     theme_publish()
  # })
    
    p +
      # scale_x_reverse(limits = rev(x_limits())) +
      coord_cartesian(
        xlim = x_limits(),
        ylim = c(0, y_limit + 2),
        expand = FALSE
      ) +
      theme_publish()})
  
  
  
  output$volcano <- renderPlot(
    {
      volcano_plot()
    },
    width = function()  input$plot_width  * 96,
    height = function() input$plot_height * 96,
    res = 96
  )
  
  observeEvent(input$list_sig, {
    x <- req(ds())
    sig <- x %>%
      filter(padj <= input$padj_thr) %>%
      arrange(padj)
    txt <- if (nrow(sig) == 0) "None" else paste(sig$Gene_Symbol, collapse = "\n")
    showModal(modalDialog(
      title = sprintf("Significant genes (padj <= %g)", input$padj_thr),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      tagList(
        div(style="margin-bottom:8px;", sprintf("N = %d", nrow(sig))),
        div(
          style="max-height:60vh; overflow:auto; border:1px solid #ddd; padding:8px;",
          tags$pre(txt)
        )
      )
    ))
  })
  
  
  
  output$selected_preview <- renderText({
    s <- sel_genes()
    if (length(s) == 0) "None" else paste(s, collapse = ", ")
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("volcano.", input$download_format)
    },
    content = function(file) {
      p <- volcano_plot()
      ext <- input$download_format
      w <- ifelse(is.null(input$plot_width),  8, input$plot_width)
      h <- ifelse(is.null(input$plot_height), 7, input$plot_height)
      if (ext == "png") {
        ggsave(file, plot = p, device = "png", width = w, height = h, dpi = 300 )
      } else if (ext == "pdf") {
        ggsave(file, plot = p, device = "pdf", width = w, height = h)
      } else {
        ggsave(file, plot = p, device = svglite::svglite, width = w, height = h)
      }
    }
  )
}

shinyApp(ui, server)



