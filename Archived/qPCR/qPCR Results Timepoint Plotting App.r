# app.R ---------------------------------------------------------------------
# Interactive Shiny app: ΔΔCq pipeline & visualisation
# --------------------------------------------------------------------------
# • Upload raw qPCR well results (.csv) and sample‑info (.tsv)
# • Inspect raw tables & computed ΔΔCq tables
# • Toggle ≥3‑rep outlier trimming
# • View and download log₂FC plots as SVG
# --------------------------------------------------------------------------

# ------------------------------- Packages ---------------------------------
# To run: install.packages(c("shiny", "tidyverse", "DT", "ggplot2"))

library(shiny)
library(tidyverse)
library(DT)

# ------------------------------ Helpers -----------------------------------
sem <- function(x) stats::sd(x) / sqrt(length(x))

make_results <- function(well_df, info_df, trim_outlier = FALSE) {
  # ----- metadata ---------------------------------------------------------
  sample_meta <- info_df %>%
    filter(Info_type == "Sample") %>%
    select(Sample = Info, Condition)

  target_meta <- info_df %>%
    filter(Info_type == "Target") %>%
    select(Target = Info_name, Info_key, Info_key_1)

  reference_genes <- target_meta %>%
    filter(Info_key  == "Reference" | Info_key_1 == "Reference") %>%
    distinct(Target) %>%
    pull()

  # ----- optional outlier trimming ---------------------------------------
  if (trim_outlier) {
    well_df <- well_df %>%
      left_join(sample_meta, by = "Sample") %>%
      group_by(Sample, Target, Condition) %>%
      group_modify(~ {
        if (nrow(.x) < 3) return(.x)
        med <- stats::median(.x$Cq)
        .x %>% mutate(dist = abs(Cq - med)) %>% slice(-which.max(dist)) %>% select(-dist)
      }) %>%
      ungroup()
  }

  # ----- replicate‑means --------------------------------------------------
  sample_target_means <- well_df %>%
    group_by(Sample, Target) %>%
    summarise(Cq = base::mean(Cq), .groups = "drop")

  ref_means <- sample_target_means %>%
    filter(Target %in% reference_genes) %>%
    group_by(Sample) %>%
    summarise(ref_mean = base::mean(Cq), .groups = "drop")

  dcq_tbl <- sample_target_means %>%
    inner_join(ref_means, by = "Sample") %>%
    mutate(dCq = Cq - ref_mean) %>%
    inner_join(sample_meta, by = "Sample")

  baseline_tbl <- dcq_tbl %>%
    filter(Condition == "Day0") %>%
    group_by(Target) %>%
    summarise(baseline_dCq = base::mean(dCq), .groups = "drop")

  results <- dcq_tbl %>%
    left_join(baseline_tbl, by = "Target") %>%
    mutate(ddCq   = dCq - baseline_dCq,
           log2FC = -ddCq)

  summary_tbl <- results %>%
    group_by(Target, Condition) %>%
    summarise(mean_log2FC = base::mean(log2FC),
              sem_log2FC  = sem(log2FC), .groups = "drop")

  list(replicates = results, summary = summary_tbl)
}

plot_log2fc <- function(res) {
  ggplot() +
    geom_point(data = res$replicates,
               aes(Condition, log2FC, group = Sample),
               alpha = 0.7, position = position_dodge(width = 0.25)) +
    geom_line(data = res$summary,
              aes(Condition, mean_log2FC, group = 1), linewidth = 1.1) +
    geom_errorbar(data = res$summary,
                  aes(Condition,
                      ymin = mean_log2FC - sem_log2FC,
                      ymax = mean_log2FC + sem_log2FC),
                  width = 0.25) +
    facet_wrap(~ Target, scales = "free_y") +
    labs(title = "Log₂ fold‑change (ΔΔCq) relative to Day 0",
         y = "log₂(FC)") +
    theme_bw() +
    theme(legend.position = "none")
}

# ------------------------------- UI ---------------------------------------
ui <- fluidPage(
  titlePanel("ΔΔCq Analysis & Plotter"),
  sidebarLayout(
    sidebarPanel(
      fileInput("well_csv", "Upload WellResultOmit.csv", accept = ".csv"),
      fileInput("info_tsv", "Upload qPCR_sample_info.tsv", accept = c(".tsv", ".txt")),
      checkboxInput("trim", "Trim furthest outlier when ≥3 tech replicates", value = FALSE),
      downloadButton("dl_plot", "Download Plot (SVG)")
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
        tabPanel("Raw Well Table", DTOutput("raw_well")),
        tabPanel("Raw Info Table", DTOutput("raw_info")),
        tabPanel("ΔΔCq Results", DTOutput("res_table")),
        tabPanel("Plot", plotOutput("plot", height = "600px"))
      )
    )
  )
)

# ----------------------------- Server -------------------------------------
server <- function(input, output, session) {
  # Reactive: load uploads --------------------------------------------------
  well_raw <- reactive({
    req(input$well_csv)
    read_csv(input$well_csv$datapath, show_col_types = FALSE) %>%
      mutate(Cq = na_if(Cq, "Undetermined"), Cq = as.numeric(Cq)) %>%
      filter(!is.na(Cq))
  })

  info_raw <- reactive({
    req(input$info_tsv)
    read_tsv(input$info_tsv$datapath, show_col_types = FALSE)
  })

  # Reactive: results ------------------------------------------------------
  res_obj <- reactive({
    make_results(well_raw(), info_raw(), trim_outlier = input$trim)
  })

  # Outputs ----------------------------------------------------------------
  output$raw_well <- renderDT({
    datatable(well_raw(), options = list(pageLength = 10))
  })

  output$raw_info <- renderDT({
    datatable(info_raw(), options = list(pageLength = 10))
  })

  output$res_table <- renderDT({
    datatable(res_obj()$replicates, options = list(pageLength = 10))
  })

  output$plot <- renderPlot({
    plot_log2fc(res_obj())
  })

  # Download handler -------------------------------------------------------
  output$dl_plot <- downloadHandler(
    filename = function() {
      paste0("log2FC_plot_", Sys.Date(), if (input$trim) "_trim" else "", ".svg")
    },
    content = function(file) {
      ggsave(file, plot_log2fc(res_obj()), device = "svg", width = 9, height = 6)
    },
    contentType = "image/svg+xml"
  )
}

# ------------------------------- Run --------------------------------------
shinyApp(ui, server)