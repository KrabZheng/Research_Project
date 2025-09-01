## app.R: Shiny app for qPCR ΔΔCq plotting and SVG download

library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# Helper: greyscale check
grey_check <- function(cols, threshold = 30) {
  lums <- sapply(cols, function(c) {
    rgb <- col2rgb(c)
    0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  })
  diffs <- combn(lums, 2, FUN = function(x) abs(diff(x)))
  if (any(diffs < threshold)) {
    warning("Palette colours differ by < ", threshold,
            " in greyscale—some may be indistinguishable in B/W.")
  }
}

# Helper: generate palette with repetition if needed
get_colors <- function(n, palette, thresh) {
  maxcols <- brewer.pal.info[palette, "maxcolors"]
  if (n <= maxcols) {
    cols <- brewer.pal(n, palette)
  } else {
    base <- brewer.pal(maxcols, palette)
    cols <- rep(base, length.out = n)
    warning("Requested ", n, " colors; palette max is ", maxcols, ". Repeating colors.")
  }
  grey_check(cols, thresh)
  cols
}

# Plot 1: original (Gene × Patient), color by Gene
plot_orig <- function(df, palette, thresh) {
  df2 <- df %>% filter(!is.na(log2FC) & Condition != "Control")
  if (nrow(df2) == 0) return(NULL)
  genes <- unique(df2$Gene)
  cols <- get_colors(length(genes), palette, thresh)
  names(cols) <- genes

  ggplot(df2, aes(x = Condition, y = log2FC, fill = Gene)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = log2FC - se_ddCt, ymax = log2FC + se_ddCt),
                  position = position_dodge(0.8), width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = cols) +
    facet_grid(Gene ~ Patient, scales = "fixed", space = "fixed",
               labeller = labeller(Gene = label_wrap_gen(15))) +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 12) +
    labs(x = "Condition", y = "Log₂ Fold-Change vs. Control") +
    theme(
      aspect.ratio = 1,
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y = element_text(angle = 0, hjust = 1),
      panel.grid.major.y = element_line(),
      legend.position = "none"
    )
}

# Plot 2: By Gene (facet by Gene), colour by Patient, x = Condition, with legend
plot_byGene <- function(df, palette, thresh) {
  df2 <- df %>% filter(!is.na(log2FC) & Condition != "Control")
  if (nrow(df2) == 0) return(NULL)
  patients <- unique(df2$Patient)
  cols <- get_colors(length(patients), palette, thresh)
  names(cols) <- patients

  ggplot(df2, aes(x = Condition, y = log2FC, fill = Patient)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = log2FC - se_ddCt, ymax = log2FC + se_ddCt),
                  position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.5) +
    scale_fill_manual(name = "Patient", values = cols) +
    facet_wrap(~Gene, scales = "fixed", labeller = label_wrap_gen(15)) +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 12) +
    labs(x = "Condition", y = "Log₂ Fold-Change vs. Control") +
    theme(
      aspect.ratio    = 1,
      strip.background= element_blank(),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(),
      legend.position = "right"
    )
}

# Plot 3: By Patient (facet by Patient), colour by Gene, x = Condition, with legend
plot_byPatient <- function(df, palette, thresh) {
  df2 <- df %>% filter(!is.na(log2FC) & Condition != "Control")
  if (nrow(df2) == 0) return(NULL)
  genes <- unique(df2$Gene)
  cols <- get_colors(length(genes), palette, thresh)
  names(cols) <- genes

  ggplot(df2, aes(x = Condition, y = log2FC, fill = Gene)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = log2FC - se_ddCt, ymax = log2FC + se_ddCt),
                  position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.5) +
    scale_fill_manual(name = "Gene", values = cols) +
    facet_wrap(~Patient, scales = "fixed") +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 12) +
    labs(x = "Condition", y = "Log₂ Fold-Change vs. Control") +
    theme(
      aspect.ratio     = 1,
      strip.background = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.text.x     = element_text(angle = 0),
      panel.grid.major.y = element_line(),
      legend.position  = "right"
    )
}

# UI
ui <- fluidPage(
  titlePanel("qPCR ΔΔCq Plotter"),
  sidebarLayout(
    sidebarPanel(
      fileInput("wellFile",  "Upload WellResult CSV",        accept=".csv"),
      fileInput("infoFile",  "Upload qPCR Sample Info TSV", accept=c(".tsv",".txt")),
      selectInput("palette", "Brewer Palette",
                  choices=rownames(brewer.pal.info), selected="YlGnBu"),
      numericInput("greyThresh","Greyscale threshold", value=30, min=0, max=255),
      hr(),
      downloadButton("downloadOrig","Download All (Gene × Patient)"),
      downloadButton("downloadGene","Download by Gene"),
      downloadButton("downloadPatient","Download by Patient")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("All (Gene × Patient)", plotOutput("qcPlot")),
        tabPanel("By Gene", plotOutput("qcPlotGene")),
        tabPanel("By Patient", plotOutput("qcPlotPatient"))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  results <- reactive({
    req(input$wellFile, input$infoFile)
    well <- read_csv(input$wellFile$datapath) %>% filter(!Omit)
    info <- read_tsv(input$infoFile$datapath)
    sample_info <- info %>% filter(Info_type=="Sample") %>%
      select(Plate_Sample=Info, SampleName=Info_name, Patient, Condition)
    target_info <- info %>% filter(Info_type=="Target") %>%
      select(Plate_Target=Info, Gene=Info_name, Role=Info_key)
    dat <- well %>%
      left_join(sample_info, by=c("Sample"="Plate_Sample")) %>%
      left_join(target_info, by=c("Target"="Plate_Target")) %>%
      select(Patient, Condition, Gene, Role, Cq)
    dat <- dat %>% group_by(Patient, Condition, Gene, Role) %>%
      mutate(dev=abs(Cq-median(Cq,na.rm=TRUE)), n_rep=sum(!is.na(Cq))) %>%
      filter(!(n_rep>=3 & dev==max(dev))) %>% ungroup() %>% select(-dev,-n_rep)
    mean_cq <- dat %>% group_by(Patient,Condition,Gene,Role) %>%
      summarise(mean_Cq=mean(Cq,na.rm=TRUE), sd_Cq=sd(Cq,na.rm=TRUE), n=sum(!is.na(Cq)), .groups="drop")
    ref_cq <- mean_cq %>% filter(Role=="Reference") %>%
      select(Patient,Condition,ref_Cq=mean_Cq,sd_ref=sd_Cq,n_ref=n)
    target_cq <- mean_cq %>% filter(Role!="Reference") %>%
      select(Patient,Condition,Gene,targ_Cq=mean_Cq,sd_targ=sd_Cq,n_targ=n)
    dq <- target_cq %>% left_join(ref_cq, by=c("Patient","Condition")) %>%
      mutate(delta_Cq=targ_Cq-ref_Cq,
             sd_delta=sqrt(sd_targ^2+sd_ref^2),
             se_delta=sd_delta/sqrt(pmin(n_targ,n_ref,na.rm=TRUE)))
    ctrl <- dq %>% filter(Condition=="Control") %>%
      select(Patient,Gene,control_delta=delta_Cq,se_control=se_delta)
    dq %>% left_join(ctrl, by=c("Patient","Gene")) %>%
      mutate(delta_delta_Cq=delta_Cq-control_delta,
             se_ddCt=sqrt(se_delta^2+se_control^2),
             log2FC=delta_delta_Cq)
  })

  output$qcPlot <- renderPlot({ plot_orig(results(), input$palette, input$greyThresh) })
  output$qcPlotGene <- renderPlot({ plot_byGene(results(), input$palette, input$greyThresh) })
  output$qcPlotPatient <- renderPlot({ plot_byPatient(results(), input$palette, input$greyThresh) })

  # Download handlers
  save_svg <- function(fun, filename) {
    function(file) {
      p <- fun(results(), input$palette, input$greyThresh)
      ggsave(file, p, device="svg", width=2048, height=2048, units="px")
    }
  }
  output$downloadOrig <- downloadHandler("qpcr_all.svg", save_svg(plot_orig, "qpcr_all.svg"))
  output$downloadGene <- downloadHandler("qpcr_gene.svg", save_svg(plot_byGene, "qpcr_gene.svg"))
  output$downloadPatient <- downloadHandler("qpcr_patient.svg", save_svg(plot_byPatient, "qpcr_patient.svg"))
}

# Launch app
shinyApp(ui, server)