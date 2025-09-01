# BCA Assay Results Automation app.R

library(shiny)
library(readxl)
library(ggplot2)
library(writexl)

ui <- fluidPage(
  titlePanel("BCA Assay Calculator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "1) Upload Excel (.xlsx)", accept = ".xlsx"),
      numericInput("sheet", "Sheet # to read", value = 1, min = 1),
      tags$hr(),
      h4("Standards"),
      textInput("std_rows",     "Rows (e.g. 1:8)",          value = "1:8"),
      numericInput("std_col1",  "Absorbance col #1 (B→2)", value = 2, min = 1),
      numericInput("std_col2",  "Absorbance col #2 (C→3)", value = 3, min = 1),
      textInput("std_conc",     "Concentrations (mg/mL)",  value = "0,0.25,0.5,1,2,4,8,16"),
      tags$hr(),
      h4("Samples"),
      textInput("smp_rows",     "Rows (e.g. 1:8)",          value = "1:8"),
      numericInput("smp_col1",  "Abs col #1 (D→4)",         value = 4, min = 1),
      numericInput("smp_col2",  "Abs col #2 (E→5)",         value = 5, min = 1),
      numericInput("smp_idcol", "SampleID col # (N→14)",    value = 14, min = 1),
      tags$hr(),
      h4("Outputs"),
      textInput("out_xlsx",   "Excel output filename",   value = "BCA_Assay_Results.xlsx"),
      textInput("out_svg",    "Plot SVG filename",       value = "BCA_Assay_Standard_Curve.svg"),
      tags$hr(),
      downloadButton("dl_xlsx", "Download Results (.xlsx)"),
      downloadButton("dl_svg",  "Download Plot (SVG)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Raw Data",    dataTableOutput("raw_tbl")),
        tabPanel("Results",     dataTableOutput("res_tbl")),
        tabPanel("Standard",    dataTableOutput("std_tbl")),
        tabPanel("Plot",        plotOutput("std_plot", height="400px"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive raw data
  raw_df <- reactive({
    req(input$file)
    read_excel(input$file$datapath, sheet = input$sheet)
  })

  # Parse and validate parameters
  params <- reactive({
    list(
      std_rows = eval(parse(text=input$std_rows)),
      std_col1 = input$std_col1,
      std_col2 = input$std_col2,
      std_conc = as.numeric(strsplit(input$std_conc, ",")[[1]]),
      smp_rows = eval(parse(text=input$smp_rows)),
      smp_col1 = input$smp_col1,
      smp_col2 = input$smp_col2,
      smp_id   = input$smp_idcol,
      out_xlsx = input$out_xlsx,
      out_svg  = input$out_svg
    )
  })

  # Build standards table
  std_df <- reactive({
    p <- params(); df <- raw_df()
    abs1 <- df[[ p$std_col1 ]][ p$std_rows ]
    abs2 <- df[[ p$std_col2 ]][ p$std_rows ]
    d <- data.frame(
      Conc = p$std_conc,
      Abs1 = abs1,
      Abs2 = abs2
    )
    d$MeanAbs  <- rowMeans(d[,c("Abs1","Abs2")])
    blank      <- d$MeanAbs[ d$Conc == 0 ]
    d$CorrAbs  <- d$MeanAbs - blank
    d
  })

  # Fit model & build results table
  res_df <- reactive({
    p <- params(); df <- raw_df(); dstd <- std_df()
    blank <- dstd$MeanAbs[ dstd$Conc == 0 ]
    fit   <- lm(CorrAbs ~ Conc, data = dstd)
    slope     <- coef(fit)["Conc"]
    intercept <- coef(fit)["(Intercept)"]

    abs1 <- df[[ p$smp_col1 ]][ p$smp_rows ]
    abs2 <- df[[ p$smp_col2 ]][ p$smp_rows ]
    out <- data.frame(
      SampleID = df[[ p$smp_id ]][ p$smp_rows ],
      MeanAbs  = rowMeans(data.frame(abs1,abs2))
    )
    out$CorrAbs         <- out$MeanAbs - blank
    out$Conc_mg_per_mL  <- (out$CorrAbs - intercept) / slope
    out
  })

  # Render tables & plot
  output$raw_tbl  <- renderDataTable({ raw_df() }, options=list(scrollX=TRUE))
  output$std_tbl  <- renderDataTable({ std_df() }, options=list(scrollX=TRUE))
  output$res_tbl  <- renderDataTable({ res_df() }, options=list(scrollX=TRUE))
  output$std_plot <- renderPlot({
    dstd <- std_df(); dout <- res_df()
    fit   <- lm(CorrAbs ~ Conc, data = dstd)
    slope     <- coef(fit)["Conc"]
    intercept <- coef(fit)["(Intercept)"]

    plot_df <- rbind(
      data.frame(Conc=dstd$Conc, CorrAbs=dstd$CorrAbs, Type="Standard"),
      data.frame(Conc=dout$Conc_mg_per_mL, CorrAbs=dout$CorrAbs, Type="Sample")
    )
    ggplot(plot_df, aes(x=Conc,y=CorrAbs,color=Type)) +
      geom_point(size=3) +
      geom_abline(slope=slope, intercept=intercept, linetype="dashed") +
      labs(
        title="BCA Assay Standard Curve with Sample Points",
        x="Protein concentration (mg/mL)", y="Blank-corrected A (562 nm)"
      ) +
      theme_minimal() +
      theme(legend.title=element_blank())
  })

  # Download handlers
  output$dl_xlsx <- downloadHandler(
    filename = function() params()$out_xlsx,
    content  = function(path) {
      write_xlsx(
        list(
          RAW             = raw_df(),
          Results         = res_df(),
          `Standard Curve`= std_df()
        ),
        path
      )
    }
  )
  output$dl_svg <- downloadHandler(
    filename = function() params()$out_svg,
    content  = function(file) {
      # re-create plot and save
      dstd <- std_df(); dout <- res_df()
      fit   <- lm(CorrAbs ~ Conc, data = dstd)
      slope     <- coef(fit)["Conc"]
      intercept <- coef(fit)["(Intercept)"]
      plot_df <- rbind(
        data.frame(Conc=dstd$Conc, CorrAbs=dstd$CorrAbs, Type="Standard"),
        data.frame(Conc=dout$Conc_mg_per_mL, CorrAbs=dout$CorrAbs, Type="Sample")
      )
      g <- ggplot(plot_df, aes(x=Conc,y=CorrAbs,color=Type)) +
        geom_point(size=3) +
        geom_abline(slope=slope, intercept=intercept, linetype="dashed") +
        theme_minimal()
      ggsave(file, plot=g, device="svg", width=8, height=6, units="in")
    }
  )
}

shinyApp(ui, server)