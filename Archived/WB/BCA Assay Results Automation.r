# === BCA Assay: automatic concentration calculation + plotting + multi-sheet export ===

# (1) Load required libraries
if (!requireNamespace("readxl", quietly = TRUE))  install.packages("readxl")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")

library(readxl)
library(ggplot2)
library(writexl)

# ————— USER PARAMETERS —————
# Path to your raw Excel file
setwd("/Users/SofTicE/Desktop/Research Project/WB/20250612 Sen Exp 65 BCA/") # Set working directory to where your file is located
file_path       <- "20250612 Sen Exp 65 BCA.xlsx"

# Standards live in these rows and columns:
std_rows        <- 1:8                # e.g. rows 1 through 8
std_col_abs1    <- 2                  # column B
std_col_abs2    <- 3                  # column C
std_conc        <- c(0,0.25,0.5,1,2,4,8,16)  # mg/mL, same length as std_rows

# Samples live in these rows and columns:
smp_rows        <- 1:5                # e.g. rows 1 through 5
smp_col_abs1    <- 4                  # column D
smp_col_abs2    <- 5                  # column E
smp_col_id      <- 14                 # column N for SampleID

# Output filenames
out_xlsx        <- "BCA_Assay_Results.xlsx"
out_plot_svg    <- "BCA_Assay_Standard_Curve.svg"

# ————— END USER PARAMETERS —————

# (2) Read raw data
raw_df <- read_excel(file_path, sheet = 1)

# (3) Process standards
std_abs1 <- raw_df[[std_col_abs1]][std_rows]
std_abs2 <- raw_df[[std_col_abs2]][std_rows]
std_df   <- data.frame(
  Conc    = std_conc,
  Abs1    = std_abs1,
  Abs2    = std_abs2
)
std_df$MeanAbs  <- rowMeans(std_df[, c("Abs1","Abs2")])
blank          <- std_df$MeanAbs[std_df$Conc == 0]
std_df$CorrAbs <- std_df$MeanAbs - blank

# (4) Linear fit
fit       <- lm(CorrAbs ~ Conc, data = std_df)
slope     <- coef(fit)["Conc"]
intercept <- coef(fit)["(Intercept)"]
r2        <- summary(fit)$r.squared
cat(sprintf(
  "Standard curve:\n  CorrAbs = %.4f × Conc + %.4f   (R² = %.3f)\n\n",
  slope, intercept, r2
))

# (5) Process samples
smp_abs1 <- raw_df[[smp_col_abs1]][smp_rows]
smp_abs2 <- raw_df[[smp_col_abs2]][smp_rows]
results_df <- data.frame(
  SampleID        = raw_df[[smp_col_id]][smp_rows],
  MeanAbs         = rowMeans(data.frame(smp_abs1, smp_abs2)),
  CorrAbs         = rowMeans(data.frame(smp_abs1, smp_abs2)) - blank
)
results_df$Conc_mg_per_mL <- (results_df$CorrAbs - intercept) / slope

# (6) Plot
plot_df <- rbind(
  data.frame(Conc = std_df$Conc,    CorrAbs = std_df$CorrAbs,    Type = "Standard"),
  data.frame(Conc = results_df$Conc_mg_per_mL, CorrAbs = results_df$CorrAbs, Type = "Sample")
)
p <- ggplot(plot_df, aes(x = Conc, y = CorrAbs, color = Type)) +
  geom_point(size = 3) +
  geom_abline(slope = slope, intercept = intercept, linetype = "dashed") +
  labs(
    title = "BCA Assay Standard Curve with Sample Points",
    x = "Protein concentration (mg/mL)",
    y = "Blank-corrected absorbance (562 nm)"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title   = element_text(hjust = 0.5)
  )

# (7) Save plot as SVG
ggsave(out_plot_svg, plot = p, width = 8, height = 6, units = "in", device = "svg")
cat("✓ Plot saved to:", out_plot_svg, "\n")

# (8) Write multi-sheet Excel
write_xlsx(
  list(
    RAW            = raw_df,
    Results        = results_df,
    `Standard Curve` = std_df
  ),
  path = out_xlsx
)
cat("✓ Data exported to:", out_xlsx, "\n")