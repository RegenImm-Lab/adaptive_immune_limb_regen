# ---  Load required packages ---
library(survival)
library(survminer)

# --- Read your data ---
F0_survival <- read.table(header = FALSE, text = "
474 0 Rag1-High
474 0 Rag1-High
474 0 Rag1-High
474 1 Rag1-High
140 1 Rag1-High
279 1 Rag1-High
126 1 Rag1-High
474 0 Rag1-High
474 0 Rag1-High
273 1 Rag1-High
230 0 Rag1-High
75 1 Rag1-High
362 1 Rag1-High
122 1 Rag1-High
407 1 Rag1-High
460 0 Rag1-High
82 1 Rag1-High
376 1 Rag1-High
376 1 Rag1-High
397 1 Rag1-High
145 1 Rag1-High
282 1 Rag1-High
460 0 Rag1-High
474 0 Rag1-Low
231 1 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
474 0 Rag1-Low
262 1 Rag1-Low
474 0 Rag1-Low
122 1 Rag1-Low
130 1 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
75 1 Rag1-Low
460 0 Rag1-Low
75 1 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
75 1 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
460 0 Rag1-Low
82 1 Rag1-Low
460 0 Rag1-Low
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
474 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
460 0 WildType
130 1 WildType
130 1 WildType
")

# --- Name the columns ---
colnames(F0_survival) <- c("time", "event", "group")

# --- Verify structure ---
str(F0_survival)
nrow(F0_survival)

# --- Convert group to factor ---
#F0_survival$group <- as.factor(F0_survival$group) this or the next line
F0_survival$group <- factor(F0_survival$group, levels = c("WildType", "Rag1-Low", "Rag1-High"))


# --- Fit the Kaplan-Meier survival model ---
fit <- survfit(Surv(time, event) ~ group, data = F0_survival)

# --- Plot with ggsurvplot ---
p = ggsurvplot(
  fit,
  data = F0_survival,
  #risk.table = TRUE, table not run in the code
  pval = TRUE,
  pval.coord = c(490.0, 0.65), # coordinates of the p value
  conf.int = FALSE,
  #legend = "right",
  legend = c(0.25, 0.25), # coordinates of the legend
  xlab = "Age (days)",
  ylab = "Percentage Survival",
  title = "Survival Curves by Group",
  legend.title = "",
  legend.labs = c("Wild Type (n=58)", "Rag1 Low edit (n=45)", "Rag1 High edit (n=23)"),
  palette = c("goldenrod3", "gray50", "steelblue3"),
  ggtheme = theme_classic(base_family = "Arial")+
    theme(
      plot.margin = margin(10, 65, 10, 10)  # top, right, bottom, left
      # ↑ this increases the right margin (80), shrinking the plot area
    ),
  font.title = c(16, "bold"),       # title font size
  font.x = c(14, "bold"),           # x-axis font size
  font.y = c(14, "bold"),           # y-axis font size
  font.tickslab = c(12, "plain"),   # tick labels font
  font.legend = c(11, "plain"),     # legend font
  risk.table.fontsize = 4,
  surv.scale = "percent" # percentage of survival
)
# ---- Add coord_cartesian to the *plot object* ----
p$plot <- p$plot + coord_cartesian(clip = "off")

# Print
p
