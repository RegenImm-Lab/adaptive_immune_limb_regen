g_3707 = c(93, 98, 20, 79, 4, 10, 70, 65, 4, 17, 2, 52, 0, 1, 25, 0, 42, 85, 83, 47, 37, 64, 17, 61, 41, 20, 1, 0, 91, 0, 0, 0, 0, 23, 0, 12, 73, 13, 11, 1, 0, 56, 0, 4, 8, 0, 0, 0, 68, 0, 60, 50, 16, 87, 1, 9, 2, 43, 63, 0, 0, 35, 2, 0, 28, 78, 0)
g_5008 = c(54, 0, 68, 68, 99, 34, 14, 87, 22, 28, 11, 89, 56, 88, 15, 60, 37, 91, 97, 24, 90, 84, 53, 47, 61, 60, 26, 41, 94, 74, 54, 6, 79, 12, 25, 94, 48, 98, 24, 23, 93, 30, 84, 95, 22, 37, 51, 97, 69, 53, 23, 83, 98, 59, 73, 87, 80, 94, 9, 60, 97, 79, 30, 94, 99, 54)

#to perform T test
t_result = t.test(g_3707, g_5008)
print(t_result)


# to make a boxplot
boxplot(g_3707, g_5008,
        names = c("gRNA 3707", "gRNA 5008"),
        main = "Knockout scores of two gRNAs to P. waltl's  Rag1",
        xlab = "",
        ylab = "Score",
        col = c("maroon4", "turquoise4"))

# to make a violin plot
install.packages("vioplot")
library(vioplot)

vioplot(g_3707, g_5008,
        names = c("gRNA 3707", "gRNA 5008"),
        col = c("goldenrod3", "steelblue3"),
        rectCol = NA,
        main = "Knockout scores of two gRNAs to P. waltl's  Rag1",
        xlab = "",
        ylab = "")
mtext(side=2, line=3, "Knockout Score", font=2, cex=1.2)
points(jitter(rep(1, length(g_3707)), amount =0.1), g_3707, col = "black", pch = 19, cex = 0.5)
points(jitter(rep(2, length(g_5008)), amount =0.1), g_5008, col = "black", pch = 19, cex = 0.5)

for(i in 1:length(g_3707)) {
  segments(x0 = jitter(1, amount = 0.1), y0 = g_3707[i], 
           x1 = jitter(2, amount = 0.1), y1 = g_5008[i], 
           col = "gray", lty = 1)  # Adjust color and line type as needed
}

install.packages("broom")
library(broom)
tidy(t_result)
  

#dotchart(g_3707, g_5008,
        #names = c("gRNA 3707", "gRNA 5008"),
        #col = c("goldenrod3", "steelblue3"),
        #main = "Comparing editing efficiency of two gRNAs to P. waltl Rag1",
        #xlab = "",
        #ylab = "Score")
  
  
#gRNA editing efficiency irrespective of the knockout score (independent sampling)
edit_3707 = c(93, 98, 24, 79, 7, 10, 72, 65, 12, 18, 2, 52, 17, 1, 25, 5, 42, 85, 83, 47, 37, 64, 18, 73, 47, 26, 1, 15, 91, 16, 0, 0, 31, 23, 0, 12, 73, 15, 11, 1, 0, 70, 0, 9, 89, 0, 0, 1, 68, 4, 60, 52, 20, 87, 1, 9, 3, 45, 73, 0, 23, 35, 2, 0,28, 78, 0)
edit_5008 = c(92, 100, 73, 68, 99, 34, 19, 87, 25, 82, 11, 90, 56, 88, 15, 60, 46, 91, 97, 85, 90, 84, 53, 47, 88, 67, 41, 41, 94, 74, 55, 6, 79, 12, 68, 94, 52, 98, 46, 23, 93, 50, 84, 95, 29, 37, 51, 97, 69, 56, 32, 85, 98, 59, 81, 87, 80, 94, 9, 60, 97, 84, 32, 94, 99, 66)
t_edit = t.test(edit_3707, edit_5008)
print(t_edit)
tidy(t_edit)

# to make a boxplot
boxplot(edit_3707, edit_5008,
        names = c("gRNA 3707", "gRNA 5008"),
        main = "Comparing editing efficiency of two gRNAs to P. waltl Rag1",
        xlab = "",
        ylab = "Indel Score",
        col = c("goldenrod3", "steelblue3"))


# to make a violin plot
vioplot(edit_3707, edit_5008,
        names = c("gRNA 3707", "gRNA 5008"),
        col = c("goldenrod3", "steelblue3"),
        main = "Editing efficiency of two gRNAs to P. waltl Rag1",
        xlab = "",
        ylab = "")
mtext(side=2, line=3, "Indel Score", font=2, cex=1.2)



#to perform a paired T test: I have removed a value from g_3707 because no corresponding value to g_5008
g_3707_pair = c(93, 98, 24, 79, 7, 10, 72, 65, 12, 18, 2, 52, 17, 1, 25, 5, 42, 85, 83, 47, 37, 64, 18, 73, 47, 26, 1, 15, 91, 16, 0, 0, 31, 0, 12, 73, 15, 11, 1, 0, 70, 0, 9, 89, 0, 0, 1, 68, 4, 60, 52, 20, 87, 1, 9, 3, 45, 73, 0, 23, 35, 2, 0, 28, 78, 0)
g_5008_pair = c(92, 100, 73, 68, 99, 34, 19, 87, 25, 82, 11, 90, 56, 88, 15, 60, 46, 91, 97, 85, 90, 84, 53, 47, 88, 67, 41, 41, 94, 74, 55, 6, 79, 12, 68, 94, 52, 98, 46, 23, 93, 50, 84, 95, 29, 37, 51, 97, 69, 56, 32, 85, 98, 59, 81, 87, 80, 94, 9, 60, 97, 84, 32, 94, 99, 66)

t_edit_paired = t.test(g_3707_pair, g_5008_pair, paired = TRUE)
print(t_edit_paired)

# to make a violin plot
install.packages("vioplot")
library(vioplot)

vioplot(g_3707_pair, g_5008_pair,
        names = c("gRNA_3707", "gRNA_5008"),
        col = c("goldenrod3", "steelblue3"),
        rectCol = NA,
        main = "Editing efficiency of two gRNAs to P. waltl's  Rag1",
        xlab = "",
        ylab = "")
mtext(side=2, line=3, "Indel Score", font=2, cex=1.2)
points(jitter(rep(1, length(g_3707_pair)), amount =0.1), g_3707_pair, col = "black", pch = 19, cex = 0.5)
points(jitter(rep(2, length(g_5008_pair)), amount =0.1), g_5008_pair, col = "black", pch = 19, cex = 0.5)

for(i in 1:length(g_3707_pair)) {
  segments(x0 = jitter(1, amount = 0.1), y0 = g_3707_pair[i], 
           x1 = jitter(2, amount = 0.1), y1 = g_5008_pair[i], 
           col = "gray", lty = 1)  # Adjust color and line type as needed
}

# to arrange the test result statistics values in a table
install.packages("broom")
library(broom)
tidy(t_edit_paired)


# to make correlation for the indels 
Indel_3707 <- c(93, 98, 24, 79, 7, 10, 72, 65, 12, 18, 2, 52, 17, 1, 25, 5, 42, 85, 83, 47, 37, 64, 18, 73, 47, 26, 1, 15, 91, 16, 0, 0, 31, 0, 12, 73, 15, 11, 1, 0, 70, 0, 9, 89, 0, 0, 1, 68, 4, 60, 52, 20, 87, 1, 9, 3, 45, 73, 0, 23, 35, 2, 0, 28, 78, 0)
Indel_5008 <- c(92, 100, 73, 68, 99, 34, 19, 87, 25, 82, 11, 90, 56, 88, 15, 60, 46, 91, 97, 85, 90, 84, 53, 47, 88, 67, 41, 41, 94, 74, 55, 6, 79, 12, 68, 94, 52, 98, 46, 23, 93, 50, 84, 95, 29, 37, 51, 97, 69, 56, 32, 85, 98, 59, 81, 87, 80, 94, 9, 60, 97, 84, 32, 94, 99, 66)

# --- Combine into a dataframe ---
df <- data.frame(Indel_3707, Indel_5008)

# --- Compute correlation ---
cor_result <- cor.test(df$Indel_3707, df$Indel_5008, method = "pearson")

# --- Print correlation summary ---
print(cor_result)

# --- Create correlation plot ---
ggplot(df, aes(x = Indel_3707, y = Indel_5008)) +
  geom_point(color = "black", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "brown", se = TRUE, linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Correlation: Indel_3707 and Indel_5008",
    x = "Indel_3707",
    y = "Indel_5008",
    caption = paste0("Pearson r = ", round(cor_result$estimate, 3),
                     ", p = ", signif(cor_result$p.value, 3))
  ) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(size = 12, hjust = 0.5)
  )









KO_3707 = c(93, 98, 20, 79, 4, 10, 70, 65, 4, 17, 2, 52, 0, 1, 25, 0, 42, 85, 83, 47, 37, 64, 17, 61, 41, 20, 1, 0, 91, 0, 0, 0, 0, 0, 12, 73, 13, 11, 1, 0, 56, 0, 4, 8, 0, 0, 0, 68, 0, 60, 50, 16, 87, 1, 9, 2, 43, 63, 0, 0, 35, 2, 0, 28, 78, 0)
KO_5008 = c(54, 0, 68, 68, 99, 34, 14, 87, 22, 28, 11, 89, 56, 88, 15, 60, 37, 91, 97, 24, 90, 84, 53, 47, 61, 60, 26, 41, 94, 74, 54, 6, 79, 12, 25, 94, 48, 98, 24, 23, 93, 30, 84, 95, 22, 37, 51, 97, 69, 53, 23, 83, 98, 59, 73, 87, 80, 94, 9, 60, 97, 79, 30, 94, 99, 54)

# --- Combine into a dataframe ---
df_KO <- data.frame(KO_3707, KO_5008)

# --- Compute correlation ---
cor_result_KO <- cor.test(df_KO$KO_3707, df_KO$KO_5008, method = "pearson")

# --- Print correlation summary ---
print(cor_result_KO)

# --- Create correlation plot ---
ggplot(df_KO, aes(x = KO_3707, y = KO_5008)) +
  geom_point(color = "black", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "brown", se = TRUE, linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Correlation of knockouts: 3707 and 5008",
    x = "3707 knockout",
    y = "5008 knockout",
    caption = paste0("Pearson r = ", round(cor_result_KO$estimate, 3),
                     ", p = ", signif(cor_result_KO$p.value, 3))
  ) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(size = 12, hjust = 0.5)
  )


# --- Combine and compute pairwise sums ---
df_Indel <- data.frame(Indel_3707, Indel_5008)
df_Indel$sum_pair <- Indel_3707 + Indel_5008

# --- Create violin plot ---
ggplot(df_Indel, aes(x = "", y = sum_pair)) +
  geom_violin(fill = "steelblue3", color = "gray30", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "goldenrod3", color = "gray30", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.05, alpha = 0.6, color = "black", size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution: Pairwise Sums of Indels",
    x = "",
    y = "Sum of paired values"
  ) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# --- Combine into dataframe and compute sums ---
df_Indel <- data.frame(
  Indel_3707 = Indel_3707,
  Indel_5008 = Indel_5008,
  Sum = Indel_3707 + Indel_5008
)

# --- Display the first few rows ---
head(df_Indel)

# --- Display all the rows ---
print(df_Indel)



# --- Combine and compute pairwise sums ---
df_KO <- data.frame(KO_3707, KO_5008)
df_KO$sum_pair <- KO_3707 + KO_5008

# --- Create violin plot ---
ggplot(df_KO, aes(x = "", y = sum_pair)) +
  geom_violin(fill = "steelblue3", color = "gray30", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "goldenrod3", color = "gray30", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.05, alpha = 0.6, color = "black", size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution: Pairwise Sums of Knockouts",
    x = "",
    y = "Sum of paired values"
  ) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )




# --- Combine into dataframe and compute sums ---
df_KO <- data.frame(
  KO_3707 = KO_3707,
  KO_5008 = KO_5008,
  Sum = KO_3707 + KO_5008
)

# --- Display the first few rows ---
head(df_KO)

print(df_KO)
# --- Optionally, save to CSV ---
write.csv(df_KO, "KO_sum_dataset.csv", row.names = FALSE)


# --- We want to put colour coding to the datapoints in the violin plot based on our conditions ---

library(tidyverse)

# --- Combine into dataframe ---
df_KO <- data.frame(
  KO_3707 = KO_3707,
  KO_5008 = KO_5008
)

# --- Compute pairwise sum ---
df_KO$sum_pair <- df_KO$KO_3707 + df_KO$KO_5008

# --- Apply the classification conditions ---
df_KO <- df_KO %>%
  mutate(
    category = case_when(
      KO_3707 >= 90 | KO_5008 >= 90 ~ "Any value ≥ 90",
      sum_pair >= 140 ~ "Sum ≥ 140, both values < 90",
      TRUE ~ "Sum < 140, both values < 90"
    )
  )

# --- Define the factor order (legend order) ---
df_KO$category <- factor(
  df_KO$category,
  levels = c(
    "Any value ≥ 90",
    "Sum ≥ 140, both values < 90",
    "Sum < 140, both values < 90"
  )
)

# --- Map colors to the new labels ---
color_map <- c(
  "Any value ≥ 90" = "steelblue4",
  "Sum ≥ 140, both values < 90" = "steelblue1",
  "Sum < 140, both values < 90" = "gray60"
)

# --- Violin plot ---
ggplot(df_KO, aes(x = "", y = sum_pair)) +
  geom_violin(fill = "gray95", color = "gray30", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "gray30", outlier.shape = NA) +
  geom_jitter(aes(color = category), width = 0.08, size = 2.5, alpha = 0.9) +
  scale_color_manual(values = color_map, drop = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    title = "        Pairwise sum of categories with conditions",
    x = "",
    y = "Sum of paired knockout scores",
    color = "Conditions"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Statistic knockout
stats_KO_3707 <- c(
  Mean   = mean(KO_3707),
  Median = median(KO_3707),
  SD     = sd(KO_3707)
)

stats_KO_5008 <- c(
  Mean   = mean(KO_5008),
  Median = median(KO_5008),
  SD     = sd(KO_5008)
)

stats_KO_3707
stats_KO_5008

KO_paired_diff <- stats_KO_5008 - stats_KO_3707

median(KO_paired_diff)   # median difference
mean(KO_paired_diff)     # mean difference
sd(KO_paired_diff)       # SD of differences


# Statistic Indel
stats_Indel_3707 <- c(
  Mean   = mean(KO_3707),
  Median = median(KO_3707),
  SD     = sd(KO_3707)
)

stats_Indel_5008 <- c(
  Mean   = mean(KO_5008),
  Median = median(KO_5008),
  SD     = sd(KO_5008)
)

stats_Indel_3707
stats_Indel_5008


Indel_paired_diff <- stats_Indel_5008 - stats_Indel_3707

median(Indel_paired_diff)   # median difference
mean(Indel_paired_diff)     # mean difference
sd(Indel_paired_diff)       # SD of differences
