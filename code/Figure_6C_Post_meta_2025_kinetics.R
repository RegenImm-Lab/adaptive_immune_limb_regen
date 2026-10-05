if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
library(tidyverse)


# --- Data from 2025 experiment ---
groups <- c("Rag1 +/+", "Rag1 +/-", "Rag1 -/-")
times  <- c("Early-bud", "Medium-bud", "Palette", "Digit")

values_list <- list(
  "Rag1 +/+" = list(
    "Early-bud" = c(13, 14),
    "Medium-bud" = c(17, 18),
    "Palette" = c(21, 22),
    "Digit" = c(26, 26)
  ),
  "Rag1 +/-" = list(
    "Early-bud" = c(14, 13, 14, 13, 13, 13, 13, 13, 13),
    "Medium-bud" = c(18, 16, 17, 17, 16, 16, 16, 16, 16),
    "Palette" = c(21, 20, 21, 21, 21, 21, 21, 20, 20),
    "Digit" = c(29, 25, 26, 26, 26, 28, 26, 23, 24)
  ),
  "Rag1 -/-" = list(
    "Early-bud" = c(13, 13, 14, 13),
    "Medium-bud" = c(16, 17, 17, 16),
    "Palette" = c(20, 22, 22, 21),
    "Digit" = c(26, 27, 26, 25)
  )
)

data <- map_dfr(groups, function(g) {
  map_dfr(times, function(t) {
    tibble(group = g, time = t, value = values_list[[g]][[t]])
  })
})

# --- Clean & Order Factors ---
data <- data %>%
  mutate(group = str_squish(group),
         time  = str_squish(time)) %>%
  mutate(group = factor(group, levels = unique(group)),
         time  = factor(time, levels = c("Early-bud", "Medium-bud", "Palette", "Digit")))

# --- Bootstrap CI function for median ---
bootstrap_median_ci <- function(x, nboot = 2000, conf = 0.95) {
  boots <- replicate(nboot, median(sample(x, replace = TRUE)))
  ci <- quantile(boots, probs = c((1 - conf)/2, 1 - (1 - conf)/2))
  tibble(median_value = median(x),
         lower = ci[1],
         upper = ci[2])
}

# --- Compute median and 95% CI for each group and time ---
median_df <- data %>%
  group_by(group, time) %>%
  summarise(bootstrap_median_ci(value), .groups = "drop")

print(median_df)

# --- Perform Kruskal–Wallis tests for each time point ---
kw_results <- data %>%
  group_by(time) %>%
  summarise(
    statistic = kruskal.test(value ~ group)$statistic,
    p_value   = kruskal.test(value ~ group)$p.value,
    .groups = "drop"
  )

print(kw_results)

# --- Merge summary and test results for overview ---
summary_table <- median_df %>%
  left_join(kw_results, by = "time") %>%
  select(time, group, median_value, lower, upper, statistic, p_value)

print(summary_table)

# --- Define your preferred colors ---
my_colors <- c("Rag1 +/+" = "goldenrod3",
               "Rag1 +/-" = "gray50",
               "Rag1 -/-" = "steelblue3")

# --- Plot with CI error bars (median time on x-axis) ---
ggplot(median_df, aes(x = median_value, y = time, group = group)) +
  geom_line(aes(color = group), linewidth = 1.2) +
  geom_point(aes(fill = group), shape = 21, color = "black", size = 3, stroke = 0.6) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = group), height = 0.15, linewidth = 0.8) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "  Median Time to Limb Regeneration Events",
    x = "Median time (days)",
    y = "Regeneration event",
    color = "Group",
    fill = "Group"
  ) +
  xlim(0, max(median_df$upper, na.rm = TRUE) + 2) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )






# --- Data from 2025 experiment but combining Hets ans WT---
groups <- c("Rag1 +/+", "Rag1 -/-")
times  <- c("Early-bud", "Medium-bud", "Palette", "Digit")

values_list <- list(
  "Rag1 +/+" = list(
    "Early-bud" = c(13, 14, 14, 13, 14, 13, 13, 13, 13, 13, 13),
    "Medium-bud" = c(17, 18, 18, 16, 17, 17, 16, 16, 16, 16, 16),
    "Palette" = c(21, 22, 21, 20, 21, 21, 21, 21, 21, 20, 20),
    "Digit" = c(26, 26, 29, 25, 26, 26, 26, 28, 26, 23, 24)
  ),
  #"Rag1 +/-" = list(
   # "Early-bud" = c(),
    #"Medium-bud" = c(),
    #"Palette" = c(),
    #"Digit" = c()
  #),
  "Rag1 -/-" = list(
    "Early-bud" = c(13, 13, 14, 13),
    "Medium-bud" = c(16, 17, 17, 16),
    "Palette" = c(20, 22, 22, 21),
    "Digit" = c(26, 27, 26, 25)
  )
)

data <- map_dfr(groups, function(g) {
  map_dfr(times, function(t) {
    tibble(group = g, time = t, value = values_list[[g]][[t]])
  })
})

# --- Clean & Order Factors ---
data <- data %>%
  mutate(group = str_squish(group),
         time  = str_squish(time)) %>%
  mutate(group = factor(group, levels = unique(group)),
         time  = factor(time, levels = c("Early-bud", "Medium-bud", "Palette", "Digit")))

# --- Bootstrap CI function for median ---
bootstrap_median_ci <- function(x, nboot = 2000, conf = 0.95) {
  boots <- replicate(nboot, median(sample(x, replace = TRUE)))
  ci <- quantile(boots, probs = c((1 - conf)/2, 1 - (1 - conf)/2))
  tibble(median_value = median(x),
         lower = ci[1],
         upper = ci[2])
}

# --- Compute median and 95% CI for each group and time ---
median_df <- data %>%
  group_by(group, time) %>%
  summarise(bootstrap_median_ci(value), .groups = "drop")

print(median_df)

# --- Perform Kruskal–Wallis tests for each time point ---
kw_results <- data %>%
  group_by(time) %>%
  summarise(
    statistic = kruskal.test(value ~ group)$statistic,
    p_value   = kruskal.test(value ~ group)$p.value,
    .groups = "drop"
  )

print(kw_results)

# --- Merge summary and test results for overview ---
summary_table <- median_df %>%
  left_join(kw_results, by = "time") %>%
  select(time, group, median_value, lower, upper, statistic, p_value)

print(summary_table)

# --- Define your preferred colors ---
my_colors <- c("Rag1 +/+" = "goldenrod3",
               #"Rag1 +/-" = "gray50",
               "Rag1 -/-" = "steelblue3")

# --- Plot with CI error bars (median time on x-axis) ---
ggplot(median_df, aes(x = median_value, y = time, group = group)) +
  geom_line(aes(color = group), linewidth = 1.2) +
  geom_point(aes(fill = group), shape = 21, color = "black", size = 3, stroke = 0.6) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = group), height = 0.15, linewidth = 0.8) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "  Median Time to Limb Regeneration Events",
    x = "Median time (days)",
    y = "Regeneration event",
    color = "Group",
    fill = "Group"
  ) +
  xlim(0, max(median_df$upper, na.rm = TRUE) + 2) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
