# =============================================================================
# Kelp Coverage Analysis
# Treatment effects on offspring growth (kelp coverage proxy)
# 2x2 factorial: Nutrient (High/Low) x Light (High/Low)
# Linear mixed model with cross as random effect
#
# kelp_coverage = fraction of ROI pixels with R:G ratio > 1.05
# Brown/gold kelp against white background — higher coverage = more growth
# =============================================================================

setwd("/Users/harrysilke/Desktop/UNI/Honours/imas_hons/kelp-photograph-pipeline")

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(scales)

# =============================================================================
# 1. Load and prepare data
# =============================================================================

raw <- read_csv("results/aggregated_20260402_123846.csv")

master <- readxl::read_excel("pilot_master_sheet.xlsx") %>%
  rename(
    cross_id  = `Cross ID...1`,
    nutrient  = `Nutrient treatment`,
    light     = `Light treatment`,
    female_id = `Parent 1 - culture ID`,
    male_id   = `Parent 2 - culture ID`
  ) %>%
  select(cross_id, nutrient, light, female_id, male_id) %>%
  mutate(
    cross_combo = paste0(female_id, " x ", male_id),
    nutrient    = factor(nutrient, levels = c("Low", "High")),
    light       = factor(light,    levels = c("Low", "High"))
  )

# kelp_coverage is identical across channels — keep one row per cross x date x ROI
df <- raw %>%
  filter(channel == "green") %>%
  mutate(
    cross_id  = as.integer(str_remove(sample_id, "C")),
    date      = factor(date_readable,
                       levels = c("10.03.2026", "17.03.2026", "27.03.2026"),
                       ordered = TRUE),
    timepoint = as.numeric(date)
  ) %>%
  left_join(master, by = "cross_id") %>%
  mutate(cross_id = factor(cross_id)) %>%
  filter(!is.na(date), !is.na(nutrient), !is.na(light))

cat("Crosses per treatment cell:\n")
df %>%
  distinct(cross_id, nutrient, light) %>%
  count(nutrient, light) %>%
  print()

cat("\nKelp coverage summary:\n")
print(summary(df$kelp_coverage))

# =============================================================================
# 2. Visualise: coverage over time by treatment
# =============================================================================

treatment_means <- df %>%
  group_by(date, nutrient, light) %>%
  summarise(
    mean_cov = mean(kelp_coverage, na.rm = TRUE),
    se       = sd(kelp_coverage, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light))

ggplot() +
  geom_line(
    data = df %>% mutate(treatment = paste0("N:", nutrient, " L:", light)),
    aes(x = date, y = kelp_coverage, group = cross_id, colour = treatment),
    alpha = 0.2
  ) +
  geom_line(
    data = treatment_means,
    aes(x = date, y = mean_cov, colour = treatment, group = treatment),
    linewidth = 1.4
  ) +
  geom_errorbar(
    data = treatment_means,
    aes(x = date, ymin = mean_cov - se, ymax = mean_cov + se, colour = treatment),
    width = 0.2, linewidth = 1
  ) +
  geom_point(
    data = treatment_means,
    aes(x = date, y = mean_cov, colour = treatment),
    size = 3
  ) +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent) +
  labs(
    title    = "Kelp coverage over time by treatment",
    subtitle = "Thin lines = individual crosses, thick lines = treatment means ± SE",
    x        = "Date",
    y        = "Kelp coverage (% of ROI)",
    colour   = "Treatment"
  ) +
  theme_bw(base_size = 13)

ggsave("results/coverage_over_time.png", width = 10, height = 6, dpi = 150)
cat("Plot saved: results/coverage_over_time.png\n")

# =============================================================================
# 3. Mixed model: coverage ~ nutrient * light * timepoint + (1 | cross_id)
# =============================================================================

m_full <- tryCatch(
  lmer(
    kelp_coverage ~ nutrient * light * timepoint + (1 + timepoint | cross_id),
    data = df, REML = TRUE
  ),
  error = function(e) {
    cat("Full model failed — simplifying to random intercept only\n")
    lmer(
      kelp_coverage ~ nutrient * light * timepoint + (1 | cross_id),
      data = df, REML = TRUE
    )
  }
)

if (isSingular(m_full)) {
  cat("Full model singular — simplifying to random intercept only\n")
  m_full <- lmer(
    kelp_coverage ~ nutrient * light * timepoint + (1 | cross_id),
    data = df, REML = TRUE
  )
}

cat("\n=== Model summary ===\n")
print(summary(m_full))

cat("\n=== ANOVA table (Type III) ===\n")
print(anova(m_full))

# =============================================================================
# 4. Post-hoc: pairwise treatment comparisons
# =============================================================================

emm <- emmeans(m_full, ~ nutrient * light | timepoint, at = list(timepoint = 1:3))

cat("\n=== Estimated marginal means by treatment and timepoint ===\n")
print(emm)

cat("\n=== Pairwise contrasts ===\n")
print(pairs(emm, adjust = "tukey"))

# =============================================================================
# 5. Trajectory plot with model estimates
# =============================================================================

emm_df <- as.data.frame(emm) %>%
  mutate(
    treatment = paste0("N:", nutrient, " L:", light),
    date = factor(
      c("10.03.2026", "17.03.2026", "27.03.2026")[timepoint],
      levels = c("10.03.2026", "17.03.2026", "27.03.2026")
    )
  )

raw_means <- df %>%
  group_by(nutrient, light, date) %>%
  summarise(
    mean_cov = mean(kelp_coverage, na.rm = TRUE),
    se       = sd(kelp_coverage, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light))

ggplot() +
  geom_point(
    data = raw_means,
    aes(x = date, y = mean_cov, colour = treatment, group = treatment),
    size = 2, alpha = 0.6, position = position_dodge(0.1)
  ) +
  geom_errorbar(
    data = raw_means,
    aes(x = date, ymin = mean_cov - se, ymax = mean_cov + se,
        colour = treatment, group = treatment),
    width = 0.15, alpha = 0.6, position = position_dodge(0.1)
  ) +
  geom_line(
    data = emm_df,
    aes(x = date, y = emmean, colour = treatment, group = treatment),
    linewidth = 1.3
  ) +
  geom_ribbon(
    data = emm_df,
    aes(x = date, ymin = lower.CL, ymax = upper.CL,
        fill = treatment, group = treatment),
    alpha = 0.15
  ) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent) +
  labs(
    title    = "Kelp coverage trajectories by treatment",
    subtitle = "Lines = model estimates (95% CI), points = raw means ± SE",
    x        = "Date",
    y        = "Kelp coverage (% of ROI)",
    colour   = "Treatment",
    fill     = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("results/coverage_trajectory.png", width = 10, height = 6, dpi = 150)
cat("Plot saved: results/coverage_trajectory.png\n")

# =============================================================================
# 6. Faceted plot: one panel per treatment
# =============================================================================

df %>%
  mutate(treatment = paste0("Nutrient: ", nutrient, "\nLight: ", light)) %>%
  ggplot(aes(x = date, y = kelp_coverage, group = cross_id)) +
  geom_line(alpha = 0.3, colour = "steelblue") +
  geom_point(alpha = 0.4, colour = "steelblue", size = 1.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               colour = "black", linewidth = 1.4) +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               colour = "black", size = 3) +
  facet_wrap(~ treatment, ncol = 2) +
  scale_y_continuous(labels = percent) +
  labs(
    title    = "Individual cross trajectories by treatment",
    subtitle = "Thin lines = individual crosses, thick line = treatment mean",
    x        = "Date",
    y        = "Kelp coverage (% of ROI)"
  ) +
  theme_bw(base_size = 13)

ggsave("results/coverage_faceted.png", width = 10, height = 8, dpi = 150)
cat("Plot saved: results/coverage_faceted.png\n")

# =============================================================================
# 7. Cross combination analysis
# =============================================================================

combo_means <- df %>%
  group_by(cross_combo, nutrient, light, date) %>%
  summarise(
    mean_cov = mean(kelp_coverage, na.rm = TRUE),
    se       = sd(kelp_coverage, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light))

ggplot(combo_means, aes(x = date, y = mean_cov,
                         colour = treatment, group = treatment)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_cov - se, ymax = mean_cov + se),
                width = 0.2, linewidth = 0.8) +
  facet_wrap(~ cross_combo, ncol = 3) +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent) +
  labs(
    title    = "Kelp coverage by cross combination and treatment",
    subtitle = "Mean ± SE per treatment per timepoint",
    x        = "Date",
    y        = "Kelp coverage (% of ROI)",
    colour   = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 30, hjust = 1)
  )

ggsave("results/coverage_cross_combo.png", width = 12, height = 7, dpi = 150)
cat("Plot saved: results/coverage_cross_combo.png\n")

# Heatmap — use most recent available timepoint per cross
final_tp <- df %>%
  group_by(cross_id, cross_combo, nutrient, light) %>%
  filter(date == max(date)) %>%
  ungroup() %>%
  group_by(cross_combo, nutrient, light) %>%
  summarise(mean_cov = mean(kelp_coverage, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(cross_combo)) %>%
  mutate(treatment = paste0("N:", nutrient, "\nL:", light))

ggplot(final_tp, aes(x = treatment, y = cross_combo, fill = mean_cov)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = percent(mean_cov, accuracy = 0.1)), size = 4) +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = median(final_tp$mean_cov, na.rm = TRUE),
    name     = "Coverage",
    labels   = percent
  ) +
  labs(
    title    = "Kelp coverage at most recent timepoint",
    subtitle = "By cross combination and treatment",
    x        = "Treatment",
    y        = "Cross combination"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

ggsave("results/coverage_heatmap.png", width = 9, height = 5, dpi = 150)
cat("Plot saved: results/coverage_heatmap.png\n")

# =============================================================================
# 8. Rate of change: 10.03 to most recent timepoint
# =============================================================================

change <- df %>%
  group_by(cross_id, date, nutrient, light, cross_combo) %>%
  summarise(kelp_coverage = mean(kelp_coverage, na.rm = TRUE), .groups = "drop") %>%
  filter(date %in% c("10.03.2026", "27.03.2026")) %>%
  pivot_wider(names_from = date, values_from = kelp_coverage) %>%
  mutate(change = `27.03.2026` - `10.03.2026`) %>%
  filter(!is.na(change))

change %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light)) %>%
  ggplot(aes(x = treatment, y = change, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent) +
  labs(
    title    = "Change in kelp coverage: 10.03 to 27.03",
    subtitle = "Positive = growth, negative = decline",
    x        = "Treatment",
    y        = "Change in coverage"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("results/coverage_change.png", width = 9, height = 6, dpi = 150)
cat("Plot saved: results/coverage_change.png\n")

# =============================================================================
# 9. Boxplots at each timepoint
# =============================================================================

df %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light)) %>%
  ggplot(aes(x = date, y = kelp_coverage, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent) +
  labs(
    title  = "Distribution of kelp coverage by treatment over time",
    x      = "Date",
    y      = "Kelp coverage (% of ROI)",
    fill   = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("results/coverage_boxplots.png", width = 10, height = 6, dpi = 150)
cat("Plot saved: results/coverage_boxplots.png\n")

# =============================================================================
# 10. Population-level performance
# =============================================================================

pop_long <- df %>%
  filter(!is.na(cross_combo)) %>%
  mutate(
    pop1 = str_extract(cross_combo, "^[^\\s]+"),
    pop2 = str_extract(cross_combo, "[^\\s]+$")
  ) %>%
  pivot_longer(cols = c(pop1, pop2), names_to = "parent", values_to = "population")

pop_long %>%
  group_by(population, date, nutrient, light) %>%
  summarise(
    mean_cov = mean(kelp_coverage, na.rm = TRUE),
    se       = sd(kelp_coverage, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light)) %>%
  ggplot(aes(x = date, y = mean_cov, colour = treatment, group = treatment)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_cov - se, ymax = mean_cov + se),
                width = 0.2, linewidth = 0.7) +
  facet_wrap(~ population, ncol = 3) +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent) +
  labs(
    title    = "Kelp coverage by population of origin",
    subtitle = "Each population appears in crosses as either parent",
    x        = "Date",
    y        = "Kelp coverage (% of ROI)",
    colour   = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 30, hjust = 1)
  )

ggsave("results/coverage_by_population.png", width = 12, height = 7, dpi = 150)
cat("Plot saved: results/coverage_by_population.png\n")

# =============================================================================
# 11. Correlation: 10.03 vs 27.03 coverage per cross
# =============================================================================

tp_wide <- df %>%
  group_by(cross_id, date, nutrient, light, cross_combo) %>%
  summarise(kelp_coverage = mean(kelp_coverage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = date, values_from = kelp_coverage) %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light)) %>%
  filter(!is.na(`10.03.2026`) & !is.na(`27.03.2026`))

ggplot(tp_wide, aes(x = `10.03.2026`, y = `27.03.2026`, colour = treatment)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  labs(
    title    = "Kelp coverage: 10.03 vs 27.03 per cross",
    subtitle = "Dashed line = no change; points above = growth",
    x        = "Coverage at 10.03.2026",
    y        = "Coverage at 27.03.2026",
    colour   = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("results/coverage_correlation.png", width = 9, height = 7, dpi = 150)
cat("Plot saved: results/coverage_correlation.png\n")

cat("\nAll plots saved to results/\n")


# Reflectance over time by treatment
raw %>%
  filter(channel == "green") %>%
  mutate(
    cross_id  = as.integer(str_remove(sample_id, "C")),
    date      = factor(date_readable,
                       levels = c("10.03.2026", "17.03.2026", "27.03.2026"),
                       ordered = TRUE)
  ) %>%
  left_join(master, by = "cross_id") %>%
  filter(!is.na(date), !is.na(nutrient), !is.na(light)) %>%
  group_by(date, nutrient, light) %>%
  summarise(
    mean_ref = mean(mean_reflectance, na.rm = TRUE),
    se       = sd(mean_reflectance, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(treatment = paste0("N:", nutrient, " L:", light)) %>%
  ggplot(aes(x = date, y = mean_ref, colour = treatment, group = treatment)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_ref - se, ymax = mean_ref + se),
                width = 0.2, linewidth = 1) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title  = "Mean green reflectance over time by treatment",
    x      = "Date",
    y      = "Mean green reflectance (%)",
    colour = "Treatment"
  ) +
  theme_bw(base_size = 13)

ggsave("results/reflectance_green_trajectory.png", width = 10, height = 6, dpi = 150)

# Reflectance heatmap for all 3 timepoints
raw %>%
  filter(channel == "green") %>%
  mutate(
    cross_id  = as.integer(str_remove(sample_id, "C")),
    date      = factor(date_readable,
                       levels = c("10.03.2026", "17.03.2026", "27.03.2026"),
                       ordered = TRUE)
  ) %>%
  left_join(master, by = "cross_id") %>%
  filter(!is.na(date), !is.na(nutrient), !is.na(light), !is.na(cross_combo)) %>%
  group_by(cross_combo, nutrient, light, date) %>%
  summarise(mean_ref = mean(mean_reflectance, na.rm = TRUE), .groups = "drop") %>%
  mutate(treatment = paste0("N:", nutrient, "\nL:", light)) %>%
  ggplot(aes(x = treatment, y = cross_combo, fill = mean_ref)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = round(mean_ref, 1)), size = 3) +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = 50,
    name     = "Reflectance (%)"
  ) +
  facet_wrap(~ date, ncol = 3) +
  labs(
    title    = "Mean green reflectance by cross combination and treatment",
    subtitle = "Faceted by timepoint",
    x        = "Treatment",
    y        = "Cross combination"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.text.x  = element_text(angle = 30, hjust = 1),
    strip.text   = element_text(size = 10)
  )

ggsave("results/reflectance_heatmap_all_timepoints.png", width = 15, height = 6, dpi = 150)
cat("Plot saved: results/reflectance_heatmap_all_timepoints.png\n")


# Coverage heatmap for all 3 timepoints
df %>%
  filter(!is.na(cross_combo)) %>%
  group_by(cross_combo, nutrient, light, date) %>%
  summarise(mean_cov = mean(kelp_coverage, na.rm = TRUE), .groups = "drop") %>%
  mutate(treatment = paste0("N:", nutrient, "\nL:", light)) %>%
  ggplot(aes(x = treatment, y = cross_combo, fill = mean_cov)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = percent(mean_cov, accuracy = 0.1)), size = 3) +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = median(df$kelp_coverage, na.rm = TRUE),
    name     = "Coverage",
    labels   = percent
  ) +
  facet_wrap(~ date, ncol = 3) +
  labs(
    title    = "Kelp coverage by cross combination and treatment",
    subtitle = "Faceted by timepoint",
    x        = "Treatment",
    y        = "Cross combination"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text  = element_text(size = 10)
  )

ggsave("results/coverage_heatmap_all_timepoints.png", width = 15, height = 6, dpi = 150)
cat("Plot saved: results/coverage_heatmap_all_timepoints.png\n")


# =============================================================================
# Kelp-only reflectance analysis
# =============================================================================

# Prepare kelp reflectance data
kelp_ref <- raw %>%
  filter(channel == "green") %>%
  mutate(
    cross_id  = as.integer(str_remove(sample_id, "C")),
    date      = factor(date_readable,
                       levels = c("10.03.2026", "17.03.2026", "27.03.2026"),
                       ordered = TRUE)
  ) %>%
  left_join(master, by = "cross_id") %>%
  filter(!is.na(date), !is.na(nutrient), !is.na(light),
         !is.na(kelp_mean_reflectance), kelp_mean_reflectance != "") %>%
  mutate(
    kelp_mean_reflectance = as.numeric(kelp_mean_reflectance),
    treatment = paste0("N:", nutrient, " L:", light)
  )

# --- 1. Trajectory plot ---
kelp_ref %>%
  group_by(date, nutrient, light, treatment) %>%
  summarise(
    mean_ref = mean(kelp_mean_reflectance, na.rm = TRUE),
    se       = sd(kelp_mean_reflectance, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  ggplot(aes(x = date, y = mean_ref, colour = treatment, group = treatment)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_ref - se, ymax = mean_ref + se),
                width = 0.2, linewidth = 1) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = "Kelp-only mean reflectance over time by treatment",
    subtitle = "Only pixels classified as kelp (R:G > 1.05)",
    x        = "Date",
    y        = "Mean reflectance (%)",
    colour   = "Treatment"
  ) +
  theme_bw(base_size = 13)

ggsave("results/kelp_reflectance_trajectory.png", width = 10, height = 6, dpi = 150)
cat("Plot saved: results/kelp_reflectance_trajectory.png\n")

# --- 2. Heatmap at most recent timepoint ---
kelp_ref %>%
  filter(!is.na(cross_combo)) %>%
  group_by(cross_id, cross_combo, nutrient, light, date) %>%
  summarise(kelp_mean_reflectance = mean(kelp_mean_reflectance, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(cross_id, cross_combo, nutrient, light) %>%
  filter(date == max(date)) %>%
  ungroup() %>%
  group_by(cross_combo, nutrient, light) %>%
  summarise(mean_ref = mean(kelp_mean_reflectance, na.rm = TRUE), .groups = "drop") %>%
  mutate(treatment = paste0("N:", nutrient, "\nL:", light)) %>%
  ggplot(aes(x = treatment, y = cross_combo, fill = mean_ref)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = round(mean_ref, 1)), size = 4) +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = 50,
    name     = "Reflectance (%)"
  ) +
  labs(
    title    = "Kelp-only reflectance at most recent timepoint",
    subtitle = "By cross combination and treatment",
    x        = "Treatment",
    y        = "Cross combination"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank())

ggsave("results/kelp_reflectance_heatmap.png", width = 9, height = 5, dpi = 150)
cat("Plot saved: results/kelp_reflectance_heatmap.png\n")

# --- 3. Coverage vs reflectance scatter ---
kelp_ref %>%
  group_by(cross_id, date, nutrient, light, treatment, cross_combo) %>%
  summarise(
    coverage  = mean(kelp_coverage, na.rm = TRUE),
    ref       = mean(kelp_mean_reflectance, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  ggplot(aes(x = coverage, y = ref, colour = treatment)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(labels = percent) +
  facet_wrap(~ date, ncol = 3) +
  labs(
    title    = "Kelp coverage vs kelp-only reflectance",
    subtitle = "Each point = one cross; faceted by timepoint",
    x        = "Kelp coverage (% of ROI)",
    y        = "Kelp-only mean reflectance (%)",
    colour   = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("results/coverage_vs_reflectance.png", width = 12, height = 5, dpi = 150)
cat("Plot saved: results/coverage_vs_reflectance.png\n")