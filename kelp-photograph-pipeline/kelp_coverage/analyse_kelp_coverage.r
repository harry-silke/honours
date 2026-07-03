#!/usr/bin/env Rscript
# Visualise kelp master_results.csv and rank individuals by growth rate.
#
# Input schema: master_results.csv from run_imagej.py. One row per pair per
# date, with left_/right_ columns per individual (left = id1, right = id2).
# If you use the measure_kelp_coverage.py per-folder output instead, change
# the column names in side_long() (there they are coverage_left,
# area_cm2_left, refl_R_left, ... rather than left_kelp_area_cm2 etc.).
#
# "Growth rate" here is RGR, the slope of ln(area) over day, the size-
# normalised standard in growth biology. Absolute slope of area over day is
# also reported. Treatment comes from the id suffix: none = 14C, B = 17C,
# C = 20C.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(lubridate); library(stringr)
})

## ---- config ----
in_csv  <- "master_results.csv"
out_dir <- "figures"
top_n   <- 15
dir.create(out_dir, showWarnings = FALSE)

## ---- load ----
df <- read_csv(in_csv, show_col_types = FALSE)

## ---- wide (per pair) -> long (per individual) ----
side_long <- function(d, side, idcol) {
  p <- paste0(side, "_")
  tibble(
    pair_key   = d$pair_key,
    date_raw   = as.character(d$date),
    individual = as.character(d[[idcol]]),
    area       = d[[paste0(p, "kelp_area_cm2")]],
    coverage   = d[[paste0(p, "kelp_coverage")]],
    refl_red   = d[[paste0(p, "kelp_mean_reflectance_red")]],
    refl_green = d[[paste0(p, "kelp_mean_reflectance_green")]],
    refl_blue  = d[[paste0(p, "kelp_mean_reflectance_blue")]]
  )
}
long <- bind_rows(side_long(df, "left", "id1"),
                  side_long(df, "right", "id2")) %>%
  filter(!is.na(individual), individual != "1000")

## ---- derive treatment, cross id, time ----
long <- long %>%
  mutate(
    suffix      = str_extract(individual, "[A-Za-z]$"),
    suffix      = ifelse(is.na(suffix), "", suffix),
    treatment   = recode(suffix, "" = "14C", "B" = "17C", "C" = "20C",
                         .default = paste0("?", suffix)),
    cross_id    = as.integer(str_replace_all(individual, "[^0-9]", "")),
    date        = dmy(date_raw),
    refl_bright = rowMeans(cbind(refl_red, refl_green, refl_blue), na.rm = TRUE)
  )
if (all(is.na(long$date)))
  stop("Could not parse 'date' as DD.MM.YYYY; check the date format.")
long <- long %>% mutate(day = as.numeric(date - min(date, na.rm = TRUE)))

## ---- per-individual growth ----
growth <- long %>%
  filter(!is.na(area), area > 0) %>%
  group_by(individual, treatment, cross_id) %>%
  summarise(
    n_time     = n(),
    day_span   = max(day) - min(day),
    area_first = area[which.min(day)],
    area_last  = area[which.max(day)],
    abs_slope  = if (n() >= 2 && max(day) > min(day))
                   unname(coef(lm(area ~ day))[2]) else NA_real_,
    rgr        = if (n() >= 2 && max(day) > min(day))
                   unname(coef(lm(log(area) ~ day))[2]) else NA_real_,
    .groups = "drop"
  ) %>%
  filter(n_time >= 2, day_span > 0, !is.na(rgr)) %>%
  arrange(desc(rgr))

top <- head(growth, top_n)
write_csv(growth, file.path(out_dir, "growth_all.csv"))
write_csv(top,    file.path(out_dir, "top_growth.csv"))
cat(sprintf("Top %d individuals by RGR (per day):\n", top_n))
print(as.data.frame(top[c("individual", "treatment", "n_time",
                          "area_first", "area_last", "rgr", "abs_slope")]))

## ---- plots ----
theme_set(theme_minimal(base_size = 12))

# 1. trajectories of the top individuals
traj <- long %>% filter(individual %in% top$individual, !is.na(area), area > 0)
p1 <- ggplot(traj, aes(day, area, group = individual, colour = treatment)) +
  geom_line() + geom_point(size = 1) +
  labs(title = sprintf("Growth trajectories: top %d by RGR", top_n),
       x = "day", y = "kelp area (cm2)")
ggsave(file.path(out_dir, "top_trajectories.png"), p1,
       width = 8, height = 5, dpi = 140)

# 2. RGR distribution by treatment
p2 <- ggplot(growth, aes(treatment, rgr, fill = treatment)) +
  geom_violin(alpha = .4, colour = NA) +
  geom_jitter(width = .12, size = .8, alpha = .5) +
  labs(title = "Relative growth rate by treatment",
       x = NULL, y = "RGR (per day)") +
  guides(fill = "none")
ggsave(file.path(out_dir, "rgr_by_treatment.png"), p2,
       width = 6, height = 4.5, dpi = 140)

# helper: mean +/- SE over time by treatment
summ <- function(d, col) {
  d %>% filter(!is.na(.data[[col]])) %>%
    group_by(treatment, day) %>%
    summarise(mean = mean(.data[[col]]),
              se = sd(.data[[col]]) / sqrt(n()), .groups = "drop")
}
ribbon <- function(s, ylab, title) {
  ggplot(s, aes(day, mean, colour = treatment, fill = treatment)) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
                alpha = .2, colour = NA) +
    geom_line() + geom_point() +
    labs(title = title, x = "day", y = ylab)
}

# 3-5. area, coverage, reflectance over time by treatment
ggsave(file.path(out_dir, "area_over_time.png"),
       ribbon(summ(long, "area"), "area (cm2)",
              "Mean kelp area over time by treatment"),
       width = 7.5, height = 5, dpi = 140)
ggsave(file.path(out_dir, "coverage_over_time.png"),
       ribbon(summ(long, "coverage"), "coverage",
              "Mean coverage over time by treatment"),
       width = 7.5, height = 5, dpi = 140)
ggsave(file.path(out_dir, "reflectance_over_time.png"),
       ribbon(summ(long, "refl_bright"), "reflectance (calibrated)",
              "Mean kelp reflectance over time by treatment"),
       width = 7.5, height = 5, dpi = 140)

cat(sprintf("\nWrote tables and figures to %s/\n", out_dir))