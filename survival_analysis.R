# =============================================================================
# Survival Analysis: Kaplan-Meier and Cox Regression
#
# Author: Nikhil Kirtipal
#
# Demonstrates:
#   - Kaplan-Meier estimation and log-rank test
#   - Cox proportional hazards regression
#   - PH assumption testing (Schoenfeld residuals)
#   - Publication-ready summary tables
#
# Data: survival::lung (NCCTG Lung Cancer, 228 patients) - built into the package
#       Loprinzi et al. (1994), J Clin Oncol 12:601-607
#
# Tested: R 4.4.2, survival 3.6-6, survminer 0.5.2, gtsummary 2.5.1,
#         gt 1.3.0, dplyr 1.2.0, tidyr 1.3.2
#
# Replace the data preparation section to use your own clinical dataset.
# =============================================================================

library(survival)
library(survminer)
library(gtsummary)
library(gt)
library(dplyr)
library(tidyr)

data("lung", package = "survival")


# ---- 0. Data preparation ----------------------------------------------------
# lung codes status as 1 = censored, 2 = dead. Survival functions expect
# 0 = censored, 1 = event, so subtract 1.
#
# ph.ecog has a single patient at level 3, which gives an unstable univariable
# model, so it is dropped. Complete cases are used throughout so that the
# unadjusted and adjusted models are fitted on the same rows.

dat <- lung %>%
  filter(!is.na(ph.ecog), ph.ecog < 3) %>%
  mutate(
    status  = status - 1,
    sex     = factor(sex, levels = 1:2, labels = c("Male", "Female")),
    ph.ecog = factor(ph.ecog, levels = 0:2,
                     labels = c("Asymptomatic", "Ambulatory", "In bed <50%"))
  ) %>%
  drop_na(time, status, sex, age, ph.ecog, wt.loss)

cat("n =", nrow(dat), " events =", sum(dat$status), "\n")

# Variable labels feed straight through to the gtsummary tables
dat <- dat %>%
  labelled::set_variable_labels(
    time     = "Follow-up time (days)",
    age      = "Age (years)",
    sex      = "Sex",
    ph.ecog  = "ECOG performance status",
    ph.karno = "Karnofsky score (physician)",
    wt.loss  = "Weight loss (lbs, 6 months)"
  )


# ---- 1. Kaplan-Meier --------------------------------------------------------

surv_obj <- Surv(time = dat$time, event = dat$status)
fit <- survfit(surv_obj ~ sex, data = dat)

print(fit)
print(surv_median(fit))

km_plot <- ggsurvplot(
  fit,
  data         = dat,
  pval         = TRUE,          # log-rank test
  conf.int     = TRUE,
  risk.table   = TRUE,
  xlab         = "Days",
  ylab         = "Survival probability",
  legend.title = "Sex",
  legend.labs  = c("Male", "Female")
)
print(km_plot)


# ---- 2. Baseline characteristics --------------------------------------------

tab1 <- dat %>%
  select(age, sex, ph.ecog, wt.loss) %>%
  tbl_summary(
    by        = sex,
    missing   = "no",
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      all_continuous()  ~ "{mean} ± {sd}"
    ),
    digits = list(
      all_categorical() ~ c(0, 1),
      all_continuous()  ~ c(1, 1)
    )
  ) %>%
  bold_labels() %>%
  add_overall(last = TRUE) %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 3))

print(tab1)


# ---- 3. Cox proportional hazards --------------------------------------------

cox <- coxph(Surv(time, status) ~ sex + age + ph.ecog + wt.loss, data = dat)
summary(cox)


# ---- 4. Check the proportional hazards assumption ---------------------------
# The Cox model assumes hazard ratios stay constant over time. This is the
# assumption the entire model rests on, and it is routinely skipped.
#
# cox.zph tests the correlation between scaled Schoenfeld residuals and time.
#   p < 0.05 means the assumption is violated for that term.
#
# If violated: stratify on the offending variable, add a time-varying
# coefficient, or move to an accelerated failure time (AFT) model.

zph <- cox.zph(cox)
print(zph)

par(mfrow = c(2, 2))
plot(zph)
par(mfrow = c(1, 1))

# Stratified alternative, if sex were to violate PH:
# cox_strat <- coxph(Surv(time, status) ~ strata(sex) + age + ph.ecog + wt.loss,
#                    data = dat)


# ---- 5. Unadjusted hazard ratios --------------------------------------------

hr_uni <- dat %>%
  select(time, status, sex, age, ph.ecog, wt.loss) %>%
  tbl_uvregression(
    method       = coxph,
    y            = Surv(time, status),
    exponentiate = TRUE,
    pvalue_fun   = ~style_pvalue(.x, digits = 3)
  ) %>%
  modify_column_merge(
    pattern = "{estimate} ({conf.low}, {conf.high})",
    rows    = !is.na(estimate)
  ) %>%
  modify_header(estimate ~ "**HR (95% CI)**") %>%
  bold_labels()

print(hr_uni)


# ---- 6. Adjusted hazard ratios ----------------------------------------------

hr_multi <- cox %>%
  tbl_regression(
    exponentiate = TRUE,
    pvalue_fun   = ~style_pvalue(.x, digits = 3)
  ) %>%
  modify_column_merge(
    pattern = "{estimate} ({conf.low}, {conf.high})",
    rows    = !is.na(estimate)
  ) %>%
  modify_header(estimate ~ "**HR (95% CI)**") %>%
  bold_labels()

print(hr_multi)


# ---- 7. Side-by-side table --------------------------------------------------

hr_table <- tbl_merge(
  tbls        = list(hr_uni, hr_multi),
  tab_spanner = c("**Unadjusted**", "**Adjusted**")
)

print(hr_table)


# ---- 8. Export --------------------------------------------------------------
# ggsave() does not work on a ggsurvplot object (it is a list, not a ggplot),
# so the plot is written through a png device instead.

dir.create("output", showWarnings = FALSE)

tab1     %>% as_gt() %>% gtsave("output/baseline_table.docx")
hr_table %>% as_gt() %>% gtsave("output/hazard_ratios.docx")

png("output/km_curve.png", width = 2400, height = 2100, res = 300)
print(km_plot)
dev.off()


# =============================================================================
# Notes
#
# Events per variable: aim for at least 10 events per model term. This example
#   has 150 events and 5 terms, which is comfortable. Fitting 9 covariates to
#   50 events is not.
#
# Missing data: coxph drops rows with NA silently. Handled up front here so
#   every model uses the same sample. Otherwise the unadjusted and adjusted
#   estimates come from different n and are not comparable.
#
# Multiple testing: p-values here are unadjusted. Correct them if screening
#   many covariates.
# =============================================================================
