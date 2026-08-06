# Survival Analysis in R: Kaplan–Meier and Cox Regression 

R workflow for time-to-event analysis: Kaplan–Meier curves, Cox proportional hazards regression, assumption testing, and summary tables.

The example workflow uses the built-in **NCCTG Lung Cancer** dataset (`survival::lung`; Loprinzi et al. 1994, *J Clin Oncol* 12:601–607), so the
repository can be run immediately without downloading external data.
![Kaplan–Meier Curve](output/km_curve.png)

---

## Features

- Kaplan–Meier survival curves with confidence intervals and risk table
- Log-rank test
- Baseline characteristics table (`gtsummary::tbl_summary`)
- Cox proportional hazards regression
- Proportional hazards assumption testing using Schoenfeld residuals (`cox.zph`)
- Unadjusted and adjusted hazard ratio tables, side by side
- Tables exported to Word with **gtsummary** and **gt**

---

## Requirements

```r
install.packages(c(
  "survival",
  "survminer",
  "gtsummary",
  "gt",
  "dplyr",
  "tidyr",
  "labelled"
))
```

Tested with:

- R 4.4.2
- survival 3.6-6
- survminer 0.5.2
- gtsummary 2.5.1
- gt 1.3.0
- dplyr 1.2.0
- tidyr 1.3.2

> **Note:** `gtsummary` 2.x is required. Some functions, such as `modify_column_merge()` and `modify_header()`, behave differently in older releases.

---

## Using your own data

Replace the preprocessing section with your own clinical dataset. The script expects:

- a follow-up time variable
- an event indicator (`0 = censored`, `1 = event`)
- one or more clinical covariates

Everything downstream can remain unchanged.

Three practical points built into the script:

- **Handle missing values up front.** `coxph()` drops rows with `NA` silently, so unadjusted and adjusted models would otherwise be fitted on different samples and would not be comparable.
- **Aim for at least 10 events per model parameter.** This example has 150 events and 5 terms.
- **Check the proportional hazards assumption** before interpreting hazard ratios. It is the assumption the entire model rests on, and it is routinely skipped.

---

## Example output

In the lung data, age is significant unadjusted (HR 1.02, p = 0.027) but not after adjustment (p = 0.177) — a clean case of confounding, and the reason the side-by-side table is worth producing.

The proportional hazards assumption holds here (global p = 0.393). Had it been violated, the options would be stratifying on the offending variable, adding a time-varying coefficient, or moving to an accelerated failure time model.

---

## Disclaimer

This repository is intended as a learning resource and starting point for survival analysis in R.

The example uses the built-in **NCCTG Lung Cancer** dataset (`survival::lung`) only to demonstrate the analysis workflow. The example results should not be interpreted as new scientific findings or clinical recommendations.

Before applying this workflow to your own research, verify data quality, assess model assumptions, and choose variables based on your study design.

---

## License

MIT
