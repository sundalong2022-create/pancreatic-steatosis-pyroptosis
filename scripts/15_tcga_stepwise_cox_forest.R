#!/usr/bin/env Rscript

# TCGA-PAAD stepwise Cox models and a PSPS-focused forest plot.
# The plot tracks the PSPS estimate as age, stage and purity are added;
# the accompanying CSV retains all terms from every model for auditability.
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(broom)
})

options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork = TRUE)
indir <- file.path(root, 'results/tcga_paad_external')
outdir <- file.path(indir, 'figures')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read_csv(file.path(indir, 'survival_input.csv'), show_col_types = FALSE) %>%
  mutate(
    time_days = as.numeric(time_days),
    event = as.integer(event),
    PSPS_core = as.numeric(PSPS_core),
    age_years = as.numeric(age_years),
    purity = as.numeric(purity),
    stage = as.character(stage)
  )

# Keep stage harmonisation identical to the primary TCGA validation script.
stage_clean <- ifelse(
  grepl('^Stage I$|^Stage IA$|^Stage IB$', dat$stage, ignore.case = TRUE), 'I',
  ifelse(grepl('^Stage III', dat$stage, ignore.case = TRUE), 'III',
  ifelse(grepl('^Stage II', dat$stage, ignore.case = TRUE), 'II',
  ifelse(grepl('^Stage IV', dat$stage, ignore.case = TRUE), 'IV', NA))))
dat$stage_group <- factor(stage_clean, levels = c('I', 'II', 'III', 'IV'))

model_specs <- list(
  Model_1_univariable = c('PSPS_core'),
  Model_2_plus_age = c('PSPS_core', 'age_years'),
  Model_3_plus_stage = c('PSPS_core', 'age_years', 'stage_group'),
  Model_4_plus_purity = c('PSPS_core', 'age_years', 'stage_group', 'purity')
)

fit_one <- function(model_name, vars) {
  keep <- complete.cases(dat[, c('time_days', 'event', vars), drop = FALSE])
  d <- dat[keep, , drop = FALSE]
  if (nrow(d) < 50L || sum(d$event, na.rm = TRUE) < 10L) {
    stop('Insufficient complete cases for ', model_name)
  }
  d$stage_group <- droplevels(d$stage_group)
  if ('stage_group' %in% vars && nlevels(d$stage_group) < 2L) {
    stop('Insufficient stage variation for ', model_name)
  }
  rhs <- paste(vars, collapse = ' + ')
  fit <- coxph(as.formula(paste0('Surv(time_days, event) ~ ', rhs)), data = d)
  tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    transmute(
      model = model_name,
      model_variables = paste(vars, collapse = ' + '),
      n_complete = nrow(d),
      events = sum(d$event, na.rm = TRUE),
      term,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p_value = p.value
    )
  list(fit = fit, data = d, terms = tidy_fit)
}

fits <- lapply(names(model_specs), function(nm) fit_one(nm, model_specs[[nm]]))
names(fits) <- names(model_specs)
all_terms <- bind_rows(lapply(fits, `[[`, 'terms'))
write_csv(all_terms, file.path(indir, 'tcga_stepwise_cox_forest_data.csv'))

# The primary forest display is PSPS across the prespecified adjustment steps.
# Full covariate estimates remain available in tcga_stepwise_cox_forest_data.csv.
psps_plot <- all_terms %>%
  filter(term == 'PSPS_core') %>%
  mutate(
    model = factor(model, levels = names(model_specs)),
    label = c(
      'Model 1: PSPS',
      'Model 2: PSPS + age',
      'Model 3: PSPS + age + stage',
      'Model 4: PSPS + age + stage + purity'
    )[as.integer(model)],
    label = factor(label, levels = rev(c(
      'Model 1: PSPS',
      'Model 2: PSPS + age',
      'Model 3: PSPS + age + stage',
      'Model 4: PSPS + age + stage + purity'
    )))
  )

p <- ggplot(psps_plot, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = 'grey45') +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = label), linewidth = 0.65, colour = '#2C7FB8') +
  geom_point(size = 2.5, colour = '#D95F02') +
  scale_x_log10() +
  theme_classic(base_size = 11) +
  labs(
    title = 'TCGA-PAAD: stepwise Cox estimates for PSPS_core',
    x = 'Hazard ratio per 1-unit increase in PSPS_core (log scale)',
    y = NULL
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(outdir, 'TCGA_PSPS_stepwise_Cox_forest.pdf'),
  p, width = 8.5, height = 4.8, useDingbats = FALSE
)
ggsave(
  file.path(outdir, 'TCGA_PSPS_stepwise_Cox_forest.png'),
  p, width = 8.5, height = 4.8, dpi = 220
)

summary_tbl <- all_terms %>%
  distinct(model, model_variables, n_complete, events) %>%
  arrange(match(model, names(model_specs)))
write_csv(summary_tbl, file.path(indir, 'tcga_stepwise_cox_model_summary.csv'))

cat('Wrote: ', file.path(indir, 'tcga_stepwise_cox_forest_data.csv'), '\n', sep = '')
cat('Wrote: ', file.path(outdir, 'TCGA_PSPS_stepwise_Cox_forest.pdf'), '\n', sep = '')
print(psps_plot %>% select(model, label, n_complete, events, HR, CI_low, CI_high, p_value))
