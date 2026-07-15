# Compare legacy env_div to new env_div
setwd("/Users/matthewgreen/Library/CloudStorage/Dropbox/Manuscipts/Chapter 3/Chapter3/")
legacy <- new.env()
newenv <- new.env()
if (file.exists("data/legacy_env_div.RData")) load("data/legacy_env_div.RData", envir = legacy)
if (file.exists("data/new_lake_analytics.RData")) load("data/new_lake_analytics.RData", envir = newenv)
# Attempt to find env_div objects
legacy_env <- if (exists("env_div", envir = legacy)) legacy$env_div else NULL
new_env <- if (exists("env_div", envir = newenv)) newenv$env_div else NULL
if (is.null(legacy_env)) { cat("legacy env_div not found\n") }
if (is.null(new_env)) { cat("new env_div not found\n") }
if (!is.null(legacy_env) && !is.null(new_env)) {
  library(dplyr)
  legacy_env <- legacy_env %>% mutate(lake_id = as.integer(lake_id), survey_date = as.character(survey_date))
  new_env <- new_env %>% mutate(lake_id = as.integer(lake_id), survey_date = as.character(survey_date))
  # Join on lake_id and survey_date
  cmp <- inner_join(legacy_env, new_env, by = c("lake_id","survey_date"), suffix = c("_legacy","_new"))
  if (nrow(cmp)==0) { cat("No matching rows between legacy and new env_div\n") }
  # Compute diffs for N1, betas.LCBD, Com.Size
  vars <- c("N1","betas.LCBD","Com.Size")
  diffs <- lapply(vars, function(v){
    a <- cmp[[paste0(v,"_legacy")]]
    b <- cmp[[paste0(v,"_new")]]
    data.frame(var=v,
               n_match = length(a),
               mean_legacy = mean(a, na.rm=TRUE),
               mean_new = mean(b, na.rm=TRUE),
               mean_diff = mean(b - a, na.rm=TRUE),
               median_diff = median(b - a, na.rm=TRUE),
               mean_abs_diff = mean(abs(b - a), na.rm=TRUE),
               ratio_mean = mean(b / (a + 1e-12), na.rm=TRUE))
  })
  diffs_df <- do.call(rbind, diffs)
  write.csv(diffs_df, file = "data/env_div_comparison_summary.csv", row.names = FALSE)
  cat("Comparison summary written to data/env_div_comparison_summary.csv\n")
  # Also write full cmp for inspection
  write.csv(cmp, file = "data/env_div_full_comparison.csv", row.names = FALSE)
  cat("Full comparison written to data/env_div_full_comparison.csv\n")
  print(diffs_df)
}
