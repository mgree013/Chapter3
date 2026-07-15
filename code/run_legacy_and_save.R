# Wrapper to run legacy Lake.analysis.R and save env_div
setwd("/Users/matthewgreen/Library/CloudStorage/Dropbox/Manuscipts/Chapter 3/Chapter3/")
# Run the legacy analysis script
source("code/Lake.analysis.R", echo = FALSE)
# Save env_div if it exists
if (exists("env_div")) {
  dir.create("data", showWarnings = FALSE)
  save(env_div, file = "data/legacy_env_div.RData")
  cat("legacy env_div saved\n")
} else {
  cat("env_div not found after sourcing Lake.analysis.R\n")
}
