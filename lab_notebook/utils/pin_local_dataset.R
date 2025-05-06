### utils/pin_local_dataset.R`

source("../lab_notebook_hub/utils/git_helpers.R")

pin_local_dataset <- function(filepath) {
  sha <- get_git_commit(filepath)
  date <- Sys.Date()
  dest <- sprintf("data/%s_%s", date, basename(filepath))
  file.copy(filepath, dest)
  
  # Log it
  log_entry <- tibble::tibble(
    date = date,
    original = filepath,
    pinned = dest,
    sha = sha
  )
  log_path <- "../lab_notebook_hub/log/data_reference_log.csv"
  if (file.exists(log_path)) {
    log <- read_csv(log_path)
    log <- bind_rows(log, log_entry)
  } else {
    log <- log_entry
  }
  write_csv(log, log_path)
  return(dest)
}
