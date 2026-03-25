#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(readr)
})

# -----------------------------------------------------------------------------
# Helper: Convert plate reader time into minutes since first read
#
# Handles two formats:
#   1) Excel numeric time (fraction of a day, e.g. 0.00104166...)
#   2) "hh:mm:ss" strings (e.g. "00:01:30")
#
# Returns numeric minutes relative to the first non-missing entry (so first = 0).
# -----------------------------------------------------------------------------
time_to_minutes <- function(raw_vec) {
  suppressWarnings({
    as_num <- as.numeric(raw_vec)
  })
  looks_like_excel_time <- !is.na(as_num) & as_num >= 0 & as_num < 1
  
  # Excel fraction-of-day -> minutes
  minutes_from_excel <- as_num * 24 * 60
  minutes_from_excel[!looks_like_excel_time] <- NA_real_
  
  # "hh:mm:ss" string -> minutes
  parsed_hms <- suppressWarnings(lubridate::hms(raw_vec))
  minutes_from_hms <- as.numeric(parsed_hms) / 60
  minutes_from_hms[is.na(parsed_hms)] <- NA_real_
  
  # Prefer numeric interpretation if present, else hms()
  combined_minutes <- ifelse(
    !is.na(minutes_from_excel),
    minutes_from_excel,
    minutes_from_hms
  )
  
  # Baseline shift so first valid = 0
  first_valid <- combined_minutes[which(!is.na(combined_minutes))[1]]
  combined_minutes <- combined_minutes - first_valid
  
  combined_minutes
}

# -----------------------------------------------------------------------------
# Helper: Produce consistent "hh:mm:ss" string for QC
#
# Works for both numeric Excel time and "hh:mm:ss" text input.
# -----------------------------------------------------------------------------
pretty_time_string <- function(raw_vec) {
  suppressWarnings({
    as_num <- as.numeric(raw_vec)
  })
  looks_like_excel_time <- !is.na(as_num) & as_num >= 0 & as_num < 1
  
  # Excel numeric -> seconds
  sec_from_excel <- as_num * 24 * 60 * 60
  sec_from_excel[!looks_like_excel_time] <- NA_real_
  
  # "hh:mm:ss" -> seconds
  parsed_hms <- suppressWarnings(lubridate::hms(raw_vec))
  sec_from_hms <- as.numeric(parsed_hms)
  sec_from_hms[is.na(parsed_hms)] <- NA_real_
  
  # choose whichever we have
  seconds_total <- ifelse(
    !is.na(sec_from_excel),
    sec_from_excel,
    sec_from_hms
  )
  
  # format to HH:MM:SS
  formatted <- ifelse(
    is.na(seconds_total),
    NA_character_,
    sprintf(
      "%02d:%02d:%02d",
      floor(seconds_total / 3600),
      floor((seconds_total %% 3600) / 60),
      floor(seconds_total %% 60)
    )
  )
  
  formatted
}

# -----------------------------------------------------------------------------
# clean_block()
#
# Takes a rectangular block (one channel: donor or acceptor).
# The first row of the block is the header row:
#   Time, Temperature, A1, A2, ..., H12
# We drop Temperature, reshape wells, add Cycle, and compute time.
#
# Returns long tibble with:
#   Cycle
#   Well
#   Time_<channel>_raw        (hh:mm:ss as string)
#   Time_<channel>_min        (numeric minutes since first read, rounded later)
#   Intensity_<channel>       (numeric)
# -----------------------------------------------------------------------------
clean_block <- function(raw_block, channel_prefix) {
  
  header_row <- raw_block %>% slice(1)
  data_rows  <- raw_block %>% slice(-1)
  
  # Promote first row to actual column names
  colnames(data_rows) <- as.character(unlist(header_row))
  
  # Drop Temperature or similar columns we don't need
  keep_cols <- colnames(data_rows)[
    !str_detect(colnames(data_rows), regex("^temp", ignore_case = TRUE))
  ]
  data_rows <- data_rows[, keep_cols, drop = FALSE]
  
  # Identify the time column (first one containing "time")
  time_col_name <- colnames(data_rows)[
    str_detect(colnames(data_rows), regex("time", ignore_case = TRUE))
  ][1]
  
  if (is.na(time_col_name)) {
    stop(paste0("No time column found in ", channel_prefix, " block."))
  }
  
  # Prepare standardized column names
  time_raw_col <- paste0("Time_", channel_prefix, "_raw")
  time_min_col <- paste0("Time_", channel_prefix, "_min")
  intensity_col <- paste0("Intensity_", channel_prefix)
  
  # Rename the time column to the standardized raw name
  data_rows <- data_rows %>%
    rename(!!time_raw_col := all_of(time_col_name))
  
  # Add Cycle index (1,2,3,... acquisition order)
  data_rows <- data_rows %>%
    mutate(Cycle = dplyr::row_number())
  
  # Compute nice hh:mm:ss (raw QC view),
  # and elapsed minutes from first read (numeric)
  data_rows <- data_rows %>%
    mutate(
      !!time_min_col := time_to_minutes(.data[[time_raw_col]]),
      !!time_raw_col := pretty_time_string(.data[[time_raw_col]])
    )
  
  # Figure out which columns correspond to wells
  meta_cols <- c("Cycle", time_raw_col, time_min_col)
  well_cols <- setdiff(colnames(data_rows), meta_cols)
  
  # wide -> long per well
  long_df <- data_rows %>%
    pivot_longer(
      cols = all_of(well_cols),
      names_to = "Well",
      values_to = intensity_col
    ) %>%
    mutate(
      # normalize well names, e.g. A01 -> A1
      Well = str_replace(Well, "^([A-Ha-h])0*([0-9]+)$", "\\1\\2"),
      
      # convert intensity to numeric (suppressing warnings for blanks etc.)
      !!intensity_col := suppressWarnings(
        as.numeric(.data[[intensity_col]])
      )
    )
  
  long_df
}

# -----------------------------------------------------------------------------
# process_plate_file()
#
# Arguments:
#   input            : Excel filename (string, must be in working dir)
#   sheet            : sheet index (1) or name ("Sheet1")
#   donor_start/end  : row range (in Excel-like 1-based rows) of donor block
#   acceptor_start/end: row range of acceptor block
#   output           : CSV filename to write
#
# What it does:
#   - Reads Excel
#   - Extracts donor + acceptor blocks
#   - Cleans them
#   - Merges by (Cycle, Well)
#   - Drops wells that are completely unused
#   - Rounds the time columns to avoid ugly floating-point tails
#   - Writes a tidy CSV with:
#       Cycle
#       Well
#       Time_donor_raw
#       Time_donor_min
#       Time_acceptor_raw
#       Time_acceptor_min
#       Intensity_donor
#       Intensity_acceptor
#
#   (No ratio column. Normalization is downstream.)
# -----------------------------------------------------------------------------
process_plate_file <- function(input,
                               sheet,
                               donor_start, donor_end,
                               acceptor_start, acceptor_end,
                               output) {
  
  # Basic safety
  if (!is.character(input)) {
    stop("'input' must be a string, e.g. \"myfile.xlsx\".")
  }
  if (!file.exists(input)) {
    stop(paste0(
      "File not found: ", input,
      "\nWorking directory is: ", getwd(),
      "\nFiles here: ", paste(list.files(), collapse = ", ")
    ))
  }
  
  # Read the full sheet with no assumed headers
  raw_all <- readxl::read_excel(
    path = input,
    sheet = sheet,
    col_names = FALSE
  ) %>%
    mutate(across(everything(), as.character))
  
  n_rows <- nrow(raw_all)
  if (donor_start < 1 || donor_end > n_rows ||
      acceptor_start < 1 || acceptor_end > n_rows) {
    stop("Specified row ranges exceed number of rows in the sheet.")
  }
  
  # Slice the donor and acceptor regions
  donor_block <- raw_all %>% slice(donor_start:donor_end)
  acceptor_block <- raw_all %>% slice(acceptor_start:acceptor_end)
  
  # Clean blocks individually
  donor_long <- clean_block(donor_block, "donor")
  acceptor_long <- clean_block(acceptor_block, "acceptor")
  
  # Merge by Cycle + Well
  merged <- full_join(donor_long, acceptor_long, by = c("Cycle", "Well")) %>%
    mutate(
      any_donor    = !is.na(Intensity_donor),
      any_acceptor = !is.na(Intensity_acceptor)
    )
  
  # Identify wells that are actually used at all
  used_wells <- merged %>%
    group_by(Well) %>%
    summarise(
      used = any(any_donor | any_acceptor, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(used) %>%
    pull(Well)
  
  # Final tidy table
  final_df <- merged %>%
    filter(Well %in% used_wells) %>%
    select(
      Cycle,
      Well,
      Time_donor_raw,
      Time_donor_min,
      Time_acceptor_raw,
      Time_acceptor_min,
      Intensity_donor,
      Intensity_acceptor
    ) %>%
    arrange(Well, Cycle) %>%
    # round the minute columns to 3 decimal places for nicer output
    mutate(
      Time_donor_min      = round(Time_donor_min, 3),
      Time_acceptor_min   = round(Time_acceptor_min, 3)
    )
  
  # Write CSV
  readr::write_csv(final_df, output)
  
  invisible(final_df)
}

# -----------------------------------------------------------------------------
# Command-line interface for batch usage
# -----------------------------------------------------------------------------
if (!interactive()) {
  option_list <- list(
    make_option(c("--input"), type = "character"),
    make_option(c("--sheet"), type = "character", default = "1"),
    make_option(c("--donor_start"), type = "integer"),
    make_option(c("--donor_end"), type = "integer"),
    make_option(c("--acceptor_start"), type = "integer"),
    make_option(c("--acceptor_end"), type = "integer"),
    make_option(c("--output"), type = "character")
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  
  # Convert sheet to numeric if it looks numeric
  sheet_val <- suppressWarnings(as.numeric(opt$sheet))
  if (!is.na(sheet_val)) {
    opt$sheet <- sheet_val
  }
  
  do.call(
    process_plate_file,
    as.list(opt)
  )
}
