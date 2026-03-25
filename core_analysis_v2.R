############################################################
## core_analysis_v2.R
##
## Core pipeline for kinetic plate-reader FRET protease assays.
## The file is a revised version "core_analysis.R", originally by Marius Krogsgaard Thomsen.
## Revised in March 2026 to account for assays with multiple substrate variants
## 
##
## Features:
##  - Step 1+2: Read + merge tidy kinetic data with plate map
##  - Step 3: Bleach/drift correction per variant using substrate controls
##  - Step 4: Reference dynamic range:
##        R_min_global = intact substrate signal
##        R_max_global = plateaued WT-like wells you trust
##  - Step 5: Fraction_cleaved
##  - Step 6: Initial rate fitting using adaptive "fraction_ceiling"
##  - Step 7: QC workbook export (xlsx)
##
## ── Key changes from core_analysis.R ──────────────────
## CHANGE 1: Variant is now substrate identity
## - The variable "variant" is assumed to be substrate identity (SubG/SubV/SubP).
## - The variable "wt_variant" is replaced by "most_cleaved_condition"
##
## CHANGE 2: Bleach-correction is per substrate variant (step 3)
## - Finds a delta_D and delta_A for each substrate variant.
##
## CHANGE 3: Trusted wells are determined by the variable "most_cleaved_condition" (step 4)
## - The variable "wt_variant" has been replaced.
## - Trusted wells used in the calculation of R_max are therefore based on "most_cleaved_condition"
##
## CHANGE 4: R_max and R_min are are per variant (step 4)
## - Finds an R_max and R_min for each substrate variant
##
## CHANGE 5: Collapsed dynamic range detection (step 4)
## - Detects substrate variants with no activity
## - Computes R_ctrl_plateau = mean R of the same substrate's control wells over the same cycle window. 
## - Declares "collapsed" when (R_max - R_ctrl_plateau) < min_signal_above_ctrl_frac * R_min
## 
## CHANGE 6: "wt_variant" is replaced by "most_cleaved_condition" in wrapper function (step 7)
## - Due to CHANGE 1
##
## CHANGE 7: Excluding NA fraction_cleaved data during regression (step 6)
## - Filters out collapsed-substrate wells.
## - Collapsed substrate variants return NA rate during regression rather when a near-zero slope.
##
## ───────────────────────────────────────────────────────────
##
## Main entry point:
##
##   run_core_pipeline_from_tidy(
##     tidy_csv_file,
##     plate_map_file,
##     S0_uM                      = 5,
##     fit_window_min             = 5,
##     fraction_ceiling           = 0.50,
##     most_cleaved_condition     = "1to1",   ## CHANGE 3/9: replaces wt_variant
##     wt_wells_for_Rmax          = NULL,
##     plateau_points             = 3,
##     plateau_span_cycles        = 10,
##     plateau_tolerance          = 0.01,
##     min_signal_above_ctrl_frac = 0.10,     ## CHANGE 5: new parameter
##     qc_output_file             = "core_output_QC.xlsx"
##   )
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(writexl)
  library(stats)
})

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

############################################################
## Step 1+2.  Read and merge kinetic data with plate map.
############################################################
read_and_merge_inputs <- function(tidy_csv_file,
                                  plate_map_file) {
  
  # Read csv with cleaned kinetic data
  kinetic_raw <- readr::read_csv(
    file = tidy_csv_file,
    col_types = cols(
      Cycle               = col_integer(),
      Well                = col_character(),
      Time_donor_raw      = col_character(),
      Time_donor_min      = col_double(),
      Time_acceptor_raw   = col_character(),
      Time_acceptor_min   = col_double(),
      Intensity_donor     = col_double(),
      Intensity_acceptor  = col_double()
    )
  )
  
  # Read csv with plate map
  plate_map <- readr::read_csv(
    file = plate_map_file,
    col_types = cols(
      Well      = col_character(),
      Variant   = col_character(),
      Condition = col_character(),
      Type      = col_character()
    )
  )
  
  # Merging the csv of the plate map into the intensity data
  merged_raw <- kinetic_raw %>%
    left_join(plate_map, by = "Well") %>%
    select(
      Well, Variant, Condition, Type,
      Cycle,
      Time_donor_raw, Time_donor_min,
      Time_acceptor_raw, Time_acceptor_min,
      Intensity_donor, Intensity_acceptor
    )
  
  list(
    plate_map   = plate_map,
    merged_raw  = merged_raw
  )
}


############################################################
## Step 3.  Per-Variant bleach/drift correction.
##
## - Substrate controls define per-cycle donor_ctrl, acceptor_ctrl.
## - ref_cycle = smallest Cycle > 1.
## - delta_D, delta_A = per cycle substrate controls substrated by intensity at reference cycle
## - Additive bleach correction - donor and acceptor intensities are corrected by substrating delta_D
## - Correction is specific for each substrate
##
## CHANGE 2 vs original: correction is now computed and applied
## per-Variant (substrate) rather than pooling all substrate_control
## wells together.
##
## Original approach:
##   ctrl_summary <- merged_raw %>%
##     filter(Type == "substrate_control") %>%
##     group_by(Cycle) %>%                     # no Variant grouping
##     summarise(donor_ctrl = mean(...), ...)
##   # single delta_D / delta_A applied to ALL sample wells
##
## v2 approach:
##   group_by(Variant, Cycle) so each substrate gets drift deltas
##   derived only from its own control wells.  The join to corrected
##   is by c("Variant", "Cycle") instead of just "Cycle".
##
############################################################
apply_bleach_correction <- function(merged_raw) {
  
  if (!any(merged_raw$Type == "substrate_control", na.rm = TRUE)) {
    stop("No wells with Type == 'substrate_control'. Cannot perform bleach/drift correction.")
  }
  
  # Defining ref_cycle - reference cycle > 1. First cycle is used for warm up of the reader
  ref_cycle <- merged_raw %>%
    filter(Type == "substrate_control", Cycle > 1) %>%
    summarise(rc = min(Cycle, na.rm = TRUE)) %>%
    pull(rc)
  
  if (length(ref_cycle) == 0 || is.na(ref_cycle) || is.infinite(ref_cycle)) {
    stop("Could not determine ref_cycle. Check substrate_control wells have Cycle > 1.")
  }
  
  ## CHANGE 2: group_by(Variant, Cycle) — original used group_by(Cycle) only.
  # Finding mean of donor and acceptor intensities across replicates. For control wells.
  ctrl_per_variant <- merged_raw %>%
    filter(Type == "substrate_control") %>%
    group_by(Variant, Cycle) %>%
    summarise(
      donor_ctrl    = safe_mean(Intensity_donor),
      acceptor_ctrl = safe_mean(Intensity_acceptor),
      .groups = "drop"
    )
  
  ## CHANGE 2: reference values kept per-Variant — original had a single
  ## donor_ref / acceptor_ref scalar shared by all wells.
  ctrl_ref <- ctrl_per_variant %>%
    filter(Cycle == ref_cycle) %>%
    select(Variant, donor_ref = donor_ctrl, acceptor_ref = acceptor_ctrl)
  
  # Finding deltas by substracting mean intensities with intensity at reference cycle
  bleach_deltas_table <- ctrl_per_variant %>%
    left_join(ctrl_ref, by = "Variant") %>%
    mutate(
      delta_D = donor_ctrl    - donor_ref,
      delta_A = acceptor_ctrl - acceptor_ref
    ) %>%
    select(Variant, Cycle, donor_ctrl, acceptor_ctrl, delta_D, delta_A)
  
  ## CHANGE 2: join on c("Variant", "Cycle") — original joined on "Cycle" only,
  ## meaning every sample well received the pooled-control delta regardless
  ## of which substrate it contained.
  # Correcting donor and acceptor intensities for photobleaching by substracting delta derived from control wells
  # R_corr - ratio between corrected intensities
  corrected <- merged_raw %>%
    filter(Type == "sample") %>%
    left_join(
      bleach_deltas_table %>% select(Variant, Cycle, delta_D, delta_A),
      by = c("Variant", "Cycle")
    ) %>%
    mutate(
      donor_corr    = Intensity_donor    - delta_D,
      acceptor_corr = Intensity_acceptor - delta_A,
      R_corr        = donor_corr / acceptor_corr
    )
  
  list(
    corrected           = corrected,
    ref_cycle           = ref_cycle,
    bleach_deltas_table = bleach_deltas_table
  )
}

############################################################
## Step 4.  Per-Variant dynamic range (R_min, R_max).
## 
## R_min:
##   - Mean R_corr of substrate_control wells at ref_cycle.
##   - One R_min for each substrate variant
## 
## R_max
## - We consider only specific trusted wells to define "fully cleaved".
##     Sources of trusted wells:
##       * If wt_wells_for_Rmax is given (e.g. c("A1","B1")), we
##         ONLY use those wells.
##       * Otherwise, we use sample wells where condition == most_cleaved_condition
##
##   - For each trusted well:
##       1. Sort by Cycle.
##       2. Scan forward: at each index i, look ahead
##          plateau_span_cycles cycles.
##          If R_corr[i] and R_corr[i + plateau_span_cycles]
##          differ by <= plateau_tolerance *relative*, we call
##          that well "plateaued" starting at cycle i.
##       3. plateau_ratio_for_that_well = mean R_corr from that
##          trigger index to the end of that well's trace.
##
##   - If ≥1 trusted well plateaued:
##       R_max_global = mean of those plateau_ratios.
##       method_used  = "plateau_detected"
##
##   - If NO trusted well plateaued:
##       Fall back to taking the last `plateau_points` cycles
##       from each trusted well, pool them, and average their
##       R_corr. method_used = "tail_points_fallback".
##
## R_ctrl_plateau
## - Finds ratio between donor and acceptor intensities for substrate control wells in that interval 
##   where R_max has been defined
## - Used to detect inactive substrates:
##      - If R_max - R_ctrl_plateau < min_signal_above_ctrl_frac * R_min substrate is inactive
##
## CHANGE 4 vs original: R_min and R_max are now computed separately
## for each Variant (substrate) instead of as single global scalars.
##
## CHANGE 3 vs original: trusted wells for R_max are selected by
##   Condition == most_cleaved_condition  (default "1to1")
## instead of
##   Variant == wt_variant               (original default "eTEVp").
## "eTEVp" does not appear in this plate's Variant column at all,
## so the original filter returned zero rows and the pipeline either
## crashed or silently fell back to using all wells mixed together.
##
## CHANGE 5 vs original: collapsed dynamic range detection.
## The original had no such check. Dividing by R_max - R_min when
## that difference is near zero produced Fraction_cleaved
## values far outside [0, 1]. v2 compares R_max to the concurrent
## control level (R_ctrl_plateau) rather than simply to R_min at
## ref_cycle, because both SubP sample and SubP control drift upward
## together due to acceptor photobleaching. Comparing to ref_cycle
## alone would overestimate the real dynamic range.
##
############################################################
compute_dynamic_range <- function(merged_raw,
                                  corrected,
                                  ref_cycle,
                                  bleach_deltas_table,
                                  most_cleaved_condition     = "1to1",
                                  wt_wells_for_Rmax          = NULL,
                                  plateau_points             = 3,
                                  plateau_span_cycles        = 10,
                                  plateau_tolerance          = 0.01,
                                  min_signal_above_ctrl_frac = 0.10) {
  
  all_variants <- sort(unique(corrected$Variant))
  
  ## CHANGE 4: R_min computed per-Variant by grouping substrate_control
  ## wells on Variant before summarising.
  ## Original computed a single R_min_global = donor_ref_mean / acceptor_ref_mean
  ## from the pooled-control summary row at ref_cycle (one scalar for all wells).
  # R_min - the mean R_corr at reference cycle
  r_min_by_variant <- merged_raw %>%
    filter(Type == "substrate_control", Cycle == ref_cycle) %>%
    mutate(R_ctrl = Intensity_donor / Intensity_acceptor) %>%
    group_by(Variant) %>%
    summarise(
      donor_ctrl_mean    = safe_mean(Intensity_donor),
      acceptor_ctrl_mean = safe_mean(Intensity_acceptor),
      R_min              = safe_mean(R_ctrl),
      .groups = "drop"
    )
  
  ## Plateau detection for one variant
  plateau_detect_one_variant <- function(var_name) {
    
    ## CHANGE 3: trusted wells for R_max selected by Condition (enzyme dilution)
    ## rather than by Variant (enzyme name). The original filtered:
    ##   filter(Variant == wt_variant)  # e.g. Variant == "eTEVp"
    ## In this plate Variant holds the substrate name, so that filter returned
    ## zero rows for every substrate. v4 filters:
    ##   filter(Variant == var_name, Condition == most_cleaved_condition)
    ## which selects the highest-enzyme-concentration wells for each substrate.
    if (!is.null(wt_wells_for_Rmax) && !is.null(wt_wells_for_Rmax[[var_name]])) {
      trusted_df <- corrected %>%
        filter(Variant == var_name, Well %in% wt_wells_for_Rmax[[var_name]])
    } else {
      trusted_df <- corrected %>%
        filter(Variant == var_name, Condition == most_cleaved_condition)
      if (nrow(trusted_df) == 0)
        trusted_df <- corrected %>% filter(Variant == var_name)
    }
    
    if (nrow(trusted_df) == 0) {
      return(list(
        R_max = NA_real_, method_used = "no_trusted_wells",
        plateau_cyc_start = NA_integer_, plateau_cyc_end = NA_integer_,
        plateau_by_well = tibble(Well=character(), Variant=character(),
                                 plateau_detected=logical(), plateau_cycle_trigger=integer(),
                                 plateau_ratio=double())
      ))
    }
    # Finding R_max - maximum ratio for intensities
    plateau_by_well <- trusted_df %>%
      group_by(Well) %>%
      arrange(Cycle, .by_group = TRUE) %>%
      group_modify(~{
        df <- .x; cycs <- df$Cycle; Rvals <- df$R_corr; n <- nrow(df)
        plat_idx <- NA_integer_
        # Scan through time to find first plateau index
        for (i in seq_len(n)) {
          j <- i + plateau_span_cycles
          if (j > n) break
          R_now <- Rvals[i]; R_fut <- Rvals[j]
          if (is.na(R_now) || is.na(R_fut)) next
          rel_chg <- if (abs(R_now) < .Machine$double.eps) abs(R_fut - R_now)
          else abs(R_fut - R_now) / abs(R_now)
          if (rel_chg <= plateau_tolerance) { plat_idx <- i; break }
        }
        if (!is.na(plat_idx)) {
          tibble(plateau_detected=TRUE,
                 plateau_cycle_trigger=cycs[plat_idx],
                 plateau_ratio=mean(Rvals[plat_idx:n], na.rm=TRUE),
                 plateau_cyc_start=cycs[plat_idx],
                 plateau_cyc_end=cycs[n])
        } else {
          tibble(plateau_detected=FALSE, plateau_cycle_trigger=NA_integer_,
                 plateau_ratio=NA_real_, plateau_cyc_start=NA_integer_,
                 plateau_cyc_end=NA_integer_)
        }
      }) %>%
      ungroup() %>%
      mutate(Variant = var_name)
    
    detected <- plateau_by_well %>% filter(plateau_detected, !is.na(plateau_ratio))
    
    if (nrow(detected) > 0) {
      R_max             <- safe_mean(detected$plateau_ratio)
      method_used       <- "plateau_detected"
      plateau_cyc_start <- min(detected$plateau_cyc_start, na.rm=TRUE)
      plateau_cyc_end   <- max(detected$plateau_cyc_end,   na.rm=TRUE)
    } else {
      tail_df <- trusted_df %>%
        group_by(Well) %>%
        mutate(rank_desc = dense_rank(desc(Cycle))) %>%
        filter(rank_desc <= plateau_points) %>%
        ungroup()
      R_max             <- safe_mean(tail_df$R_corr)
      method_used       <- "tail_points_fallback"
      plateau_cyc_start <- min(tail_df$Cycle, na.rm=TRUE)
      plateau_cyc_end   <- max(tail_df$Cycle, na.rm=TRUE)
    }
    
    list(R_max=R_max, method_used=method_used,
         plateau_cyc_start=plateau_cyc_start, plateau_cyc_end=plateau_cyc_end,
         plateau_by_well=plateau_by_well)
  }
  
  rmax_results <- setNames(lapply(all_variants, plateau_detect_one_variant), all_variants)
  
  ## CHANGE 5: compute R_ctrl_plateau — the mean R of each substrate's own
  ## substrate_control wells over the same cycle window used to estimate R_max.
  ## The original had no equivalent; it only compared R_max to R_min (the control
  ## at ref_cycle, i.e. t≈0). 
  ## By comparing to R_ctrl_plateau (concurrent control), the test measures the
  ## genuine cleavage-driven signal above the drifting baseline.
  compute_ctrl_plateau <- function(var_name) {
    res <- rmax_results[[var_name]]
    cs  <- res$plateau_cyc_start
    ce  <- res$plateau_cyc_end
    if (is.na(cs) || is.na(ce)) return(NA_real_)
    merged_raw %>%
      filter(Variant==var_name, Type=="substrate_control", Cycle>=cs, Cycle<=ce) %>%
      mutate(R_ctrl = Intensity_donor / Intensity_acceptor) %>%
      pull(R_ctrl) %>%
      safe_mean()
  }
  
  refs_by_variant <- r_min_by_variant %>%
    mutate(
      R_max          = sapply(Variant, function(v) rmax_results[[v]]$R_max),
      method_used    = sapply(Variant, function(v) rmax_results[[v]]$method_used),
      R_ctrl_plateau = sapply(Variant, compute_ctrl_plateau)
    ) %>%
    mutate(
      ## CHANGE 5: v4 flags the substrate as collapsed and propagates NA downstream,
      ## allowing the pipeline to complete and report all substrates.
      ## The threshold compares the sample plateau to the concurrent control (R_ctrl_plateau).
      signal_above_ctrl = R_max - R_ctrl_plateau,
      min_signal_thresh = min_signal_above_ctrl_frac * R_min,
      collapsed         = is.na(signal_above_ctrl) | signal_above_ctrl < min_signal_thresh,
      method_used       = ifelse(collapsed, "collapsed_dynamic_range", method_used)
    )
  
  wt_plateau_by_well <- bind_rows(lapply(all_variants, function(v) rmax_results[[v]]$plateau_by_well))
  
  list(refs_by_variant = refs_by_variant, wt_plateau_by_well = wt_plateau_by_well)
}

############################################################
## Step 5.  Fraction_cleaved — per Variant, no clamping.
##
## Fraction_cleaved = (R_corr - R_min) / (R_max - R_min)
## - 0 - no substrate cleaved
## - 1 - all substrate cleaved
##
## Notes:
## - Values may fall slightly outside [0,1] due to noise or imperfect
##   dynamic-range estimates; we preserve them for transparency.
## 
## CHANGE 4: R_min and R_max looked up per-Variant via a join on
## refs_by_variant rather than using single global scalars.
## Original: Fraction_cleaved_raw = (R_corr - R_min_global) / (R_max_global - R_min_global)
## v4:       same formula but R_min / R_max are substrate-specific.
##
##
## CHANGE 5: collapsed substrates receive Fraction_cleaved = NA
## instead of a division-by-near-zero result. The original had no
## such guard.
############################################################
compute_fraction_cleaved <- function(corrected, refs_by_variant) {
  corrected %>%
    ## CHANGE 4: join on Variant to get per-substrate R_min / R_max.
    ## Original passed R_min_global and R_max_global as scalar arguments.
    left_join(refs_by_variant %>% select(Variant, R_min, R_max, collapsed), by="Variant") %>%
    mutate(
      Fraction_cleaved = ifelse(
        collapsed | is.na(R_max) | is.na(R_min) | (R_max == R_min),
        NA_real_,
        (R_corr - R_min) / (R_max - R_min)
      )
    ) %>%
    select(Well, Variant, Condition, Type, Cycle, Time_donor_min, Fraction_cleaved, R_corr)
}

############################################################
## Step 6.  Initial cleavage rates per well.
##
## Adaptive biochemical window:
##   - Remove Cycle == 1.
##   - Use all points with Fraction_cleaved <= fraction_ceiling
##     (default 0.30).
##   - If that set has <2 points, fall back to using all
##     remaining post-warmup points.
##   - Fit Fraction_cleaved ~ Time_donor_min.
##   - slope * S0_uM gives rate in µM/min.
##   - Also record n_points_used and which points were used.
##
## CHANGE 7 vs original: the valid-point filter now explicitly
## excludes NA Fraction_cleaved values:
##   filter(Cycle != 1, !is.na(Fraction_cleaved))
## Original only filtered Cycle != 1. If collapsed substrates produce NA
## Fraction_cleaved (Change 5), lm() would receive a data frame with
## zero usable rows; the early-return guard (nrow(df_valid) < 2) now
## correctly catches this and returns NA rate rather than an error or
## artefactual near-zero slope.
##
## Everything else (fraction_ceiling adaptive window, slope * S0_uM
## conversion, fit-point flagging) is unchanged from the original.
############################################################
compute_initial_rates <- function(fraction_table,
                                  S0_uM            = 5,
                                  fraction_ceiling = 0.30,
                                  fit_window_min   = 5) {
  # Creating table to do regression on
  per_well_list <- fraction_table %>%
    group_by(Well, Variant, Condition) %>%
    group_split()
  
  compute_rate_for_group <- function(df_well) {
    
    empty_sum <- tibble(
      Well=df_well$Well[1], Variant=df_well$Variant[1], Condition=df_well$Condition[1],
      slope_fraction_per_min=NA_real_, Initial_rate_uM_per_min=NA_real_, n_points_used=0L
    )
    empty_fp <- df_well %>%
      mutate(used_in_fit=FALSE, fraction_ceiling=fraction_ceiling, fit_window_min=fit_window_min) %>%
      select(Well, Variant, Condition, Cycle, Time_donor_min, Fraction_cleaved,
             used_in_fit, fraction_ceiling, fit_window_min)
    
    ## CHANGE 7: added !is.na(Fraction_cleaved) to the filter.
    ## Original: df_well_valid <- df_well %>% filter(Cycle != 1) %>% arrange(...)
    df_valid <- df_well %>% filter(Cycle != 1, !is.na(Fraction_cleaved)) %>% arrange(Time_donor_min)
    if (nrow(df_valid) < 2) return(list(summary_row=empty_sum, fit_points=empty_fp))
    
    df_pref <- df_valid %>% filter(Fraction_cleaved <= fraction_ceiling) # filtering datapoints below fraction ceiling
    df_fit  <- if (nrow(df_pref) >= 2) df_pref else df_valid  # if there are two or less datapoints below fraction ceiling regression is performed on all datapoints
    
    # Performing regression and finding slope
    slope <- NA_real_
    fit   <- tryCatch(lm(Fraction_cleaved ~ Time_donor_min, data=df_fit), error=function(e) NULL)
    if (!is.null(fit)) {
      cand <- unname(coef(fit)["Time_donor_min"])
      if (!is.na(cand)) slope <- cand
    }
    
    fit_key <- paste(df_fit$Cycle, df_fit$Time_donor_min)
    df_flagged <- df_well %>%
      mutate(
        used_in_fit      = paste(Cycle, Time_donor_min) %in% fit_key,
        fraction_ceiling = fraction_ceiling,
        fit_window_min   = fit_window_min
      ) %>%
      select(Well, Variant, Condition, Cycle, Time_donor_min, Fraction_cleaved,
             used_in_fit, fraction_ceiling, fit_window_min)
    
    list(
      summary_row = tibble(
        Well=df_well$Well[1], Variant=df_well$Variant[1], Condition=df_well$Condition[1],
        slope_fraction_per_min=slope, Initial_rate_uM_per_min=slope*S0_uM, # finding initial rate
        n_points_used=nrow(df_fit)
      ),
      fit_points = df_flagged
    )
  }
  
  all_results     <- per_well_list %>% map(compute_rate_for_group)
  rates_table     <- all_results %>% map("summary_row") %>% bind_rows()
  rate_fit_points <- all_results %>% map("fit_points")  %>% bind_rows()
  list(rates_table=rates_table, rate_fit_points=rate_fit_points)
}

############################################################
## Step 7.  Write QC workbook.
############################################################
write_qc_workbook <- function(plate_map, corrected, refs, fraction_table,
                              rates_table, bleach_deltas_table, rate_fit_points,
                              qc_output_file = "core_output_QC.xlsx") {
  
  bleach_corrected_raw <- corrected %>%
    select(Well, Variant, Condition, Type, Cycle, Time_donor_min,
           Intensity_donor, Intensity_acceptor, donor_corr, acceptor_corr, R_corr)
  
  writexl::write_xlsx(
    list(
      plate_map            = plate_map,
      bleach_corrected_raw = bleach_corrected_raw,
      bleach_deltas        = bleach_deltas_table,
      refs                 = refs,
      fraction_table       = fraction_table,
      rates_table          = rates_table,
      rate_fit_points      = rate_fit_points
    ),
    path = qc_output_file
  )
}

############################################################
## Step 8.  Main wrapper.
##
## Orchestrates everything and returns:
##   plate_map
##   merged_raw
##   corrected
##   bleach_deltas_table
##   refs
##   fraction_table
##   rates_table
##   rate_fit_points
##   wt_plateau_by_well
##
## Also writes qc_output_file.
## 
## CHANGE 7 vs original: parameter wt_variant (default "eTEVp") has
## been replaced by most_cleaved_condition (default "1to1").
## wt_variant was passed to compute_dynamic_range to select trusted
## wells by Variant == wt_variant. In this plate Variant = substrate,
## so "eTEVp" never matched. most_cleaved_condition selects by
## Condition == "1to1" (highest enzyme dilution, most cleavage).
##
## CHANGE 5: new parameter min_signal_above_ctrl_frac (default 0.10)
## passed through to compute_dynamic_range for collapsed detection.
## No equivalent existed in the original.
##
## All other parameters (S0_uM, fit_window_min, fraction_ceiling,
## wt_wells_for_Rmax, plateau_points, plateau_span_cycles,
## plateau_tolerance, qc_output_file) are unchanged.
############################################################
run_core_pipeline_from_tidy <- function(
    tidy_csv_file,
    plate_map_file,
    S0_uM                      = 5,
    fit_window_min             = 5,
    fraction_ceiling           = 0.50,
    most_cleaved_condition     = "1to1",    ## CHANGE 9: replaces wt_variant = "eTEVp"
    wt_wells_for_Rmax          = NULL,
    plateau_points             = 3,
    plateau_span_cycles        = 10,
    plateau_tolerance          = 0.01,
    min_signal_above_ctrl_frac = 0.10,      ## CHANGE 5: new parameter, no original equivalent
    qc_output_file             = "core_output_QC.xlsx"
) {
  
  inputs     <- read_and_merge_inputs(tidy_csv_file, plate_map_file)
  plate_map  <- inputs$plate_map
  merged_raw <- inputs$merged_raw
  
  bleach_out          <- apply_bleach_correction(merged_raw)
  corrected           <- bleach_out$corrected
  ref_cycle           <- bleach_out$ref_cycle
  bleach_deltas_table <- bleach_out$bleach_deltas_table
  
  dyn_range <- compute_dynamic_range(
    merged_raw                 = merged_raw,
    corrected                  = corrected,
    ref_cycle                  = ref_cycle,
    bleach_deltas_table        = bleach_deltas_table,
    most_cleaved_condition     = most_cleaved_condition,
    wt_wells_for_Rmax          = wt_wells_for_Rmax,
    plateau_points             = plateau_points,
    plateau_span_cycles        = plateau_span_cycles,
    plateau_tolerance          = plateau_tolerance,
    min_signal_above_ctrl_frac = min_signal_above_ctrl_frac
  )
  refs_by_variant <- dyn_range$refs_by_variant
  
  fraction_table <- compute_fraction_cleaved(corrected, refs_by_variant)
  
  rate_out        <- compute_initial_rates(fraction_table, S0_uM, fraction_ceiling, fit_window_min)
  rates_table     <- rate_out$rates_table
  rate_fit_points <- rate_out$rate_fit_points
  
  refs <- refs_by_variant %>%
    mutate(
      ref_cycle                  = ref_cycle,
      most_cleaved_condition     = most_cleaved_condition,
      plateau_points             = plateau_points,
      plateau_span_cycles        = plateau_span_cycles,
      plateau_tolerance          = plateau_tolerance,
      min_signal_above_ctrl_frac = min_signal_above_ctrl_frac,
      S0_uM                      = S0_uM,
      fit_window_min             = fit_window_min,
      fraction_ceiling           = fraction_ceiling
    )
  
  write_qc_workbook(
    plate_map=plate_map, corrected=corrected, refs=refs,
    fraction_table=fraction_table, rates_table=rates_table,
    bleach_deltas_table=bleach_deltas_table, rate_fit_points=rate_fit_points,
    qc_output_file=qc_output_file
  )
  
  list(
    plate_map           = plate_map,
    merged_raw          = merged_raw,
    corrected           = corrected,
    bleach_deltas_table = bleach_deltas_table,
    refs                = refs,
    fraction_table      = fraction_table,
    rates_table         = rates_table,
    rate_fit_points     = rate_fit_points,
    wt_plateau_by_well  = dyn_range$wt_plateau_by_well
  )
}

############################################################
## End of core_analysis_v4.R
############################################################
