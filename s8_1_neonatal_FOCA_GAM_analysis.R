## ================================================================
## Script: neonatal_FOCA_GAM_analysis.R
## Purpose:
##   Model age-related trajectories of internetwork FOCA (Adversarial Connectivity)
##   during the early postnatal period using mass univariate GAMs.
##
## Description:
##   - For each pair among 20 functional networks, fit a GAM:
##       FOCA_ij ~ s(PMA, k=3) + Sex + Motion + ScanBirthInterval
##   - Extract F-statistic, p-value, and Z-value (directional)
##   - Save results as .mat for each network pair
##
## Note:
##   All paths are relative to the current working directory.
##
## ================================================================

## --- Load required packages ---
library(R.matlab)
library(mgcv)
library(nlme)
library(visreg)

## --- Define relative paths ---
ResultPath <- "./outputs/AgeEffects_FOCA"
AtlasLoading_Folder <- "./data/dHCP/FOCA_AgeInputs"

## --- Define network names ---
network_names <- c("DMN Par","DMN Ant","DMN Dor","DMN Ret",
                   "Vis Lat","Vis Str","Vis V5","Vis V1",
                   "FP","DAN","PreMot","Lang","Sal","AMN",
                   "PMN","Hand","Face","Foot","Aud","SCAN")

## --- Create result folder if not exists ---
EffectResultDir <- file.path(ResultPath, "EffectResult")
if (!dir.exists(EffectResultDir)) dir.create(EffectResultDir, recursive = TRUE)

cat("=== Starting Neonatal FOCA GAM Analysis ===\n")

## --- Main analysis loop ---
for (i in seq_along(network_names)) {
  for (j in seq_along(network_names)) {

    ## Skip upper triangle and diagonal
    if (i <= j) next

    cat("Processing: Network_", i, "_", j, "\n")

    ## --- Step 1. Load FOCA (Adversarial Connectivity) data ---
    data_path <- file.path(AtlasLoading_Folder, "ResultLynch20",
                           paste0("Individual_AC_", i, "_", j, "_DeleteSubject.mat"))
    if (!file.exists(data_path)) {
      warning(paste("Missing file:", data_path))
      next
    }
    Measure_All <- readMat(data_path)
    Measure <- Measure_All$Adversarial.Connectivity

    ## Remove zero columns (subjects with all zeros)
    valid_cols <- which(colSums(Measure) != 0)
    Data_Valid <- as.numeric(Measure[, valid_cols])

    ## --- Step 2. Load covariates ---
    behavior_path <- file.path(AtlasLoading_Folder, "ResultLynch20",
                               paste0("AgeScore_", i, "_", j, ".mat"))
    if (!file.exists(behavior_path)) {
      warning(paste("Missing file:", behavior_path))
      next
    }
    Behavior_All <- readMat(behavior_path)
    Behavior <- data.frame(
      scan_age = as.numeric(Behavior_All$AgeScore[[1]]),
      sex = as.factor(Behavior_All$AgeScore[[2]]),
      mFD = as.numeric(Behavior_All$AgeScore[[3]]),
      TimeInterval = as.numeric(Behavior_All$AgeScore[[4]])
    )

    ## --- Step 3. Fit GAM model ---
    Gam_Model <- gam(Data_Valid ~ s(scan_age, bs = "cs", k = 3) +
                       sex + TimeInterval + mFD,
                     data = Behavior, method = "REML", na.action = na.omit)

    ## --- Step 4. Extract statistics ---
    Gam_Summary <- summary(Gam_Model)
    Gam_F_Vector <- Gam_Summary$s.table["s(scan_age)", "F"]
    Gam_P_Vector <- Gam_Summary$s.table["s(scan_age)", "p-value"]

    ## Compute Z value (direction determined by linear model)
    lm_Model <- lm(Data_Valid ~ scan_age + sex + mFD, data = Behavior)
    Age_T <- summary(lm_Model)$coefficients["scan_age", "t value"]
    Gam_Z_Vector <- qnorm(Gam_P_Vector / 2, lower.tail = FALSE)
    if (Age_T < 0) Gam_Z_Vector <- -Gam_Z_Vector

    ## --- Step 5. Save results ---
    output_path <- file.path(EffectResultDir,
                             paste0("Network_", i, "_", j, "_AgeEffect.mat"))
    writeMat(output_path,
             Gam_P_Vector = Gam_P_Vector,
             Gam_F_Vector = Gam_F_Vector,
             Gam_Z_Vector = Gam_Z_Vector)
  }
}

cat("=== GAM analysis complete. Results saved to:", EffectResultDir, "===\n")
