# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT --------------------------------------
# ============================================================================ #
# Dieses Skript stellt Funktionen zur Berechnung verschiedener statistischer 
# Kennzahlen fuer die Lasso Regression bereit. Dazu gehoeren:
# - Bias, Root Mean Square Error (RMSE), Median, IQR (Erstes Quartil (Q1); Drittes Quartil (Q3),
#   Range (Minimum, Maximum)
# - Ueberdeckungswahrscheinlichkeit (Coverage Probability)
# - Trefferquote (Hit Rate), Falschalarmquote (False Alarm Rate)
# - Lambda
#
# Am Ende des Skripts werden die berechneten Ergebnisse in einem DataFrame 
# zusammengefuehrt.
#
# Zusaetzlich werden die Ergebnisse fall- und kennzahlenspezifisch in neun 
# separate Word-Dokumente exportiert:
# - 3 Dokumente (fuer Fall A-C) - Guete von Punktschaetzern:       Mean, Bias, RMSE, Median, IQR, Range und Lambda
# - 3 Dokumente (fuer Fall A-C) - Guete von Konfidenzintervallen:  Coverage Probabilities (CPs)
# - 3 Dokumente (fuer Fall A-C) - Guete von Variablenselektionn:   Hit Rate (HR) und False Alarm Rate (FAR)
# ============================================================================ #

# ============================================================================ #
# BENOETIGTE PAKETE INSTALLIEREN UND LADEN -----------------------------------
# ============================================================================ #

# install.packages("dplyr")
# install.packages("officer")
# install.packages("flextable")

library(dplyr)  
library(officer)
library(flextable)

# ============================================================================ #
# FUNKTIONEN DEFINIEREN ------------------------------------------------------
# ============================================================================ #

# Funktion zur Rundung auf 2 Nachkommastellen
# Diese Funktion wird verwendet, um alle berechneten Werte mit einer Genauigkeit 
# von zwei Nachkommastellen zu formatieren.
format_value <- function(x, digits = 2) {
  round(x, digits)
}

# Funktion zur Berechnung des Bias
# Bias misst die durchschnittliche Differenz zwischen den geschaetzten Werten und 
# den wahren Werten. Ein Bias von null zeigt eine unverzerrte Schaetzung an.
bias <- function(true, est) {
  format_value(mean(est - true))
}

# Funktion zur Berechnung des Root Mean Square Error (RMSE)
# RMSE misst die durchschnittliche quadratische Abweichung der Schaetzwerte von den wahren 
# Werten. Es ist ein Mass fuer die Genauigkeit der Schaetzungen.
rmse <- function(true, est) {
  format_value(sqrt(mean((est - true)^2)))
}

# Funktion zur Berechnung von Median, IQR (Q1; Q3), Min und Max sowie Range
calc_stats <- function(data) {
  q1 <- format_value(quantile(data, 0.25))    # Q1
  q3 <- format_value(quantile(data, 0.75))    # Q3
  median_value <- format_value(median(data))  # Median
  min_value <- format_value(min(data))        # Minimum
  max_value <- format_value(max(data))        # Maximum
  
  # Median und IQR formatieren
  median_iqr <- paste0(median_value, " [", q1, ";", q3, "]")
  
  # Range formatieren
  range_value <- paste0(min_value, ",", max_value)
  
  return(c(Median_IQR = median_iqr, Range = range_value))
}

# Funktion zur Berechnung der Coverage Probability (Geschaetzte Ueberdeckungswahrscheinlichkeit)
# Diese Funktion berechnet den Anteil der Faelle, bei denen die wahren Werte 
# innerhalb der Konfidenzintervalle liegen.
coverage_prob <- function(ci) {
  format_value(mean(ci == TRUE) * 100)
}

# Funktion zur Berechnung der Hit Rate fuer Effektstaerken von relevanten Dummy-Variablen
# Diese Funktion berechnet den Anteil der Faelle, bei denen 0 bei relevanten 
# Variablen nicht vom KI eingeschlossen wird
hit_rate_effect <- function(true_val, ci_lower, ci_upper) {
  # Nur relevante Faelle (true_val != 0) betrachten
  relevant <- true_val != 0
  # Berechnung, ob 0 ausserhalb des Konfidenzintervalls liegt = HIT
  hits <- (ci_lower[relevant] > 0) | (ci_upper[relevant] < 0)
  # Berechnung der Hit Rate fuer relevante Faelle
  format_value(mean(hits) * 100)
}
# Falls 0 nicht eingeschlossen wird: Variable richtigerweise als richtig erkannt (TP/H)
# Falls 0 eingeschlossen wird: Variable faelschlicherweise als falsch erkannt (FN)
# TP: True Positive, H: Hit, FN: False Negative

# Funktion zur Berechnung der False Alarm Rate fuer Effektstaerken irrelevanter Dummy-Variablen
# Diese Funktion berechnet den Anteil der Faelle, bei denen 0 nicht vom KI bei irrelevanten
# Variablen nicht vom KI eingeschlossen wird
false_alarm_rate_effect <- function(true_val, ci_lower, ci_upper) {
  # Nur irrelevante Faelle (true_val == 0) betrachten
  irrelevant <- true_val == 0
  # Berechnung, ob 0 ausserhalb des Konfidenzintervalls liegt = FALSE ALARM
  false_alarms <- (ci_lower[irrelevant] > 0) | (ci_upper[irrelevant] < 0)
  # Berechnung der False Alarm Rate
  format_value(mean(false_alarms) * 100)
}
# Falls 0 nicht eingeschlossen wird: Variable faelschlicherweise als richtig erkannt (FP/FA)
# Falls 0 eingeschlossen wird: Variable richtigerweise als falsch erkannt (TN)
# TN: True Negative, FP: False Positive, FA: False Alarm

# ============================================================================ #
# DATEN LADEN ----------------------------------------------------------------
# ============================================================================ #
load("4.2_RESULTS_SIMSTUD_LASSO.RData")   # simulierte Daten laden
metric_data_lasso <- final.results.lasso  # Daten in das Objekt 'metric_data_lasso' speichern

# Ueberpruefen, ob fehlende Werte in den geladenen Daten vorhanden sind.
any(is.na(metric_data_lasso))

# ============================================================================ #
# PARAMETER UND LEERE LISTE ERSTELLEN ----------------------------------------
# ============================================================================ #
# Definition der verschiedenen Kategorien von CASE.ID, REF.ID und n.total.
case_ids <- c("A", "B", "C1", "C2")   # Fall-IDs
ref_ids <- c("REFA", "REFB", "REFC")  # Referenz-IDs
n_total_vals <- c(50, 500)            # Stichprobengroessen

# Erstellung einer leeren Liste zur Speicherung der Ergebnisse
results_list <- list()

# ============================================================================ #
# SCHLEIFEN UEBER ALLE FAKTOR-KOMBINATIONEN -----------------------------------
# ============================================================================ #
# In dieser Schleife werden die Daten fuer jede Kombination von CASE.ID, REF.ID und n.total 
# gefiltert und die entsprechenden statistischen Kennzahlen berechnet.
for (case_id in case_ids) {
  for (ref_id in ref_ids) {
    for (n_total in n_total_vals) {
      
      # Filtere die Daten fuer die aktuelle Kombination von CASE.ID, REF.ID und n.total
      filtered_data_lasso <- subset(metric_data_lasso, CASE.ID == case_id & REF.ID == ref_id & n.total == n_total)
      
      # Berechnung der Mittelwerte der geschaetzten Achsenabschnitt- und Effektstaerken
      mean_est_int <- format_value(mean(filtered_data_lasso$EST.INT))
      mean_est_slope_one <- format_value(mean(filtered_data_lasso$EST.SLOPE.ONE))
      mean_est_slope_two <- format_value(mean(filtered_data_lasso$EST.SLOPE.TWO))
      
      # Berechnung von Bias fuer Achsenabschnitt und Effektstaerken
      bias_int <- bias(filtered_data_lasso$TRUE.INT, filtered_data_lasso$EST.INT)
      bias_slope_one <- bias(filtered_data_lasso$TRUE.SLOPE.ONE, filtered_data_lasso$EST.SLOPE.ONE)
      bias_slope_two <- bias(filtered_data_lasso$TRUE.SLOPE.TWO, filtered_data_lasso$EST.SLOPE.TWO)
      
      # Berechnung von RMSE fuer Achsenabschnitt und Effektstaerken
      rmse_int <- rmse(filtered_data_lasso$TRUE.INT, filtered_data_lasso$EST.INT)
      rmse_slope_one <- rmse(filtered_data_lasso$TRUE.SLOPE.ONE, filtered_data_lasso$EST.SLOPE.ONE)
      rmse_slope_two <- rmse(filtered_data_lasso$TRUE.SLOPE.TWO, filtered_data_lasso$EST.SLOPE.TWO)
      
      # Berechnung von Median und IQR fuer Achsenabschnitt und Effektstaerken
      stats_int <- calc_stats(filtered_data_lasso$EST.INT)
      stats_slope_one <- calc_stats(filtered_data_lasso$EST.SLOPE.ONE)
      stats_slope_two <- calc_stats(filtered_data_lasso$EST.SLOPE.TWO)
      
      # Coverage Probability fuer Achsenabschnitt und Effektstaerken mit Bootstrap-Konfidenzintervallen
      coverage_prob_int <- coverage_prob(filtered_data_lasso$BOOT.CI.INT)
      coverage_prob_slope_one <- coverage_prob(filtered_data_lasso$BOOT.CI.SLOPE.ONE)
      coverage_prob_slope_two <- coverage_prob(filtered_data_lasso$BOOT.CI.SLOPE.TWO)
      
      # Hit Rate und False Alarm Rate fuer Slope One mit Bootstrap-Konfidenzintervallen
      hit_rate_slope_one <- hit_rate_effect(filtered_data_lasso$TRUE.SLOPE.ONE, 
                                            filtered_data_lasso$BOOT.CI.SLOPE.ONE.LOWER, 
                                            filtered_data_lasso$BOOT.CI.SLOPE.ONE.UPPER)
      false_alarm_rate_slope_one <- false_alarm_rate_effect(filtered_data_lasso$TRUE.SLOPE.ONE, 
                                                            filtered_data_lasso$BOOT.CI.SLOPE.ONE.LOWER, 
                                                            filtered_data_lasso$BOOT.CI.SLOPE.ONE.UPPER)
      
      # Hit Rate und False Alarm Rate fuer Slope Two mit Bootstrap-Konfidenzintervallen
      hit_rate_slope_two <- hit_rate_effect(filtered_data_lasso$TRUE.SLOPE.TWO, 
                                            filtered_data_lasso$BOOT.CI.SLOPE.TWO.LOWER, 
                                            filtered_data_lasso$BOOT.CI.SLOPE.TWO.UPPER)
      false_alarm_rate_slope_two <- false_alarm_rate_effect(filtered_data_lasso$TRUE.SLOPE.TWO, 
                                                            filtered_data_lasso$BOOT.CI.SLOPE.TWO.LOWER, 
                                                            filtered_data_lasso$BOOT.CI.SLOPE.TWO.UPPER)
      
      # Mittleres Lambda
      mean_lambda <- format_value(mean(filtered_data_lasso$BEST.LAMBDA))
      
      # Speichern der Ergebnisse fuer die aktuelle Kombination in einer Liste
      results_list_lasso[[paste(case_id, ref_id, n_total, sep = "_")]] <- list(
        
        Case = case_id, Reference = ref_id, n = n_total,
        
        # Guete von Punktschaetzern
        M_INT = mean_est_int, M_SlopeOne = mean_est_slope_one, M_SlopeTwo = mean_est_slope_two,
        Bias_INT = bias_int, Bias_SlopeOne = bias_slope_one, Bias_SlopeTwo = bias_slope_two,
        RMSE_INT = rmse_int, RMSE_SlopeOne = rmse_slope_one, RMSE_SlopeTwo = rmse_slope_two,
        Median_IQR_INT = stats_int["Median_IQR"], Range_INT = stats_int["Range"],
        Median_IQR_SlopeOne = stats_slope_one["Median_IQR"], Range_SlopeOne = stats_slope_one["Range"],
        Median_IQR_SlopeTwo = stats_slope_two["Median_IQR"], Range_SlopeTwo = stats_slope_two["Range"],
        
        # Guete von Konfidenzintervallen
        CP_INT_BOOT = coverage_prob_int, CP_SlopeOne_BOOT = coverage_prob_slope_one, CP_SlopeTwo_BOOT = coverage_prob_slope_two,
        
        # Guete von Variablenselektion
        HR_SlopeOne_BOOT = hit_rate_slope_one, FAR_SlopeOne_BOOT = false_alarm_rate_slope_one,
        HR_SlopeTwo_BOOT = hit_rate_slope_two, FAR_SlopeTwo_BOOT = false_alarm_rate_slope_two,
        
        # Mittleres Lambda
        M_Lambda = mean_lambda
      )
    }
  }
}

# ============================================================================ #
# ERGEBNISSE SPEICHERN -------------------------------------------------------
# ============================================================================ #

# Zusammenfuehren der Ergebnisse in ein DataFrame
lasso.analysis.results <- do.call(rbind, lapply(results_list_lasso, as.data.frame))

# Ergebnisse anzeigen
View(lasso.analysis.results)

# Speichern der Ergebnisse in einer RData-Datei
save(lasso.analysis.results, file = "6.1_METRIC_ANALYSIS_LASSO.RData")

# ============================================================================ #
# ERGEBNISSE ZUSAETZLICH NACH FALL EXPORTIEREN --------------------------------
# ============================================================================ #

# Ordner fuer die Ergebnisse
# (mit der Voraussetzung, dass der Ordner bereits im selben Verzeichnis existiert)
output_dir_lasso <- "6.2_Tables_LassoRegression"

# Fuer jeden Fall die Ergebnisse in separate Dokumente exportieren
for (case_id in case_ids) {
  
  # Daten nach dem aktuellen Fall filtern
  case_data_lasso <- subset(lasso.analysis.results, Case == case_id)
  
  # Erstellung von DataFrames fuer die verschiedenen Kennzahlen:
  
  # Guete von Punktschaetzern und weitere deskriptive Angaben
  mean_bias_rmse_median_iqr_range_df_lasso <- case_data_lasso[, c("Case", "Reference", "n", 
                                                                  "M_INT", "M_SlopeOne", "M_SlopeTwo", 
                                                                  "Bias_INT", "Bias_SlopeOne", "Bias_SlopeTwo", 
                                                                  "RMSE_INT", "RMSE_SlopeOne", "RMSE_SlopeTwo",
                                                                  "Median_IQR_INT", "Range_INT",
                                                                  "Median_IQR_SlopeOne", "Range_SlopeOne",
                                                                  "Median_IQR_SlopeTwo", "Range_SlopeTwo",
                                                                  "M_Lambda")]
  
  # Guete von Konfidenzintervallen
  cp_df_lasso <- case_data_lasso[, c("Case", "Reference", "n",
                                     "CP_INT_BOOT", "CP_SlopeOne_BOOT", "CP_SlopeTwo_BOOT")]
  
  # Guete von Variablenselektion
  hr_far_df_lasso <- case_data_lasso[, c("Case", "Reference", "n",
                                         "HR_SlopeOne_BOOT", "FAR_SlopeOne_BOOT", 
                                         "HR_SlopeTwo_BOOT", "FAR_SlopeTwo_BOOT")]
  
  # Export in drei separate Word-Dokumente
  
  # Dokument fuer Mean, Bias, RMSE, Median, IQR, Range und Lambda
  doc1 <- read_docx()
  doc1 <- body_add_flextable(doc1, flextable(mean_bias_rmse_median_iqr_range_df_lasso))
  print(doc1, target = file.path(output_dir_lasso, 
                                 paste0("Case_", case_id, "_A_Mean_Bias_RMSE_Median_IQR_Range_Lambda_Lasso.docx")))
  
  # Dokument fuer Coverage Probabilities
  doc2 <- read_docx()
  doc2 <- body_add_flextable(doc2, flextable(cp_df_lasso))
  print(doc2, target = file.path(output_dir_lasso, 
                                 paste0("Case_", case_id, "_B_Coverage_Probabilities_Lasso.docx")))
  
  # Dokument fuer Hit Rate und False Alarm Rate
  doc3 <- read_docx()
  doc3 <- body_add_flextable(doc3, flextable(hr_far_df_lasso))
  print(doc3, target = file.path(output_dir_lasso,
                                 paste0("Case_", case_id, "_C_HitRate_FalseAlarmRate_Lasso.docx")))
}

# ============================================================================ #
# ENDE DES SKRIPTS -----------------------------------------------------------
# ============================================================================ #