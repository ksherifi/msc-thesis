# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT --------------------------------------
# ============================================================================ #
# Dieses Skript erstellt Histogramme zur Visualisierung der (asymptotischen) Normalverteilung
# der Punktschaetzer fuer beide Regressionsverfahren in den Faellen A und B sowie 
# den Sub-Faellen C1 und C2. Dabei erfolgt die Visualisierung der (asymptotischen) 
# Normalverteilung  separat fuer die Parameter Achsenabschnitt, Effektstaerke D1 
# und Effektstaerke D2.
# 
# Die Simulationsergebnisse werden zuerst fuer beide Regresssionsverfahren
# geladen und in einem kombinierten Datensatz zusammengefuehrt. Anschliessend werden die
# Daten in ein langes Format transformiert, um die Parameterwerte uebersichtlich darzustellen.
# Fuer jede Referenzkategorie und jeden Fall werden die wahren Parameterwerte definiert und
# als vertikale Referenzlinien in die Histogramme eingefuegt.
# 
# Die Histogramme zeigen die Verteilungen der Punktschaetzer um den wahren Parameterwert und erlauben
# den direkten Vergleich der Methoden fuer verschiedene Stichprobengroessen (n=50 und n=500).
# Die x-Achsen fuer jede Facette sind individuell auf den Bereich der wahren Werte ±1 skaliert,
# um die Ergebnisse besser vergleichbar zu machen. Die Ausgabe erfolgt als JPEG-Dateien fuer
# jeden Fall und wird im Verzeichnis "9.1_Histograms" gespeichert.
# ============================================================================ #

# ============================================================================ #
# BENOETIGTE PAKETE INSTALLIEREN UND LADEN -----------------------------------
# ============================================================================ #

# install.packages("ggplot2")   # Zum Erstellen und Gestalten der Histogramme
# install.packages("dplyr")     # Zum Bearbeiten und Filtern der Daten
# install.packages("tidyr")     # Zum Umformen der Daten in ein langes Format
# install.packages("ggh4x")     # Erweitert ggplot2, z.B. fuer separate y-Achsen pro Facette
# install.packages("purrr")     # Fuer das Erstellen dynamischer Achsenskalierungen
# install.packages("tibble")    # Erleichtert die Arbeit mit Tabellenformaten fuer die wahren Werte

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(purrr)
library(tibble)

# ============================================================================ #
# DATEN LADEN UND VORBEREITEN ------------------------------------------------
# ============================================================================ #

# Laden der gespeicherten Simulationsergebnisse fuer multiple Regression und Lasso
load("4.1_RESULTS_SIMSTUD_LINREG.RData")  # Ergebnisse der multiplen Regression laden
load("4.2_RESULTS_SIMSTUD_LASSO.RData")   # Ergebnisse der Lasso Regression laden

# Zuweisung der geladenen Daten an spezifische Objekte fuer weitere Verarbeitung
data_linreg <- final.results.linreg
data_lasso <- final.results.lasso

# ============================================================================ #
# DATEN AUFBEITEN UND KOMBINIEREN --------------------------------------------
# ============================================================================ #

# Auswahl der relevanten Spalten fuer Histogramme zur multiplen Regression
data_histograms_linreg <- data_linreg[, c("n.total", "CASE.ID", "REF.ID", 
                                          "EST.INT", "EST.SLOPE.ONE", "EST.SLOPE.TWO", "METHOD")]

# Auswahl der relevanten Spalten fuer Histogramme zur Lasso Regression
data_histograms_lasso <- data_lasso[, c("n.total", "CASE.ID", "REF.ID", 
                                        "EST.INT", "EST.SLOPE.ONE", "EST.SLOPE.TWO", "METHOD")]

# ============================================================================ #
# KOMBINIEREN BEIDER DATENSAETZE IN EINEM DATAFRAME --------------------------
# ============================================================================ #

# Zusammenfuehren der beiden Datensaetze in einen einzigen DataFrame
combined_data <- bind_rows(data_histograms_linreg, data_histograms_lasso)

# ============================================================================ #
# ANPASSEN DER 'METHOD'-SPALTE FUER KONSISTENTE BEZEICHNUNGEN ----------------
# ============================================================================ #

# Umkodierung der 'METHOD'-Spalte in Faktoren mit konsistenten und klaren Labels
combined_data <- combined_data %>%
  mutate(METHOD = factor(METHOD, 
                         levels = c("LinReg", "Lasso"),
                         labels = c("Multiple Regression", "Lasso Regression")))

# ============================================================================ #
# DATEN TRANSFORMIEREN IN LONG-FORMAT ----------------------------------------
# ============================================================================ #

# Umstrukturierung der kombinierten Daten in ein langes Format, um die Parameter als separate Zeilen zu haben
combined_long <- combined_data %>%
  pivot_longer(
    cols = starts_with("EST"),  # Auswahl der EST-Variablen (EST.INT, EST.SLOPE.ONE, EST.SLOPE.TWO)
    names_to = "Parameter",     # Name der neuen Spalte fuer Parameter
    values_to = "value"         # Name der neuen Spalte fuer Werte
  ) %>%
  mutate(
    # Umkodierung der 'Parameter'-Spalte in Faktoren mit verstaendlichen Labels
    Parameter = factor(Parameter, 
                       levels = c("EST.INT", "EST.SLOPE.ONE", "EST.SLOPE.TWO"),
                       labels = c("Achsenabschnitt", "Effektstärke D1", "Effektstärke D2")),
    # Umkodierung der 'REF.ID'-Spalte in Faktoren mit verstaendlichen Labels
    REF.ID = factor(REF.ID, 
                    levels = c("REFA", "REFB", "REFC"),
                    labels = c("Referenz A", "Referenz B", "Referenz C")),
    # Umkodierung von 'n.total' in einen Faktor
    n.total = factor(n.total)
  )

# ============================================================================ #
# DEFINITION DER WAHREN WERTE FUER ALLE FAELLE -------------------------------
# ============================================================================ #

# Wahre Werte fuer jeden Fall und jede Referenzkategorie
true_values_all_cases <- list(
  
  # Fall A: Wahre Parameterwerte fuer alle Referenzen
  A = list(
    "Referenz A" = c("Achsenabschnitt" = 5, "Effektstärke D1" = 0, "Effektstärke D2" = 0),
    "Referenz B" = c("Achsenabschnitt" = 5, "Effektstärke D1" = 0, "Effektstärke D2" = 0),
    "Referenz C" = c("Achsenabschnitt" = 5, "Effektstärke D1" = 0, "Effektstärke D2" = 0)
  ),
  
  # Fall B: Wahre Parameterwerte fuer alle Referenzen
  B = list(
    "Referenz A" = c("Achsenabschnitt" = 5, "Effektstärke D1" = 5, "Effektstärke D2" = -5),
    "Referenz B" = c("Achsenabschnitt" = 10, "Effektstärke D1" = -5, "Effektstärke D2" = -10),
    "Referenz C" = c("Achsenabschnitt" = 0, "Effektstärke D1" = 5, "Effektstärke D2" = 10)
  ),
  
  # Sub-Fall C1: Wahre Parameterwerte fuer alle Referenzen
  C1 = list(
    "Referenz A" = c("Achsenabschnitt" = 10, "Effektstärke D1" = 0, "Effektstärke D2" = -5),
    "Referenz B" = c("Achsenabschnitt" = 10, "Effektstärke D1" = 0, "Effektstärke D2" = -5),
    "Referenz C" = c("Achsenabschnitt" = 5, "Effektstärke D1" = 5, "Effektstärke D2" = 5)
  ),
  
  # Sub-Fall C2: Wahre Parameterwerte fuer alle Referenzen
  C2 = list(
    "Referenz A" = c("Achsenabschnitt" = 10, "Effektstärke D1" = 0, "Effektstärke D2" = -10),
    "Referenz B" = c("Achsenabschnitt" = 10, "Effektstärke D1" = 0, "Effektstärke D2" = -10),
    "Referenz C" = c("Achsenabschnitt" = 0, "Effektstärke D1" = 10, "Effektstärke D2" = 10)
  )
  
)

# Umwandlung der Liste der wahren Werte in ein langes DataFrame zur einfachen Verknuepfung mit den geschaetzten Werten
true_vals_long_all <- imap_dfr(true_values_all_cases, function(case_values, case_id) {
  tibble(
    REF.ID = rep(names(case_values), each = 3),
    Parameter = rep(names(case_values[[1]]), times = 3),
    true_value = unlist(case_values)
  )
}, .id = "CASE.ID")

# ============================================================================ #
# VERBINDEN DER GESCHAETZTEN UND WAHREN WERTE --------------------------------
# ============================================================================ #

# Verknuepfung der geschaetzten Werte mit den wahren Werten anhand von CASE.ID, REF.ID und Parameter
data_with_true <- combined_long %>%
  left_join(true_vals_long_all, by = c("CASE.ID", "REF.ID", "Parameter")) %>%
  # Erstellung einer neuen Variable, die Methode und Stichprobengroesse kombiniert, um verschiedene Fuellfarben zu ermoeglichen
  mutate(METHOD_SAMPLE = interaction(METHOD, n.total, sep = "_"))

# ============================================================================ #
# DEFINIEREN DER BENUTZERSPEZIFISCHEN FARBEN FUER DIE HISTOGRAMME ------------
# ============================================================================ #

# Festlegen der Farben fuer die verschiedenen Methoden und Stichprobengroessen
colors_custom <- scale_fill_manual(
  values = c(
    "Multiple Regression_50" = "red3",
    "Multiple Regression_500" = "darkred",
    "Lasso Regression_50" = "dodgerblue",
    "Lasso Regression_500" = "darkblue"
  ),
  # Festlegen der Beschriftungen fuer die Legende unter Verwendung von Ausdrucksfunktionen fuer mathematische Symbole
  labels = c(
    "Multiple Regression_50" = expression(italic(n) * " = 50"),
    "Multiple Regression_500" = expression(italic(n) * " = 500"),
    "Lasso Regression_50" = expression(italic(n) * " = 50"),
    "Lasso Regression_500" = expression(italic(n) * " = 500")
  ),
  name = "Stichprobengrösse"  # Titel der Legende
)

# ============================================================================ #
# DEFINIEREN DER HISTOGRAMM-ERSTELLEN FUNKTION ------------  *+*+*KI-HELP*+*+*
# ============================================================================ #

# Funktion zur Erstellung und Speicherung von Histogrammen fuer spezifische Methoden und Faelle
create_histogram <- function(data_case, method_filter, case_filter, output_dir = "9.1_Histograms") {
  
  # Filtern der Daten basierend auf der gewaehlten Methode und dem Fall
  data_filtered <- data_case %>%
    filter(METHOD == method_filter, CASE.ID == case_filter)
  
  # Extrahieren der wahren Werte fuer den aktuellen Fall und Methode zur Darstellung als vertikale Linien
  true_vals_case <- data_filtered %>%
    select(REF.ID, Parameter, true_value) %>%
    distinct() %>%
    arrange(REF.ID, Parameter) %>%
    mutate(Facet_ID = paste(REF.ID, Parameter, sep = "_")) %>%
    mutate(
      lower = true_value - 1,  # Untere Grenze fuer die x-Achse
      upper = true_value + 1   # Obere Grenze fuer die x-Achse
    )
  
  # Erstellung einer Liste von Skalen fuer jede Facette, um individuelle x-Achsen-Beschraenkungen zu ermoeglichen
  facet_scales_x <- setNames(
    map2(true_vals_case$lower, true_vals_case$upper, ~ scale_x_continuous(limits = c(.x, .y))),
    true_vals_case$Facet_ID
  )
  
  # Erstellung des Histogramms mit ggplot2
  plot <- ggplot(data_filtered, aes(x = value, fill = METHOD_SAMPLE)) +
    geom_histogram(aes(y = ..count..), alpha = 0.6, position = "identity", bins = 30) +  # Histogramm der geschaetzten Werte
    geom_vline(aes(xintercept = true_value), color = "darkred", linetype = "solid", size = 0.3) +  # Linie fuer den wahren Wert
    facet_wrap2(~ REF.ID + Parameter, nrow = 3, scales = "free_x", axes = "x") +  # Facettierung nach Referenz und Parameter
    facetted_pos_scales(x = facet_scales_x) +  # Anwenden der individuellen Skalen
    colors_custom +  # Anwenden der benutzerdefinierten Farben
    labs(
      x = "Verteilung der Punktschätzer um den Wahren Parameterwert",  # Beschriftung der x-Achse
      y = "Häufigkeit"  # Beschriftung der y-Achse
    ) +
    theme_minimal(base_size = 16) +  # Verwendung eines minimalen Themas mit Basis-Schriftgroesse 16
    theme(
      axis.title.x = element_text(face = "bold", size = 20, margin = margin(t = 25)),  # Formatierung der x-Achsentitel
      axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 25)),  # Formatierung der y-Achsentitel
      
      axis.text.x = element_text(face = "bold", size = 13),  # Formatierung der x-Achsentexte
      axis.text.y = element_text(face = "bold", size = 13),  # Formatierung der y-Achsentexte
      
      axis.ticks.x = element_line(size = 0.8),  # Formatierung der x-Achsen-Ticks
      axis.ticks.y = element_line(size = 0.8),  # Formatierung der y-Achsen-Ticks
      axis.minor.ticks.x.bottom = element_line(color = "black", size = 0.4),  # Formatierung der kleinen x-Ticks
      axis.minor.ticks.length.x.bottom = unit(2, "mm"),  # Laenge der kleinen x-Ticks
      axis.minor.ticks.y.left = element_line(color = "black", size = 0.4),  # Formatierung der kleinen y-Ticks
      axis.minor.ticks.length.y.left = unit(2, "mm"),  # Laenge der kleinen y-Ticks
      
      legend.margin = margin(10, 10, 10, 10),  # Abstand der Legende
      legend.background = element_rect(color = "black", fill = "white", size = 0.5),  # Hintergrund der Legende
      legend.text = element_text(face = "plain", size = 15, family = "Arial"),  # Formatierung der Legendentexte
      legend.title = element_text(face = "bold", size = 17, family = "Arial"),  # Formatierung des Legendetitels
      legend.position = "top",  # Position der Legende
      
      panel.border = element_rect(colour = "black", fill = NA, size = 0.8),  # Rahmen um die Panels
      panel.spacing = unit(0.45, "cm"),  # Abstand zwischen den Panels
      
      plot.margin = margin(2, 2, 2, 2),  # Raender des gesamten Plots
      
      strip.background = element_rect(fill = "grey90", colour = "black", size = 0.8),  # Hintergrund der Facetten-Titel
      strip.clip = "on",  # Zuschneiden der Facetten-Titel
      strip.text = element_text(face = "bold", size = 15)  # Formatierung der Facetten-Titel
    )
  
  # Anpassen der y-Achsen-Grenzen basierend auf der gewaehlten Methode
  if (method_filter == "Multiple Regression") {
    plot <- plot + scale_y_continuous(limits = c(0, 400))
  } else if (method_filter == "Lasso Regression") {
    plot <- plot + scale_y_continuous(limits = c(0, 850))
  }
  
  # Anzeigen des Plots in der R-Konsole
  print(plot)
  
  # Speichern des erstellten Plots als JPEG-Datei im angegebenen Ausgabeordner
  ggsave(filename = file.path(output_dir, 
                              paste0("Histogramm_", method_filter, "_Fall_", case_filter, ".jpeg")),
         plot = plot, 
         width = 9.5, 
         height = 10.5,
         dpi = 600)
}

# ============================================================================ #
# ERSTELLEN UND SPEICHERN DER HISTOGRAMME ------------------------------------
# ============================================================================ #

# Definieren aller verfuegbaren Faelle und Methoden fuer die Schleife
all_cases <- names(true_values_all_cases)  # "A", "B", "C1", "C2"
all_methods <- c("Multiple Regression", "Lasso Regression")  # Methoden zur Analyse

# Durchlaufen aller Kombinationen von Faellen und Methoden, um Histogramme zu erstellen
for (current_case in all_cases) {
  for (current_method in all_methods) {
    create_histogram(
      data_case = data_with_true,       # Datensatz mit geschaetzten und wahren Werten
      method_filter = current_method,   # Aktuelle Methode (Multiple Regression oder Lasso Regression)
      case_filter = current_case,       # Aktueller Fall (A, B, C1, C2)
      output_dir = "9.1_Histograms"     # Verzeichnis zum Speichern der Histogramme
    )
  }
}

# ============================================================================ #
# ENDE DES SKRIPTS -----------------------------------------------------------
# ============================================================================ #