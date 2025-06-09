# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT --------------------------------------
# ============================================================================ #
# Dieses Skript erstellt Boxplots zur Visualisierung der Erwartungstreue der Punktschaetzer
# fuer beide Regressionsverfahren in den Faellen A und B sowie den Sub-Faellen C1 und C2.
# Dabei erfolgt die Visualisierung der Erwartungstreue  separat fuer die Parameter
# Achsenabschnitt, Effektstaerke D1 und Effektstaerke D2.
#
# Die Simulationsergebnisse werden zuerst fuer beide Regressionsverfahren
# geladen und in einem kombinierten Datensatz zusammengefuehrt. Anschliessend werden die
# Daten in ein langes Format transformiert, um die Parameterwerte uebersichtlich darzustellen.
# Fuer jede Referenzkategorie und jeden Fall werden die wahren Parameterwerte definiert und
# als horizontale Referenzlinien in die Boxplots eingefuegt.
#
# Die Boxplots zeigen die Abweichungen der geschaetzten Werte vom wahren Wert und erlauben
# den direkten Vergleich der Methoden fuer verschiedene Stichprobengroessen (n=50 und n=500).
# Die y-Achsen fuer jede Facette sind individuell auf den Bereich der wahren Werte ±1.5 skaliert,
# um die Ergebnisse besser vergleichbar zu machen. Die Ausgabe erfolgt als JPEG-Dateien fuer
# jeden Fall und wird im Verzeichnis "8.1_Boxplots" gespeichert.
# ============================================================================ #

# ============================================================================ #
# BENOETIGTE PAKETE INSTALLIEREN UND LADEN -----------------------------------
# ============================================================================ #

# install.packages("ggplot2")   # Zum Erstellen und Gestalten der Boxplots
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

# Auswahl der relevanten Spalten fuer Boxplots zur multiplen Regression
data_boxplots_linreg <- data_linreg[, c("n.total", "CASE.ID", "REF.ID", 
                                        "EST.INT", "EST.SLOPE.ONE", "EST.SLOPE.TWO", "METHOD")]

# Auswahl der relevanten Spalten fuer Boxplots zur Lasso Regression
data_boxplots_lasso <- data_lasso[, c("n.total", "CASE.ID", "REF.ID", 
                                      "EST.INT", "EST.SLOPE.ONE", "EST.SLOPE.TWO", "METHOD")]

# ============================================================================ #
# KOMBINIEREN BEIDER DATENSAETZE IN EINEM DATAFRAME --------------------------
# ============================================================================ #

# Zusammenfuehren der beiden Datensaetze in einen einzigen DataFrame
combined_data <- bind_rows(data_boxplots_linreg, data_boxplots_lasso)

# ============================================================================ #
# ANPASSEN DER 'METHOD'-SPALTE FUER KONSISTENTE BEZEICHNUNGEN ----------------
# ============================================================================ #

# Umkodierung der 'METHOD'-Spalte in Faktoren mit klaren und verstaendlichen Labels
combined_data <- combined_data %>%
  mutate(METHOD = factor(METHOD,
                         levels = c("LinReg", "Lasso"),
                         labels = c("Multiple Regression", "Lasso Regression")))

# ============================================================================ #
# DATEN TRANSFORMIEREN IN LONG-FORMAT --------------------   *+*+*KI-HELP*+*+*
# ============================================================================ #

# Umstrukturierung der kombinierten Daten in ein langes Format, um die Parameter als separate Zeilen zu haben
combined_long <- combined_data %>%
  pivot_longer(
    cols = starts_with("EST"),  # Auswahl der EST-Variablen (EST.INT, EST.SLOPE.ONE, EST.SLOPE.TWO)
    names_to = "Parameter",     # Name der neuen Spalte fuer Parameter
    values_to = "value"         # Name der neuen Spalte fuer Werte
  ) %>%
  mutate(
    # Umkodierung der 'Parameter'-Spalte in klar benannte Faktoren mit verstaendlichen Labels
    Parameter = factor(Parameter, 
                       levels = c("EST.INT", "EST.SLOPE.ONE", "EST.SLOPE.TWO"),
                       labels = c("Achsenabschnitt", "Effektstärke D1", "Effektstärke D2")),
    # Umkodierung der 'REF.ID'-Spalte in klar benannte Faktoren mit verstaendlichen Labels
    REF.ID = factor(REF.ID, 
                    levels = c("REFA", "REFB", "REFC"),
                    labels = c("Referenz A", "Referenz B", "Referenz C")),
    # Umkodierung von 'n.total' in einen Faktor zur Kategorisierung nach Stichprobengroesse
    n.total = factor(n.total)
  )

# ============================================================================ #
# DEFINITION DER WAHREN WERTE FUER ALLE FAELLE -------------------------------
# ============================================================================ #

# Wahre Werte fuer jeden Fall und jede Referenzkategorie
true_values_all <- list(
  
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

# ============================================================================ #
# DEFINIEREN DER FUNKTION ZUR ERSTELLUNG VON BOXPLOTS -----  *+*+*KI-HELP*+*+*
# ============================================================================ #

# Definition einer Funktion zur Erstellung von Boxplots fuer jeden Fall
create_boxplot <- function(case_id, true_values_case, data_combined_long, output_dir = "8.1_Boxplots") {
  
  # Filtern der Daten fuer den aktuellen Fall
  data_case <- data_combined_long %>%
    filter(CASE.ID == case_id)
  
  # Transformieren der wahren Werte in ein langes Format
  true_vals_long <- true_values_case %>%
    enframe(name = "REF.ID", value = "parameters") %>%    # Umwandlung der Liste in ein DataFrame mit REF.ID und parameters
    unnest_wider(parameters) %>%                          # Erweiterung der Parameter-Spalten (Achsenabschnitt, Effektstaerke D1, Effektstaerke D2)
    pivot_longer(
      cols = -REF.ID,                                     # Ausschluss von REF.ID aus der Pivotierung
      names_to = "Parameter",                             # Name der neuen Spalte fuer Parameter
      values_to = "value"                                 # Name der neuen Spalte fuer Werte
    )
  
  # Berechnen der y-Achsen-Grenzen fuer jede Facette basierend auf den wahren Werten ±1.5
  facet_scales <- true_vals_long %>%
    mutate(
      lower = value - 1.5,  # Untere Grenze fuer die y-Achse
      upper = value + 1.5   # Obere Grenze fuer die y-Achse
    ) %>%
    group_by(REF.ID, Parameter) %>%                          
    summarize(lower = first(lower), upper = first(upper), .groups = "drop") %>%       # Zusammenfassung der Grenzen pro Facette
    mutate(scale = map2(lower, upper, ~ scale_y_continuous(limits = c(.x, .y)))) %>%  # Erstellung individueller Skalen
    pull(scale)                                                                       # Extraktion der Skalen als Liste
  
  # Definieren der Farbpalette fuer die 'METHOD'-Kategorie
  colors_method <- scale_fill_manual(
    values = c("Multiple Regression" = "red3",
               "Lasso Regression" = "dodgerblue"),
    name = "Methode"
  )
  
  # Aufbau der ggplot2-Grafik
  plot <- ggplot(data = data_case, 
                 mapping = aes(x = n.total, y = value, fill = METHOD)) +
    
    # Hinzufuegen der Boxplots
    geom_boxplot(
      position = position_dodge(width = 0.8),  # Positionierung der Boxplots nebeneinander fuer verschiedene Methoden
      outlier.colour = "black",                # Farbe der Ausreisserpunkte
      outlier.shape = 42,                      # Form der Ausreisserpunkte
      outlier.size = 5,                        # Groesse der Ausreisserpunkte
      outlier.stroke = 0.5,                    # Strichstaerke der Ausreisserpunkte
      staplewidth = 1,                         # Breite der "Staples" (Verbindungen der Boxen)
      width = 0.6                              # Breite der Boxplots
    ) +
    
    # Hinzufuegen der horizontalen Linien fuer die wahren Parameterwerte
    geom_hline(
      data = true_vals_long, 
      aes(yintercept = value), 
      color = "darkred",        # Farbe der Linien
      linetype = "solid",       # Linientyp (durchgezogen)
      size = 0.3                # Linienstaerke
    ) +
    
    # Facettierung nach Referenz-ID und Parameter
    facet_wrap2(
      facets = ~ REF.ID + Parameter, 
      nrow = 3,                     # Anzahl der Zeilen in der Facettierung
      scales = "free_y",            # Freie Skalierung der y-Achsen fuer jede Facette
      axes = "all"                  # Anzeige aller Achsen
    ) +
    
    # Anwendung der individuellen Facet-Skalen
    facetted_pos_scales(y = facet_scales) +
    
    # Anwendung der benutzerdefinierten Farbskala
    colors_method +
    
    # Beschriftung der Achsen
    labs(
      x = expression(bold("Stichprobengrösse (" * bolditalic(n) * " = 50, 500)")),  # Dynamischer Achsentitel mit spezifischer Notation
      y = "Abweichung der Punktschätzer vom Wahren Parameterwert"                   # Achsentitel fuer die y-Achse
    ) +
    
    # Anwendung eines minimalistischen Themas mit Anpassungen
    theme_minimal(base_size = 16, base_family = "Arial") +
    theme(
      axis.title.x = element_text(face = "bold", size = 20, margin = margin(t = 25)),  # Formatierung der x-Achsentitel
      axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 25)),  # Formatierung der y-Achsentitel
      
      axis.text.x = element_text(face = "bold", size = 16),  # Formatierung der x-Achsentexte
      axis.text.y = element_text(face = "bold", size = 16),  # Formatierung der y-Achsentexte
      
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
  
  # Anzeigen der Abbildung in der R-Konsole
  print(plot)
  
  # Speichern der Abbildung als JPEG-Datei im angegebenen Ausgabeverzeichnis
  ggsave(filename = file.path(output_dir, paste0("Boxplot_Fall_", case_id, ".jpeg")),
         plot = plot,
         width = 9.5, 
         height = 10.5, 
         dpi = 600)
  
}

# ============================================================================ #
# ERSTELLEN UND SPEICHERN DER BOXPLOTS FUER ALLE FAELLE ----------------------
# ============================================================================ #

# Iterieren ueber alle definierten Faelle und Boxplots erstellen
for(case_id in names(true_values_all)) {
  create_boxplot(
    case_id = case_id, 
    true_values_case = true_values_all[[case_id]],
    data_combined_long = combined_long,
    output_dir = "8.1_Boxplots"  # Verzeichnis zum Speichern der Boxplots
  )
}

# ============================================================================ #
# ENDE DES SKRIPTS ----------------------------------------------------------- 
# ============================================================================ #