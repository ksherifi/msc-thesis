# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT --------------------------------------
# ============================================================================ #
# Dieses Skript erstellt Visualisierungen fuer drei verschiedene Kategorien (A, B, C)
# mit Dummy-Variablen und speichert die Plots als Bilddateien (JPEG).
# 
# Die Visualisierungen zeigen Punkte fuer jede Kategorie in verschiedenen Panels.
# Es wird sichergestellt, dass alle Panels die gleiche Groesse haben, und dass die 
# Beschriftungen auf den Achsen korrekt und konsistent sind.
# 
# Die Legenden und Titel aendern sich entsprechend der Referenzkategorie (A, B, C).
# 
# Die Ausgabe erfolgt als JPEG-Dateien fuer jede Referezkatgorie und wird im Verzeichnis 
# "10.1_Figures_Methodology" gespeichert.
# ============================================================================ #

# ============================================================================ #
# BENOETIGTE PACKAGES INSTALLIEREN UND LADEN ---------------------------------
# ============================================================================ #

# install.packages("cowplot")
# install.packages("ggplot2")

library(cowplot)  # Erstellen und Anordnen von mehreren Plots in einem Raster
library(ggplot2)  # Erstellen von Grafiken

# ============================================================================ #
# ABBILDUNGEN MIT REFERENZKATEGORIE A ----------------------------------------
# ============================================================================ #
# ============================================================================ #
# FUNKTION ZUR ERSTELLUNG EINZELNER PLOTS -----------------  *+*+*KI-HELP*+*+*
# ============================================================================ #
# Die Funktion `create_plot` wird verwendet, um einzelne Plots fuer jede Kombination
# von Werten zu erstellen. Diese Funktion generiert einen einfachen Scatterplot
# fuer die drei Kategorien A, B und C und formatiert den Plot entsprechend.
# Die Funktion erhaelt folgende Parameter:
# - A_value, B_value, C_value: Werte fuer die drei Kategorien A, B und C
# - plot_title: Titel fuer den Plot (z.B. "1.1", "1.2", etc.)
# - y_axis: Logischer Wert, der bestimmt, ob die y-Achse angezeigt werden soll

create_plot <- function(A_value, B_value, C_value, plot_title, y_axis) {
  # Erstellung eines Dataframes mit den Kategorien A, B, C und deren Werten
  data <- data.frame(
    Kategorie = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    Wert = c(A_value, B_value, C_value)
  )
  
  # Erstellung des Scatterplots
  p <- ggplot(data, aes(x = Kategorie, y = Wert, color = Kategorie)) +
    geom_point(size = 5) +  # Zeichnen der Punkte mit einer Groesse von 5
    scale_color_manual(values = c("A" = "red", "B" = "blue", "C" = "green")) +  # Manuelle Zuweisung der Farben zu den Kategorien
    theme_minimal() +  # Einfaches Layout ohne unnoetige Elemente
    theme(legend.position = "none",  # Die Legende wird in jedem Einzelplot deaktiviert
          axis.title.x = element_blank(), # Entfernen des Titels der x-Achse
          axis.title.y = element_blank(), # Entfernen des Titels der y-Achse
          axis.text.x = element_blank(), # Entfernen der Beschriftungen auf der x-Achse
          axis.text.y = element_text(face = "bold", size = 15),  # Formatierung der y-Achsen-Beschriftung
          plot.title = element_text(hjust = 0.5, face = "bold", size = 15),  # Formatierung und Zentrierung des Titels
          panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # Schwarzer Rahmen um jeden Panel
          panel.spacing = unit(0.6, "cm")) +  # Abstand zwischen den Panels
    labs(title = plot_title) +  # Setzen des Titels fuer den Plot
    scale_y_continuous(limits = c(0, 10), breaks = c(0, 5, 10))  # Begrenzung der y-Achse zwischen 0 und 10
  
  # Falls der Plot nicht in der ersten Spalte ist, werden die y-Achsen-Beschriftungen entfernt
  if (!y_axis) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  # Entfernen der y-Achsenticks und -beschriftungen
  }
  
  return(p)  # Rueckgabe des erstellten Plots
}

# ============================================================================ #
# DATEN FUER DIE KATEGORIEN DEFINIEREN ---------------------------------------
# ============================================================================ #
# Definition der Werte fuer jede Kombination von Kategorien A, B und C
# Die Werte werden in einer Liste `value_sets` gespeichert, wobei jede Liste drei Werte
# enthaelt, die jeweils fuer A, B und C stehen.

value_sets <- list(
  c(10, 10, 10), c(10, 10, 5), c(10, 10, 0),
  c(10, 5, 10), c(10, 5, 5), c(10, 5, 0),
  c(10, 0, 10), c(10, 0, 5), c(10, 0, 0),
  c(5, 10, 10), c(5, 10, 5), c(5, 10, 0),
  c(5, 5, 10), c(5, 5, 5), c(5, 5, 0),
  c(5, 0, 10), c(5, 0, 5), c(5, 0, 0),
  c(0, 10, 10), c(0, 10, 5), c(0, 10, 0),
  c(0, 5, 10), c(0, 5, 5), c(0, 5, 0),
  c(0, 0, 10), c(0, 0, 5), c(0, 0, 0)
)

# Titel fuer jeden Plot, die in der Visualisierung angezeigt werden
plot_titles_A <- paste0("1.", 1:27)
plot_titles_B <- paste0("2.", 1:27)
plot_titles_C <- paste0("3.", 1:27)

# ============================================================================ #
# VORAUSSETZUNGEN ZUR GENERIERUNG DER PLOTS FUER REFERENZ A  *+*+*KI-HELP*+*+*
# ============================================================================ #
# Erstellung der 27 Plots fuer die Visualisierung von Referenz A.
plots_A <- list()  # Leere Liste zur Speicherung der erstellten Plots
for (i in 1:27) {
  # Wenn der Plot in der ersten Spalte (jede 9. Position) ist, wird die y-Achse angezeigt
  y_axis <- (i %% 9 == 1)
  
  # Erstellen des Plots und Speichern in der Liste
  plots_A[[i]] <- create_plot(value_sets[[i]][1], value_sets[[i]][2], value_sets[[i]][3], plot_titles_A[i], y_axis)
}

# ============================================================================ #
# LEGENDENPLOT FUER REFERENZ A ----------------------------  *+*+*KI-HELP*+*+*
# ============================================================================ #
# Erstellung eines Dummy-Datensatzes fuer die Legende, um die Farben der Kategorien A, B und C
# zu zeigen. Der Dummy-Datensatz wird in der Legende verwendet, ohne dass zusaetzliche Punkte
# auf der Plotflaeche gezeichnet werden.

legend_plot_A <- ggplot(data.frame(Kategorie = factor(c("A", "B", "C"), levels = c("A", "B", "C"))),
                        aes(x = Kategorie, y = Kategorie, color = Kategorie)) +
  geom_point(size = 0) +  # Keine Punkte auf der Plotflaeche zeichnen
  scale_color_manual(values = c("A" = "red", "B" = "blue", "C" = "green"),
                     labels = c(expression(bolditalic(A) * " (Referenzkategorie)"),
                                expression(bolditalic(B) * " (Dummy-Variable 1)"),
                                expression(bolditalic(C) * " (Dummy-Variable 2)"))) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  # Nur Legende anzeigen, ohne Punkte
  theme_void() +  # Leeres Layout fuer die Legende (ohne Achsen oder Titel)
  theme(legend.position = "top",  # Legende oben platzieren
        legend.title = element_text(face = "bold", size = 17),  # Formatierung des Legendentitels
        legend.text = element_text(size = 15, face = "plain"),  # Formatierung der Legendentexte
        legend.box.background = element_rect(color = "black", fill = "white", size = 0.5),  # Rahmen um die Legende
        legend.margin = margin(10, 10, 10, 10),  # Abstand um die Legende herum
        plot.margin = margin(10, 10, 10, 10)) +  # Rand um den gesamten Legendenplot
  labs(color = "Kategorien")  # Legendentitel

# Anordnung der Plots und Speichern des Plots fuer Referenz A
final_plot_A <- plot_grid(plotlist = plots_A, ncol = 9, align = "hv", rel_widths = rep(1, 9))
combined_plot_A <- plot_grid(legend_plot_A, final_plot_A, ncol = 1, rel_heights = c(0.15, 1))

# Speichern des Plots fuer Referenz A
ggsave(filename = file.path("10.1_Figures_Methodology", "METHODIK_Visualisierungen_ReferenzA.jpeg"),
       plot = combined_plot_A,
       width = 9.5, 
       height = 10.5, 
       dpi = 600)

# ============================================================================ #
# ABBILDUNGEN MIT REFERENZKATEGORIE B ----------------------------------------
# ============================================================================ #
# Erstellung der Plots fuer Referenz B mit neuen Titeln und einer angepassten Legende

plots_B <- list()  # Leere Liste zur Speicherung der erstellten Plots fuer Referenz B
for (i in 1:27) {
  y_axis <- (i %% 9 == 1)  # Y-Achse nur fuer die erste Spalte
  plots_B[[i]] <- create_plot(value_sets[[i]][1], value_sets[[i]][2], value_sets[[i]][3], plot_titles_B[i], y_axis)
}

# ============================================================================ #
# LEGENDENPLOT FUER REFERENZ B -----------------------------------------------
# ============================================================================ #

# Anpassung der Farben: B ist die Referenzkategorie (blau), A (rot), C (gruen)
legend_plot_B <- ggplot(data.frame(Kategorie = factor(c("B", "A", "C"), levels = c("B", "A", "C"))),  # Reihenfolge B, A, C
                        aes(x = Kategorie, y = Kategorie, color = Kategorie)) +
  geom_point(size = 0) +  # Keine Punkte auf der Plotflaeche zeichnen
  scale_color_manual(values = c("B" = "blue", "A" = "red", "C" = "green"),  # Farben entsprechend aendern
                     labels = c(expression(bolditalic(B) * " (Referenzkategorie)"),
                                expression(bolditalic(A) * " (Dummy-Variable 1)"),
                                expression(bolditalic(C) * " (Dummy-Variable 2)"))) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  # Nur Legende anzeigen, ohne Punkte
  theme_void() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 17),
        legend.text = element_text(size = 15, face = "plain"),
        legend.box.background = element_rect(color = "black", fill = "white", size = 0.5),
        legend.margin = margin(10, 10, 10, 10),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(color = "Kategorien")

# Anordnung der Plots und Speichern des Plots fuer Referenz B
final_plot_B <- plot_grid(plotlist = plots_B, ncol = 9, align = "hv", rel_widths = rep(1, 9))
combined_plot_B <- plot_grid(legend_plot_B, final_plot_B, ncol = 1, rel_heights = c(0.15, 1))

ggsave(filename = file.path("10.1_Figures_Methodology", "METHODIK_Visualisierungen_ReferenzB.jpeg"),
       plot = combined_plot_B,
       width = 9.5, 
       height = 10.5, 
       dpi = 600)

# ============================================================================ #
# ABBILDUNGEN MIT REFERENZKATEGORIE C ----------------------------------------
# ============================================================================ #
# Erstellung der Plots fuer Referenz C mit neuen Titeln und einer angepassten Legende

plots_C <- list()  # Leere Liste zur Speicherung der erstellten Plots fuer Referenz C
for (i in 1:27) {
  y_axis <- (i %% 9 == 1)  # Y-Achse nur fuer die erste Spalte
  plots_C[[i]] <- create_plot(value_sets[[i]][1], value_sets[[i]][2], value_sets[[i]][3], plot_titles_C[i], y_axis)
}

# ============================================================================ #
# LEGENDENPLOT FUER REFERENZ C -----------------------------------------------
# ============================================================================ #
# Anpassung der Farben: C ist die Referenzkategorie (gruen), A (rot), B (blau)
legend_plot_C <- ggplot(data.frame(Kategorie = factor(c("C", "A", "B"), levels = c("C", "A", "B"))),  # Reihenfolge C, A, B
                        aes(x = Kategorie, y = Kategorie, color = Kategorie)) +
  geom_point(size = 0) +  # Keine Punkte auf der Plotflaeche zeichnen
  scale_color_manual(values = c("C" = "green", "A" = "red", "B" = "blue"),  # Farben entsprechend aendern
                     labels = c(expression(bolditalic(C) * " (Referenzkategorie)"),
                                expression(bolditalic(A) * " (Dummy-Variable 1)"),
                                expression(bolditalic(B) * " (Dummy-Variable 2)"))) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  # Nur Legende anzeigen, ohne Punkte
  theme_void() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 17),
        legend.text = element_text(size = 15, face = "plain"),
        legend.box.background = element_rect(color = "black", fill = "white", size = 0.5),
        legend.margin = margin(10, 10, 10, 10),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(color = "Kategorien")

# Anordnung der Plots und Speichern des Plots fuer Referenz C
final_plot_C <- plot_grid(plotlist = plots_C, ncol = 9, align = "hv", rel_widths = rep(1, 9))
combined_plot_C <- plot_grid(legend_plot_C, final_plot_C, ncol = 1, rel_heights = c(0.15, 1))

ggsave(filename = file.path("10.1_Figures_Methodology", "METHODIK_Visualisierungen_ReferenzC.jpeg"),
       plot = combined_plot_C,
       width = 9.5, 
       height = 10.5, 
       dpi = 600)

# ============================================================================ #
# ENDE DES SKRIPTS -----------------------------------------------------------
# ============================================================================ #