# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT
# ============================================================================ #
# 
# Dieses Skript fuehrt mehrfache Simulationsdurchlaeufe fuer beide 
# Regressionsverfahren durch. Es werden sowohl multiple als auch Lasso 
# Regressionsmodelle auf Datensaetze mit verschiedenen Stichprobengroessen und 
# Regressionsparametern angewendet.
# 
# Ziel dieses Skripts ist es, alle Daten zu simulieren, die fuer die statistische 
# Analyse der Modellstruktur und der Guetekriterien benoetigt werden. Dabei werden 
# die generierten Datensaetze anhand vorgegebener Parameter (Regressionsparameter 
# und Stichprobengroessen) simuliert.
#
# Zwischenergebnisse werden in regelmaessigen Abstaenden gespeichert, um den Fortschritt 
# der Simulation zu sichern und Datenverluste bei laengeren Berechnungen zu vermeiden.
# 
# Die Ergebnisse werden am Ende in verschiedenen .RData-Dateien abgelegt, die sowohl 
# die einzelnen Ergebnisse der multiplen und Lasso Regressionsmodelle als auch die 
# aggregierten Endergebnisse enthalten.
# 
# ---------------------------------------------------------------------------- #
# ============================================================================ #

# ============================================================================ #
# R-SKRIPTE MIT ONESIM-FUNKTIONEN ZUR MULTIPLEN UND LASSO REGRESSION LADEN ---
# ============================================================================ #

source("2_FINAL_ONESIM_LINREG.R")
source("3_FINAL_ONESIM_LASSOREG.R")

# ============================================================================ #
# BENOETIGTE PAKETE INSTALLIEREN UND LADEN -----------------------------------
# ============================================================================ #

# Installieren der benoetigten Pakete (falls nicht bereits installiert)
# install.packages("glmnet")
# install.packages("boot")

# Laden der Pakete
library(glmnet)  # Fuer Lasso Regression
library(boot)  # Fuer Bootstrapping

# ============================================================================ #
# FUNKTION FÜR VOLLSTÄNDIGE SIMULATION ---------------------------------------
# ============================================================================ #
# 
# Funktion zur Durchfuehrung der Simulation ueber alle Faelle und Bedingungen hinweg.
# Parameter:
# - num.datasets: Anzahl der Datensaetze, die erstellt werden sollen (= Anzahl Simulationsdurchlaeufe)
# - n.total: Vektor der Stichprobengroessen, die simuliert werden sollen
# - betas: Vektor der wahren Regressionsparameter, die in der Simulation verwendet werden
# - save.interval: Wie haeufig die Ergebnisse zwischengespeichert werden sollen
# - linreg.file: Dateiname fuer die Speicherung der Ergebnisse zur multiplen Regression
# - lasso.file: Dateiname fuer die Speicherung der Ergebnisse zur Lasso Regression

sim.stud.all.cases <- function(num.datasets, 
                               n.total, 
                               betas,
                               save.interval = 1, 
                               linreg.file = "4.1_RESULTS_SIMSTUD_LINREG.RData", 
                               lasso.file = "4.2_RESULTS_SIMSTUD_LASSO.RData"){
  
  # Erstellen eines Datenrahmens, der alle Kombinationen von Iterationen (num.datasets) 
  # und Stichprobengroessen (n.total) enthaelt. Jeder Durchlauf der Simulation wird 
  # auf einer dieser Kombinationen basieren.
  df.combinations <- expand.grid(i = 1:num.datasets, n.total = n.total)
  
  # Leere Datenrahmen fuer die Endergebnisse der multiplen und Lasso Regression.
  # Hier werden die Ergebnisse aus allen Simulationen gesammelt.
  final.results.linreg <- data.frame()
  final.results.lasso <- data.frame()
  
  # Startzeitpunkt der gesamten Simulation wird aufgezeichnet, um die Gesamtdauer zu berechnen.
  start.time <- Sys.time()
  print(paste("Simulation gestartet um:", start.time))
  
  # Schleife ueber alle Kombinationen aus Iterationen und Stichprobengroessen.
  for (i in 1:nrow(df.combinations)){
    iteration.start.time <- Sys.time()  # Startzeitpunkt der aktuellen Iteration
    
    # Durchfuehrung von einem Simulationsdurchlauf fuer die multiple Regression 
    # fuer die aktuelle Stichprobengroesse.
    interim.results.linreg <- one.sim.linreg.all.cases(n.total = df.combinations[i,]$n.total,
                                                       betas = betas)
    
    # Durchfuehrung von einem Simulationsdurchlauf fuer die Lasso Regression 
    # fuer die aktuelle Stichprobengroesse.
    interim.results.lasso <- one.sim.lasso.all.cases(n.total = df.combinations[i,]$n.total,
                                                     betas = betas)
    
    # Anzeige, welche Iteration abgeschlossen wurde, um den Fortschritt zu ueberwachen.
    print(paste("Iteration:", i, "abgeschlossen"))
    
    # Die Ergebnisse der aktuellen Simulation werden mit den bisherigen Ergebnissen 
    # zusammengefuehrt. Die Funktion 'rbind' fuegt den neuen Datensatz (eine Zeile) 
    # an den bestehenden Datenrahmen an.
    final.results.linreg <- rbind(final.results.linreg, cbind(df.combinations[i,], 
                                                              interim.results.linreg))
    final.results.lasso <- rbind(final.results.lasso, cbind(df.combinations[i,], 
                                                            interim.results.lasso))
    
    # Zwischenspeicherung der Ergebnisse nach jeder 'save.interval' Iteration.
    # Die Modulo-Bedingung 'i %% save.interval == 0' prueft, ob die aktuelle Iteration 
    # ein Vielfaches von 'save.interval' ist. Falls ja, werden die Ergebnisse gespeichert.
    if (i %% save.interval == 0) {
      # Speichern der Ergebnisse der multiplen Regression in der angegebenen Datei.
      save(final.results.linreg, file = linreg.file)
      # Speichern der Ergebnisse der Lasso Regression in der angegebenen Datei.
      save(final.results.lasso, file = lasso.file)
      
      # Anzeige, dass die Zwischenspeicherung erfolgt ist.
      print(paste("Zwischenspeicherung nach", i, "Iterationen abgeschlossen"))
      
    }
    
    # Endzeitpunkt der aktuellen Iteration wird aufgezeichnet, um die Dauer der Iteration 
    # zu berechnen und auszugeben.
    iteration.end.time <- Sys.time()
    print(paste("Dauer der Iteration", i, ":", iteration.end.time - iteration.start.time))
  }
  
  # Endzeitpunkt der gesamten Simulation wird aufgezeichnet.
  end.time <- Sys.time()
  print(paste("Simulation beendet um:", end.time))
  
  # Die Gesamtdauer der Simulation wird berechnet und ausgegeben.
  print(paste("Gesamtdauer der Simulation:", end.time - start.time))
  
  # Rueckgabe der Endergebnisse der linearen und Lasso Regression als Liste. 
  # Dies wird hauptsaechlich verwendet, falls ein Fehler auftritt oder wenn eine 
  # explizite Rueckgabe benoetigt wird.
  return(list(linreg.results = final.results.linreg, lasso.results = final.results.lasso))
}

# ============================================================================ #
# START DER SIMULATION -------------------------------------------------------
# ============================================================================ #
# set.seed: Stellt sicher, dass die Ergebnisse reproduzierbar sind, indem ein 
# Startwert fuer den Zufallszahlengenerator gesetzt wird.
set.seed(123)

# Aufruf der Simulationsfunktion mit folgenden Parametern:
# - 1000 Simulationen
# - Stichprobengroessen 50 und 500
# - Beta-Werte (-10, -5, 0, 5, 10)
# - Zwischenspeicherung nach jeder zweiten Iteration
data.sim.studs.all.cases <- sim.stud.all.cases(num.datasets = 1000,
                                               n.total = c(50, 500),
                                               betas = c(-10, -5, 0, 5, 10),
                                               save.interval = 2)

# Die Ergebnisse der Simulation für beide Regressionsverfahren werden zusätzlich 
# als .RData-Datei abgespeichert.
save(data.sim.studs.all.cases, file = "4.3_RESULTS_SIMSTUD_ALL.RData")

# ============================================================================ #
# ENDE DES SKRIPTS -----------------------------------------------------------
# ============================================================================ #