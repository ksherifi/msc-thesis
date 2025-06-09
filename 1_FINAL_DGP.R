# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT --------------------------------------
# ============================================================================ #
# Dieses R-Skript enthaelt drei Funktionen (dgp.case.a, dgp.case.b, dgp.case.c), 
# die jeweils Datensaetze fuer unterschiedliche Szenarien (Fall A, Fall B und Fall C) 
# generieren. Diese Daten dienen zur Simulation einer abhaengigen Variable 
# (YTRUE) basierend auf kategorialen unabhaengigen Variablen (A, B, C) und 
# verschiedenen Regressionsparametern (betas).
# 
# Jede Funktion beinhaltet folgende Schritte:
# 1. Zufaellige Ziehung einer kategorialen Variable mit drei Kategorien (A, B, C).
# 2. Dummykodierung der Kategorien, wobei eine Kategorie (A) als Referenz 
#    festgelegt wird.
# 3. Generierung eines zufaelligen Fehlerterms.
# 4. Simulation der abhaengigen Variable (YTRUE) basierend auf den 
#    Dummykodierungen, den Regressionsparametern (betas) und dem Fehlerterm.
# 5. Rueckgabe eines DataFrames mit den generierten Daten.
# 
# Das Skript stellt den ersten Teil einer dreiteiligen Simulationsstudie dar.
# ============================================================================ #

# ============================================================================ #
# DATENGENERIERUNGSFUNKTIONEN DEFINIEREN ------------------------------------- #
# ============================================================================ #

# ============================================================================ #
# FUNKTION FUER FALL A -------------------------------------------------------- #
# ============================================================================ #
# Simuliert Daten fuer eine multiple lineare Regression mit einer kategorialen 
# Variablen, die drei Auspraegungen (A, B, C) umfasst, wobei A als Referenzkategorie dient. 
# Die abhaengige Variable YTRUE wird auf Grundlage von Dummy-Variablen und den 
# Regressionsparametern (betas) sowie einem Fehlerterm erzeugt.
# 
# Parameter:
# - n.total: Anzahl der zu generierenden Beobachtungen.
# - betas: Vektor der Regressionsparameter
# 
# Rueckgabe:
# - DataFrame mit den generierten Dummy-Variablen und der abhaengigen Variable YTRUE.
# ============================================================================ #
dgp.case.a <- function(n.total, betas){
  
  # Ziehung der kategorialen Variable mit den Kategorien A, B und C
  cat.var.a <- as.factor(sample(x = LETTERS[1:3], size = n.total, 
                                replace = TRUE, prob = c(0.33, 0.33, 0.33)))
  
  # Dummykodierung: A wird als Referenz festgelegt, B und C werden als Dummy-Variablen kodiert
  d.code.a <- model.matrix(~ cat.var.a - 1)
  
  # Generierung des Fehlerterms (normalverteilte Zufallszahlen)
  err.a <- rnorm(n = n.total, mean = 0, sd = 1) 
  
  # Berechnung der abhaengigen Variable YTRUE basierend auf den Regressionsparametern 
  # und dem Fehlerterm
  y.true.a <- betas[4] + betas[3] * d.code.a[, "cat.var.aB"] + 
    betas[3] * d.code.a[, "cat.var.aC"] + err.a
  
  # Erstellung eines DataFrames mit den generierten Variablen und Rueckgabe der Daten
  data.dgp.case.a <- data.frame(CASE.ID = "A", A = d.code.a[, "cat.var.aA"],
                                B = d.code.a[, "cat.var.aB"], C = d.code.a[, "cat.var.aC"], 
                                ERR = err.a, YTRUE = y.true.a)
  
  return(data.dgp.case.a)
}

# Beispielhafte Nutzung der Funktion fuer Fall A
# set.seed(123)
# data.dgp.case.a <- dgp.case.a(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.dgp.case.a, file = "1.1_RESULTS_DGP_A.RData")

# ============================================================================ #
# FUNKTION FUER FALL B -------------------------------------------------------- #
# ============================================================================ #
# Aehnlich wie Fall A, jedoch werden unterschiedliche Effektstaerken fuer die 
# Dummy-Variablen B und C verwendet.
# 
# Parameter:
# - n.total: Anzahl der zu generierenden Beobachtungen.
# - betas: Vektor der Regressionsparameter
# 
# Rueckgabe:
# - DataFrame mit den generierten Dummy-Variablen und der abhaengigen Variable YTRUE.
# ============================================================================ #
dgp.case.b <- function(n.total, betas){
  
  # Ziehung der kategorialen Variable mit den Kategorien A, B und C
  cat.var.b <- as.factor(sample(x = LETTERS[1:3], size = n.total, 
                                replace = TRUE, prob = c(0.33, 0.33, 0.33)))
  
  # Dummykodierung: A wird als Referenz festgelegt
  d.code.b <- model.matrix(~ cat.var.b - 1)
  
  # Generierung des Fehlerterms (normalverteilte Zufallszahlen)
  err.b <- rnorm(n = n.total, mean = 0, sd = 1)
  
  # Berechnung der abhaengigen Variable YTRUE basierend auf den Regressionsparametern 
  # und dem Fehlerterm
  y.true.b <- betas[4] + betas[4] * d.code.b[, "cat.var.bB"] + 
    betas[2] * d.code.b[, "cat.var.bC"] + err.b
  
  # Erstellung eines DataFrames mit den generierten Variablen und Rueckgabe der Daten
  data.dgp.case.b <- data.frame(CASE.ID = "B", A = d.code.b[, "cat.var.bA"], 
                                B = d.code.b[, "cat.var.bB"], C = d.code.b[, "cat.var.bC"], 
                                ERR = err.b, YTRUE = y.true.b)
  
  return(data.dgp.case.b)
}

# Beispielhafte Nutzung der Funktion fuer Fall B
# set.seed(456)
# data.dgp.case.b <- dgp.case.b(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.dgp.case.b, file = "1.2_RESULTS_DGP_B.RData")

# ============================================================================ #
# FUNKTION FUER FALL C -------------------------------------------------------- #
# ============================================================================ #
# Erweiterte Simulation mit zwei unterschiedlichen Sub-Faellen (C1 und C2), 
# bei denen die abhaengige Variable YTRUE unterschiedlich berechnet wird.
# 
# Parameter:
# - n.total: Anzahl der zu generierenden Beobachtungen.
# - betas: Vektor der Regressionsparameter
# 
# Rueckgabe:
# - DataFrame mit den generierten Dummy-Variablen und der abhaengigen Variable YTRUE fuer beide Sub-Faelle.
# ============================================================================ #
dgp.case.c <- function(n.total, betas){
  
  # Ziehung der kategorialen Variable mit den Kategorien A, B und C
  cat.var.c <- as.factor(sample(x = LETTERS[1:3], size = n.total, 
                                replace = TRUE, prob = c(0.33, 0.33, 0.33)))
  
  # Dummykodierung: A wird als Referenz festgelegt
  d.code.c <- model.matrix(~ cat.var.c - 1)
  
  # Generierung des Fehlerterms (normalverteilte Zufallszahlen)
  err.c <- rnorm(n = n.total, mean = 0, sd = 1)
  
  # Szenario 1: Berechnung der abhaengigen Variable YTRUE fuer Fall C1 basierend 
  # auf den Regressionsparametern und dem Fehlerterm
  y.true.cone <- betas[5] + betas[3] * d.code.c[, "cat.var.cB"] + 
    betas[2] * d.code.c[, "cat.var.cC"] + err.c
  
  # Szenario 2: Berechnung der abhaengigen Variable YTRUE fuer Fall C2 basierend 
  # auf den Regressionsparametern 
  y.true.ctwo <- betas[5] + betas[3] * d.code.c[, "cat.var.cB"] + 
    betas[1] * d.code.c[, "cat.var.cC"] + err.c
  
  # Erstellung von DataFrames mit den generierten Variablen fuer beide Szenarien
  data.cone <- data.frame(CASE.ID = "C1", A = d.code.c[, "cat.var.cA"], 
                          B = d.code.c[, "cat.var.cB"], C = d.code.c[, "cat.var.cC"], 
                          ERR = err.c, YTRUE = y.true.cone)
  
  data.ctwo <- data.frame(CASE.ID = "C2", A = d.code.c[, "cat.var.cA"], 
                          B = d.code.c[, "cat.var.cB"], C = d.code.c[, "cat.var.cC"], 
                          ERR = err.c, YTRUE = y.true.ctwo)
  
  # Erstellung eines gemeinsamen DataFrames und Rueckgabe der Daten
  data.dgp.case.c <- rbind(data.cone, data.ctwo)
  
  return(data.dgp.case.c)
}

# Beispielhafte Nutzung der Funktion fuer Fall C
# set.seed(789)
# data.dgp.case.c <- dgp.case.c(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.dgp.case.c, file = "1.3_RESULTS_DGP_C.RData")

# ============================================================================ #
# ENDE DES SKRIPTS ----------------------------------------------------------- #
# ============================================================================ #