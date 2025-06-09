# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT --------------------------------------
# ============================================================================ #
# Dieses Skript dient zur Schaetzung linearer Regressionsmodelle fuer drei verschiedene 
# Faelle (Fall A, B und C).
#
# Fuer jeden Fall wird ein lineares Regressionsmodell geschaetzt, die 
# Regressionsparameter werden extrahiert, herkoemmliche 95%-Konfidenzintervalle und 
# 95%-Bootstrap-Konfidenzintervalle werden berechnet. Ausserdem wird jeweils geprüft, 
# ob die Konfidenzintervalle den wahren Wert eines Parameters enthalten.
#
# Struktur des Skripts:
# 1. Datenvorbereitung fuer die Regression.
# 2. Durchfuehrung der Multiplen Regression fuer Fall A, B und C.
# 3. Schaetzung der Regressionsmodelle und Berechnung herkömmlicher und Bootstrap-
#    Konfidenzintervalle.
# 4. Ueberpruefung ob wahrer Regressionsparameter in den Konfidenzintervallen liegt
#
# Verwendete Methoden:
# - Lineare Regression mit Dummy-Variablen für kategoriale Praediktoren
# - Bootstrapping zur Berechnung der Konfidenzintervalle
#
# Pakete:
# - 'boot'-Paket für das Bootstrapping
#
# ---------------------------------------------------------------------------- #
# ============================================================================ #

# ============================================================================ #
# R-SKRIPT MIT DGP-FUNKTIONEN LADEN ------------------------------------------
# ============================================================================ #

source("1_FINAL_DGP.R") # Laedt die DGP-Funktionen

# ============================================================================ #
# BENOETIGTE PAKETE INSTALLIEREN UND LADEN -----------------------------------
# ============================================================================ #

# Installieren der benoetigten Pakete (falls noch nicht installiert)
# install.packages("boot")

# Laden der erforderlichen Pakete
library(boot)

# =========================== CASE A - LINEARE REGRESSION ================= ####

A.one.sim.linreg <- function(n.total, betas) {
  
  # Abrufen der simulierten Daten fuer Fall A
  data.dgp.case.a <- dgp.case.a(n.total = n.total, betas = betas)
  
  # A.A ----
  # Erstellung der Designmatrix fuer die multiple Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.a)
  y <- data.dgp.case.a$YTRUE
  
  # Multiples lineares Modell anpassen
  A.refA.linreg <- lm(y ~ x)
  
  # Geschaetzte Regresssionsparametern und Konfidenzintervalle berechnen
  A.refA.est.int.linreg <- as.numeric(coef(A.refA.linreg)[1])
  A.refA.ci.int.lower <- confint(A.refA.linreg, level = 0.95)[1, 1]
  A.refA.ci.int.upper <- confint(A.refA.linreg, level = 0.95)[1, 2]
  A.refA.cis.int.linreg <- A.refA.ci.int.lower <= betas[4] & A.refA.ci.int.upper >= betas[4]
  
  A.refA.est.slope.one.linreg <- as.numeric(coef(A.refA.linreg)[2])
  A.refA.ci.slope.one.lower <- confint(A.refA.linreg, level = 0.95)[2, 1]
  A.refA.ci.slope.one.upper <- confint(A.refA.linreg, level = 0.95)[2, 2]
  A.refA.cis.slope.one.linreg <- A.refA.ci.slope.one.lower <= betas[3] & A.refA.ci.slope.one.upper >= betas[3]
  
  A.refA.est.slope.two.linreg <- as.numeric(coef(A.refA.linreg)[3])
  A.refA.ci.slope.two.lower <- confint(A.refA.linreg, level = 0.95)[3, 1]
  A.refA.ci.slope.two.upper <- confint(A.refA.linreg, level = 0.95)[3, 2]
  A.refA.cis.slope.two.linreg <- A.refA.ci.slope.two.lower <= betas[3] & A.refA.ci.slope.two.upper >= betas[3]
  
  # Bootstrapping fuer Konfidenzintervalle
  A.refA.bootstrap <- function(data, indices) {
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ B + C - 1, data = d)
    y <- d$YTRUE
    A.refA.lm.boot <- lm(y ~ x, data = d)
    return(coef(A.refA.lm.boot))
  }
  
  # Bootstrap-Schaetzungen mit 1000 Wiederholungen
  A.refA.boot.results <- boot(data = data.dgp.case.a, statistic = A.refA.bootstrap, R = 1000)
  
  # Bootstrap-Konfidenzintervalle extrahieren
  A.refA.boot.ci.int.lower <- boot.ci(A.refA.boot.results, type = "perc", index = 1)$percent[4]
  A.refA.boot.ci.int.upper <- boot.ci(A.refA.boot.results, type = "perc", index = 1)$percent[5]
  A.refA.boot.cis.int <- A.refA.boot.ci.int.lower <= betas[4] & A.refA.boot.ci.int.upper >= betas[4]
  
  A.refA.boot.ci.slope.one.lower <- boot.ci(A.refA.boot.results, type = "perc", index = 2)$percent[4]
  A.refA.boot.ci.slope.one.upper <- boot.ci(A.refA.boot.results, type = "perc", index = 2)$percent[5]
  A.refA.boot.cis.slope.one <- A.refA.boot.ci.slope.one.lower <= betas[3] & A.refA.boot.ci.slope.one.upper >= betas[3]
  
  A.refA.boot.ci.slope.two.lower <- boot.ci(A.refA.boot.results, type = "perc", index = 3)$percent[4]
  A.refA.boot.ci.slope.two.upper <- boot.ci(A.refA.boot.results, type = "perc", index = 3)$percent[5]
  A.refA.boot.cis.slope.two <- A.refA.boot.ci.slope.two.lower <= betas[3] & A.refA.boot.ci.slope.two.upper >= betas[3]
  
  # Erstellung des Ergebnisdatensatzes mit Schaetzungen und Konfidenzintervallen
  A.refA.data.linreg <- data.frame(
    CASE.ID = "A",
    REF.ID = "REFA",
    TRUE.INT = betas[4],
    EST.INT = A.refA.est.int.linreg,
    CI.INT.LOWER = A.refA.ci.int.lower,
    CI.INT.UPPER = A.refA.ci.int.upper,
    BOOT.CI.INT.LOWER = A.refA.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = A.refA.boot.ci.int.upper,
    CI.INT = A.refA.cis.int.linreg,
    BOOT.CI.INT = A.refA.boot.cis.int,
    TRUE.SLOPE.ONE = betas[3],
    EST.SLOPE.ONE = A.refA.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = A.refA.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = A.refA.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = A.refA.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = A.refA.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = A.refA.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = A.refA.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[3],
    EST.SLOPE.TWO = A.refA.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = A.refA.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = A.refA.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = A.refA.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = A.refA.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = A.refA.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = A.refA.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  # A.B ----
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.a)
  y <- data.dgp.case.a$YTRUE
  
  A.refB.linreg <- lm(y ~ x)
  
  A.refB.est.int.linreg <- as.numeric(coef(A.refB.linreg)[1])
  A.refB.ci.int.lower <- confint(A.refB.linreg, level = 0.95)[1, 1]
  A.refB.ci.int.upper <- confint(A.refB.linreg, level = 0.95)[1, 2]
  A.refB.cis.int.linreg <- A.refB.ci.int.lower <= betas[4] & A.refB.ci.int.upper >= betas[4]
  
  A.refB.est.slope.one.linreg <- as.numeric(coef(A.refB.linreg)[2])
  A.refB.ci.slope.one.lower <- confint(A.refB.linreg, level = 0.95)[2, 1]
  A.refB.ci.slope.one.upper <- confint(A.refB.linreg, level = 0.95)[2, 2]
  A.refB.cis.slope.one.linreg <- A.refB.ci.slope.one.lower <= betas[3] & A.refB.ci.slope.one.upper >= betas[3]
  
  A.refB.est.slope.two.linreg <- as.numeric(coef(A.refB.linreg)[3])
  A.refB.ci.slope.two.lower <- confint(A.refB.linreg, level = 0.95)[3, 1]
  A.refB.ci.slope.two.upper <- confint(A.refB.linreg, level = 0.95)[3, 2]
  A.refB.cis.slope.two.linreg <- A.refB.ci.slope.two.lower <= betas[3] & A.refB.ci.slope.two.upper >= betas[3]
  
  A.refB.bootstrap <- function(data, indices) {
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + C - 1, data = d)
    y <- d$YTRUE
    A.refB.lm.boot <- lm(y ~ x, data = d)
    return(coef(A.refB.lm.boot))
  }
  
  A.refB.boot.results <- boot(data = data.dgp.case.a, statistic = A.refB.bootstrap, R = 1000)
  
  A.refB.boot.ci.int.lower <- boot.ci(A.refB.boot.results, type = "perc", index = 1)$percent[4]
  A.refB.boot.ci.int.upper <- boot.ci(A.refB.boot.results, type = "perc", index = 1)$percent[5]
  A.refB.boot.cis.int <- A.refB.boot.ci.int.lower <= betas[4] & A.refB.boot.ci.int.upper >= betas[4]
  
  A.refB.boot.ci.slope.one.lower <- boot.ci(A.refB.boot.results, type = "perc", index = 2)$percent[4]
  A.refB.boot.ci.slope.one.upper <- boot.ci(A.refB.boot.results, type = "perc", index = 2)$percent[5]
  A.refB.boot.cis.slope.one <- A.refB.boot.ci.slope.one.lower <= betas[3] & A.refB.boot.ci.slope.one.upper >= betas[3]
  
  A.refB.boot.ci.slope.two.lower <- boot.ci(A.refB.boot.results, type = "perc", index = 3)$percent[4]
  A.refB.boot.ci.slope.two.upper <- boot.ci(A.refB.boot.results, type = "perc", index = 3)$percent[5]
  A.refB.boot.cis.slope.two <- A.refB.boot.ci.slope.two.lower <= betas[3] & A.refB.boot.ci.slope.two.upper >= betas[3]
  
  A.refB.data.linreg <- data.frame(
    CASE.ID = "A",
    REF.ID = "REFB",
    TRUE.INT = betas[4],
    EST.INT = A.refB.est.int.linreg,
    CI.INT.LOWER = A.refB.ci.int.lower,
    CI.INT.UPPER = A.refB.ci.int.upper,
    BOOT.CI.INT.LOWER = A.refB.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = A.refB.boot.ci.int.upper,
    CI.INT = A.refB.cis.int.linreg,
    BOOT.CI.INT = A.refB.boot.cis.int,
    TRUE.SLOPE.ONE = betas[3],
    EST.SLOPE.ONE = A.refB.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = A.refB.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = A.refB.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = A.refB.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = A.refB.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = A.refB.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = A.refB.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[3],
    EST.SLOPE.TWO = A.refB.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = A.refB.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = A.refB.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = A.refB.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = A.refB.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = A.refB.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = A.refB.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  # A.C ----
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.a)
  y <- data.dgp.case.a$YTRUE
  
  # Modelanpassung fuer Fall A, Referenz C
  A.refC.linreg <- lm(y ~ x)
  
  # Geschaetzte Regresssionsparametern und KI's fuer Fall A, Referenz C
  A.refC.est.int.linreg <- as.numeric(coef(A.refC.linreg)[1])
  A.refC.ci.int.lower <- confint(A.refC.linreg, level = 0.95)[1, 1]
  A.refC.ci.int.upper <- confint(A.refC.linreg, level = 0.95)[1, 2]
  A.refC.cis.int.linreg <- A.refC.ci.int.lower <= betas[4] & A.refC.ci.int.upper >= betas[4]
  
  A.refC.est.slope.one.linreg <- as.numeric(coef(A.refC.linreg)[2])
  A.refC.ci.slope.one.lower <- confint(A.refC.linreg, level = 0.95)[2, 1]
  A.refC.ci.slope.one.upper <- confint(A.refC.linreg, level = 0.95)[2, 2]
  A.refC.cis.slope.one.linreg <- A.refC.ci.slope.one.lower <= betas[3] & A.refC.ci.slope.one.upper >= betas[3]
  
  A.refC.est.slope.two.linreg <- as.numeric(coef(A.refC.linreg)[3])
  A.refC.ci.slope.two.lower <- confint(A.refC.linreg, level = 0.95)[3, 1]
  A.refC.ci.slope.two.upper <- confint(A.refC.linreg, level = 0.95)[3, 2]
  A.refC.cis.slope.two.linreg <- A.refC.ci.slope.two.lower <= betas[3] & A.refC.ci.slope.two.upper >= betas[3]
  
  # Bootstrap-Konfidenzintervalle fuer Fall A, Referenz C
  A.refC.bootstrap <- function(data, indices) {
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + B - 1, data = d)
    y <- d$YTRUE
    A.refC.lm.boot <- lm(y ~ x, data = d)
    return(coef(A.refC.lm.boot))
  }
  
  A.refC.boot.results <- boot(data = data.dgp.case.a, statistic = A.refC.bootstrap, R = 1000)
  
  A.refC.boot.ci.int.lower <- boot.ci(A.refC.boot.results, type = "perc", index = 1)$percent[4]
  A.refC.boot.ci.int.upper <- boot.ci(A.refC.boot.results, type = "perc", index = 1)$percent[5]
  A.refC.boot.cis.int <- A.refC.boot.ci.int.lower <= betas[4] & A.refC.boot.ci.int.upper >= betas[4]
  
  A.refC.boot.ci.slope.one.lower <- boot.ci(A.refC.boot.results, type = "perc", index = 2)$percent[4]
  A.refC.boot.ci.slope.one.upper <- boot.ci(A.refC.boot.results, type = "perc", index = 2)$percent[5]
  A.refC.boot.cis.slope.one <- A.refC.boot.ci.slope.one.lower <= betas[3] & A.refC.boot.ci.slope.one.upper >= betas[3]
  
  A.refC.boot.ci.slope.two.lower <- boot.ci(A.refC.boot.results, type = "perc", index = 3)$percent[4]
  A.refC.boot.ci.slope.two.upper <- boot.ci(A.refC.boot.results, type = "perc", index = 3)$percent[5]
  A.refC.boot.cis.slope.two <- A.refC.boot.ci.slope.two.lower <= betas[3] & A.refC.boot.ci.slope.two.upper >= betas[3]
  
  A.refC.data.linreg <- data.frame(
    CASE.ID = "A",
    REF.ID = "REFC",
    TRUE.INT = betas[4],
    EST.INT = A.refC.est.int.linreg,
    CI.INT.LOWER = A.refC.ci.int.lower,
    CI.INT.UPPER = A.refC.ci.int.upper,
    BOOT.CI.INT.LOWER = A.refC.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = A.refC.boot.ci.int.upper,
    CI.INT = A.refC.cis.int.linreg,
    BOOT.CI.INT = A.refC.boot.cis.int,
    TRUE.SLOPE.ONE = betas[3],
    EST.SLOPE.ONE = A.refC.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = A.refC.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = A.refC.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = A.refC.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = A.refC.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = A.refC.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = A.refC.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[3],
    EST.SLOPE.TWO = A.refC.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = A.refC.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = A.refC.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = A.refC.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = A.refC.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = A.refC.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = A.refC.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  # FINAL DATA FRAME ----
  A.data.all <- rbind(A.refA.data.linreg, A.refB.data.linreg, A.refC.data.linreg)
  
  return(A.data.all)
}

# set.seed(123)
# data.one.sim.linreg.case.a <- A.one.sim.linreg(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.linreg.case.a, file = "2.1_RESULTS_ONESIM_LINREG_A.RData")

# =========================== CASE B - LINEARE REGRESSION ================= ####
B.one.sim.linreg <- function(n.total, betas){
  
  # DGP-Funktion fuer Fall B abrufen
  data.dgp.case.b <- dgp.case.b(n.total = n.total, betas = betas)
  
  # B.A ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.b)
  y <- data.dgp.case.b$YTRUE
  
  # Modelanpassung fuer Fall B, Referenz A 
  B.refA.linreg <- lm(y ~ x)
  
  # Geschaetzte Regresssionsparametern und KI's fuer Fall B, Referenz A
  B.refA.est.int.linreg <- as.numeric(coef(B.refA.linreg)[1])
  B.refA.ci.int.lower <- confint(B.refA.linreg, level = 0.95)[1, 1]
  B.refA.ci.int.upper <- confint(B.refA.linreg, level = 0.95)[1, 2]
  B.refA.cis.int.linreg <- B.refA.ci.int.lower <= betas[4] & B.refA.ci.int.upper >= betas[4]
  
  B.refA.est.slope.one.linreg <- as.numeric(coef(B.refA.linreg)[2])
  B.refA.ci.slope.one.lower <- confint(B.refA.linreg, level = 0.95)[2, 1]
  B.refA.ci.slope.one.upper <- confint(B.refA.linreg, level = 0.95)[2, 2]
  B.refA.cis.slope.one.linreg <- B.refA.ci.slope.one.lower <= betas[4] & B.refA.ci.slope.one.upper >= betas[4]
  
  B.refA.est.slope.two.linreg <- as.numeric(coef(B.refA.linreg)[3])
  B.refA.ci.slope.two.lower <- confint(B.refA.linreg, level = 0.95)[3, 1]
  B.refA.ci.slope.two.upper <- confint(B.refA.linreg, level = 0.95)[3, 2]
  B.refA.cis.slope.two.linreg <- B.refA.ci.slope.two.lower <= betas[2] & B.refA.ci.slope.two.upper >= betas[2]
  
  # Bootstrap-Konfidenzintervalle fuer Fall B, Referenz A
  B.refA.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ B + C - 1, data = d)
    y <- d$YTRUE
    B.refA.lm.boot <- lm(y ~ x, data = d)
    return(coef(B.refA.lm.boot))
  }
  
  B.refA.boot.results <- boot(data = data.dgp.case.b, statistic = B.refA.bootstrap, R = 1000)
  
  # Konfidenzintervalle extrahieren fuer jeden geschaetzten Regresssionsparametern
  B.refA.boot.ci.int.lower <- boot.ci(B.refA.boot.results, type = "perc", index = 1)$percent[4]
  B.refA.boot.ci.int.upper <- boot.ci(B.refA.boot.results, type = "perc", index = 1)$percent[5]
  B.refA.boot.cis.int <- B.refA.boot.ci.int.lower <= betas[4] & B.refA.boot.ci.int.upper >= betas[4]
  
  B.refA.boot.ci.slope.one.lower <- boot.ci(B.refA.boot.results, type = "perc", index = 2)$percent[4]
  B.refA.boot.ci.slope.one.upper <- boot.ci(B.refA.boot.results, type = "perc", index = 2)$percent[5]
  B.refA.boot.cis.slope.one <- B.refA.boot.ci.slope.one.lower <= betas[4] & B.refA.boot.ci.slope.one.upper >= betas[4]
  
  B.refA.boot.ci.slope.two.lower <- boot.ci(B.refA.boot.results, type = "perc", index = 3)$percent[4]
  B.refA.boot.ci.slope.two.upper <- boot.ci(B.refA.boot.results, type = "perc", index = 3)$percent[5]
  B.refA.boot.cis.slope.two <- B.refA.boot.ci.slope.two.lower <= betas[2] & B.refA.boot.ci.slope.two.upper >= betas[2]
  
  B.refA.data.linreg <- data.frame(
    CASE.ID = "B",
    REF.ID = "REFA",
    TRUE.INT = betas[4],
    EST.INT = B.refA.est.int.linreg,
    CI.INT.LOWER = B.refA.ci.int.lower,
    CI.INT.UPPER = B.refA.ci.int.upper,
    BOOT.CI.INT.LOWER = B.refA.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = B.refA.boot.ci.int.upper,
    CI.INT = B.refA.cis.int.linreg,
    BOOT.CI.INT = B.refA.boot.cis.int,
    TRUE.SLOPE.ONE = betas[4],
    EST.SLOPE.ONE = B.refA.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = B.refA.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = B.refA.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = B.refA.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = B.refA.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = B.refA.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = B.refA.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[2],
    EST.SLOPE.TWO = B.refA.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = B.refA.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = B.refA.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = B.refA.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = B.refA.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = B.refA.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = B.refA.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  # B.B ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.b)
  y <- data.dgp.case.b$YTRUE
  
  # Modelanpassung fuer Fall B, Referenz B
  B.refB.linreg <- lm(y ~ x)
  
  # Geschaetzte Regresssionsparametern und KI's fuer Fall B, Referenz B
  B.refB.est.int.linreg <- as.numeric(coef(B.refB.linreg)[1])
  B.refB.ci.int.lower <- confint(B.refB.linreg, level = 0.95)[1, 1]
  B.refB.ci.int.upper <- confint(B.refB.linreg, level = 0.95)[1, 2]
  B.refB.cis.int.linreg <- B.refB.ci.int.lower <= betas[5] & B.refB.ci.int.upper >= betas[5]
  
  B.refB.est.slope.one.linreg <- as.numeric(coef(B.refB.linreg)[2])
  B.refB.ci.slope.one.lower <- confint(B.refB.linreg, level = 0.95)[2, 1]
  B.refB.ci.slope.one.upper <- confint(B.refB.linreg, level = 0.95)[2, 2]
  B.refB.cis.slope.one.linreg <- B.refB.ci.slope.one.lower <= betas[2] & B.refB.ci.slope.one.upper >= betas[2]
  
  B.refB.est.slope.two.linreg <- as.numeric(coef(B.refB.linreg)[3])
  B.refB.ci.slope.two.lower <- confint(B.refB.linreg, level = 0.95)[3, 1]
  B.refB.ci.slope.two.upper <- confint(B.refB.linreg, level = 0.95)[3, 2]
  B.refB.cis.slope.two.linreg <- B.refB.ci.slope.two.lower <= betas[1] & B.refB.ci.slope.two.upper >= betas[1]
  
  # Bootstrap-Konfidenzintervalle fuer Fall B, Referenz B
  B.refB.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + C - 1, data = d)
    y <- d$YTRUE
    B.refB.lm.boot <- lm(y ~ x, data = d)
    return(coef(B.refB.lm.boot))
  }
  
  B.refB.boot.results <- boot(data = data.dgp.case.b, statistic = B.refB.bootstrap, R = 1000)
  
  B.refB.boot.ci.int.lower <- boot.ci(B.refB.boot.results, type = "perc", index = 1)$percent[4]
  B.refB.boot.ci.int.upper <- boot.ci(B.refB.boot.results, type = "perc", index = 1)$percent[5]
  B.refB.boot.cis.int <- B.refB.boot.ci.int.lower <= betas[5] & B.refB.boot.ci.int.upper >= betas[5]
  
  B.refB.boot.ci.slope.one.lower <- boot.ci(B.refB.boot.results, type = "perc", index = 2)$percent[4]
  B.refB.boot.ci.slope.one.upper <- boot.ci(B.refB.boot.results, type = "perc", index = 2)$percent[5]
  B.refB.boot.cis.slope.one <- B.refB.boot.ci.slope.one.lower <= betas[2] & B.refB.boot.ci.slope.one.upper >= betas[2]
  
  B.refB.boot.ci.slope.two.lower <- boot.ci(B.refB.boot.results, type = "perc", index = 3)$percent[4]
  B.refB.boot.ci.slope.two.upper <- boot.ci(B.refB.boot.results, type = "perc", index = 3)$percent[5]
  B.refB.boot.cis.slope.two <- B.refB.boot.ci.slope.two.lower <= betas[1] & B.refB.boot.ci.slope.two.upper >= betas[1]
  
  B.refB.data.linreg <- data.frame(
    CASE.ID = "B",
    REF.ID = "REFB",
    TRUE.INT = betas[5],
    EST.INT = B.refB.est.int.linreg,
    CI.INT.LOWER = B.refB.ci.int.lower,
    CI.INT.UPPER = B.refB.ci.int.upper,
    BOOT.CI.INT.LOWER = B.refB.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = B.refB.boot.ci.int.upper,
    CI.INT = B.refB.cis.int.linreg,
    BOOT.CI.INT = B.refB.boot.cis.int,
    TRUE.SLOPE.ONE = betas[2],
    EST.SLOPE.ONE = B.refB.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = B.refB.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = B.refB.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = B.refB.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = B.refB.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = B.refB.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = B.refB.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[1],
    EST.SLOPE.TWO = B.refB.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = B.refB.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = B.refB.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = B.refB.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = B.refB.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = B.refB.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = B.refB.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  # B.C ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.b)
  y <- data.dgp.case.b$YTRUE
  
  # Modelanpassung fuer Fall B, Referenz C
  B.refC.linreg <- lm(y ~ x)
  
  # Geschaetzte Regresssionsparametern und KI's fuer Fall B, Referenz C
  B.refC.est.int.linreg <- as.numeric(coef(B.refC.linreg)[1])
  B.refC.ci.int.lower <- confint(B.refC.linreg, level = 0.95)[1, 1]
  B.refC.ci.int.upper <- confint(B.refC.linreg, level = 0.95)[1, 2]
  B.refC.cis.int.linreg <- B.refC.ci.int.lower <= betas[3] & B.refC.ci.int.upper >= betas[3]
  
  B.refC.est.slope.one.linreg <- as.numeric(coef(B.refC.linreg)[2])
  B.refC.ci.slope.one.lower <- confint(B.refC.linreg, level = 0.95)[2, 1]
  B.refC.ci.slope.one.upper <- confint(B.refC.linreg, level = 0.95)[2, 2]
  B.refC.cis.slope.one.linreg <- B.refC.ci.slope.one.lower <= betas[4] & B.refC.ci.slope.one.upper >= betas[4]
  
  B.refC.est.slope.two.linreg <- as.numeric(coef(B.refC.linreg)[3])
  B.refC.ci.slope.two.lower <- confint(B.refC.linreg, level = 0.95)[3, 1]
  B.refC.ci.slope.two.upper <- confint(B.refC.linreg, level = 0.95)[3, 2]
  B.refC.cis.slope.two.linreg <- B.refC.ci.slope.two.lower <= betas[5] & B.refC.ci.slope.two.upper >= betas[5]
  
  # Bootstrap-Konfidenzintervalle fuer Fall B, Referenz C
  B.refC.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + B - 1, data = d)
    y <- d$YTRUE
    B.refC.lm.boot <- lm(y ~ x, data = d)
    return(coef(B.refC.lm.boot))
  }
  
  B.refC.boot.results <- boot(data = data.dgp.case.b, statistic = B.refC.bootstrap, R = 1000)
  
  B.refC.boot.ci.int.lower <- boot.ci(B.refC.boot.results, type = "perc", index = 1)$percent[4]
  B.refC.boot.ci.int.upper <- boot.ci(B.refC.boot.results, type = "perc", index = 1)$percent[5]
  B.refC.boot.cis.int <- B.refC.boot.ci.int.lower <= betas[3] & B.refC.boot.ci.int.upper >= betas[3]
  
  B.refC.boot.ci.slope.one.lower <- boot.ci(B.refC.boot.results, type = "perc", index = 2)$percent[4]
  B.refC.boot.ci.slope.one.upper <- boot.ci(B.refC.boot.results, type = "perc", index = 2)$percent[5]
  B.refC.boot.cis.slope.one <- B.refC.boot.ci.slope.one.lower <= betas[4] & B.refC.boot.ci.slope.one.upper >= betas[4]
  
  B.refC.boot.ci.slope.two.lower <- boot.ci(B.refC.boot.results, type = "perc", index = 3)$percent[4]
  B.refC.boot.ci.slope.two.upper <- boot.ci(B.refC.boot.results, type = "perc", index = 3)$percent[5]
  B.refC.boot.cis.slope.two <- B.refC.boot.ci.slope.two.lower <= betas[5] & B.refC.boot.ci.slope.two.upper >= betas[5]
  
  B.refC.data.linreg <- data.frame(
    CASE.ID = "B",
    REF.ID = "REFC",
    TRUE.INT = betas[3],
    EST.INT = B.refC.est.int.linreg,
    CI.INT.LOWER = B.refC.ci.int.lower,
    CI.INT.UPPER = B.refC.ci.int.upper,
    BOOT.CI.INT.LOWER = B.refC.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = B.refC.boot.ci.int.upper,
    CI.INT = B.refC.cis.int.linreg,
    BOOT.CI.INT = B.refC.boot.cis.int,
    TRUE.SLOPE.ONE = betas[4],
    EST.SLOPE.ONE = B.refC.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = B.refC.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = B.refC.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = B.refC.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = B.refC.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = B.refC.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = B.refC.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[5],
    EST.SLOPE.TWO = B.refC.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = B.refC.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = B.refC.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = B.refC.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = B.refC.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = B.refC.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = B.refC.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  # FINAL DATA FRAME ----
  B.data.all <- rbind(B.refA.data.linreg, B.refB.data.linreg, B.refC.data.linreg)
  
  return(B.data.all)
}

# set.seed(456)
# data.one.sim.linreg.case.b <- B.one.sim.linreg(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.linreg.case.b, file = "2.2_RESULTS_ONESIM_LINREG_B.RData")

# =========================== CASE C - LINEARE REGRESSION ================= ####
C.one.sim.linreg <- function(n.total, betas){
  
  # DGP-Funktion fuer Fall C abrufen
  data.dgp.case.c <- dgp.case.c(n.total = n.total, betas = betas)
  
  # Daten nach 'Case-ID' aufteilen zur weiteren Verarbeitung
  data.dgp.case.cone <- subset(data.dgp.case.c, CASE.ID == "C1")
  data.dgp.case.ctwo <- subset(data.dgp.case.c, CASE.ID == "C2")
  
  # SUBCASE C1 ----
  ## C1.A ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.cone)
  y <- data.dgp.case.cone$YTRUE
  
  # Modellanpassung
  Cone.refA.linreg <- lm(y ~ x)
  
  # Geschaetzte Regresssionsparametern und klassische Konfidenzintervalle
  Cone.refA.est.int.linreg <- as.numeric(coef(Cone.refA.linreg)[1])
  Cone.refA.ci.int.lower <- confint(Cone.refA.linreg, level = 0.95)[1, 1]
  Cone.refA.ci.int.upper <- confint(Cone.refA.linreg, level = 0.95)[1, 2]
  Cone.refA.cis.int.linreg <- Cone.refA.ci.int.lower <= betas[5] & Cone.refA.ci.int.upper >= betas[5]
  
  Cone.refA.est.slope.one.linreg <- as.numeric(coef(Cone.refA.linreg)[2])
  Cone.refA.ci.slope.one.lower <- confint(Cone.refA.linreg, level = 0.95)[2, 1]
  Cone.refA.ci.slope.one.upper <- confint(Cone.refA.linreg, level = 0.95)[2, 2]
  Cone.refA.cis.slope.one.linreg <- Cone.refA.ci.slope.one.lower <= betas[3] & Cone.refA.ci.slope.one.upper >= betas[3]
  
  Cone.refA.est.slope.two.linreg <- as.numeric(coef(Cone.refA.linreg)[3])
  Cone.refA.ci.slope.two.lower <- confint(Cone.refA.linreg, level = 0.95)[3, 1]
  Cone.refA.ci.slope.two.upper <- confint(Cone.refA.linreg, level = 0.95)[3, 2]
  Cone.refA.cis.slope.two.linreg <- Cone.refA.ci.slope.two.lower <= betas[2] & Cone.refA.ci.slope.two.upper >= betas[2]
  
  # Bootstrap-Konfidenzintervalle
  Cone.refA.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ B + C - 1, data = d)
    y <- d$YTRUE
    Cone.refA.lm.boot <- lm(y ~ x, data = d)
    return(coef(Cone.refA.lm.boot))
  }
  
  Cone.refA.boot.results <- boot(data = data.dgp.case.cone, statistic = Cone.refA.bootstrap, R = 1000)
  
  Cone.refA.boot.ci.int.lower <- boot.ci(Cone.refA.boot.results, type = "perc", index = 1)$percent[4]
  Cone.refA.boot.ci.int.upper <- boot.ci(Cone.refA.boot.results, type = "perc", index = 1)$percent[5]
  Cone.refA.boot.cis.int <- Cone.refA.boot.ci.int.lower <= betas[5] & Cone.refA.boot.ci.int.upper >= betas[5]
  
  Cone.refA.boot.ci.slope.one.lower <- boot.ci(Cone.refA.boot.results, type = "perc", index = 2)$percent[4]
  Cone.refA.boot.ci.slope.one.upper <- boot.ci(Cone.refA.boot.results, type = "perc", index = 2)$percent[5]
  Cone.refA.boot.cis.slope.one <- Cone.refA.boot.ci.slope.one.lower <= betas[3] & Cone.refA.boot.ci.slope.one.upper >= betas[3]
  
  Cone.refA.boot.ci.slope.two.lower <- boot.ci(Cone.refA.boot.results, type = "perc", index = 3)$percent[4]
  Cone.refA.boot.ci.slope.two.upper <- boot.ci(Cone.refA.boot.results, type = "perc", index = 3)$percent[5]
  Cone.refA.boot.cis.slope.two <- Cone.refA.boot.ci.slope.two.lower <= betas[2] & Cone.refA.boot.ci.slope.two.upper >= betas[2]
  
  # Ergebnisdatensatz fuer C1.A
  Cone.refA.data.linreg <- data.frame(
    CASE.ID = "C1",
    REF.ID = "REFA",
    TRUE.INT = betas[5],
    EST.INT = Cone.refA.est.int.linreg,
    CI.INT.LOWER = Cone.refA.ci.int.lower,
    CI.INT.UPPER = Cone.refA.ci.int.upper,
    BOOT.CI.INT.LOWER = Cone.refA.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = Cone.refA.boot.ci.int.upper,
    CI.INT = Cone.refA.cis.int.linreg,
    BOOT.CI.INT = Cone.refA.boot.cis.int,
    TRUE.SLOPE.ONE = betas[3],
    EST.SLOPE.ONE = Cone.refA.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = Cone.refA.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = Cone.refA.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = Cone.refA.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = Cone.refA.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = Cone.refA.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = Cone.refA.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[2],
    EST.SLOPE.TWO = Cone.refA.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = Cone.refA.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = Cone.refA.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = Cone.refA.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = Cone.refA.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = Cone.refA.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = Cone.refA.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  ## C1.B ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.cone)
  y <- data.dgp.case.cone$YTRUE
  
  Cone.refB.linreg <- lm(y ~ x)
  
  Cone.refB.est.int.linreg <- as.numeric(coef(Cone.refB.linreg)[1])
  Cone.refB.ci.int.lower <- confint(Cone.refB.linreg, level = 0.95)[1, 1]
  Cone.refB.ci.int.upper <- confint(Cone.refB.linreg, level = 0.95)[1, 2]
  Cone.refB.cis.int.linreg <- Cone.refB.ci.int.lower <= betas[5] & Cone.refB.ci.int.upper >= betas[5]
  
  Cone.refB.est.slope.one.linreg <- as.numeric(coef(Cone.refB.linreg)[2])
  Cone.refB.ci.slope.one.lower <- confint(Cone.refB.linreg, level = 0.95)[2, 1]
  Cone.refB.ci.slope.one.upper <- confint(Cone.refB.linreg, level = 0.95)[2, 2]
  Cone.refB.cis.slope.one.linreg <- Cone.refB.ci.slope.one.lower <= betas[3] & Cone.refB.ci.slope.one.upper >= betas[3]
  
  Cone.refB.est.slope.two.linreg <- as.numeric(coef(Cone.refB.linreg)[3])
  Cone.refB.ci.slope.two.lower <- confint(Cone.refB.linreg, level = 0.95)[3, 1]
  Cone.refB.ci.slope.two.upper <- confint(Cone.refB.linreg, level = 0.95)[3, 2]
  Cone.refB.cis.slope.two.linreg <- Cone.refB.ci.slope.two.lower <= betas[2] & Cone.refB.ci.slope.two.upper >= betas[2]
  
  # Bootstrap-Konfidenzintervalle
  Cone.refB.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + C - 1, data = d)
    y <- d$YTRUE
    Cone.refB.lm.boot <- lm(y ~ x, data = d)
    return(coef(Cone.refB.lm.boot))
  }
  
  Cone.refB.boot.results <- boot(data = data.dgp.case.cone, statistic = Cone.refB.bootstrap, R = 1000)
  
  Cone.refB.boot.ci.int.lower <- boot.ci(Cone.refB.boot.results, type = "perc", index = 1)$percent[4]
  Cone.refB.boot.ci.int.upper <- boot.ci(Cone.refB.boot.results, type = "perc", index = 1)$percent[5]
  Cone.refB.boot.cis.int <- Cone.refB.boot.ci.int.lower <= betas[5] & Cone.refB.boot.ci.int.upper >= betas[5]
  
  Cone.refB.boot.ci.slope.one.lower <- boot.ci(Cone.refB.boot.results, type = "perc", index = 2)$percent[4]
  Cone.refB.boot.ci.slope.one.upper <- boot.ci(Cone.refB.boot.results, type = "perc", index = 2)$percent[5]
  Cone.refB.boot.cis.slope.one <- Cone.refB.boot.ci.slope.one.lower <= betas[3] & Cone.refB.boot.ci.slope.one.upper >= betas[3]
  
  Cone.refB.boot.ci.slope.two.lower <- boot.ci(Cone.refB.boot.results, type = "perc", index = 3)$percent[4]
  Cone.refB.boot.ci.slope.two.upper <- boot.ci(Cone.refB.boot.results, type = "perc", index = 3)$percent[5]
  Cone.refB.boot.cis.slope.two <- Cone.refB.boot.ci.slope.two.lower <= betas[2] & Cone.refB.boot.ci.slope.two.upper >= betas[2]
  
  # Ergebnisdatensatz fuer C1.B
  Cone.refB.data.linreg <- data.frame(
    CASE.ID = "C1",
    REF.ID = "REFB",
    TRUE.INT = betas[5],
    EST.INT = Cone.refB.est.int.linreg,
    CI.INT.LOWER = Cone.refB.ci.int.lower,
    CI.INT.UPPER = Cone.refB.ci.int.upper,
    BOOT.CI.INT.LOWER = Cone.refB.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = Cone.refB.boot.ci.int.upper,
    CI.INT = Cone.refB.cis.int.linreg,
    BOOT.CI.INT = Cone.refB.boot.cis.int,
    TRUE.SLOPE.ONE = betas[3],
    EST.SLOPE.ONE = Cone.refB.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = Cone.refB.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = Cone.refB.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = Cone.refB.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = Cone.refB.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = Cone.refB.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = Cone.refB.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[2],
    EST.SLOPE.TWO = Cone.refB.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = Cone.refB.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = Cone.refB.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = Cone.refB.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = Cone.refB.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = Cone.refB.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = Cone.refB.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  ## C1.C ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.cone)
  y <- data.dgp.case.cone$YTRUE
  
  Cone.refC.linreg <- lm(y ~ x)
  
  Cone.refC.est.int.linreg <- as.numeric(coef(Cone.refC.linreg)[1])
  Cone.refC.ci.int.lower <- confint(Cone.refC.linreg, level = 0.95)[1, 1]
  Cone.refC.ci.int.upper <- confint(Cone.refC.linreg, level = 0.95)[1, 2]
  Cone.refC.cis.int.linreg <- Cone.refC.ci.int.lower <= betas[4] & Cone.refC.ci.int.upper >= betas[4]
  
  Cone.refC.est.slope.one.linreg <- as.numeric(coef(Cone.refC.linreg)[2])
  Cone.refC.ci.slope.one.lower <- confint(Cone.refC.linreg, level = 0.95)[2, 1]
  Cone.refC.ci.slope.one.upper <- confint(Cone.refC.linreg, level = 0.95)[2, 2]
  Cone.refC.cis.slope.one.linreg <- Cone.refC.ci.slope.one.lower <= betas[4] & Cone.refC.ci.slope.one.upper >= betas[4]
  
  Cone.refC.est.slope.two.linreg <- as.numeric(coef(Cone.refC.linreg)[3])
  Cone.refC.ci.slope.two.lower <- confint(Cone.refC.linreg, level = 0.95)[3, 1]
  Cone.refC.ci.slope.two.upper <- confint(Cone.refC.linreg, level = 0.95)[3, 2]
  Cone.refC.cis.slope.two.linreg <- Cone.refC.ci.slope.two.lower <= betas[4] & Cone.refC.ci.slope.two.upper >= betas[4]
  
  # Bootstrap-Konfidenzintervalle
  Cone.refC.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + B - 1, data = d)
    y <- d$YTRUE
    Cone.refC.lm.boot <- lm(y ~ x, data = d)
    return(coef(Cone.refC.lm.boot))
  }
  
  Cone.refC.boot.results <- boot(data = data.dgp.case.cone, statistic = Cone.refC.bootstrap, R = 1000)
  
  Cone.refC.boot.ci.int.lower <- boot.ci(Cone.refC.boot.results, type = "perc", index = 1)$percent[4]
  Cone.refC.boot.ci.int.upper <- boot.ci(Cone.refC.boot.results, type = "perc", index = 1)$percent[5]
  Cone.refC.boot.cis.int <- Cone.refC.boot.ci.int.lower <= betas[4] & Cone.refC.boot.ci.int.upper >= betas[4]
  
  Cone.refC.boot.ci.slope.one.lower <- boot.ci(Cone.refC.boot.results, type = "perc", index = 2)$percent[4]
  Cone.refC.boot.ci.slope.one.upper <- boot.ci(Cone.refC.boot.results, type = "perc", index = 2)$percent[5]
  Cone.refC.boot.cis.slope.one <- Cone.refC.boot.ci.slope.one.lower <= betas[4] & Cone.refC.boot.ci.slope.one.upper >= betas[4]
  
  Cone.refC.boot.ci.slope.two.lower <- boot.ci(Cone.refC.boot.results, type = "perc", index = 3)$percent[4]
  Cone.refC.boot.ci.slope.two.upper <- boot.ci(Cone.refC.boot.results, type = "perc", index = 3)$percent[5]
  Cone.refC.boot.cis.slope.two <- Cone.refC.boot.ci.slope.two.lower <= betas[4] & Cone.refC.boot.ci.slope.two.upper >= betas[4]
  
  # Ergebnisdatensatz fuer C1.C
  Cone.refC.data.linreg <- data.frame(
    CASE.ID = "C1",
    REF.ID = "REFC",
    TRUE.INT = betas[4],
    EST.INT = Cone.refC.est.int.linreg,
    CI.INT.LOWER = Cone.refC.ci.int.lower,
    CI.INT.UPPER = Cone.refC.ci.int.upper,
    BOOT.CI.INT.LOWER = Cone.refC.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = Cone.refC.boot.ci.int.upper,
    CI.INT = Cone.refC.cis.int.linreg,
    BOOT.CI.INT = Cone.refC.boot.cis.int,
    TRUE.SLOPE.ONE = betas[4],
    EST.SLOPE.ONE = Cone.refC.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = Cone.refC.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = Cone.refC.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = Cone.refC.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = Cone.refC.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = Cone.refC.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = Cone.refC.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[4],
    EST.SLOPE.TWO = Cone.refC.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = Cone.refC.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = Cone.refC.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = Cone.refC.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = Cone.refC.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = Cone.refC.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = Cone.refC.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  

  
  # SUBCASE C2 ----
  ## C2.A ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.ctwo)
  y <- data.dgp.case.ctwo$YTRUE
  
  Ctwo.refA.linreg <- lm(y ~ x)
  
  Ctwo.refA.est.int.linreg <- as.numeric(coef(Ctwo.refA.linreg)[1])
  Ctwo.refA.ci.int.lower <- confint(Ctwo.refA.linreg, level = 0.95)[1, 1]
  Ctwo.refA.ci.int.upper <- confint(Ctwo.refA.linreg, level = 0.95)[1, 2]
  Ctwo.refA.cis.int.linreg <- Ctwo.refA.ci.int.lower <= betas[5] & Ctwo.refA.ci.int.upper >= betas[5]
  
  Ctwo.refA.est.slope.one.linreg <- as.numeric(coef(Ctwo.refA.linreg)[2])
  Ctwo.refA.ci.slope.one.lower <- confint(Ctwo.refA.linreg, level = 0.95)[2, 1]
  Ctwo.refA.ci.slope.one.upper <- confint(Ctwo.refA.linreg, level = 0.95)[2, 2]
  Ctwo.refA.cis.slope.one.linreg <- Ctwo.refA.ci.slope.one.lower <= betas[3] & Ctwo.refA.ci.slope.one.upper >= betas[3]
  
  Ctwo.refA.est.slope.two.linreg <- as.numeric(coef(Ctwo.refA.linreg)[3])
  Ctwo.refA.ci.slope.two.lower <- confint(Ctwo.refA.linreg, level = 0.95)[3, 1]
  Ctwo.refA.ci.slope.two.upper <- confint(Ctwo.refA.linreg, level = 0.95)[3, 2]
  Ctwo.refA.cis.slope.two.linreg <- Ctwo.refA.ci.slope.two.lower <= betas[1] & Ctwo.refA.ci.slope.two.upper >= betas[1]
  
  # Bootstrap-Konfidenzintervalle
  Ctwo.refA.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ B + C - 1, data = d)
    y <- d$YTRUE
    Ctwo.refA.lm.boot <- lm(y ~ x, data = d)
    return(coef(Ctwo.refA.lm.boot))
  }
  
  Ctwo.refA.boot.results <- boot(data = data.dgp.case.ctwo, statistic = Ctwo.refA.bootstrap, R = 1000)
  
  Ctwo.refA.boot.ci.int.lower <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 1)$percent[4]
  Ctwo.refA.boot.ci.int.upper <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 1)$percent[5]
  Ctwo.refA.boot.cis.int <- Ctwo.refA.boot.ci.int.lower <= betas[5] & Ctwo.refA.boot.ci.int.upper >= betas[5]
  
  Ctwo.refA.boot.ci.slope.one.lower <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 2)$percent[4]
  Ctwo.refA.boot.ci.slope.one.upper <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 2)$percent[5]
  Ctwo.refA.boot.cis.slope.one <- Ctwo.refA.boot.ci.slope.one.lower <= betas[3] & Ctwo.refA.boot.ci.slope.one.upper >= betas[3]
  
  Ctwo.refA.boot.ci.slope.two.lower <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 3)$percent[4]
  Ctwo.refA.boot.ci.slope.two.upper <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 3)$percent[5]
  Ctwo.refA.boot.cis.slope.two <- Ctwo.refA.boot.ci.slope.two.lower <= betas[1] & Ctwo.refA.boot.ci.slope.two.upper >= betas[1]
  
  # Ergebnisdatensatz fuer C2.A
  Ctwo.refA.data.linreg <- data.frame(
    CASE.ID = "C2",
    REF.ID = "REFA",
    TRUE.INT = betas[5],
    EST.INT = Ctwo.refA.est.int.linreg,
    CI.INT.LOWER = Ctwo.refA.ci.int.lower,
    CI.INT.UPPER = Ctwo.refA.ci.int.upper,
    BOOT.CI.INT.LOWER = Ctwo.refA.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = Ctwo.refA.boot.ci.int.upper,
    CI.INT = Ctwo.refA.cis.int.linreg,
    BOOT.CI.INT = Ctwo.refA.boot.cis.int,
    TRUE.SLOPE.ONE = betas[3],
    EST.SLOPE.ONE = Ctwo.refA.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = Ctwo.refA.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = Ctwo.refA.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = Ctwo.refA.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = Ctwo.refA.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = Ctwo.refA.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = Ctwo.refA.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[1],
    EST.SLOPE.TWO = Ctwo.refA.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = Ctwo.refA.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = Ctwo.refA.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = Ctwo.refA.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = Ctwo.refA.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = Ctwo.refA.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = Ctwo.refA.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  ## C2.B ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.ctwo)
  y <- data.dgp.case.ctwo$YTRUE
  
  Ctwo.refB.linreg <- lm(y ~ x)
  
  Ctwo.refB.est.int.linreg <- as.numeric(coef(Ctwo.refB.linreg)[1])
  Ctwo.refB.ci.int.lower <- confint(Ctwo.refB.linreg, level = 0.95)[1, 1]
  Ctwo.refB.ci.int.upper <- confint(Ctwo.refB.linreg, level = 0.95)[1, 2]
  Ctwo.refB.cis.int.linreg <- Ctwo.refB.ci.int.lower <= betas[5] & Ctwo.refB.ci.int.upper >= betas[5]
  
  Ctwo.refB.est.slope.one.linreg <- as.numeric(coef(Ctwo.refB.linreg)[2])
  Ctwo.refB.ci.slope.one.lower <- confint(Ctwo.refB.linreg, level = 0.95)[2, 1]
  Ctwo.refB.ci.slope.one.upper <- confint(Ctwo.refB.linreg, level = 0.95)[2, 2]
  Ctwo.refB.cis.slope.one.linreg <- Ctwo.refB.ci.slope.one.lower <= betas[3] & Ctwo.refB.ci.slope.one.upper >= betas[3]
  
  Ctwo.refB.est.slope.two.linreg <- as.numeric(coef(Ctwo.refB.linreg)[3])
  Ctwo.refB.ci.slope.two.lower <- confint(Ctwo.refB.linreg, level = 0.95)[3, 1]
  Ctwo.refB.ci.slope.two.upper <- confint(Ctwo.refB.linreg, level = 0.95)[3, 2]
  Ctwo.refB.cis.slope.two.linreg <- Ctwo.refB.ci.slope.two.lower <= betas[1] & Ctwo.refB.ci.slope.two.upper >= betas[1]
  
  # Bootstrap-Konfidenzintervalle
  Ctwo.refB.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + C - 1, data = d)
    y <- d$YTRUE
    Ctwo.refB.lm.boot <- lm(y ~ x, data = d)
    return(coef(Ctwo.refB.lm.boot))
  }
  
  Ctwo.refB.boot.results <- boot(data = data.dgp.case.ctwo, statistic = Ctwo.refB.bootstrap, R = 1000)
  
  Ctwo.refB.boot.ci.int.lower <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 1)$percent[4]
  Ctwo.refB.boot.ci.int.upper <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 1)$percent[5]
  Ctwo.refB.boot.cis.int <- Ctwo.refB.boot.ci.int.lower <= betas[5] & Ctwo.refB.boot.ci.int.upper >= betas[5]
  
  Ctwo.refB.boot.ci.slope.one.lower <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 2)$percent[4]
  Ctwo.refB.boot.ci.slope.one.upper <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 2)$percent[5]
  Ctwo.refB.boot.cis.slope.one <- Ctwo.refB.boot.ci.slope.one.lower <= betas[3] & Ctwo.refB.boot.ci.slope.one.upper >= betas[3]
  
  Ctwo.refB.boot.ci.slope.two.lower <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 3)$percent[4]
  Ctwo.refB.boot.ci.slope.two.upper <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 3)$percent[5]
  Ctwo.refB.boot.cis.slope.two <- Ctwo.refB.boot.ci.slope.two.lower <= betas[1] & Ctwo.refB.boot.ci.slope.two.upper >= betas[1]
  
  # Ergebnisdatensatz fuer C2.B
  Ctwo.refB.data.linreg <- data.frame(
    CASE.ID = "C2",
    REF.ID = "REFB",
    TRUE.INT = betas[5],
    EST.INT = Ctwo.refB.est.int.linreg,
    CI.INT.LOWER = Ctwo.refB.ci.int.lower,
    CI.INT.UPPER = Ctwo.refB.ci.int.upper,
    BOOT.CI.INT.LOWER = Ctwo.refB.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = Ctwo.refB.boot.ci.int.upper,
    CI.INT = Ctwo.refB.cis.int.linreg,
    BOOT.CI.INT = Ctwo.refB.boot.cis.int,
    TRUE.SLOPE.ONE = betas[3],
    EST.SLOPE.ONE = Ctwo.refB.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = Ctwo.refB.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = Ctwo.refB.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = Ctwo.refB.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = Ctwo.refB.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = Ctwo.refB.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = Ctwo.refB.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[1],
    EST.SLOPE.TWO = Ctwo.refB.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = Ctwo.refB.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = Ctwo.refB.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = Ctwo.refB.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = Ctwo.refB.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = Ctwo.refB.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = Ctwo.refB.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  ## C2.C ----
  # Vorbereitung der Daten fuer Lineare Regression
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.ctwo)
  y <- data.dgp.case.ctwo$YTRUE
  
  Ctwo.refC.linreg <- lm(y ~ x)
  
  Ctwo.refC.est.int.linreg <- as.numeric(coef(Ctwo.refC.linreg)[1])
  Ctwo.refC.ci.int.lower <- confint(Ctwo.refC.linreg, level = 0.95)[1, 1]
  Ctwo.refC.ci.int.upper <- confint(Ctwo.refC.linreg, level = 0.95)[1, 2]
  Ctwo.refC.cis.int.linreg <- Ctwo.refC.ci.int.lower <= betas[3] & Ctwo.refC.ci.int.upper >= betas[3]
  
  Ctwo.refC.est.slope.one.linreg <- as.numeric(coef(Ctwo.refC.linreg)[2])
  Ctwo.refC.ci.slope.one.lower <- confint(Ctwo.refC.linreg, level = 0.95)[2, 1]
  Ctwo.refC.ci.slope.one.upper <- confint(Ctwo.refC.linreg, level = 0.95)[2, 2]
  Ctwo.refC.cis.slope.one.linreg <- Ctwo.refC.ci.slope.one.lower <= betas[5] & Ctwo.refC.ci.slope.one.upper >= betas[5]
  
  Ctwo.refC.est.slope.two.linreg <- as.numeric(coef(Ctwo.refC.linreg)[3])
  Ctwo.refC.ci.slope.two.lower <- confint(Ctwo.refC.linreg, level = 0.95)[3, 1]
  Ctwo.refC.ci.slope.two.upper <- confint(Ctwo.refC.linreg, level = 0.95)[3, 2]
  Ctwo.refC.cis.slope.two.linreg <- Ctwo.refC.ci.slope.two.lower <= betas[5] & Ctwo.refC.ci.slope.two.upper >= betas[5]
  
  # Bootstrap-Konfidenzintervalle
  Ctwo.refC.bootstrap <- function(data, indices){
    d <- data[indices, ]
    x <- model.matrix(YTRUE ~ A + B - 1, data = d)
    y <- d$YTRUE
    Ctwo.refC.lm.boot <- lm(y ~ x, data = d)
    return(coef(Ctwo.refC.lm.boot))
  }
  
  Ctwo.refC.boot.results <- boot(data = data.dgp.case.ctwo, statistic = Ctwo.refC.bootstrap, R = 1000)
  
  Ctwo.refC.boot.ci.int.lower <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 1)$percent[4]
  Ctwo.refC.boot.ci.int.upper <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 1)$percent[5]
  Ctwo.refC.boot.cis.int <- Ctwo.refC.boot.ci.int.lower <= betas[3] & Ctwo.refC.boot.ci.int.upper >= betas[3]
  
  Ctwo.refC.boot.ci.slope.one.lower <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 2)$percent[4]
  Ctwo.refC.boot.ci.slope.one.upper <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 2)$percent[5]
  Ctwo.refC.boot.cis.slope.one <- Ctwo.refC.boot.ci.slope.one.lower <= betas[5] & Ctwo.refC.boot.ci.slope.one.upper >= betas[5]
  
  Ctwo.refC.boot.ci.slope.two.lower <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 3)$percent[4]
  Ctwo.refC.boot.ci.slope.two.upper <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 3)$percent[5]
  Ctwo.refC.boot.cis.slope.two <- Ctwo.refC.boot.ci.slope.two.lower <= betas[5] & Ctwo.refC.boot.ci.slope.two.upper >= betas[5]
  
  # Ergebnisdatensatz fuer C2.C
  Ctwo.refC.data.linreg <- data.frame(
    CASE.ID = "C2",
    REF.ID = "REFC",
    TRUE.INT = betas[3],
    EST.INT = Ctwo.refC.est.int.linreg,
    CI.INT.LOWER = Ctwo.refC.ci.int.lower,
    CI.INT.UPPER = Ctwo.refC.ci.int.upper,
    BOOT.CI.INT.LOWER = Ctwo.refC.boot.ci.int.lower,
    BOOT.CI.INT.UPPER = Ctwo.refC.boot.ci.int.upper,
    CI.INT = Ctwo.refC.cis.int.linreg,
    BOOT.CI.INT = Ctwo.refC.boot.cis.int,
    TRUE.SLOPE.ONE = betas[5],
    EST.SLOPE.ONE = Ctwo.refC.est.slope.one.linreg,
    CI.SLOPE.ONE.LOWER = Ctwo.refC.ci.slope.one.lower,
    CI.SLOPE.ONE.UPPER = Ctwo.refC.ci.slope.one.upper,
    BOOT.CI.SLOPE.ONE.LOWER = Ctwo.refC.boot.ci.slope.one.lower,
    BOOT.CI.SLOPE.ONE.UPPER = Ctwo.refC.boot.ci.slope.one.upper,
    CI.SLOPE.ONE = Ctwo.refC.cis.slope.one.linreg,
    BOOT.CI.SLOPE.ONE = Ctwo.refC.boot.cis.slope.one,
    TRUE.SLOPE.TWO = betas[5],
    EST.SLOPE.TWO = Ctwo.refC.est.slope.two.linreg,
    CI.SLOPE.TWO.LOWER = Ctwo.refC.ci.slope.two.lower,
    CI.SLOPE.TWO.UPPER = Ctwo.refC.ci.slope.two.upper,
    BOOT.CI.SLOPE.TWO.LOWER = Ctwo.refC.boot.ci.slope.two.lower,
    BOOT.CI.SLOPE.TWO.UPPER = Ctwo.refC.boot.ci.slope.two.upper,
    CI.SLOPE.TWO = Ctwo.refC.cis.slope.two.linreg,
    BOOT.CI.SLOPE.TWO = Ctwo.refC.boot.cis.slope.two,
    METHOD = "LinReg"
  )
  
  # FINAL DATA FRAME ----
  
  C.data.all <- rbind(Cone.refA.data.linreg, Cone.refB.data.linreg,Cone.refC.data.linreg,
                      Ctwo.refA.data.linreg, Ctwo.refB.data.linreg,Ctwo.refC.data.linreg)
  
  return(C.data.all)
}

# set.seed(789)
# data.one.sim.linreg.case.c <- C.one.sim.linreg(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.linreg.case.c, file = "2.3_RESULTS_ONESIM_LINREG_C.RData")

# ========================= ALL CASES - LINEARE REGRESSION ================ ####
one.sim.linreg.all.cases <- function(n.total, betas){
  
  A.one.sim.data <- A.one.sim.linreg(n.total = n.total, betas = betas)
  B.one.sim.data <- B.one.sim.linreg(n.total = n.total, betas = betas)
  C.one.sim.data <- C.one.sim.linreg(n.total = n.total, betas = betas)
  
  all.cases.one.sim.data <- rbind(A.one.sim.data,
                                  B.one.sim.data,
                                  C.one.sim.data)
  
  return(all.cases.one.sim.data)
}

# set.seed(123)
# data.one.sim.linreg.all.cases <- one.sim.linreg.all.cases(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.linreg.all.cases, file = "2.4_RESULTS_ONESIM_LINREG_ALL_CASES.RData")

# ============================================================================ #
# ENDE DES SKRIPTS -----------------------------------------------------------
# ============================================================================ #