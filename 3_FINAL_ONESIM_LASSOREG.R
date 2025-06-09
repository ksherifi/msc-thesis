# ============================================================================ #
# ALLGEMEINE INFORMATIONEN ZUM R-SKRIPT
# ============================================================================ #
# Dieses Skript fuehrt Lasso Regressionen fuer mehrere Szenarien (Fall A, B und C) durch.
# Es beinhaltet die Schaetzung von Regressionsparametern, die Verwendung von Kreuzvalidierung 
# zur Auswahl des optimalen Lambda-Werts und die Berechnung von Bootstrap-
# Konfidenzintervallen, um die Robustheit der Schaetzungen zu pruefen.
# 
# Uebersicht:
# 1. Datenvorbereitung fuer die Regression.
# 2. Durchfuehrung der Lasso Regression fuer Fall A, B und C.
# 3. Schaetzung der Regressionsmodelle und Berechnung der Bootstrap-Konfidenzintervalle.
# 4. Ueberpruefung ob wahrer Regressionsparameter in den Konfidenzintervallen liegt
# 
# Verwendete Pakete:
# - 'glmnet' fuer Lasso Regression.
# - 'boot' fuer Bootstrapping.
# 
# ---------------------------------------------------------------------------- #
# ============================================================================ #

# ============================================================================ #
# R-SKRIPT MIT DGP-FUNKTIONEN LADEN ------------------------------------------
# ============================================================================ #

source("1_FINAL_DGP.R")  # Laedt die DGP-Funktionen

# ============================================================================ #
# BENOETIGTE PAKETE INSTALLIEREN UND LADEN -----------------------------------
# ============================================================================ #

# Installieren der benoetigten Pakete (falls nicht bereits installiert)
# install.packages("glmnet")
# install.packages("boot")

# Laden der Pakete
library(glmnet)  # Fuer Lasso Regression
library(boot)  # Fuer Bootstrapping

# =========================== CASE A - LASSO REGRESSION ================== ####
A.one.sim.lasso <- function(n.total, betas){
  
  # DGP-Funktion fuer Fall A abrufen
  data.dgp.case.a <- dgp.case.a(n.total = n.total, betas = betas)
  
  # A.A ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.a)
  y <- data.dgp.case.a$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz A (mit Achsenabschnitt)
  A.refA.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  A.refA.best.lambda <- A.refA.cv.fit$lambda.min
  A.refA.lasso.fit <- glmnet(x, y, alpha = 1, lambda = A.refA.best.lambda,
                             standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  A.refA.coef <- as.numeric(coef(A.refA.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall A, Referenz A
  A.refA.est.int.lasso <- A.refA.coef[1]  # Achsenabschnitt
  A.refA.est.slope.one.lasso <- A.refA.coef[2]  # Parameter fuer B
  A.refA.est.slope.two.lasso <- A.refA.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  A.refA.bootstrap <- function(data, indices){
    d <- data[indices,]
    A.refA.xboot <- model.matrix(YTRUE ~ B + C - 1, data = d)
    A.refA.yboot <- d$YTRUE
    A.refA.lasso.mod <- glmnet(A.refA.xboot, A.refA.yboot, alpha = 1, 
                               lambda = A.refA.best.lambda, standardize = FALSE)
    return(as.numeric(coef(A.refA.lasso.mod)))
  }
  
  # Bootstrapping durchfuehren
  A.refA.boot.results <- boot(data = data.dgp.case.a, statistic = A.refA.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  A.refA.ci.int <- boot.ci(A.refA.boot.results, type = "perc", index = 1)$percent[4:5]
  A.refA.ci.slope.one <- boot.ci(A.refA.boot.results, type = "perc", index = 2)$percent[4:5]
  A.refA.ci.slope.two <- boot.ci(A.refA.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefung, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  A.refA.cis.int.lasso <- A.refA.ci.int[1] <= betas[4] & A.refA.ci.int[2] >= betas[4]
  A.refA.cis.slope.one.lasso <- A.refA.ci.slope.one[1] <= betas[3] & A.refA.ci.slope.one[2] >= betas[3]
  A.refA.cis.slope.two.lasso <- A.refA.ci.slope.two[1] <= betas[3] & A.refA.ci.slope.two[2] >= betas[3]
  
  # Konfidenzintervalle numerisch abspeichern
  A.refA.ci.int.lower <- A.refA.ci.int[1]
  A.refA.ci.int.upper <- A.refA.ci.int[2]
  A.refA.ci.slope.one.lower <- A.refA.ci.slope.one[1]
  A.refA.ci.slope.one.upper <- A.refA.ci.slope.one[2]
  A.refA.ci.slope.two.lower <- A.refA.ci.slope.two[1]
  A.refA.ci.slope.two.upper <- A.refA.ci.slope.two[2]
  
  # Ergebnisdatensatz erstellen
  A.refA.data.lasso <- data.frame(CASE.ID = "A",
                                  REF.ID = "REFA",
                                  TRUE.INT = betas[4],
                                  EST.INT = A.refA.est.int.lasso,
                                  BOOT.CI.INT = A.refA.cis.int.lasso,
                                  BOOT.CI.INT.LOWER = A.refA.ci.int.lower,
                                  BOOT.CI.INT.UPPER = A.refA.ci.int.upper,
                                  TRUE.SLOPE.ONE = betas[3],
                                  EST.SLOPE.ONE = A.refA.est.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE = A.refA.cis.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE.LOWER = A.refA.ci.slope.one.lower,
                                  BOOT.CI.SLOPE.ONE.UPPER = A.refA.ci.slope.one.upper,
                                  TRUE.SLOPE.TWO = betas[3],
                                  EST.SLOPE.TWO = A.refA.est.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO = A.refA.cis.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO.LOWER = A.refA.ci.slope.two.lower,
                                  BOOT.CI.SLOPE.TWO.UPPER = A.refA.ci.slope.two.upper,
                                  BEST.LAMBDA = A.refA.best.lambda,
                                  METHOD = "Lasso")
  
  # A.B ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.a)
  y <- data.dgp.case.a$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz B (mit Achsenabschnitt)
  A.refB.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  A.refB.best.lambda <- A.refB.cv.fit$lambda.min
  A.refB.lasso.fit <- glmnet(x, y, alpha = 1, lambda = A.refB.best.lambda, 
                             standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  A.refB.coef <- as.numeric(coef(A.refB.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall A, Referenz B
  A.refB.est.int.lasso <- A.refB.coef[1]  # Achsenabschnitt
  A.refB.est.slope.one.lasso <- A.refB.coef[2]  # Parameter fuer A
  A.refB.est.slope.two.lasso <- A.refB.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  A.refB.bootstrap <- function(data, indices){
    d <- data[indices,]
    A.refB.xboot <- model.matrix(YTRUE ~ A + C - 1, data = d)
    A.refB.yboot <- d$YTRUE
    A.refB.lasso.mod <- glmnet(A.refB.xboot, A.refB.yboot, alpha = 1, 
                               lambda = A.refB.best.lambda, standardize = FALSE)
    return(as.numeric(coef(A.refB.lasso.mod)))
  }
  
  A.refB.boot.results <- boot(data = data.dgp.case.a, statistic = A.refB.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  A.refB.ci.int <- boot.ci(A.refB.boot.results, type = "perc", index = 1)$percent[4:5]
  A.refB.ci.slope.one <- boot.ci(A.refB.boot.results, type = "perc", index = 2)$percent[4:5]
  A.refB.ci.slope.two <- boot.ci(A.refB.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  A.refB.cis.int.lasso <- A.refB.ci.int[1] <= betas[4] & A.refB.ci.int[2] >= betas[4]
  A.refB.cis.slope.one.lasso <- A.refB.ci.slope.one[1] <= betas[3] & A.refB.ci.slope.one[2] >= betas[3]
  A.refB.cis.slope.two.lasso <- A.refB.ci.slope.two[1] <= betas[3] & A.refB.ci.slope.two[2] >= betas[3]
  
  # Konfidenzintervalle numerisch abspeichern
  A.refB.ci.int.lower <- A.refB.ci.int[1]
  A.refB.ci.int.upper <- A.refB.ci.int[2]
  A.refB.ci.slope.one.lower <- A.refB.ci.slope.one[1]
  A.refB.ci.slope.one.upper <- A.refB.ci.slope.one[2]
  A.refB.ci.slope.two.lower <- A.refB.ci.slope.two[1]
  A.refB.ci.slope.two.upper <- A.refB.ci.slope.two[2]
  
  # Ergebnisdatensatz mit den numerischen Konfidenzintervallen
  A.refB.data.lasso <- data.frame(CASE.ID = "A",
                                  REF.ID = "REFB",
                                  TRUE.INT = betas[4],
                                  EST.INT = A.refB.est.int.lasso,
                                  BOOT.CI.INT = A.refB.cis.int.lasso,
                                  BOOT.CI.INT.LOWER = A.refB.ci.int.lower,
                                  BOOT.CI.INT.UPPER = A.refB.ci.int.upper,
                                  TRUE.SLOPE.ONE = betas[3],
                                  EST.SLOPE.ONE = A.refB.est.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE = A.refB.cis.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE.LOWER = A.refB.ci.slope.one.lower,
                                  BOOT.CI.SLOPE.ONE.UPPER = A.refB.ci.slope.one.upper,
                                  TRUE.SLOPE.TWO = betas[3],
                                  EST.SLOPE.TWO = A.refB.est.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO = A.refB.cis.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO.LOWER = A.refB.ci.slope.two.lower,
                                  BOOT.CI.SLOPE.TWO.UPPER = A.refB.ci.slope.two.upper,
                                  BEST.LAMBDA = A.refB.best.lambda,
                                  METHOD = "Lasso")
  
  # A.C ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.a)
  y <- data.dgp.case.a$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz B (mit Achsenabschnitt)
  A.refC.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  A.refC.best.lambda <- A.refC.cv.fit$lambda.min
  A.refC.lasso.fit <- glmnet(x, y, alpha = 1, lambda = A.refC.best.lambda, 
                             standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  A.refC.coef <- as.numeric(coef(A.refC.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall A, Referenz B
  A.refC.est.int.lasso <- A.refC.coef[1]  # Achsenabschnitt
  A.refC.est.slope.one.lasso <- A.refC.coef[2]  # Parameter fuer A
  A.refC.est.slope.two.lasso <- A.refC.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  A.refC.bootstrap <- function(data, indices){
    d <- data[indices,]
    A.refC.xboot <- model.matrix(YTRUE ~ A + B - 1, data = d)
    A.refC.yboot <- d$YTRUE
    A.refC.lasso.mod <- glmnet(A.refC.xboot, A.refC.yboot, alpha = 1, 
                               lambda = A.refC.best.lambda, standardize = FALSE)
    return(as.numeric(coef(A.refC.lasso.mod)))
  }
  
  A.refC.boot.results <- boot(data = data.dgp.case.a, statistic = A.refC.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  A.refC.ci.int <- boot.ci(A.refC.boot.results, type = "perc", index = 1)$percent[4:5]
  A.refC.ci.slope.one <- boot.ci(A.refC.boot.results, type = "perc", index = 2)$percent[4:5]
  A.refC.ci.slope.two <- boot.ci(A.refC.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  A.refC.cis.int.lasso <- A.refC.ci.int[1] <= betas[4] & A.refC.ci.int[2] >= betas[4]
  A.refC.cis.slope.one.lasso <- A.refC.ci.slope.one[1] <= betas[3] & A.refC.ci.slope.one[2] >= betas[3]
  A.refC.cis.slope.two.lasso <- A.refC.ci.slope.two[1] <= betas[3] & A.refC.ci.slope.two[2] >= betas[3]
  
  # Konfidenzintervalle numerisch abspeichern
  A.refC.ci.int.lower <- A.refC.ci.int[1]
  A.refC.ci.int.upper <- A.refC.ci.int[2]
  A.refC.ci.slope.one.lower <- A.refC.ci.slope.one[1]
  A.refC.ci.slope.one.upper <- A.refC.ci.slope.one[2]
  A.refC.ci.slope.two.lower <- A.refC.ci.slope.two[1]
  A.refC.ci.slope.two.upper <- A.refC.ci.slope.two[2]
  
  # Ergebnisdatensatz mit den numerischen Konfidenzintervallen
  A.refC.data.lasso <- data.frame(CASE.ID = "A",
                                  REF.ID = "REFC",
                                  TRUE.INT = betas[4],
                                  EST.INT = A.refC.est.int.lasso,
                                  BOOT.CI.INT = A.refC.cis.int.lasso,
                                  BOOT.CI.INT.LOWER = A.refC.ci.int.lower,
                                  BOOT.CI.INT.UPPER = A.refC.ci.int.upper,
                                  TRUE.SLOPE.ONE = betas[3],
                                  EST.SLOPE.ONE = A.refC.est.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE = A.refC.cis.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE.LOWER = A.refC.ci.slope.one.lower,
                                  BOOT.CI.SLOPE.ONE.UPPER = A.refC.ci.slope.one.upper,
                                  TRUE.SLOPE.TWO = betas[3],
                                  EST.SLOPE.TWO = A.refC.est.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO = A.refC.cis.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO.LOWER = A.refC.ci.slope.two.lower,
                                  BOOT.CI.SLOPE.TWO.UPPER = A.refC.ci.slope.two.upper,
                                  BEST.LAMBDA = A.refC.best.lambda,
                                  METHOD = "Lasso")
  
  # FINAL DATA FRAME ----
  A.data.all <- rbind(A.refA.data.lasso, A.refB.data.lasso, A.refC.data.lasso)
  
  return(A.data.all)
}

# set.seed(123)
# data.one.sim.lassoreg.case.a <- A.one.sim.lasso(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.lassoreg.case.a, file = "3.1_RESULTS_ONESIM_LASSO_A.RData")

# =========================== CASE B - LASSO REGRESSION ================== ####
B.one.sim.lasso <- function(n.total, betas){
  
  # DGP-Funktion fuer Fall B abrufen
  data.dgp.case.b <- dgp.case.b(n.total = n.total, betas = betas)
  
  # B.A. ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.b)
  y <- data.dgp.case.b$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz A (mit Achsenabschnitt)
  B.refA.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  B.refA.best.lambda <- B.refA.cv.fit$lambda.min
  B.refA.lasso.fit <- glmnet(x, y, alpha = 1, lambda = B.refA.best.lambda, 
                             standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  B.refA.coef <- as.numeric(coef(B.refA.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall B, Referenz A
  B.refA.est.int.lasso <- B.refA.coef[1]  # Achsenabschnitt
  B.refA.est.slope.one.lasso <- B.refA.coef[2]  # Parameter fuer B
  B.refA.est.slope.two.lasso <- B.refA.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  B.refA.bootstrap <- function(data, indices){
    d <- data[indices,]
    B.refA.xboot <- model.matrix(YTRUE ~ B + C - 1, data = d)
    B.refA.yboot <- d$YTRUE
    B.refA.lasso.mod <- glmnet(B.refA.xboot, B.refA.yboot, alpha = 1, 
                               lambda = B.refA.best.lambda, standardize = FALSE)
    return(as.numeric(coef(B.refA.lasso.mod)))
  }
  
  B.refA.boot.results <- boot(data = data.dgp.case.b, statistic = B.refA.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  B.refA.ci.int <- boot.ci(B.refA.boot.results, type = "perc", index = 1)$percent[4:5]
  B.refA.ci.slope.one <- boot.ci(B.refA.boot.results, type = "perc", index = 2)$percent[4:5]
  B.refA.ci.slope.two <- boot.ci(B.refA.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  B.refA.cis.int.lasso <- B.refA.ci.int[1] <= betas[4] & B.refA.ci.int[2] >= betas[4]
  B.refA.cis.slope.one.lasso <- B.refA.ci.slope.one[1] <= betas[4] & B.refA.ci.slope.one[2] >= betas[4]
  B.refA.cis.slope.two.lasso <- B.refA.ci.slope.two[1] <= betas[2] & B.refA.ci.slope.two[2] >= betas[2]
  
  # Konfidenzintervalle numerisch abspeichern
  B.refA.ci.int.lower <- B.refA.ci.int[1]
  B.refA.ci.int.upper <- B.refA.ci.int[2]
  B.refA.ci.slope.one.lower <- B.refA.ci.slope.one[1]
  B.refA.ci.slope.one.upper <- B.refA.ci.slope.one[2]
  B.refA.ci.slope.two.lower <- B.refA.ci.slope.two[1]
  B.refA.ci.slope.two.upper <- B.refA.ci.slope.two[2]
  
  # Ergebnisdatensatz mit den numerischen Konfidenzintervallen
  B.refA.data.lasso <- data.frame(CASE.ID = "B",
                                  REF.ID = "REFA",
                                  TRUE.INT = betas[4],
                                  EST.INT = B.refA.est.int.lasso,
                                  BOOT.CI.INT = B.refA.cis.int.lasso,
                                  BOOT.CI.INT.LOWER = B.refA.ci.int.lower,
                                  BOOT.CI.INT.UPPER = B.refA.ci.int.upper,
                                  TRUE.SLOPE.ONE = betas[4],
                                  EST.SLOPE.ONE = B.refA.est.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE = B.refA.cis.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE.LOWER = B.refA.ci.slope.one.lower,
                                  BOOT.CI.SLOPE.ONE.UPPER = B.refA.ci.slope.one.upper,
                                  TRUE.SLOPE.TWO = betas[2],
                                  EST.SLOPE.TWO = B.refA.est.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO = B.refA.cis.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO.LOWER = B.refA.ci.slope.two.lower,
                                  BOOT.CI.SLOPE.TWO.UPPER = B.refA.ci.slope.two.upper,
                                  BEST.LAMBDA = B.refA.best.lambda,
                                  METHOD = "Lasso")
  
  # B.B ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.b)
  y <- data.dgp.case.b$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz B (mit Achsenabschnitt)
  B.refB.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  B.refB.best.lambda <- B.refB.cv.fit$lambda.min
  B.refB.lasso.fit <- glmnet(x, y, alpha = 1, lambda = B.refB.best.lambda,
                             standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  B.refB.coef <- as.numeric(coef(B.refB.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall B, Referenz B
  B.refB.est.int.lasso <- B.refB.coef[1]  # Achsenabschnitt
  B.refB.est.slope.one.lasso <- B.refB.coef[2]  # Parameter fuer A
  B.refB.est.slope.two.lasso <- B.refB.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  B.refB.bootstrap <- function(data, indices){
    d <- data[indices,]
    B.refB.xboot <- model.matrix(YTRUE ~ A + C - 1, data = d)
    B.refB.yboot <- d$YTRUE
    B.refB.lasso.mod <- glmnet(B.refB.xboot, B.refB.yboot, alpha = 1, 
                               lambda = B.refB.best.lambda, standardize = FALSE)
    return(as.numeric(coef(B.refB.lasso.mod)))
  }
  
  B.refB.boot.results <- boot(data = data.dgp.case.b, statistic = B.refB.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  B.refB.ci.int <- boot.ci(B.refB.boot.results, type = "perc", index = 1)$percent[4:5]
  B.refB.ci.slope.one <- boot.ci(B.refB.boot.results, type = "perc", index = 2)$percent[4:5]
  B.refB.ci.slope.two <- boot.ci(B.refB.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  B.refB.cis.int.lasso <- B.refB.ci.int[1] <= betas[5] & B.refB.ci.int[2] >= betas[5]
  B.refB.cis.slope.one.lasso <- B.refB.ci.slope.one[1] <= betas[2] & B.refB.ci.slope.one[2] >= betas[2]
  B.refB.cis.slope.two.lasso <- B.refB.ci.slope.two[1] <= betas[1] & B.refB.ci.slope.two[2] >= betas[1]
  
  # Konfidenzintervalle numerisch abspeichern
  B.refB.ci.int.lower <- B.refB.ci.int[1]
  B.refB.ci.int.upper <- B.refB.ci.int[2]
  B.refB.ci.slope.one.lower <- B.refB.ci.slope.one[1]
  B.refB.ci.slope.one.upper <- B.refB.ci.slope.one[2]
  B.refB.ci.slope.two.lower <- B.refB.ci.slope.two[1]
  B.refB.ci.slope.two.upper <- B.refB.ci.slope.two[2]
  
  # Ergebnisdatensatz mit den numerischen Konfidenzintervallen
  B.refB.data.lasso <- data.frame(CASE.ID = "B",
                                  REF.ID = "REFB",
                                  TRUE.INT = betas[5],
                                  EST.INT = B.refB.est.int.lasso,
                                  BOOT.CI.INT = B.refB.cis.int.lasso,
                                  BOOT.CI.INT.LOWER = B.refB.ci.int.lower,
                                  BOOT.CI.INT.UPPER = B.refB.ci.int.upper,
                                  TRUE.SLOPE.ONE = betas[2],
                                  EST.SLOPE.ONE = B.refB.est.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE = B.refB.cis.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE.LOWER = B.refB.ci.slope.one.lower,
                                  BOOT.CI.SLOPE.ONE.UPPER = B.refB.ci.slope.one.upper,
                                  TRUE.SLOPE.TWO = betas[1],
                                  EST.SLOPE.TWO = B.refB.est.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO = B.refB.cis.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO.LOWER = B.refB.ci.slope.two.lower,
                                  BOOT.CI.SLOPE.TWO.UPPER = B.refB.ci.slope.two.upper,
                                  BEST.LAMBDA = B.refB.best.lambda,
                                  METHOD = "Lasso")
  
  # B.C ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.b)
  y <- data.dgp.case.b$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz C (mit Achsenabschnitt)
  B.refC.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  B.refC.best.lambda <- B.refC.cv.fit$lambda.min
  B.refC.lasso.fit <- glmnet(x, y, alpha = 1, lambda = B.refC.best.lambda,
                             standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  B.refC.coef <- as.numeric(coef(B.refC.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall B, Referenz C
  B.refC.est.int.lasso <- B.refC.coef[1]  # Achsenabschnitt
  B.refC.est.slope.one.lasso <- B.refC.coef[2]  # Parameter fuer A
  B.refC.est.slope.two.lasso <- B.refC.coef[3]  # Parameter fuer B
  
  # Bootstrap fuer Konfidenzintervalle
  B.refC.bootstrap <- function(data, indices){
    d <- data[indices,]
    B.refC.xboot <- model.matrix(YTRUE ~ A + B - 1, data = d)
    B.refC.yboot <- d$YTRUE
    B.refC.lasso.mod <- glmnet(B.refC.xboot, B.refC.yboot, alpha = 1, 
                               lambda = B.refC.best.lambda, standardize = FALSE)
    return(as.numeric(coef(B.refC.lasso.mod)))
  }
  
  B.refC.boot.results <- boot(data = data.dgp.case.b, statistic = B.refC.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  B.refC.ci.int <- boot.ci(B.refC.boot.results, type = "perc", index = 1)$percent[4:5]
  B.refC.ci.slope.one <- boot.ci(B.refC.boot.results, type = "perc", index = 2)$percent[4:5]
  B.refC.ci.slope.two <- boot.ci(B.refC.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  B.refC.cis.int.lasso <- B.refC.ci.int[1] <= betas[3] & B.refC.ci.int[2] >= betas[3]
  B.refC.cis.slope.one.lasso <- B.refC.ci.slope.one[1] <= betas[4] & B.refC.ci.slope.one[2] >= betas[4]
  B.refC.cis.slope.two.lasso <- B.refC.ci.slope.two[1] <= betas[5] & B.refC.ci.slope.two[2] >= betas[5]
  
  # Konfidenzintervalle numerisch abspeichern
  B.refC.ci.int.lower <- B.refC.ci.int[1]
  B.refC.ci.int.upper <- B.refC.ci.int[2]
  B.refC.ci.slope.one.lower <- B.refC.ci.slope.one[1]
  B.refC.ci.slope.one.upper <- B.refC.ci.slope.one[2]
  B.refC.ci.slope.two.lower <- B.refC.ci.slope.two[1]
  B.refC.ci.slope.two.upper <- B.refC.ci.slope.two[2]
  
  # Ergebnisdatensatz mit den numerischen Konfidenzintervallen
  B.refC.data.lasso <- data.frame(CASE.ID = "B",
                                  REF.ID = "REFC",
                                  TRUE.INT = betas[3],
                                  EST.INT = B.refC.est.int.lasso,
                                  BOOT.CI.INT = B.refC.cis.int.lasso,
                                  BOOT.CI.INT.LOWER = B.refC.ci.int.lower,
                                  BOOT.CI.INT.UPPER = B.refC.ci.int.upper,
                                  TRUE.SLOPE.ONE = betas[4],
                                  EST.SLOPE.ONE = B.refC.est.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE = B.refC.cis.slope.one.lasso,
                                  BOOT.CI.SLOPE.ONE.LOWER = B.refC.ci.slope.one.lower,
                                  BOOT.CI.SLOPE.ONE.UPPER = B.refC.ci.slope.one.upper,
                                  TRUE.SLOPE.TWO = betas[5],
                                  EST.SLOPE.TWO = B.refC.est.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO = B.refC.cis.slope.two.lasso,
                                  BOOT.CI.SLOPE.TWO.LOWER = B.refC.ci.slope.two.lower,
                                  BOOT.CI.SLOPE.TWO.UPPER = B.refC.ci.slope.two.upper,
                                  BEST.LAMBDA = B.refC.best.lambda,
                                  METHOD = "Lasso")
  
  # FINAL DATA FRAME ----
  B.data.all <- rbind(B.refA.data.lasso, B.refB.data.lasso, B.refC.data.lasso)
  
  return(B.data.all)
}

# set.seed(456)
# data.one.sim.lassoreg.case.b <- B.one.sim.lasso(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.lassoreg.case.b, file = "3.2_RESULTS_ONESIM_LASSO_B.RData")


# =========================== CASE C - LASSO REGRESSION ================== ####
C.one.sim.lasso <- function(n.total, betas){
  
  data.dgp.case.c <- dgp.case.c(n.total = n.total, betas = betas)
  
  # Subsets fuer die unterschiedlichen Subfaelle (C1, C2) von Fall C
  data.dgp.case.cone <- subset(data.dgp.case.c, CASE.ID == "C1")
  data.dgp.case.ctwo <- subset(data.dgp.case.c, CASE.ID == "C2")
  
  # SUBCASE C1 ----
  ## C1.A ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.cone)
  y <- data.dgp.case.cone$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz A (mit Achsenabschnitt)
  Cone.refA.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  Cone.refA.best.lambda <- Cone.refA.cv.fit$lambda.min
  Cone.refA.lasso.fit <- glmnet(x, y, alpha = 1, lambda = Cone.refA.best.lambda,
                                standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  Cone.refA.coef <- as.numeric(coef(Cone.refA.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall C1, Referenz A
  Cone.refA.est.int.lasso <- Cone.refA.coef[1]  # Achsenabschnitt
  Cone.refA.est.slope.one.lasso <- Cone.refA.coef[2]  # Parameter fuer B
  Cone.refA.est.slope.two.lasso <- Cone.refA.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  Cone.refA.bootstrap <- function(data, indices){
    d <- data[indices,]
    Cone.refA.xboot <- model.matrix(YTRUE ~ B + C - 1, data = d)
    Cone.refA.yboot <- d$YTRUE
    Cone.refA.lasso.mod <- glmnet(Cone.refA.xboot, Cone.refA.yboot, alpha = 1, 
                                  lambda = Cone.refA.best.lambda, standardize = FALSE)
    return(as.numeric(coef(Cone.refA.lasso.mod)))
  }
  
  Cone.refA.boot.results <- boot(data = data.dgp.case.cone, statistic = Cone.refA.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  Cone.refA.ci.int <- boot.ci(Cone.refA.boot.results, type = "perc", index = 1)$percent[4:5]
  Cone.refA.ci.slope.one <- boot.ci(Cone.refA.boot.results, type = "perc", index = 2)$percent[4:5]
  Cone.refA.ci.slope.two <- boot.ci(Cone.refA.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  Cone.refA.cis.int.lasso <- Cone.refA.ci.int[1] <= betas[5] & Cone.refA.ci.int[2] >= betas[5]
  Cone.refA.cis.slope.one.lasso <- Cone.refA.ci.slope.one[1] <= betas[3] & Cone.refA.ci.slope.one[2] >= betas[3]
  Cone.refA.cis.slope.two.lasso <- Cone.refA.ci.slope.two[1] <= betas[2] & Cone.refA.ci.slope.two[2] >= betas[2]
  
  # Konfidenzintervalle numerisch abspeichern
  Cone.refA.ci.int.lower <- Cone.refA.ci.int[1]
  Cone.refA.ci.int.upper <- Cone.refA.ci.int[2]
  Cone.refA.ci.slope.one.lower <- Cone.refA.ci.slope.one[1]
  Cone.refA.ci.slope.one.upper <- Cone.refA.ci.slope.one[2]
  Cone.refA.ci.slope.two.lower <- Cone.refA.ci.slope.two[1]
  Cone.refA.ci.slope.two.upper <- Cone.refA.ci.slope.two[2]
  
  # Ergebnis speichern
  Cone.refA.data.lasso <- data.frame(CASE.ID = "C1",
                                     REF.ID = "REFA",
                                     TRUE.INT = betas[5],
                                     EST.INT = Cone.refA.est.int.lasso,
                                     BOOT.CI.INT = Cone.refA.cis.int.lasso,
                                     BOOT.CI.INT.LOWER = Cone.refA.ci.int.lower,
                                     BOOT.CI.INT.UPPER = Cone.refA.ci.int.upper,
                                     TRUE.SLOPE.ONE = betas[3],
                                     EST.SLOPE.ONE = Cone.refA.est.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE = Cone.refA.cis.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE.LOWER = Cone.refA.ci.slope.one.lower,
                                     BOOT.CI.SLOPE.ONE.UPPER = Cone.refA.ci.slope.one.upper,
                                     TRUE.SLOPE.TWO = betas[2],
                                     EST.SLOPE.TWO = Cone.refA.est.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO = Cone.refA.cis.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO.LOWER = Cone.refA.ci.slope.two.lower,
                                     BOOT.CI.SLOPE.TWO.UPPER = Cone.refA.ci.slope.two.upper,
                                     BEST.LAMBDA = Cone.refA.best.lambda,
                                     METHOD = "Lasso")
  
  ## C1.B ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.cone)
  y <- data.dgp.case.cone$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz B (mit Achsenabschnitt)
  Cone.refB.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  Cone.refB.best.lambda <- Cone.refB.cv.fit$lambda.min
  Cone.refB.lasso.fit <- glmnet(x, y, alpha = 1, lambda = Cone.refB.best.lambda,
                                standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  Cone.refB.coef <- as.numeric(coef(Cone.refB.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall C1, Referenz B
  Cone.refB.est.int.lasso <- Cone.refB.coef[1]  # Achsenabschnitt
  Cone.refB.est.slope.one.lasso <- Cone.refB.coef[2]  # Parameter fuer A
  Cone.refB.est.slope.two.lasso <- Cone.refB.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  Cone.refB.bootstrap <- function(data, indices){
    d <- data[indices,]
    Cone.refB.xboot <- model.matrix(YTRUE ~ A + C - 1, data = d)
    Cone.refB.yboot <- d$YTRUE
    Cone.refB.lasso.mod <- glmnet(Cone.refB.xboot, Cone.refB.yboot, alpha = 1, 
                                  lambda = Cone.refB.best.lambda, standardize = FALSE)
    return(as.numeric(coef(Cone.refB.lasso.mod)))
  }
  
  Cone.refB.boot.results <- boot(data = data.dgp.case.cone, statistic = Cone.refB.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  Cone.refB.ci.int <- boot.ci(Cone.refB.boot.results, type = "perc", index = 1)$percent[4:5]
  Cone.refB.ci.slope.one <- boot.ci(Cone.refB.boot.results, type = "perc", index = 2)$percent[4:5]
  Cone.refB.ci.slope.two <- boot.ci(Cone.refB.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  Cone.refB.cis.int.lasso <- Cone.refB.ci.int[1] <= betas[5] & Cone.refB.ci.int[2] >= betas[5]
  Cone.refB.cis.slope.one.lasso <- Cone.refB.ci.slope.one[1] <= betas[3] & Cone.refB.ci.slope.one[2] >= betas[3]
  Cone.refB.cis.slope.two.lasso <- Cone.refB.ci.slope.two[1] <= betas[2] & Cone.refB.ci.slope.two[2] >= betas[2]
  
  # Konfidenzintervalle numerisch abspeichern
  Cone.refB.ci.int.lower <- Cone.refB.ci.int[1]
  Cone.refB.ci.int.upper <- Cone.refB.ci.int[2]
  Cone.refB.ci.slope.one.lower <- Cone.refB.ci.slope.one[1]
  Cone.refB.ci.slope.one.upper <- Cone.refB.ci.slope.one[2]
  Cone.refB.ci.slope.two.lower <- Cone.refB.ci.slope.two[1]
  Cone.refB.ci.slope.two.upper <- Cone.refB.ci.slope.two[2]
  
  # Ergebnis speichern
  Cone.refB.data.lasso <- data.frame(CASE.ID = "C1",
                                     REF.ID = "REFB",
                                     TRUE.INT = betas[5],
                                     EST.INT = Cone.refB.est.int.lasso,
                                     BOOT.CI.INT = Cone.refB.cis.int.lasso,
                                     BOOT.CI.INT.LOWER = Cone.refB.ci.int.lower,
                                     BOOT.CI.INT.UPPER = Cone.refB.ci.int.upper,
                                     TRUE.SLOPE.ONE = betas[3],
                                     EST.SLOPE.ONE = Cone.refB.est.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE = Cone.refB.cis.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE.LOWER = Cone.refB.ci.slope.one.lower,
                                     BOOT.CI.SLOPE.ONE.UPPER = Cone.refB.ci.slope.one.upper,
                                     TRUE.SLOPE.TWO = betas[2],
                                     EST.SLOPE.TWO = Cone.refB.est.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO = Cone.refB.cis.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO.LOWER = Cone.refB.ci.slope.two.lower,
                                     BOOT.CI.SLOPE.TWO.UPPER = Cone.refB.ci.slope.two.upper,
                                     BEST.LAMBDA = Cone.refB.best.lambda,
                                     METHOD = "Lasso")
  
  ## C1.C ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.cone)
  y <- data.dgp.case.cone$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz C (mit Achsenabschnitt)
  Cone.refC.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  Cone.refC.best.lambda <- Cone.refC.cv.fit$lambda.min
  Cone.refC.lasso.fit <- glmnet(x, y, alpha = 1, lambda = Cone.refC.best.lambda,
                                standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  Cone.refC.coef <- as.numeric(coef(Cone.refC.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall C1, Referenz C
  Cone.refC.est.int.lasso <- Cone.refC.coef[1]  # Achsenabschnitt
  Cone.refC.est.slope.one.lasso <- Cone.refC.coef[2]  # Parameter fuer A
  Cone.refC.est.slope.two.lasso <- Cone.refC.coef[3]  # Parameter fuer B
  
  # Bootstrap fuer Konfidenzintervalle
  Cone.refC.bootstrap <- function(data, indices){
    d <- data[indices,]
    Cone.refC.xboot <- model.matrix(YTRUE ~ A + B - 1, data = d)
    Cone.refC.yboot <- d$YTRUE
    Cone.refC.lasso.mod <- glmnet(Cone.refC.xboot, Cone.refC.yboot, alpha = 1, 
                                  lambda = Cone.refC.best.lambda, standardize = FALSE)
    return(as.numeric(coef(Cone.refC.lasso.mod)))
  }
  
  Cone.refC.boot.results <- boot(data = data.dgp.case.cone, statistic = Cone.refC.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  Cone.refC.ci.int <- boot.ci(Cone.refC.boot.results, type = "perc", index = 1)$percent[4:5]
  Cone.refC.ci.slope.one <- boot.ci(Cone.refC.boot.results, type = "perc", index = 2)$percent[4:5]
  Cone.refC.ci.slope.two <- boot.ci(Cone.refC.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  Cone.refC.cis.int.lasso <- Cone.refC.ci.int[1] <= betas[4] & Cone.refC.ci.int[2] >= betas[4]
  Cone.refC.cis.slope.one.lasso <- Cone.refC.ci.slope.one[1] <= betas[4] & Cone.refC.ci.slope.one[2] >= betas[4]
  Cone.refC.cis.slope.two.lasso <- Cone.refC.ci.slope.two[1] <= betas[4] & Cone.refC.ci.slope.two[2] >= betas[4]
  
  # Konfidenzintervalle numerisch abspeichern
  Cone.refC.ci.int.lower <- Cone.refC.ci.int[1]
  Cone.refC.ci.int.upper <- Cone.refC.ci.int[2]
  Cone.refC.ci.slope.one.lower <- Cone.refC.ci.slope.one[1]
  Cone.refC.ci.slope.one.upper <- Cone.refC.ci.slope.one[2]
  Cone.refC.ci.slope.two.lower <- Cone.refC.ci.slope.two[1]
  Cone.refC.ci.slope.two.upper <- Cone.refC.ci.slope.two[2]
  
  # Ergebnis speichern
  Cone.refC.data.lasso <- data.frame(CASE.ID = "C1",
                                     REF.ID = "REFC",
                                     TRUE.INT = betas[4],
                                     EST.INT = Cone.refC.est.int.lasso,
                                     BOOT.CI.INT = Cone.refC.cis.int.lasso,
                                     BOOT.CI.INT.LOWER = Cone.refC.ci.int.lower,
                                     BOOT.CI.INT.UPPER = Cone.refC.ci.int.upper,
                                     TRUE.SLOPE.ONE = betas[4],
                                     EST.SLOPE.ONE = Cone.refC.est.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE = Cone.refC.cis.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE.LOWER = Cone.refC.ci.slope.one.lower,
                                     BOOT.CI.SLOPE.ONE.UPPER = Cone.refC.ci.slope.one.upper,
                                     TRUE.SLOPE.TWO = betas[4],
                                     EST.SLOPE.TWO = Cone.refC.est.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO = Cone.refC.cis.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO.LOWER = Cone.refC.ci.slope.two.lower,
                                     BOOT.CI.SLOPE.TWO.UPPER = Cone.refC.ci.slope.two.upper,
                                     BEST.LAMBDA = Cone.refC.best.lambda,
                                     METHOD = "Lasso")
  
  # SUBCASE C2 ----
  ## C2.A ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ B + C - 1, data = data.dgp.case.ctwo)
  y <- data.dgp.case.ctwo$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz A (mit Achsenabschnitt)
  Ctwo.refA.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  Ctwo.refA.best.lambda <- Ctwo.refA.cv.fit$lambda.min
  Ctwo.refA.lasso.fit <- glmnet(x, y, alpha = 1, lambda = Ctwo.refA.best.lambda, 
                                standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  Ctwo.refA.coef <- as.numeric(coef(Ctwo.refA.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall C2, Referenz A
  Ctwo.refA.est.int.lasso <- Ctwo.refA.coef[1]  # Achsenabschnitt
  Ctwo.refA.est.slope.one.lasso <- Ctwo.refA.coef[2]  # Parameter fuer B
  Ctwo.refA.est.slope.two.lasso <- Ctwo.refA.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  Ctwo.refA.bootstrap <- function(data, indices){
    d <- data[indices,]
    Ctwo.refA.xboot <- model.matrix(YTRUE ~ B + C - 1, data = d)
    Ctwo.refA.yboot <- d$YTRUE
    Ctwo.refA.lasso.mod <- glmnet(Ctwo.refA.xboot, Ctwo.refA.yboot, alpha = 1, 
                                  lambda = Ctwo.refA.best.lambda, standardize = FALSE)
    return(as.numeric(coef(Ctwo.refA.lasso.mod)))
  }
  
  Ctwo.refA.boot.results <- boot(data = data.dgp.case.ctwo, statistic = Ctwo.refA.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  Ctwo.refA.ci.int <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 1)$percent[4:5]
  Ctwo.refA.ci.slope.one <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 2)$percent[4:5]
  Ctwo.refA.ci.slope.two <- boot.ci(Ctwo.refA.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  Ctwo.refA.cis.int.lasso <- Ctwo.refA.ci.int[1] <= betas[5] & Ctwo.refA.ci.int[2] >= betas[5]
  Ctwo.refA.cis.slope.one.lasso <- Ctwo.refA.ci.slope.one[1] <= betas[3] & Ctwo.refA.ci.slope.one[2] >= betas[3]
  Ctwo.refA.cis.slope.two.lasso <- Ctwo.refA.ci.slope.two[1] <= betas[1] & Ctwo.refA.ci.slope.two[2] >= betas[1]
  
  # Konfidenzintervalle numerisch abspeichern
  Ctwo.refA.ci.int.lower <- Ctwo.refA.ci.int[1]
  Ctwo.refA.ci.int.upper <- Ctwo.refA.ci.int[2]
  Ctwo.refA.ci.slope.one.lower <- Ctwo.refA.ci.slope.one[1]
  Ctwo.refA.ci.slope.one.upper <- Ctwo.refA.ci.slope.one[2]
  Ctwo.refA.ci.slope.two.lower <- Ctwo.refA.ci.slope.two[1]
  Ctwo.refA.ci.slope.two.upper <- Ctwo.refA.ci.slope.two[2]
  
  # Ergebnis speichern
  Ctwo.refA.data.lasso <- data.frame(CASE.ID = "C2",
                                     REF.ID = "REFA",
                                     TRUE.INT = betas[5],
                                     EST.INT = Ctwo.refA.est.int.lasso,
                                     BOOT.CI.INT = Ctwo.refA.cis.int.lasso,
                                     BOOT.CI.INT.LOWER = Ctwo.refA.ci.int.lower,
                                     BOOT.CI.INT.UPPER = Ctwo.refA.ci.int.upper,
                                     TRUE.SLOPE.ONE = betas[3],
                                     EST.SLOPE.ONE = Ctwo.refA.est.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE = Ctwo.refA.cis.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE.LOWER = Ctwo.refA.ci.slope.one.lower,
                                     BOOT.CI.SLOPE.ONE.UPPER = Ctwo.refA.ci.slope.one.upper,
                                     TRUE.SLOPE.TWO = betas[1],
                                     EST.SLOPE.TWO = Ctwo.refA.est.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO = Ctwo.refA.cis.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO.LOWER = Ctwo.refA.ci.slope.two.lower,
                                     BOOT.CI.SLOPE.TWO.UPPER = Ctwo.refA.ci.slope.two.upper,
                                     BEST.LAMBDA = Ctwo.refA.best.lambda,
                                     METHOD = "Lasso")
  
  ## C2.B ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + C - 1, data = data.dgp.case.ctwo)
  y <- data.dgp.case.ctwo$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz B (mit Achsenabschnitt)
  Ctwo.refB.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  Ctwo.refB.best.lambda <- Ctwo.refB.cv.fit$lambda.min
  Ctwo.refB.lasso.fit <- glmnet(x, y, alpha = 1, lambda = Ctwo.refB.best.lambda, 
                                standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  Ctwo.refB.coef <- as.numeric(coef(Ctwo.refB.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall C2, Referenz B
  Ctwo.refB.est.int.lasso <- Ctwo.refB.coef[1]  # Achsenabschnitt
  Ctwo.refB.est.slope.one.lasso <- Ctwo.refB.coef[2]  # Parameter fuer A
  Ctwo.refB.est.slope.two.lasso <- Ctwo.refB.coef[3]  # Parameter fuer C
  
  # Bootstrap fuer Konfidenzintervalle
  Ctwo.refB.bootstrap <- function(data, indices){
    d <- data[indices,]
    Ctwo.refB.xboot <- model.matrix(YTRUE ~ A + C - 1, data = d)
    Ctwo.refB.yboot <- d$YTRUE
    Ctwo.refB.lasso.mod <- glmnet(Ctwo.refB.xboot, Ctwo.refB.yboot, alpha = 1, 
                                  lambda = Ctwo.refB.best.lambda, standardize = FALSE)
    return(as.numeric(coef(Ctwo.refB.lasso.mod)))
  }
  
  Ctwo.refB.boot.results <- boot(data = data.dgp.case.ctwo, statistic = Ctwo.refB.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  Ctwo.refB.ci.int <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 1)$percent[4:5]
  Ctwo.refB.ci.slope.one <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 2)$percent[4:5]
  Ctwo.refB.ci.slope.two <- boot.ci(Ctwo.refB.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  Ctwo.refB.cis.int.lasso <- Ctwo.refB.ci.int[1] <= betas[5] & Ctwo.refB.ci.int[2] >= betas[5]
  Ctwo.refB.cis.slope.one.lasso <- Ctwo.refB.ci.slope.one[1] <= betas[3] & Ctwo.refB.ci.slope.one[2] >= betas[3]
  Ctwo.refB.cis.slope.two.lasso <- Ctwo.refB.ci.slope.two[1] <= betas[1] & Ctwo.refB.ci.slope.two[2] >= betas[1]
  
  # Konfidenzintervalle numerisch abspeichern
  Ctwo.refB.ci.int.lower <- Ctwo.refB.ci.int[1]
  Ctwo.refB.ci.int.upper <- Ctwo.refB.ci.int[2]
  Ctwo.refB.ci.slope.one.lower <- Ctwo.refB.ci.slope.one[1]
  Ctwo.refB.ci.slope.one.upper <- Ctwo.refB.ci.slope.one[2]
  Ctwo.refB.ci.slope.two.lower <- Ctwo.refB.ci.slope.two[1]
  Ctwo.refB.ci.slope.two.upper <- Ctwo.refB.ci.slope.two[2]
  
  # Ergebnis-Datenframe fuer Referenz B in Fall C2
  Ctwo.refB.data.lasso <- data.frame(CASE.ID = "C2",
                                     REF.ID = "REFB",
                                     TRUE.INT = betas[5],
                                     EST.INT = Ctwo.refB.est.int.lasso,
                                     BOOT.CI.INT = Ctwo.refB.cis.int.lasso,
                                     BOOT.CI.INT.LOWER = Ctwo.refB.ci.int.lower,
                                     BOOT.CI.INT.UPPER = Ctwo.refB.ci.int.upper,
                                     TRUE.SLOPE.ONE = betas[3],
                                     EST.SLOPE.ONE = Ctwo.refB.est.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE = Ctwo.refB.cis.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE.LOWER = Ctwo.refB.ci.slope.one.lower,
                                     BOOT.CI.SLOPE.ONE.UPPER = Ctwo.refB.ci.slope.one.upper,
                                     TRUE.SLOPE.TWO = betas[1],
                                     EST.SLOPE.TWO = Ctwo.refB.est.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO = Ctwo.refB.cis.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO.LOWER = Ctwo.refB.ci.slope.two.lower,
                                     BOOT.CI.SLOPE.TWO.UPPER = Ctwo.refB.ci.slope.two.upper,
                                     BEST.LAMBDA = Ctwo.refB.best.lambda,
                                     METHOD = "Lasso")
  
  ## C2.C ----
  # Vorbereitung der Daten fuer Lasso Regression
  x <- model.matrix(YTRUE ~ A + B - 1, data = data.dgp.case.ctwo)
  y <- data.dgp.case.ctwo$YTRUE
  
  # Lasso Regression mit Kreuzvalidierung fuer Referenz C (mit Achsenabschnitt)
  Ctwo.refC.cv.fit <- cv.glmnet(x, y, alpha = 1, standardize = FALSE)
  Ctwo.refC.best.lambda <- Ctwo.refC.cv.fit$lambda.min
  Ctwo.refC.lasso.fit <- glmnet(x, y, alpha = 1, lambda = Ctwo.refC.best.lambda, 
                                standardize = FALSE)
  
  # Extraktion der Regressionsparametern fuer das lambda, das den minimalen Fehler liefert
  Ctwo.refC.coef <- as.numeric(coef(Ctwo.refC.lasso.fit, s = "lambda.min"))
  
  # Geschaetzte Regressionsparametern fuer Fall C2, Referenz C
  Ctwo.refC.est.int.lasso <- Ctwo.refC.coef[1]  # Achsenabschnitt
  Ctwo.refC.est.slope.one.lasso <- Ctwo.refC.coef[2]  # Parameter fuer A
  Ctwo.refC.est.slope.two.lasso <- Ctwo.refC.coef[3]  # Parameter fuer B
  
  # Bootstrap fuer Konfidenzintervalle
  Ctwo.refC.bootstrap <- function(data, indices){
    d <- data[indices,]
    Ctwo.refC.xboot <- model.matrix(YTRUE ~ A + B - 1, data = d)
    Ctwo.refC.yboot <- d$YTRUE
    Ctwo.refC.lasso.mod <- glmnet(Ctwo.refC.xboot, Ctwo.refC.yboot, alpha = 1, 
                                  lambda = Ctwo.refC.best.lambda, standardize = FALSE)
    return(as.numeric(coef(Ctwo.refC.lasso.mod)))
  }
  
  Ctwo.refC.boot.results <- boot(data = data.dgp.case.ctwo, statistic = Ctwo.refC.bootstrap, R = 1000)
  
  # Konfidenzintervalle berechnen
  Ctwo.refC.ci.int <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 1)$percent[4:5]
  Ctwo.refC.ci.slope.one <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 2)$percent[4:5]
  Ctwo.refC.ci.slope.two <- boot.ci(Ctwo.refC.boot.results, type = "perc", index = 3)$percent[4:5]
  
  # Ueberpruefen, ob die wahren Werte innerhalb der Konfidenzintervalle liegen
  Ctwo.refC.cis.int.lasso <- Ctwo.refC.ci.int[1] <= betas[3] & Ctwo.refC.ci.int[2] >= betas[3]
  Ctwo.refC.cis.slope.one.lasso <- Ctwo.refC.ci.slope.one[1] <= betas[5] & Ctwo.refC.ci.slope.one[2] >= betas[5]
  Ctwo.refC.cis.slope.two.lasso <- Ctwo.refC.ci.slope.two[1] <= betas[5] & Ctwo.refC.ci.slope.two[2] >= betas[5]
  
  # Konfidenzintervalle numerisch abspeichern
  Ctwo.refC.ci.int.lower <- Ctwo.refC.ci.int[1]
  Ctwo.refC.ci.int.upper <- Ctwo.refC.ci.int[2]
  Ctwo.refC.ci.slope.one.lower <- Ctwo.refC.ci.slope.one[1]
  Ctwo.refC.ci.slope.one.upper <- Ctwo.refC.ci.slope.one[2]
  Ctwo.refC.ci.slope.two.lower <- Ctwo.refC.ci.slope.two[1]
  Ctwo.refC.ci.slope.two.upper <- Ctwo.refC.ci.slope.two[2]
  
  # Ergebnis speichern
  Ctwo.refC.data.lasso <- data.frame(CASE.ID = "C2",
                                     REF.ID = "REFC",
                                     TRUE.INT = betas[3],
                                     EST.INT = Ctwo.refC.est.int.lasso,
                                     BOOT.CI.INT = Ctwo.refC.cis.int.lasso,
                                     BOOT.CI.INT.LOWER = Ctwo.refC.ci.int.lower,
                                     BOOT.CI.INT.UPPER = Ctwo.refC.ci.int.upper,
                                     TRUE.SLOPE.ONE = betas[5],
                                     EST.SLOPE.ONE = Ctwo.refC.est.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE = Ctwo.refC.cis.slope.one.lasso,
                                     BOOT.CI.SLOPE.ONE.LOWER = Ctwo.refC.ci.slope.one.lower,
                                     BOOT.CI.SLOPE.ONE.UPPER = Ctwo.refC.ci.slope.one.upper,
                                     TRUE.SLOPE.TWO = betas[5],
                                     EST.SLOPE.TWO = Ctwo.refC.est.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO = Ctwo.refC.cis.slope.two.lasso,
                                     BOOT.CI.SLOPE.TWO.LOWER = Ctwo.refC.ci.slope.two.lower,
                                     BOOT.CI.SLOPE.TWO.UPPER = Ctwo.refC.ci.slope.two.upper,
                                     BEST.LAMBDA = Ctwo.refC.best.lambda,
                                     METHOD = "Lasso")
  
  # FINAL DATA FRAME ----
  C.data.all <- rbind(Cone.refA.data.lasso, Cone.refB.data.lasso, Cone.refC.data.lasso,
                      Ctwo.refA.data.lasso, Ctwo.refB.data.lasso, Ctwo.refC.data.lasso)
  
  return(C.data.all)
}

# set.seed(789)
# data.one.sim.lassoreg.case.c <- C.one.sim.lasso(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.lassoreg.case.c, file = "3.3_RESULTS_ONESIM_LASSO_C.RData")

# ========================= ALL CASES - LASSO REGRESSION ================ ####

one.sim.lasso.all.cases <- function(n.total, betas){
  
  A.one.sim.data <- A.one.sim.lasso(n.total = n.total, betas = betas)
  B.one.sim.data <- B.one.sim.lasso(n.total = n.total, betas = betas)
  C.one.sim.data <- C.one.sim.lasso(n.total = n.total, betas = betas)
  
  all.cases.one.sim.data <- rbind(A.one.sim.data, 
                                  B.one.sim.data, 
                                  C.one.sim.data)
  
  return(all.cases.one.sim.data)
}

# set.seed(123)
# data.one.sim.lasso.all.cases <- one.sim.lasso.all.cases(n.total = 100, betas = c(-10, -5, 0, 5, 10))
# save(data.one.sim.lasso.all.cases, file = "3.4_RESULTS_ONESIM_LASSO_ALL_CASES.RData")

# ============================================================================ #
# ENDE DES SKRIPTS -----------------------------------------------------------
# ============================================================================ #