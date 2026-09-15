library(quantmod)
library(rugarch)
library(tseries)
library(FinTS)

##MXNUSD

ticker    <- "MXN=X"
fecha_fin <- Sys.Date()
fecha_ini <- fecha_fin - 5 * 365  # aprox. 5 años de historia

mxn <- getSymbols(
  Symbols     = ticker,
  src         = "yahoo",
  from        = fecha_ini,
  to          = fecha_fin,
  auto.assign = FALSE
)

# Nos quedamos con el precio de cierre y eliminamos NAs (días sin cotización)
precio <- na.omit(Cl(mxn))
colnames(precio) <- "Cierre"
retornos_mxn <- na.omit(diff(log(precio)) * 100)
colnames(retornos_mxn) <- "retorno"

plot(retornos, main = "Retornos diarios - Brent", col = "darkred")

Box.test(retornos^2, lag = 12, type = "Ljung-Box")

# Prueba ARCH-LM de Engle
FinTS::ArchTest(retornos, lags = 12)

spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

fit_garch_mxn <- ugarchfit(spec = spec_garch, data = retornos_mxn, solver = "hybrid")

show(fit_garch_mxn)          # resumen completo: coeficientes, criterios de info, pruebas
coef(fit_garch_mxn)          # coeficientes estimados
infocriteria(fit_garch_mxn)  # AIC, BIC, etc.

coef(fit_garch_mxn)["alpha1"]+coef(fit_garch_mxn)["beta1"]
coef(fit_garch_mxn)["omega"]/(1-(coef(fit_garch_mxn)["alpha1"]+coef(fit_garch_mxn)["beta1"]))

# Samsung

ticker    <- "005930.KS"
fecha_fin <- Sys.Date()
fecha_ini <- fecha_fin - 5 * 365  # aprox. 5 años de historia

sa <- getSymbols(
  Symbols     = ticker,
  src         = "yahoo",
  from        = fecha_ini,
  to          = fecha_fin,
  auto.assign = FALSE
)

# Nos quedamos con el precio de cierre y eliminamos NAs (días sin cotización)
precio <- na.omit(Cl(sa))
colnames(precio) <- "Cierre"
retornos_sa <- na.omit(diff(log(precio)) * 100)
colnames(retornos_sa) <- "retorno"

plot(retornos_sa, main = "Retornos diarios - sa", col = "darkred")

Box.test(retornos_sa^2, lag = 12, type = "Ljung-Box")

# Prueba ARCH-LM de Engle
FinTS::ArchTest(retornos_sa, lags = 12)

spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

fit_garch_sa <- ugarchfit(spec = spec_garch, data = retornos_sa, solver = "hybrid")

show(fit_garch_sa)          # resumen completo: coeficientes, criterios de info, pruebas
coef(fit_garch_sa)          # coeficientes estimados
infocriteria(fit_garch_sa)  # AIC, BIC, etc.


# café
ticker    <- "KC=F"
fecha_fin <- Sys.Date()
fecha_ini <- fecha_fin - 5 * 365  # aprox. 5 años de historia

caf <- getSymbols(
  Symbols     = ticker,
  src         = "yahoo",
  from        = fecha_ini,
  to          = fecha_fin,
  auto.assign = FALSE
)

# Nos quedamos con el precio de cierre y eliminamos NAs (días sin cotización)
precio <- na.omit(Cl(caf))
colnames(precio) <- "Cierre"
retornos_caf <- na.omit(diff(log(precio)) * 100)
colnames(retornos_caf) <- "retorno"

plot(retornos_caf, main = "Retornos diarios - caf", col = "darkred")

Box.test(retornos_caf^2, lag = 12, type = "Ljung-Box")

# Prueba ARCH-LM de Engle
FinTS::ArchTest(retornos_caf, lags = 12)

spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

fit_garch_caf <- ugarchfit(spec = spec_garch, data = retornos_caf, solver = "hybrid")

show(fit_garch_caf)          # resumen completo: coeficientes, criterios de info, pruebas
coef(fit_garch_caf)          # coeficientes estimados
infocriteria(fit_garch_caf)  # AIC, BIC, etc.


