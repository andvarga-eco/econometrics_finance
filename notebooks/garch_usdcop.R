paquetes <- c("quantmod", "rugarch", "tseries", "FinTS")
faltantes <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]
if (length(faltantes) > 0) install.packages(faltantes)

invisible(lapply(paquetes, library, character.only = TRUE))

# ---------------------------------------------------------------------------
# 2. DESCARGA DE DATOS DESDE YAHOO FINANCE
# ---------------------------------------------------------------------------
# Ticker del petróleo Brent (ICE Brent Crude Oil futures) en Yahoo Finance: "BZ=F"

ticker    <- "usdcop=x"
fecha_fin <- Sys.Date()
fecha_ini <- fecha_fin - 5 * 365  # aprox. 5 años de historia

usdcop<- getSymbols(
  Symbols     = ticker,
  src         = "yahoo",
  from        = fecha_ini,
  to          = fecha_fin,
  auto.assign = FALSE
)

# Nos quedamos con el precio de cierre y eliminamos NAs (días sin cotización)
precio <- na.omit(Cl(usdcop))
colnames(precio) <- "Cierre"

plot(precio, main = "Precio de cierre - USDCOP", col = "steelblue")

# ---------------------------------------------------------------------------
# 3. RETORNOS Y ANÁLISIS EXPLORATORIO
# ---------------------------------------------------------------------------
# Log-retornos en porcentaje (estabiliza la escala numérica para la optimización)

retornos <- na.omit(diff(log(precio)) * 100)
colnames(retornos) <- "retorno"

plot(retornos, main = "Retornos diarios - USDCOP", col = "darkred")
hist(retornos, breaks = 60, main = "Distribución de los retornos", col = "gray80",
     xlab = "Retorno (%)")

# Estadísticos descriptivos
summary(retornos)
cat("SD:",sd(retornos),"\n")
skewness <- function(x) {
  x <- x - mean(x)
  (mean(x^3)) / (mean(x^2))^1.5
}
cat("Asimetría:", skewness(as.numeric(retornos)), "\n")
moments::kurtosis(retornos)

# Estacionariedad de la serie de retornos (Dickey-Fuller aumentado)
adf.test(retornos)

# Autocorrelación de los retornos y de los retornos al cuadrado
par(mfrow = c(2, 2))
acf(retornos,      main = "ACF retornos")
pacf(retornos,     main = "PACF retornos")
acf(retornos^2,    main = "ACF retornos^2")
pacf(retornos^2,   main = "PACF retornos^2")
par(mfrow = c(1, 1))

# Aparente efecto día de la semana. Dummies por día

wday_num<-.indexwday(retornos)

days_matrix <- model.matrix(~ factor(wday_num) - 1)
# Name the columns appropriately
colnames(days_matrix) <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri")[unique(wday_num) + 1]

days_matrix<-days_matrix[,-4]

# ---------------------------------------------------------------------------
# 4. PRUEBAS DE EFECTOS ARCH (justifican el uso de un modelo GARCH)
# ---------------------------------------------------------------------------

# Ljung-Box sobre los retornos al cuadrado: H0 = no hay autocorrelación
# (si se rechaza, hay evidencia de heterocedasticidad condicional -> ARCH)
Box.test(retornos^2, lag = 12, type = "Ljung-Box")

# Prueba ARCH-LM de Engle
FinTS::ArchTest(retornos, lags = 12)

# ---------------------------------------------------------------------------
# 5. ESPECIFICACIÓN Y ESTIMACIÓN DEL MODELO GARCH(1,1)
# ---------------------------------------------------------------------------
# Media: ARMA(0,0) con constante | Varianza: GARCH(1,1)
# Distribución: t-Student estandarizada (los retornos financieros suelen
# tener colas más pesadas que la normal)

spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

fit_garch <- ugarchfit(spec = spec_garch, data = retornos, solver = "hybrid")

show(fit_garch)          # resumen completo: coeficientes, criterios de info, pruebas
coef(fit_garch)          # coeficientes estimados
infocriteria(fit_garch)  # AIC, BIC, etc.

coef(fit_garch)["alpha1"]+coef(fit_garch)["beta1"]
coef(fit_garch)["omega"]/(1-(coef(fit_garch)["alpha1"]+coef(fit_garch)["beta1"]))
# ---------------------------------------------------------------------------
# 6. DIAGNÓSTICO DEL MODELO
# ---------------------------------------------------------------------------
# Se evalúa si, tras estimar el GARCH, los residuales estandarizados se
# comportan como ruido blanco y si desaparece la heterocedasticidad
# condicional (si no desaparece, el modelo no capturó bien la dinámica).

res_std <- residuals(fit_garch, standardize = TRUE)

# 6.1 Autocorrelación remanente en los residuales estandarizados
Box.test(res_std,    lag = 10, type = "Ljung-Box")
# 6.2 Autocorrelación remanente en los residuales al cuadrado (efectos ARCH residuales)
Box.test(res_std^2,  lag = 10, type = "Ljung-Box")
# 6.3 Prueba ARCH-LM sobre los residuales estandarizados
FinTS::ArchTest(res_std, lags = 10)

# 6.4 Normalidad / bondad de ajuste de la distribución supuesta
jarque.bera.test(as.numeric(res_std))
qqnorm(as.numeric(res_std), main = "QQ-plot de residuales estandarizados")
qqline(as.numeric(res_std), col = "red")

# 6.5 Prueba de sesgo de signo (Engle-Ng): detecta asimetrías/efecto apalancamiento
#     no capturadas por un GARCH simétrico
signbias(fit_garch)

# 6.6 Estabilidad de los parámetros en el tiempo (test de Nyblom)
nyblom(fit_garch)

# 6.7 Panel de gráficos de diagnóstico integrados de rugarch
#     (series con volatilidad, ACF de residuales, QQ-plot, News Impact Curve, etc.)
# En modo interactivo (RStudio) recorre las opciones con Enter:
plot(fit_garch, which = "all")
# Para un gráfico puntual, usar por ejemplo:
# plot(fit_garch, which = 9)   # QQ-plot de residuales estandarizados
# plot(fit_garch, which = 11)  # ACF de residuales estandarizados al cuadrado

vol<-sigma(fit_garch)
plot(vol)
# ---------------------------------------------------------------------------
# 7. PRONÓSTICO DE VOLATILIDAD
# ---------------------------------------------------------------------------
tail(sigma(fit_garch))
pronostico <- ugarchforecast(fit_garch, n.ahead = 30)
pronostico
plot(pronostico, which = 1)  # pronóstico de la serie (retornos)
plot(pronostico, which = 3)  # pronóstico de la volatilidad (sigma) condicional

# Las pruebas de residuales sugieren que aún hay correlación de orden 5 en la ecuación de la media
# Incorporamos regresores externos

spec_garch_dow <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE, external.regressors=days_matrix),
  distribution.model = "std"
)

fit_garch_dow <- ugarchfit(spec = spec_garch_dow, data = retornos, solver = "hybrid")

show(fit_garch_dow)          # resumen completo: coeficientes, criterios de info, pruebas
coef(fit_garch_dow)          # coeficientes estimados
infocriteria(fit_garch_dow)  # AIC, BIC, etc.

coef(fit_garch_dow)["alpha1"]+coef(fit_garch_dow)["beta1"]
coef(fit_garch_dow)["omega"]/(1-(coef(fit_garch_dow)["alpha1"]+coef(fit_garch_dow)["beta1"]))

res_std_dow <- residuals(fit_garch_dow, standardize = TRUE)

Box.test(res_std_dow,    lag = 10, type = "Ljung-Box")
Box.test(res_std_dow^2,  lag = 10, type = "Ljung-Box")
FinTS::ArchTest(res_std_dow, lags = 10)

jarque.bera.test(as.numeric(res_std_dow))
qqnorm(as.numeric(res_std_dow), main = "QQ-plot de residuales estandarizados")
qqline(as.numeric(res_std_dow), col = "red")

signbias(fit_garch_dow)
nyblom(fit_garch_dow)

plot(fit_garch_dow, which = "all")

vol_dow<-sigma(fit_garch_dow)
plot(vol_dow)
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
tail(sigma(fit_garch_dow))
pronostico_dow <- ugarchforecast(fit_garch_dow, n.ahead = 30)
pronostico_dow
plot(pronostico_dow, which = 1)  # pronóstico de la serie (retornos)
plot(pronostico_dow, which = 3)  # pronóstico de la volatilidad (sigma) condicional
