#############################################################################
# ESTIMACIÓN DE UN MODELO GARCH PARA EL PRECIO DEL PETRÓLEO BRENT
# Datos diarios, últimos 5 años, descargados directamente de Yahoo Finance
#
# Flujo de trabajo:
#   1. Paquetes y descarga de datos
#   2. Construcción de retornos y análisis exploratorio
#   3. Pruebas de efectos ARCH (justificación del modelo GARCH)
#   4. Especificación y estimación del modelo GARCH(1,1)
#   5. Diagnóstico del modelo (residuales estandarizados)
#   6. Pronóstico de volatilidad
#   7. (Opcional) Modelos asimétricos GJR-GARCH / EGARCH y comparación
#############################################################################

# ---------------------------------------------------------------------------
# 1. PAQUETES
# ---------------------------------------------------------------------------

paquetes <- c("quantmod", "rugarch", "tseries", "FinTS")
faltantes <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]
if (length(faltantes) > 0) install.packages(faltantes)

invisible(lapply(paquetes, library, character.only = TRUE))

# ---------------------------------------------------------------------------
# 2. DESCARGA DE DATOS DESDE YAHOO FINANCE
# ---------------------------------------------------------------------------
# Ticker del petróleo Brent (ICE Brent Crude Oil futures) en Yahoo Finance: "BZ=F"

ticker    <- "BZ=F"
fecha_fin <- Sys.Date()
fecha_ini <- fecha_fin - 5 * 365  # aprox. 5 años de historia

brent <- getSymbols(
  Symbols     = ticker,
  src         = "yahoo",
  from        = fecha_ini,
  to          = fecha_fin,
  auto.assign = FALSE
)

# Nos quedamos con el precio de cierre y eliminamos NAs (días sin cotización)
precio <- na.omit(Cl(brent))
colnames(precio) <- "Cierre"

plot(precio, main = "Precio de cierre - Petróleo Brent (BZ=F)", col = "steelblue")

# ---------------------------------------------------------------------------
# 3. RETORNOS Y ANÁLISIS EXPLORATORIO
# ---------------------------------------------------------------------------
# Log-retornos en porcentaje (estabiliza la escala numérica para la optimización)

retornos <- na.omit(diff(log(precio)) * 100)
colnames(retornos) <- "retorno"

plot(retornos, main = "Retornos diarios - Brent", col = "darkred")
hist(retornos, breaks = 60, main = "Distribución de los retornos", col = "gray80",
     xlab = "Retorno (%)")

# Estadísticos descriptivos
summary(retornos)
sd(retornos)
skewness <- function(x) {
  x <- x - mean(x)
  (mean(x^3)) / (mean(x^2))^1.5
}
cat("Asimetría:", skewness(as.numeric(retornos)), "\n")

# Estacionariedad de la serie de retornos (Dickey-Fuller aumentado)
adf.test(retornos)

# Autocorrelación de los retornos y de los retornos al cuadrado
par(mfrow = c(2, 2))
acf(retornos,      main = "ACF retornos")
pacf(retornos,     main = "PACF retornos")
acf(retornos^2,    main = "ACF retornos^2")
pacf(retornos^2,   main = "PACF retornos^2")
par(mfrow = c(1, 1))

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

pronostico <- ugarchforecast(fit_garch, n.ahead = 20)
pronostico
plot(pronostico, which = 1)  # pronóstico de la serie (retornos)
plot(pronostico, which = 3)  # pronóstico de la volatilidad (sigma) condicional

# ---------------------------------------------------------------------------
# 8. (OPCIONAL) MODELOS ASIMÉTRICOS Y COMPARACIÓN
# ---------------------------------------------------------------------------
# El precio del petróleo suele mostrar asimetría en la volatilidad (choques
# negativos de precio elevan la volatilidad más que choques positivos).
# Si la prueba de sesgo de signo (paso 6.5) resulta significativa, conviene
# comparar el GARCH simétrico con especificaciones asimétricas.

spec_gjr <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
fit_gjr <- ugarchfit(spec = spec_gjr, data = retornos, solver = "hybrid")

spec_egarch <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
fit_egarch <- ugarchfit(spec = spec_egarch, data = retornos, solver = "hybrid")

comparacion <- data.frame(
  modelo = c("sGARCH(1,1)", "gjrGARCH(1,1)", "eGARCH(1,1)"),
  AIC = c(infocriteria(fit_garch)[1], infocriteria(fit_gjr)[1], infocriteria(fit_egarch)[1]),
  BIC = c(infocriteria(fit_garch)[2], infocriteria(fit_gjr)[2], infocriteria(fit_egarch)[2])
)
comparacion  # el modelo con menor AIC/BIC es preferible
