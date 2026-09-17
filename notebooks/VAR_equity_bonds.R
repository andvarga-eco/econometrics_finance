################################################################################
# VAR y cobertura dinámica para un portafolio equity-renta fija en Colombia
#
# Réplica en R del ejercicio: COLCAP (equity) + TES (renta fija) [+ TRM opcional]
#
# Autor: (adaptar)
# Requiere: quantmod, vars, urca, tseries, dplyr, readr, zoo, ggplot2
################################################################################

# ---- 0. Paquetes ------------------------------------------------------------
paquetes <- c("quantmod", "vars", "urca", "tseries", "dplyr",
              "readr", "zoo", "ggplot2", "tidyr")
instalar <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]
if (length(instalar) > 0) install.packages(instalar)
invisible(lapply(paquetes, library, character.only = TRUE))

# ==============================================================================
# 1. DESCARGA DE DATOS
# ==============================================================================

# --- 1a. Equity: ETF iShares MSCI COLCAP (ICOLCAP.CL) vía Yahoo Finance ------
#     (Si prefieres el índice COLCAP "puro", descárgalo manualmente de
#      SUAMECA/Investing.com como CSV y sáltate este bloque -> ver 1c)

fecha_inicio <- "2020-01-01"
fecha_fin    <- Sys.Date()

getSymbols("ICOLCAP.CL", src = "yahoo",
           from = fecha_inicio, to = fecha_fin, auto.assign = TRUE)

equity_px <- Cl(ICOLCAP.CL)                # precio de cierre
colnames(equity_px) <- "px_equity"

# --- 1b. Renta fija: TES ------------------------------------------------------
# El Banco de la República (SUAMECA) no tiene una API simple tipo quantmod,
# así que lo más práctico es:
#   1. Ir a https://suameca.banrep.gov.co/estadisticas-economicas/
#   2. Buscar "TES tasa fija" al plazo de referencia (p. ej. 10 años) o el
#      índice de precios de TES si está disponible
#   3. Descargar el CSV y ajustar la ruta/columnas abajo
#
# Se espera un CSV con columnas: fecha (YYYY-MM-DD), tasa (o precio) del TES

ruta_tes <- "https://raw.githubusercontent.com/andvarga-eco/econometrics_finance/refs/heads/main/notebooks/1_Tasa_cero_cupon_iqy.csv"   # <-- ajustar ruta al archivo descargado

tes_raw <- read_csv(ruta_tes, locale = locale(decimal_mark = ","))
# Ajusta nombres de columnas según el CSV real de SUAMECA:
names(tes_raw)[5]<-"fecha"
tes_raw <- tes_raw|>arrange(fecha)|>filter(Denominación=="Pesos colombianos" & Plazo=="Cinco años")
tes_raw<-tes_raw|>select(c(fecha,'Tasa (%)'))
tes_raw<-tes_raw|>filter(fecha>'2020-01-01')


# Conversión de la tasa cero cupón a precio

plazo_dias <- 1825   # plazo constante en días (1825 = 5 años); ajustar al nodo de curva usado
names(tes_raw)[2]<-"z_tasa"
tes_raw$z_tasa<-as.numeric(tes_raw$z_tasa)
tes_raw <- tes_raw |>
  mutate(
    precio_tes = 100 / (1 + z_tasa)^(plazo_dias / 365),
    ret_bond   = log(precio_tes / lag(precio_tes)) * 100  # rendimiento log, en % (misma escala que equity)
  )

# ==============================================================================
# 2. CONSTRUCCIÓN DEL PANEL DE RENDIMIENTOS
# ==============================================================================

equity_df <- data.frame(fecha = index(equity_px),
                         px_equity = as.numeric(equity_px)) %>%
  mutate(ret_equity = log(px_equity / lag(px_equity)) * 100) %>%
  select(fecha, ret_equity)

tes_df <- tes_raw %>%
  select(fecha, ret_bond)

# Unión por fecha común (alinea calendarios BVC vs. mercado de deuda)
panel <- inner_join(equity_df, tes_df, by = "fecha") %>%
  filter(!is.na(ret_equity), !is.na(ret_bond)) %>%
  arrange(fecha)

cat("Observaciones en el panel alineado:", nrow(panel), "\n")
head(panel)

# ==============================================================================
# 3. PRUEBAS DE ESTACIONARIEDAD (ADF)
# ==============================================================================

adf_equity <- ur.df(panel$ret_equity, type = "drift", selectlags = "AIC")
adf_bond   <- ur.df(panel$ret_bond,   type = "drift", selectlags = "AIC")

cat("\n--- ADF: ret_equity ---\n"); print(summary(adf_equity))
cat("\n--- ADF: ret_bond ---\n");   print(summary(adf_bond))
# Si el estadístico tau es más negativo que el valor crítico al 5%, se rechaza
# la raíz unitaria -> la serie es estacionaria y se puede usar el VAR en niveles.

# ==============================================================================
# 4. SELECCIÓN DE REZAGOS Y ESTIMACIÓN DEL VAR
# ==============================================================================

datos_var <- panel %>% select(ret_equity, ret_bond)

seleccion <- VARselect(datos_var, lag.max = 10, type = "const")
print(seleccion$selection)
p_optimo <- seleccion$selection["AIC(n)"]

modelo_var <- VAR(datos_var, p = p_optimo, type = "const")
summary(modelo_var)

serial.test(modelo_var,lags.bg = 5,type="BG")
serial.test(modelo_var,lags.pt=10,type="PT.asymptotic")

# ==============================================================================
# 5. CAUSALIDAD DE GRANGER
# ==============================================================================

causalidad_eq_a_bond <- causality(modelo_var, cause = "ret_equity")
causalidad_bond_a_eq <- causality(modelo_var, cause = "ret_bond")

cat("\n--- Granger: equity -> bond ---\n")
print(causalidad_eq_a_bond$Granger)
cat("\n--- Granger: bond -> equity ---\n")
print(causalidad_bond_a_eq$Granger)

# ==============================================================================
# 6. FUNCIONES IMPULSO-RESPUESTA (IRF) Y FEVD
# ==============================================================================

irf_resultado <- irf(modelo_var, impulse = "ret_equity", response = "ret_bond",
                      n.ahead = 10, boot = TRUE, ci = 0.95)
plot(irf_resultado, main = "Respuesta de bonos ante shock en equity")

fevd_resultado <- fevd(modelo_var, n.ahead = 10)
print(fevd_resultado)
plot(fevd_resultado)

# ==============================================================================
# 7. COBERTURA ESTÁTICA (mínima varianza, todo el periodo)
# ==============================================================================

Sigma_hat <- summary(modelo_var)$covres   # matriz de covarianza de residuos
cov_eq_bond <- Sigma_hat["ret_equity", "ret_bond"]
var_bond    <- Sigma_hat["ret_bond", "ret_bond"]
var_eq      <- Sigma_hat["ret_equity", "ret_equity"]

h_static <- -cov_eq_bond / var_bond
cat(sprintf("\nRatio de cobertura estático h* = %.4f\n", h_static))

var_hedged <- var_eq + h_static^2 * var_bond + 2 * h_static * cov_eq_bond
reduccion  <- 1 - var_hedged / var_eq
cat(sprintf("Reducción de varianza con cobertura estática: %.1f%%\n", reduccion * 100))

# ==============================================================================
# 8. COBERTURA DINÁMICA (VAR en ventana móvil)
# ==============================================================================

ventana <- 100   # días por ventana
paso    <- 20    # frecuencia de reestimación

n <- nrow(datos_var)
resultados_h <- data.frame()

for (inicio in seq(1, n - ventana, by = paso)) {
  fin <- inicio + ventana - 1
  sub_datos <- datos_var[inicio:fin, ]

  modelo_sub <- tryCatch(
    VAR(sub_datos, p = 1, type = "const"),
    error = function(e) NULL
  )
  if (is.null(modelo_sub)) next

  S <- summary(modelo_sub)$covres
  h_t <- -S["ret_equity", "ret_bond"] / S["ret_bond", "ret_bond"]

  resultados_h <- rbind(resultados_h,
                         data.frame(fecha_fin = panel$fecha[fin], h_t = h_t))
}

print(resultados_h)

# Gráfico del ratio de cobertura dinámico en el tiempo
ggplot(resultados_h, aes(x = fecha_fin, y = h_t)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = h_static, linetype = "dashed", color = "firebrick") +
  labs(title = "Ratio de cobertura dinámico (VAR, ventana móvil) vs. estático",
       subtitle = "Línea roja punteada = h* estático (todo el periodo)",
       x = "Fecha", y = expression(h[t])) +
  theme_minimal()

# ==============================================================================
# 9. COMPARACIÓN DE PnL: SIN COBERTURA vs. ESTÁTICA vs. DINÁMICA
# ==============================================================================

posicion_equity <- 1e6   # COP o USD, según la serie usada

pnl_sin <- c(); pnl_est <- c(); pnl_din <- c()

for (i in seq_len(nrow(resultados_h))) {
  idx_fin <- which(panel$fecha == resultados_h$fecha_fin[i])
  idx_ini <- idx_fin + 1
  idx_last <- min(idx_fin + paso, n)
  if (idx_ini > idx_last) next

  bloque <- datos_var[idx_ini:idx_last, ]
  h_t <- resultados_h$h_t[i]

  pnl_sin <- c(pnl_sin, posicion_equity * bloque$ret_equity / 100)
  pnl_est <- c(pnl_est, posicion_equity * bloque$ret_equity / 100 +
                 h_static * posicion_equity * bloque$ret_bond / 100)
  pnl_din <- c(pnl_din, posicion_equity * bloque$ret_equity / 100 +
                 h_t * posicion_equity * bloque$ret_bond / 100)
}

cat(sprintf("\nVolatilidad PnL SIN cobertura:     %s\n", format(sd(pnl_sin), big.mark = ",")))
cat(sprintf("Volatilidad PnL cobertura ESTÁTICA: %s (reducción %.1f%%)\n",
            format(sd(pnl_est), big.mark = ","), (1 - sd(pnl_est)/sd(pnl_sin)) * 100))
cat(sprintf("Volatilidad PnL cobertura DINÁMICA: %s (reducción %.1f%%)\n",
            format(sd(pnl_din), big.mark = ","), (1 - sd(pnl_din)/sd(pnl_sin)) * 100))

################################################################################
# NOTAS:
# - Ajusta 'ruta_tes' (y opcionalmente 'ruta_trm') a los archivos descargados
#   manualmente de SUAMECA (https://suameca.banrep.gov.co/).
# - Si usas el índice COLCAP puro en vez del ETF ICOLCAP.CL, reemplaza el
#   bloque 1a por la lectura de tu CSV descargado de Investing.com o SUAMECA,
#   con el mismo formato fecha/precio.
# - Para incluir la TRM como tercera variable, descomenta el bloque 1c,
#   agrégala a 'datos_var' y repite el pipeline (VARselect, VAR, IRF, FEVD).
# - La duración modificada usada para convertir tasa -> rendimiento de precio
#   de TES es una aproximación; ajústala al plazo/duración real del título.
################################################################################