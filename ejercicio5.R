# --- DEFINICIÓN DE FUNCIONES ---

#' Función 1: Realiza la simulación y devuelve los intervalos
#'
#' Esta función ejecuta el bucle de simulación y devuelve
#' una matriz con los límites inferiores y superiores.
#'
simular_intervalos <- function(n, confianza, metodo, lambda_true, num_sims) {
  
  # Matriz para guardar los límites
  limites <- matrix(NA, nrow = num_sims, ncol = 2)
  
  for (i in 1:num_sims) {
    muestra <- rpois(n, lambda = lambda_true)
    media_muestral <- mean(muestra)
    
    if (metodo == "asintotico") {
      # Método de Wald (asintótico)
      z <- qnorm(1 - (1 - confianza) / 2)
      
      if (media_muestral == 0) {
        # El método de Wald falla en 0. Usamos poisson.test(0) 
        # como un reemplazo más robusto que un IC de [0,0].
        ic <- poisson.test(0, T = n, conf.level = confianza)$conf.int
        limites[i, 1] <- ic[1]
        limites[i, 2] <- ic[2]
      } else {
        error_est <- z * sqrt(media_muestral / n)
        limites[i, 1] <- media_muestral - error_est
        limites[i, 2] <- media_muestral + error_est
      }
      
    } else {
      # Método Exacto (usado como "General" para n pequeño)
      ic <- poisson.test(sum(muestra), T = n, conf.level = confianza)$conf.int
      limites[i, 1] <- ic[1]
      limites[i, 2] <- ic[2]
    }
  }
  
  return(limites) # Devuelve la matriz de resultados
}


#' Función 2: Grafica los intervalos de confianza
#'
#' Toma la matriz de límites, el valor verdadero y un título,
#' y genera una gráfica.
#'
graficar_intervalos <- function(limites, lambda_true, titulo) {
  
  num_sims <- nrow(limites)
  cobertura <- 0
  colores <- rep(NA, num_sims)
  
  # Calcular cobertura y colores para la gráfica
  for (i in 1:num_sims) {
    # Revisar NAs (aunque no deberían ocurrir con este código)
    if (!is.na(limites[i, 1]) && !is.na(limites[i, 2])) {
      if (limites[i, 1] <= lambda_true && limites[i, 2] >= lambda_true) {
        cobertura <- cobertura + 1
        colores[i] <- "black" # Cubre
      } else {
        colores[i] <- "red"   # No cubre
      }
    } else {
      colores[i] <- "gray" # En caso de NAs
    }
  }
  
  # Calcular porcentaje de cobertura
  cobertura_pct <- round((cobertura / num_sims) * 100)
  
  # Graficar
  plot(1:num_sims, ylim = range(c(limites, lambda_true), na.rm = TRUE), type = "n",
       xlab = "Muestra", ylab = expression(lambda), main = titulo)
  segments(x0 = 1:num_sims, y0 = limites[, 1],
           x1 = 1:num_sims, y1 = limites[, 2], col = colores)
  abline(h = lambda_true, col = "blue", lwd = 2) # Valor verdadero
  legend("topright", legend = paste("Cobertura:", cobertura_pct, "%"), bty = "n")
}


#' Función 3: Analiza y reporta las longitudes de los intervalos
#'
#' Toma una matriz de límites y un título, e imprime en
#' la consola el análisis de longitudes.
#'
analizar_longitudes <- function(limites, titulo) {
  
  # Calcular las longitudes de todos los intervalos
  # na.rm = TRUE por si acaso
  longitudes <- limites[, 2] - limites[, 1]
  
  # Encontrar el índice (número de muestra) del más corto y más largo
  idx_min <- which.min(longitudes)
  idx_max <- which.max(longitudes)
  
  # Obtener los intervalos y longitudes correspondientes
  intervalo_min <- limites[idx_min, ]
  longitud_min <- longitudes[idx_min]
  
  intervalo_max <- limites[idx_max, ]
  longitud_max <- longitudes[idx_max]
  
  # Imprimir los resultados en la consola
  cat("\n Análisis de Longitud para:", titulo, "---\n")
  cat("Intervalo MÁS CORTO (Muestra #", idx_min, "):\n")
  cat("  [", round(intervalo_min[1], 4), ", ", round(intervalo_min[2], 4), "]\n")
  cat("  Longitud:", round(longitud_min, 4), "\n\n")
  
  cat("Intervalo MÁS LARGO (Muestra #", idx_max, "):\n")
  cat("  [", round(intervalo_max[1], 4), ", ", round(intervalo_max[2], 4), "]\n")
  cat("  Longitud:", round(longitud_max, 4), "\n")
}


# Configuración inicial
set.seed(123) # Para reproducibilidad
lambda_true <- 8
num_sims <- 100


# Ejecutar las simulaciones y gráficas

# --- (a) n=15, Conf=80%, Método General ---
titulo_a <- "(a) n=15, 80% Conf (General)"
limites_a <- simular_intervalos(n = 15, confianza = 0.80, metodo = "exacto",
                                lambda_true = lambda_true, num_sims = num_sims)
graficar_intervalos(limites = limites_a, lambda_true = lambda_true, titulo = titulo_a)
analizar_longitudes(limites = limites_a, titulo = titulo_a)
# --- (b) n=15, Conf=95%, Método General ---
titulo_b <- "(b) n=15, 95% Conf (General)"
limites_b <- simular_intervalos(n = 15, confianza = 0.95, metodo = "exacto",
                                lambda_true = lambda_true, num_sims = num_sims)
graficar_intervalos(limites = limites_b, lambda_true = lambda_true, titulo = titulo_b)
analizar_longitudes(limites = limites_b, titulo = titulo_b)
# --- (c) n=150, Conf=80%, Método Asintótico ---
titulo_c <- "(c) n=150, 80% Conf (Asintótico)"
limites_c <- simular_intervalos(n = 150, confianza = 0.80, metodo = "asintotico",
                                lambda_true = lambda_true, num_sims = num_sims)
graficar_intervalos(limites = limites_c, lambda_true = lambda_true, titulo = titulo_c)
analizar_longitudes(limites = limites_c, titulo = titulo_c)
# --- (d) n=150, Conf=95%, Método Asintótico ---
titulo_d <- "(d) n=150, 95% Conf (Asintótico)"
limites_d <- simular_intervalos(n = 150, confianza = 0.95, metodo = "asintotico",
                                lambda_true = lambda_true, num_sims = num_sims)
graficar_intervalos(limites = limites_d, lambda_true = lambda_true, titulo = titulo_d)
analizar_longitudes(limites = limites_d, titulo = titulo_d)

