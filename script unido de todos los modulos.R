# =============================================================================
# SCRIPT MAESTRO: Pipeline de Expresión Diferencial RNA-seq (GSE160299)
# Autor: Sandra
# Universidad: UAM
# Descripción: Orquestador para ejecutar todos los módulos individuales del proyecto
# =============================================================================

# Limpiar el entorno de trabajo antes de empezar
rm(list = ls())

# Mensaje de inicio
cat("====================================================\n")
cat("Iniciando el Pipeline Automatizado de RNA-seq...\n")
cat("====================================================\n\n")

# Tiempo inicial
inicio_pipeline <- Sys.time()

# 1. Lista ordenada de los 12 scripts individuales a ejecutar
scripts_pipeline <- c(
  "0_configuración.R",
  "01_paquetes.R",
  "02_Plot helper.R",
  "03_GEO download.R",
  "04_ Metadata parsing.R",
  "05_Count matrix parsing.R",
  "06_ Metadata and library QC.R",
  "07_DESeq2 model.R",
  "08_GEO2R-like QC plots.R",
  "09_ Differential expression.R",
  "10_Heatmaps and gene profiles.R",
  "11_Run manifest and session info.R"
)

# 2. Bucle para ejecutar cada script con manejo de errores básico
for (script in scripts_pipeline) {
  if (file.exists(script)) {
    cat(paste0("\n[", format(Sys.time(), "%H:%M:%S"), "] Ejecutando: ", script, "...\n"))
    
    # Ejecuta el script individual en el entorno global
    source(script, local = FALSE, echo = FALSE)
    
    cat(paste0("[OK] Finalizó con éxito: ", script, "\n"))
  } else {
    stop(paste0("ERROR CRÍTICO: No se encontró el archivo '", script, 
                "'. Verifica que esté en esta misma carpeta y que el nombre coincida exactamente."))
  }
}

# 3. Mensaje final de éxito
final_pipeline <- Sys.time()
tiempo_total <- round(difftime(final_pipeline, inicio_pipeline, units = "mins"), 2)

cat("\n====================================================\n")
cat("¡PIPELINE FINALIZADO CON ÉXITO!\n")
cat("Los 12 módulos se ejecutaron de principio a fin.\n")
cat("Las gráficas y tablas están listas en la carpeta de resultados.\n")
cat("Tiempo total de ejecución: ", tiempo_total, " minutos.\n")
cat("====================================================\n")