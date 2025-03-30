

# librerias ---------------------------------------------------------------
library(readxl)
library(tidyverse)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("SparseArray")
library("SummarizedExperiment")

# descarga datos ----------------------------------------------------------

## importamos ----
# comprobamos si el documento excelt tiene una o varias hojas
excel_sheets("GastricCancer_NMR.xlsx")

# importamos los datos
datos_raw <- read_xlsx(path = "GastricCancer_NMR.xlsx", sheet = "Data")

# importamos la hoja ¨peak¨
metadata <- read_xlsx(path = "GastricCancer_NMR.xlsx", sheet = "Peak")

## descripcion de datos----
### datos ----
str(datos_raw)
colnames(datos_raw)
# los datos tienen 140 filas (samples) y 153 columnas incluyendo 4 columnas con metadatos de los 
# samples y 149 metabolitos

# tipos de samples: 17 QC y 123 muestras experimentales
table(datos_raw$SampleType)

# clasificados en cuatro clases
datos_raw %>% 
  count(SampleType, Class)

### metadatos metabolitos----
str(metadata)
# 149 filas para cada metabolito (features) y 5 columnas con informacion de cada metabolito

# objeto SummarizedExperiment ---------------------------------------------
## assaydata ----
# matriz con medidas de las observaciones
data.metabolomics <- datos_raw[, c(2,5:153)] %>% 
  pivot_longer(cols = -c(1)) %>% 
  pivot_wider(names_from = SampleID) %>% 
  column_to_rownames(var = "name") %>% 
  as.matrix()

## colData ----
# dataframe con informacion (en las columnas) de las muestras (en las filas)
metadatos_muestras <- datos_raw[, 2:4] %>% 
  column_to_rownames(var = "SampleID")
# head(metadatos_observaciones)

# revisamos que las filas de colData se llaman igual que las columnas de data.assay
table(colnames(data.metabolomics) == rownames(metadatos_muestras))

## rowData ----
metadatos_metabolitos <- metadata[, 2:5] %>% 
  # cambiamos nombres
  dplyr::rename(nombre_metabolito = Label) %>% 
  column_to_rownames(var = "Name")
  
## summarizedExperiment ----
metabo_sumexp <- ßSummarizedExperiment(# assays es una lista en la que cada elemento es un experimento (conjunto de datos) distintos
                     assays = list(metabolomics = data.metabolomics),
                     colData = metadatos_muestras,
                     rowData = metadatos_metabolitos)

# exploracion datos -------------------------------------------------------



# outputs -----------------------------------------------------------------

## rda summarizedexperiment----
save(metabo_sumexp, file = "summarized_experiment.rda")

# input for markdown metadatos

## datos txt ----

