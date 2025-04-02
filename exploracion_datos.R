

# librerias ---------------------------------------------------------------
library(readxl)
library(tidyverse)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("SummarizedExperiment")
# BiocManager::install("SparseArray")
library("SummarizedExperiment")


# estilo graficos ---------------------------------------------------------
my_theme <- 
  theme_bw() +
  theme(legend.position = "bottom")

theme_set(my_theme)

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
  dplyr::count(SampleType, Class)

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
metabo_sumexp <- SummarizedExperiment(
  # assays es una lista en la que cada elemento es un experimento (conjunto de datos) distintos
  assays = list(raw_data = data.metabolomics),
  colData = metadatos_muestras,
  rowData = metadatos_metabolitos)

# exploracion datos -------------------------------------------------------

## heatmap - overview ----
# sample 111 y 105 - outliers
# QC las unicas que cluster
# algunos metabolitos faltan en muchas muestras
library(pheatmap)
colors_heat <- list(Class = c("BN" = "blue", 
                              "GC" = "darkred", 
                              "HE" = "lightblue", 
                              "QC" = "grey"),
                    SampleType = c("Sample" = "lightblue",
                                   "QC" = "grey"
                                   )
                    )

pheatmap(
  log10(assay(metabo_sumexp)),                
  scale = "row",  # escalar los datos por fila (metabolito)  
  cluster_rows = FALSE,
  # clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "euclidean",  
  clustering_method = "complete",
  show_rownames = TRUE,         
  show_colnames = TRUE,
  annotation_col = as.data.frame(colData(metabo_sumexp)[1]), # Class para anotar muestras/columnas
  annotation_colors = colors_heat,
  # annotation_row = as.data.frame(rowData(metabo_sumexp)[1]), # nombre_metabolito para anotar metabolitos/filas
  main = "Heatmap de los metabolitos por muestras", 
  fontsize = 7,
  annotation_legend = TRUE)


## QC tutorial ----
# included metabolites
metabolites_inclusion <- rowData(metabo_sumexp) %>% 
  as.data.frame() %>% 
  # missing in less than 10% of all samples AND variation of detection <20%
  filter(Perc_missing < 10 & QC_RSD < 20) %>%  
  pull(nombre_metabolito)

# represent exclusion
rowData(metabo_sumexp) %>% 
  as.data.frame() %>% 
  # add exclusion
  mutate(inclusion = ifelse(Perc_missing < 10 & QC_RSD < 20, T, F)) %>% 
  # filter(Perc_missing != 0) %>% 
  # arrange(desc(Perc_missing))
  ggplot(aes(x = fct_reorder(nombre_metabolito, -Perc_missing), # orden decreciente % missing
  )) +
  geom_col(aes(y = QC_RSD/2,
               color = inclusion,
  ), position = "dodge", # escalamos a la mitad del porcentaje
  fill = "darkred",
  alpha = 0.5) +
  geom_col(aes(y = Perc_missing,
               color = inclusion), 
           position = "dodge", 
           fill = "grey",
           alpha = 0.7) +
  scale_color_manual(values = c("white", "black")) +
  scale_y_continuous(name = "Perc. missing", 
                     limits = c(0,70),
                     breaks = seq(0,70, 10),
                     sec.axis = sec_axis(~ .*2, name = "QC-RSD",
                                         breaks = seq(0,160,20)
                     )) +  # duplicamos la escala para ver cada eje bien
  geom_hline(yintercept = c(10), linetype = "dotted") +
  labs(title = "Metabolites % of missing detection (left-y) and variation of measurement (QC-RSC; right-y) ",
       subtitle = paste0("Colored bar represents included metabolites due to %missing<10% and QC-RSD<20 (n=", length(metabolites_inclusion), ")")
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line.y.right = element_line(color = "darkred"))

## missing readouts in QC ----
# dataset_exploration <- 
dataset_exploration <- assay(metabo_sumexp) %>% 
  as.data.frame() %>% 
  rownames_to_column("var") %>% 
  pivot_longer(cols = -c(1)) %>% 
  left_join(., as.data.frame(rowData(metabo_sumexp)[1]) %>% 
              rownames_to_column("var")) %>% 
  left_join(., as.data.frame(colData(metabo_sumexp)[2]) %>% 
              rownames_to_column(var = "name")
              )

# in all samples - might be relevant if some are i.e., tumor-specific
missing_overview <- dataset_exploration %>% 
  group_by(nombre_metabolito, Class) %>% 
  summarise(n_missing = sum(is.na(value)),
            percent_missing = mean(is.na(value))*100) %>% 
  ungroup()

head(missing_overview)

# in QC - could not be measured properly
missing_overview %>% 
  filter(Class == "QC" & n_missing != 0) %>% 
  ggplot(aes(x = fct_reorder(nombre_metabolito, -percent_missing), # orden decreciente % missing
             y = percent_missing)) +
  geom_col(aes(color = ifelse(nombre_metabolito %in% metabolites_inclusion, T, F)),
           fill = "grey80") +
  geom_text(aes(label = n_missing, 
                y = percent_missing/2),
            color = "white"
  ) +
  scale_color_manual(values = c("white", "black")) +
  labs(title = "%missing in QC samples for (missing) metabolites",
       subtitle = "Number on the bar indicates the number of QC samples with the missing readout\nBorde negro indica los metabolitos incluidos en estudio segun criterios del tutorial",
       y = "%QC samples missing",
       x = element_blank()) +
  scale_y_continuous(limits = c(0,60),
                     breaks = seq(0,60, 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## boxplots ----
dataset_exploration %>%
  ggplot(aes(x = nombre_metabolito,
             y = value)) +
  geom_line(aes(group = name),
            color = "grey",
            alpha = .5) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5))

## POMA ----
# BiocManager::install("POMA")
library(POMA)
library(ggtext)
library(magrittr)

### dataset analysis ----
# subset with QC conditions - 52 analytes
(metabo_sumexp_qc <- metabo_sumexp[rowData(metabo_sumexp)$Perc_missing < 10 & rowData(metabo_sumexp)$QC_RSD < 20])

### missing data imputation ----
assays(metabo_sumexp_qc)[["imputed_data"]] <- assay(metabo_sumexp_qc %>% 
                                                      PomaImpute(method = "knn", 
                                                                 zeros_as_na = TRUE, 
                                                                 remove_na = TRUE, 
                                                                 cutoff = 10)) # ya hemos preseleccionados metabolitos con <10% missing, por eso outcome de "0 features removed"

# ver objeto con raw e imputed data
metabo_sumexp_qc

# convertir a factor - necesario Poma
colData(metabo_sumexp_qc)$SampleType <- as.factor(colData(metabo_sumexp_qc)$SampleType)

# antes NAs
table(is.na(assay(metabo_sumexp_qc, "raw_data")))
# ahora no
table(is.na(assay(metabo_sumexp_qc, "imputed_data")))

### normalization ----
# log (la que use para primera inspeccion)
assays(metabo_sumexp_qc)[["normalized_log"]] <- assay(SummarizedExperiment(
  assays = list(imputed_data = assays(metabo_sumexp_qc)[["imputed_data"]]), 
  colData = colData(metabo_sumexp_qc), 
  rowData = rowData(metabo_sumexp_qc)) %>% 
    PomaNorm(method = "log_scaling"))


# log-pareto: valor/sqrt(sd(metabolito)) - conserva mejor la diferencia entre metabolitos con mucha vs poca varianza
assays(metabo_sumexp_qc)[["normalized_logpareto"]] <- assay(SummarizedExperiment(
  assays = list(imputed_data = assays(metabo_sumexp_qc)[["imputed_data"]]), 
  colData = colData(metabo_sumexp_qc), 
  rowData = rowData(metabo_sumexp_qc)) %>% 
    PomaNorm(method = "log_pareto"))

# escala original
# box samples por metabolito
SummarizedExperiment(assays = list(imputed = assay(metabo_sumexp_qc, "imputed_data")),
                     colData = colData(metabo_sumexp_qc), 
                     rowData = rowData(metabo_sumexp_qc)) %>% 
  PomaBoxplots(.,
               x = "samples",
               theme_params = list(axistext = "y")
               )

# distribucion metabolitos
SummarizedExperiment(assays = list(imputed = assay(metabo_sumexp_qc, "imputed_data")),
                     colData = colData(metabo_sumexp_qc), 
                     rowData = rowData(metabo_sumexp_qc)) %>% 
  PomaDensity(.,
              x = "features", theme_params = list(legend_title = FALSE)) +
  theme(legend.position = "none")

# log-pareto
SummarizedExperiment(assays = list(imputed = assay(metabo_sumexp_qc, "normalized_logpareto")),
                     colData = colData(metabo_sumexp_qc), 
                     rowData = rowData(metabo_sumexp_qc)) %>% 
  PomaBoxplots(.,
               x = "samples",
               theme_params = list(axistext = "y")
  )

# distribucion
SummarizedExperiment(assays = list(imputed = assay(metabo_sumexp_qc, "normalized_logpareto")),
                     colData = colData(metabo_sumexp_qc), 
                     rowData = rowData(metabo_sumexp_qc)) %>% 
  PomaDensity(.,
              x = "features", theme_params = list(legend_title = FALSE)) +
  theme(legend.position = "none")

# log - mejor a simple vista por menos outliers
SummarizedExperiment(assays = list(imputed = assay(metabo_sumexp_qc, "normalized_log")),
                     colData = colData(metabo_sumexp_qc), 
                     rowData = rowData(metabo_sumexp_qc)) %>% 
  PomaBoxplots(.,
               x = "samples",
               theme_params = list(axistext = "y")
  )

# distribucion
SummarizedExperiment(assays = list(imputed = assay(metabo_sumexp_qc, "normalized_log")),
                     colData = colData(metabo_sumexp_qc), 
                     rowData = rowData(metabo_sumexp_qc)) %>% 
  PomaDensity(.,
              x = "features", 
              theme_params = list(legend_title = FALSE)) +
  theme(legend.position = "none")

### outliers ----
metabo_sumexp_norm_log <- SummarizedExperiment(assays = list(imputed = assay(metabo_sumexp_qc, "normalized_log")),
                                               colData = colData(metabo_sumexp_qc), 
                                               rowData = rowData(metabo_sumexp_qc))

PomaOutliers(metabo_sumexp_norm_log)$polygon_plot

PomaOutliers(metabo_sumexp_norm_log)$outliers

# que muestra es
colData(metabo_sumexp_norm_log)["sample_42", ]

### preprocessed analysis dataset---- 
# 52 metabolitos en 139 muestras
pre_processed <- PomaOutliers(metabo_sumexp_norm_log)$data

# add name of metabolites
rownames(pre_processed) <- rowData(metabo_sumexp_qc)[[1]]

pca_qc_samples <- PomaPCA(pre_processed,
                          scale = FALSE,
                          labels = FALSE,
                          ellipse = FALSE)

pca_qc_samples$biplot

# analisis datos -----
# remove QC
pre_processed_analysis <- pre_processed[,colData(pre_processed)$SampleType == "Sample"]

colData(pre_processed_analysis)$SampleType <- droplevels(colData(pre_processed_analysis)$SampleType)
colData(pre_processed_analysis)$Class <- droplevels(colData(pre_processed_analysis)$Class)

table(colData(pre_processed_analysis)$Class)

## PCA-Samples----
pca_samples <- PomaPCA(pre_processed_analysis,
                       outcome = "Class",
                       scale = FALSE,
                       labels = FALSE,
                       ellipse = FALSE)

pca_samples$factors_plot

# metabolitos posicionando los puntos
pca_samples$biplot

# pc2 drives GC downwards - check most contributing metabolites
pca_samples$loadings %>% 
  filter(PC2<0) %>% 
  arrange(PC2) %>% 
  print(n = Inf)

metabolites_gcpc2 <- pca_samples$loadings %>% 
  filter(PC2<0) %>% 
  arrange(PC2) %>% pull(feature)

## heatmap ----
pre_processed_analysis_he_gc <- pre_processed[,colData(pre_processed)$Class %in% c("GC", "HE")]

colData(pre_processed_analysis_he_gc)$Class <- droplevels(colData(pre_processed_analysis_he_gc)$Class)

pre_processed_analysis_he_gc[rownames(pre_processed_analysis_he_gc) %in% metabolites_gcpc2] %>% 
  PomaHeatmap(covs = c("Class"), # covariates to plot (e.g., treatment, sex, etc)
              feature_names = TRUE
  )

## limma ----
limma_results <- pre_processed_analysis_he_gc %>% 
  PomaLimma(contrast = "HE-GC", 
            outcome = "Class",
            covs = NULL,
            adjust = "fdr",
            block = NULL)

# volcano
limma_results %>% 
  dplyr::select(feature, log2FC, pvalue) %>% 
  PomaVolcano(labels = TRUE)

# sacar los expresados diferencial y significativamente
metabolitos_limma <- limma_results %>% 
  filter(pvalue<=0.05) %>% 
  pull(feature)

# heatmap
pre_processed_analysis_he_gc[rownames(pre_processed_analysis_he_gc) %in% metabolitos_limma] %>% 
  PomaHeatmap(covs = c("Class"), # covariates to plot (e.g., treatment, sex, etc)
              feature_names = TRUE
              )

# overlap PCA y limma
metabolitos_limma[!(metabolitos_limma %in% metabolites_gcpc2)]

metabolites_gcpc2[!(metabolites_gcpc2 %in% metabolitos_limma)]

# outputs -----------------------------------------------------------------

## rda summarizedexperiment----
save(metabo_sumexp, file = "summarized_experiment.rda")

## datos txt ----
write.table(x = assay(metabo_sumexp),
            file = "datos.txt",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
