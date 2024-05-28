library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # solo para mm10
library(BSgenome.Mmusculus.UCSC.mm10)
#library(EnsDb.Mmusculus.v75) # solo para mm9
#library(EnsDb.Hsapiens.v86) # hg38?
library(ggplot2)
library(patchwork)

datos <- '8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5'
metadatos <- '8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_singlecell.csv'
fragmentos <- '8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz'
genoma <- 'mm9'
#genoma <- 'hg38'

genes_interesantes <- c('Ndufa4l2', 'Cox4i2', 'Vtn', 'Notch3', 'Pdgfrb', 'Higd1b')

datos_scrna_seq <- '../HoltzBrainOutput/int/exprMat_filtered_log.Rds'

#tfs_de_scenic <- unique(readRDS('../HoltzBrainOutput/tfs_descubiertos.Rds')$regulones.TF)
tfs_de_scenic <- c ('Tbx2', 'Zic1', 'Hic1', 'Foxd1', 'Mef2c')

tipos_celulares <- read.csv ('../HoltzBrainOutput/int/cels', header = FALSE)

if (!file.exists('counts.Rds')) {
    counts <- Read10X_h5(datos)
    saveRDS (counts, 'counts.Rds')
} else {
    counts <- readRDS ('counts.Rds')
}
metadata <- read.csv (metadatos, header = TRUE, row.names = 1)
brain_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(':','-'),
    genome = genoma,
    fragments = fragmentos,
    min.cells = 1
)

brain <- CreateSeuratObject(
    counts = brain_assay,
    assay = 'peaks',
    project = 'ATAC',
    meta.data = metadata
)

if (!file.exists('annotations.Rds')) {
    annotations <- GetGRangesFromEnsDb (ensdb = EnsDb.Mmusculus.v79)
    saveRDS (annotations, 'annotations.Rds')
} else {
    annotations <- readRDS ('annotations.Rds')
}

#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- genoma

Annotation(brain) <- annotations

########################
## CONTROL DE CALIDAD ##
########################

# Calcular la intensidad de señal de los nucleosomas.
brain <- NucleosomeSignal (object = brain)

brain$nucleosome_group <- ifelse(brain$nucleosome_signal > 4, 'NS>4', 'NS<4')

# NS<4: células bien secuenciadas
# NS>4: células mal secuenciadas
# En NS<4 deberían verse tres picos: uno corresponde a zonas sin nucleosomas (eucromatina que mapea a TSSs), el segundo a lecturas de nucleosoma/mononucleosoma, que flanquea regiones abiertas, y el tercer pico a lecturas dinucleosomicas. Los tres picos indican calidad
FragmentHistogram(object=brain, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
## Ndufa4l2 Mmus
#FragmentHistogram(object=brain, group.by = 'nucleosome_group', region = 'chr10-127345000-127356000')
## Cox4i2 Mmus
#FragmentHistogram(object=brain, group.by = 'nucleosome_group', region = 'chr2-152590000-152608000')
## Higd1b Mmus
#FragmentHistogram(object=brain, group.by = 'nucleosome_group', region = 'chr11-102716675-102738866')

# Otra métrica de calidad: comprobar enriquecimiento de la integración de la Tn5 en los TSSs: nº de sitios de integración de la Tn5 alrededor del TSS normalizado por el nº de sitios de integración en sus flancos 2kb alrededor. 10-15 aceptable, >15 ideal
if (!file.exists('brain_with_tss.Rds')) {
    brain <- TSSEnrichment(brain, fast = FALSE)
    saveRDS (brain, 'brain_with_tss.Rds')
} else {
    brain <- readRDS ('brain_with_tss.Rds')
}

brain$high.tss <- ifelse (brain$TSS.enrichment > 2, 'High', 'Low')
TSSPlot (brain, group.by = 'high.tss') #+NoLegend()

brain$pct_reads_in_peaks <- brain$peak_region_fragments / brain$passed_filters * 100
brain$blacklist_ratio <- brain$blacklist_region_fragments / brain$peak_region_fragments

# Último paso de control de calidad: ratio de lecturas que mapean en picos, ratio de lecturas que caen en regiones prohibidas
VlnPlot(
    object = brain,
    features = c(
        'pct_reads_in_peaks',
        'peak_region_fragments',
        'TSS.enrichment'
    ),
    pt.size = 0.1,
    ncol = 3
)
VlnPlot(
    object = brain,
    features = c(
        'blacklist_ratio',
        'nucleosome_signal'
    ),
    pt.size = 0.1,
    ncol = 2
)

# Limpiamos los datos que no sean de calidad
brain <- subset(
    x = brain,
    # Células que tengan entre 3000-100000 lecturas que mapeen en picos
    subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    # Ratio de picos de >40%
    pct_reads_in_peaks > 40 &
    # Porcentaje de lecturas en regiones ruidosas <0.025%
    blacklist_ratio < 0.025 &
    # Señal nucleosómica <4
    nucleosome_signal < 4 &
    # Puntuación TSS >2
    TSS.enrichment > 2
)

# Normalizado y reducción lineal de dimensionalidad.
# TF-IDF normaliza células y picos, dando importancia a picos poco comunes
brain <- RunTFIDF(brain)

#brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- FindTopFeatures(brain, min.cutoff = 'q75') # 25% de los mejores picos

# SVD: hermana de tSNE y UMAP. Reduce la dimensionalidad en la matriz TF-IDF usando solo los picos del paso previo. Combinar TF-IDF con SVD recibe el nombre de Latent Semantic Indexing (LSI)
if (!file.exists('data_svd.Rds')) {
    brain <- RunSVD(object = brain)
    saveRDS (brain, 'data_svd.Rds')
} else {
    brain <- readRDS ('data_svd.Rds')
}

# La primera componente de LSI captura la profundidad del secuenciado (variación técnica) y no varianza biológica. Por ello, debe eliminarse del análisis posterior. Se encuentra inversamente muy correlacionada con el nº de lecturas total en las células, así que habremos de omitirla
DepthCor(brain)

DimPlot(
    object = brain,
    label = TRUE,
    reduction = 'lsi',
    dims = c(2,3)
)

brain <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindClusters(
  object = brain,
  algorithm = 3,
  resolution = 1.2,
  verbose = TRUE
)
brain <- RunUMAP(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)

# Clustering UMAP. 
DimPlot(object = brain, label = TRUE) #+ NoLegend()

# Generar una matriz de actividad génica. Extraemos coordenadas de los genes extendidas 2 kb upstream (para incluir al promotor) y se cuentan en cada célula el nº de lecturas que mapean en esa región. Esto pretende emular como si hubiéramos secuenciado el transcriptoma celular permitiéndonos tratar los datos de ATAC-seq como si fuera RNA-seq.
if (!file.exists('gene_activities.Rds')) {
    gene.activities <- GeneActivity(brain)
    saveRDS (gene.activities, 'gene_activities.Rds')
} else {
    gene.activities <- readRDS ('gene_activities.Rds')
}

# add the gene activity matrix to the Seurat object as a new assay
brain[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = brain,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)

# show cell types with at least 50 cells
idents.plot <- names(which(table(Idents(brain)) > 50))

brain <- RegionStats(
    brain,
    genome = BSgenome.Mmusculus.UCSC.mm10
)

brain <- LinkPeaks(
    object = brain,
    peak.assay = 'ATAC',
    expression.assay = 'SCT',
    genes.use = c ('Ndufa4l2', 'Cox4i2', 'Tbx2', 'Zic1', 'Hic1', 'Foxd1')
)

# Nº integraciones de la Tn5 
CoveragePlot(
  object = brain,
  region = c("Ndufa4l2"),
  expression.assay = 'SCT',
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = c("Higd1b"),
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = c("Cox4i2"),
  expression.assay = 'SCT',
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = c("Pdgfrb"),
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = 'Vtn',
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = 'Notch3',
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = tfs_de_scenic[1],
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = tfs_de_scenic[2],
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = tfs_de_scenic[3],
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = tfs_de_scenic[4],
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = tfs_de_scenic[5],
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)

# CAMBIO
#DefaultAssay(brain) <- 'RNA'

FeaturePlot(
  object = brain,
  features = 'Ndufa4l2',
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = "Vtn",
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = "Higd1b",
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = c('Pdgfrb'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = c('Vtn'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = c('Cox4i2'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = tfs_de_scenic[1],
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = tfs_de_scenic[2],
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = tfs_de_scenic[3],
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = tfs_de_scenic[4],
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
FeaturePlot(
  object = brain,
  features = tfs_de_scenic[5],
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)
## INTEGRAR SCRNA-SEQ ##

# Load the pre-processed scRNA-seq data
sc_rna_seq <- readRDS(datos_scrna_seq)
#sc_rna_seq <- SCopeLoomR::get_dgem(SCopeLoomR::open_loom(datos_scrna_seq))
#sc_rna_seq <- UpdateSeuratObject(sc_rna_seq)

if (!file.exists('scRNA_seq_data.Rds')) {
    sc_rna_seq <- CreateSeuratObject(
        #t(sc_rna_seq),
        sc_rna_seq,
        meta.data = tipos_celulares
    )
    sc_rna_seq <- FindVariableFeatures(
      object = sc_rna_seq,
      nfeatures = 5000
    )
    sc_rna_seq <- NormalizeData(
      object = sc_rna_seq,
      assay = 'RNA',
      normalization.method = 'LogNormalize',
      scale.factor = median(sc_rna_seq$nCount_RNA)
    )
    saveRDS (sc_rna_seq, 'scRNA_seq_data.Rds')
} else {
    sc_rna_seq <- readRDS ('scRNA_seq_data.Rds')
}

if (!file.exists('transfer_anchors.Rds')) {
    transfer.anchors <- FindTransferAnchors(
        reference = sc_rna_seq,
        query = brain,
        features = VariableFeatures(object = sc_rna_seq),
        reference.assay = "RNA",
        query.assay = "RNA",
        reduction = "cca",
        dims = 1:30
    )
    saveRDS (transfer.anchors, 'transfer_anchors.Rds')
} else {
    transfer.anchors <- readRDS ('transfer_anchors.Rds')
}

#DefaultAssay(sc_rna_seq) <- 'RNA'
#DefaultAssay(brain) <- 'RNA' # Cambiamos a RNA para hacer ambos objetos compatibles. Antes era 'peaks'
#transfer.anchors <- FindTransferAnchors(
#  reference = sc_rna_seq,
#  features = genes_rna_seq,
#  query = brain,
#  reduction = 'cca',
#  dims = 1:30
#)

if (!file.exists('predicted_labels.Rds')) {
    predicted.labels <- TransferData(
        anchorset = transfer.anchors,
        refdata = sc_rna_seq$orig.ident,
        weight.reduction = brain[['lsi']],
        dims = 2:30
    )
    saveRDS (predicted.labels, 'predicted_labels.Rds')
} else {
    predicted.labels <- readRDS ('predicted_labels.Rds')
}

brain <- AddMetaData(object = brain, metadata = predicted.labels)

#sc_rna_seq <- RunUMAP(
#    object = sc_rna_seq,
#    reduction = 'lsi',
#    dims = 2:30
#)

# AQUI
#plot1 <- DimPlot(sc_rna_seq, group.by = 'subclass', label = TRUE, repel = TRUE) + ggtitle('scRNA-seq')
#plot1

plot2 <- DimPlot(brain, group.by = 'predicted.id', label = TRUE, repel = TRUE) + ggtitle('scATAC-seq predicted labels')
plot2

RidgePlot (brain, features = c('Ndufa4l2', 'Cox4i2', 'Pdgfrb', 'Vtn', 'Notch3'))
RidgePlot (brain, features = c('Tbx2', 'Zic1', 'Hic1', 'Foxd1', 'Mef2c'))

CoveragePlot(
  object = brain,
  region = c("Ndufa4l2"),
  idents = unique(predicted.labels),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = c("Cox4i2"),
  idents = predicted.labels,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = c("Pdgfrb"),
  idents = predicted.labels,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = c("Vtn"),
  idents = predicted.labels,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)
CoveragePlot(
  object = brain,
  region = c("Notch3"),
  idents = predicted.labels,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)

#switch back to working with peaks instead of gene activities
DefaultAssay(brain) <- 'peaks'

# ENCONTRAR PICOS DE MARCADORES SEGÚN CLUSTERING PROPIO
if (!file.exists('da_peaks_marcadores_pericitos.Rds')) {
    da_peaks <- FindMarkers(
        object = brain,
        ident.1 = 20, # El cluster de pericitos (más probable)
        test.use = 'LR',
        latent.vars = 'nCount_peaks'
    )
    saveRDS (da_peaks, 'da_peaks_marcadores_pericitos.Rds')
} else {
    da_peaks <- readRDS ('da_peak_marcadores_pericitos.Rds')
}

# ENCONTRAR PICOS DE MARCADORES SEGÚN LA PREDICCIÓN DE
# QUÉ CÉLULAS SON PERICITOS (SCRNA-SEQ)
# OJO: TARDA MUCHO
#Idents(brain) <- "predicted.id"
#if (!file.exists('da_peaks_marcadores_pericitos_predecidos.Rds')) {
#    da_peaks_predicted <- FindMarkers(
#        object = brain,
#        ident.1 = c("Brain.Mural3","Brain.Mural4","Brain.Mural5","Brain.Mural6"),
#        #ident.2 = c("Brain.EC3","Brain.EC4","Brain.EC5","Brain.EC6"),
#        #ident.3 = c("Brain.Pdgfra1","Brain.Pdgfra2"),
#        test.use = 'LR',
#        latent.vars = 'nCount_peaks'
#    )
#    saveRDS (da_peaks_predicted, 'da_peaks_marcadores_pericitos_predecidos.Rds')
#} else {
#    da_peaks_predicted <- readRDS ('da_peak_marcadores_pericitos_predecidos.Rds')
#}

plot1 <- VlnPlot(
  object = brain,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  #idents = c("L4","L5 IT","L2/3 IT")
)
plot2 <- FeaturePlot(
  object = brain,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  max.cutoff = 'q95'
)
plot1 | plot2

open_l23 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_l456 <- rownames(da_peaks[da_peaks$avg_log2FC < 3, ])
closest_l23 <- ClosestFeature(brain, open_l23)
closest_l456 <- ClosestFeature(brain, open_l456)
