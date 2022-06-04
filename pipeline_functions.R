
# Run these installs once, then comment all the install lines
# Uncomment to install

#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#remotes::install_github("mojaveazure/seurat-disk")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#install.packages('ape')
#install.packages('SeuratDisk')

#BiocManager::install('Seurat')
#BiocManager::install("SingleCellExperiment")
#BiocManager::install('scRNAseq')
#BiocManager::install('scater')
#BiocManager::install('flexmix')
#BiocManager::install('splines')
#BiocManager::install('biomaRt')
#BiocManager::install('miQC')
#BiocManager::install('limma')
#BiocManager::install("glmGamPoi")
#devtools::install_github("satijalab/sctransform", ref = "develop")
#devtools::install_github('satijalab/seurat-data')

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scRNAseq)
  library(scater)
  library(flexmix)
  library(splines)
  library(biomaRt)
  library(miQC)
  library(Seurat)
  library(SeuratData)
  library(SeuratDisk)
  library(patchwork)
  library(dplyr)
  library(rstudioapi)
  library(gridExtra)
})

######### Global constants
# for object loading
min_cells_parameter = 3
min_features_parameter = 200
variable_features_parameter = 2000

# for QC
nFeature_RNA_min_parameter = 200
nFeature_RNA_max_parameter = 6000
percent.mt_parameter = 5



########### helper function for finding similar gene names in an entire list searchpattern is a simple regex
#e.g.: get_matches(rownames(seurat.obj.name), "[Mm][Tt]-")
get_matches = function(searchlist, searchpattern){
  searchlist[vapply(searchlist, function(x) all(grepl(searchpattern, x)), logical(1))]
}

############################## LOADING OF DATASETS

load_dataset_from_location = function(location){
  dataset = Read10X(data.dir = location)
  return(dataset)
}

# Use this one for the Vidal datasets
alt_load_dataset_from_location = function(location){
  dataset = read.table(file=location, header=T, row.names=1, sep=',', as.is=T)
  return(dataset)
}

create_seurat_object_from_data = function(dataset){
  seurat.obj = CreateSeuratObject(counts = dataset, min.cells = min_cells_parameter, min.features = min_features_parameter)
  return(seurat.obj)
}

data_location_to_seurat_object = function(location, alt_load_function = FALSE){
  data = NULL
  print(paste('Starting to load object at location: ', location, sep=''))
  if (alt_load_function){
    data = alt_load_dataset_from_location(location)
  }else{
    data = load_dataset_from_location(location)
  }
  data = create_seurat_object_from_data(data)
  return(data)
}

##############################  quality controls

perform_custom_QC = function(seurat.obj, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter){
  seurat.obj[['percent.mt']] = PercentageFeatureSet(seurat.obj, pattern = "^[Mm][Tt]-", assay='RNA') # mark the percentage of (reads mapping to) mitochondrial genes per cell
  seurat.obj = subset(seurat.obj, subset = nFeature_RNA > nFeature_RNA_min_parameter & nFeature_RNA < nFeature_RNA_max_parameter & percent.mt < percent.mt_parameter)
  return(seurat.obj)
}

perform_QC = function(seurat.obj){
  seurat.obj[['percent.mt']] = PercentageFeatureSet(seurat.obj, pattern = "^[Mm][Tt]-", assay='RNA') # mark the percentage of (reads mapping to) mitochondrial genes per cell
  seurat.obj = subset(seurat.obj, subset = nFeature_RNA > nFeature_RNA_min_parameter & nFeature_RNA < nFeature_RNA_max_parameter & percent.mt < percent.mt_parameter)
  return(seurat.obj)
}

create_qc_plots = function(seurat.obj){
  seurat.obj[['percent.mt']] = PercentageFeatureSet(seurat.obj, pattern = "^[Mm][Tt]-")
  violin_feature_count_percent = VlnPlot(seurat.obj, features = c('nFeature_RNA','nCount_RNA', 'percent.mt'), ncol = 3)
  count_mt_plot = FeatureScatter(seurat.obj, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
  count_feature_plot = FeatureScatter(seurat.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
  return(list(violin_feature_count_percent, count_mt_plot, count_feature_plot))
}
############################## Clustering

perform_normalization = function(seurat.obj){
  seurat.obj = SCTransform(seurat.obj, verbose=FALSE, vars.to.regress = 'percent.mt', vst.flavor = 'v2')
  return(seurat.obj)
}

perform_clustering = function(seurat.obj, dimensionality, resolution, reduction_method = 'UMAP'){
  seurat.obj = RunPCA(seurat.obj, features = VariableFeatures(seurat.obj), verbose=F)
  seurat.obj = FindNeighbors(object = seurat.obj, dims = 1:dimensionality, verbose=F)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution, verbose=F)
  if (reduction_method == 'tSNE'){# in tsne we use resolution as perplexity
    seurat.obj = RunTSNE(seurat.obj, seed.use = magic_seed_number, perplexity = resolution, dims = 1:dimensionality, tsne.method='Rtsne')
  } else if (reduction_method == 'UMAP'){
    seurat.obj = RunUMAP(seurat.obj, seed.use = magic_seed_number, dims = 1:dimensionality, verbose=F)
  } else {
    print(paste('reduction_method argument (', reduction_method,') in perform_clustering not recognized, using UMAP.', sep=''))
    seurat.obj = RunUMAP(seurat.obj, seed.use = magic_seed_number, dims = 1:dimensionality, verbose=F)
  } 
  return(seurat.obj)
}

calculate_all_marker_genes = function(seurat.obj){
  all_marker_genes = FindAllMarkers(seurat.obj,
                                           only.pos=T,
                                           min.pct=0.25,
                                           logfc.threshold=0.25,
                                            verbose=F) %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
  return(all_marker_genes)
}

plot_clusters = function(seurat.obj){
  clusterplot = DimPlot(seurat.obj, pt.size=1, label=T) + ggplot2::theme_dark()
  return(clusterplot)
}

calculate_and_save_clusters_and_heatmap = function(dataset, path, name, epicardial_factors, all_marker_genes){
  top10_genes_per_cluster = all_marker_genes %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
  png(width=1200, height=1200)
  clusterpath = paste(path, 'clusters_', name, '.png', sep='')
  generated_clusters = plot_clusters(dataset)
  ggsave(clusterpath)
  dev.off()
  png(width=1200, height=1200)
  heatmappath = paste(path, 'heatmap_', name, '.png', sep='')
  generated_heatmap = DoHeatmap(dataset, features = top10_genes_per_cluster$gene)
  ggsave(heatmappath)
  dev.off()
  png(width=1200, height=1200)
  featurepath = paste(path, 'epi.activity_', name, '.png', sep='')
  cluster_feature_plot = FeaturePlot(dataset, features = epicardial_factors, cols=c('lightgrey','red'), pt.size=1, order=T)+generated_clusters
  ggsave(featurepath)
  dev.off()
  return(list(generated_clusters, generated_heatmap, cluster_feature_plot))
}

calculate_and_save_violin_plots = function(dataset, path, name, epicardial_factors){
  png(width=1200, height=1200)
  violinpath = paste(path, 'violin_plots_', name, '.png', sep='')
  violinplot = VlnPlot(dataset, features = epicardial_factors, stack=T, flip=T) +
    theme(legend.position = 'none') +ggtitle(name)
  ggsave(violinpath)
  dev.off()
  
  tablepath = paste(path, 'clustertable_', name, '.png', sep='')
  png(tablepath, width=1200, height=1200)
  clustertable = table(dataset@meta.data$seurat_clusters)
  tableplot = plot(clustertable)
  ggsave(tablepath, tableplot)
  dev.off()
  
  return(violinplot)
}

############################### entire analysis pipeline

# loads the objects in the location list and integrates them, returning a single object
integrate_pipeline = function(location.list, assumed_dimensionality, clustering_resolution, reduction_method = 'UMAP', alt_load_function=F){
  dataset.list = lapply(X=location.list, FUN=function(dataset_location){
    print('Starting loading of dataset at location:')
    print(dataset_location)
    if (alt_load_function){
      dataset = alt_load_dataset_from_location(dataset_location)
    } else {
      dataset = load_dataset_from_location(dataset_location)
    }
    print('Starting converting of dataset...')
    dataset = create_seurat_object_from_data(dataset)
    print('Finished loading and converting of dataset')
    return(dataset)
  })
  print('Starting QC of datasets')
  dataset.list = lapply(X=dataset.list, FUN = perform_QC)
  print('Finished QC of datasets')
  
  print('Starting normalization of datasets')
  dataset.list = lapply(X=dataset.list, FUN = perform_normalization)
  print('Finished normalization of datasets')
  
  print('Starting integration of datasets')
  features = SelectIntegrationFeatures(dataset.list, nfeatures=3000)
  dataset.list = PrepSCTIntegration(dataset.list, anchor.features = features)
  integration.anchors = FindIntegrationAnchors(dataset.list, normalization.method = 'SCT', anchor.features = features)
  dataset.combined = IntegrateData(anchorset = integration.anchors, normalization.method = 'SCT')
  print('Finished integration of datasets')
  return(dataset.combined)
}

# analyses the object given, clusters it and returns the resulting object
dataset_analysis_pipeline = function(dataset, assumed_dimensionality, clustering_resolution, reduction_method = 'UMAP'){
  print('Starting QC of datasets')
  dataset = perform_QC(dataset)
  print('Finished QC of datasets')
  
  print('Starting normalization of datasets')
  dataset = perform_normalization(dataset)
  print('Finished normalization of datasets')
  
  print('Starting clustering of data')
  dataset = perform_clustering(dataset, assumed_dimensionality, clustering_resolution, reduction_method = reduction_method)
  print('Finished clustering')
  return(dataset)
}

# loads and analyses a single object at a given location.
analysis_pipeline = function(dataset_location, assumed_dimensionality, clustering_resolution, reduction_method = 'UMAP', alt_load_function=FALSE){
  print('Starting loading of dataset at location:')
  print(dataset_location)
  if (alt_load_function){
    dataset = alt_load_dataset_from_location(dataset_location)
  } else {
    dataset = load_dataset_from_location(dataset_location)
  }
  print('Starting converting of dataset...')
  dataset = create_seurat_object_from_data(dataset)
  print('Finished loading and converting of dataset')
  dataset = dataset_analysis_pipeline(dataset, assumed_dimensionality, clustering_resolution, reduction_method)
  return(dataset)
}

label_cells = function(dataset, sourcename, age, dpi=0, stagename){
  dataset$source = sourcename
  dataset$age = age
  dataset$dpi = dpi
  dataset$stage = stagename
  return(dataset)
}

initialize = function(seed){
  magic_seed_number = seed
  set.seed(magic_seed_number)
  
  setwd(dirname(getActiveDocumentContext()$path))
  print(getwd())
}