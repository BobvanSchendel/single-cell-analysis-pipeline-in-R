library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
print(getwd())

source('pipeline_functions.R')

seed = 1992

initialize(seed)

# epicardial selection has to be performed on the following datasets


##################### Loading datasets as seurat objects and finding dimensionality

output_directory = 'seurat_data_objects'
prenatal_directory = paste(output_directory, 'prenatal_objects', sep='/')
postnatal_healthy_directory = paste(output_directory, 'postnatal_healthy_objects', sep='/')
postnatal_diseased_directory = paste(output_directory, 'postnatal_diseased_objects', sep='/')

#prenatal datasets
jackson_e12.5 = LoadH5Seurat(paste(prenatal_directory, 'QC_jackson_e12.5.h5Seurat', sep='/'))

#postnatal healthy datasets
forte_healthy_d0   = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_forte_healthy_d0.h5Seurat'  , sep='/'))
forte_healthy_d7 = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_forte_healthy_d7.h5Seurat', sep='/'))

wang_healthy_p8m1 = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_wang_healthy_p8m1.h5Seurat', sep='/'))
wang_healthy_p8m3 = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_wang_healthy_p8m3.h5Seurat', sep='/'))

vidal_young = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_vidal_young.h5Seurat', sep='/'))
vidal_old = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_vidal_old.h5Seurat', sep='/'))

#postnatal diseased datasets
forte_diseased_d1   = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d1.h5Seurat', sep='/'))
forte_diseased_d3 = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d3.h5Seurat', sep='/'))
forte_diseased_d5   = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d5.h5Seurat', sep='/'))
forte_diseased_d7   = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d7.h5Seurat', sep='/'))
forte_diseased_d14  = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d14.h5Seurat', sep='/'))
forte_diseased_d28  = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d28.h5Seurat', sep='/'))

wang_diseased_p8m1 = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_wang_diseased_p8m1.h5Seurat', sep='/'))
wang_diseased_p8m3 = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_wang_diseased_p8m3.h5Seurat', sep='/'))



jackson_e12.5 = RunPCA(jackson_e12.5, features = VariableFeatures(jackson_e12.5), verbose=F)
jackson_e12.5.elbow = ElbowPlot(jackson_e12.5)
jackson_e12.5.dimensionality = 5
# dim 5

forte_healthy_d0 = RunPCA(forte_healthy_d0, features = VariableFeatures(forte_healthy_d0), verbose=F)
forte_healthy_d7 = RunPCA(forte_healthy_d7, features = VariableFeatures(forte_healthy_d7), verbose=F)
forte_healthy_d0.elbow = ElbowPlot(forte_healthy_d0)
forte_healthy_d0.dimensionality = 5
# dim 5

forte_healthy_d7.elbow = ElbowPlot(forte_healthy_d7)
forte_healthy_d7.dimensionality = 5
# dim 5




wang_healthy_p8m1 = RunPCA(wang_healthy_p8m1, features = VariableFeatures(wang_healthy_p8m1), verbose=F)
wang_healthy_p8m3 = RunPCA(wang_healthy_p8m3, features = VariableFeatures(wang_healthy_p8m3), verbose=F)
wang_healthy_p8m1.elbow = ElbowPlot(wang_healthy_p8m1)
wang_healthy_p8m1.dimensionality = 5
# dim 5

wang_healthy_p8m3.elbow = ElbowPlot(wang_healthy_p8m3)
wang_healthy_p8m3.dimensionality = 10
# dim 10



vidal_young = RunPCA(vidal_young, features = VariableFeatures(vidal_young), verbose=F)
vidal_old = RunPCA(vidal_old, features = VariableFeatures(vidal_old), verbose=F)
vidal_young.elbow = ElbowPlot(vidal_young)
vidal_young.dimensionality = 10
# dim 10

vidal_old.elbow = ElbowPlot(vidal_old)
vidal_old.dimensionality = 10
# dim 10

forte_diseased_d1 = RunPCA(forte_diseased_d1, features = VariableFeatures(forte_diseased_d1), verbose=F)
forte_diseased_d3 = RunPCA(forte_diseased_d3, features = VariableFeatures(forte_diseased_d3), verbose=F)
forte_diseased_d5 = RunPCA(forte_diseased_d5, features = VariableFeatures(forte_diseased_d5), verbose=F)
forte_diseased_d7 = RunPCA(forte_diseased_d7, features = VariableFeatures(forte_diseased_d7), verbose=F)
forte_diseased_d14 = RunPCA(forte_diseased_d14, features = VariableFeatures(forte_diseased_d14), verbose=F)
forte_diseased_d28 = RunPCA(forte_diseased_d28, features = VariableFeatures(forte_diseased_d28), verbose=F)
forte_diseased_d1.elbow = ElbowPlot(forte_diseased_d1)
forte_diseased_d1.dimensionality = 4
# dim 4

forte_diseased_d3.elbow = ElbowPlot(forte_diseased_d3)
forte_diseased_d3.dimensionality = 8
# dim 8

forte_diseased_d5.elbow = ElbowPlot(forte_diseased_d5)
forte_diseased_d5.dimensionality = 8
# dim 8

forte_diseased_d7.elbow = ElbowPlot(forte_diseased_d7)
forte_diseased_d7.dimensionality = 8
# dim 8

forte_diseased_d14.elbow = ElbowPlot(forte_diseased_d14)
forte_diseased_d14.dimensionality = 5
# dim 5

forte_diseased_d28.elbow = ElbowPlot(forte_diseased_d28)
forte_diseased_d28.dimensionality = 4
# dim 4


wang_diseased_p8m1 = RunPCA(wang_diseased_p8m1, features = VariableFeatures(wang_diseased_p8m1), verbose=F)
wang_diseased_p8m3 = RunPCA(wang_diseased_p8m3, features = VariableFeatures(wang_diseased_p8m3), verbose=F)
wang_diseased_p8m1.elbow = ElbowPlot(wang_diseased_p8m1)
wang_diseased_p8m1.dimensionality = 5
# dim 5

wang_diseased_p8m3.elbow = ElbowPlot(wang_diseased_p8m3)
wang_diseased_p8m3.dimensionality = 6
# dim 6

############################# Clustering and epicardial marker gene selection

#prenatal
jackson_e12.5.resolution = 0.5
jackson_e12.5 = perform_clustering(jackson_e12.5, dimensionality = jackson_e12.5.dimensionality, resolution = jackson_e12.5.resolution)

#postnatal healthy
forte_healthy_d0.resolution = 0.5
forte_healthy_d0 = perform_clustering(forte_healthy_d0, dimensionality = forte_healthy_d0.dimensionality, resolution = forte_healthy_d0.resolution)
forte_healthy_d7.resolution = 0.5
forte_healthy_d7 = perform_clustering(forte_healthy_d7, dimensionality = forte_healthy_d0.dimensionality, resolution = forte_healthy_d7.resolution)

wang_healthy_p8m1.resolution = 0.5
wang_healthy_p8m1 = perform_clustering(wang_healthy_p8m1, dimensionality = wang_healthy_p8m1.dimensionality, resolution = wang_healthy_p8m1.resolution)
wang_healthy_p8m3.resolution = 0.5
wang_healthy_p8m3 = perform_clustering(wang_healthy_p8m3, dimensionality = wang_healthy_p8m3.dimensionality, resolution = wang_healthy_p8m3.resolution)

vidal_young.resolution = 0.5
vidal_young = perform_clustering(vidal_young, dimensionality = vidal_young.dimensionality, resolution = vidal_young.resolution)
vidal_old.resolution = 0.5
vidal_old = perform_clustering(vidal_old, dimensionality = vidal_old.dimensionality, resolution = vidal_old.resolution)

#postnatal diseased
forte_diseased_d1.resolution = 0.5
forte_diseased_d1 = perform_clustering(forte_diseased_d1, dimensionality = forte_diseased_d1.dimensionality, resolution = forte_diseased_d1.resolution)
forte_diseased_d3.resolution = 0.5
forte_diseased_d3 = perform_clustering(forte_diseased_d3, dimensionality = forte_diseased_d3.dimensionality, resolution = forte_diseased_d3.resolution)
forte_diseased_d5.resolution = 0.5
forte_diseased_d5 = perform_clustering(forte_diseased_d5, dimensionality = forte_diseased_d5.dimensionality, resolution = forte_diseased_d5.resolution)
forte_diseased_d7.resolution = 0.5
forte_diseased_d7 = perform_clustering(forte_diseased_d7, dimensionality = forte_diseased_d7.dimensionality, resolution = forte_diseased_d7.resolution)
forte_diseased_d14.resolution = 0.5
forte_diseased_d14 = perform_clustering(forte_diseased_d14, dimensionality = forte_diseased_d14.dimensionality, resolution = forte_diseased_d14.resolution)
forte_diseased_d28.resolution = 0.5
forte_diseased_d28 = perform_clustering(forte_diseased_d28, dimensionality = forte_diseased_d28.dimensionality, resolution = forte_diseased_d28.resolution)

wang_diseased_p8m1.resolution = 0.5
wang_diseased_p8m1 = perform_clustering(wang_diseased_p8m1, dimensionality = wang_diseased_p8m1.dimensionality, resolution = wang_diseased_p8m1.resolution)
wang_diseased_p8m3.resolution = 0.5
wang_diseased_p8m3 = perform_clustering(wang_diseased_p8m3, dimensionality = wang_diseased_p8m3.dimensionality, resolution = wang_diseased_p8m3.resolution)



combined_epicardial_marker_genes = c('Wt1','Upk3b','Cebpb','Krt18','Msln','Clu','Dmkn')
combined_epicardial_marker_genes_without_Upk3b = c('Wt1','Cebpb','Krt18','Msln','Clu','Dmkn')

# IMPORTANT: Tcf21 is an inverse marker gene for epicardial cells.

DefaultAssay(jackson_e12.5) = 'RNA'
DefaultAssay(forte_healthy_d0) = 'RNA'
DefaultAssay(forte_healthy_d7) = 'RNA'
DefaultAssay(wang_healthy_p8m1) = 'RNA'
DefaultAssay(wang_healthy_p8m3) = 'RNA'
DefaultAssay(vidal_young) = 'RNA'
DefaultAssay(vidal_old) = 'RNA'
DefaultAssay(forte_diseased_d1) = 'RNA'
DefaultAssay(forte_diseased_d3) = 'RNA'
DefaultAssay(forte_diseased_d5) = 'RNA'
DefaultAssay(forte_diseased_d7) = 'RNA'
DefaultAssay(forte_diseased_d14) = 'RNA'
DefaultAssay(forte_diseased_d28) = 'RNA'
DefaultAssay(wang_diseased_p8m1) = 'RNA'
DefaultAssay(wang_diseased_p8m3) = 'RNA'

# these plots are saved in the plots/epicardial selection folder
jackson_e12.5.selectionplots = DotPlot(jackson_e12.5, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(jackson_e12.5, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
forte_healthy_d0.selectionplots = DotPlot(forte_healthy_d0, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_healthy_d0, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
forte_healthy_d7.selectionplots = DotPlot(forte_healthy_d7, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_healthy_d7, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
wang_healthy_p8m1.selectionplots = DotPlot(wang_healthy_p8m1, features = combined_epicardial_marker_genes_without_Upk3b, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(wang_healthy_p8m1, features = combined_epicardial_marker_genes_without_Upk3b, pt.size=1, stack=TRUE, flip=FALSE)
wang_healthy_p8m3.selectionplots = DotPlot(wang_healthy_p8m3, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(wang_healthy_p8m3, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
vidal_young.selectionplots = DotPlot(vidal_young, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(vidal_young, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
vidal_old.selectionplots = DotPlot(vidal_old, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(vidal_old, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)

forte_diseased_d1.selectionplots = DotPlot(forte_diseased_d1, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_diseased_d1, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
forte_diseased_d3.selectionplots = DotPlot(forte_diseased_d3, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_diseased_d3, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
forte_diseased_d5.selectionplots = DotPlot(forte_diseased_d5, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_diseased_d5, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
forte_diseased_d7.selectionplots = DotPlot(forte_diseased_d7, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_diseased_d7, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
forte_diseased_d14.selectionplots = DotPlot(forte_diseased_d14, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_diseased_d14, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
forte_diseased_d28.selectionplots = DotPlot(forte_diseased_d28, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(forte_diseased_d28, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
wang_diseased_p8m1.selectionplots = DotPlot(wang_diseased_p8m1, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(wang_diseased_p8m1, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)
wang_diseased_p8m3.selectionplots = DotPlot(wang_diseased_p8m3, features = combined_epicardial_marker_genes, cols=c('black', 'green'), dot.scale=10) +
  VlnPlot(wang_diseased_p8m3, features = combined_epicardial_marker_genes, pt.size=1, stack=TRUE, flip=FALSE)



#epicardial cluster selection and saving object to disk

jackson_e12.5_epicardial = subset(x = jackson_e12.5, idents=c(2, 12, 13))

forte_healthy_d0_epicardial = subset(x = forte_healthy_d0, idents=c(4, 7, 8, 10))
forte_healthy_d7_epicardial = subset(x = forte_healthy_d7, idents=c(8, 11))

wang_healthy_p8m1_epicardial = subset(x = wang_healthy_p8m1, idents=c(1, 2, 3, 6))
wang_healthy_p8m3_epicardial = subset(x = wang_healthy_p8m3, idents=c(0, 2, 3))

vidal_young_epicardial = subset(x = vidal_young, idents=c(1, 6, 8, 10))
vidal_old_epicardial = subset(x = vidal_old, idents=c(1, 6, 7, 10))

forte_diseased_d1_epicardial = subset(x = forte_diseased_d1, idents=c(1, 3, 6, 11))
forte_diseased_d3_epicardial = subset(x = forte_diseased_d3, idents=c(6, 7, 14))
forte_diseased_d5_epicardial = subset(x = forte_diseased_d5, idents=c(2, 6, 8))
forte_diseased_d7_epicardial = subset(x = forte_diseased_d7, idents=c(2, 7, 12, 13))
forte_diseased_d14_epicardial = subset(x = forte_diseased_d14, idents=c(4, 10, 15))
forte_diseased_d28_epicardial = subset(x = forte_diseased_d28, idents=c(4, 9))

wang_diseased_p8m1_epicardial = subset(x = forte_diseased_d28, idents=c(0, 1, 4, 5))
wang_diseased_p8m3_epicardial = subset(x = forte_diseased_d28, idents=c(4, 5))
