library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
print(getwd())

source('pipeline_functions.R')

seed = 1992

initialize(seed)

################### Loading of the datasets

output_directory = 'seurat_data_objects'
prenatal_directory = paste(output_directory, 'prenatal_objects', sep='/')
postnatal_healthy_directory = paste(output_directory, 'postnatal_healthy_objects', sep='/')
postnatal_diseased_directory = paste(output_directory, 'postnatal_diseased_objects', sep='/')

#prenatal datasets
quijada_e12.5 = LoadH5Seurat(paste(prenatal_directory, 'QC_quijada_e12.5.h5Seurat', sep='/'))
quijada_e16.5 = LoadH5Seurat(paste(prenatal_directory, 'QC_quijada_e16.5.h5Seurat', sep='/'))

jackson_e12.5 = LoadH5Seurat(paste(prenatal_directory, 'QC_jackson_e12.5_epicardial.h5Seurat', sep='/'))

#postnatal healthy datasets
yuan = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_yuan.h5Seurat', sep='/'))

forte_healthy_d0   = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_forte_healthy_d0_epicardial.h5Seurat'  , sep='/'))
forte_healthy_d7 = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_forte_healthy_d7_epicardial.h5Seurat', sep='/'))

wang_healthy_p8m1 = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_wang_healthy_p8m1_epicardial.h5Seurat', sep='/'))
wang_healthy_p8m3 = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_wang_healthy_p8m3_epicardial.h5Seurat', sep='/'))

vidal_young = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_vidal_young_epicardial.h5Seurat', sep='/'))
vidal_old = LoadH5Seurat(paste(postnatal_healthy_directory, 'QC_vidal_old_epicardial.h5Seurat', sep='/'))

#postnatal diseased datasets
forte_diseased_d1   = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d1_epicardial.h5Seurat', sep='/'))
forte_diseased_d3 = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d3_epicardial.h5Seurat', sep='/'))
forte_diseased_d5   = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d5_epicardial.h5Seurat', sep='/'))
forte_diseased_d7   = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d7_epicardial.h5Seurat', sep='/'))
forte_diseased_d14  = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d14_epicardial.h5Seurat', sep='/'))
forte_diseased_d28  = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_forte_diseased_d28_epicardial.h5Seurat', sep='/'))

wang_diseased_p8m1 = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_wang_diseased_p8m1_epicardial.h5Seurat', sep='/'))
wang_diseased_p8m3 = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_wang_diseased_p8m3_epicardial.h5Seurat', sep='/'))

hesse_epdc    = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_hesse_epdc.h5Seurat', sep='/'))
hesse_MI_epdc = LoadH5Seurat(paste(postnatal_diseased_directory, 'QC_hesse_MI_epdc.h5Seurat', sep='/'))

integration_function = function(datasets, k.weight=0){
  integration_features = SelectIntegrationFeatures(datasets, nfeatures=3000)
  datasets = PrepSCTIntegration(datasets, anchor.features = integration_features)
  integration_anchors = FindIntegrationAnchors(datasets, normalization.method = 'SCT', anchor.features = integration_features)
  if (k.weight == 0){
    integrated_datasets = IntegrateData(anchorset = integration_anchors, normalization.method = 'SCT')
  }else{
    integrated_datasets = IntegrateData(anchorset = integration_anchors, normalization.method = 'SCT', k.weight = k.weight)
  }
  return(integrated_datasets)
}

quijada_complete = integration_function(list(quijada_e12.5, quijada_e16.5))

jackson_integrated = jackson_e12.5

yuan_integrated = yuan

forte_integrated = integration_function(list(forte_healthy_d0,
                                           forte_healthy_d7,
                                           forte_diseased_d1,
                                           forte_diseased_d3,
                                           forte_diseased_d5,
                                           forte_diseased_d7,
                                           forte_diseased_d14,
                                           forte_diseased_d28))

DefaultAssay(wang_healthy_p8m1) = 'SCT'
DefaultAssay(wang_healthy_p8m3) = 'SCT'
DefaultAssay(wang_diseased_p8m1) = 'SCT'
DefaultAssay(wang_diseased_p8m3) = 'SCT'

wang_integrated = integration_function(list(wang_healthy_p8m1,
                                          wang_healthy_p8m3,
                                          wang_diseased_p8m1,
                                          wang_diseased_p8m3))

vidal_integrated = integration_function(list(vidal_young,
                                           vidal_old))

hesse_integrated = integration_function(list(hesse_epdc, hesse_MI_epdc))

all_integrated = integration_function(list(quijda_complete,
                                           jackson_integrated,
                                           yuan_integrated,
                                           forte_integrated,
                                           wang_integrated,
                                           vidal_integrated,
                                           hesse_integrated))

SaveH5Seurat(all_integrated, paste(prenatal_directory, 'all_integrated.h5Seurat', sep='/'), overwrite=TRUE)

















































