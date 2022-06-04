library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
print(getwd())

source('pipeline_functions.R')

seed = 1992

initialize(seed)

##################### Loading datasets as seurat objects

output_directory = 'seurat_data_objects'
prenatal_directory = paste(output_directory, 'prenatal_objects', sep='/')
postnatal_healthy_directory = paste(output_directory, 'postnatal_healthy_objects', sep='/')
postnatal_diseased_directory = paste(output_directory, 'postnatal_diseased_objects', sep='/')

#prenatal datasets
quijada_e12.5 = LoadH5Seurat(paste(prenatal_directory, 'quijada_e12.5.h5Seurat', sep='/'))
quijada_e16.5 = LoadH5Seurat(paste(prenatal_directory, 'quijada_e16.5.h5Seurat', sep='/'))
jackson_e12.5 = LoadH5Seurat(paste(prenatal_directory, 'jackson_e12.5.h5Seurat', sep='/'))
jackson_e12.5_CKO = LoadH5Seurat(paste(prenatal_directory, 'jackson_e12.5_CKO.h5Seurat', sep='/'))

#postnatal healthy datasets
yuan = LoadH5Seurat(paste(postnatal_healthy_directory, 'yuan.h5Seurat', sep='/'))

forte_healthy_d7_1 = LoadH5Seurat(paste(postnatal_healthy_directory, 'forte_healthy_d7_1.h5Seurat', sep='/'))
forte_healthy_d7_2 = LoadH5Seurat(paste(postnatal_healthy_directory, 'forte_healthy_d7_2.h5Seurat', sep='/'))
forte_healthy_d7_3 = LoadH5Seurat(paste(postnatal_healthy_directory, 'forte_healthy_d7_3.h5Seurat', sep='/'))
forte_healthy_d7_4 = LoadH5Seurat(paste(postnatal_healthy_directory, 'forte_healthy_d7_4.h5Seurat', sep='/'))
forte_healthy_d0   = LoadH5Seurat(paste(postnatal_healthy_directory, 'forte_healthy_d0.h5Seurat'  , sep='/'))
forte_healthy_d7_5 = LoadH5Seurat(paste(postnatal_healthy_directory, 'forte_healthy_d7_5.h5Seurat', sep='/'))

wang_healthy_p8m1 = LoadH5Seurat(paste(postnatal_healthy_directory, 'wang_healthy_p8m1.h5Seurat', sep='/'))
wang_healthy_p8m3 = LoadH5Seurat(paste(postnatal_healthy_directory, 'wang_healthy_p8m3.h5Seurat', sep='/'))

vidal_y1 = LoadH5Seurat(paste(postnatal_healthy_directory, 'vidal_y1.h5Seurat', sep='/'))
vidal_y2 = LoadH5Seurat(paste(postnatal_healthy_directory, 'vidal_y2.h5Seurat', sep='/'))
vidal_y3 = LoadH5Seurat(paste(postnatal_healthy_directory, 'vidal_y3.h5Seurat', sep='/'))

vidal_o1 = LoadH5Seurat(paste(postnatal_healthy_directory, 'vidal_o1.h5Seurat', sep='/'))
vidal_o2 = LoadH5Seurat(paste(postnatal_healthy_directory, 'vidal_o2.h5Seurat', sep='/'))
vidal_o3 = LoadH5Seurat(paste(postnatal_healthy_directory, 'vidal_o3.h5Seurat', sep='/'))

#postnatal diseased datasets
forte_diseased_d1   = LoadH5Seurat(paste(postnatal_diseased_directory, 'forte_diseased_d1.h5Seurat', sep='/'))
forte_diseased_d3_1 = LoadH5Seurat(paste(postnatal_diseased_directory, 'forte_diseased_d3_1.h5Seurat', sep='/'))
forte_diseased_d5   = LoadH5Seurat(paste(postnatal_diseased_directory, 'forte_diseased_d5.h5Seurat', sep='/'))
forte_diseased_d7   = LoadH5Seurat(paste(postnatal_diseased_directory, 'forte_diseased_d7.h5Seurat', sep='/'))
forte_diseased_d14  = LoadH5Seurat(paste(postnatal_diseased_directory, 'forte_diseased_d14.h5Seurat', sep='/'))
forte_diseased_d28  = LoadH5Seurat(paste(postnatal_diseased_directory, 'forte_diseased_d28.h5Seurat', sep='/'))
forte_diseased_d3_2 = LoadH5Seurat(paste(postnatal_diseased_directory, 'forte_diseased_d3_2.h5Seurat', sep='/'))

wang_diseased_p8m1 = LoadH5Seurat(paste(postnatal_diseased_directory, 'wang_diseased_p8m1.h5Seurat', sep='/'))
wang_diseased_p8m3 = LoadH5Seurat(paste(postnatal_diseased_directory, 'wang_diseased_p8m3.h5Seurat', sep='/'))

hesse_epdc1    = LoadH5Seurat(paste(postnatal_diseased_directory, 'hesse_epdc1.h5Seurat', sep='/'))
hesse_epdc2    = LoadH5Seurat(paste(postnatal_diseased_directory, 'hesse_epdc2.h5Seurat', sep='/'))
hesse_MI1_epdc = LoadH5Seurat(paste(postnatal_diseased_directory, 'hesse_MI1_epdc.h5Seurat', sep='/'))
hesse_MI2_epdc = LoadH5Seurat(paste(postnatal_diseased_directory, 'hesse_MI2_epdc.h5Seurat', sep='/'))
hesse_MI3_epdc = LoadH5Seurat(paste(postnatal_diseased_directory, 'hesse_MI3_epdc.h5Seurat', sep='/'))

#################### merging datasets that can be merged: Vidal old, vidal young, 

forte_healthy_d7 = merge(forte_healthy_d7_1, c(
  forte_healthy_d7_2,
  forte_healthy_d7_3,
  forte_healthy_d7_4,
  forte_healthy_d7_5
))

vidal_young = merge(vidal_y1, c(vidal_y2, vidal_y3))
vidal_old = merge(vidal_o1, c(vidal_o2, vidal_o3))

forte_diseased_d3 = merge(forte_diseased_d3_1, forte_diseased_d3_2)

hesse_epdc = merge(hesse_epdc1, hesse_epdc2)
hesse_MI_epdc = merge(hesse_MI1_epdc, c(hesse_MI2_epdc, hesse_MI3_epdc))

#################### generating Quality Control plots for all datasets before QC

# default QC parameters
nFeature_RNA_min_parameter = 200
nFeature_RNA_max_parameter = 6000
percent.mt_parameter = 10

# prenatal
quijada_e12.5.plots.before = create_qc_plots(quijada_e12.5)
quijada_e12.5 = perform_custom_QC(quijada_e12.5, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
quijada_e12.5 = perform_normalization(quijada_e12.5)
quijada_e12.5 = label_cells(quijada_e12.5, 'quijada', 'e12.5', stagename='prenatal')
quijada_e12.5.plots.after = create_qc_plots(quijada_e12.5)

quijada_e16.5.plots.before = create_qc_plots(quijada_e16.5)
quijada_e16.5 = perform_custom_QC(quijada_e16.5, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
quijada_e16.5 = perform_normalization(quijada_e16.5)
quijada_e16.5 = label_cells(quijada_e16.5, 'quijada', 'e16.5', stagename='prenatal')
quijada_e16.5.plots.after = create_qc_plots(quijada_e16.5)

jackson_e12.5.plots.before = create_qc_plots(jackson_e12.5)
jackson_e12.5 = perform_custom_QC(jackson_e12.5, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
jackson_e12.5 = perform_normalization(jackson_e12.5)
jackson_e12.5 = label_cells(jackson_e12.5, 'jackson-weaver', 'e12.5', stagename='prenatal')
jackson_e12.5.plots.after = create_qc_plots(jackson_e12.5)

SaveH5Seurat(quijada_e12.5, paste(prenatal_directory, 'QC_quijada_e12.5.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(quijada_e16.5, paste(prenatal_directory, 'QC_quijada_e16.5.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(jackson_e12.5, paste(prenatal_directory, 'QC_jackson_e12.5.h5Seurat', sep='/'), overwrite=TRUE)

# postnatal healthy

# For YUAN we need to filter the cells based on their barcodes first. this is because there are both WT and CKO cells in the dataset.
# according to the author, the barcodes ending in *-1, *-2 and *-3 come from WT samples

# filter the barcodes here
all_barcodes = colnames(yuan)
all_barcodes_length = lapply(all_barcodes, FUN=function(x) 18 == nchar(x))

all(unlist(all_barcodes_length)) # so they are all length 18

all_barcode_ends = lapply(all_barcodes, FUN=function(x) substring(x, first = 18, last = 18))

all_barcode_ends

all_WT_cells = lapply(all_barcode_ends, FUN=function(x) (x == '1')|(x == '2')|(x == '3'))

only_WT_barcodes = all_barcodes[unlist(all_WT_cells)]

yuan = subset(x = yuan, cells=only_WT_barcodes)

yuan.plots.before = create_qc_plots(yuan)
yuan = perform_custom_QC(yuan, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
yuan = perform_normalization(yuan)
yuan = label_cells(yuan, 'yuan', 120, stagename='postnatal healthy')
yuan.plots.after = create_qc_plots(yuan)

forte_healthy_d0.plots.before = create_qc_plots(forte_healthy_d0)
forte_healthy_d0 = perform_custom_QC(forte_healthy_d0, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_healthy_d0 = perform_normalization(forte_healthy_d0)
forte_healthy_d0 = label_cells(forte_healthy_d0, 'forte', 70, stagename='postnatal healthy')
forte_healthy_d0.plots.after = create_qc_plots(forte_healthy_d0)

forte_healthy_d7.plots.before = create_qc_plots(forte_healthy_d7)
forte_healthy_d7 = perform_custom_QC(forte_healthy_d7, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_healthy_d7 = perform_normalization(forte_healthy_d7)
forte_healthy_d7 = label_cells(forte_healthy_d7, 'forte', 77, stagename='postnatal healthy')
forte_healthy_d7.plots.after = create_qc_plots(forte_healthy_d7)

wang_healthy_p8m1.plots.before = create_qc_plots(wang_healthy_p8m1)
wang_healthy_p8m1 = perform_custom_QC(wang_healthy_p8m1, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
wang_healthy_p8m1 = perform_normalization(wang_healthy_p8m1)
wang_healthy_p8m1 = label_cells(wang_healthy_p8m1, 'wang', 9, stagename='postnatal healthy')
wang_healthy_p8m1.plots.after = create_qc_plots(wang_healthy_p8m1)

wang_healthy_p8m3.plots.before = create_qc_plots(wang_healthy_p8m3)
wang_healthy_p8m3 = perform_custom_QC(wang_healthy_p8m3, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
wang_healthy_p8m3 = perform_normalization(wang_healthy_p8m3)
wang_healthy_p8m3 = label_cells(wang_healthy_p8m3, 'wang', 11, stagename='postnatal healthy')
wang_healthy_p8m3.plots.after = create_qc_plots(wang_healthy_p8m3)

vidal_young.plots.before = create_qc_plots(vidal_young)
vidal_young = perform_custom_QC(vidal_young, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
vidal_young = perform_normalization(vidal_young)
vidal_young = label_cells(vidal_young, 'vidal', 84, stagename='postnatal healthy')
vidal_young.plots.after = create_qc_plots(vidal_young)

vidal_old.plots.before = create_qc_plots(vidal_young)
vidal_old = perform_custom_QC(vidal_old, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
vidal_old = perform_normalization(vidal_old)
vidal_old = label_cells(vidal_old, 'vidal', 600, stagename='postnatal healthy')
vidal_old.plots.after = create_qc_plots(vidal_old)

SaveH5Seurat(yuan, paste(postnatal_healthy_directory, 'QC_yuan.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(forte_healthy_d0, paste(postnatal_healthy_directory, 'QC_forte_healthy_d0.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_healthy_d7, paste(postnatal_healthy_directory, 'QC_forte_healthy_d7.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(wang_healthy_p8m1, paste(postnatal_healthy_directory, 'QC_wang_healthy_p8m1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(wang_healthy_p8m3, paste(postnatal_healthy_directory, 'QC_wang_healthy_p8m3.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(vidal_young, paste(postnatal_healthy_directory, 'QC_vidal_young.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(vidal_old, paste(postnatal_healthy_directory, 'QC_vidal_old.h5Seurat', sep='/'), overwrite=TRUE)

# postnatal diseased

forte_diseased_d1.plots.before = create_qc_plots(forte_diseased_d1)
forte_diseased_d1 = perform_custom_QC(forte_diseased_d1, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_diseased_d1 = perform_normalization(forte_diseased_d1)
forte_diseased_d1 = label_cells(forte_diseased_d1, 'forte', 71, 1, stagename='postnatal diseased')
forte_diseased_d1.plots.after = create_qc_plots(forte_diseased_d1)

forte_diseased_d3.plots.before = create_qc_plots(forte_diseased_d3)
forte_diseased_d3 = perform_custom_QC(forte_diseased_d3, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_diseased_d3 = perform_normalization(forte_diseased_d3)
forte_diseased_d3 = label_cells(forte_diseased_d3, 'forte', 73, 3, stagename='postnatal diseased')
forte_diseased_d3.plots.after = create_qc_plots(forte_diseased_d3)

forte_diseased_d5.plots.before = create_qc_plots(forte_diseased_d5)
forte_diseased_d5 = perform_custom_QC(forte_diseased_d5, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_diseased_d5 = perform_normalization(forte_diseased_d5)
forte_diseased_d5 = label_cells(forte_diseased_d5, 'forte', 75, 5, stagename='postnatal diseased')
forte_diseased_d5.plots.after = create_qc_plots(forte_diseased_d5)

forte_diseased_d7.plots.before = create_qc_plots(forte_diseased_d7)
forte_diseased_d7 = perform_custom_QC(forte_diseased_d7, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_diseased_d7 = perform_normalization(forte_diseased_d7)
forte_diseased_d7 = label_cells(forte_diseased_d7, 'forte', 77, 7, stagename='postnatal diseased')
forte_diseased_d7.plots.after = create_qc_plots(forte_diseased_d7)

forte_diseased_d14.plots.before = create_qc_plots(forte_diseased_d14)
forte_diseased_d14 = perform_custom_QC(forte_diseased_d14, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_diseased_d14 = perform_normalization(forte_diseased_d14)
forte_diseased_d14 = label_cells(forte_diseased_d14, 'forte', 84, 14, stagename='postnatal diseased')
forte_diseased_d14.plots.after = create_qc_plots(forte_diseased_d14)

forte_diseased_d28.plots.before = create_qc_plots(forte_diseased_d28)
forte_diseased_d28 = perform_custom_QC(forte_diseased_d28, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
forte_diseased_d28 = perform_normalization(forte_diseased_d28)
forte_diseased_d28 = label_cells(forte_diseased_d28, 'forte', 98, 28, stagename='postnatal diseased')
forte_diseased_d28.plots.after = create_qc_plots(forte_diseased_d28)

wang_diseased_p8m1.plots.before = create_qc_plots(wang_diseased_p8m1)
wang_diseased_p8m1 = perform_custom_QC(wang_diseased_p8m1, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
wang_diseased_p8m1 = perform_normalization(wang_diseased_p8m1)
wang_diseased_p8m1 = label_cells(wang_diseased_p8m1, 'wang', 9, 1, stagename='postnatal diseased')
wang_diseased_p8m1.plots.after = create_qc_plots(wang_diseased_p8m1)

wang_diseased_p8m3.plots.before = create_qc_plots(wang_diseased_p8m3)
wang_diseased_p8m3 = perform_custom_QC(wang_diseased_p8m3, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
wang_diseased_p8m3 = perform_normalization(wang_diseased_p8m3)
wang_diseased_p8m3 = label_cells(wang_diseased_p8m3, 'wang', 11, 3, stagename='postnatal diseased')
wang_diseased_p8m3.plots.after = create_qc_plots(wang_diseased_p8m3)

hesse_epdc.plots.before = create_qc_plots(hesse_epdc)
hesse_epdc = perform_custom_QC(hesse_epdc, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
hesse_epdc = perform_normalization(hesse_epdc)
hesse_epdc = label_cells(hesse_epdc, 'hesse', 75, 5, stagename='postnatal diseased')
hesse_epdc.plots.after = create_qc_plots(hesse_epdc)

hesse_MI_epdc.plots.before = create_qc_plots(hesse_MI_epdc)
hesse_MI_epdc = perform_custom_QC(hesse_MI_epdc, nFeature_RNA_min_parameter, nFeature_RNA_max_parameter, percent.mt_parameter)
hesse_MI_epdc = perform_normalization(hesse_MI_epdc)
hesse_MI_epdc = label_cells(hesse_MI_epdc, 'hesse', 75, 5, stagename='postnatal diseased')
hesse_MI_epdc.plots.after = create_qc_plots(hesse_MI_epdc)

SaveH5Seurat(forte_diseased_d1, paste(postnatal_diseased_directory, 'QC_forte_diseased_d1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d3, paste(postnatal_diseased_directory, 'QC_forte_diseased_d3.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d5, paste(postnatal_diseased_directory, 'QC_forte_diseased_d5.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d7, paste(postnatal_diseased_directory, 'QC_forte_diseased_d7.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d14, paste(postnatal_diseased_directory, 'QC_forte_diseased_d14.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d28, paste(postnatal_diseased_directory, 'QC_forte_diseased_d28.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(wang_diseased_p8m1, paste(postnatal_diseased_directory, 'QC_wang_diseased_p8m1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(wang_diseased_p8m3, paste(postnatal_diseased_directory, 'QC_wang_diseased_p8m3.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(hesse_epdc, paste(postnatal_diseased_directory, 'QC_hesse_epdc.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(hesse_MI_epdc, paste(postnatal_diseased_directory, 'QC_hesse_MI_epdc.h5Seurat', sep='/'), overwrite=TRUE)


