library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
print(getwd())

source('pipeline_functions.R')

seed = 1992

initialize(seed)
library(Matrix)
################################ Locations of datasets

#prenatal datasets
quijada_e12.5 = 'quijada/quijada_e12.5'
quijada_e16.5 = 'quijada/quijada_e16.5'

jackson_e12.5 = 'jackson/jackson_e12.5'
jackson_e12.5_CKO = 'jackson/jackson_e12.5_CKO'

#postnatal healthy datasets
yuan_barcodes = 'yuan/barcodes.csv'
yuan_counts = 'yuan/counts.mtx'
yuan_genes = 'yuan/genes.txt'

forte_healthy_d7_1 = 'forte/17004-10999'
forte_healthy_d7_2 = 'forte/17004-13609'
forte_healthy_d7_3 = 'forte/17008-11003'
forte_healthy_d7_4 = 'forte/17008-13613'
forte_healthy_d0   = 'forte/17010-12037'
forte_healthy_d7_5 = 'forte/18001-03499'

#wang_healthy_p1m1 = 'wang/postnatal_healthy/wangp1m1' # don't use. doesn't seem to work
#wang_healthy_p1m3 = 'wang/postnatal_healthy/wangp1m3' # don't use. doesn't seem to work

wang_healthy_p8m1 = 'wang/postnatal_healthy/wangp8m1' # don't use. doesn't seem to work
wang_healthy_p8m3 = 'wang/postnatal_healthy/wangp8m3'

vidal_y1 = 'vidal/vidalw12/Y1.csv'
vidal_y2 = 'vidal/vidalw12/Y2.csv'
vidal_y3 = 'vidal/vidalw12/Y3.csv'

vidal_o1 = 'vidal/vidalm18/O1.csv'
vidal_o2 = 'vidal/vidalm18/O2.csv'
vidal_o3 = 'vidal/vidalm18/O3.csv'

#hesse_sham1 =     'hesse/Sham1_CF' # Don't use, contains only fibroblast/stromal cells
#hesse_sham3 =     'hesse/Sham3_CF' # Don't use, contains only fibroblast/stromal cells
#hesse_sham4 =     'hesse/Sham4_CF' # Don't use, contains only fibroblast/stromal cells

#postnatal diseased datasets
forte_diseased_d1   =  'forte/17013-12040'
forte_diseased_d3_1 =  'forte/17014-12041'
forte_diseased_d5   =  'forte/17015-12042'
forte_diseased_d7   =  'forte/17016-12043'
forte_diseased_d14  =  'forte/17017-15642'
forte_diseased_d28  =  'forte/17018-15643'
forte_diseased_d3_2 =  'forte/18002-03500'

#wang_diseased_p1m1 = 'wang/postnatal_diseased/wangp1m1' # Don't use. is too early in the postnatal cycle to classify as postnatal.
#wang_diseased_p1m3 = 'wang/postnatal_diseased/wangp1m3' # Don't use. is too early in the postnatal cycle to classify as postnatal.
wang_diseased_p8m1 = 'wang/postnatal_diseased/wangp8m1'
wang_diseased_p8m3 = 'wang/postnatal_diseased/wangp8m3'

#hesse_acf1 =      'hesse/aCF1_TdTomato' # Don't use, contains only fibroblast/stromal cells
#hesse_acf2 =      'hesse/aCF2_TdTomato' # Don't use, contains only fibroblast/stromal cells
hesse_epdc1 =     'hesse/EPDC1_TdTomato' # NOTE: is ALSO post-MI
hesse_epdc2 =     'hesse/EPDC2_TdTomato' # NOTE: is ALSO post-MI
#hesse_MI1_acf =   'hesse/MI1_aCF' # Don't use, contains only fibroblast/stromal cells
hesse_MI1_epdc =  'hesse/MI1_EPDC'
#hesse_MI2_acf =   'hesse/MI2_aCF' # Don't use, contains only fibroblast/stromal cells
hesse_MI2_epdc =  'hesse/MI2_EPDC'
#hesse_MI3_acf =   'hesse/MI3_aCF' # Don't use, contains only fibroblast/stromal cells
hesse_MI3_epdc =  'hesse/MI3_EPDC'

# loading and converting of datasets

#prenatal datasets
quijada_e12.5 = data_location_to_seurat_object(quijada_e12.5)
quijada_e16.5 = data_location_to_seurat_object(quijada_e16.5)

jackson_e12.5 = data_location_to_seurat_object(jackson_e12.5)
jackson_e12.5_CKO = data_location_to_seurat_object(jackson_e12.5_CKO)

#postnatal healthy
yuan_counts = readMM(yuan_counts)
yuan_cells = read.csv(yuan_barcodes)
colnames(yuan_counts) = yuan_cells[[1]]

yuan_genes = scan(yuan_genes, what = character(), sep = '\n')
rownames(yuan_counts) = yuan_genes
yuan = CreateSeuratObject(yuan_counts)

forte_healthy_d7_1 = data_location_to_seurat_object(forte_healthy_d7_1)
forte_healthy_d7_2 = data_location_to_seurat_object(forte_healthy_d7_2)
forte_healthy_d7_3 = data_location_to_seurat_object(forte_healthy_d7_3)
forte_healthy_d7_4 = data_location_to_seurat_object(forte_healthy_d7_4)
forte_healthy_d0   = data_location_to_seurat_object(forte_healthy_d0)
forte_healthy_d7_5 = data_location_to_seurat_object(forte_healthy_d7_5)

wang_healthy_p8m1 = data_location_to_seurat_object(wang_healthy_p8m1)
wang_healthy_p8m3 = data_location_to_seurat_object(wang_healthy_p8m3)

vidal_y1 = data_location_to_seurat_object(vidal_y1, alt_load_function=TRUE)
vidal_y2 = data_location_to_seurat_object(vidal_y2, alt_load_function=TRUE)
vidal_y3 = data_location_to_seurat_object(vidal_y3, alt_load_function=TRUE)

vidal_o1 = data_location_to_seurat_object(vidal_o1, alt_load_function=TRUE)
vidal_o2 = data_location_to_seurat_object(vidal_o2, alt_load_function=TRUE)
vidal_o3 = data_location_to_seurat_object(vidal_o3, alt_load_function=TRUE)

#postnatal diseased datasets
forte_diseased_d1   = data_location_to_seurat_object(forte_diseased_d1)
forte_diseased_d3_1 = data_location_to_seurat_object(forte_diseased_d3_1)
forte_diseased_d5   = data_location_to_seurat_object(forte_diseased_d5)
forte_diseased_d7   = data_location_to_seurat_object(forte_diseased_d7)
forte_diseased_d14  = data_location_to_seurat_object(forte_diseased_d14)
forte_diseased_d28  = data_location_to_seurat_object(forte_diseased_d28)
forte_diseased_d3_2 = data_location_to_seurat_object(forte_diseased_d3_2)

#wang_diseased_p1m1 = data_location_to_seurat_object(wang_diseased_p1m1)
#wang_diseased_p1m3 = data_location_to_seurat_object(wang_diseased_p1m3)
wang_diseased_p8m1 = data_location_to_seurat_object(wang_diseased_p8m1)
wang_diseased_p8m3 = data_location_to_seurat_object(wang_diseased_p8m3)

hesse_epdc1 =     data_location_to_seurat_object(hesse_epdc1)
hesse_epdc2 =     data_location_to_seurat_object(hesse_epdc2)
hesse_MI1_epdc =  data_location_to_seurat_object(hesse_MI1_epdc)
hesse_MI2_epdc =  data_location_to_seurat_object(hesse_MI2_epdc)
hesse_MI3_epdc =  data_location_to_seurat_object(hesse_MI3_epdc)
print('Finished converting data into seurat objects...\nStarting to ')

# Merging of objects from the same experimental run (for which there can be assumed to be no batch effects)
vidal_young = merge(vidal_y1, c(vidal_y2, vidal_y3), add.cell.ids = c('y1','y2','y3'))
vidal_old = merge(vidal_o1, c(vidal_o2, vidal_o3), add.cell.ids = c('o1','o2','o3'))

# saving of converted data objects
output_directory = 'seurat_data_objects'
prenatal_directory = paste(output_directory, 'prenatal_objects', sep='/')
postnatal_healthy_directory = paste(output_directory, 'postnatal_healthy_objects', sep='/')
postnatal_diseased_directory = paste(output_directory, 'postnatal_diseased_objects', sep='/')
if (!dir.exists(prenatal_directory)){
  dir.create(prenatal_directory)
}
if (!dir.exists(postnatal_healthy_directory)){
  dir.create(postnatal_healthy_directory)
}
if (!dir.exists(postnatal_diseased_directory)){
  dir.create(postnatal_diseased_directory)
}

#prenatal datasets
SaveH5Seurat(quijada_e12.5, paste(prenatal_directory, 'quijada_e12.5.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(quijada_e16.5, paste(prenatal_directory, 'quijada_e16.5.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(jackson_e12.5, paste(prenatal_directory, 'jackson_e12.5.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(jackson_e12.5_CKO, paste(prenatal_directory, 'jackson_e12.5_CKO.h5Seurat', sep='/'), overwrite=TRUE)


#postnatal healthy datasets
SaveH5Seurat(yuan, paste(postnatal_healthy_directory, 'yuan.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(forte_healthy_d7_1, paste(postnatal_healthy_directory, 'forte_healthy_d7_1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_healthy_d7_2, paste(postnatal_healthy_directory, 'forte_healthy_d7_2.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_healthy_d7_3, paste(postnatal_healthy_directory, 'forte_healthy_d7_3.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_healthy_d7_4, paste(postnatal_healthy_directory, 'forte_healthy_d7_4.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_healthy_d0, paste(postnatal_healthy_directory, 'forte_healthy_d0.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_healthy_d7_5, paste(postnatal_healthy_directory, 'forte_healthy_d7_5.h5Seurat', sep='/'), overwrite=TRUE)

#SaveH5Seurat(wang_healthy_p1m1, paste(postnatal_healthy_directory, 'wang_healthy_p1m1.h5Seurat', sep='/'), overwrite=TRUE)
#SaveH5Seurat(wang_healthy_p1m3, paste(postnatal_healthy_directory, 'wang_healthy_p1m3.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(wang_healthy_p8m1, paste(postnatal_healthy_directory, 'wang_healthy_p8m1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(wang_healthy_p8m3, paste(postnatal_healthy_directory, 'wang_healthy_p8m3.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(vidal_y1, paste(postnatal_healthy_directory, 'vidal_y1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(vidal_y2, paste(postnatal_healthy_directory, 'vidal_y2.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(vidal_y3, paste(postnatal_healthy_directory, 'vidal_y3.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(vidal_o1, paste(postnatal_healthy_directory, 'vidal_o1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(vidal_o2, paste(postnatal_healthy_directory, 'vidal_o2.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(vidal_o3, paste(postnatal_healthy_directory, 'vidal_o3.h5Seurat', sep='/'), overwrite=TRUE)

#postnatal diseased datasets
SaveH5Seurat(forte_diseased_d1, paste(postnatal_diseased_directory, 'forte_diseased_d1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d3_1, paste(postnatal_diseased_directory, 'forte_diseased_d3_1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d5, paste(postnatal_diseased_directory, 'forte_diseased_d5.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d7, paste(postnatal_diseased_directory, 'forte_diseased_d7.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d14, paste(postnatal_diseased_directory, 'forte_diseased_d14.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d28, paste(postnatal_diseased_directory, 'forte_diseased_d28.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(forte_diseased_d3_2, paste(postnatal_diseased_directory, 'forte_diseased_d3_2.h5Seurat', sep='/'), overwrite=TRUE)

#SaveH5Seurat(wang_diseased_p1m1, paste(postnatal_diseased_directory, 'wang_diseased_p1m1.h5Seurat', sep='/'), overwrite=TRUE)
#SaveH5Seurat(wang_diseased_p1m3, paste(postnatal_diseased_directory, 'wang_diseased_p1m3.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(wang_diseased_p8m1, paste(postnatal_diseased_directory, 'wang_diseased_p8m1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(wang_diseased_p8m3, paste(postnatal_diseased_directory, 'wang_diseased_p8m3.h5Seurat', sep='/'), overwrite=TRUE)

SaveH5Seurat(hesse_epdc1, paste(postnatal_diseased_directory, 'hesse_epdc1.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(hesse_epdc2, paste(postnatal_diseased_directory, 'hesse_epdc2.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(hesse_MI1_epdc, paste(postnatal_diseased_directory, 'hesse_MI1_epdc.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(hesse_MI2_epdc, paste(postnatal_diseased_directory, 'hesse_MI2_epdc.h5Seurat', sep='/'), overwrite=TRUE)
SaveH5Seurat(hesse_MI3_epdc, paste(postnatal_diseased_directory, 'hesse_MI3_epdc.h5Seurat', sep='/'), overwrite=TRUE)









