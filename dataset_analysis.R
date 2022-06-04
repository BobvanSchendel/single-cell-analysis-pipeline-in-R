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

all_integrated = LoadH5Seurat(paste(output_directory, 'post_all_integrated.h5Seurat', sep='/'))

all_integrated

print('printing number of NA cells:')
print(sum(is.na(all_integrated@assays$RNA@data)))

all_integrated@assays$integrated@data

# load list of neurotrophic factors

DefaultAssay(all_integrated) = 'SCT'

DimPlot(all_integrated, group.by='stage', label=T) + labs(title='Integrated datasets')

# Identifying possible marker genes differentially expressed AND also being a neurotrophic factor
neurotrophic_factors = read.csv('neurotrophic_factors.csv', header=F, fileEncoding='UTF-8-BOM')
neurotrophic_factors = unlist(neurotrophic_factors)

# some neurotrophic factors are not present in the sequenced objects so we have to filter them out or FindAllMarkers throws an error
sequenced_neurotrophic_factors = intersect(neurotrophic_factors, rownames(all_integrated))

Idents(all_integrated) = 'stage'
all_integrated = PrepSCTFindMarkers(all_integrated, assay = 'SCT')

SaveH5Seurat(all_integrated, paste(output_directory, 'post_all_integrated.h5Seurat', sep='/'), overwrite=T)
all_stages_markers = FindAllMarkers(object=all_integrated, features=sequenced_neurotrophic_factors, logfc.threshold = 0.25)

all_stages_markers
# p_val : p_val (unadjusted)
# avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.
# pct.1 : The percentage of cells where the gene is detected in the first group
# pct.2 : The percentage of cells where the gene is detected in the second group
# p_val_adj : Adjusted p-value, based on bonferroni correction using all genes in the dataset.

# considering FindAllMarkers now uses the 'stage' groups (prenatal/healthy/diseased) it compares each cluster to the other 2 clusters in the following way:
# pct.1 corresponds to cells in column 'cluster'
# pct.2 corresponds to cells in ALL the other clusters (so in our case, in the other 2 clusters)
# sometimes the p-value is 0, this is simply a case of the value being too small to be represented properly in R so it gets the value 0, these are thus considered significant


stage_dimplot = DimPlot(all_integrated, group.by='stage', label=T, shuffle=T) + labs(title='Integrated cells grouped by stage', size=element_text(size=5))
source_dimplot = DimPlot(all_integrated, group.by='source', label=T, shuffle=T) + labs(title='Integrated cells grouped by source', size=element_text(size=5))
age_dimplot = DimPlot(all_integrated, group.by='age',
                      order = c('e12.5','e16.5','9','11','70','71','73','75','77','84','98','120','600'), label=T) + labs(title='Integrated cells grouped by age', size=element_text(size=5))
age_split_dimplot = DimPlot(all_integrated, group.by='age', split.by='stage',
                      order = c('e12.5','e16.5','9','11','70','71','73','75','77','84','98','120','600'), label=T) + labs(title='Integrated cells grouped by age', size=element_text(size=5)) & NoLegend()
prenatal_dimplot = DimPlot(all_integrated, group.by='stage', label=F, cols=c('gray', 'gray', 'red'), shuffle=T) + labs(title='Prenatal cells', size=element_text(size=5)) & NoLegend()
postnatal_healthy = DimPlot(all_integrated, group.by='stage', label=F, cols=c('gray', 'red', 'gray'), shuffle=T) + labs(title='Postnatal healthy cells') & NoLegend()
postnatal_diseased = DimPlot(all_integrated, group.by='stage', label=F, cols=c('red', 'gray', 'gray'), shuffle=T) + labs(title='Postnatal diseased cells') & NoLegend()
