remotes::install_github("HenrikBengtsson/matrixStats", ref="develop")
BiocManager::install("BiocSingular")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scran")
BiocManager::install("batchelor")

library(Seurat)
library(SingleCellExperiment)
library(scran)
library(batchelor)
library(biomaRt)

setwd("~/scRNAseq-subtype-mapping/")

#NOT sample-corrected data for Differential expression between ADNC groups
p1 <- synapser::synGet('syn31618644')
AllenAD <- readRDS(p1$path)
head(x=AllenAD[[]])
table(AllenAD$sex)
length(unique(AllenAD$`Donor ID`))
AllenAD

counts <- GetAssayData(object = AllenAD, slot = "counts")
metadata <- AllenAD@meta.data
metadata$TAG <- rownames(metadata)



###########################################################################################
# 

#replace ensembl identifiers with gene short names
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='https://www.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}

Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}

ensembl <- row.names(counts)
ensembl <- as.data.frame(ensembl)
ensembl$gene_short_name <- Make.Gene.Symb(ensembl$ensembl)
gene_short_name <- ensembl$gene_short_name
ensembl$gene_short_name2 <- ifelse(ensembl$gene_short_name=="", ensembl$ensembl, ensembl$gene_short_name)

gene_short_name2 <- ensembl$gene_short_name2
rownames(counts) <- ensembl$gene_short_name2

#delete genes that did not get translated by biomaRt (virtually all are RNA genes, pseudogenes, and have low expression)
counts <- counts[-which(grepl("^ENSG",counts@Dimnames[[1]])),]
counts <- counts[-which(grepl("^[0-9]",counts@Dimnames[[1]])),]

dim(counts)

#remake seurat object, which will fix column headers with spaces in the metadata and make rownames all unique
allen <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 3)
dim(allen)
#calculate size factors for DE analysis later:
libsizes <- colSums(counts)
allen <- AddMetaData(allen, metadata=libsizes, col.name='libsizes')
allen$size.factor <- allen$libsizes/mean(allen$libsizes)

metadata <- allen@meta.data
meta_bydonor <- subset(metadata, select = c(Donor.ID, Age.at.death, sex, ADNC, Cognitive.status, Braak.stage, Thal.phase, CERAD.score, APOE4.status, Lewy.body.disease.pathology, LATE.NC.stage,Microinfarct.pathology, PMI, ethnicity))
meta_bydonor <- meta_bydonor[!duplicated(meta_bydonor),]
table(meta_bydonor$ADNC)


#need to pare the dataset down (monocle has memory issues with cells over ~25000 or 30000, even with very large EC2 instance sizes)
#data is overrrepresented by individuals in the High ADNC group (Not AD=9, Low=12, Int=21, High=42)
#after excluding individuals with high microinfarct path and high pmi:
#for female samples:
#meta_bydonor2 <- subset(meta_bydonor, meta_bydonor$sex=='female')
meta_bydonor2 <- subset(meta_bydonor, meta_bydonor$sex=='male')
table(meta_bydonor2$ADNC)
meta_bydonor2 <- subset(meta_bydonor2, meta_bydonor2$ADNC!='Reference')
#meta_bydonor2 <- subset(meta_bydonor2, meta_bydonor2$Microinfarct.pathology=='0 to 3 microinfarcts')
#meta_bydonor2 <- subset(meta_bydonor2, meta_bydonor2$PMI!='8.7 to 11.4 hours')
table(meta_bydonor2$ADNC)


# table(MonRun$Donor.ID)
# table(MonRun$ADNC)
# DonorIDs <- unique(MonRun$Donor.ID)
# metadata_monrun <- metadata[metadata$Donor.ID %in% DonorIDs,]
# metadata_monrun$include <- 'Yes'
# metada_monrun <- subset(metadata_monrun, select=c(Donor.ID, include))
# metadata_all <- merge(metadata, metadata_monrun, all=TRUE)

 #randomly select 10 Highs

#Intermed <- subset(meta_bydonor2, meta_bydonor2$ADNC=='Intermediate')
High <- subset(meta_bydonor2, meta_bydonor2$ADNC=='High')
#sample1 <- sample(nrow(Intermed), 10)
#InterSample <- Intermed[sample1,]
sample2 <- sample(nrow(High), 10)
HighSample <- High[sample2,]

#InterSample$Donor.ID<-as.character(InterSample$Donor.ID)
#InterSample$Donor.ID
HighSample$Donor.ID<-as.character(HighSample$Donor.ID)
HighSample$Donor.ID

#collect IDs 
NoAD <- subset(meta_bydonor2, meta_bydonor2$ADNC=='Not AD')
Low <- subset(meta_bydonor2, meta_bydonor2$ADNC=='Low')
InterSample <- subset(meta_bydonor2, meta_bydonor2$ADNC=='Intermediate')
#final data set
metadata_small <- rbind(NoAD, Low, InterSample, HighSample)

metadata_small$include <- "Yes"
metada_small <- subset(metadata_small, select=c(Donor.ID, include))
metadata_all <- merge(metadata, metadata_small, all=TRUE)

#reorder metadata rows (same as 'TAG' column variable) to match columns of counts matrix:
colnames_counts <- as.data.frame(colnames(allen))
rownames(metadata_all)<-metadata_all$TAG
names(colnames_counts)[names(colnames_counts) == "colnames(allen)"] <- "columnnames"
metadata_all <- metadata_all[order(match(metadata_all$TAG, colnames_counts$columnnames)),]
#add inclusion variable to metadata of seurat object
inclusions <- metadata_all$include
allen <- AddMetaData(allen, metadata=inclusions, col.name='include')

#table(meta_bydonor2$ADNC)

#For all samples: remove ADNC Reference group and non-AD individuals; 
#For female dataset, also remove individuals with high PMI (8.7-11.4 hours), and individuals with high microinfarcts (> 0-3 microinfarcts) to shrink data size for monocle2
#allen <- subset(allen, ADNC!='Reference')
#allen <- subset(allen, ADNC!='Not AD')

#allen <- subset(allen, Microinfarct.pathology=='0 to 3 microinfarcts')
#allen <- subset(allen, PMI!='8.7 to 11.4 hours')
allen <- subset(allen, include=='Yes')


table(allen$ADNC)
dim(allen)

allen$ADNC<-as.character(allen$ADNC)
allen$ADNC<-as.factor(allen$ADNC)
length(unique(allen$Donor.ID))
table(allen$Donor.ID)
dim(allen)
table(allen$ADNC)

saveRDS(allen, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_Seurat_notADincl.RDS") 
saveRDS(allen, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_Seurat_notADincl.RDS") 
#save to synapse:
#file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_Seurat.RDS', parentId='syn31924906')
#file <- synapser::synStore(file)


#batch correct by sex to shrink sample size:
#AstroFemale <- subset(allen, sex=='female')
#AstroFemale <- subset(allen, sex=='male')
#length(unique(AstroFemale$Donor.ID))
#dim(AstroFemale)

counts <- GetAssayData(object = allen, slot = "counts")
metadata <- allen@meta.data

AstroF <- CreateSeuratObject(counts = counts, project = "AllenSCRNAseq", min.cells=3)
AstroF <- AddMetaData(AstroF, metadata)
head(x=AstroF[[]])
dim(AstroF)
AstroF$Donor.ID<-as.character(AstroF$Donor.ID)
AstroF$Donor.ID<-as.factor(AstroF$Donor.ID)

length(unique(AstroF$Donor.ID))

saveRDS(AstroF, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_Seurat.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_Seurat.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_Seurat_notADincl.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

saveRDS(AstroF, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_Seurat_notADincl.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_Seurat_notADincl.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

#AstroF <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_Seurat_notADincl.RDS")
astroFobj <- synapser::synGet('syn31987861')
AstroF <- readRDS(astroFobj$path)

AstroF[["percent.mt"]] <- PercentageFeatureSet(AstroF, pattern = "^MT-")
VlnPlot(AstroF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(AstroF, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AstroF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


length(unique(AstroF$Donor.ID))
#AstroM <- subset(AllenAD, sex=='male')


#To visualize the cells, run Seurat pipeline to normalize data and perform dimension reduction
AstroF <- SCTransform(AstroF, do.scale = TRUE, verbose = TRUE, conserve.memory=TRUE)

# AstroF <- NormalizeData(AstroF, normalization.method = "LogNormalize", scale.factor = 10000)
# AstroF <- FindVariableFeatures(AstroF, selection.method = "vst", nfeatures = 2000)
# 
# top10 <- head(VariableFeatures(AstroF), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(AstroF)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# all.genes <- rownames(AstroF)
# AstroF <- ScaleData(AstroF, features = all.genes)




AstroF <- RunPCA(AstroF, npcs = 30)
AstroF <- RunUMAP(AstroF, dims = 1:30)
#for cleaner clusters, decrease resolution in FindNeighbors
AstroF <- FindNeighbors(AstroF, reduction = "umap", dims = 1:2)
AstroF <- FindClusters(AstroF, resolution = 0.5, verbose = TRUE)
DimPlot(AstroF, group.by = "Supertype", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(AstroF, reduction="umap", label = TRUE, repel=TRUE)
DimPlot(AstroF, group.by = "ADNC", reduction="umap", label = TRUE, repel=TRUE)
DimPlot(AstroF, group.by = "Donor.ID", reduction="umap", label = TRUE, repel=TRUE)



#For batch correction: extract counts and metadata from seurat object, and make a singlecellexperiment object
counts <- GetAssayData(object = AstroF, slot = "counts")
#counts <- GetAssayData(object = glialcellsM, slot = "counts")
counts <- as.matrix(counts)
metadata <- AstroF@meta.data

#create single single experiment object
sce <- SingleCellExperiment(list(counts=counts), 
                            colData=metadata)
sce
head(assay(sce[0:20,0:20]))
head(colData(sce))

clusters <- quickCluster(sce, min.size=100)
sce <- computeSumFactors(sce, cluster=clusters)
#sce <- computeSumFactors(sce, cluster=clusters, BPPARAM=MulticoreParam(8))
summary(sizeFactors(sce))

sce <- logNormCounts(sce)
head(logcounts(sce[0:20,0:20]))
head(rownames(logcounts(sce)))

dim(logcounts(sce))


sce$Donor.ID <- as.character(sce$Donor.ID)
sce$Donor.ID <- as.factor(sce$Donor.ID)

#sce4 <- fastMNN(sce, batch = sce$ros_ids, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
sce4 <- fastMNN(sce, batch = sce$Donor.ID, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)

#sce4 <- fastMNN(sce, batch = sce$ros_ids, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce5 <- fastMNN(sce, batch = sce$batch, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)
#sce4 <- fastMNN(sce, batch = sce$batch, subset.row=chosen.hvgs, correct.all=TRUE, cos.norm=TRUE, auto.merge=TRUE, prop.k=0.05)


corrected_sample <- assay(sce4)
#corrected_batch <- assay(sce5)


saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected.RDS")
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

#save female metadata
metadata <- AstroF@meta.data
metadata$TAG <- rownames(metadata)
write.csv(metadata, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata.csv", row.names = FALSE)
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata.csv', parentId='syn31924906')
file <- synapser::synStore(file)




#save randomly selected subset of patients with 'Not_AD' subjects
saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected2.RDS")
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected2.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

#save female metadata
metadata <- AstroF@meta.data
metadata$TAG <- rownames(metadata)
write.csv(metadata, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata2.csv", row.names = FALSE)
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata2.csv', parentId='syn31924906')
file <- synapser::synStore(file)




#save male sample-corrected matrix
saveRDS(corrected_sample, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_samplecorrected.RDS")
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_samplecorrected.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

#save male metadata
metadata <- AstroF@meta.data
metadata$TAG <- rownames(metadata)
write.csv(metadata, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_metadata.csv", row.names = FALSE)
#save to synapse
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_metadata.csv', parentId='syn31924906')
file <- synapser::synStore(file)

