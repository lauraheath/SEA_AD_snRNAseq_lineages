install.packages("devtools")
library(devtools)
install_github("lhe17/nebula")
install.packages('Seurat')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(nebula)
library(biomaRt)
library(Seurat)


#Calculate differential expression by ADNC status, using nebula (runs a linear mixed effect model to account for having
#multiple cells per donor)

#Upload celltype-specific Seurat object with metadata, both sexes included (astrocytes)
p1 <- synapser::synGet('syn31618644')
AllenAD <- readRDS(p1$path)
#look at metadata
head(x=AllenAD[[]])
table(AllenAD$sex)
#how many donors?
length(unique(AllenAD$`Donor ID`))
#how many genes/cells?
dim(AllenAD)

#extract the counts matrix and save the metadata file for later
counts <- GetAssayData(object = AllenAD, slot = "counts")
seaMeta <- AllenAD@meta.data
#create a column containing barcodes
seaMeta$TAG <- rownames(seaMeta)


#replace ensembl identifiers with gene short names with following functions
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

#extract ensembl gene names from counts matrix
ensembl <- row.names(counts)
ensembl <- as.data.frame(ensembl)
#run function to get short gene names
ensembl$gene_short_name <- Make.Gene.Symb(ensembl$ensembl)
gene_short_name <- ensembl$gene_short_name
#for ensembl genes that did not receive a gene short name from biomart, replace with the ensembl ID
ensembl$gene_short_name <- ifelse(ensembl$gene_short_name=="", ensembl$ensembl, ensembl$gene_short_name)

rownames(counts) <- ensembl$gene_short_name

#delete genes that did not get translated by biomaRt (virtually all are RNA genes, pseudogenes, and have low expression)
counts2 <- counts[-which(grepl("^ENSG",counts@Dimnames[[1]])),]
counts3 <- counts2[-which(grepl("^[0-9]",counts2@Dimnames[[1]])),]

dim(counts3)

#recreate seurat object, which will fix column headers with spaces in the metadata and make rownames unique (seurat will list changed headers in red)
seaAD <- CreateSeuratObject(counts = counts3, meta.data = seaMeta, min.cells = 3)
dim(seaAD)
#calculate library size factors for DE analysis and add them to the seurat object:
libsizes <- colSums(counts3)
seaAD <- AddMetaData(seaAD, metadata=libsizes, col.name='libsizes')
seaAD$size.factor <- seaAD$libsizes/mean(seaAD$libsizes)


table(seaAD$ADNC)
#remove reference donors
seaAD <- subset(seaAD, ADNC!='Reference')
#allen <- subset(allen, ADNC!='Not AD')

seaAD$ADNC <- as.character(seaAD$ADNC)
seaAD$ADNC <- as.factor(seaAD$ADNC)

#differential testing comparisons:
#test1: Not AD vs any AD
#test2: Not AD vs Low/Int AD
#test3: Low/Int AD vs High AD

#create separate objects for each comparison
test1 <- seaAD
test1$diag <- ifelse(test1$ADNC=='Not AD', 0, 1)
test2 <- subset(seaAD, ADNC!='High')
test2$diag <- ifelse(test2$ADNC=='Not AD', 0, 1)
test3 <- subset(seaAD, ADNC!='Not AD')
test3$diag <- ifelse(test3$ADNC!='High', 0, 1)

#look at numbers of cells per group, double check that 'diag' is appropriately labeled
table(test1$diag, test1$ADNC)
table(test2$diag, test2$ADNC)
table(test3$diag, test3$ADNC)
table(seaAD$ADNC)


#order the metadata by donor ID, then reorder the counts matrix to match:
meta <- test1@meta.data
mat <- GetAssayData(object=test1, slot="counts")
meta <- meta[order(meta$Donor.ID),]
col.order <- meta$TAG
mat <- mat[,col.order]
#set design matrix (adjust for sex) and run nebula (this takes a while, be patient)
modmat = model.matrix(~diag+sex, data=meta)
model1 <- nebula(mat, meta$Donor.ID, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the comparison, and save
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_diag, "fdr")
model1sum <- cbind(adjustp, model1sum)
model1sum$model <- 'NotAD_vs_All'
head(model1sum)


#model test2:
#order the metadata by subject ID, then reorder the counts matrix to match:
meta <- test2@meta.data
mat <- GetAssayData(object=test2, slot="counts")
meta <- meta[order(meta$Donor.ID),]
col.order <- meta$TAG
mat <- mat[,col.order]
#set design matrix and run nebulaa
modmat = model.matrix(~diag+sex, data=meta)
model2 <- nebula(mat, meta$Donor.ID, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the comparison, and save for now
model2sum <- model2$summary
adjustp <- p.adjust(model2sum$p_diag, "fdr")
model2sum <- cbind(adjustp, model2sum)
model2sum$model <- 'NotADvsLowInt'
head(model2sum)

#model test3
meta <- test3@meta.data
mat <- GetAssayData(object=test3, slot="counts")
meta <- meta[order(meta$Donor.ID),]
col.order <- meta$TAG
mat <- mat[,col.order]
modmat = model.matrix(~diag+sex, data=meta)
model3 <- nebula(mat, meta$Donor.ID, pred=modmat, offset=meta$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the comparison, and save for now
model3sum <- model3$summary
adjustp <- p.adjust(model3sum$p_diag, "fdr")
model3sum <- cbind(adjustp, model3sum)
model3sum$model <- 'LowInt_vs_High'
head(model3sum)

#combine the three sets of DE gene results
DEgenes <- rbind(model1sum, model2sum)
DEgenes <- rbind(DEgenes, model3sum)


write.csv(DEgenes, file='~/SEA_AD_snRNAseq_lineages/data_objects/AllenAstro_DE_genesALL.csv', row.names=FALSE)

#save in celltype-specific folder in project SEA_AD_lineage_analysis (syn35549553)
file <- synapser::File(path='~/SEA_AD_snRNAseq_lineages/data_objects/SEA_AD_CELLTYPEHERE_DEgenes.csv', parentId='FILL IN CELLTYPE FOLDER ID HERE')
file <- synapser::synStore(file)

#how many genes were significant?
pvalue05 <- DEgenes[(DEgenes$p_diag<0.05),]
length(unique(pvalue05$gene))
