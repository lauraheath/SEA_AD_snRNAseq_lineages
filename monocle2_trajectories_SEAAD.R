install.packages('lme4')

library(monocle)
library(dplyr)
library(ggplot2)
library(lme4)
library(nebula)
library(biomaRt)
library(readxl)
synapser::synLogin()


##############################################################################




#upload sample-corrected counts matrix
#counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected.RDS")
#metadata <- read.csv(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata.csv")

counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_samplecorrected.RDS")
metadata <- read.csv(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_metadata.csv")

#counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_samplecorrected_withrefs.RDS")
#metadata <- read.csv(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_metadata_withrefs.csv")

#or upload from synapse:
#female sample-corrected matrix
#countsObj1 <- synapser::synGet('syn31925053')
#counts <- readRDS(countsObj1$path)
#counts <- as.matrix(counts)
#metobj1 <- synapser::synGet('syn31925214')
#metadata <- read.csv(metobj1$path)
#dim(counts)
#dim(metadata)

#OR
#male sample-corrected matrix
#countsObj1 <- synapser::synGet('syn32110604')
#counts <- readRDS(countsObj1$path)
#counts <- as.matrix(counts)
#metobj1 <- synapser::synGet('syn31925214')
#metadata <- read.csv(metobj1$path)
#dim(counts)
#dim(metadata)



#DEgenes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AllenAstro_DE_genesALL.csv')
#or upload from synapse:
#deobj1 <- synapser::synGet('syn31925276')
DEgenes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AllenAstro_DE_genesALL.csv')
pvalue05 <- DEgenes[(DEgenes$p_diag<0.05),]
length(unique(pvalue05$gene))

#corrected sample with NOT AD patients:
counts <- readRDS(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_samplecorrected2.RDS")
metadata <- read.csv(file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_metadata2.csv")

#DEgenes <- read.csv(deobj1$path)
#DEgenes05 <- subset(DEgenes, DEgenes$p_diag<0.05)
#length(unique(DEgenes05$gene))


genes2 <- c()
#for (gene in unique(c(as.vector(DEgenes05$gene)))){
for (gene in unique(c(as.vector(pvalue05$gene)))){
  #for (gene in unique(c(as.vector(wilcox_DEgenes$genes)))){
  
  if (gene %in% rownames(counts)){
    genes2 <- c(genes2,which(rownames(counts)==gene))
  }
}

counts2 <- counts[genes2,]
dim(counts2)


gene_short_name <- rownames(counts2)
dim(counts2)
gene_short_name <- as.data.frame(gene_short_name)
rownames(gene_short_name)<-gene_short_name$gene_short_name

temp <- counts2
temp2 <- metadata

rownames(temp) <- NULL
colnames(temp) <- NULL



pd <- new("AnnotatedDataFrame", data = temp2)
fd <- new("AnnotatedDataFrame", data = gene_short_name)



MonRun <- newCellDataSet(as.matrix(temp),
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = tobit())

#visualize the cells
MonRun <- reduceDimension(MonRun, max_components=2, reduction_method = 'tSNE', norm_method = 'none', pseudo_expr = 0, verbose = TRUE)
MonRun <- clusterCells(MonRun)
plot_cell_clusters(MonRun, 1, 2, color = "Supertype")
plot_cell_clusters(MonRun, 1, 2, color = "Donor.ID")
plot_cell_clusters(MonRun, 1, 2, color = "ADNC")


MonRun$test1 <- ifelse(MonRun$Donor.ID=='H21.33.020' | MonRun$Donor.ID=='H20.33.038', 1, 0)
table(MonRun$test1, MonRun$Donor.ID)

plot_cell_clusters(MonRun, 1, 2, color = "test1")

MonRun$test1 <- ifelse(MonRun$Donor.ID=='H21.33.001','test','bgrnd')
table(MonRun$test1, MonRun$Donor.ID)
plot_cell_clusters(MonRun, 1, 2, color = "test1")

#reduce dimensions again with DDRTree for the trajectory inference
MonRun <- reduceDimension(MonRun, max_components=2, reduction_method = 'DDRTree', norm_method = 'none', pseudo_expr = 0, verbose = TRUE)

#or upload this monocle object with reduced dimensions from synapse:
#monobj <- synapser::synGet('syn31925811')
#MonRun <- readRDS(monobj$path)

#or upload this monocle object with reduced dimensions from synapse:
monobj <- synapser::synGet('syn32151522')
MonRun <- readRDS(monobj$path)

#save the Monocle object with reduced dimensions:
#saveRDS(MonRun, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_reduceddim.RDS')
#file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_reduceddim.RDS', parentId='syn31924906')
#file <- synapser::synStore(file)


MonRun <- orderCells(MonRun)
#saveRDS(MonRun, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_ordered.RDS')

#with Not AD patients:
saveRDS(MonRun, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_ordered2.RDS')
#file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_ordered.RDS', parentId='syn31924906')
#file <- synapser::synStore(file)

# MonRun <- orderCells(MonRun)
saveRDS(MonRun, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_MonRun_ordered2.RDS')
# MonRun <- readRDS(file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_MonRun_ordered_withrefs.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_MonRun_ordered.RDS', parentId='syn31924906')
file <- synapser::synStore(file)


MonRun <- readRDS(file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_MonRun_ordered.RDS')
#MonRunObj <- synapser::synGet('syn31957036')
#MonRun <- readRDS(MonRunObj$path)

####### For non-batch-corrected matrix, run RunMonocleTobit. For batch-corrected matrix, run RunMonocleTobit2 ###

#MonRun <- RunMonocleTobit2(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroF_samplecorrected_tree.tiff',height=125,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroM_samplecorrected_tree.tiff',height=125,width=200,units='mm',res=300)
#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroM_samplecorrectedwithREFS_tree.tiff',height=125,width=200,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "ADNC",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis") + guides(colour=guide_legend(override.aes=list(size=3)))
g
dev.off()

#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroM_state_withREFS_tree.tiff',height=125,width=200,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State") + guides(colour=guide_legend(override.aes=list(size=3)))
g
#dev.off()

g<- plot_cell_trajectory(MonRun,color_by = "Supertype",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Supertype")
g


plot_cell_trajectory(MonRun,color_by = "Donor.ID",show_branch_points=F,use_color_gradient = F,cell_size = 1)
plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)

#determine which state has highest proportion of control cells compared to AD cells
df <- list()
df$State <- MonRun$State
df$ADNC <- MonRun$ADNC
df <- as.data.frame(df)
table(df$State)
proptable <- with(df, table(ADNC, State)) %>% prop.table(margin = 2)
proptable

table(MonRun$State)

#If necessary, reset root state to the state with the highest proportion of control cells

#Astro F: set state 7 as root state
MonRun <- orderCells(MonRun, root_state = 7)
#Astro M: set state 7 as root state (homeostatic state)
#MonRun <- orderCells(MonRun, root_state = 7)


plot_cell_trajectory(MonRun,color_by = "Pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1)

plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)


#relabel states:
#AstroF 
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 7] <- 1
MonRun$State2[MonRun$State == 4] <- 2
MonRun$State2[MonRun$State == 5] <- 3
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 3] <- 5
MonRun$State2[MonRun$State == 2] <- 6
MonRun$State2[MonRun$State == 1] <- 7

#AstroM
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 7] <- 1
MonRun$State2[MonRun$State == 3] <- 2
MonRun$State2[MonRun$State == 5] <- 3
MonRun$State2[MonRun$State == 2] <- 4
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 9] <- 6
MonRun$State2[MonRun$State == 1] <- 7
MonRun$State2[MonRun$State == 6] <- 8
 
#AstroF with NOT AD patients
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 7] <- 1
MonRun$State2[MonRun$State == 1] <- 2
MonRun$State2[MonRun$State == 2] <- 3
MonRun$State2[MonRun$State == 3] <- 4
MonRun$State2[MonRun$State == 4] <- 5
MonRun$State2[MonRun$State == 5] <- 5
MonRun$State2[MonRun$State == 6] <- 6

#AstroM with NOT AD patients
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 1] <- 1
MonRun$State2[MonRun$State == 11] <- 2
MonRun$State2[MonRun$State == 2] <- 3
MonRun$State2[MonRun$State == 10] <- 3
MonRun$State2[MonRun$State == 3] <- 3
MonRun$State2[MonRun$State == 9] <- 4
MonRun$State2[MonRun$State == 4] <- 5
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 5] <- 6
MonRun$State2[MonRun$State == 6] <- 7
MonRun$State2[MonRun$State == 7] <- 8



MonRun$State2 <- as.character(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)
plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1)
df <- list()
df$State <- MonRun$State2
df$ADNC <- MonRun$ADNC
df <- as.data.frame(df)
table(df$State)
proptable <- with(df, table(ADNC, State)) %>% prop.table(margin = 2)
proptable
plot_cell_trajectory(MonRun,color_by = "Supertype",show_branch_points=F,use_color_gradient = F,cell_size = 1)

#save the MonRun object with saved reassigned States and ordered cells:
saveRDS(MonRun, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_ordered_notADincl.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_MonRun_ordered_notADincl.RDS', parentId='syn31924906')
file <- synapser::synStore(file)

saveRDS(MonRun, file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_MonRun_ordered.RDS')
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_MonRun_ordered.RDS', parentId='syn31924906')
file <- synapser::synStore(file)




#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroF_state_tree.tiff',height=125,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroM_state_tree.tiff',height=125,width=200,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State") + guides(colour=guide_legend(override.aes=list(size=3)))
g
dev.off()


#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroF_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroM_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=Braak.stage, y=scale(Pseudotime,center=F),fill=Braak.stage)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()

#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenastroF_bargraph_Thal.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AllenastroM_bargraph_Thal.tiff',height=85,width=100,units='mm',res=300)

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=Thal.phase, y=scale(Pseudotime,center=F),fill=Thal.phase)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Thal phase",y="Pseudotime",x="Thal phase")
g
dev.off()


#change levels of CERAD variables
MonRun$CERAD.score <- factor(MonRun$CERAD.score, levels = c("Absent", "Sparse", "Moderate", "Frequent"))
#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroF_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroM_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=CERAD.score, y=scale(Pseudotime,center=F),fill=CERAD.score)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD.score")
g
dev.off()

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=Lewy.body.disease.pathology, y=scale(Pseudotime,center=F),fill=Lewy.body.disease.pathology)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Lewy body disease path",y="Pseudotime",x="Lewy body disease path")
g

head(pData(MonRun))
x <- list()
x$cellID <- MonRun$TAG
x$Donor.ID <- MonRun$Donor.ID
x$Pseudotime <- MonRun$Pseudotime
x$State <- MonRun$State2
x$Supertype <- MonRun$Supertype
x$Years.of.education <- MonRun$Years.of.education
x$Age.at.death <- MonRun$Age.at.death
x$sex <- MonRun$sex
x$ethnicity <- MonRun$ethnicity
x$Cognitive.status <- MonRun$Cognitive.status
x$disease <- MonRun$disease
x$ADNC <- MonRun$ADNC
x$Braak.stage <- MonRun$Braak.stage
x$Thal.phase <- MonRun$Thal.phase
x$CERAD.score <- MonRun$CERAD.score
x$APOE4.statu <- MonRun$APOE4.status
x$Lewybodypath <- MonRun$Lewy.body.disease.pathology
x$Late.NC.stage <- MonRun$LATE.NC.stage
x$Specimen.ID <- MonRun$Specimen.ID
x$PMI <- MonRun$PMI
x$size.factor <- MonRun$size.factor


#females: rename and create a scaled pseudotime variable
Fvariables <- as.data.frame(x)
Fvariables$pseudotime_sc <- scale(Fvariables$Pseudotime, center=F)

#save variables file for later
#write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_allgenes_pstimeStates.csv", row.names=FALSE)
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_allgenes_pstimeStates.csv', parentId='syn31924906')
file <- synapser::synStore(file)

write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_allgenes_pstimeStates.csv", row.names=FALSE)
file <- synapser::File(path='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_allgenes_pstimeStates.csv', parentId='syn31924906')
file <- synapser::synStore(file)

write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_withReferences_pstimeStates.csv", row.names=FALSE)
write.csv(Fvariables, file="~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_withNotAD_pstimeStates.csv", row.names=FALSE)

#run logistic regression comparing pseudotiem between cases and controls only
#casecontrolF <- subset(Fvariables, Fvariables$simpleDiagnosis=='AD'|Fvariables$simpleDiagnosis=='Cont')

Fvariables$diag1 <- ifelse(Fvariables$disease=='dementia', 1, 0)


summary(glm(diag1 ~ pseudotime_sc,Fvariables,family='binomial'))

#Not AD vs low
NvsL <- subset(Fvariables, Fvariables$ADNC=='Not AD' | Fvariables$ADNC=='Low')
NvsL$diag <- ifelse(NvsL$ADNC=='Not AD', 0, 1)
summary(glm(diag ~ pseudotime_sc,NvsL,family='binomial'))
#Not AD vs intermediate
NvsI <- subset(Fvariables, Fvariables$ADNC=='Not AD' | Fvariables$ADNC=='Intermediate')
NvsI$diag <- ifelse(NvsI$ADNC=='Not AD', 0, 1)
summary(glm(diag ~ pseudotime_sc,NvsI,family='binomial'))
#Not AD vs High
NvsH <- subset(Fvariables, Fvariables$ADNC=='Not AD' | Fvariables$ADNC=='High')
NvsH$diag <- ifelse(NvsH$ADNC=='Not AD', 0, 1)
summary(glm(diag ~ pseudotime_sc,NvsH,family='binomial'))

#low vs intermediate
LvsI <- subset(Fvariables, Fvariables$ADNC=='Low' | Fvariables$ADNC=='Intermediate')
LvsI$diag <- ifelse(LvsI$ADNC=='Low', 0, 1)
summary(glm(diag ~ pseudotime_sc,LvsI,family='binomial'))

#Low vs High
LvsH <- subset(Fvariables, Fvariables$ADNC=='Low' | Fvariables$ADNC=='High')
LvsH$diag <- ifelse(LvsH$ADNC=='Low', 0, 1)
summary(glm(diag ~ pseudotime_sc,LvsH,family='binomial'))

#Intermediate vs High
IvsH <- subset(Fvariables, Fvariables$ADNC=='Intermediate' | Fvariables$ADNC=='High')
IvsH$diag <- ifelse(IvsH$ADNC=='Intermediate', 0, 1)
summary(glm(diag ~ pseudotime_sc,IvsH,family='binomial'))

#summary(lmer(diag3 ~ pseudotime_sc + pmi + (pseudotime_sc|ros_ids), Fvariables))
#summary(lmer(diag3 ~ pseudotime_sc + (1|ros_ids), Fvariables))

theme_set(theme_gray())
Fvariables$ADNC <- ordered(Fvariables$ADNC, levels = c("Not AD", "Low", "Intermediate", "High"))
#tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroF_bargraph_ADNC.tiff',height=170,width=200,units='mm',res=300)
tiff(file='~/scRNAseq-subtype-mapping/figures/AllenAstroM_bargraph_ADNC.tiff.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(Fvariables,aes(x=ADNC,
                             y=pseudotime_sc,
                             color=ADNC)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(4)[1:4]) 
g
dev.off()




#Fvariables$Diagnosis2 <- ifelse(Fvariables$Diagnosis=='Cont', 0,
#                                 ifelse(Fvariables$Diagnosis=='Early', 1, 2))
# Fvariables$Diagnosis2 <- factor(Fvariables$Diagnosis2,levels = c(0:2))
# Fvariables$braaksc <- factor(Fvariables$braaksc,levels = c(1:6))
# Fvariables$ceradsc <- factor(Fvariables$ceradsc,levels = c(1:4))
# Fvariables$cogdx <- factor(Fvariables$cogdx,levels = c(1:6))
# 
# 
# #run proportional odds logistic regression for neuropath/cognitive endpoints:
# diagfit <- MASS::polr(Diagnosis2 ~ pseudotime_sc,Fvariables)
# braakfit <- MASS::polr(braaksc ~ pseudotime_sc,Fvariables)
# ceradfit <- MASS::polr(ceradsc ~ pseudotime_sc,Fvariables)
# cogdxfit <- MASS::polr(cogdx ~ pseudotime_sc,Fvariables)
# 
# cat('diagnosis p-value: ',pt(abs(summary(diagfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
# cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
# cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
# cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')
# 
# 
# 
# #BEAM (branch point genes)
# plot_cell_trajectory(MonRun, color_by = "State")
# BEAM_res <- BEAM(MonRun, branch_point = , cores = 1)
# BEAM_res <- BEAM_res[order(BEAM_res$qval),]
# BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
# 
# #by patient


table(MonRun$Donor.ID)

MonRun2 <- MonRun


MonRun2$Donor.ID2 <- ifelse(MonRun2$Donor.ID=='H21.33.001',1,0)
MonRun2$Donor.ID2 <- as.factor(MonRun2$Donor.ID2)
plot_cell_trajectory(MonRun2,color_by = "Donor.ID2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)


MonRun2$Donor.ID2 <- ifelse(MonRun2$Donor.ID=='H20.33.033',1,0)
MonRun2$Donor.ID2 <- as.factor(MonRun2$Donor.ID2)
plot_cell_trajectory(MonRun2,color_by = "Donor.ID2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)


table(MonRun2$PMI)
MonRun2$PMI2 <- ifelse(MonRun2$PMI=='8.7 to 11.4 hours',1,0)
MonRun2$PMI2 <- as.factor(MonRun2$PMI2)
plot_cell_trajectory(MonRun2,color_by = "PMI2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# 
# MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS14', 'ROS14', 'OTHER')
# plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# 
# MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS15', 'ROS15', 'OTHER')
# plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# 
# MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS16', 'ROS16', 'OTHER')
# plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# 
# MonRun2$ros_ids2 <- ifelse(MonRun2$ros_ids=='ROS17', 'ROS17', 'OTHER')
# plot_cell_trajectory(MonRun2,color_by = "ros_ids2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)



##### find genes that change expression as a function of pseudotime #####


#female seurat object
p <- synapser::synGet('syn31987861')
astro <- readRDS(p$path)

head(x=astro[[]])
counts <- GetAssayData(object=astro, slot="counts")
counts <- as.matrix(counts)
dim(counts)


#upload dataframe with pseudotimes and states from monocle2_micro_trajectories.R
pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroF_allgenes_pstimeStates.csv')

#order the metadata by subject ID, then reorder the counts matrix to match:
rownames(pstimes) <- pstimes$cellID
pstimes <- pstimes[order(pstimes$Donor.ID),]
col.order <- pstimes$cellID
counts <- counts[,col.order]


###### Find genes that change as a function of pseudotime
#perform nebula on pseudotime ~ expression
subjectID <- pstimes$Donor.ID
modmat = model.matrix(~Pseudotime, data=pstimes)
model1 <- nebula(counts, subjectID, pred=modmat, offset=pstimes$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_Pseudotime, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)


write.csv(model1sum, file="~/scRNAseq-subtype-mapping/DEresults/AllenAstroF_genes_changing_by_pstime.csv", row.names=FALSE)






#male seurat object
p2 <- synapser::synGet('syn32109687')
astro <- readRDS(p2$path)

head(x=astro[[]])
counts <- GetAssayData(object=astro, slot="counts")
counts <- as.matrix(counts)
dim(counts)

#upload dataframe with pseudotimes and states from monocle2_micro_trajectories.R
pstimes <- read.csv(file='~/scRNAseq-subtype-mapping/data_objects/AllenAstroM_allgenes_pstimeStates.csv')

#order the metadata by subject ID, then reorder the counts matrix to match:
rownames(pstimes) <- pstimes$cellID
pstimes <- pstimes[order(pstimes$Donor.ID),]
col.order <- pstimes$cellID
counts <- counts[,col.order]


###### Find genes that change as a function of pseudotime
#perform nebula on pseudotime ~ expression
subjectID <- pstimes$Donor.ID
modmat = model.matrix(~Pseudotime, data=pstimes)
model1 <- nebula(counts, subjectID, pred=modmat, offset=pstimes$size.factor)
#get dataframe of results, perform FDR to adjust for multiple comparisons, label the state, and save for now
model1sum <- model1$summary
adjustp <- p.adjust(model1sum$p_Pseudotime, "fdr")
model1sum <- cbind(adjustp, model1sum)
head(model1sum)

write.csv(model1sum, file="~/scRNAseq-subtype-mapping/DEresults/AllenAstroM_genes_changing_by_pstime.csv", row.names=FALSE)
