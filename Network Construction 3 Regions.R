#Read in the three regions s score spreadsheets
NACdata <- read.csv(file="BXD_NA_Sscores_All_Data_mm9.csv", header=TRUE, sep=",") 
PFCdata <- read.csv(file="BXD_PFC_Sscores_All_Data_mm9.csv", header=TRUE, sep=",") 
VTAdata <- read.csv(file="BXD_VTA_Sscores_All_Data_mm9.csv", header=TRUE, sep=",")
#filtering data by removing empty spaces and then those of a q value above .2
NACdata_F <- NACdata[complete.cases(NACdata$emp.qvalue),]
PFCdata_F <- PFCdata[complete.cases(PFCdata$emp.qvalue),]
VTAdata_F <- VTAdata[complete.cases(VTAdata$emp.qvalue),] 
NACdata_F <- NACdata_F[NACdata_F$emp.qvalue<=.2,]
PFCdata_F <- PFCdata_F[PFCdata_F$emp.qvalue<=.2,]
VTAdata_F <- VTAdata_F[VTAdata_F$emp.qvalue<=.2,]  
NAC_Genelist <- data.frame(probset_id=NACdata_F$Probeset,entrez_id=NACdata_F$Entrez) 
PFC_Genelist <- data.frame(probset_id=PFCdata_F$Probeset,entrez_id=PFCdata_F$Entrez) 
VTA_Genelist <- data.frame(probset_id=VTAdata_F$Probeset,entrez_id=VTAdata_F$Entrez)
Genelist <- rbind(NAC_Genelist,PFC_Genelist,VTA_Genelist) 
Genelist <- unique(Genelist) 
#-------------------------------------------------------------------------------------
# Pulling the data for relevant genes
NACdata_Genelist<-NACdata[NACdata$Probeset %in% Genelist$probset_id,] 
PFCdata_Genelist<-PFCdata[PFCdata$Probeset %in% Genelist$probset_id,] 
VTAdata_Genelist<-VTAdata[VTAdata$Probeset %in% Genelist$probset_id,] 
save(Genelist, 
     NACdata_Genelist, PFCdata_Genelist, VTAdata_Genelist, 	file="Genelist_filtered_3regions.Rdata")
#-------------------------------------------------------------------------------------
#Opening WGCNA and writing data to excel spreadsheets for easy access
library(WGCNA) 
options(stringsAsFactors = FALSE) 
disableWGCNAThreads() 
lnames = load(file = "Genelist_filtered_3regions.Rdata")
write.csv(PFCdata_Genelist, "PFCdata_Genelist.csv")
write.csv(NACdata_Genelist, "NACdata_Genelist.csv")
write.csv(NACdata_Genelist, "VTAdata_Genelist.csv")
#-------------------------------------------------------------------------------------
# Data frames as listed below for relevant information
datExpr0 = as.data.frame(t(PFCdata_Genelist[, -c(1:29)]))
datExpr1 = as.data.frame(t(NACdata_Genelist[, -c(1:21)]))
datExpr2 = as.data.frame(t(VTAdata_Genelist[, -c(1:21)]))
names(datExpr0) = PFCdata_Genelist$Probeset
rownames(datExpr0) = names(PFCdata_Genelist)[-c(1:29)]
datExpr_0<-datExpr0[-(1:2),]
names(datExpr1) = NACdata_Genelist$Probeset
rownames(datExpr1) = names(NACdata_Genelist)[-c(1:21)]
datExpr_1<-datExpr1[-(1:2),]
names(datExpr2) = VTAdata_Genelist$Probeset
rownames(datExpr2) = names(VTAdata_Genelist)[-c(1:21)]
datExpr_2<-datExpr2[-(1:2),]
#-----------------------------------------------------------------------------------------
sampleTree0=flashClust::flashClust(dist(datExpr_0), method = "average")
par(cex = 0.6); 
par(mar = c(0,4,2,0)) 
plot(sampleTree0, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr_0, powerVector = powers, networkType = "unsigned", verbose = 0)
par(mfrow = c(1,2)) 
cex1 = 0.9; 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red") 
abline(h=0.87,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
#------------------------------------------------------------------------------------
# For NAC
sampleTree1 = flashClust::flashClust(dist(datExpr_1), method = "average"); 
par(cex = 0.6); 
par(mar = c(0,4,2,0)) 
plot(sampleTree1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
par(cex = 0.6); 
par(mar = c(0,4,2,0)) 
powers1 = c(c(1:10), seq(from = 12, to=20, by=2))
sft_1 = pickSoftThreshold(datExpr_1, powerVector = powers1, networkType = "unsigned", verbose = 0)
par(mfrow = c(1,2)) 
cex1 = 0.9; 
plot(sft_1$fitIndices[,1], -sign(sft_1$fitIndices[,3])*sft_1$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"))
text(sft_1$fitIndices[,1], -sign(sft_1$fitIndices[,3])*sft_1$fitIndices[,2], labels=powers1,cex=cex1,col="red") 
abline(h=0.84,col="red")
plot(sft_1$fitIndices[,1], sft_1$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity")) 
text(sft_1$fitIndices[,1], sft_1$fitIndices[,5], labels=powers1, cex=cex1, col="red")
#----------------------------------------------------------------------------------------
#For VTA
sampleTree2 = flashClust::flashClust(dist(datExpr_2), method = "average"); 
par(cex = 0.6); 
par(mar = c(0,4,2,0)) 
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
par(cex = 0.6); 
par(mar = c(0,4,2,0)) 
powers2 = c(c(1:10), seq(from = 12, to=20, by=2))
sft_2 = pickSoftThreshold(datExpr_2, powerVector = powers2, networkType = "unsigned", verbose = 0)
par(mfrow = c(1,2)) 
cex1 = 0.9; 
plot(sft_2$fitIndices[,1], -sign(sft_2$fitIndices[,3])*sft_2$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"))
text(sft_2$fitIndices[,1], -sign(sft_2$fitIndices[,3])*sft_2$fitIndices[,2], labels=powers2,cex=cex1,col="red") 
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers2, cex=cex1, col="red")





#----------------------------------------------------------------------------------------
#Network Construction
net0 = blockwiseModules(datExpr_0, maxBlockSize= 10000, power = 12, TOMType = "unsigned", minModuleSize = 20, deepSplit = 3, reassignThreshold = 0, detectcutheight=.995, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "PFCanalysisTOM", verbose = 0) 
mergedColors0 = labels2colors(net0$colors) 
net1 = blockwiseModules(datExpr_1, maxBlockSize= 10000, power = 9, TOMType = "unsigned", minModuleSize = 20, deepSplit = 3, reassignThreshold = 0, detectcutheight=.995, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "NACanalysisTOM", verbose = 0) 
mergedColors1= labels2colors(net1$colors)
net2 = blockwiseModules(datExpr_0, maxBlockSize= 10000, power = 9, TOMType = "unsigned", minModuleSize = 20, deepSplit = 3, reassignThreshold = 0, detectcutheight=.995, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "VTAanalysisTOM", verbose = 0) 
mergedColors2= labels2colors(net2$colors)
#----------------------------------------------------------------------------------------
#Plotting Dendrograms
plotDendroAndColors(net0$dendrograms[[1]], mergedColors0[net0$blockGenes[[1]]], "Deep Split 3 PFC", dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)
plotDendroAndColors(net1$dendrograms[[1]], mergedColors1[net1$blockGenes[[1]]], "Deep Split 3 NAC", dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)
plotDendroAndColors(net2$dendrograms[[1]], mergedColors2[net2$blockGenes[[1]]], "Deep Split 3 VTA", dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)
table(mergedColors0)
table(mergedColors1)
table(mergedColors2)
# Saving module colors and labels from each dendrogram. PFC
#----------------------------------------------------------------------------------------
moduleLabels0 = net0$colors 
moduleColors0 = labels2colors(net0$colors) 
MEs0 = net0$MEs; 
geneTree0 = net0$dendrograms[[1]]; 
save(datExpr_0, MEs0, moduleLabels0, moduleColors0,
     geneTree0, file = "PFCdata_autonetwork.RData")
#--------------------------------------------------------------------------------------
# Saving module names and colors for NAC
moduleLabels1 = net1$colors 
moduleColors1 = labels2colors(net1$colors) 
MEs1 = net1$MEs; 
geneTree1 = net1$dendrograms[[1]]; 
save(datExpr_1, MEs1, moduleLabels1, moduleColors1,
     geneTree1, file = "NACdata_autonetwork.RData")
#-------------------------------------------------------------------------------------
# For VTA
moduleLabels2 = net2$colors 
moduleColors2 = labels2colors(net2$colors) 
MEs2 = net2$MEs; 
geneTree2 = net2$dendrograms[[1]]; 
save(datExpr_2, MEs2, moduleLabels2, moduleColors2,
     geneTree2, file = "VTAdata_autonetwork.RData")
#-------------------------------------------------------------------------------------
# Loading region data and mouse annotation file, labeling probes PFC
lnames_0 <- load(file = "PFCdata_autonetwork.RData")
annot<-read.csv(file = "Mouse430_2.na35.annot.csv") 
probes_0 = names(datExpr_0) 
probes2annot_0 = match(probes_0, annot$Probe.Set.ID) 
allLLIDs_0 <- annot$Entrez.Gene[probes2annot_0]
#---------------------------------------------------------------------------------------
#Loading names saved in last NAC session and mouse file
lnames_1 <- load(file = "NACdata_autonetwork.RData")
annot<-read.csv(file = "Mouse430_2.na35.annot.csv") 
probes_1 = names(datExpr_1) 
probes2annot_1 = match(probes_1, annot$Probe.Set.ID) 
allLLIDs_1 <- annot$Entrez.Gene[probes2annot_1]
#---------------------------------------------------------------------------------------
# VTA
lnames_2 <- load(file = "VTAdata_autonetwork.RData")
annot<-read.csv(file = "Mouse430_2.na35.annot.csv") 
probes_2 = names(datExpr_2) 
probes2annot_2 = match(probes_2, annot$Probe.Set.ID) 
allLLIDs_2 <- annot$Entrez.Gene[probes2annot_2]
#---------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# sort modules for external softwares PFC
intModules0 <- unique(moduleColors0)
#intModules0 =c("black","blue","brown","cyan","green","greenyellow","grey","magenta","pink","purple","red","salmon","tan","turquoise","yellow")
for (module in intModules0)
{# Select module probes
  modGenes0 = (moduleColors0==module)
  # Get their entrez ID codes
  modLLIDs0 = allLLIDs_0[modGenes0];
  # Write them into a file
  fileName0 = paste("entrezIDs-PFC-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs0), file = fileName0,
              row.names = FALSE, col.names = FALSE)}

fileName_0= paste("LocusLinkIDs-all-PFC.txt", sep="");
write.table(as.data.frame(allLLIDs_0), file = fileName_0,
            row.names = FALSE, col.names = FALSE)
#-------------------------------------------------------------------------------------------------
# NAC
intModules1 <- unique(moduleColors1)
for (module in intModules1)
{ modGenes1 = (moduleColors1==module)
modLLIDs1 = allLLIDs_1[modGenes1];
fileName1 = paste("entrezIDs-NAC-", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs1), file = fileName1,
            row.names = FALSE, col.names = FALSE)}
fileName_1= paste("LocusLinkIDs-all-NAC.txt", sep="");
write.table(as.data.frame(allLLIDs_1), file = fileName_1,
            row.names = FALSE, col.names = FALSE)
#----------------------------------------------------------------------------------------------------
# VTA
intModules2 <- unique(moduleColors2)
for (module in intModules2)
{modGenes2= (moduleColors2==module)
modLLIDs2 = allLLIDs_2[modGenes2];
fileName2 = paste("entrezIDs-VTA", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs2), file = fileName2,
            row.names = FALSE, col.names = FALSE)}
fileName_2= paste("LocusLinkIDs-all-VTA.txt", sep="");
write.table(as.data.frame(allLLIDs_2), file = fileName_2,
            row.names = FALSE, col.names = FALSE)
#--------------------------------------------------------------------------------------------------
