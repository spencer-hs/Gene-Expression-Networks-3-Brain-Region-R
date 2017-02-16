###NETWORK CONSTRUCTION OUTLINE MUST BE RAN FIRST, SEVERAL VARIABLES ARE USED
#IN SAME DIRECTORY AS BELOW

setwd("C:\\Users\\Godjira\\Dropbox\\Miles & Heather\\Heather's WGCNA\\WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE);
#__________________________________________________
PFC_Sscores_exp=read.csv("PFCdata_Genelist.csv")
NAC_Sscores_exp=read.csv("NACdata_Genelist.csv")
VTA_Sscores_exp=read.csv("VTAdata_Genelist.csv")
#__________________________________________________
datExpr_PFC = as.data.frame(t(PFC_Sscores_exp[, -c(1:29)]))
names(datExpr_PFC) = PFC_Sscores_exp$Probeset
rownames(datExpr_PFC) = names(PFC_Sscores_exp)[-c(1:29)]
datExpr_PFC<-datExpr_PFC[-(1:2),]
#___________________________________________________
datExpr_NAC = as.data.frame(t(NAC_Sscores_exp[, -c(1:21)]))
names(datExpr_NAC) = NAC_Sscores_exp$Probeset
rownames(datExpr_NAC) = names(NAC_Sscores_exp)[-c(1:21)]
datExpr_NAC<-datExpr_NAC[-(1:2),]
#____________________________________________________
datExpr_VTA = as.data.frame(t(VTA_Sscores_exp[, -c(1:21)]))
names(datExpr_VTA) = VTA_Sscores_exp$Probeset
rownames(datExpr_VTA) = names(VTA_Sscores_exp)[-c(1:21)]
datExpr_VTA<-datExpr_VTA[-(1:2),]
#_____________________________________________________
#Observing Data for average for detection of outliers, make sure no missing values
meanExpressionPFC=apply( datExpr_PFC,1,mean, na.rm=T)
meanExpressionNAC=apply( datExpr_NAC,1,mean, na.rm=T)
meanExpressionVTA=apply( datExpr_VTA,1,mean, na.rm=T)

NumberMissingPFC=apply( is.na(data.frame(datExpr_PFC)),1, sum)
NumberMissingNAC=apply( is.na(data.frame(datExpr_NAC)),1, sum)
NumberMissingVTA=apply( is.na(data.frame(datExpr_VTA)),1, sum)
NumberMissingPFC
NumberMissingNAC
NumberMissingVTA

sizeGrWindow(9, 5)
barplot(meanExpressionPFC,
        xlab = "Strains", ylab = "Mean expression",
        main ="Mean S Scores expression in PFC across strains",
        names.arg = c(1:27), cex.names = 0.7)
sizeGrWindow(9, 5)
barplot(meanExpressionNAC,
        xlab = "Strains", ylab = "Mean expression",
        main ="Mean S Scores expression in NAC across strains",
        names.arg = c(1:35), cex.names = 0.7)
sizeGrWindow(9, 5)
barplot(meanExpressionVTA,
        xlab = "Strains", ylab = "Mean expression",
        main ="Mean S Scores expression in VTA across strains",
        names.arg = c(1:35), cex.names = 0.7)
#_____________________________________________________________________
# here we define the adjacency matrix using soft thresholding with beta=6
ADJPFC=abs(cor(datExpr_PFC,use="p"))^6
k_PFC=softConnectivity(datE=datExpr_PFC,type="unsigned", power=12)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k_PFC)
scaleFreePlot(k_PFC, main="Check scale free topology\n")

ADJNAC=abs(cor(datExpr_NAC,use="p"))^6
k_NAC=softConnectivity(datE=datExpr_NAC,type="unsigned", power=9)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k_NAC)
scaleFreePlot(k_NAC, main="Check scale free topology\n")

ADJVTA=abs(cor(datExpr_VTA,use="p"))^6
k_VTA=softConnectivity(datE=datExpr_VTA,type="unsigned", power=9)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k_VTA)
scaleFreePlot(k_VTA, main="Check scale free topology\n")
#___________________________________________________________
#Dissimilarity ?
dissADJ_PFC=1-ADJPFC
dissADJ_NAC=1-ADJNAC
dissADJ_VTA=1-ADJVTA

# Plot the dendrogram with module colors
sizeGrWindow(10,5);
#Use plotdendro from WGCNA construction, run that prior so as to have net0,net1, net2
plotDendroAndColors(net0$dendrograms[[1]], mergedColors0[net0$blockGenes[[1]]], "Deep Split 3 PFC", dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)
plotDendroAndColors(net1$dendrograms[[1]], mergedColors1[net1$blockGenes[[1]]], "Deep Split 3 NAC", dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)
plotDendroAndColors(net2$dendrograms[[1]], mergedColors2[net2$blockGenes[[1]]], "Deep Split 3 VTA", dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)

datME_PFC=moduleEigengenes(datExpr_PFC, mergedColors0)$eigengenes
signif(cor(datME_PFC, use="p"), 2)

datME_NAC=moduleEigengenes(datExpr_NAC, mergedColors1)$eigengenes
signif(cor(datME_NAC, use="p"), 2)

datME_VTA=moduleEigengenes(datExpr_VTA, mergedColors2)$eigengenes
signif(cor(datME_VTA, use="p"), 2)

#___________________________________________________
# define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes
dissimME_PFC=(1-t(cor(datME_PFC, method="p")))/2
hclustdatME_PFC=hclust(as.dist(dissimME_PFC), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME_PFC, main="Clustering tree based of the module eigengenes")

dissimME_NAC=(1-t(cor(datME_NAC, method="p")))/2
hclustdatME_NAC=hclust(as.dist(dissimME_NAC), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME_NAC, main="Clustering tree based of the module eigengenes")

dissimME_VTA=(1-t(cor(datME_VTA, method="p")))/2
hclustdatME_VTA=hclust(as.dist(dissimME_VTA), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME_VTA, main="Clustering tree based of the module eigengenes")
#_________________________________________________________________________
#OPTIONAL HEAT MAP of module
sizeGrWindow(12,12)
par(mfrow=c(3,1), mar=c(0, 2, 4, 0))

# for the module
which.module="blue";
plotMat(t(scale(datExpr_PFC[,mergedColors0==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
#______________________________________________________________________

Alldegrees_PFC=intramodularConnectivity(ADJPFC, mergedColors0)
head(AlldegreesPFC)
write.csv(Alldegrees_PFC, "PFCconnectivitydegrees_Genelist.csv")

Alldegrees_NAC=intramodularConnectivity(ADJNAC,mergedColors1)
head(AlldegreesNAC)
write.csv(Alldegrees_NAC, "NACconnectivitydegrees_Genelist.csv")

Alldegrees_VTA=intramodularConnectivity(ADJVTA, mergedColors2)
head(AlldegreesVTA)
write.csv(Alldegrees_VTA, "VTAconnectivitydegrees_Genelist.csv")
#_______________________________________________________________
#Module Membership
datKME_PFC=signedKME(datExpr_PFC, datME_PFC, outputColumnName="MM.")
write.csv(datKME_PFC, "PFCmodulemembership_Genelist.csv")

datKME_NAC=signedKME(datExpr_NAC, datME_NAC, outputColumnName="MM.")
write.csv(datKME_NAC, "NACmodulemembership_Genelist.csv")

datKME_VTA=signedKME(datExpr_VTA, datME_VTA, outputColumnName="MM.")
write.csv(datKME_VTA, "VTAmodulemembership_Genelist.csv")
#__________________________________________________________________
#Optional for interest: module membership measures and intramodular connectivity

sizeGrWindow(8,6)
par(mfrow=c(2,2))

which.color="yellow";
restrictGenes=mergedColors0==which.color
verboseScatterplot(Alldegrees_PFC$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color="blue";
restrictGenes=mergedColors0==which.color
verboseScatterplot(Alldegrees_PFC$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")


