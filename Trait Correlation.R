lnames_0 = load(file = "PFC_Anxiety_dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames_0
# Load network data saved in the second part.
lnames_0 <- load(file = "PFCdata_autonetwork.RData")
lnames_0
nGenespfc = ncol(datExpr0);
nSamplespfc = nrow(datExpr0);
# Recalculate MEs with color labels
MEs_0 = moduleEigengenes(datExpr_0, moduleColors0)$eigengenes
MEs0 = orderMEs(MEs_0)
moduleTraitCor0 = cor(MEs0, datTraits0);
moduleTraitPvalue0 = corPvalueStudent(moduleTraitCor0, nSamplespfc);
textMatrix0 = paste(signif(moduleTraitCor0, 2), "\n(",
                    signif(moduleTraitPvalue0, 1), ")", sep = "");
dim(textMatrix0) = dim(moduleTraitCor0)
par(mar = c(3.5,6,0.89,0.2));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor0,
               xLabels = names(datTraits0),
               yLabels = names(MEs0),
               ySymbols = names(MEs0),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix0,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-anxiety relationships-PFC"))
#-------------------------------------------------------------------------
#NAC
lnames_1 = load(file = "NAC_Anxiety_dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames_1
# Load network data saved in the second part.
lnames_1 <- load(file = "NACdata_autonetwork.RData")
lnames_1
nGenesnac = ncol(datExpr_1_corr);
nSamplesnac = nrow(datExpr_1_corr);
# Recalculate MEs with color labels
MEs_1 = moduleEigengenes(datExpr_1_corr, moduleColors1)$eigengenes
MEs1 = orderMEs(MEs_1)
moduleTraitCor1 = cor(MEs1, datTraits1);
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamplesnac);
textMatrix1 = paste(signif(moduleTraitCor1, 2), "\n(",
                    signif(moduleTraitPvalue1, 1), ")", sep = "");
dim(textMatrix1) = dim(moduleTraitCor1)
par(mar = c(1.5,6,0.89,0.2));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor1,
                                   xLabels = names(datTraits1),
                                   yLabels = names(MEs1),
                                   ySymbols = names(MEs1),
                                   colorLabels = FALSE,
                                   colors = greenWhiteRed(50),
                                   textMatrix = textMatrix1,
                                   setStdMargins = FALSE,
                                   cex.text = 0.5,
                                   zlim = c(-1,1),
                                   main = paste("Module-anxiety relationships-NAC"))
##-------------
##VTA
lnames_2 = load(file = "VTA_Anxiety_dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames_2
# Load network data saved in the second part.
lnames_2 <- load(file = "VTAdata_autonetwork.RData")
lnames_2
nGenesvta = ncol(datExpr_2_corr);
nSamplesvta = nrow(datExpr_2_corr);
# Recalculate MEs with color labels
MEs_2 = moduleEigengenes(datExpr_2_corr, moduleColors2)$eigengenes
MEs2 = orderMEs(MEs_2)
moduleTraitCor2 = cor(MEs2, datTraits2);
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamplesvta);
textMatrix2 = paste(signif(moduleTraitCor2, 2), "\n(",
                    signif(moduleTraitPvalue2, 1), ")", sep = "");
dim(textMatrix2) = dim(moduleTraitCor2)
par(mar = c(1.5,6,0.89,0.2));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(datTraits2),
               yLabels = names(MEs2),
               ySymbols = names(MEs2),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-anxiety relationships-VTA"))