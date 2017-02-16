###PFC
traitData = read.csv("Ethanol BXD phenotypes_tester.csv",header=TRUE, sep=",");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(2)];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
pfcSamples = rownames(datExpr_0);
strains= allTraits$Strain
traitRows0 = match(pfcSamples, strains);
datTraits0 = allTraits[traitRows0, -1];
rownames(datTraits0) = allTraits[traitRows0, 1];
save(datExpr_0, datTraits0, file = "PFC_Anxiety_dataInput.RData")
#--------------------------------------------------------
##NAC
traitData = read.csv("Ethanol BXD phenotypes_tester.csv",header=TRUE, sep=",");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(2)];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
datExpr_1_corr=datExpr_1[-c(27,30), ]
##datExpr_1_corr=datExpr_1[-c(19,22), ]
nacSamples = rownames(datExpr_1_corr);

strains= allTraits$Strain
traitRows1 = match(nacSamples, strains);
datTraits1 = allTraits[traitRows1, -1];

rownames(datTraits1) = allTraits[traitRows1, 1];
save(datExpr_1, datTraits1, file = "NAC_Anxiety_dataInput.RData")
##------------------------------------------------------
##VTA
traitData = read.csv("Ethanol BXD phenotypes_tester.csv",header=TRUE, sep=",");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(2)];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
datExpr_2_corr=datExpr_2[-c(27,30), ]
nacSamples = rownames(datExpr_2_corr);

strains= allTraits$Strain
traitRows2 = match(nacSamples, strains);
datTraits2 = allTraits[traitRows2, -1];

rownames(datTraits2) = allTraits[traitRows2, 1];
save(datExpr_2, datTraits2, file = "VTA_Anxiety_dataInput.RData")
##------------------------------------------------------