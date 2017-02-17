# Gene-Expression-Networks-3-Brain-Region-R

Construction of Gene Expression Network across 3 brain regions in the presence of Ethanol using WGCNA in R

The goal was to construct gene expression networks in B6 mice across three brain regions after acute exposure to ethanol. The script walks through the process of inputting each data set, filtering out at a determined significance level, binding the 3 regions together into a master gene list and removing redundancies before extracting the probeset_ids for each row in the data set.

From there, extraneous information (to this analysis) in the dataframe is removed, and a sampleTree is constructed and pickSoftThreshold employed to use R and WGCNA built in functionality to decide the best power for the analysis. Once completed for each brain region, the network construction function of WGCNA, blockwiseModules is begun, taking in the dataframe of each region, the block size determined by the memory and strength of the computing machine, the power selected previously, and other default setting as suggested by the tutorial methods of WGCNA.

Dendrograms are then plotted to show the users the options for the deepsplit and merge cut height. A table is given to specify the number of genes in each color, colors then set here to the name of a module. The annotation data, provided by Affymetrix, is loaded into R, and the names of the members in each module are matched to probeset_ID from the original data frame and against the given annotation data.

The unique Entrez IDs are then saved under a file named for the module it is a member of. Annotation data is saved into another file with the module name in the title.



