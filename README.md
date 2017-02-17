# Gene-Expression-Networks-3-Brain-Region-R

CONSTRUCTION OF GENE EXPRESSION NETWORKS ACROSS 3 BRAIN REGIONS USING R

The goal was to construct gene expression networks in B6 mice across three brain regions after acute exposure to ethanol. The script walks through the process of inputting each data set, filtering out at a determined significance level, binding the 3 regions together into a master gene list and removing redundancies before extracting the probeset_ids for each row in the data set.

From there, extraneous information (to this analysis) in the dataframe is removed, and a sampleTree is constructed and pickSoftThreshold employed to use R and WGCNA built in functionality to decide the best power for the analysis. Once completed for each brain region, the network construction function of WGCNA, blockwiseModules is begun, taking in the dataframe of each region, the block size determined by the memory and strength of the computing machine, the power selected previously, and other default setting as suggested by the tutorial methods of WGCNA.

Dendrograms are then plotted to show the users the options for the deepsplit and merge cut height. A table is given to specify the number of genes in each color, colors then set here to the name of a module. The annotation data, provided by Affymetrix, is loaded into R, and the names of the members in each module are matched to probeset_ID from the original data frame and against the given annotation data.

The unique Entrez IDs are then saved under a file named for the module it is a member of. Annotation data is saved into another file with the module name in the title.

CONNECTIVITY MEASURES BETWEEN MODULES

After network construction is run, it was of interest to know the similarity and dissimilarity between the modules in each region. The gene lists per brain region were loaded and extraneous information removed. The data was then observed for an average, to ensure there were no outliers, and no missing values. The average was plotted for better visualization. 

An adjacency matrix was constructed for each region, observing the same powers that were used in the construction of their networks. A histogram and Scale Free Topology were plotted. Dissimilarity matrixes were calculated by subtracting the adjacency matrix from one. 

The dendrograms from network construction were then run again so that the user had the variables net0, net1, net2 to proceed with. Module Eigengenes were then defined with a set significance and correlation. 

A dissimilarity measure was again defined, but now between the module Eigengenes, so as to keep track of the sign of the correlation. The Eigengenes were then plotted to a dendrogram, and a heatmap was constructed for specific modules to show its similarity to others. 

The measures of connectivity were written to excel files and a module membership list was constructed and also written to excel files. 

Module membership and intramodular connectivity measures could then be shown as a scatterplot. 

TRAIT INPUT AND CORRELATION


