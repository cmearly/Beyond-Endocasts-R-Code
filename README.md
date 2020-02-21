# Beyond-Endocasts-R-Code
R code used to run the analyses published in Early et al. 2020 Beyond endocasts

"BayesModelS_v24.R" - This is the version of BayesModelS (Nunn and Zhu, 2014) used in the analyses in this publication.


The following R scripts should be used, in order and as indicated by annotations within the script, to run analyses on optic tectum and optic lobe size.

"Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R" - This is the code for predicting optic tectum volume in extinct birds based on optic lobe surface area using BayesModelS (Nunn and Zhu, 2014).
"Early et al 2020 - Optic Lobe Optic Tectum - Step 2 Fossil Merging.R" - This is the code for incorporating fossils into datasets for use in the previous script.
"Early et al 2020 - Optic Lobe Optic Tectum - Step 3 phyANCOVAs.R" - This is the code for using the R package evomap (Smaers and Mongle, 2014) to test for differences in the optic tectum volumes of fossils predicted in the previous scripts with BayesModelS.


The following R scripts should be used, in order and as indicated by annotations within the script, to run analyses on hyperpallium and Wulst size.

"Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R" - This is the code for predicting hyperpallium volume in extinct birds based on Wulst surface area using BayesModelS (Nunn and Zhu, 2014).
"Early et al 2020 - Wulst Hyperpallium - Step 2 Fossil Merging.R" - This is the code for incorporating fossils into datasets for use in the previous script.
"Early et al 2020 - Wulst Hyperpallium - Step 3 phyANCOVAs.R" - This is the code for using the R package evomap (Smaers and Mongle, 2014) to test for differences in the hyperpallium volumes of fossils predicted in the previous scripts with BayesModelS.

The following files are called in the above analyses:
"opticlobes4.csv" and "wulst3.csv" have all data for extant birds.
"allWulstandhyperpallium.csv," "allhyperpallium.csv," "allopticlobesandtectum.csv," and "alloptictectum.csv" have all data for extant birds plus the values for extinct birds predicted with BayesModelS (Nunn and Zhu, 2014).
"alloptictectumandhyperpalliumtaxa.tre" is the tree used in all analyses. It was compiled following the methods described in the manuscript.


Citations
Nunn CL, Zhu L. 2014. Phylogenetic prediction to identify "evolutionary singularities." In Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Ed. LZ Garamszegi. Springer: Berlin/Heidelberg, Germany. pp. 481-514.
Smaers JB, Mongle CS. 2014. evomap: R package for the evolutionary mapping of continuous traits. https://github.com/JeroenSmaers/evomap.
