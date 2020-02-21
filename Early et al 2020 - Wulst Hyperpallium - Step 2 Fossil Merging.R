
##NOTE: Must be run in same session as "Early et al 2020 - Wulst Hyperpallium - 
#Step 1 Fossil Predictions.R" right after running the PGLS in that script

##NOTE: All trees called here were prepared first by adding fossils in Mesquite.




##########################################################################
############__Preparing trees and datasets for BayesModelS__##############
##########################################################################




##Dinornis##

#check data from running "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantWulstlogged.csv and add Dinornis: 
#Dinornis_robustus FMNH_PA_35 7.356037191
#then save as DinornisWulstlogged.csv

#Upload and check Dinornis dataset
Dinornis.data<-read.csv("DinornisWulstlogged.csv", header=TRUE, row.names=1)
head(Dinornis.data)
sapply(Dinornis.data, class)

#Bring in tree from Mesquite
Dinornis.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Dinornis.tree<-treedata(Dinornis.tree, Dinornis.data,sort=TRUE,warnings=TRUE)$phy
Dinornis.data<-treedata(Dinornis.tree, Dinornis.data,sort=TRUE,warnings=TRUE)$data
plot(Dinornis.tree,cex=0.5)

#make sure dataset is in proper format
Dinornis.data<-as.data.frame(Dinornis.data)
Dinornis.data

#assign specimen numbers to vector
specimens.Dinornis.data<-data.frame(Dinornis.data$Specimen_No)
specimens.Dinornis.data

#output the sorted csv file for use in BayesModelS
write.csv(Dinornis.data, "Dinornis_Wulst_data.csv")

#output the sorted tree file for use in BayesModelS
write.tree(Dinornis.tree, file="Dinornis_Wulst_tree.tre")

#Make tree block with Dinornis by opening new .tre file in Notepad
#and copying the tree twice, then save as "Dinornis_Wulst_block.tre"


##Llallawavis##

#check data from running "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantWulstlogged.csv and add Llallawavis: 
#Llallawavis_scagliai MMP_5050 7.40278683
#then save as LlallawavisWulstlogged.csv

#Upload and check Llallawavis dataset
Llallawavis.data<-read.csv("LlallawavisWulstlogged.csv", header=TRUE, row.names=1)
head(Llallawavis.data)
sapply(Llallawavis.data, class)

#Bring in tree from Mesquite
Llallawavis.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Llallawavis.tree<-treedata(Llallawavis.tree, Llallawavis.data,sort=TRUE,warnings=TRUE)$phy
Llallawavis.data<-treedata(Llallawavis.tree, Llallawavis.data,sort=TRUE,warnings=TRUE)$data
plot(Llallawavis.tree,cex=0.5)

#make sure dataset is in proper format
Llallawavis.data<-as.data.frame(Llallawavis.data)
Llallawavis.data

#assign specimen numbers to vector
specimens.Llallawavis.data<-data.frame(Llallawavis.data$Specimen_No)
specimens.Llallawavis.data

#output the sorted csv file for use in BayesModelS
write.csv(Llallawavis.data, "Llallawavis_Wulst_data.csv")

#output the sorted tree file for use in BayesModelS
write.tree(Llallawavis.tree, file="Llallawavis_Wulst_tree.tre")

#Make tree block with Llallawavis by opening new .tre file in Notepad
#and copying the tree twice, then save as "Llallawavis_Wulst_block.tre"


##Miocene galliform##

#check data from running "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantWulstlogged.csv and add Miocene galliform: 
#Miocene_galliform AMNH_8269 4.851053224
#then save as MiocenegalliformWulstlogged.csv

#Upload and check Miocene galliform dataset
Miocenegalliform.data<-read.csv("MiocenegalliformWulstlogged.csv", header=TRUE, row.names=1)
head(Miocenegalliform.data)
sapply(Miocenegalliform.data, class)

#Bring in tree from Mesquite
Miocenegalliform.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Miocenegalliform.tree<-treedata(Miocenegalliform.tree, Miocenegalliform.data,sort=TRUE,warnings=TRUE)$phy
Miocenegalliform.data<-treedata(Miocenegalliform.tree, Miocenegalliform.data,sort=TRUE,warnings=TRUE)$data
plot(Miocenegalliform.tree,cex=0.5)

#make sure dataset is in proper format
Miocenegalliform.data<-as.data.frame(Miocenegalliform.data)
Miocenegalliform.data

#assign specimen numbers to vector
specimens.Miocenegalliform.data<-data.frame(Miocenegalliform.data$Specimen_No)
specimens.Miocenegalliform.data

#output the sorted csv file for use in BayesModelS
write.csv(Miocenegalliform.data, "Miocenegalliform_Wulst_data.csv")

#output the sorted tree file for use in BayesModelS
write.tree(Miocenegalliform.tree, file="Miocenegalliform_Wulst_tree.tre")

#Make tree block with Miocene galliform by opening new .tre file in Notepad
#and copying the tree twice, then save as "Miocenegalliform_Wulst_block.tre"


##Paraptenodytes##

#check data from running "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantWulstlogged.csv and add Paraptenodytes: 
#Paraptenodytes_antarcticus AMNH_3338 6.775651438
#then save as ParaptenodytesWulstlogged.csv

#Upload and check Paraptenodytes dataset
Paraptenodytes.data<-read.csv("ParaptenodytesWulstlogged.csv", header=TRUE, row.names=1)
head(Paraptenodytes.data)
sapply(Paraptenodytes.data, class)

#Bring in tree from Mesquite
Paraptenodytes.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Paraptenodytes.tree<-treedata(Paraptenodytes.tree, Paraptenodytes.data,sort=TRUE,warnings=TRUE)$phy
Paraptenodytes.data<-treedata(Paraptenodytes.tree, Paraptenodytes.data,sort=TRUE,warnings=TRUE)$data
plot(Paraptenodytes.tree,cex=0.5)

#make sure dataset is in proper format
Paraptenodytes.data<-as.data.frame(Paraptenodytes.data)
Paraptenodytes.data

#assign specimen numbers to vector
specimens.Paraptenodytes.data<-data.frame(Paraptenodytes.data$Specimen_No)
specimens.Paraptenodytes.data

#output the sorted csv file for use in BayesModelS
write.csv(Paraptenodytes.data, "Paraptenodytes_Wulst_data.csv")

#output the sorted tree file for use in BayesModelS
write.tree(Paraptenodytes.tree, file="Paraptenodytes_Wulst_tree.tre")

#Make tree block with Paraptenodytes by opening new .tre file in Notepad
#and copying the tree twice, then save as "Paraptenodytes_Wulst_block.tre"


##Psilopterus##

#check data from running "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantWulstlogged.csv and add Psilopterus: 
#Psilopterus_lemoinei AMNH_9257 6.752620782
#then save as PsilopterusWulstlogged.csv

#Upload and check Psilopterus dataset
Psilopterus.data<-read.csv("PsilopterusWulstlogged.csv", header=TRUE, row.names=1)
head(Psilopterus.data)
sapply(Psilopterus.data, class)

#Bring in tree from Mesquite
Psilopterus.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Psilopterus.tree<-treedata(Psilopterus.tree, Psilopterus.data,sort=TRUE,warnings=TRUE)$phy
Psilopterus.data<-treedata(Psilopterus.tree, Psilopterus.data,sort=TRUE,warnings=TRUE)$data
plot(Psilopterus.tree,cex=0.5)

#make sure dataset is in proper format
Psilopterus.data<-as.data.frame(Psilopterus.data)
Psilopterus.data

#assign specimen numbers to vector
specimens.Psilopterus.data<-data.frame(Psilopterus.data$Specimen_No)
specimens.Psilopterus.data

#output the sorted csv file for use in BayesModelS
write.csv(Psilopterus.data, "Psilopterus_Wulst_data.csv")

#output the sorted tree file for use in BayesModelS
write.tree(Psilopterus.tree, file="Psilopterus_Wulst_tree.tre")

#Make tree block with Psilopterus by opening new .tre file in Notepad
#and copying the tree twice, then save as "Psilopterus_Wulst_block.tre"