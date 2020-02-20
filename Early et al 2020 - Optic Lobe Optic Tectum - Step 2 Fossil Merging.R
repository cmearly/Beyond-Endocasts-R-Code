
##NOTE: Must be run in same session as "Early et al 2020 - Optic Lobe Optic Tectum - 
#Step 1 Fossil Predictions.R" right after running the PGLS in that script

##NOTE: All trees called here were prepared first by adding fossils in Mesquite.




##########################################################################
############__Preparing trees and datasets for BayesModelS__##############
##########################################################################



##Archaeopteryx##


#check data from running "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open "S_extantopticlogged.csv" and add Archaeopteryx: 
#Archaeopteryx_lithographica NHMUK_PV_OR_37001 4.74673051
#then save as Archaeopteryxopticlogged.csv

#Upload and check Archaeopteryx dataset
Archaeopteryx.data<-read.csv("Archaeopteryxopticlogged.csv", header=TRUE, row.names=1)
head(Archaeopteryx.data)
sapply(Archaeopteryx.data, class)

#Bring in tree from Mesquite
Archaeopteryx.tree <- read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Archaeopteryx.tree<-treedata(Archaeopteryx.tree, Archaeopteryx.data,sort=TRUE,warnings=TRUE)$phy
Archaeopteryx.data<-treedata(Archaeopteryx.tree, Archaeopteryx.data,sort=TRUE,warnings=TRUE)$data

#make sure dataset is in proper format
Archaeopteryx.data<-as.data.frame(Archaeopteryx.data)
Archaeopteryx.data

plot(Archaeopteryx.tree, cex=0.5)

#assign specimen numbers to vector
specimens.Archaeopteryx.data<-data.frame(Archaeopteryx.data$Specimen_No)
specimens.Archaeopteryx.data

#output the sorted csv file for use in BayesModelS
write.csv(Archaeopteryx.data, "Archaeopteryx_data.csv")

#output the sorted tree file for use in BayesModels
write.tree(Archaeopteryx.tree, file="Archaeopteryx_tree.tre")

#Make tree block with Archaeopteryx by opening new .tre file in Notepad
#and copying the tree twice, then save as "Archaeopteryx_block.tre"


##Dinornis##

#check data from running "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantopticlogged.csv and add Dinornis: 
#Dinornis_robustus FMNH_PA_35 5.269341594
#then save as Dinornisopticlogged.csv

#Upload and check Dinornis dataset
Dinornis.data<-read.csv("Dinornisopticlogged.csv", header=TRUE, row.names=1)
head(Dinornis.data)
sapply(Dinornis.data, class)

#Bring in tree from Mesquite
Dinornis.tree <- read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Dinornis.tree<-treedata(Dinornis.tree, Dinornis.data,sort=TRUE,warnings=TRUE)$phy
Dinornis.data<-treedata(Dinornis.tree, Dinornis.data,sort=TRUE,warnings=TRUE)$data

#make sure dataset is in proper format
Dinornis.data<-as.data.frame(Dinornis.data)
Dinornis.data

#assign specimen numbers to vector
specimens.Dinornis.data<-data.frame(Dinornis.data$Specimen_No)
specimens.Dinornis.data

#output the sorted csv file for use in BayesModels
write.csv(Dinornis.data, "Dinornis_data.csv")

#output the sorted tree file for use in BayesModels
write.tree(Dinornis.tree, file="Dinornis_tree.tre")

#Make tree block with Dinornis by opening new .tre file in Notepad
#and copying the tree twice, then save as "Dinornis_block.tre"


##Lithornis##

#check data from running "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantopticlogged.csv and add Lithornis: 
#Lithornis_plebius USNM_336534 5.150860897
#then save as Lithornisopticlogged.csv

#Upload and check Lithornis dataset
Lithornis.data<-read.csv("Lithornisopticlogged.csv", header=TRUE, row.names=1)
head(Lithornis.data)
sapply(Lithornis.data, class)

#Bring in tree from Mesquite
Lithornis.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Lithornis.tree<-treedata(Lithornis.tree, Lithornis.data,sort=TRUE,warnings=TRUE)$phy
Lithornis.data<-treedata(Lithornis.tree, Lithornis.data,sort=TRUE,warnings=TRUE)$data

#make sure dataset is in proper format
Lithornis.data<-as.data.frame(Lithornis.data)
Lithornis.data

#assign specimen numbers to vector
specimens.Lithornis.data<-data.frame(Lithornis.data$Specimen_No)
specimens.Lithornis.data

#output the sorted csv file for use in BayesModels
write.csv(Lithornis.data, "Lithornis_data.csv")

#output the sorted tree file for use in BayesModels
write.tree(Lithornis.tree, file="Lithornis_tree.tre")

#Make tree block with Lithornis by opening new .tre file in Notepad
#and copying the tree twice, then save as "Lithornis_block.tre"


##Llallawavis##

#check data from running "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantopticlogged.csv and add Llallawavis: 
#Llallawavis_scagliai MMP_5050 6.498658585
#then save as Llallawavisopticlogged.csv

#Upload and check Llallawavis dataset
Llallawavis.data<-read.csv("Llallawavisopticlogged.csv", header=TRUE, row.names=1)
head(Llallawavis.data)
sapply(Llallawavis.data, class)

#Bring in tree from Mesquite
Llallawavis.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Llallawavis.tree<-treedata(Llallawavis.tree, Llallawavis.data,sort=TRUE,warnings=TRUE)$phy
Llallawavis.data<-treedata(Llallawavis.tree, Llallawavis.data,sort=TRUE,warnings=TRUE)$data

#make sure dataset is in proper format
Llallawavis.data<-as.data.frame(Llallawavis.data)
Llallawavis.data

#assign specimen numbers to vector
specimens.Llallawavis.data<-data.frame(Llallawavis.data$Specimen_No)
specimens.Llallawavis.data

#output the sorted csv file for use in BayesModels
write.csv(Llallawavis.data, "Llallawavis_data.csv")

#output the sorted tree file for use in BayesModels
write.tree(Llallawavis.tree, file="Llallawavis_tree.tre")

#Make tree block with Llallawavis by opening new .tre file in Notepad
#and copying the tree twice, then save as "Llallawavis_block.tre"


##Miocene galliform##

#check data from running "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantopticlogged.csv and add Miocene galliform: 
#Miocene_galliform AMNH_8269 5.152736527
#then save as Miocenegalliformopticlogged.csv

#Upload and check Miocene galliform dataset
Miocenegalliform.data<-read.csv("Miocenegalliformopticlogged.csv", header=TRUE, row.names=1)
head(Miocenegalliform.data)
sapply(Miocenegalliform.data, class)

#Bring in tree from Mesquite
Miocenegalliform.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Miocenegalliform.tree<-treedata(Miocenegalliform.tree, Miocenegalliform.data,sort=TRUE,warnings=TRUE)$phy
Miocenegalliform.data<-treedata(Miocenegalliform.tree, Miocenegalliform.data,sort=TRUE,warnings=TRUE)$data

#make sure dataset is in proper format
Miocenegalliform.data<-as.data.frame(Miocenegalliform.data)
Miocenegalliform.data

#assign specimen numbers to vector
specimens.Miocenegalliform.data<-data.frame(Miocenegalliform.data$Specimen_No)
specimens.Miocenegalliform.data

#output the sorted csv file for use in BayesModels
write.csv(Miocenegalliform.data, "Miocenegalliform_data.csv")

#output the sorted tree file for use in BayesModels
write.tree(Miocenegalliform.tree, file="Miocenegalliform_tree.tre")

#Make tree block with Miocene galliform by opening new .tre file in Notepad
#and copying the tree twice, then save as "Miocenegalliform_block.tre"


##Paraptenodytes##

#check data from running "Ch 2 Optic Lobe Optic Tectum Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantopticlogged.csv and add Paraptenodytes: 
#Paraptenodytes_antarcticus AMNH_3338 6.118145655
#then save as Paraptenodytesopticlogged.csv

#Upload and check Paraptenodytes dataset
Paraptenodytes.data<-read.csv("Paraptenodytesopticlogged.csv", header=TRUE, row.names=1)
head(Paraptenodytes.data)
sapply(Paraptenodytes.data, class)

#Bring in tree from Mesquite
Paraptenodytes.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Paraptenodytes.tree<-treedata(Paraptenodytes.tree, Paraptenodytes.data,sort=TRUE,warnings=TRUE)$phy
Paraptenodytes.data<-treedata(Paraptenodytes.tree, Paraptenodytes.data,sort=TRUE,warnings=TRUE)$data

#make sure dataset is in proper format
Paraptenodytes.data<-as.data.frame(Paraptenodytes.data)
Paraptenodytes.data

#assign specimen numbers to vector
specimens.Paraptenodytes.data<-data.frame(Paraptenodytes.data$Specimen_No)
specimens.Paraptenodytes.data

#output the sorted csv file for use in BayesModels
write.csv(Paraptenodytes.data, "Paraptenodytes_data.csv")

#output the sorted tree file for use in BayesModels
write.tree(Paraptenodytes.tree, file="Paraptenodytes_tree.tre")

#Make tree block with Paraptenodytes by opening new .tre file in Notepad
#and copying the tree twice, then save as "Paraptenodytes_block.tre"


##Psilopterus##

#check data from running "Ch 2 Optic Lobe Optic Tectum Fossil Predictions.R"
head(sort.data)
nrow(sort.data)

#Open S_extantopticlogged.csv and add Psilopterus: 
#Psilopterus_lemoinei AMNH_9257 6.11916491
#then save as Psilopterusopticlogged.csv

#Upload and check Psilopterus dataset
Psilopterus.data<-read.csv("Psilopterusopticlogged.csv", header=TRUE, row.names=1)
head(Psilopterus.data)
sapply(Psilopterus.data, class)

#Bring in tree from Mesquite
Psilopterus.tree <-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#Sort tree and data to same order
Psilopterus.tree<-treedata(Psilopterus.tree, Psilopterus.data,sort=TRUE,warnings=TRUE)$phy
Psilopterus.data<-treedata(Psilopterus.tree, Psilopterus.data,sort=TRUE,warnings=TRUE)$data

#make sure dataset is in proper format
Psilopterus.data<-as.data.frame(Psilopterus.data)
Psilopterus.data

#assign specimen numbers to vector
specimens.Psilopterus.data<-data.frame(Psilopterus.data$Specimen_No)
specimens.Psilopterus.data

#output the sorted csv file for use in BayesModels
write.csv(Psilopterus.data, "Psilopterus_data.csv")

#output the sorted tree file for use in BayesModels
write.tree(Psilopterus.tree, file="Psilopterus_tree.tre")

#Make tree block with Psilopterus by opening new .tre file in Notepad
#and copying the tree twice, then save as "Psilopterus_block.tre"