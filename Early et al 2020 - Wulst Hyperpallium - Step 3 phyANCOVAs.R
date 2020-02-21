##### Phylogenetic ANCOVAs for Hyperpallium and Wulst of Fossils #####

require(phytools)
require(geiger)
require(nlme)
require(evomap)
library(caper)
library(nlme)
library(car)

### !!! Don't forget to set your working directory! !!! ###



##########################################################################
#######################__Preparing for analysis__######################### 
##############__hyperpallium volume vs. brain rest volume__###############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and hyperpallium
#volumes, including the predicted hyperpallium volumes from the BayesModelS
#code done in "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R"
data<-read.csv("allhyperpallium.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"brain_rest_vol"; Y<-"hyperpallium_vol"

#limit the dataset to your X and Y variables of interest
data<-data[,c(which(colnames(data)==X),which(colnames(data)==Y)),drop=F]
data<-na.omit(data) 

#log transform the data (note: R defaults to ln)
data<-log(data)

#bring in the tree file
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#sort the tree and the data for each other
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-treedata(tree,data,sort=T,warnings=T)$data   #match the data to the tree
plot(tree,cex=0.5)

#make the sorted data a dataframe
data<-as.data.frame(data)

#if desired, change the column names of the data for repeatability
colnames(data)<-c("Independent","Dependent")     #Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets.


##Explore the data##

#do a simple linear regression between variables
fit_H_B<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_H_B)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
H_vs_B_Data<-data

#get variables
Xvar<-H_vs_B_Data[,which(colnames(H_vs_B_Data)==paste("Independent")), drop=F]
Yvar<-H_vs_B_Data[,which(colnames(H_vs_B_Data)==paste("Dependent")), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)
model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML")
model
summary(model)

PGLS_HvsBR_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_HvsBR_Lambda

PGLS_HvsBR_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_HvsBR_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species.



##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset :__#######
##############__hyperpallium volume vs. brain rest volume__###############
##########################################################################



##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:93,95:97)
length(Others)
Dinornis<-c(94)
length(Dinornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Dinornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Dinornis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Llallawavis vs. all others

#assign group identities
Others<-c(1:34,36:97)
length(Others)
Llallawavis<-c(35)
length(Llallawavis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Llallawavis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Llallawavis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Miocene_galliform vs. all others

#assign group identities
Others<-c(1:81,83:97)
length(Others)
Miocene_galliform<-c(82)
length(Miocene_galliform)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Miocene_galliform)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Miocene_galliform]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Paraptenodytes vs. all others

#assign group identities
Others<-c(1:53,55:97)
length(Others)
Paraptenodytes<-c(54)
length(Paraptenodytes)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Paraptenodytes)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Paraptenodytes]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant:  
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Psilopterus vs. all others

#assign group identities
Others<-c(1:33,35:97)
length(Others)
Psilopterus<-c(34)
length(Psilopterus)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Psilopterus)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Psilopterus]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)






##########################################################################
#######################__Preparing for analysis__######################### 
########__Wulst surface area vs. endocast rest surface area__########
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and hyperpallium
#volumes, including the predicted hyperpallium volumes from the BayesModelS
#code done in "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allWulstandhyperpallium.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"endo_rest_SA"; Y<-"Wulst_SA"

#limit the dataset to your X and Y variables of interest
data<-data[,c(which(colnames(data)==X),which(colnames(data)==Y)),drop=F]
data<-na.omit(data) 

#log transform the data (note: R defaults to ln)
data<-log(data)

#bring in the tree file
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#sort the tree and the data for each other
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-treedata(tree,data,sort=T,warnings=T)$data   #match the data to the tree

#check the tree
plot(tree,cex=0.5)

#make the sorted data a dataframe
data<-as.data.frame(data)

#if desired, change the column names of the data for repeatability
colnames(data)<-c("Independent","Dependent")     #Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets.


##Explore the data##

#do a simple linear regression between variables
fit_W_E<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_W_E)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
W_vs_E_Data<-data

#get variables
Xvar<-W_vs_E_Data[,which(colnames(W_vs_E_Data)==paste("Independent")), drop=F]
Yvar<-W_vs_E_Data[,which(colnames(W_vs_E_Data)==paste("Dependent")), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)
model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML")
model
summary(model)

PGLS_WvsER_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_WvsER_Lambda

PGLS_WvsER_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_WvsER_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species.




##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__######## 
########__Wulst surface area vs. endocast rest surface area__########
##########################################################################



##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:26,28:30)
length(Others)
Dinornis<-c(27)
length(Dinornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Dinornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Dinornis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Llallawavis vs. all others

#assign group identities
Others<-c(1:6,8:30)
length(Others)
Llallawavis<-c(7)
length(Llallawavis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Llallawavis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Llallawavis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Miocene_galliform vs. all others

#assign group identities
Others<-c(1:20,22:30)
length(Others)
Miocene_galliform<-c(21)
length(Miocene_galliform)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Miocene_galliform)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Miocene_galliform]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Paraptenodytes vs. all others

#assign group identities
Others<-c(1:12,14:30)
length(Others)
Paraptenodytes<-c(13)
length(Paraptenodytes)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Paraptenodytes)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Paraptenodytes]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant:  
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Psilopterus vs. all others

#assign group identities
Others<-c(1:5,7:30)
length(Others)
Psilopterus<-c(6)
length(Psilopterus)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Psilopterus)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Psilopterus]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)






##########################################################################
#######################__Preparing for analysis__######################### 
#############__Wulst surface area vs. foramen magnum area__###############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and hyperpallium
#volumes, including the predicted hyperpallium volumes from the BayesModelS
#code done in "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allWulstandhyperpallium.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_SA"; Y<-"Wulst_SA"

#limit the dataset to your X and Y variables of interest
data<-data[,c(which(colnames(data)==X),which(colnames(data)==Y)),drop=F]
data<-na.omit(data) 

#log transform the data (note: R defaults to ln)
data<-log(data)

#bring in the tree file
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#sort the tree and the data for each other
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-treedata(tree,data,sort=T,warnings=T)$data   #match the data to the tree

#check the tree
plot(tree,cex=0.5)

#make the sorted data a dataframe
data<-as.data.frame(data)

#if desired, change the column names of the data for repeatability
colnames(data)<-c("Independent","Dependent")     #Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets.


##Explore the data##

#do a simple linear regression between variables
fit_W_FM<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_W_FM)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
W_vs_FM_Data<-data

#get variables
Xvar<-W_vs_FM_Data[,which(colnames(W_vs_FM_Data)==paste("Independent")), drop=F]
Yvar<-W_vs_FM_Data[,which(colnames(W_vs_FM_Data)==paste("Dependent")), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)
model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML")
model
summary(model)

PGLS_WvsFM_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_WvsFM_Lambda

PGLS_WvsFM_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_WvsFM_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species.




##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__########
#############__Wulst surface area vs. foramen magnum area__###############
##########################################################################



##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:26,28:30)
length(Others)
Dinornis<-c(27)
length(Dinornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Dinornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Dinornis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Llallawavis vs. all others

#assign group identities
Others<-c(1:6,8:30)
length(Others)
Llallawavis<-c(7)
length(Llallawavis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Llallawavis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Llallawavis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Miocene_galliform vs. all others

#assign group identities
Others<-c(1:20,22:30)
length(Others)
Miocene_galliform<-c(21)
length(Miocene_galliform)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Miocene_galliform)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Miocene_galliform]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Paraptenodytes vs. all others

#assign group identities
Others<-c(1:12,14:30)
length(Others)
Paraptenodytes<-c(13)
length(Paraptenodytes)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Paraptenodytes)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Paraptenodytes]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant:  
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Psilopterus vs. all others

#assign group identities
Others<-c(1:5,7:30)
length(Others)
Psilopterus<-c(6)
length(Psilopterus)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Psilopterus)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Psilopterus]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)






##########################################################################
#######################__Preparing for analysis__######################### 
#############__Wulst surface area vs. foramen magnum width__##############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModelS
#code done in "Ch 2 Optic Lobe Optic Tectum Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allWulstandhyperpallium.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_w"; Y<-"Wulst_SA"

#limit the dataset to your X and Y variables of interest
data<-data[,c(which(colnames(data)==X),which(colnames(data)==Y)),drop=F]
data<-na.omit(data) 

#log transform the data (note: R defaults to ln)
data<-log(data)

#bring in the tree file
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#sort the tree and the data for each other
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-treedata(tree,data,sort=T,warnings=T)$data   #match the data to the tree

#check the tree
plot(tree,cex=0.5)

#make the sorted data a dataframe
data<-as.data.frame(data)

#if desired, change the column names of the data for repeatability
colnames(data)<-c("Independent","Dependent")     #Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets.


##Explore the data##

#do a simple linear regression between variables
fit_W_FMw<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_W_FMw)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
W_vs_FMw_Data<-data

#get variables
Xvar<-W_vs_FMw_Data[,which(colnames(W_vs_FMw_Data)==paste("Independent")), drop=F]
Yvar<-W_vs_FMw_Data[,which(colnames(W_vs_FMw_Data)==paste("Dependent")), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)
model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML")
model
summary(model)

PGLS_WvsFMw_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_WvsFMw_Lambda

PGLS_WvsFMw_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_WvsFMw_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species.




##########################################################################
#########__Running phyANCOVA code with dataset minus Dinornis:__########## 
#############__Wulst surface area vs. foramen magnum width__##############
##########################################################################



##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:26,28:30)
length(Others)
Dinornis<-c(27)
length(Dinornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Dinornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Dinornis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Llallawavis vs. all others

#assign group identities
Others<-c(1:6,8:30)
length(Others)
Llallawavis<-c(7)
length(Llallawavis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Llallawavis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Llallawavis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Miocene_galliform vs. all others

#assign group identities
Others<-c(1:20,22:30)
length(Others)
Miocene_galliform<-c(21)
length(Miocene_galliform)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Miocene_galliform)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Miocene_galliform]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Paraptenodytes vs. all others

#assign group identities
Others<-c(1:12,14:30)
length(Others)
Paraptenodytes<-c(13)
length(Paraptenodytes)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Paraptenodytes)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Paraptenodytes]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant:  
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Psilopterus vs. all others

#assign group identities
Others<-c(1:5,7:30)
length(Others)
Psilopterus<-c(6)
length(Psilopterus)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Psilopterus)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Psilopterus]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)






##########################################################################
#######################__Preparing for analysis__######################### 
#############__hyperpallium volume vs. foramen magnum area__##############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and hyperpallium
#volumes, including the predicted hyperpallium volumes from the BayesModelS
#code done in "Early et al 2020 - Wulst Hyperpallium - Step 1 Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allWulstandhyperpallium.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_SA"; Y<-"hyperpallium_vol"

#limit the dataset to your X and Y variables of interest
data<-data[,c(which(colnames(data)==X),which(colnames(data)==Y)),drop=F]
data<-na.omit(data) 

#log transform the data (note: R defaults to ln)
data<-log(data)

#bring in the tree file
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#sort the tree and the data for each other
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-treedata(tree,data,sort=T,warnings=T)$data   #match the data to the tree

#check the tree
plot(tree,cex=0.5)

#make the sorted data a dataframe
data<-as.data.frame(data)

#if desired, change the column names of the data for repeatability
colnames(data)<-c("Independent","Dependent")     #Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets.


##Explore the data##

#do a simple linear regression between variables
fit_H_FM<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_H_FM)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
H_vs_FM_Data<-data

#get variables
Xvar<-H_vs_FM_Data[,which(colnames(H_vs_FM_Data)==paste("Independent")), drop=F]
Yvar<-H_vs_FM_Data[,which(colnames(H_vs_FM_Data)==paste("Dependent")), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)
model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML")
model
summary(model)

PGLS_HvsFM_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_HvsFM_Lambda

PGLS_HvsFM_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_HvsFM_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species.




##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__######## 
#############__hyperpallium volume vs. foramen magnum area__##############
##########################################################################



##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:26,28:30)
length(Others)
Dinornis<-c(27)
length(Dinornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Dinornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Dinornis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Llallawavis vs. all others

#assign group identities
Others<-c(1:6,8:30)
length(Others)
Llallawavis<-c(7)
length(Llallawavis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Llallawavis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Llallawavis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Miocene_galliform vs. all others

#assign group identities
Others<-c(1:20,22:30)
length(Others)
Miocene_galliform<-c(21)
length(Miocene_galliform)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Miocene_galliform)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Miocene_galliform]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Paraptenodytes vs. all others

#assign group identities
Others<-c(1:12,14:30)
length(Others)
Paraptenodytes<-c(13)
length(Paraptenodytes)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Paraptenodytes)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Paraptenodytes]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant:  
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Psilopterus vs. all others

#assign group identities
Others<-c(1:5,7:30)
length(Others)
Psilopterus<-c(6)
length(Psilopterus)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Psilopterus)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Psilopterus]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)






##########################################################################
#######################__Preparing for analysis__######################### 
#############__hyperpallium volume vs. foramen magnum width__#############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModelS
#code done in "Ch 2 Optic Lobe Optic Tectum Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allWulstandhyperpallium.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_w"; Y<-"hyperpallium_vol"

#limit the dataset to your X and Y variables of interest
data<-data[,c(which(colnames(data)==X),which(colnames(data)==Y)),drop=F]
data<-na.omit(data) 

#log transform the data (note: R defaults to ln)
data<-log(data)

#bring in the tree file
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#sort the tree and the data for each other
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-treedata(tree,data,sort=T,warnings=T)$data   #match the data to the tree

#check the tree
plot(tree,cex=0.5)

#make the sorted data a dataframe
data<-as.data.frame(data)

#if desired, change the column names of the data for repeatability
colnames(data)<-c("Independent","Dependent")     #Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets.


##Explore the data##

#do a simple linear regression between variables
fit_H_FMw<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_H_FMw)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
H_vs_FMw_Data<-data

#get variables
Xvar<-H_vs_FMw_Data[,which(colnames(H_vs_FMw_Data)==paste("Independent")), drop=F]
Yvar<-H_vs_FMw_Data[,which(colnames(H_vs_FMw_Data)==paste("Dependent")), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)
model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML")
model
summary(model)

PGLS_HvsFMw_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_HvsFMw_Lambda

PGLS_HvsFMw_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_HvsFMw_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species.




##########################################################################
#########__Running phyANCOVA code with dataset minus Dinornis:__########## 
#############__hyperpallium volume vs. foramen magnum width__#############
##########################################################################



##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:26,28:30)
length(Others)
Dinornis<-c(27)
length(Dinornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Dinornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Dinornis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Llallawavis vs. all others

#assign group identities
Others<-c(1:6,8:30)
length(Others)
Llallawavis<-c(7)
length(Llallawavis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Llallawavis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Llallawavis]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Miocene_galliform vs. all others

#assign group identities
Others<-c(1:20,22:30)
length(Others)
Miocene_galliform<-c(21)
length(Miocene_galliform)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Miocene_galliform)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Miocene_galliform]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Paraptenodytes vs. all others

#assign group identities
Others<-c(1:12,14:30)
length(Others)
Paraptenodytes<-c(13)
length(Paraptenodytes)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Paraptenodytes)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Paraptenodytes]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant:  
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Psilopterus vs. all others

#assign group identities
Others<-c(1:5,7:30)
length(Others)
Psilopterus<-c(6)
length(Psilopterus)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Psilopterus)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Psilopterus]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)




####### Adjusting for rate of false discovery with Benjamini & Hochberg #######
####### hyperpallium and Wulst with all fossils #######

HvBR <- c(0.65,0.89,0.69,0.90,0.75)
WvER <- c(0.15,0.40,0.92,0.91,0.78)
WvFMa <- c(0.92,0.64,0.98,0.98,0.98)
WvFMw <- c(0.55,0.10,0.74,0.32,0.48)
HvFMa <- c(0.95,0.62,0.72,0.86,0.93)
HvFMw <- c(0.73,0.16,1.0,0.36,0.48)

adjHvBR <- p.adjust(HvBR, method="BH", n = length(HvBR))
write.csv(adjHvBR)

adjWvER <- p.adjust(WvER, method="BH", n = length(WvER))
write.csv(adjWvER)

adjWvFMa <- p.adjust(WvFMa, method="BH", n = length(WvFMa))
write.csv(adjWvFMa)

adjWvFMw <- p.adjust(WvFMw, method="BH", n = length(WvFMw))
write.csv(adjWvFMw)

adjHvFMa <- p.adjust(HvFMa, method="BH", n = length(HvFMa))
write.csv(adjHvFMa)

adjHvFMw <- p.adjust(HvFMw, method="BH", n = length(HvFMw))
write.csv(adjHvFMw)
