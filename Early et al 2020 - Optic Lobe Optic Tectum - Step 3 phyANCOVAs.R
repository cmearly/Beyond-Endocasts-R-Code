##### Phylogenetic ANCOVAs for Optic Tectum and Optic Lobe of Fossils #####

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
##############__optic tectum volume vs. brain rest volume__###############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModelS
#code run in "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R"
data<-read.csv("alloptictectum.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"brain_rest_vol"; Y<-"tectum_vol"

#limit the dataset to your X and Y variables of interest
data<-data[,c(which(colnames(data)==X),which(colnames(data)==Y)),drop=F]
data<-na.omit(data) 

#log transform the data
data<-log(data)

#bring in the tree file
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")

#sort the tree and the data for each other
tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-treedata(tree,data,sort=T,warnings=T)$data   #match the data to the tree

#make the sorted data a dataframe
data<-as.data.frame(data)

#if desired, change the column names of the data for repeatability
colnames(data)<-c("Independent","Dependent")     

#check the tree
plot(tree, cex=0.5)


##Explore the data##

#do a simple linear regression between variables
fit_OT_B<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_OT_B)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
OT_vs_B_Data<-data

#get variables
Xvar<-OT_vs_B_Data[,which(colnames(OT_vs_B_Data)==paste("Independent")), drop=F]
Yvar<-OT_vs_B_Data[,which(colnames(OT_vs_B_Data)==paste("Dependent")), drop=F]
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

#check phylogenetic signal in regression

PGLS_OTvsBR_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_OTvsBR_Lambda

PGLS_OTvsBR_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_OTvsBR_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to use to assign group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species



##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset :__#######
##############__optic tectum volume vs. brain rest volume__###############
##########################################################################


##phyANCOVA of Archaeopteryx vs. all others 

#assign group identities
Others<-c(1:177)
length(Others)
Archaeopteryx<-c(178)
length(Archaeopteryx)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Archaeopteryx)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Archaeopteryx]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:173,175:178)
length(Others)
Dinornis<-c(174)
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


##phyANCOVA of Lithornis vs. all others

#assign group identities
Others<-c(1:170,172:178)
length(Others)
Lithornis<-c(171)
length(Lithornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Lithornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Lithornis]<-"B" 
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
Others<-c(1:73,75:178)
length(Others)
Llallawavis<-c(74)
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


##phyANCOVA of Miocene galliform vs. all others

#assign group identities
Others<-c(1:150,152:178)
length(Others)
Miocene_galliform<-c(151)
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
Others<-c(1:100,102:178)
length(Others)
Paraptenodytes<-c(101)
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
Others<-c(1:72,74:178)
length(Others)
Psilopterus<-c(73)
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
########__optic lobe surface area vs. endocast rest surface area__########
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModelS
#code done in "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R", 
#as well as endocast and foramen magnum values
data<-read.csv("allopticlobesandtectum.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"endo_rest_SA"; Y<-"optic_lobes_SA"

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
colnames(data)<-c("Independent","Dependent")     


##Explore the data##

#do a simple linear regression between variables
fit_OL_E<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_OL_E)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
OL_vs_E_Data<-data

#get variables
Xvar<-OL_vs_E_Data[,which(colnames(OL_vs_E_Data)==paste("Independent")), drop=F]
Yvar<-OL_vs_E_Data[,which(colnames(OL_vs_E_Data)==paste("Dependent")), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)

#make optic lobe a vector
OL<-as.vector(data)[,2]

#Pagel's lambda for optic lobe
lambda_OL<-phylosig(tree,OL,method="lambda", test=TRUE, nsim=9999)
lambda_OL

#make endo rest a vector
ER<-as.vector(data)[,2]

#Pagel's lambda for endo rest
lambda_ER<-phylosig(tree,ER,method="lambda", test=TRUE, nsim=9999)
lambda_ER

model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML", bounds=list(lambda=c(0.001,1)))
model
summary(model)

#check phylogenetic signal in regression

PGLS_OLvsER_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_OLvsER_Lambda

PGLS_OLvsER_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_OLvsER_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])


#Now we can begin to run the phyANCOVA code for all the extinct + extant species




##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__######## 
########__optic lobe surface area vs. endocast rest surface area__########
##########################################################################



##phyANCOVA of Archaeopteryx vs. all others 

#assign group identities
Others<-c(1:40)
length(Others)
Archaeopteryx<-c(41)
length(Archaeopteryx)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Archaeopteryx)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Archaeopteryx]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:36,38:41)
length(Others)
Dinornis<-c(37)
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


##phyANCOVA of Lithornis vs. all others

#assign group identities
Others<-c(1:34,36:41)
length(Others)
Lithornis<-c(35)
length(Lithornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Lithornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Lithornis]<-"B" 
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
Others<-c(1:6,8:41)
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


##phyANCOVA of Miocene galliform vs. all others

#assign group identities
Others<-c(1:28,30:41)
length(Others)
Miocene_galliform<-c(29)
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
Others<-c(1:17,19:41)
length(Others)
Paraptenodytes<-c(18)
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
Others<-c(1:5,7:41)
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
###########__optic lobe surface area vs. foramen magnum area__############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModels
#code done in "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allopticlobesandtectum.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_SA"; Y<-"optic_lobes_SA"

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
fit_OL_FM<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_OL_FM)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
OL_vs_FM_Data<-data

#get variables
Xvar<-OL_vs_FM_Data[,which(colnames(OL_vs_FM_Data)==paste("Independent")), drop=F]
Yvar<-OL_vs_FM_Data[,which(colnames(OL_vs_FM_Data)==paste("Dependent")), drop=F]
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

#check phylogenetic signal of regression

PGLS_OLvsFM_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_OLvsFM_Lambda

PGLS_OLvsFM_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_OLvsFM_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])





##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__########
###########__optic lobe surface area vs. foramen magnum area__############
##########################################################################



##phyANCOVA of Archaeopteryx vs. all others 

#assign group identities
Others<-c(1:40)
length(Others)
Archaeopteryx<-c(41)
length(Archaeopteryx)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Archaeopteryx)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Archaeopteryx]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:36,38:41)
length(Others)
Dinornis<-c(37)
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


##phyANCOVA of Lithornis vs. all others

#assign group identities
Others<-c(1:34,36:41)
length(Others)
Lithornis<-c(35)
length(Lithornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Lithornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Lithornis]<-"B" 
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
Others<-c(1:6,8:41)
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


##phyANCOVA of Miocene galliform vs. all others

#assign group identities
Others<-c(1:28,30:41)
length(Others)
Miocene_galliform<-c(29)
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
Others<-c(1:17,19:41)
length(Others)
Paraptenodytes<-c(18)
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
Others<-c(1:5,7:39)
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
##########__optic lobe surface area vs. foramen magnum width__############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModels
#code done in "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allopticlobesandtectum.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_w"; Y<-"optic_lobes_SA"

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
fit_OL_FMw<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_OL_FMw)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
OL_vs_FMw_Data<-data

#get variables
Xvar<-OL_vs_FMw_Data[,which(colnames(OL_vs_FMw_Data)==paste("Independent")), drop=F]
Yvar<-OL_vs_FMw_Data[,which(colnames(OL_vs_FMw_Data)==paste("Dependent")), drop=F]
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
slope_max<-summary(model)$coefficients[2,1]+((qt(0.975,summary(model)$df[2]))*summary(model)$coefficients[2,2])
slope_min<-summary(model)$coefficients[2,1]-((qt(0.975,summary(model)$df[2]))*summary(model)$coefficients[2,2])
slope_max
slope_min

#check phylogenetic signal of regression

PGLS_OLvsFMw_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_OLvsFMw_Lambda

PGLS_OLvsFMw_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_OLvsFMw_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])






##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__########
##########__optic lobe surface area vs. foramen magnum width__############
##########################################################################



##phyANCOVA of Archaeopteryx vs. all others

#assign group identities
Others<-c(1:40)
length(Others)
Archaeopteryx<-c(41)
length(Archaeopteryx)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Archaeopteryx)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Archaeopteryx]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:36,38:41)
length(Others)
Dinornis<-c(37)
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


##phyANCOVA of Lithornis vs. all others

#assign group identities
Others<-c(1:34,36:39)
length(Others)
Lithornis<-c(35)
length(Lithornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Lithornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Lithornis]<-"B" 
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
Others<-c(1:6,8:39)
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


##phyANCOVA of Miocene galliform vs. all others

#assign group identities
Others<-c(1:28,30:39)
length(Others)
Miocene_galliform<-c(29)
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
Others<-c(1:17,19:39)
length(Others)
Paraptenodytes<-c(18)
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
Others<-c(1:5,7:39)
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
#############__optic tectum volume vs. foramen magnum area__##############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModels
#code done in "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allopticlobesandtectum.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_SA"; Y<-"optic_tectum_vol"

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
fit_OT_FM<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_OT_FM)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
OT_vs_FM_Data<-data

#get variables
Xvar<-OT_vs_FM_Data[,which(colnames(OT_vs_FM_Data)==paste("Independent")), drop=F]
Yvar<-OT_vs_FM_Data[,which(colnames(OT_vs_FM_Data)==paste("Dependent")), drop=F]
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

#check phylogenetic signal of regressions

PGLS_OTvsFM_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_OTvsFM_Lambda

PGLS_OTvsFM_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_OTvsFM_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))
#copy this output into Excel to search for assigning group identities

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])






##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__######## 
#############__optic tectum volume vs. foramen magnum area__##############
##########################################################################



##phyANCOVA of Archaeopteryx vs. all others

#assign group identities
Others<-c(1:40)
length(Others)
Archaeopteryx<-c(41)
length(Archaeopteryx)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Archaeopteryx)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Archaeopteryx]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:36,38:41)
length(Others)
Dinornis<-c(37)
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


##phyANCOVA of Lithornis vs. all others

#assign group identities
Others<-c(1:34,36:41)
length(Others)
Lithornis<-c(35)
length(Lithornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Lithornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Lithornis]<-"B" 
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
Others<-c(1:6,8:41)
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


##phyANCOVA of Miocene galliform vs. all others

#assign group identities
Others<-c(1:28,30:41)
length(Others)
Miocene_galliform<-c(29)
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
Others<-c(1:17,19:41)
length(Others)
Paraptenodytes<-c(18)
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
Others<-c(1:5,7:41)
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
#############__optic tectum volume vs. foramen magnum width__#############
##########################################################################



##Data and tree preparation##

#bring in the dataset; this will be brain rest volumes and optic tectum
#volumes, including the predicted optic tectum volumes from the BayesModels
#code done in "Early et al 2020 - Optic Lobe Optic Tectum - Step 1 Fossil Predictions.R", as well
#as endocast and foramen magnum values
data<-read.csv("allopticlobesandtectum.csv", row.names=1, header=TRUE)

#assign your X and Y variables
X<-"FM_w"; Y<-"optic_tectum_vol"

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
colnames(data)<-c("Independent","Dependent")    


##Explore the data##

#do a simple linear regression between variables
fit_OT_FMw<-lm(Dependent~Independent, data=data)

#test for outliers
outlierTest(fit_OT_FMw)


##Run a PGLS##

head(data)
row.names(data)

#make dummy data object
OT_vs_FMw_Data<-data

#get variables
Xvar<-OT_vs_FMw_Data[,which(colnames(OT_vs_FMw_Data)==paste("Independent")), drop=F]
Yvar<-OT_vs_FMw_Data[,which(colnames(OT_vs_FMw_Data)==paste("Dependent")), drop=F]
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

#check phylogenetic signal of regression

PGLS_OTvsFMw_Lambda<-phylosig(tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_OTvsFMw_Lambda

PGLS_OTvsFMw_BlomK<-phylosig(tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_OTvsFMw_BlomK


##Figure out tip labels for tree to assign group identities

#assign n to equal the number of tips in the tree
n<-length(tree$tip.label)

#instead, run this code to get a numbered list of the tip labels  
write.csv(as.matrix(tree$tip.label))

#rescale the tree to the PGLS you just ran
tree<-rescale(tree,"lambda",model$param[[2]])






##########################################################################
#######__Running phyANCOVA code with extant and extinct dataset:__######## 
#############__optic tectum volume vs. foramen magnum width__#############
##########################################################################



##phyANCOVA of Archaeopteryx vs. all others

#assign group identities
Others<-c(1:40)
length(Others)
Archaeopteryx<-c(41)
length(Archaeopteryx)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Archaeopteryx)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Archaeopteryx]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

#set up model
Model<-model.matrix(as.formula(Dependent~Independent),data)

#set up testing for differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#test for differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)


##phyANCOVA of Dinornis vs. all others

#assign group identities
Others<-c(1:36,38:41)
length(Others)
Dinornis<-c(37)
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


##phyANCOVA of Lithornis vs. all others

#assign group identities
Others<-c(1:34,36:41)
length(Others)
Lithornis<-c(35)
length(Lithornis)

#check to make sure your group identities add up to the number of tips
n == ((length(Others)) + (length(Lithornis)))

#set up the groups for testing in differences in intercept holding slope constant
grpI<-rep("A",length(rownames(data))) 
grpI[Lithornis]<-"B" 
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
Others<-c(1:6,8:41)
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


##phyANCOVA of Miocene galliform vs. all others

#assign group identities
Others<-c(1:28,30:41)
length(Others)
Miocene_galliform<-c(29)
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
Others<-c(1:17,19:41)
length(Others)
Paraptenodytes<-c(18)
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
Others<-c(1:5,7:41)
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




####### Adjusting for rate of false discovery with Benjamini & Hochberg procedure #######
####### optic tectum and lobe with all fossils #######

#make lists of the p-values of phyANCOVAs of all fossils and keep them in the same order
#for each combination of variables, then bring them into R
#for example, these are the p-values for Archaeopteryx, Lithornis, Dinornis, 
#Miocene galliform, Paraptenodytes, Psilopterus, and Llallawavis, in that order
OTvBR <- c(0.61,0.086,0.0005,0.16,0.95,0.79,0.82)
OLvER <- c(0.38,0.36,0.016,0.90,0.57,0.46,0.37)
OLvFMa <- c(0.78,0.99,0.012,0.83,0.63,0.43,0.39)
OLvFMw <- c(0.22,0.69,0.022,0.73,0.60,0.22,0.42)
OTvFMa <- c(0.98,0.97,0.010,0.79,0.59,0.59,0.58)
OTvFMw <- c(0.25,0.70,0.027,0.71,0.55,0.33,0.59)

#perform the adjustment, then copy the list of the adjusted p-values for reporting

adjOTvBR <- p.adjust(OTvBR, method="BH", n = length(OTvBR))
write.csv(adjOTvBR)

adjOLvER <- p.adjust(OLvER, method="BH", n = length(OLvER))
write.csv(adjOLvER)

adjOLvFMa <- p.adjust(OLvFMa, method="BH", n = length(OLvFMa))
write.csv(adjOLvFMa)

adjOLvFMw <- p.adjust(OLvFMw, method="BH", n = length(OLvFMw))
write.csv(adjOLvFMw)

adjOTvFMa <- p.adjust(OTvFMa, method="BH", n = length(OTvFMa))
write.csv(adjOTvFMa)

adjOTvFMw <- p.adjust(OTvFMw, method="BH", n = length(OTvFMw))
write.csv(adjOTvFMw)


