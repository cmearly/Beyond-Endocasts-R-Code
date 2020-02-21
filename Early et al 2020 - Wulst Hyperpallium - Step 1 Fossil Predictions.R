##### Predicting Hyperpallium Volume of Fossils #####

library(ape)
library(caper)
library(car)
library(evomap)
library(geiger)
library(phytools)
require(phytools)
library(stats)


###!!!Don't forget to set your working directory!!!###




##########################################################################
#######################__Preparing for analysis__######################### 
##########################################################################



##Data Input##

#bring in dataset
data <- read.csv("wulst3.csv", row.names=1, header=T)

#make sure data table is in proper format
data<-as.data.frame(data)

#check how R classifies each column
sapply(data,class)

#check columns and first 6 rows of data
head(data)

#limit data to variables of interest
W_H.dat<-na.omit(data[,c(1,2,3)])
head(W_H.dat)
nrow(W_H.dat)

#log transform
W_H.dat[,2:3]<-log(W_H.dat[,2:3])


##Regression diagnostics##

#view relationships of all variables
scatterplotMatrix(data[,-1])

#simple linear regression between variables
fit_W_H<-lm(Hyperpallium_Vol~Wulst_SA, data=W_H.dat)

#test for outliers
outlierTest(fit_W_H)
qqPlot(fit_W_H) 
leveragePlots(fit_W_H) 

##test for non-normality
sresid_W_H<-studres(fit_W_H)
hist(sresid_W_H, freq=FALSE, breaks=7, main="distribution of studentized residuals, H vs W", col=rainbow(22))
xfit_W_H<-seq(min(sresid_W_H), max(sresid_W_H))
yfit_W_H<-dnorm(xfit_W_H)
lines(xfit_W_H, yfit_W_H, lwd=3, col="deeppink")

##check non-constant error variance / evaluate homoscedasticity
spreadLevelPlot(fit_W_H)

##check for nonlinearity
crPlots(fit_W_H)


##Tree manipulation##

#bring in phylogeny
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")
tree

#visually inspect tree
is.rooted(tree)
is.binary(tree)
is.ultrametric(tree)

#check that data and tree represent same taxa
name.check(W_H.dat,tree)

#prunes the tree to only include the members in the dataset and sorts it
sort.tree<-treedata(tree,W_H.dat,sort=TRUE,warnings=TRUE)$phy
sort.data<-treedata(tree,W_H.dat,sort=TRUE,warnings=TRUE)$data
sort.data<-as.data.frame(sort.data)
sort.tree
sort.tree<-force.ultrametric(sort.tree)
is.ultrametric(sort.tree)
plot(sort.tree,cex=0.5)

#make sure the tree has everyone you expected it to!

#export this new tree if you want to add fossils to it
write.tree(sort.tree, file="S_Wulsttree.tre", append=FALSE, digits=10, tree.names=FALSE)

#export the sorted dataset if you want to add fossils to it for phylogenetic prediction
write.csv(sort.data, "S_extantWulstlogged.csv")


##Check phylogenetic signal in each variable##

#first, make sure that variables are numeric
sapply(sort.data,class)

#make variables numeric
sort.data$Wulst_SA<-as.numeric(as.character(sort.data$Wulst_SA))
sort.data$Hyperpallium_Vol<-as.numeric(as.character(sort.data$Hyperpallium_Vol))

#check again to make sure conversion to numeric worked
sapply(sort.data,class)

#make Wulst a vector
W<-as.vector(sort.data)[,2]

#Pagel's lambda for Wulst
lambda_W<-phylosig(sort.tree,W,method="lambda", test=TRUE, nsim=9999)
lambda_W

#Blomberg's K for Wulst
BlomK_W<-phylosig(sort.tree,W,method="K", test=TRUE, nsim=9999)
BlomK_W

#make Wulst a vector
H<-as.vector(sort.data)[,3]

#Pagel's lambda for hyperpallium
lambda_H<-phylosig(sort.tree,H,method="lambda", test=TRUE, nsim=9999)
lambda_H

#Blomberg's K for hyperpallium
BlomK_H<-phylosig(sort.tree,H,method="K", test=TRUE, nsim=9999)
BlomK_H

#The values for lambda_H and BlomK_H will be drawn upon if you use BayesModelS
#to predict hyperpallium values for fossils




##########################################################################
#####################__Running a PGLS with caper__########################
##########################################################################


#check the data
head(sort.data)
row.names(sort.data)

#assign variables of interest
Y<-"Hyperpallium_Vol"
X<-"Wulst_SA"

#make dummy data object
H_vs_W_Data<-sort.data

#get variables
Yvar<-H_vs_W_Data[,which(colnames(H_vs_W_Data)==paste(Y)), drop=F]
Xvar<-H_vs_W_Data[,which(colnames(H_vs_W_Data)==paste(X)), drop=F]
data2<-cbind(Xvar, Yvar)
data2

#PGLS
data_pgls<-cbind(data2, rownames(data2))
tree_pgls<-sort.tree
colnames(data_pgls)[which(colnames(data_pgls)=="rownames(data2)")]<-"Species"
data.pgls <- comparative.data(tree_pgls, data_pgls, vcv=TRUE, names.col=Species)
model <- pgls(Yvar[,1] ~ Xvar[,1], data=data.pgls, lambda="ML")
model
summary(model)

#check residuals for phylogenetic signal
PGLS_HvsW_Lambda<-phylosig(sort.tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_HvsW_Lambda
PGLS_HvsW_BlomK<-phylosig(sort.tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_HvsW_BlomK




##########################################################################
#####__Predicting hyperpallium volume from Wulst surface area__######
##########################################################################



## !!! FIRST, run the code in "Early et al 2020 - Wulst Hyperpallium - Step 2 Fossil Merging.R" !!! ##

#load BayesModelS
source('BayesModelS_v24.R', chdir = TRUE)


##Dinornis##

#bring in the tree with just the extant sample of interest + Dinornis
treeDataAll <- read.tree("Dinornis_Wulst_block.tre")

#bring in the dataset of log-transformed Wulst surface areas and 
#hyperpallium volumes with Wulst surface area of Dinornis added
data1<-read.csv("Dinornis_Wulst_data.csv", header=TRUE)

#check to make sure data is in proper format
data1<-as.data.frame(data1)
colnames(data1)
data1

#identify factorName for BayesModelS
factorName<-c(data1$Specimen_No)

#identify the specimen for which you want to predict the missing value for BayesModelS
missingList<-c("Dinornis_robustus") 
missingList

#assign column names
colnames(data1)

#identify the formula for BayesModelS
formula<-"Hyperpallium_Vol~Wulst_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_H$lambda, kappaValue = BlomK_H$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Dinornis.H<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Dinornis
#hyperpallium volume
Dinornis.H.25q<-Dinornis.H[3]
Dinornis.H.Mean<-Dinornis.H[5]
Dinornis.H.75q<-Dinornis.H[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean hyperpallium volume for Dinornis and save it in the csv
#to be used for phyANCOVA
Dinornis.H.Mean


##Llallawavis##

#bring in the tree with just the extant sample of interest + Llallawavis
treeDataAll <- read.tree("Llallawavis_Wulst_block.tre")

#bring in the dataset of log-transformed Wulst surface areas and 
#hyperpallium volumes with Wulst surface area of Llallawavis added
data1<-read.csv("Llallawavis_Wulst_data.csv", header=TRUE)

#check to make sure data is in proper format
data1<-as.data.frame(data1)
colnames(data1)
data1

#identify factorName for BayesModelS
factorName<-c(data1$Specimen_No)

#identify the specimen for which you want to predict the missing value for BayesModelS
missingList<-c("Llallawavis_scagliai") 
missingList

#assign column names
colnames(data1)

#identify the formula for BayesModelS
formula<-"Hyperpallium_Vol~Wulst_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_H$lambda, kappaValue = BlomK_H$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Llallawavis.H<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Llallawavis
#hyperpallium volume
Llallawavis.H.25q<-Llallawavis.H[3]
Llallawavis.H.Mean<-Llallawavis.H[5]
Llallawavis.H.75q<-Llallawavis.H[6]

#this allows you to visualize MCMC convergence in a given run. Save this plH
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean hyperpallium volume for Llallawavis and save it in the csv
#to be used for phyANCOVA
Llallawavis.H.Mean


##Miocene galliform##

#bring in the tree with just the extant sample of interest + Miocene galliform
treeDataAll <- read.tree("Miocenegalliform_Wulst_block.tre")

#bring in the dataset of log-transformed Wulst surface areas and 
#hyperpallium volumes with Wulst surface area of Miocene galliform added
data1<-read.csv("Miocenegalliform_Wulst_data.csv", header=TRUE)

#check to make sure data is in proper format
data1<-as.data.frame(data1)
colnames(data1)
data1

#identify factorName for BayesModelS
factorName<-c(data1$Specimen_No)

#identify the specimen for which you want to predict the missing value for BayesModelS
missingList<-c("Miocene_galliform") 
missingList

#assign column names
colnames(data1)

#identify the formula for BayesModelS
formula<-"Hyperpallium_Vol~Wulst_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_H$lambda, kappaValue = BlomK_H$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Miocenegalliform.H<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Miocene galliform
#hyperpallium volume
Miocenegalliform.H.25q<-Miocenegalliform.H[3]
Miocenegalliform.H.Mean<-Miocenegalliform.H[5]
Miocenegalliform.H.75q<-Miocenegalliform.H[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean hyperpallium volume for Miocene galliform and save it in the csv
#to be used for phyANCOVA
Miocenegalliform.H.Mean


##Paraptenodytes##

#bring in the tree with just the extant sample of interest + Paraptenodytes
treeDataAll <- read.tree("Paraptenodytes_Wulst_block.tre")

#bring in the dataset of log-transformed Wulst surface areas and 
#hyperpallium volumes with Wulst surface area of Paraptenodytes added
data1<-read.csv("Paraptenodytes_Wulst_data.csv", header=TRUE)

#check to make sure data is in proper format
data1<-as.data.frame(data1)
colnames(data1)
data1

#identify factorName for BayesModelS
factorName<-c(data1$Specimen_No)

#identify the specimen for which you want to predict the missing value for BayesModelS
missingList<-c("Paraptenodytes_antarcticus") 
missingList

#assign column names
colnames(data1)

#identify the formula for BayesModelS
formula<-"Hyperpallium_Vol~Wulst_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_H$lambda, kappaValue = BlomK_H$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Paraptenodytes.H<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Paraptenodytes
#hyperpallium volume
Paraptenodytes.H.25q<-Paraptenodytes.H[3]
Paraptenodytes.H.Mean<-Paraptenodytes.H[5]
Paraptenodytes.H.75q<-Paraptenodytes.H[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean hyperpallium volume for Paraptenodytes and save it in the csv
#to be used for phyANCOVA
Paraptenodytes.H.Mean


##Psilopterus##

#bring in the tree with just the extant sample of interest + Psilopterus
treeDataAll <- read.tree("Psilopterus_Wulst_block.tre")

#bring in the dataset of log-transformed Wulst surface areas and 
#hyperpallium volumes with Wulst surface area of Psilopterus added
data1<-read.csv("Psilopterus_Wulst_data.csv", header=TRUE)

#check to make sure data is in proper format
data1<-as.data.frame(data1)
colnames(data1)
data1

#identify factorName for BayesModelS
factorName<-c(data1$Specimen_No)

#identify the specimen for which you want to predict the missing value for BayesModelS
missingList<-c("Psilopterus_lemoinei") 
missingList

#assign column names
colnames(data1)

#identify the formula for BayesModelS
formula<-"Hyperpallium_Vol~Wulst_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_H$lambda, kappaValue = BlomK_H$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Psilopterus.H<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Psilopterus
#hyperpallium volume
Psilopterus.H.25q<-Psilopterus.H[3]
Psilopterus.H.Mean<-Psilopterus.H[5]
Psilopterus.H.75q<-Psilopterus.H[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean hyperpallium volume for Psilopterus and save it in the csv
#to be used for phyANCOVA
Psilopterus.H.Mean

