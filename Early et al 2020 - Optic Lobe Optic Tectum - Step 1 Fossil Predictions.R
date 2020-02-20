##### Predicting Optic Tectum Volume of Fossils #####

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
data <- read.csv("opticlobes4.csv", row.names=1, header=T)

#make sure data table is in proper format
data<-as.data.frame(data)

#check how R classifies each column
sapply(data,class)

#check columns and first 6 rows of data
head(data)

#limit data to variables of interest
OL_OT.dat<-na.omit(data[,c(1,2,3)])
head(OL_OT.dat)
nrow(OL_OT.dat)

#log transform
OL_OT.dat[,2:3]<-log(OL_OT.dat[,2:3])


##Regression diagnostics##

#view relationships of all variables
scatterplotMatrix(data[,-1])

#simple linear regression between variables
fit_OL_OT<-lm(optic_tectum_Vol~optic_lobes_SA, data=OL_OT.dat)

#test for outliers
outlierTest(fit_OL_OT)
qqPlot(fit_OL_OT)
leveragePlots(fit_OL_OT)

#test for non-normality
sresid_OL_OT<-studres(fit_OL_OT)
hist(sresid_OL_OT, freq=FALSE, breaks=7, main="distribution of studentized residuals, OT vs OL", col=rainbow(22))
xfit_OL_OT<-seq(min(sresid_OL_OT), max(sresid_OL_OT))
yfit_OL_OT<-dnorm(xfit_OL_OT)
lines(xfit_OL_OT, yfit_OL_OT, lwd=3, col="deeppink")

#check non-constant error variance / evaluate homoscedasticity
spreadLevelPlot(fit_OL_OT)

#check for nonlinearity
crPlots(fit_OL_OT)


##Tree manipulation##

#bring in phylogeny
tree<-read.nexus("alloptictectumandhyperpalliumtaxa.tre")
tree

#inspect tree
is.rooted(tree)
is.binary(tree)
is.ultrametric(tree)

#check that data and tree represent same taxa
name.check(OL_OT.dat,tree)

#prune the tree to only include the members in the dataset and sort it
sort.tree<-treedata(tree,OL_OT.dat,sort=TRUE,warnings=TRUE)$phy
sort.data<-treedata(tree,OL_OT.dat,sort=TRUE,warnings=TRUE)$data
sort.data<-as.data.frame(sort.data)
sort.tree
sort.tree<-force.ultrametric(sort.tree)
is.ultrametric(sort.tree)
plot(sort.tree,cex=0.5)

#make sure the tree has everyone you expected it to!

#export this new tree if you ever want to add fossils to it
write.tree(sort.tree, file="S_optictree.tre", append=FALSE, digits=10, tree.names=FALSE)

#export the sorted dataset if you want to add fossils to it for phylogenetic prediction
write.csv(sort.data, "S_extantopticlogged.csv")


##Check phylogenetic signal in each variable##

##first, make sure that variables are numeric
sapply(sort.data,class)

#make variables numeric
sort.data$optic_lobes_SA<-as.numeric(as.character(sort.data$optic_lobes_SA))
sort.data$optic_tectum_Vol<-as.numeric(as.character(sort.data$optic_tectum_Vol))

#check again to make sure conversion to numeric worked
sapply(sort.data,class)

#make optic lobe a vector
OL<-as.vector(sort.data)[,2]

#Pagel's lambda for optic lobe
lambda_OL<-phylosig(sort.tree,OL,method="lambda", test=TRUE, nsim=9999)
lambda_OL

#Blomberg's K for optic lobe
BlomK_OL<-phylosig(sort.tree,OL,method="K", test=TRUE, nsim=9999)
BlomK_OL

#make optic tectum a vector
OT<-as.vector(sort.data)[,3]

#Pagel's lambda for optic tectum
lambda_OT<-phylosig(sort.tree,OT,method="lambda", test=TRUE, nsim=9999)
lambda_OT

#Blomberg's K for optic tectum
BlomK_OT<-phylosig(sort.tree,OT,method="K", test=TRUE, nsim=9999)
BlomK_OT

#The values for lambda_OT and BlomK_OT will be drawn upon if you use BayesModelS
#to predict optic tectum values for fossils




##########################################################################
#####################__Running a PGLS with caper__########################
##########################################################################


#check the data
head(sort.data)
row.names(sort.data)

#assign variables of interest
Y<-"optic_tectum_Vol"
X<-"optic_lobes_SA"

#make dummy data object
OT_vs_OL_Data<-sort.data

#get variables
Yvar<-OT_vs_OL_Data[,which(colnames(OT_vs_OL_Data)==paste(Y)), drop=F]
Xvar<-OT_vs_OL_Data[,which(colnames(OT_vs_OL_Data)==paste(X)), drop=F]
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
PGLS_OTvsOL_Lambda<-phylosig(sort.tree, model$residuals, method="lambda", test=TRUE, nsim=999)
PGLS_OTvsOL_Lambda
PGLS_OTvsOL_BlomK<-phylosig(sort.tree, model$residuals, method="K", test=TRUE, nsim=999)
PGLS_OTvsOL_BlomK




##########################################################################
#####__Predicting optic tectum volume from optic lobe surface area__######
##########################################################################



## !!! FIRST, run the code in "Early et al 2020 - Optic Lobe Optic Tectum - Step 2 Fossil Merging.R" !!! ##

#load BayesModelS
source('BayesModelS_v24.R', chdir = TRUE)


##Archaeopteryx##

#bring in the tree with just the extant sample of interest + Archaeopteryx
treeDataAll <- read.tree("Archaeopteryx_block.tre")

#bring in the dataset of log-transformed optic lobe surface areas and 
#optic tectum volumes with optic lobe surface area of Archaeopteryx added
data1<-read.csv("Archaeopteryx_data.csv", header=TRUE)

#check to make sure data is in proper format
data1<-as.data.frame(data1)
colnames(data1)
data1

#identify factorName for BayesModelS
factorName<-c(data1$Specimen_No)

#identify the specimen for which you want to predict the missing value for BayesModelS
missingList<-c("Archaeopteryx_lithographica") 
missingList

#assign column names
colnames(data1)

#identify the formula for BayesModelS
formula<-"optic_tectum_Vol~optic_lobes_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_OT$lambda, kappaValue = BlomK_OT$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Archaeopteryx.OT<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Archaeopteryx
#optic tectum volume
Archaeopteryx.OT.25q<-Archaeopteryx.OT[3]
Archaeopteryx.OT.Mean<-Archaeopteryx.OT[5]
Archaeopteryx.OT.75q<-Archaeopteryx.OT[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean optic tectum volume for Archaeopteryx and save it in the csv
#to be used for phyANCOVA
Archaeopteryx.OT.Mean


##Dinornis##

#bring in the tree with just the extant sample of interest + Dinornis
treeDataAll <- read.tree("Dinornis_block.tre")

#bring in the dataset of log-transformed optic lobe surface areas and 
#optic tectum volumes with optic lobe surface area of Dinornis added
data1<-read.csv("Dinornis_data.csv", header=TRUE)

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
formula<-"optic_tectum_Vol~optic_lobes_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_OT$lambda, kappaValue = BlomK_OT$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Dinornis.OT<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Dinornis
#optic tectum volume
Dinornis.OT.25q<-Dinornis.OT[3]
Dinornis.OT.Mean<-Dinornis.OT[5]
Dinornis.OT.75q<-Dinornis.OT[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean optic tectum volume for Dinornis and save it in the csv
#to be used for phyANCOVA
Dinornis.OT.Mean


##Lithornis##

#bring in the tree with just the extant sample of interest + Lithornis
treeDataAll <- read.tree("Lithornis_block.tre")

#bring in the dataset of log-transformed optic lobe surface areas and 
#optic tectum volumes with optic lobe surface area of Lithornis added
data1<-read.csv("Lithornis_data.csv", header=TRUE)

#check to make sure data is in proper format
data1<-as.data.frame(data1)
colnames(data1)
data1

#identify factorName for BayesModelS
factorName<-c(data1$Specimen_No)

#identify the specimen for which you want to predict the missing value for BayesModelS
missingList<-c("Lithornis_plebius") 
missingList

#assign column names
colnames(data1)

#identify the formula for BayesModelS
formula<-"optic_tectum_Vol~optic_lobes_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_OT$lambda, kappaValue = BlomK_OT$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Lithornis.OT<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Lithornis
#optic tectum volume
Lithornis.OT.25q<-Lithornis.OT[3]
Lithornis.OT.Mean<-Lithornis.OT[5]
Lithornis.OT.75q<-Lithornis.OT[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean optic tectum volume for Lithornis and save it in the csv
#to be used for phyANCOVA
Lithornis.OT.Mean


##Llallawavis##

#bring in the tree with just the extant sample of interest + Llallawavis
treeDataAll <- read.tree("Llallawavis_block.tre")

#bring in the dataset of log-transformed optic lobe surface areas and 
#optic tectum volumes with optic lobe surface area of Llallawavis added
data1<-read.csv("Llallawavis_data.csv", header=TRUE)

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
formula<-"optic_tectum_Vol~optic_lobes_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_OT$lambda, kappaValue = BlomK_OT$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Llallawavis.OT<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Llallawavis
#optic tectum volume
Llallawavis.OT.25q<-Llallawavis.OT[3]
Llallawavis.OT.Mean<-Llallawavis.OT[5]
Llallawavis.OT.75q<-Llallawavis.OT[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean optic tectum volume for Llallawavis and save it in the csv
#to be used for phyANCOVA
Llallawavis.OT.Mean


##Miocene galliform##

#bring in the tree with just the extant sample of interest + Miocene galliform
treeDataAll <- read.tree("Miocenegalliform_block.tre")

#bring in the dataset of log-transformed optic lobe surface areas and 
#optic tectum volumes with optic lobe surface area of Miocene galliform added
data1<-read.csv("Miocenegalliform_data.csv", header=TRUE)

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
formula<-"optic_tectum_Vol~optic_lobes_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_OT$lambda, kappaValue = BlomK_OT$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Miocenegalliform.OT<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Miocene galliform
#optic tectum volume
Miocenegalliform.OT.25q<-Miocenegalliform.OT[3]
Miocenegalliform.OT.Mean<-Miocenegalliform.OT[5]
Miocenegalliform.OT.75q<-Miocenegalliform.OT[6]

#this allows you to visualize MCMC convergence in a given run
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean optic tectum volume for the Miocene galliform and save it in the csv
#to be used for phyANCOVA
Miocenegalliform.OT.Mean


##Paraptenodytes##

#bring in the tree with just the extant sample of interest + Paraptenodytes
treeDataAll <- read.tree("Paraptenodytes_block.tre")

#bring in the dataset of log-transformed optic lobe surface areas and 
#optic tectum volumes with optic lobe surface area of Paraptenodytes added
data1<-read.csv("Paraptenodytes_data.csv", header=TRUE)

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
formula<-"optic_tectum_Vol~optic_lobes_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_OT$lambda, kappaValue = BlomK_OT$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Paraptenodytes.OT<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Paraptenodytes
#optic tectum volume
Paraptenodytes.OT.25q<-Paraptenodytes.OT[3]
Paraptenodytes.OT.Mean<-Paraptenodytes.OT[5]
Paraptenodytes.OT.75q<-Paraptenodytes.OT[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean optic tectum volume for Paraptenodytes and save it in the csv
#to be used for phyANCOVA
Paraptenodytes.OT.Mean


##Psilopterus##

#bring in the tree with just the extant sample of interest + Psilopterus
treeDataAll <- read.tree("Psilopterus_block.tre")

#bring in the dataset of log-transformed optic lobe surface areas and 
#optic tectum volumes with optic lobe surface area of Psilopterus added
data1<-read.csv("Psilopterus_data.csv", header=TRUE)

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
formula<-"optic_tectum_Vol~optic_lobes_SA"

#run BayesModelS! This will take a while to run through. The code will output
#various diagnostics that will indicate if you had set it up to run for enough
#iterations to yield robust results.
bmselection = blm(formula, data1, treeDataAll, factorName=factorName, missingList = missingList, currentValue = 0, nposterior = 200000, burnin = 20000, thin = 50, varSelection = "lambda", lambdaValue = lambda_OT$lambda, kappaValue = BlomK_OT$K)
modelChecking(bmselection, missingList)
analysis(bmselection)
predict(bmselection, missingList)
exp(predict(bmselection, missingList))
Psilopterus.OT<-predict(bmselection, missingList)

#get R to store the 25th quartile, mean, and 75th quartile values for Psilopterus
#optic tectum volume
Psilopterus.OT.25q<-Psilopterus.OT[3]
Psilopterus.OT.Mean<-Psilopterus.OT[5]
Psilopterus.OT.75q<-Psilopterus.OT[6]

#this allows you to visualize MCMC convergence in a given run. Save this plot
#for supplementals.
outputPosterior = read.csv("result.csv")
plot(outputPosterior$lkhoodSample, xlab="Iteration", ylab="Likelihood", cex=0.8, pch=1)

#call the mean optic tectum volume for Psilopterus and save it in the csv
#to be used for phyANCOVA
Psilopterus.OT.Mean

