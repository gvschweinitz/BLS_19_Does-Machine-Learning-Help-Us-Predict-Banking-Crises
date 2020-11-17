# This script sets all parameters that are common across CV, IS and REC estimations
# To be called at the start of every code.
# Gregor von Schweinitz, 2020-11-17

############## Preamble ###############################################################

options(warn=1)
# options(warn=0)	 # default / factory setting: print list of warnings after execution
# options(warn=1)  # print warnings as they occur
# options(warn=2)  # treat warnings as errors
# 
# to show traceback, the following has to be activated:
#   Menu: Debug --> OnError --> Error Inspector


######################################
# PACKAGES AND R-CODES 
packages <- c("abind","data.table","doParallel","foreach","e1071","kknn","nnet","optimx","png","plyr","pROC","randomForest","reldist","rpart","rpart.plot","tictoc")
for (p in packages){
  check = require(p,character.only=TRUE)
  # INSTALLATION NEEDS TO BE DONE MANUALLY...
}

programs <- c("calc.train.bs.R","do.tables.R","estimate_cv.R","libraries_new.R","train_bs.R","train_bs_par.R")
for (p in programs){
  source(paste(wd_master,"/",p,sep=""))
}
rm(check,p,packages,programs)

#####################################
# LOAD DATASET AND CHOOSE VARIABLES
col.crisisdata <- "Crisis"  # set what name you expect for the crisis variable column (in the data)
col.crisis <- "c0"  # choose how you want the crisis column to be named (what our code expects)

# Select explanatory variables, grouped by model setup
usevars.list_in <- list(c("Var1","Var2","Var3"),
                        c("Var4", "Var5", "Var6", "Var7"),   
                        c("Var8", "Var9", "Var10")
                        )                     

usevars.list <- c(usevars.list_in,list(unique(unlist(usevars.list_in))))
cols.expl <- unique(unlist(usevars.list))                 # define columns with explanatory variables
usevars <- c("Country","Date", cols.expl,col.crisisdata)
ind.list = 1:length(usevars.list)

# Create dataset restricted to variables used in the following
inputData <-read.csv(paste(wd_master,"/",data_file,sep=""), header=TRUE, sep=",",dec=".",stringsAsFactors=FALSE)

# Controls for recursive window and for crossvalidation
recursivestart <- 2005.5    # choose start and end date for recursive exercise  ("2005Q3")
recursiveend <- 2016.5      # note: evaluation set could be empty if recursivestart > end of crisis database - 3 years

####################################
# Models to be estimated
methodz.ind <- c("logit", "trees", "knn", "rf", "svm", "nen")  
methodz.ind <- c("logit","trees")
combi <- expand.grid(ind.list,methodz.ind) # This is where the dataset selection "ind.list" enters into methodz.
methodz <- paste(combi[,2],combi[,1],sep=".")

####################################
# Output selection
doAUC <- TRUE
doBPS <- TRUE
doFmeas <- TRUE
doAUC.pvals <- TRUE
do.shift <- TRUE    # shift bootstrap distribution of Ur and BPS by their median

####################################
# Parameter settings
parameters = list()
# general
parameters$cols.expl <- usevars.list  # dataset will be selected later on, based on methodz
parameters$col.resp <- "pre"
parameters$col.crisis <- col.crisis
# parameters$methodz.agg = methodz # use all models in methodz
parameters$seed.value <- 1234567890
parameters$standardizeList = TRUE  # TRUE / FALSE. Standardize RHS variables before inputting to methods. 
parameters$drop.horizon = TRUE
parameters$horiz = c(5,12)
parameters$posthoriz = 0
parameters$svm.doProb <- TRUE
parameters$svmk <- "radial"

# evaluation criteria and output
parameters$optimizethreshold <- TRUE
parameters$u.alessidetken <- TRUE
parameters$mu <- 0.5
parameters$decimals <- 3

###################################
# More detailed data preparation
data.global = inputData[ , usevars]
colnames(data.global)[dim(data.global)[2]] <- col.crisis
data.global <- data.global[complete.cases(data.global[,1:(dim(data.global)[2]-1)]),]
data.global <- as.data.table(data.global,key=c("Country","Date"))
data.global <- calculate.pre(data.global, parameters$horiz, parameters$posthoriz) # adds pre, latepre, post
data.global$c0_NA<-0
data.global$c0_NA[is.na(data.global$c0)]=1   # add column indicating when crisis variable is NA
