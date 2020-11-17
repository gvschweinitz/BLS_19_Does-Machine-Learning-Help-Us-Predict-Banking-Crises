####################################
# Setting some names
datasets.eval <- "data.db.train"
resnames <- "res.cv"
posvecs.eval <- "pos.eval.cv"  # posvecs.eval is used to select from resnames (in script_tables)
resnames.bs <- NULL

####################################
# Bootstrap settings (turn off)
do.bootstrap <- FALSE
R.BS <- 1
blocklength <- 12           # blocklength for autocorrelation
doAUC.pvals <- FALSE        # make sure that there is no extra computation for AUC pvals
alphavec <- c(0.05,0.95)    # quantiles to be plotted as confidence bands  c(0.05,0.95)  c(0.1,0.9)

####################################
#  hyperparameter adaptations
parameters$optimizationMode <- TRUE
methodz <- methodz[!grepl("logit",methodz)]
# parameters$methodz.agg <- methodz

params <- setNames(as.list(rep(NA,length(methodz))),methodz)
P.method <- setNames(rep(0,length(methodz)),methodz)

### Parameter Grids
if (length(grep("knn",methodz))>0){
  para.vec <- seq(from=1, to=121, by =2)
  posm <- grep("knn",methodz)
  for (pm in posm){
    params[[methodz[pm]]] <- matrix(para.vec,nrow=length(para.vec),ncol=1,dimnames=list(c(1:length(para.vec)),"knnk"))
    P.method[methodz[pm]] <- length(para.vec)
  }
}
if (length(grep("trees",methodz))>0){
  posm <- grep("trees",methodz)
  para.vec = seq(from=0, to=0.3, length.out=100)
  for (pm in posm){
    params[[methodz[pm]]] <- matrix(para.vec,nrow=length(para.vec),ncol=1,dimnames=list(c(1:length(para.vec)),"treescp"))
    P.method[methodz[pm]] <- length(para.vec)
  }
}
if (length(grep("rf",methodz))>0){
  posm <- grep("rf",methodz)
  for (pm in posm){
    dset <- as.numeric(strsplit(methodz[pm],".",fixed=TRUE)[[1]][2])
    K <- length(parameters$cols.expl[[dset]])
    para.vec = expand.grid(nodesize=c(1:16),rfmtry=1:(K-1))
    params[[methodz[pm]]] <- para.vec
    P.method[methodz[pm]] <- dim(para.vec)[1]
  }
}
if (length(grep("svm",methodz))>0){
  posm <- grep("svm",methodz)
  z.svmg <- 2^(seq(-15,3,2))
  z.svmc <- c(seq(from=0.01, to=0.06, length.out=10),2^(seq(-3,9,2)))
  para.vec = expand.grid(svmg=z.svmg,svmc=z.svmc)
  for (pm in posm){
    params[[methodz[pm]]] <- para.vec
    P.method[methodz[pm]] <- dim(para.vec)[1]
  }
}
if (length(grep("nen",methodz))>0){
  posm <- grep("nen",methodz)
  z.nendecay <- 10^c(-8:1)
  for (pm in posm){
    dset <- as.numeric(strsplit(methodz[pm],".",fixed=TRUE)[[1]][2])
    K <- length(parameters$cols.expl[[dset]])
    # ub <- dim(data.db[rowSums(data.db[,c("latepre","post","c0_NA","c0")],na.rm=TRUE)==0,])[1]/10
    # ub <- floor((ub-1)/(K+2))
    z.nensize <- 1:K
    para.vec = expand.grid(nendecay=z.nendecay,nensize=z.nensize)
    params[[methodz[pm]]] <- para.vec
    P.method[methodz[pm]] <- dim(para.vec)[1]
  }
}
P <- max(P.method)

stat.select = "U_r" # optimize hyperparameters on: 1: U_r, 2: AUC

####################################
# Adaptation of data
print(paste("Running on in-sample data. Adapting data and setting parameters$drop.horizon = TRUE."))
data.global = data.global[Date < recursivestart,]
parameters$drop.horizon = FALSE

############## Parameter-Specific Processing Block ###############################################################
#   P = length(para.vec)
maxNumParameters = 3
statnames <- c("U_r",if(doAUC){"AUC"},if (doBPS){"BPS"})
resultnames <- c("p1","p2","p3", statnames)
num.stats = length(resultnames)  # parameter values, U_r, AUC

num.methods = length(methodz)

results.paras.array <- array(dim=c(P,num.stats,num.methods,length(resnames)),
                             dimnames=list(1:P,resultnames,methodz,resnames))
# results.paras.array: parameterization, stats, methodz, datasets
print(dimnames(results.paras.array))  
