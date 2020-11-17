# This script runs the full cross-validation exercise.
# The main output ("outtable_05.csv") is not saved in the main working directory, but in the output working directory
# Note that this csv-file needs to be copied into the main working directory for further use in in-sample and recursive estimations
# Gregor von Schweinitz, 2020-11-17

sink.N = sink.number()
if (sink.N>0){for (j in 1 : sink.N){sink()}} # close existing sinks

######################################
# Set working directory and input files
wd <- getwd()
wd_master <- wd
output_directory = "/baseline_mu05"
wd_output_detail <- paste(wd_master,output_directory,sep="")
if (!dir.exists(wd_output_detail)){dir.create(wd_output_detail,recursive=TRUE)}

data_file <- "data_random.csv"                  # data file with all variables
paramtable_file <- "outtable_05.csv"
filename <- "Data_master"     # Output filename for workspace of script_master
filename <- paste(wd_output_detail,"/",filename,sep="")

######################################
# Set up specifictions of estimation. Switch here between CV, IS and REC estimation
do.is <- FALSE
do.recursive <- FALSE
do.crossval <- TRUE
if (sum(c(do.is,do.recursive,do.crossval))!=1){stop("can only run one exercise at a time")}

source("script_prep_gen.R")
tic("start.master")
if (do.is){source("script_prep_is.R")}
if (do.recursive){source("script_prep_recursive.R")}
if (do.crossval){source("script_prep_cv.R")}

do.parallel <- TRUE
if (do.parallel){
  cl <- makeCluster(12,type="PSOCK")
  registerDoParallel(cl)
}
set.seed(parameters$seed.value)

######################################
# Run estimation.
if (do.is | do.recursive){
  source("script_estimate.R")
  source("script_tables.R")
}
if (do.crossval){
  methodz.results.all <- methodz
  results.paras.array <- array(dim=c(P,num.stats,num.methods,length(resnames)),
                               dimnames=list(1:P,resultnames,methodz.results.all,resnames))
  res.optimize <- "res.cv"
  for (p in 1:P){
    tic("p-run")
    # THE FOLLOWING IS CORRECT AS LONG AS params ARE SORTED SUCH THAT EVERY METHOD HAS IDENTICAL PARAMS IN IDENTICAL POSITIONS, INDEPENDENT OF DATASET
    methodz <- names(P.method)[P.method>=p]
    num.methods <- length(methodz)
    print(paste("Starting Iteration ", p, " of ",P, " on ",paste(methodz,collapse=", "), "."), sep="")
    if (length(grep("knn",methodz))>0){
      parameters$knnk = params[[methodz[grep("knn",methodz)[1]]]][p,"knnk"]
      print(paste("knn:",params[[methodz[grep("knn",methodz)[1]]]][p,]))
    }
    if (length(grep("trees",methodz))>0) {
      parameters$treescp = params[[methodz[grep("trees",methodz)[1]]]][p,"treescp"]
      print(paste("trees:",params[[methodz[grep("trees",methodz)[1]]]][p,]))
    }
    if (length(grep("rf",methodz))>0){
      parameters$rfmtry = params[[methodz[grep("rf",methodz)[1]]]][p,"rfmtry"]
      parameters$nodesize <- params[[methodz[grep("rf",methodz)[1]]]][p,"nodesize"]
      print(paste("rf:",params[[methodz[grep("rf",methodz)[1]]]][p,]))
    }
    if (length(grep("svm",methodz))>0){
      parameters$svmg <- params[[methodz[grep("svm",methodz)[1]]]][p,"svmg"]
      parameters$svmc <- params[[methodz[grep("svm",methodz)[1]]]][p,"svmc"]
      print(paste("svm:",params[[methodz[grep("svm",methodz)[1]]]][p,]))
    }
    if (length(grep("nen",methodz))>0){
      parameters$nensize <- params[[methodz[grep("nen",methodz)[1]]]][p,"nensize"]
      parameters$nendecay <- params[[methodz[grep("nen",methodz)[1]]]][p,"nendecay"]
      print(paste("nen:",params[[methodz[grep("nen",methodz)[1]]]][p,]))
    }
    if (do.parallel){
      acomb <- function(...) abind(...,along=3)
      
      packages.export <- c("data.table","e1071","kknn","nnet","optimx","plyr","pROC","randomForest","rpart")
      
      output <- foreach(block=1:blocklength,.combine=acomb,
                        .packages=packages.export) %dopar% {
                          parameters$startDate.index <- block
                          out.est <- estimate_cv(data.global,methodz,parameters,blocklength)
                          res.cv <- out.est$res.cv
                          data.db.train <- out.est$data.db.train
                          pos.eval.cv <- out.est$pos.eval.cv
                          out.table <- do.tables(data.eval=data.db.train,resmat=res.cv,pos.eval=pos.eval.cv,methodz=methodz,parameters=parameters,
                                                 doAUC=doAUC,doBPS=doBPS,doFmeas=doFmeas)
                          return(restable.array=drop(out.table))
                        }
      temp.results.all <- output
    }else{
      for (block in 1:blocklength){
        print(paste("starting block",block))
        parameters$startDate.index <- block
        out.est <- estimate_cv(data.db,methodz,parameters,blocklength)
        res.cv <- out.est$res.cv
        data.db.train <- out.est$data.db.train
        pos.eval.cv <- out.est$pos.eval.cv
        out.table <- do.tables(data.eval=data.db.train,resmat=res.cv,pos.eval=pos.eval.cv,methodz=methodz,parameters=parameters,
                               doAUC=doAUC,doBPS=doBPS,doFmeas=doFmeas)
        restable.array <- out.table
        if (block==1){
          res.cv.all <- array(dim=c(dim(res.cv),blocklength),dimnames=c(dimnames(res.cv),list(1:12)))
          temp.results.all <- array(dim=c(dim(drop(restable.array)),blocklength),dimnames=c(dimnames(drop(restable.array)),list(1:12)))
        }
        ndim.temp <- length(dim(temp.results.all))
        res.cv.all[,,block] <- res.cv
        if (ndim.temp==2){
          temp.results.all[,block] = drop(restable.array)   # restable.array: methodz, stats
        }else if (ndim.temp==3){
          temp.results.all[,,block] = drop(restable.array)  
        }
      }
    }
    ndim.temp <- length(dim(temp.results.all))
    temp.results <- apply(temp.results.all,c(1:(ndim.temp-1)),mean,na.rm=TRUE)
    ndim.temp <- length(dim(temp.results))
    
    for (m in methodz){
      n.params <- length(params[[m]][p,])
      results.paras.array[p,1:n.params,m,] = as.matrix(params[[m]][p,])
      if (num.methods == 1){
        results.paras.array[p,statnames,m,] = temp.results[statnames]
      }else{
        results.paras.array[p,statnames,m,] = temp.results[m,statnames]
      }
    }
    toc()
  }

  methodz.results <- methodz.results.all
  
  ############## Final Processing Block ###############################################################
  # Optimization and Graphs
  #     stat.select = "AUC" # 1: U_r, 2: AUC # moved to controls block
  #     res.select <- 2
  num.methods <- length(methodz.results)    
  output.table <- matrix(nrow=num.methods,ncol=maxNumParameters*2+1+length(statnames))
  colnames(output.table) <- c("method",paste("p",1:maxNumParameters,"_val",sep=""),paste("p",1:maxNumParameters,"_name",sep=""),statnames)
  for (mm in 1:num.methods){
    res.mm <- results.paras.array[!is.na(results.paras.array[,stat.select,mm,1]),,mm,res.optimize]
    if (length(dim(res.mm))==2){
      recursive.maxpos = which.max(res.mm[,stat.select])
      pval <- res.mm[recursive.maxpos,1:maxNumParameters]
      statval <- res.mm[recursive.maxpos,statnames]
    }else{
      # logit case
      pval <- rep(NA,maxNumParameters)
      statval <- res.mm[statnames]
    }
    
    pnames <- colnames(params[[mm]])
    pnames <- c(pnames,rep(NA,maxNumParameters-length(pnames)))
    output.table[mm,] <- c(methodz.results[mm],pval,pnames,statval)
  }
  write.csv(output.table,paste(wd_output_detail,"/",paramtable_file,sep=""),row.names=FALSE)
}

if (do.parallel){stopCluster(cl)}
save.image(paste(filename,".RData",sep=""))  # save final workspace after running full script_master
toc()