train_bs_par <- function(data.db.train,data.db.test,parameters,methodz,do.bootstrap=TRUE,R.BS=1000,blocklength=12,sample.countries = FALSE,min_cr = 5){
  # ------------------------------------
  # CALL
  #   train_bs(data.db.train,data.db.test,parameters,methodz,R.BS,blocklength)
  #   train_bs(data.db.train,data.db.test,parameters,methodz,R.BS,blocklength,sample.countries)
  #   train_bs(data.db.train,data.db.test,parameters,methodz,R.BS,blocklength,sample.countries,min_cr)
  # ------------------------------------
  # DESCRIPTION
  #   performs a bootstrap on data from data.db.train in order to provide a bootstrapped forecast on data.db.test
  #   In case of in-sample estimation, data.db.train (for bootstrapping) and data.db.test (for forecasting) should be identical
  #   Bootstrap specificities:
  #     - block bootstrap
  #     - cross sectional blocks (cannot be excluded at the current moment)
  #     - panel out-of-sample (dropping whole countries)
  #   Note that this function transforms input data.tables into data.frames for estimation purposes
  # ------------------------------------
  # INPUT
  #   data.db.train:  in-sample estimation data.table with country, year, quarter, c0, pre, latepre, post, and early warning indicators
  #   data.db.test:   out-of-sample prediction data.table with country, year, quarter, c0, pre, latepre, post, and early warning indicators
  #   cols.string:    column names (pre and estimation variables)
  #   parameters:     parameters from preparation scripts
  #   methodz:        vector with methods to be estimated
  #   do.bootstrap:   Boolean if bootstrap on training data should be performed or not
  #   R.BS:           number of bootstrap draws
  #   blocklength:    number of adjacent periods in a block; RECOMMENDED: EARLY WARNING HORIZON. Reduced if only smaller blocks are available. 
  #   sample.countries:   TRUE for panel out-of-sample bootstrap, FALSE otherwise (NOT RECOMMENDED FOR RECURSIVE ESTIMATION). 
  #   min_cr:         minimum number of pre-crisis observations in bootstrap sample. Needed for valid estimations
  # ------------------------------------
  # OUTPUT
  #   list containing
  #     - data.out:   T x output x R.BS array, where 
  #                     T is the number of observations in data.db.test
  #                     output is the number of outputs (methodz x c(Prob, PrioP, AUC, MSE))
  # Gregor von Schweinitz 2020-11-17
  
  if (!do.bootstrap){
    return(train_bs(data.db.train,data.db.test,parameters,methodz,do.bootstrap=do.bootstrap))
  }
  
  # have parallel framework only for bootstrap estimations
  b <- blocklength
  # Define all possible start and end-dates for bootstrap blocks
  # Account for the fact that blocks might need to be shorter than b periods in order to give every observation an equal chance to be drawn
  # In order to manage shorter blocks at the beginning, we give the additional (b-1) starting observations "fictitious earlier starting dates" 
  time.vec <- sort(unique(data.db.train[,Date]))
  time.add <- seq(-(b-1)/4,0,1/4)
  time.vec <- c(time.vec[1]+time.add,time.vec[2:length(time.vec)])
  sample.master <- data.frame(pos1=NA,pos2=NA,Date=NA,Country=NA)
  for (coun in unique(data.db.train$Country)){
    pos <- which(data.db.train$Country==coun)
    pos1 <- c(rep(pos[1],b-1),pos)
    pos2 <- c(pos,rep(pos[length(pos)],b-1))
    time.coun <- data.db.train[pos,Date]
    time.coun1 <- data.db.train[pos1,Date]
    time.coun2 <- data.db.train[pos2,Date]
    time.add1 <- c(time.add[1:min(b,length(pos1))],rep(0,max(0,length(pos1)-b)))
    # Account for missing observations, as latepre, c0, post are not included in data.db.train
    pos.na <- which(time.coun1+b/4-0.25<time.coun2)
    for (k in pos.na){pos2[k] <- pos[max(which(time.coun<=time.coun1[k]+b/4-0.25))]}
    time.coun2 <- data.db.train[pos2,Date]
    sample.master <- rbind(sample.master,cbind(pos1=pos1,pos2=pos2,Date=time.coun1+time.add1,Country=coun))
  }
  sample.master <- sample.master[2:dim(sample.master)[1],]
  rm(time.add,coun,pos,pos1,pos2,time.coun,time.coun1,time.coun2,time.add1,pos.na)
  
  output.names <- c("PriorP",
                    paste("Prob(", methodz, ")", sep=""),
                    paste("OT(", methodz, ")", sep=""))
  data.out.r <- array(dim=c(dim(data.db.test)[1],length(output.names)),
                    dimnames=list(rownames(data.db.test),output.names))
  
  methodz.agg <- parameters$methodz.agg  # adjusted, if a method yields no results
  
  Nobs <- dim(data.db.train)[1]
  data.test <- cbind(obs=Nobs+1:dim(data.db.test)[1],data.db.test)
  
  
  # FOREACH LOOP OVER BOOTSTRAP DRAWS
  
  acomb <- function(...) abind(...,along=3)
  packages.export <- c("data.table","e1071","kknn","nnet","optimx","plyr","pROC","randomForest","rpart")
  vars.export <- c("columns","calc.train.bs","calculate.threshold")
  print("starting parallel loop")
  # data.out <- foreach(r=1:R.BS,.combine=acomb,.packages=packages.export,.export=vars.export) %do% {
  data.out <- foreach(r=1:R.BS,.combine=acomb,.packages=packages.export,.export=vars.export) %dopar% {
    if(r%%10==0) print(r)
    #### sample data ####
    PriorP <- 0
    while (PriorP<min_cr/Nobs){
      # randomly drawn countries (with replacement)
      if (sample.countries) {
        countries.bs <- sample(unique(data.db.train$Country),replace=TRUE)
      }else{
        countries.bs <- unique(data.db.train$Country)
      }
      pos.bs1 <- unlist(lapply(countries.bs,grep,sample.master$Country))
      sample.bs <- sample.master[pos.bs1,]      # restrict sample.master to selected countries
      # randomly drawn starting times
      time.bs <- sample(time.vec,replace=TRUE)
      pos.bs2 <- unlist(lapply(time.bs,grep,sample.bs$Date))
      # get blocks from data.train
      pos.bs <- unlist(apply(sample.bs[pos.bs2,],1,function(x) seq(as.numeric(x["pos1"]),as.numeric(x["pos2"]))))
      data.train <- data.db.train[pos.bs[1:Nobs],]
      PriorP <- mean(data.train[,pre])
      rm(countries.bs,pos.bs1,sample.bs,time.bs,pos.bs2,pos.bs)
    }
    realizations <- data.train[,pre]
    
    #### Standardization ####
    if (parameters$standardizeList){
      sdvars <- unique(unlist(parameters$cols.expl))
      colStats <- data.train[,lapply(.SD, function(x) c(mean(x),sd(x))),.SDcols = sdvars]
      data.train.norm = data.train
      data.test.norm = data.test
      for (sdcol in sdvars){
        data.train.norm[,(sdcol) := scale(data.train[,sdcol,with=FALSE],colStats[1,sdcol,with=FALSE],colStats[2,sdcol,with=FALSE])]
        data.test.norm[,(sdcol) := scale(data.test[,sdcol,with=FALSE],colStats[1,sdcol,with=FALSE],colStats[2,sdcol,with=FALSE])]
      }
    }else { # no standardization:
      data.train.norm = data.train
      data.test.norm = data.test
    }
    
    for (method in methodz){
      # SET PARAMETERS BY METHOD
      method.calc <- strsplit(method,".",fixed=TRUE)[[1]][1]
      dnum <- as.numeric(strsplit(method,".",fixed=TRUE)[[1]][2])
      
      usevars <- c("Country","Date",parameters$cols.expl[[dnum]],"c0",parameters$col.resp,"latepre","post")
      data.train.use <- as.data.frame(data.train.norm[,usevars,with=FALSE])
      data.test.use <- as.data.frame(data.test.norm[,usevars,with=FALSE])
      # parameters$logitvariables <- parameters$cols.expl[[dnum]]
      # parameters$datavariables <- parameters$cols.expl[[dnum]]
      cols.string <- columns(parameters$col.resp, parameters$cols.expl[[dnum]])
      
      if (!parameters$optimizationMode){
        # Set hyperparameters of method based on paramtable (which comes from an exogenous csv-file specifying hyperparameters for all estimation methods)
        namecols <- grep("_name",colnames(parameters$paramtable))
        valcols <- grep("_val",colnames(parameters$paramtable))
        for (col in 1:length(namecols)){
          if (!is.na(parameters$paramtable[method,namecols[col]])){
            parameters[[parameters$paramtable[method,namecols[col]]]] <- parameters$paramtable[method,valcols[col]]
          }
        }
      }
      
      train.out <- NULL
      try(train.out <- calc.train.bs(method.calc, cols.string, data.train.use, data.test.use, parameters),silent=FALSE)
      if (!is.null(train.out)){
        data.out.r[,paste("Prob(", method, ")", sep="")] <- train.out$temp.test
        threshold <- calculate.threshold(realizations, train.out$temp.train, parameters$mu, optimizethreshold = parameters$optimizethreshold, evaluate = FALSE, PriorP=PriorP)$threshold
        data.out.r[,paste("OT(", method, ")", sep="")] <- threshold
      }else{
        data.out.r[,paste("Prob(", method, ")", sep="")] <- NA
        data.out.r[,paste("OT(", method, ")", sep="")] <- mean(realizations)
      }
    }
    
    return(data.out.r)
  }
  
  return(list(data.out=data.out))
  
}