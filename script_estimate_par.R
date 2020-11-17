# Script performing in- and out-of-sample estimation, possibly on parallel processors
# Gregor von Schweinitz, 2020-11-17

##########################################################
# IN-SAMPLE ESTIMATION
data.db.train <- data.global
data.db.train <- data.db.train[c0_NA!=1,] # drop if c0 = NA (useless for training sample, note: strictly speaking it could happen that pre = 0/1 but crisis = NA,
# in practice however, this is not the case in our data set)
if (parameters$drop.horizon) data.db.train <- drop.c0(data.db.train, parameters$horiz) 
# drop.c0 drops observations at the end of the given sample where it is not possible to know if they are pre-crisis or not 
data.db.train <- data.db.train[complete.cases(data.db.train),]  # drop if rows contain NA (RHS variables), note: don't move this line before drop.c0 (if we had c0 information but NA in RHS vars this would lead to wrong elimination by drop.c0)
data.db.train <- data.db.train[latepre == 0,] # For training data, remove late-pre-, c0- and post-crisis periods
data.db.train <- data.db.train[c0 == 0,] # For training data, remove late-pre-, c0- and post-crisis periods
data.db.train <- data.db.train[post == 0,] # For training data, remove late-pre-, c0- and post-crisis periods

# construct "test sample"
data.db.test <- data.global
pos.eval.is <- as.numeric(row.names(match_df(data.db.test[,1:2],data.db.train[,1:2]))) # match rows in data.db.test and data.db.train

############# PRINTING
print("")
print("Number of observations contained in full sample training data:")
print("All obs: ")
print(dim(data.db.train)[1])
print("Pre-crisis periods:")
print(sum(data.db.train$pre))
print("Corresponding to approximate number of crisis events:")
print(sum(data.db.train$pre)/(parameters$horiz[2]-parameters$horiz[1]+1))
print("")

tic("est.is")
# generate in-sample (full sample) results
if (do.is){
  out.is <- train_bs(data.db.train,data.db.test,parameters,methodz,do.bootstrap=FALSE)
  # parameters <- out.is$parameters
  res.is <- drop(out.is$data.out)
  modelz.is <- out.is$modelz
  if (do.bootstrap){
    out.bs.is <- train_bs_par(data.db.train,data.db.test,parameters,methodz,do.bootstrap=TRUE,
                              R.BS=R.BS,blocklength=blocklength,sample.countries=sample.countries,min_cr=min_cr)
    res.bs.is <- out.bs.is$data.out
    rm(out.bs.is)
  }
  
  rm(out.is)
}
toc()


##########################################################
# RECURSIVE OUT-OF-SAMPLE ESTIMATION
if (do.recursive){
  endDate.lastCrisis = max(data.global[c0==1,Date], na.rm=TRUE)
  if (endDate.lastCrisis - recursiveend >= parameters$horiz[2]/4){
    # Accounting for "unknown crisis status" at the end of the database is only necessary if recursive analysis includes the last periods
    recursiveDropCriterion = FALSE
  }else{
    recursiveDropCriterion = TRUE
  } 
  
  pos.recursive <- data.global$Date>=recursivestart & data.global$Date<=recursiveend
  time.vec.rec <- sort(unique(data.global[pos.recursive,Date]))
  
  print(paste("Starting recursive exercise from ", min(time.vec.rec), " with methods: ", paste(methodz, collapse=", "), sep=""))
  
  # create arrays to be filled with recursive point estimates (res.recursive) and recursive bootstrap results (res.bs.recursive)
  output.names <- c("PriorP",
                    paste("Prob(", methodz, ")", sep=""),
                    paste("OT(", methodz, ")", sep=""))
    
  res.recursive <- array(dim=c(sum(pos.recursive),length(output.names)),
                         dimnames=list(rownames(data.global[pos.recursive,]),
                                       output.names))
  res.bs.recursive <- array(dim=c(sum(pos.recursive),length(output.names),R.BS),
                            dimnames=list(rownames(data.global[pos.recursive,]),
                                          output.names,
                                          1:R.BS))
  cnt <- 0  
  
  
  # CONSTRUCT DATASETS needed for graphs and tables
  data.db.recursive.plot = data.global[pos.recursive,]  # data set against which to plot recursive predictions
  
  # for final evaluation (for tables) use drop.c0 to cut-off observations at end of sample and remove: post, latepre and c0
  data.db.recursive.eval <- data.db.recursive.plot[c0_NA==0,] # drop if c0 = NA (otherwise drop.c0 will generate error if c0 = NA at end of sample, which is likely)
  if (dim(data.db.recursive.eval)[1] == 0){
    print("Attention: Your recursive OOS window does not contain any crisis information. (Crisis/c0 is NA everywhere in window).")
  }else{
    if (parameters$drop.horizon & recursiveDropCriterion) data.db.recursive.eval <- drop.c0(data.db.recursive.eval, parameters$horiz, training = FALSE) 
    # drop.c0 drops observations at the end of the given sample where it is not possible to know if they are pre-crisis or not 
    data.db.recursive.eval <- data.db.recursive.eval[complete.cases(data.db.recursive.eval),]  # drop if rows contain NA (RHS variables), note: don't move this line before drop.c0 (if we had c0 information but NA in RHS vars this would lead to wrong elimination by drop.c0)
    data.db.recursive.eval <- data.db.recursive.eval[post==0,] # For training or evaluation data, remove late-pre, c0 and post-crisis periods
    data.db.recursive.eval <- data.db.recursive.eval[c0==0, ]
    data.db.recursive.eval <- data.db.recursive.eval[latepre==0,]
    print("")
    print("Number of observations contained in evaluation data set for recursive :")
    print("All obs: ")
    print(dim(data.db.recursive.eval)[1])
    print("Pre-crisis periods:")
    print(sum(data.db.recursive.eval$pre))
    print("Corresponding to approximate number of crisis events:")
    print(sum(data.db.recursive.eval$pre)/(parameters$horiz[2]-parameters$horiz[1]+1))
    print("")
  }
  
  # for each row in data.db.recursive.eval, find corresponding row in data.db.recursive.plot
  pos.eval.recursive <- as.numeric(row.names(match_df(data.db.recursive.plot[,1:2],data.db.recursive.eval[,1:2]))) 
  
  
  # start looping through dates in recursive OOS window
  tic("est")
  for (i in 1:length(time.vec.rec)) {
    tic("est.rec.pre")
    # Define local training data set for recursive (changes at each iteration)
    data.oos.train <- data.global[Date<time.vec.rec[i],]
    data.oos.train <- data.oos.train[c0_NA==0,] # drop if c0 = NA (useless for training sample, note: strictly speaking it could happen that pre = 0/1 but crisis = NA,
    if (parameters$drop.horizon) data.oos.train <- drop.c0(data.oos.train, parameters$horiz)
    data.oos.train <- data.oos.train[complete.cases(data.oos.train),]
    data.oos.train <- data.oos.train[post==0,]
    data.oos.train <- data.oos.train[c0==0, ]
    data.oos.train <- data.oos.train[latepre==0,]
    
    # Define local "test" data set for which train_bs will make prediction based on training data (changes at each iteration)
    pos.oos.test <- data.db.recursive.plot$Date == time.vec.rec[i]
    data.oos.test <- data.db.recursive.plot[pos.oos.test,]
    print(paste("recursive ", time.vec.rec[i],", training #obs: ", sum(data.oos.train$pre), " pre-crisis observations from ", nrow(data.oos.train), " total",sep=""))
    
    # generate point estimates (bootstrap = FALSE)
    toc()
    tic("est.rec")
    if (dim(data.oos.test)[1]>0){
      out.recursive <- train_bs(data.oos.train,data.oos.test,parameters,methodz,do.bootstrap=FALSE)
      res.temp <- drop(out.recursive$data.out)
      res.recursive[pos.oos.test,] <- res.temp
      # generate bootstrap estimates with bootstrap on training data set
      if (do.bootstrap){
        out.temp <- train_bs_par(data.oos.train,data.oos.test,parameters,methodz,do.bootstrap=TRUE,
                                 R.BS=R.BS,blocklength=blocklength,sample.countries=sample.countries,min_cr=min_cr)
        res.bs.recursive[pos.oos.test,,] <- out.temp$data.out
        rm(out.temp)
      }
      toc()
    }
  }
  toc()
}
