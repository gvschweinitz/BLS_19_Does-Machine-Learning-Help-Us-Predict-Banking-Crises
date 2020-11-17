estimate_cv <- function(data.db,methods.all,parameters,blocklength){
  # ------------------------------------
  # CALL
  #   estimate_cv(data.db,methods.all,parameters,blocklength)
  # ------------------------------------
  # DESCRIPTION
  #   Produces estimated out-of-sample probabilities based on cross-validation. 
  #   The data are split into different blocks of a given length, containing observations from all countries. 
  #   Every block is used as test sample once, with all other blocks forming the training sample. That is, every observation is put into the cross-validation (out-of-sample) dataset once.
  #   Output collects estimated probabilities coming from different estimations.
  # ------------------------------------
  # INPUT
  #   data.db:        data.frame with country, year, quarter, c0, pre, latepre, post, and early warning indicators
  #   methods.all:    Set of models to be estimated
  #   parameters:     structure containing control parameters.
  #   blocklength:    length of blocks (number of periods)
  # ------------------------------------
  # OUTPUT
  #   list with:      - Output (Prior probability, out-of-sample estimated probability and optimal threshold) for all observations
  #                   - training data
  #                   - estimation model
  #                   - position in data and output to be used for evaluation (all positions)
  #   
  #   Gregor von Schweinitz, 2020-11-16


  ##########################################################
  posstart <- parameters$startDate.index
  if (posstart>blocklength){stop("CV starting date too late")}
  
  # CROSSVALIDATED ESTIMATION
  # construct "training sample" (i.e. full sample for estimation)
  data.db.train <- data.db
  data.db.train <- data.db.train[data.db.train$c0_NA!=1,] # drop if c0 = NA (useless for training sample, note: strictly speaking it could happen that pre = 0/1 but crisis = NA,
  # in practice however, this is not the case in our data set)
  if (parameters$drop.horizon) data.db.train <- drop.c0(data.db.train, parameters$horiz,parameters$col.crisis) 
  # drop.c0 drops observations at the end of the given sample where it is not possible to know if they are pre-crisis or not 
  data.db.train <- data.db.train[complete.cases(data.db.train),]  # drop if rows contain NA (RHS variables), note: don't move this line before drop.c0 (if we had c0 information but NA in RHS vars this would lead to wrong elimination by drop.c0)
  data.db.train <- data.db.train[post==0,] # For training data, remove late-pre-, c0- and post-crisis periods
  data.db.train <- data.db.train[c0==0, ]
  data.db.train <- data.db.train[latepre==0,]
  
  # Prepare output variables
  output.names <- c("PriorP",
                    paste("Prob(", methods.all, ")", sep=""),
                    paste("OT(", methods.all, ")", sep=""))
  res.cv <- array(dim=c(dim(data.db.train)[1],length(output.names)),
                  dimnames=list(rownames(data.db.train),output.names))
  pos.eval.cv <- 1:dim(res.cv)[1]
  
  
  # Construction of crossvalidation folds
  time.vec <- sort(unique(data.db.train[,Date]))
  numfolds <- ceiling(length(time.vec)/blocklength)
  fold <- c(rep(numfolds,posstart-1),rep(1:numfolds,each=blocklength))
  fold <- fold[1:length(time.vec)]
  
  for (f in 1:numfolds){
    time.vec.f <- time.vec[fold==f]
    data.cv.train <- data.db.train[!Date%in%time.vec.f,]
    data.cv.test <- data.db.train[Date%in%time.vec.f,]
    out.cv <- train_bs(data.cv.train,data.cv.test,parameters,methodz,do.bootstrap=FALSE)
    res.cv[data.db.train$Date%in%time.vec.f,] <- drop(out.cv$data.out)
  }
  
  
  return(list(res.cv=res.cv,data.db.train=data.db.train,pos.eval.cv=pos.eval.cv))
}