do.tables <- function(data.eval,resmat,pos.eval,methodz,parameters,doAUC=TRUE,doBPS=TRUE,doFmeas=TRUE){
  # ------------------------------------
  # CALL
  #   do.tables(data.eval,resmat,pos.eval,methodz,parameters,doAUC=TRUE,doBPS=TRUE,doFmeas=TRUE)
  # ------------------------------------
  # DESCRIPTION
  # Function creates result tables for cross-validation estimations that are needed to determine the optimal hyperparameters
  # Result tables contain for every method
  #   - entries of the confusion matrix and derived statistics
  #   - Aggregate statistics: U_r, AUC, BPS, Fmeas (depending on selection in script_prep_gen.R)
  # ------------------------------------
  # INPUT
  #   data.eval:      estimation data.table with country, year, quarter, c0, pre, latepre, post, and early warning indicators
  #   resmat:         estimation results from estimate_cv.R
  #   methodz:        methodz for which results are to be calculated
  #   parameters:     parameters from preparation scripts
  #   doAUC:          Boolean if AUC is to be calculated
  #   doBPS:          Boolean if BPS is to be calculated
  #   doFmeas:        Boolean if F-measure is to be calculated
  # ------------------------------------
  # OUTPUT
  #   restable.array  array of dimensions methodz x results, containing all relevant result statistics for every method
  
  # Gregor von Schweinitz, 2020-11-17
  
  ##########################################################
  # RESULT TABLES
  
  # select dataset (in-sample or recursive)
  initialized = 0
  mu <- parameters$mu
  methodz.results <- methodz
  data.train <- data.eval
  resmat <- resmat[pos.eval,]
  merged = cbind(data.train,resmat)
  
  # Input checking:
  if (dim(data.train)[1] != dim(resmat)[1]){  # jb added
    stop('Error in script_tables.R: Number of predictions (resmat) has to match number of obs in evaluation data set (data.train).')
  }
  realizations = merged[ ,pre]
  PriorP = mean(data.db.train[,pre],na.rm=TRUE)  # used for evaluation of predictions
  PriorP.OT = merged[,PriorP]

  for (i in 1:length(methodz.results)){
    method <- methodz.results[i]
    
    # compute signals (for given method and threshold)
    if (parameters$optimizethreshold){OT <- merged[[paste("OT(",method,")",sep="")]]}
    probabilityPredictions = merged[[paste("Prob(",method,")",sep="")]]
    predictions <- (merged[[paste("Prob(",method,")",sep="")]] >= OT)*1
    
    mergedNew = cbind(merged[,.(Country,Date)], PriorP, probabilityPredictions, predictions, realizations)  # for checking
    res <- calculate.threshold(realizations, predictions, mu, optimizethreshold = FALSE, evaluate = TRUE, PriorP=PriorP)$matrix.performance
    res[,"Threshold"] <- mean(OT)
    
    # compute AUC
    if (doAUC){ # check if there are both tranquil periods and pre-crisis periods
      checkVec = c(0,1)
      if (all(checkVec %in% unique(realizations[!is.na(probabilityPredictions)]))){
        rocest <- roc(as.factor(realizations), probabilityPredictions,direction="<",quiet=TRUE)
        AUC <- rocest$auc
      }else{
        AUC <- NA
      }
      res <- cbind(res,AUC=AUC)  # add AUC to results
    }
    if (doBPS){
      BPS <- mean((probabilityPredictions-realizations)^2)
      res <- cbind(res,BPS=BPS)
    }
    if (doFmeas){
      Fmeas <- res[,"TP"]/(res[,"TP"]+mean(res[,c("FN","FP")]))
      res <- cbind(res,Fmeas=Fmeas)
    }
    res <- cbind(mu=mu,opt=parameters$optimizethreshold*1,type=parameters$u.alessidetken*1,res)
    if (initialized == 0){
      restable.array <- array(dim=c(length(methodz.results),length(res)),
                              dimnames=list(methodz.results,colnames(res)))
      initialized = 1
    }
    restable.array[i,] <- res  # dims: methodz, stats, in-sample/OOS 
  }
  
  return(restable.array)
}
