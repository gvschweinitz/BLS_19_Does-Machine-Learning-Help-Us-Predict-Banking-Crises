# Script creating result tables for in-sample and out-of-sample estimations
# Result tables for cross-validation exercises (not as extensive) are created using the function do.tables.R
# Result tables are computed and saved into the output folder specified by wd_output_detail
# Result tables contain for every method
#   - entries of the confusion matrix and derived statistics
#   - Aggregate statistics: U_r, AUC, BPS, Fmeas (depending on selection in script_prep_gen.R)
#   - Confidence bands and corresponding probability quantiles for aggregate statistics. Probability quantiles are shifted such that confidence bands contain the point estimate
#   - p-values of comparisons of aggregate statistics across methods (separate files)
# Gregor von Schweinitz, 2020-11-17

##########################################################
# RESULT TABLES
# select dataset (in-sample or recursive)
initialized = 0
mu <- parameters$mu
methodz.results <- methodz
alpha.width <- alphavec[2] - alphavec[1]
alpha.comp <- alpha.width/2
for (d in 1:length(datasets.eval)){
  # if (!optimizationMode){print(datasets.eval[d])}

  # observations and predictions for all countries
  data.train <- get(datasets.eval[d])               # observations
  resmat = get(resnames[d])[get(posvecs.eval[d]),]  # predictions (corresponding to observations)
  # if (do.bootstrap){resmat_all.bs <- get(resnames.bs[d])[get(posvecs.eval[d]),,]}
  merged = cbind(data.train,resmat)

  # Input checking:
  if (dim(data.train)[1] != dim(resmat)[1]){  # jb added
    stop('Error in script_tables.R: Number of predictions (resmat) has to match number of obs in evaluation data set (data.train).')
  }

  pos <- 0   # initialize for loop over method (and muvec): position in final output table
  realizations = merged[ ,pre]
  PriorP = mean(data.db.train[,pre],na.rm=TRUE)  # used for evaluation of predictions
  PriorP.OT = merged[,PriorP]
    
  
  if (!parameters$optimizethreshold){
    # for each mu calculate threshold (vector, depending on PriorP)
    if (parameters$u.alessidetken){   # AD usefulness
      OT <- (1 - mu) * PriorP.OT / (mu * (1 - PriorP.OT) + (1 - mu) * PriorP.OT) # optimal threshold for AD usefulness
    }else{  # Schweinitz, Sarlin usefulness
      OT <- 1-mu
    }
  }

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
        if (rocest$auc==1){
          AUC <- rep(1,3)
        }else{
          AUC <- ci.auc(rocest,conf.level=alphavec[2]-alphavec[1],partial.auc=c(0.5,1),partial.auc.focus="spec",partial.auc.correct=TRUE)
        }
      }else{
        AUC <- rep(NA,3)
      }
      cnames <- c(colnames(res),"AUC",paste("conf_",alphavec,"(AUC)",sep=""))
      res <- cbind(res,matrix(AUC[c(2,1,3)],nrow=1,ncol=3))  # add AUC to results
      colnames(res) <- cnames
    }
    if (doBPS){
      BPS <- mean((probabilityPredictions-realizations)^2)
      res <- cbind(res,BPS=BPS)
    }
    if (doFmeas){
      Fmeas <- res[,"TP"]/(res[,"TP"]+mean(res[,c("FN","FP")]))
      res <- cbind(res,Fmeas=Fmeas)
    }
    
    if (do.bootstrap){
      tempdata <- get(resnames.bs[d])[get(posvecs.eval[d]),grep(paste("(",method,")",sep=""),dimnames(get(resnames.bs[d]))[[2]],fixed=TRUE),]
      res.bs <- apply(tempdata,3,
                      function(x){
                        probabilityPredictions <- x[,paste("Prob(",method,")",sep="")]
                        predictions <- probabilityPredictions>=x[,paste("OT(",method,")",sep="")]
                        res <- calculate.threshold(realizations, predictions, mu, optimizethreshold = FALSE, evaluate = TRUE, PriorP=PriorP)$matrix.performance
                        # Don't do AUC for bootstrap draws
                        if (doAUC){res <- cbind(res,NA,NA,NA)}
                        if (doBPS){
                          BPS <- mean((probabilityPredictions-realizations)^2)
                          res <- cbind(res,BPS=BPS)
                        }
                        if (doFmeas){
                          Fmeas <- res[,"TP"]/(res[,"TP"]+mean(res[,c("FN","FP")]))
                          res <- cbind(res,Fmeas=Fmeas)
                        }
                        return(res)})
        rownames(res.bs) <- colnames(res)
        cnames <- c(colnames(res),"q.U","probs.U.low","probs.U.high",
                    paste("conf_",alphavec,"(U_r)",sep=""),"sd_U_r")
        q.U <- mean(res.bs["U_r",]<res[,"U_r"])
        probs.U <- c(max(q.U-alpha.comp+1-max(q.U+alpha.comp,1),0),min(q.U+alpha.comp-min(q.U-alpha.comp,0),1))
        res <- cbind(res,q.U,matrix(probs.U,nrow=1),
                     min(quantile(res.bs["U_r",],probs=probs.U[1],na.rm=TRUE),res[,"U_r"]),
                     max(quantile(res.bs["U_r",],probs=probs.U[2],na.rm=TRUE),res[,"U_r"]),
                     sd(res.bs["U_r",],na.rm=TRUE)
                     )
        if (doBPS){
          cnames <- c(cnames,"q.BPS","probs.BPS.low","probs.BPS.high",
                      paste("conf_",alphavec,"(BPS)",sep=""),"sd_BPS")
          q.BPS <- mean(res.bs["BPS",]<res[,"BPS"])
          probs.BPS <- c(max(q.BPS-alpha.comp+1-max(q.BPS+alpha.comp,1),0),min(q.BPS+alpha.comp-min(q.BPS-alpha.comp,0),1))
          res <- cbind(res,q.BPS,matrix(probs.BPS,nrow=1),
                       min(quantile(res.bs["BPS",],probs=probs.BPS[1],na.rm=TRUE),res[,"BPS"]),
                       max(quantile(res.bs["BPS",],probs=probs.BPS[2],na.rm=TRUE),res[,"BPS"]),
                       sd(res.bs["BPS",],na.rm=TRUE))
        }
        
        if (doFmeas){
          cnames <- c(cnames,"q.Fmeas","probs.Fmeas.low","probs.Fmeas.high",
                      paste("conf_",alphavec,"(Fmeas)",sep=""),"sd_Fmeas")
          q.Fmeas <- mean(res.bs["Fmeas",]<res[,"Fmeas"])
          probs.Fmeas <- c(max(q.Fmeas-alpha.comp+1-max(q.Fmeas+alpha.comp,1),0),min(q.Fmeas+alpha.comp-min(q.Fmeas-alpha.comp,0),1))
          res <- cbind(res,q.Fmeas,matrix(probs.Fmeas,nrow=1),
                       min(quantile(res.bs["Fmeas",],probs=probs.Fmeas[1],na.rm=TRUE),res[,"Fmeas"]),
                       max(quantile(res.bs["Fmeas",],probs=probs.Fmeas[2],na.rm=TRUE),res[,"Fmeas"]),
                       sd(res.bs["Fmeas",],na.rm=TRUE))
        }
        colnames(res) <- cnames
        
        if (i==1){
          res.bs.array <- array(NA,dim=c(length(methodz.results),dim(res.bs)),
                              dimnames = list(methodz.results,rownames(res.bs),colnames(res.bs)))
        }
        res.bs.array[i,,] <- res.bs
      }

    # format output
    res <- cbind(mu=mu,opt=parameters$optimizethreshold*1,type=parameters$u.alessidetken*1,res)
    #initialize restable.array
    if (initialized == 0){
      restable.array <- array(dim=c(length(methodz.results),length(res),length(datasets.eval)),
                              dimnames=list(methodz.results,colnames(res),resnames))
      initialized = 1
    }
    restable.array[i,,d] <- res  # dims: methodz, stats, in-sample/OOS 
  }
  
  if (doAUC & doAUC.pvals){
    # Compare AUCs of different methods, based on point estimates
    ptable.AUC <- matrix(NA,length(methodz.results),length(methodz.results),dimnames=list(methodz.results,methodz.results))
    diag(ptable.AUC) <- 1
    for (m in 2:length(methodz.results)){
      for (m2 in 1:(m-1)){
        method <- methodz.results[m]
        method2 <- methodz.results[m2]
        if (restable.array[method,"AUC",d]==1 & restable.array[method2,"AUC",d]==1){
          ptable.AUC[m2,m] <- 1
          ptable.AUC[m,m2] <- 1
        }else{
          ptable.AUC[m,m2] <- roc.test(as.factor(realizations),merged[[paste("Prob(",method,")",sep="")]],merged[[paste("Prob(",method2,")",sep="")]],
                                       alternative="greater",boot.n=500,progress="none",quiet=TRUE)$p.value
          ptable.AUC[m2,m] <- 1-ptable.AUC[m,m2]
        }
      }
    }
    pAUC.csv <- paste(wd_output_detail,"/pvals_AUC_",resnames[d],"_",mu,".csv",sep="")
    write.csv(ptable.AUC,file=pAUC.csv)
    # }
  }

  if (do.bootstrap & length(methodz.results)>1){
    # Compare U_r, BPS and Fmeas of different methods, based on bootstrap results
    ptable.Ur <- matrix(1,length(methodz.results),length(methodz.results),dimnames=list(methodz.results,methodz.results))
    if (doBPS){ptable.BPS <- ptable.Ur}
    if (doFmeas){ptable.Fmeas <- ptable.Ur}
    for (m in 2:length(methodz.results)){
      for (m2 in 1:(m-1)){
        ptable.Ur[m,m2] <- mean(res.bs.array[m,"U_r",]<res.bs.array[m2,"U_r",])
        ptable.Ur[m2,m] <- 1-ptable.Ur[m,m2]
        if (doBPS){
          ptable.BPS[m,m2] <- mean(res.bs.array[m,"BPS",]>res.bs.array[m2,"BPS",])
          ptable.BPS[m2,m] <- 1-ptable.BPS[m,m2]
        }
        if (doFmeas){
          ptable.Fmeas[m,m2] <- mean(res.bs.array[m,"Fmeas",]<res.bs.array[m2,"Fmeas",])
          ptable.Fmeas[m2,m] <- 1-ptable.Fmeas[m,m2]
        }
      }
    }
    pU.csv <- paste(wd_output_detail,"/pvals_Ur_",resnames[d],"_",mu,".csv",sep="")
    write.csv(ptable.Ur,file=pU.csv)

    if (doBPS){
      pBPS.csv <- paste(wd_output_detail,"/pvals_BPS_",resnames[d],"_",mu,".csv",sep="")
      write.csv(ptable.BPS,file=pBPS.csv)
    }
    if (doFmeas){
      pFmeas.csv <- paste(wd_output_detail,"/pvals_Fmeas_",resnames[d],"_",mu,".csv",sep="")
      write.csv(ptable.Fmeas,file=pFmeas.csv)
    }
  }
  
  if (!do.crossval){
    # write Table in .csv (one table for every dataset, country - containing all mu and methods)
    filename.csv <- paste(wd_output_detail,"/",resnames[d],".csv",sep="")
    write.csv(restable.array[,,d],file=filename.csv)
  }
} # end of datasets loop
