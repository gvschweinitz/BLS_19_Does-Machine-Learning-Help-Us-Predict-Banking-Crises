## classification evaluation criteria functions
calculate.threshold <- function(realizations, predictions, mu, optimizethreshold = FALSE, evaluate = FALSE, PriorP=NULL) {
  # ------------------------------------
  # CALL
  #   calculate.threshold(realizations, predictions, mu)
  #   calculate.threshold(realizations, predictions, mu, optimizethreshold = FALSE)
  #   calculate.threshold(realizations, predictions, mu, optimizethreshold = FALSE, evaluate = FALSE)
  #   calculate.threshold(realizations, predictions, mu, optimizethreshold = FALSE, evaluate = FALSE, PriorP=NULL)
  # ------------------------------------
  # DESCRIPTION
  #   Calculates confusion matrix and statistics from realizations and predicted probabilities
  # ------------------------------------
  # INPUT
  #   realizations:   binary pre-crisis variable
  #   predictions:    probability or signal of being in a pre-crisis state
  #   mu:             relative weight of missing a crisis. Note that reasonable values are different for usefulness functions of Alessi-Detken and Sarlin
  #   optimizethreshold:     TRUE/FALSE if threshold is to be optimized. If FALSE, only one threshold applied to predictions
  #   evaluate:       Type of threshold to be applied to the data (if optimizethreshold = FALSE). TRUE: threshold=0.5; FALSE: long-run optimal threshold (Sarlin and von Schweinitz) 
  #   PriorP:         Prior crisis probability. If missing, set to the prior probability in realizations. Set this value if thresholds are applied to out-of-sample data
  # ------------------------------------
  # OUTPUT
  #   list with entries:
  #         - matrix.performance: confusion matrix and statistics for all tested thresholds
  #         - threshold:          selected threshold
  # Gregor von Schweinitz 2020-11-17

  
  if (length(realizations) != length(predictions)){
    stop('Error in libraries_new.R : calculate.threshold. Number of predictions has to match number of obs in evaluation data set.')
  }
  
  if (is.null(PriorP)){
    # compute P
    PriorP <- mean(realizations) 
  }else if (length(PriorP)>1){
    stop('Error in libraries_new.R : calculate.threshold. PriorP needs to be a single number')
  }

	if (optimizethreshold == TRUE && evaluate == FALSE) {
		# Define number of grid points for threshold checking
		num.gridpoints <- 100
		# Construct probability grid to be used in loop, consisting of 100 points between minimum probability and maximum probability of method
		prob.grid <- seq(min(predictions, na.rm = TRUE), max(predictions, na.rm = TRUE), length.out = (num.gridpoints + 1))
	} else if (evaluate == FALSE) {	
		num.gridpoints <- 0
		if (parameters$u.alessidetken) {
			prob.grid <- (1 - mu) * PriorP / (mu * (1 - PriorP) + (1 - mu) * PriorP)
		} else {
		  prob.grid <- 1 - mu
		}
	} else if (evaluate == TRUE) {	# if "evaluate" is supplied, then this function is used only to evaluate predictions
		num.gridpoints <- 0
		prob.grid <- 0.5	# could be any number between 0 and 1
	}
	# Construct array
  X <- array(0, dim=c((num.gridpoints + 1), 8))
  colnames(X) <- c("Threshold","TP","FP","TN","FN","FPrate","FNrate","U_r")
  N <- length(realizations)
  X[,"Threshold"] <- prob.grid
  X[,"TN"] <- sapply(prob.grid,FUN=function(x){sum(predictions[realizations==0]<x)})
  X[,"TP"] <- sapply(prob.grid,FUN=function(x){sum(predictions[realizations==1]>=x)})
  X[,"FN"] <- sapply(prob.grid,FUN=function(x){sum(predictions[realizations==1]<x)})
  X[,"FP"] <- sapply(prob.grid,FUN=function(x){sum(predictions[realizations==0]>=x)})
  X[,"FPrate"] <- X[,"FP"]/(X[,"FP"]+X[,"TN"])
  X[,"FNrate"] <- X[,"FN"]/(X[,"FN"]+X[,"TP"])
  if (parameters$u.alessidetken == TRUE) {
    minU <- min(mu, (1 - mu))
    U <- minU - (mu * X[,"FNrate"] + (1 - mu) * X[,"FPrate"])
    X[,"U_r"] <- U / minU
  } else {
    minU <- min(PriorP*mu,(1 - mu)*(1-PriorP)) # with Priorp, simplified equation
    U <- minU - (mu*X[,"FN"]/N+(1 - mu)*X[,"FP"]/N)
    X[,"U_r"] <- U / minU
  }
	threshold <- X[which.max(X[, dim(X)[2]]), 1]
	if (evaluate) threshold <- 0.5
  return(list("matrix.performance" = X, "threshold" = threshold))
}


## function for creating pre-crisis indicator for classification 
calculate.pre <- function(data, horizon, posthorizon) {
  # ------------------------------------
  # CALL
  #   calculate.pre(data, horizon, posthorizon)
  # ------------------------------------
  # DESCRIPTION
  #   Extends the dataset, including columns "pre", "latepre" and "post".
  #   The code makes use of the fact that data is a data.table
  # ------------------------------------
  # INPUT
  #   data:           data.table containing at minimum the columns Country, Date and c0
  #   horizon:        1x2 vector containing lower and upper horizon of pre-crisis window
  #   posthorizon:    upper horizon of post-crisis window. Can be zero
  # ------------------------------------
  # OUTPUT
  #   data            Original data set with three additional columns pre, latepre and post
  # Gregor von Schweinitz 2020-11-17
  
  pre.vec <- -horizon[2]:-horizon[1]
  latepre.vec <- (-horizon[1]+1):0
  post.vec <- 0:posthorizon
  
  data[,pre := 0]
  data[,latepre := 0]
  data[,post := 0]
  
  country.list <- unique(data$Country)
  
  for (country.name in country.list) {
    data.temp <- data[which(data$Country == country.name),]
    data.temp <- data.temp[CJ(Country,Date, unique=TRUE)]
    N.c <- dim(data.temp)[1]
    # only do if crises exist in current data
    if (sum(data.temp[, c0] > 0, na.rm = T)) {
      crisis.indices <- which(data.temp[, c0] == 1)
      pre.vec.c <- sort(unique(pmax(pmin(unlist(lapply(crisis.indices,FUN = function(x){x+pre.vec})),N.c),1)))
      latepre.vec.c <- sort(unique(pmax(pmin(unlist(lapply(crisis.indices,FUN = function(x){x+latepre.vec})),N.c),1)))
      post.vec.c <- sort(unique(pmax(pmin(unlist(lapply(crisis.indices,FUN = function(x){x+post.vec})),N.c),1)))
      data.temp[setdiff(post.vec.c,crisis.indices),post := 1]
      data.temp[setdiff(latepre.vec.c,crisis.indices),latepre := 1]
      data.temp[setdiff(pre.vec.c,c(latepre.vec.c,crisis.indices)),pre:=1]
      data[Country == country.name,] <- data.temp[complete.cases(data.temp),]
    }
  }
  return(data)
}

# function to drop data from end of supplied data equal to forecast horizon, if no "C0" present
drop.c0 <- function(data, horizon, training = TRUE) {	
  # ------------------------------------
  # CALL
  #   drop.c0(data, horizon, training = TRUE)
  # ------------------------------------
  # DESCRIPTION
  #   Drops observations at the end of the given sample where it is not possible to know if they are pre-crisis or not in pseudo-real-time.
  # ------------------------------------
  # INPUT
  #   data:           data.table containing at minimum the columns Country, Date, c0 and pre
  #   horizon:        1x2 vector containing lower and upper horizon h2 of pre-crisis window. Up to h2 periods are removed at the end of the dataset.
  #   training:       Boolean to indicate if the data should be used for training. If (after dropping) there are no crisis left in the data set, throws an error
  # ------------------------------------
  # OUTPUT
  #   data            data set after removing observations with unknown pre-crisis status
  # Gregor von Schweinitz 2020-11-17
  
  maxDate <- max(data$Date)
  data[,keep:=1]
  data[is.na(c0),keep:=0]
	country.list <- unique(data$Country)
	for (country.name in country.list) {
		data.temp <- data[Country == country.name,]
		data.temp <- data.temp[complete.cases(data.temp[, c0]), ] # remove possible NA in C0
		
		# helper
		last <- data.temp[Date > maxDate - horizon[2]/4,]
		if (sum(last[, c0]) == 0) {		# drop based on crisis column
		  data[.(country.name,last[,Date]),keep:=0]
		} else if (sum(last[, c0]) > 0) {
			pos <- which(last[, c0] == 1)	# find position of last actual crisis
			pos <- pos[length(pos)]		
			if (pos < dim(last)[1]) {	# if post-crisis observations exist, drop these (regardless of post-horizon)
				last <- last[(pos + 1):dim(last)[1], ]
				data[.(country.name,last[,Date]),keep:=0]
			}
		}
	}
	data <- data[keep==1,]
	data[,keep:=NULL]
	
  if (training){ # only perform check on pre-crisis periods > 0 if we are using drop.c0 on training data
  	if (nrow(data[pre==1,]) == 0) stop("Dropping of unknown events according to forecasting window caused no pre-crisis observations in training data")	# error message if horizon is too large
  }
	return(data)
}

			
columns <- function(col.resp, cols.expl) {	
  # helper function for defining output column and input formulas for all methods
  # create formula cols.string
  return(paste(col.resp,"~",paste(cols.expl,collapse=" + ")))
}
