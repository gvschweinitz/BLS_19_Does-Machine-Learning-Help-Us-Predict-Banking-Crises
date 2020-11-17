calc.train.bs <- function(method, cols.string, data.train, data.test, parameters) {
  # ------------------------------------
  # CALL
  #   calc.train.bs(method, cols.string, data.train, data.test, parameters)
  # ------------------------------------
  # DESCRIPTION
  #   Estimates a probability model. Uses estimated model to produce predictions on training and test data set for a given method
  # ------------------------------------
  # INPUT
  #   method:         method to be estimated
  #   cols.string:    column names (pre and estimation variables)
  #   data.train:     in-sample estimation data.frame with country, year, quarter, c0, pre, latepre, post, and early warning indicators
  #   data.test:      out-of-sample prediction data.frame with country, year, quarter, c0, pre, latepre, post, and early warning indicators
  #   parameters:     parameters from preparation scripts
  # ------------------------------------
  # OUTPUT
  #   list with:      - predicted probabilities on data.train
  #                   - predicted probabilities on data.test
  #                   - estimated model
  #   
  #   Gregor von Schweinitz, 2020-11-16
  #   Partially based on CrisisModeler by Peter Sarlin and Markus Holopainen
  
  # Create response column name, required for logit and knn, as well as explanatory column names
  col.resp <- strsplit(cols.string, " ~")[[1]][1]
  cols.expl <- strsplit(cols.string, "pre ~ ")[[1]][2]	
  cols.expl <- strsplit(cols.expl, " . ")[[1]]
  
  if (method == "logit") {
    model <- glm(cols.string, family = binomial(link="logit"), data = data.train)
    temp.train <- predict(model, newdata = data.train, type = "response")
    temp.test <- predict(model, newdata = data.test, type = "response")
  }
  if (method == "trees") {
    model <- prune(rpart(cols.string, method = "class", data = data.train),cp = parameters$treescp)
    temp.train <- predict(model, newdata = data.train)[,2]
    temp.test <- predict(model, newdata = data.test)[,2]
  }
  if (method == "knn") {
    data.train[, col.resp] <- as.factor(data.train[, col.resp])
    data.test[, col.resp] <- as.factor(data.test[, col.resp])
    data.all <- rbind(data.train,data.test) # combine data.train and data.test for joint prediction
    T.train <- dim(data.train)[1]
    T.all <- dim(data.all)[1]
    model <- kknn(formula(cols.string), data.train, data.all, k = parameters$knnk, distance = 2, kernel = "optimal")
    temp.train <- model$prob[1:T.train, 2]
    temp.test <- model$prob[(T.train+1):T.all, 2]
  }
  if (method == "rf") {
    data.train[, col.resp] <- as.factor(data.train[, col.resp]) # if classification, temporary switch response to factor, changes mode in randomforest()
    model <- randomForest(formula(cols.string), data = data.train, ntree = 1000, mtry = parameters$rfmtry, maxnodes = NULL, nodesize=parameters$nodesize)
    temp.train <- predict(model, newdata = data.train, type = "prob")[, 2]
    temp.test <- predict(model, newdata = data.test, type = "prob")[, 2]
  }	
  if (method == "nen") {
    # note - use sink() to suppress automatic output. Command different for non-Windows PCs
    sink("NUL") 
    model <- nnet(formula(cols.string), data = data.train, size = parameters$nensize, maxit = 10000, decay = parameters$nendecay)
    sink()
    temp.train <- as.numeric(predict(model, newdata = data.train))
    temp.test <- as.numeric(predict(model, newdata = data.test))
  }	
  if (method == "svm") {
    # SVM requires that inputs are only numerical columns, give explanatory and response only
    numeric.cols <- c(cols.expl, col.resp)
    model <- svm(formula(cols.string), data = data.train[, numeric.cols], gamma = parameters$svmg, cost = parameters$svmc, kernel = parameters$svmk, cross = 0, type = "C-classification", probability = parameters$svm.doProb)
    attributeSVM = "probabilities"
    temp.train <- attr(predict(model, data.train[, numeric.cols], probability = parameters$svm.doProb), attributeSVM)[, 2]				
    temp.test <- attr(predict(model, data.test[, numeric.cols], probability = parameters$svm.doProb), attributeSVM)[, 2]
  }
  return(list("temp.train" = temp.train, "temp.test" = temp.test, "model" = model))
  
}