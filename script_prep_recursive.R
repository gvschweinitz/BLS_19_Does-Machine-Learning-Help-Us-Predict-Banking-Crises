####################################
# Setting some names
datasets.eval <- "data.db.recursive.eval"
resnames <- "res.recursive"
posvecs.eval <- "pos.eval.recursive"  # posvecs.eval is used to select from resnames (in script_tables)
resnames.bs <- "res.bs.recursive"

##################################
# Bootstrap settings

do.bootstrap <- TRUE
R.BS <- 500                 # number of bootstrap draws (500)
blocklength <- 12           # blocklength for autocorrelation
sample.countries <- FALSE   # TRUE if countries are to be randomly excluded, FALSE otherwise (our preference: FALSE)
min_cr <- 5                 # minimum number of pre-crisis periods in the bootstrap
alphavec <- c(0.05,0.95)    # quantiles to be plotted as confidence bands  c(0.05,0.95)  c(0.1,0.9)


##################################
# parameter settings
parameters$optimizationMode <- FALSE
parameters$paramtable <- read.csv(paramtable_file,row.names="method")
for (col in grep("_val",colnames(parameters$paramtable))){
  parameters$paramtable[,col] <- as.numeric(as.character(parameters$paramtable[,col]))
}
for (col in grep("_name",colnames(parameters$paramtable))){
  parameters$paramtable[,col] <- as.character(parameters$paramtable[,col])
}
