# This script runs the full in-sample estimation.
# Note that this code needs outtable_05.csv (output from cross-validation) in the working directory to set hyperparameters of machine-learning methods
# Gregor von Schweinitz, 2020-11-17

sink.N = sink.number()
if (sink.N>0){for (j in 1 : sink.N){sink()}} # close sinks

######################################
# Set working directory and input files
wd <- getwd()
wd_master <- wd
output_directory = "/baseline_mu05/is"
wd_output_detail <- paste(wd_master,output_directory,sep="")
if (!dir.exists(wd_output_detail)){dir.create(wd_output_detail,recursive=TRUE)}

data_file <- "data_random.csv"                  # data file with all variables
paramtable_file <- "outtable_05.csv"
filename <- "Data_master"     # Output filename for workspace of script_master
filename <- paste(wd_output_detail,"/",filename,sep="")

######################################
# Set up specifictions of estimation. Switch here between CV, IS and REC estimation
do.is <- TRUE
do.recursive <- FALSE
do.crossval <- FALSE
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
  source("script_estimate_par.R")
  source("script_tables.R")
}
if (do.crossval){
  stop("Run script_master_mu05_cv.R for this case")
}

if (do.parallel){stopCluster(cl)}
save.image(paste(filename,".RData",sep=""))  # save final workspace after running full script_master
toc()