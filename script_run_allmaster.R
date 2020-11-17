# This script runs cross-validation, in-sample and out-of-sample estimations
# Gregor von Schweinitz, 2020-11-17

cat("\014")     # clear console
rm(list = ls()) # clear workspace

list.masterfiles_cv <- c("mu05_cv")
list.masterfiles_est <- c("mu05_is","mu05_rec")
ls.keep <- c(ls(),"ls.keep")

for (masterfile in list.masterfiles_cv){
  # For safety reasons, the output of the cross-validation ("outtable_05.csv") is not saved in the main working directory, but in the output working directory
  # For usage, this needs to be copied into the main working directory
  print("#########################")
  print("#########################")
  print("")
  print("")
  print("")
  print(masterfile)
  print("")
  print("")
  print("")
  print("#########################")
  print("#########################")
  source(paste("script_master_",masterfile,".R",sep=""))
  rm(list=setdiff(ls(),ls.keep))
}

for (masterfile in list.masterfiles_est){
  print("#########################")
  print("#########################")
  print("")
  print("")
  print("")
  print(masterfile)
  print("")
  print("")
  print("")
  print("#########################")
  print("#########################")
  source(paste("script_master_",masterfile,".R",sep=""))
  rm(list=setdiff(ls(),ls.keep))
}