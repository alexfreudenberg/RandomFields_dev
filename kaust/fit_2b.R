i <- as.numeric(commandArgs(trailingOnly=TRUE))
sink(paste0("~/kaust/2b_log-",i,".txt"))


library(RandomFields)
RFoptions(pivot=PIVOT_NONE, cores=3, useGPU=TRUE, sub_optimiser="soma", optimiser="optim", print=7)

  training <- read.csv(paste0("~/kaust/Sub-competition_2b/dataset",i,"_training.csv"))
  #str(training)

  #n_training <- 2e4
  #sampled_rows_training <- sample(seq_len(nrow(training)), n_training)
  #system.time(vario <- RFvariogram(data=training[sampled_rows_training,]))
  #plot(vario)
  
  model <- ~ RMspheric(var=var1, scale=scale1) + RMspheric(var=var2, scale=scale2) + RMdeclare(scaleplus=scaleplus)
  param <- list(var1=NA, scale1=NA, var2=NA,  scaleplus=NA, scale2=~scale1+scaleplus)
  lower <- list(var1=1e-1, var2=1e-1, scale1=1e-6, scaleplus=0.01)
  upper <- list(var1=20, var2=20, scale1=1, scaleplus=1)
  
  #fit1 <- RFfit(model, param=param, lower=lower, upper=upper, data=training[sampled_rows_training,], method=NULL)
  
  #fit1_parameter <- unname(fit1@self@param["value",])
  #users.guess <- list(var1=fit1_parameter[1], scale1=fit1_parameter[2], var2=fit1_parameter[3], scaleplus=fit1_parameter[5])
  
  
  len <- length(training$x)
  split_index <- 10 * floor(10*training$x) + floor(10*training$y)
  split_df <- split(training[,1:2],split_index)
  split_val <- split(training[,3],split_index)
  evalsubset <- c(1, 30, 61, 99)
  start <- Sys.time()
  fit2 <- RFfit(#model, param=param, lower=lower, upper=upper, 
    RMmatern(nu=NA, var=NA, scale=NA) + RMnugget(),
    x= split_df[evalsubset],
    data=split_val[evalsubset], 
    sub.methods="plain",
    optim.control=list(trace=6,lmm=3))
  
  cat(difftime(Sys.time(), start, unit="secs"), file=paste0("kaust/2b_time-", i, ".txt"))
  saveRDS(fit2, file=paste0("kaust/2b_fit-", i, ".RDS"))
