#!/usr/bin/env Rscript
i = as.numeric(commandArgs(trailingOnly=TRUE))
sink(paste0("~/kaust/log_dataset",i,".txt"))

library(RandomFields)
RFoptions(pivot=PIVOT_NONE, critical=-1, reoptimise=F, optimiser="optim", cores = 1,printlevel = 7)


  df<-read.csv(paste0("~/kaust/Sub-competition_1a/dataset",i,"_training.csv"))
  #df <- df[1:18e3,]  
  len <- length(df$x)
  split_index <- sapply(1:9, function(n)rep(n,len/9))
  split_df <- split(df[,1:2],split_index)
  split_val <- split(df[,3],split_index)
  print(system.time(model <- RFfit(RMmatern(nu=NA,scale=NA,var=NA) + RMnugget(var=NA),
                                   x= split_df[1:4],data=split_val[1:4], sub.methods="plain",optim.control=list(trace=6,maxit=10,lmm=3),spC=F)))
  saveRDS(model, file = paste0("~/kaust/models/dataset_update_",i,".RDS"))
  cat(timestamp())
  cat(rep("\n",10))
sink(NULL)
