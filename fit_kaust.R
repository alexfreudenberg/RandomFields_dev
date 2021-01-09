library(RandomFields)
library(optimParallel)
RFoptions(pivot=PIVOT_NONE,printlevel=6, optimiser="lbfgs", cores = 12,printlevel = 6)
cl <- makeCluster(9)     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster


for(i in 1:16){
  df<-read.csv(paste0("~/Documents/kaust/Sub-competition_1a/dataset",i,"_training.csv"))
  df <- df[1:18e3,]  
  len <- length(df$x)
  split_index <- sapply(1:9, function(n)rep(n,len/9))
  split_df <- split(df[,1:2],split_index)
  split_val <- split(df[,3],split_index)
  print(system.time(model <- RFfit(RMmatern(nu=NA,scale=NA,var=NA) + RMnugget(var=NA),
                                   x= split_df,data=split_val, sub.methods="plain",optim.control=list(trace=6,maxit=10,lmm=3),spC=F)))
  saveRDS(model, file = paste0("kaust/models/dataset_",i,".RDS"))
  cat(timestamp())
  cat(rep("\n",10))
}
