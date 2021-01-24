library(RandomFields)
RFoptions(pivot=PIVOT_NONE, critical=-1, reoptimise=F, optimiser="optimParallel", cores = 12,printlevel = 7)
model <- RMmatern(nu=1,var=1,scale=1)
RFoptions(cores=12, useGPU=T,pivot=PIVOT_NONE)
set.seed(0)
repet <- 1
pts <- 4e3
x <- runif(n=pts, min=-1, max=1)
y <- runif(n=pts, min=-1, max=1)
dta <- RFsimulate(model, x=x, y=y, n=repet, spC = FALSE)

library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

system.time(fit <- RFfit(RMmatern(nu=NA,var=NA,scale=NA)+RMnugget(var=NA),x=x,y=y,data=dta))


for(i in 12:16){
  df<-read.csv(paste0("~/kaust/Sub-competition_1a/dataset",i,"_training.csv"))
  #df <- df[1:18e3,]  
  len <- length(df$x)
  split_index <- sapply(1:9, function(n)rep(n,len/9))
  split_df <- split(df[,1:2],split_index)
  split_val <- split(df[,3],split_index)
  print(system.time(model <- RFfit(RMmatern(nu=NA,scale=NA,var=NA) + RMnugget(var=NA),
                                   x= split_df[1:4],data=split_val[1:4], sub.methods="plain",optim.control=list(trace=6,maxit=10,lmm=3),spC=F)))
  saveRDS(model, file = paste0("~/kaust/models/dataset_",i,".RDS"))
  cat(timestamp())
  cat(rep("\n",10))
}
