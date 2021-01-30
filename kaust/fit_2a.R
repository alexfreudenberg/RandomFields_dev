library(RandomFields)
RFoptions(pivot=PIVOT_NONE, cores=8, useGPU=TRUE, sub_optimiser="soma", optimiser="optim", print=7)

i <- 1
training <- read.csv(paste0("~/kaust/Sub-competition_2a/dataset",i,"_training.csv"))
str(training)

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
split_index <- sapply(1:9, function(n)rep(n,len/9))
split_df <- split(training[,1:2],split_index)
split_val <- split(training[,3],split_index)

fit2 <- RFfit(#model, param=param, lower=lower, upper=upper, 
              RMmatern(nu=NA, var=NA, scale=NA),
              x= split_df[1:4],
              data=split_val[1:4], 
              sub.methods="plain",
              optim.control=list(trace=6,lmm=3))

saveRDS(fit2, file="kaust/2a_fit.RDS")




###############
if(FALSE) {
  i <- 2
  training <- read.csv(paste0("~/kaust/Sub-competition_2a/dataset",i,"_training.csv"))
  n_training <- 1e3
  sampled_rows_training <- sample(seq_len(nrow(training)), n_training)
  system.time(vario <- RFvariogram(data=training[sampled_rows_training,]))
  plot(vario)
  
  
  model <- ~ RMspheric(var=var1, scale=scale1) + RMdampedcos(lambda = lambda, var= var2, scale=scale2)
  param <- list(var1=NA, scale1=NA, lambda=NA, var2=NA, scale2=NA)
  lower <- list(var1=1, scale1=0.1, lambda=0.1, var2=0.1, scale2=0.01)
  upper <- list(var1=100, scale1=1, lambda=0.1, var2=100, scale2=1)
  x <- seq(0, 0.7, 0.01)
  p_guess <- list(var1=20, scale1=0.3, lambda=0.2, var2=0.1, scale2=0.05)
  vario <- RFvariogram(model, param=p_guess, x=x)
  plot(x,vario, type="l")
  
  fit1 <- RFfit(model, param=param, lower=lower, upper=upper, data=training[sampled_rows_training,], method=NULL)
}