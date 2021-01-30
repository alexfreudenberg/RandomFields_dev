library(data.table)
library(RandomFields)
RFoptions(pivot=PIVOT_NONE, cores=8)

RFoptions(useGPU=TRUE)
alpha <- 10
teamname <- "Mannheim"

for(i in 1:16) {
  message(paste0("Preparing Dataset ", i, "\n"))
  training <- fread(paste0("~/kaust/Sub-competition_1a/dataset",i,"_training.csv"))
  testing <- fread(paste0("~/kaust/Sub-competition_1b/dataset",i,"_testing.csv"))
  testing$values <- NA_real_
  fit <- readRDS(paste0("kaust/dataset_update2_", i, ".RDS"))
  
  #str(training)
  #str(testing)
  #n_training <- 1e4
  #n_testing <- 2e4
  #sampled_rows_training <- sample(seq_len(nrow(training)), n_training)
  #sampled_rows_testing <- sample(seq_len(nrow(testing)), n_testing)
  
  
  for(k in 1:alpha) {
    message(paste0("Dataset ", i, ", loop ", k, "of", alpha, "loops.\n"))
    training_selection <- (((k-1)/alpha - 0.1 * 1/alpha) < training$x)  &   ((k/alpha + 0.1 * 1/alpha) > training$x)
    testing_selection <- (((k-1)/alpha) <= testing$x)  &   ((k/alpha) > testing$x)
    
    n_training <- sum(training_selection)
    n_testing <- sum(testing_selection)
    sampled_rows_training <- seq_len(nrow(training))[training_selection]
    sampled_rows_testing <- seq_len(nrow(testing))[testing_selection]
                                                   
    locations <- rbind(training[sampled_rows_training,.(x,y)], testing[sampled_rows_testing,.(x,y)])
    conditional_known <- training[sampled_rows_training,]$values
    
    known_rows <- seq_len(n_training)
    prediction_rows <- seq_len(nrow(locations))[-known_rows]
    
    system.time(comatrix <- RFcovmatrix(fit, x=data.frame(locations)))
    
    system.time(pred <- as.vector(comatrix[prediction_rows, known_rows] %*%
                        solvex(comatrix[known_rows, known_rows],
                               conditional_known)))
    
    testing[testing_selection==TRUE,values:=pred]
  }
    
  colnames(testing)[3] <- "predicted values"
  write.csv(x=testing, file=paste0("kaust/delivery/", teamname, "-1b-", i,".csv"), row.names=FALSE, quote=FALSE)
}

if(FALSE) {
  RFoptions(solve_method = "cholesky", printlevel=1,useGPU=T)
  set.seed(1)
  n <- 10000
  x <- 1:n
  y <- runif(n)
  
  
  ## FIRST EXAMPLE: full rank matrix
  M <- exp(-as.matrix(dist(x) / n)) 
  b0 <- matrix(nr=n, runif(n * 5))
  b <- M %*% b0 + runif(n)
  
  ## Without pivoting:
  RFoptions(pivot=PIVOT_NONE) ## (default)
  print(system.time(z <- solvex(M, b)))
  print(range(b - M %*% z))
  stopifnot(all(abs((b - M %*% z)) < 2e-11))
}