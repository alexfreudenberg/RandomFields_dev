teamname <- "Mannheim"

# TODO: Add Hardware Description
hardware <- "TODO"

options(scipen = 999)

params_1a <- data.frame(t(sapply(
  1:16, 
  function(i) {
    fit <- readRDS(paste0("kaust/dataset_update2_", i, ".RDS"))
    as.character(round(fit$ml$param["value",],6))
  }
)))


# TODO: add times
params_1a$time <- NA_real_

colnames(params_1a) <- c("sigma squared", "beta", "nu", "tau squared", "time in seconds")

write.csv(x=params_1a, file=paste0("kaust/delivery/", teamname, "-1a.csv"), row.names=FALSE, quote=FALSE)
cat(hardware, file=paste0("kaust/delivery/", teamname, "-1a-hardware.txt"))
