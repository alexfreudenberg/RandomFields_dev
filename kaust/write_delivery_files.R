teamname <- "TeamMS"

# TODO: Add Hardware Description
hardware <- "We trained 8 models in parallel at a time on an workstation equipped with an Intel Xeon E5-1630, 64 GB memory and an Nvidia RTX 2080 Ti GPU."

options(scipen = 999)
########################### Competition 1a #####################################
params_1a <- data.frame(t(sapply(
  1:16, 
  function(i) {
    fit <- readRDS(paste0("kaust/dataset_update2_", i, ".RDS"))
    as.character(round(fit$ml$param["value",],6))
  }
)))


params_1a$time <- c(23840.735, 19766.53,16084.675,11374.494, 21186.153,18986.247,19884.400,
17390.325, 19219.877,19274.823,13827.237,21870.02,13163.113,12717.330,17289.590,13270.578)

colnames(params_1a) <- c("sigma squared", "beta", "nu", "tau squared", "time in seconds")

write.csv(x=params_1a, file=paste0("kaust/delivery/", teamname, "-1a.csv"), row.names=FALSE, quote=FALSE)
cat(hardware, file=paste0("kaust/delivery/", teamname, "-1a-hardware.txt"))
cat(hardware, file=paste0("kaust/delivery/", teamname, "-1b-hardware.txt"))


########################### Competition 2a #####################################
cat(hardware, file=paste0("kaust/delivery/", teamname, "-2a-hardware.txt"))
cat(hardware, file=paste0("kaust/delivery/", teamname, "-2b-hardware.txt"))