library(data.table)
library(RandomFields)
RFoptions(pivot=PIVOT_NONE)

RFoptions(useGPU=FALSE)
fit <- dataset_1 

training <- fread("dataset1_training.csv")
testing <- fread("dataset1_testing.csv")

str(training)
str(testing)

n_training <- 10000
n_testing <- 10000
sampled_rows <- sample(seq_len(nrow(training)), n_training)
locations <- rbind(training[sampled_rows,.(x,y)], testing[n_testing])
conditional_known <- training[sampled_rows,]$values

known_rows <- seq_len(n_training)
prediction_rows <- seq_len(nrow(locations))[-known_rows]

system.time(comatrix <- RFcovmatrix(fit, x=data.frame(locations)))

system.time(pred <- as.vector(comatrix[prediction_rows, known_rows] %*%
                    solvex(comatrix[known_rows, known_rows],
                           conditional_known)))




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