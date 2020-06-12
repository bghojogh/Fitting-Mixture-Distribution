rm(list=ls(all=TRUE)) # remove history

#===================== function:
fit_several_Poissons <- function(x, K, max.iter=100, eps=1e-7, print_steps=FALSE){ 
  # x: data
  # K: number of distributions
  n <- length(x)
  
  # initialize parameters:
  lambda_new <- runif(K, min(x), max(x))
  w_new <- runif(K, 0, 1)
  lambda_of_iterations <- matrix(0, nrow=max.iter, ncol=K)
  lambda_of_iterations[1,] <- lambda_new
  w_of_iterations <- matrix(0, nrow=max.iter, ncol=K)
  w_of_iterations[1,] <- w_new
  
  iteration_index = 0; convg<-FALSE
  while (iteration_index <= max.iter && !convg){ 
    # update parameters:
    lambda <- lambda_new 
    w <- w_new
    if (print_steps == TRUE) {
      cat("Iteration: ", iteration_index, "\n")
      cat("lambda = ", lambda, "\n") 
      cat("pi = ", Pi, "\n")
      print("=======")           
    }
    cat("Iteration: ", iteration_index, "\n")
    
    # e-step  
    gamma = matrix(0, nrow=n, ncol=K)
    for (i in 1:n) {
      g = exp(-lambda) * lambda^x[i] / factorial(x[i])  # vector of g's for all k's
      for (k in 1:K) {
        gamma[i,k] = (w[k] * g[k]) / (sum(w * g))
      }
    }
    
    # m-step    
    lambda_new <- vector(length = K)
    w_new <- vector(length = K)
    for (k in 1:K) {
      lambda_new[k] = sum(gamma[,k]*x) / sum(gamma[,k])
      w_new[k] = (sum(gamma[,k])) / (sum(gamma))
    }
    print(w_new)
    
    # save the obtained parameters in every iteration:
    for (k in 1:K) {
      lambda_of_iterations[iteration_index,k] <- lambda_new[k]
      w_of_iterations[iteration_index,k] <- w_new[k]
    }
    
    #convg <- abs(lambda_new-lambda)<eps & abs(w_new-w)<eps
    iteration_index <- iteration_index + 1   
  }   
  return(list("lambda"=lambda_new, "w"=w_new, "gamma"=gamma, "w_of_iterations"=w_of_iterations, "lambda_of_iterations"=lambda_of_iterations)) 
} 


#===================== data:
frequencies <- c(162, 267, 271, 185, 111, 61, 120, 210, 215, 136, 73, 43, 14, 160, 230, 243, 104, 36, 15, 10, 0)
n <- length(frequencies)
x <- vector(length = 0)
for (i in 1:length(frequencies)) {
  x <- c(x, rep(i-1,frequencies[i]))
}
print(x)

#===================== plot data:
barplot(names.arg = 0:(n-1), frequencies, xlab="x", ylab="frequency")

#===================== call fit function and print results:
optimum_parameters <- fit_several_Poissons(x=x, K=3) 
w <- optimum_parameters$w
lambda <- optimum_parameters$lambda
K <- length(optimum_parameters$w)
w_of_iterations <- optimum_parameters$w_of_iterations
lambda_of_iterations <- optimum_parameters$lambda_of_iterations
cat("w = ", w, "\n")
cat("lambda = ", lambda, "\n")

#===================== plot the progress of fitted weights:
# first an empty plot just for showing axes:
plot(1,type='n',xlim=c(1,length(w_of_iterations[,1])),ylim=c(min(w_of_iterations),max(w_of_iterations)),xlab='iteration', ylab='w')
#colors = c("blue", "red", "darkgreen")
colors = c("darkgreen", "red", "blue")
sumbols = c(1, 22, 17)
for (k in 1:K) {
  lines(x=1:length(w_of_iterations[,k]), y=w_of_iterations[,k], type="o", pch=sumbols[k], col=colors[k], lwd=2)
}
grid (NULL,NULL, lty = 6, col = "cornsilk2") 

#===================== plot the progress of fitted lambdas:
# first an empty plot just for showing axes:
plot(1,type='n',xlim=c(1,length(lambda_of_iterations[,1])),ylim=c(min(lambda_of_iterations),max(lambda_of_iterations)),xlab='iteration', ylab=expression(lambda))
#colors = c("blue", "red", "darkgreen")
colors = c("darkgreen", "red", "blue")
sumbols = c(1, 22, 17)
for (k in 1:K) {
  lines(x=1:length(lambda_of_iterations[,k]), y=lambda_of_iterations[,k], type="o", pch=sumbols[k], col=colors[k], lwd=2)
}
grid (NULL,NULL, lty = 6, col = "cornsilk2") 

#===================== plot total fitted distribution:
# first an empty plot just for showing axes:
plot(1,type='n',xlim=c(0,20),ylim=c(0,0.35),xlab="x", ylab="mass")
f <- vector(length = n)
for (y in 0:(n-1)) {
  for (k in 1:K) {
    g = exp(-lambda[k]) * lambda[k]^y / factorial(y)
    f[y+1] <- f[y+1] + (w[k] * g)
  }
}
grid (NULL,NULL, lty = 6, col = "cornsilk2") 
lines(x=0:(n-1), y=f, type="o", col="mediumorchid1", lwd=5, pch=22, ylab="mass", xlab="x")

#===================== plot K fitted distributions:
#colors = c("blue", "red", "darkgreen")
colors = c("darkgreen", "red", "blue")
for (k in 1:K) {
  g <- vector(length = n)
  for (y in 0:(n-1)) {
    g[y+1] = exp(-lambda[k]) * lambda[k]^y / factorial(y)
  }
  lines(x=0:(n-1), y=g, type="o", pch=22, col=colors[k], lwd=5)
}

#===================== plot only one fitted distribution:
lambda_onlyOne = mean(x)
f_onlyOne <- vector(length = n)
for (y in 0:(n-1)) {
  f_onlyOne[y+1] = exp(-lambda_onlyOne) * lambda_onlyOne^y / factorial(y)
}
lines(x=0:(n-1), y=f_onlyOne, type="o", pch=22, col="burlywood", lwd=5)


# save workspace:
save(list=ls(all=TRUE), file="./workspace.RData")

# load workspace:
load(file="./workspace.RData")




