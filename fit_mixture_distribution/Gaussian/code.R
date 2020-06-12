rm(list=ls(all=TRUE)) # remove history

#===================== function:
fit_several_Gaussians <- function(x, K, max.iter=150, eps=1e-7, print_steps=FALSE){ 
  # x: data
  # K: number of distributions
  n <- length(x)
  
  # initialize parameters:
  mu_new <- runif(K, min(x), max(x))
  sigma_new <- runif(K, 0, (max(x) - min(x))/6)
  w_new <- runif(K, 0, 1)
  mu_of_iterations <- matrix(0, nrow=max.iter, ncol=K)
  mu_of_iterations[1,] <- mu_new
  sigma_of_iterations <- matrix(0, nrow=max.iter, ncol=K)
  sigma_of_iterations[1,] <- sigma_new
  w_of_iterations <- matrix(0, nrow=max.iter, ncol=K)
  w_of_iterations[1,] <- w_new
  
  iteration_index = 0; convg <- FALSE
  while (iteration_index <= max.iter && !convg){ 
    # update parameters:
    mu <- mu_new 
    sigma <- sigma_new 
    w <- w_new
    if (print_steps == TRUE) {
      cat("Iteration: ", iteration_index, "\n")
      cat("mu = ", mu, "\n")
      cat("sigma = ", sigma, "\n")
      cat("w = ", w, "\n")
      print("=======")           
    }
    cat("Iteration: ", iteration_index, "\n")
    
    # e-step  
    gamma = matrix(0, nrow=n, ncol=K)
    for (i in 1:n) {
      g = 1/(sqrt(2*pi*sigma^2)) * exp(-0.5*((x[i]-mu)^2/(sigma^2)))  # vector of g's for all k's
      for (k in 1:K) {
        gamma[i,k] = (w[k] * g[k]) / (sum(w * g))
      }
    }
    
    # m-step    
    mu_new <- vector(length = K)
    sigma_new <- vector(length = K)
    w_new <- vector(length = K)
    for (k in 1:K) {
      mu_new[k] = sum(gamma[,k]*x) / sum(gamma[,k])
      sigma_new[k] = sqrt(sum(gamma[,k]*(x-mu_new[k])^2) / sum(gamma[,k]))
      w_new[k] = (sum(gamma[,k])) / (sum(gamma))
    }
    
    # save the obtained parameters in every iteration:
    for (k in 1:K) {
      mu_of_iterations[iteration_index,k] <- mu_new[k]
      sigma_of_iterations[iteration_index,k] <- sigma_new[k]
      w_of_iterations[iteration_index,k] <- w_new[k]
    }
    
    #convg <- abs(mu_new-mu)<eps & abs(w_new-w)<eps
    iteration_index <- iteration_index + 1   
  }   
  return(list("mu"=mu_new, "sigma"=sigma_new, "w"=w_new, "gamma"=gamma, "w_of_iterations"=w_of_iterations, "mu_of_iterations"=mu_of_iterations, "sigma_of_iterations"=sigma_of_iterations)) 
} 


#===================== data:
n1 <- 700; mu1 <- -10; sigma1 <- 1.2
n2 <- 1000; mu2 <- 0;   sigma2 <- 2
n3 <- 500; mu3 <- 5;   sigma3 <- 5
n <- n1 + n2 + n3
x1 <- rnorm(n1, mean = mu1, sd = sigma1)
x2 <- rnorm(n2, mean = mu2, sd = sigma2)
x3 <- rnorm(n3, mean = mu3, sd = sigma3)
x <- c(x1, x2, x3)
print(x)

#===================== plot the original densities of data:
colors = c("blue", "red", "darkgreen")
x_sweep_list <- seq(-25, 20, by = 0.1) 
# plot first density:
mu <- mu1; sigma <- sigma1
g_original <- vector(length = length(x_sweep_list))
index <- 0
for (x_sweep in x_sweep_list) {
  index <- index + 1
  g_original[index] <- 1/(sqrt(2*pi*sigma^2)) * exp(-0.5*((x_sweep-mu)^2/(sigma^2)))
}
max_g1_density <- max(g_original)
plot(x=x_sweep_list, y=g_original, type="l", col=colors[1], lwd=5, ylab="density", xlab="x", xlim=c(-20,20))
# plot second density:
mu <- mu2; sigma <- sigma2
g_original <- vector(length = length(x_sweep_list))
index <- 0
for (x_sweep in x_sweep_list) {
  index <- index + 1
  g_original[index] <- 1/(sqrt(2*pi*sigma^2)) * exp(-0.5*((x_sweep-mu)^2/(sigma^2)))
}
max_g2_density <- max(g_original)
lines(x=x_sweep_list, y=g_original, type="l", col=colors[2], lwd=5, ylab="density", xlab="x")
# plot third density:
mu <- mu3; sigma <- sigma3
g_original <- vector(length = length(x_sweep_list))
index <- 0
for (x_sweep in x_sweep_list) {
  index <- index + 1
  g_original[index] <- 1/(sqrt(2*pi*sigma^2)) * exp(-0.5*((x_sweep-mu)^2/(sigma^2)))
}
max_g3_density <- max(g_original)
lines(x=x_sweep_list, y=g_original, type="l", col=colors[3], lwd=5, ylab="density", xlab="x")

#===================== plot the total original density of data:
g_original <- vector(length = length(x_sweep_list))
index <- 0
for (x_sweep in x_sweep_list) {
  index <- index + 1
  temp <- 0
  for (k in 1:3) {
    if (k==1) {
      mu <- mu1; sigma <- sigma1; 
    } else if (k==2) {
      mu <- mu2; sigma <- sigma2; 
    } else if (k==3) {
      mu <- mu3; sigma <- sigma3; 
    }
    temp <- temp + 1/(sqrt(2*pi*sigma^2)) * exp(-0.5*((x_sweep-mu)^2/(sigma^2)))
  }
  g_original[index] <- temp
}
max_g1_density <- max(g_original)
lines(x=x_sweep_list, y=g_original, type="l", col="blue", lwd=2, ylab="density", xlab="x", xlim=c(-20,20))

#===================== plot the original data:
# par(new=T)
# max_densities <- max(c(max_g1_density, max_g2_density, max_g3_density))
# plot(x=x1, y=runif(length(x1), 0, max_densities/20), col="blue", lwd=2, xlab="", ylab="", ylim=c(0,max_densities), axes=FALSE)
# par(new=T)
# plot(x=x2, y=runif(length(x2), 0, max_densities/20), col="red", lwd=2, xlab="", ylab="", ylim=c(0,max_densities), axes=FALSE)
# par(new=T)
# plot(x=x3, y=runif(length(x3), 0, max_densities/20), col="darkgreen", lwd=2, xlab="", ylab="", ylim=c(0,max_densities), axes=FALSE)
# 


#===================== call fit function and print results:
optimum_parameters <- fit_several_Gaussians(x=x, K=3) 
w <- optimum_parameters$w
mu <- optimum_parameters$mu
sigma <- optimum_parameters$sigma
K <- length(optimum_parameters$w)
w_of_iterations <- optimum_parameters$w_of_iterations
mu_of_iterations <- optimum_parameters$mu_of_iterations
sigma_of_iterations <- optimum_parameters$sigma_of_iterations
cat("w = ", w, "\n")
cat("mu = ", mu, "\n")
cat("sigma = ", sigma, "\n")

#===================== plot the progress of fitted weights:
# first an empty plot just for showing axes:
plot(1,type='n',xlim=c(1,length(w_of_iterations[,1])),ylim=c(min(w_of_iterations),max(w_of_iterations)),xlab='iteration', ylab='w')
colors = c("blue", "red", "darkgreen")
sumbols = c(1, 22, 17)
for (k in 1:K) {
  lines(x=1:length(w_of_iterations[,k]), y=w_of_iterations[,k], type="o", pch=sumbols[k], col=colors[k], lwd=2)
}
grid (NULL,NULL, lty = 6, col = "cornsilk2") 

#===================== plot the progress of fitted mu's:
# first an empty plot just for showing axes:
plot(1,type='n',xlim=c(1,length(mu_of_iterations[,1])),ylim=c(min(mu_of_iterations),max(mu_of_iterations)),xlab='iteration', ylab=expression(mu))
colors = c("blue", "red", "darkgreen")
sumbols = c(1, 22, 17)
for (k in 1:K) {
  lines(x=1:length(mu_of_iterations[,k]), y=mu_of_iterations[,k], type="o", pch=sumbols[k], col=colors[k], lwd=2)
}
grid (NULL,NULL, lty = 6, col = "cornsilk2") 

#===================== plot the progress of fitted sigma's:
# first an empty plot just for showing axes:
plot(1,type='n',xlim=c(1,length(sigma_of_iterations[,1])),ylim=c(min(sigma_of_iterations),max(sigma_of_iterations)),xlab='iteration', ylab=expression(sigma))
colors = c("blue", "red", "darkgreen")
sumbols = c(1, 22, 17)
for (k in 1:K) {
  lines(x=1:length(sigma_of_iterations[,k]), y=sigma_of_iterations[,k], type="o", pch=sumbols[k], col=colors[k], lwd=2)
}
grid (NULL,NULL, lty = 6, col = "cornsilk2") 

#===================== plot total fitted distribution:
# first an empty plot just for showing axes:
plot(1,type='n',xlim=c(-20,20),ylim=c(0,0.35),xlab="x", ylab="density")
g_original <- vector(length = length(x_sweep_list))
index <- 0
for (x_sweep in x_sweep_list) {
  index <- index + 1
  temp <- 0
  for (k in 1:K) {
    temp <- temp + (w[k] * (1/(sqrt(2*pi*sigma[k]^2)) * exp(-0.5*((x_sweep-mu[k])^2/(sigma[k]^2)))))
  }
  g_original[index] <- temp
}
max_g1_density <- max(g_original)
lines(x=x_sweep_list, y=g_original, type="l", col="mediumorchid1", lty=3, lwd=5, ylab="density", xlab="x", xlim=c(-20,20))


#===================== plot K fitted distributions:
colors = c("blue", "red", "darkgreen")
for (k in 1:K) {
  g <- vector(length = length(x_sweep_list))
  index <- 0
  for (x_sweep in x_sweep_list) {
    index <- index + 1
    g[index] <- 1/(sqrt(2*pi*sigma[k]^2)) * exp(-0.5*((x_sweep-mu[k])^2/(sigma[k]^2)))
  }
  lines(x=x_sweep_list, y=g, type="l", col=colors[k], lwd=5, ylab="density", xlab="x", xlim=c(-20,20))
}

#===================== plot only one fitted distribution:
mu_onlyOne = mean(x)
sigma_onlyOne = sd(x)
f_onlyOne <- vector(length = length(x_sweep_list))
index <- 0
for (x_sweep in x_sweep_list) {
  index <- index + 1
  f_onlyOne[index] <- 1/(sqrt(2*pi*sigma_onlyOne^2)) * exp(-0.5*((x_sweep-mu_onlyOne)^2/(sigma_onlyOne^2)))
}
lines(x=x_sweep_list, y=f_onlyOne, type="l", col="burlywood", lty=2, lwd=5, ylab="density", xlab="x", xlim=c(-20,20))


# save workspace:
save(list=ls(all=TRUE), file="./workspace.RData")

# load workspace:
load(file="./workspace.RData")




