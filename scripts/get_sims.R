

double_X <- function(n=250){
  num <- n/2
  x <- seq(-2.114743,2.114743,2*20^(1/3)/(num-1))
  y1 <- 50 + sin(pi* 10 * x) * x**4
  y2 <- -50 - sin(pi * 10 * x) * x**4
  X <- matrix(c(x,x),ncol=1)
  Y <- matrix(c(y1,y2),ncol=1)
  rotation_mat <- matrix(c(cos(pi/4),sin(pi/4),
                           -sin(pi/4),cos(pi/4)),
                         nrow=2)
  cbind(X,Y) %*% rotation_mat
}

generate_three_gaussians <- function(n=250){
  m <- ceiling(n/3)
  rbind(MASS::mvrnorm(m,mu=c(20,0),Sigma = 10*diag(2)),
        MASS::mvrnorm(m,mu=c(-10,0),Sigma = diag(2)),
        MASS::mvrnorm(m,mu=c(-5,0),Sigma = diag(2)))
}
two_lines <- function(n=250){
  m <- ceiling(n/2)
  x <- rnorm(n,mean=0,sd=1)
  y <- 2 * x + c(rep(-5,m),rep(5,n-m))
  cbind(x,y)
}

circles <- function(n=250){
  phi <- seq(0,2*pi,4*pi/n)
  x <- c(cos(phi),0.5*cos(phi)) + rnorm(n=length(phi),sd=0.02)
  y <- c(sin(phi),0.5*sin(phi)) + rnorm(n=length(phi),sd=0.02)
  cbind(x,y)
}


trefoil <- function(n=250){
  phi <- seq(0,2*pi,2*pi/n)
  x <- sin(phi) + 2*sin(2*phi) + rnorm(n=length(phi),sd=0.1)
  y <- cos(phi) - 2*cos(2*phi) + rnorm(n=length(phi),sd=0.1)
  cbind(x,y)
}

high_dim_clusters <- function(){
  mean_vals <- as.matrix(expand.grid(
    v1 <- c(5,0),
    v2 <- c(5,0),
    v3 <- c(5,0),
    v4 <- c(5,0)
  ))
  df <- c()
  for(i in 1:16){
    tmp <- MASS::mvrnorm(50,
                         mu=mean_vals[i,],
                         Sigma=diag(4))
    df <- rbind(df,tmp)
  }
  df
}
