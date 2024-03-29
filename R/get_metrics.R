
dist <- function(X)
{
  X2 <- apply (X, 1, function(z)sum(z^2))
  dist2 <- outer(X2, X2, "+")- 2*X%*%t(X)
  dist2 <- ifelse(dist2<0,0,dist2)
  return (sqrt(dist2))
}
#' @export
get_knn_res <- function(X1,num_groups,group_sizes){
  groups <- list()
  for(i in 1:num_groups){
    tmp <- matrix(1,nrow=group_sizes,ncol=group_sizes)
    groups[[i]] <- tmp
  }
  knn_mat <- Reduce(adiag,groups)
  dist_mat <- dist(X1)
  rank_mat <- apply(dist_mat,1,rank)
  knn1 <- ifelse(rank_mat==2,1,0)
  percent <- sum(sum(knn1 * knn_mat))/as.numeric(num_groups*group_sizes)
  return(percent)
}
#' @export
spearman_rho <- function(X1,X2){
  stats::cor(as.vector(dist(X1)),
             as.vector(dist(X2)),
             method="spearman")
}
#' @export
lcmc <- function(X1,X2,k){
  #'
  #'   X1 is the input matrix
  #'   X2 is the output matrix
  #'   k is the number of neighbors to consider
  #'
  n <- dim(X1)[1]
  Do  <- dist(X1)
  Do2 <- dist(X2)
  Daux  <- apply(Do,2,sort)[k+1,]
  Daux2 <- apply(Do2,2,sort)[k+1,]
  Inb  <- ifelse(Do>c(Daux),0,1)
  Inb2 <- ifelse(Do2>c(Daux2),0,1)
  percent.neighbors.same <- (sum(Inb*Inb2)-n)/(k*n) - k/(n-1)
  return(percent.neighbors.same)
}
#' @export
kmax <- function(X1,X2){
  vals <- apply(matrix(c(1:(nrow(X1)-2)),
                       nrow=1),
                1,
                function(x) lcmc(X1,X2,x))
  max(vals)
}
#' @export
qnx <- function(X1,X2,ks=NULL,kt=1){
  #'
  #'   X1 is the input matrix
  #'   X2 is the output matrix
  #'   k is the number of neighbors to consider
  #'
  n <- dim(X1)[1]
  if(is.null(ks)){
    ks = n-1
  }
  Do <- dist(X1)
  Do2 <- dist(X2)
  Daux <- apply(Do,2,sort)[ks+1,]
  Daux2 <- apply(Do2,2,sort)[ks+1,]

  rank1 <- apply(Do,2,rank)
  rank2 <- apply(Do2,2,rank)

  abs_diff <- abs(rank1-rank2)
  wt <- ifelse(abs_diff>(kt+1),0,1)

  rank1t <- apply(t(Do),2,rank)
  rank2t <- apply(t(Do2),2,rank)
  abs_difft <- abs(rank1t-rank2t)
  wtt <- ifelse(abs_difft>(kt+1),0,1)
  Inb <- ifelse(Do>Daux,0,1)
  Inb2 <- ifelse(Do2>Daux2,0,1)
  percent.neighbors.same <- (sum(Inb*Inb2*(wt+wtt)))/(2*ks*n)
  return(percent.neighbors.same)
}
#' @export
spectral <- function(X,Y,k){

  #' Since we consistently want to see a spectrum of
  #' how many nearest neighbors are kept, the focus
  #' of this function is to capture 1 to k of the
  #' percent of nearest neighbors kept.

  n <- dim(X)[1]

  distances1 <- dist(X)
  distances2 <- dist(Y)

  rank1 <- apply(distances1,2,sort)
  rank2 <- apply(distances2,2,sort)

  kept <- matrix(0,nrow=(k-1),ncol=1)

  for(i in 1:(k-1)){
    neighborhood1 <- ifelse(distances1>rank1[(i+1),],0,1)
    neighborhood2 <- ifelse(distances2>rank2[(i+1),],0,1)
    percent_neighbors_same <- (sum(neighborhood1*neighborhood2)-n)/(i*n)
    kept[i] <- percent_neighbors_same
  }
  return(sum(kept))
}

#' @export
local_error <- function(high,low){
  N = dim(high)[1]
  high_distances = dist(high)
  low_distances = dist(low)

  high_ranks = apply(high_distances,2,rank)-1
  low_ranks = apply(low_distances,2,rank)-1

  errors <- c()
  for(i in 1:(N-1)){
    err <- ifelse(high_ranks==i,1,0)*(high_distances - low_distances)**2
    mean_err <- mean(err[err!=0])
    errors <- c(errors,mean_err)
  }
  return(sum(cumsum(errors)))
}

#' @export
entropy_and_mi <- function(X1,X2){
  n = dim(X1)[1]
  mat <- coRanking::coranking(dist(X1),dist(X2))
  pij <- mat/(n-1)

  H <- - sum(sum(pij[pij!=0] * log(pij[pij!=0])))

  divisor <- pij %*% pij
  divisor <- ifelse(divisor!=0,divisor,1)
  pij_mi <- pij/divisor
  MI <- sum(sum(pij[pij!=0] * log(pij_mi[pij_mi!=0])))
  output <- list("entropy"=H,
                 "mutual_info"=MI)
  return(output)
}
#' @export
my_procrust_dist <- function(original,output,num_rotations = 10){
  phis <- seq(0,2*pi,2*pi/num_rotations)
  distances <- c()
  for(i in 1:num_rotations){
    rot_mat <- matrix(c(cos(phis[i]),sin(phis[i]),
                        -sin(phis[i]),cos(phis[i])),
                      byrow = TRUE,ncol=2,nrow=2)
    rotated <- output %*% rot_mat
    mu3 <- CVXR::Variable(1)
    ttr3 <- CVXR::Variable(1)
    mu4 <- CVXR::Variable(1)
    ttr4 <- CVXR::Variable(1)
    objective <- mean((original[,1] - ttr3 * rotated[,1]  - mu3)**2 +(original[,2] - ttr4 * rotated[,2]  - mu4)**2)
    prob <- Problem(Minimize(objective))
    result <- solve(prob)
    a <- result$getValue(mu3)
    b <- result$getValue(ttr3)
    c <- result$getValue(mu4)
    d <- result$getValue(ttr4)

    rotated[,1] <- b*rotated[,1] - a
    rotated[,2] <- d*rotated[,2] - c

    new_dist <- CVXR::norm(original-rotated,"2")

    distances <- c(distances,new_dist)
  }
  return(min(distances))

}
