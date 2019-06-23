
evaluate_output <- function(high,low){
  n = dim(high)[1]
  kmax_val <- kmax(high,low)
  qnx_val <- qnx(high,low)
  em <- entropy_and_mi(high,low)
  h <- em$entropy
  mi <- em$mutual_info
  le <- local_error(high,low)
  rk <- tryCatch({range_kept(high,low,n)},error=function(e){return(rankge_kept(high,low,(n-1)))})
  spear <- spearman_rho(high,low)
  knn1 <- get_knn_res(low,num_groups = 16,group_sizes = 50)
  data_row <- matrix(c(kmax_val,qnx_val,h,mi,le,rk,spear,knn1),nrow=1)
  return(data.frame(data_row))
}
