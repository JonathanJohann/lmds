

simulate <- function(method,dataset,seed=123123,data_seed=123123){
  if(method=="Sammon") {
    X <- get_data(dataset=dataset,seed=data_seed)
    Y <- stats::sammon(dist(X))$points
    results <- evaluate(X,Y)
    return(results)
  } else if(method == "lmds") {
    X <- get_data(dataset=dataset,seed=data_seed)
    grid <- get_grid(method=method,dataset=dataset)
    results <- c()
    for(i in 1:dim(grid)[1]){
      Y <- lmds(X=X,k=grid[i,]$k,tau=grid[i,]$tau)$X
      results <- rbind(results,cbind(grid[i,],evaluate(X,Y)))
      print(paste0("Finished -- ",100.0*i/dim(grid)[1]))
    }
    return(results)
  } else if(method == "umap") {
    X <- get_data(dataset=dataset,seed=data_seed)
    grid <- get_grid(method=method,dataset=dataset)
    results <- c()
    for(i in 1:dim(grid)[1]){
      custom.settings = umap.defaults
      custom.settings$n_neighbors=grid[i,]$k
      custom.settings$local_connectivity=grid[i,]$local_connectivity
      custom.settings$bandwidth=grid[i,]$bandwidth
      Y <- umap(X,custom.settings)$layout
      results <- rbind(results,cbind(grid[i,],evaluate(X,Y)))
      print(paste0("Finished -- ",100.0*i/dim(grid)[1]))
    }
    return(results)
  } else if(method == "tsne") {
    X <- get_data(dataset=dataset,seed=data_seed)
    grid <- get_grid(method=method,dataset=dataset)
    results <- c()
    for(i in 1:dim(grid)[1]){
      Y <- Rtsne(X,perplexity=grid[i,]$perplexity,
                 pca_center=grid[i,]$pca_center,
                 pca_scale=grid[i,]$pca_scale)
      results <- rbind(results,cbind(grid[i,],evaluate(X,Y)))
      print(paste0("Finished -- ",100.0*i/dim(grid)[1]))
    }
    return(results)
  } else {
    stop("Please provide a valid method.")
  }
}

get_data <- function(dataset,seed=123123) {
  set.seed(seed)
  df <- switch(dataset,
               "Xs"=double_X(),
               "three gaussians"=three_gaussians(),
               "two lines"=two_lines(),
               "circles"=circles(),
               "trefoil"=trefoil(),
               "clusters"=high_dim_clusters())
  df
}

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
  if(dim(high)[2]==dim(low)[2]){
    bayes <- my_procrust_dist(high,low,300)
  }
  else{
    bayes <- get_knn_res(low,num_groups = 16,group_sizes = 50)
  }

  data_row <- matrix(c(kmax_val,qnx_val,h,mi,le,rk,spear,bayes),nrow=1)
  return(data.frame(data_row))
}

get_grid <- function(method = "tsne", dataset = "two lines") {
  if(method == "tsne") {
    if(dataset == "Xs") {
      grid = expand.grid(
        perplexity = seq(5,120,10),
        pca_center = c(TRUE,FALSE),
        pca_scale = c(TRUE,FALSE)
      )
      return(grid)

    } else if(dataset == "three gaussians"){

      grid = expand.grid(
        perplexity = seq(5,45,5),
        pca_center = c(TRUE,FALSE),
        pca_scale = c(TRUE,FALSE)
      )
      return(grid)
    } else if(dataset == "two lines"){
      grid = expand.grid(
        perplexity = seq(5,45,5),
        pca_center = c(TRUE,FALSE),
        pca_scale = c(TRUE,FALSE)
      )
      return(grid)
    } else if(dataset == "circles"){
      grid = expand.grid(
        perplexity = seq(5,80,5),
        pca_center = c(TRUE,FALSE),
        pca_scale = c(TRUE,FALSE)
      )
      return(grid)

    } else if(dataset == "trefoil"){
      grid= expand.grid(
        perplexity = seq(5,80,5),
        pca_center = c(TRUE,FALSE),
        pca_scale = c(TRUE,FALSE)
      )
      return(grid)
    } else if(dataset == "clusters"){
      grid = expand.grid(
        perplexity = seq(5,250,20),
        pca_center = c(TRUE,FALSE),
        pca_scale = c(TRUE,FALSE)
      )
      return(grid)
    }else {
      stop("datasets are\
           (1) Xs\
           (2) three gaussians\
           (3) two lines\
           (4) trefoil\
           (5) circles\
           (6) clusters")
    }
  } else if(method == "umap") {
    if(dataset == "Xs") {
      grid <- expand.grid(
        k=c(5,10,20,50,100,250,300,350),
        bandwidth=seq(1,11,2),
        local_connectivity=seq(1,11,2)
      )
      return(grid)
    } else if(dataset == "three gaussians"){
      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        bandwidth=seq(1,11,2),
        local_connectivity=seq(1,11,2)
      )
      return(grid)
    } else if(dataset == "two lines"){

      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        bandwidth=seq(1,11,2),
        local_connectivity=seq(1,11,2)
      )
      return(grid)
    } else if(dataset == "circles"){
      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        bandwidth=seq(1,11,2),
        local_connectivity=seq(1,11,2)
      )
      return(grid)
    } else if(dataset == "trefoil"){
      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        bandwidth=seq(1,11,2),
        local_connectivity=seq(1,11,2)
      )
      return(grid)
    } else if(dataset == "clusters"){
      grid <- expand.grid(
        k=c(5,10,40,80,120,160,240),
        bandwidth=seq(1,11,2),
        local_connectivity=seq(1,11,2)
      )
      return(grid)
    } else {
      stop("datasets are\
           (1) Xs\
           (2) three gaussians\
           (3) two lines\
           (4) trefoil\
           (5) circles\
           (6) clusters")
    }
  } else if(method == "lmds") {
    if(dataset == "Xs") {
      grid <- expand.grid(
        k=c(5,10,20,50,100,250,300,350),
        tau=10^seq(-3,0,0.5)
      )
      return(grid)
    } else if(dataset == "three gaussians"){
      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        tau=10^seq(-3,0,0.5)
      )
      return(grid)
    } else if(dataset == "two lines"){
      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        tau=10^seq(-3,0,0.5)
      )
      return(grid)
    } else if(dataset == "circles"){
      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        tau=10^seq(-3,0,0.5)
      )
      return(grid)
    } else if(dataset == "trefoil"){
      grid <- expand.grid(
        k=c(5,10,20,30,40,50,75,90,100),
        tau=10^seq(-3,0,0.5)
      )
      return(grid)
    } else if(dataset == "clusters"){
      grid <- expand.grid(
        k=c(5,10,40,80,120,160,240),
        tau=10^seq(-3,0,0.5)
      )
      return(grid)
    } else {

      stop("datasets are\
           (1) Xs\
           (2) three gaussians\
           (3) two lines\
           (4) trefoil\
           (5) circles\
           (6) clusters")
    }
  } else {
    stop("Choose: (1) tsne (2) umap (3) lmds")
  }
}
