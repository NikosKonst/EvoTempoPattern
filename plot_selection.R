## This R file contains functions used for plot selection ###
## get_quantile_matrix: uses the pseudotime assignment to bin cells 
## inputs: 
# pseudotime: a vector of the pseudotime values calculated for each cell 
# n_permutation: number of permutations to perform, used to generate a null distribution of pseudotime assignments 
# intervals: number of intervals to bin 

## output: list containing the binned assignment for the original data, as well as the permuted sets
get_quantile_matrix <- function(pseudotime, n_permutation = 100, intervals = 20){
  quantiles <-  seq(from = 0, to =1, length = intervals + 1)   
  ##add a small bit, otherwise using < for the upper quantile will always omit the cells with pseudotime of 1
  quantiles[length(quantiles)] <- 1.01
  pseudotime <-  pseudotime/(max(pseudotime)-min(pseudotime))
  
  quantile_assignment <-  matrix(0, ncol = intervals, nrow = length(pseudotime))
  for(q in 1:intervals){
    lower <- quantiles[q]
    upper <- quantiles[(q+1)]
    count <- which(pseudotime >= lower & pseudotime<upper)
    quantile_assignment[count, q] <- 1
  }
  
  permuted_intervals <- list()
  for(i in 1:n_permutation){
    permuted_pseudotime <- pseudotime
    set.seed(i)
    permuted_pseudotime <- sample(permuted_pseudotime)
    permuted_quantiles <- matrix(0, ncol = intervals, nrow = length(pseudotime))
    for(q in 1:intervals){
      
      lower <- quantiles[q]
      upper <- quantiles[(q+1)]
      count <- which(permuted_pseudotime >= lower & permuted_pseudotime<upper)
      permuted_quantiles[count, q] <- 1
    }
    permuted_intervals[[i]]<- permuted_quantiles
  
  return(list(original = quantile_assignment, permuted = permuted_intervals))
}}

#####Reach0Filter: determines whether or not a gene's expression tends to 0#####
##inputs: 
#sample_genes: a vector of genes to test
#quantile assignment: the interval assignment for each cell, ie the output of get_quantile_matrix, entry 'original'
#count_matrix: this is a count matrix. Should be genes x cells, with genes as rownames. Cells should be in the same order as the quantile_assignment matrix
#logbase: the base to use for log normalization, default is 10 (but can be adapted)
#approach0_max: if the mean over any interval < approach0_max, this qualifies as reaching 0. Default is 0.1. 
#min_cells_for_expression: requires at least this proportion of cells to be nonzero, otherwise the mean is 0. Default is 0.05
#numCores: number of cores to use for parallelization, default is 6 (can be set based on the user)

###output: a vector containing 'pass' or 'fail' for each gene. 
reach0filter = function(sample_genes, quantile_assignment, count_matrix, logbase = 10, approach0_max = 0.1,min_cells_for_expression = 0.05, numCores = 6){
  registerDoParallel(cores = numCores)
  pvalues <- foreach(i = 1:length(sample_genes), .combine = 'c') %dopar% {
    
    gene = sample_genes[i]
    median = c()
    for(q in 1:ncol(quantile_assignment)){
      cells = which(quantile_assignment[, q]!=0)
      exp = log(count_matrix[gene, cells]+1, logbase)
      cutoff = floor(length(exp)*min_cells_for_expression)
      #exp = exp[exp!=0]
      median[q] = ifelse(length(exp)>cutoff, mean(exp), 0)
    }
    ifelse(sum(median<approach0_max)>0, 'pass', 'fail')
  }
  
  stopImplicitCluster()
  return(pvalues)
}

######Gene exp: examines gene expression and requires a minimum mean expression and minimum proportion of non-zero genes ##########
##inputs: 
#sample_genes: a vector of genes to test
#count_matrix: this is a count matrix. Should be genes x cells, with genes as rownames. Cells should be in the same order as the quantile_assignment matrix
#min_cells_for_expression: requires at least this proportion of cells to be nonzero, otherwise the median is 0. Default is 0.05
#numCores: number of cores to use for parallelization, default is 6 (because in the file you sent me, you used 6)
#min_exp: this is the minimum mean expression needed for the gene to pass the filter. Can be adapted, tried between the values:  log (2.25, base=10)= 0.35, log (2.45, base=10)= 2.45, log (2.8, base=10)=0.447

###output: a vector containing 'pass' or 'fail' for each gene. 

gene_exp_filter <- function(sample_genes,count_matrix, logbase = 10, numCores = 6,min_cells_for_expression = 0.05, min_exp = log(2.45, base = 10) ){
  registerDoParallel(cores = numCores)
  pvalues <- foreach(i = 1:length(sample_genes), .combine = 'c') %dopar% {
    gene <- sample_genes[i]
    counts <- count_matrix[gene, ]
    exp <- log(counts+1, base = logbase)
    #median_0 <- median(exp)
    mean_0 <- mean(exp)
    cell_count <- length(counts)
    exp <- exp[exp!=0]
    
    ## for nonzero gene expression, we adjusted this value to require 5% of cells to be nonzero for the gene, otherwise the median expression was 0 
    if(length(exp)!=0){
      
      mean_non0 <- ifelse(length(exp)> cell_count*min_cells_for_expression, mean(exp), 0)
    } else{
      sd <- 0
      mean_non0 <- 0
    } 
    
    
    ifelse(mean_non0>= min_exp, 'pass', 'fail')
  }
  
  stopImplicitCluster()
  return(pvalues)
}



