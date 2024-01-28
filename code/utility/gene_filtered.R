# Function for removing 0s and low variation genes
gene_filtered <- function(exp, log2_trans = F, zero_thres = 1){
  
  if(log2_trans){
    exp <- log2(exp + 1)
  }
  
  # calculate 0 percentage
  gene_zeros <- rowSums(exp == 0)/ncol(exp)
  
  # filter out all 0 genes
  gene_filtered <- exp[gene_zeros < zero_thres,]
  
  return(gene_filtered)
  
}

# Get most variable genes
get_most_variable_gene <- function(exp, quantile = 0.5){
  
  var <- rowVars(exp)
  exp <- exp[var >= quantile(var,quantile),]
  
  return(exp)
}