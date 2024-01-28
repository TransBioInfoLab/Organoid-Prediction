# Function for calculate drug response
corr_coef <- function(exp, pheno, drug_var, parallel = T, ...){
  
  if(parallel) doParallel::registerDoParallel(10)
  
  plyr::adply(
    exp,
    .margins = 1,
    .fun = function(ge){
      
      get_corr_coef(
        gene = ge, 
        pheno = pheno,
        drug_var = drug_var,
        ...
      )
      
    }, .id = "gene", .parallel = parallel
  )
}

get_corr_coef <- function(gene, pheno, drug_var, method = "pearson", ...){
  
  # fit cox regression
  corr_mod <- cor.test(
    gene, pheno[[drug_var]],
    method = method,
    ...
  )
  
  # get summary and coefficient
  corr_coef <-  broom::tidy(corr_mod)
  
  return(corr_coef)
}


# Function for gene-gene correlation
get_corr_coef_gene <- function(gene1, gene2, method = "spearman", ...){
  
  coef_df <- cor.test(gene1, gene2, method = method, ...) %>% 
    broom::tidy() %>%
    clean_names()
  
  return(coef_df)
}



corr_coef_gene <- function(exp1, exp2, parallel = T, ...){
  
  if(parallel) doParallel::registerDoParallel(10)
  
  g <- rownames(exp1)
  exp2 <- exp2[g,]
  
  re <- plyr::ldply(
    1:nrow(exp1),
    .fun = function(ge){
      
      suppressWarnings({
        get_corr_coef(
          gene1 = exp1[ge,], 
          gene2 = exp2[ge,],
          ...
        )
      })
      
    }, .parallel = parallel
  )
  
  re <- re %>%
    mutate(gene = g, .before = 1)
  
}