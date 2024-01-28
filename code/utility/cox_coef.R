# Function for cox regression

# Function for check warnings
check_warning <- function(formula, data,...){
  
  mod <- tryCatch({
    coxph(
      formula,
      data = data,
      ...
    )
  }, warning = function(w) return(NULL))
  
  if(!is.null(mod) && class(mod) == "coxph"){
    w <- F
  } else {w <- T}
  
  return(w)
}


# Function for cox regression analysis
get_cox_coef <- function(gene, pheno, time_var, event_var, adjust_var = NULL, ...){
  
  # make survival formula
  fo <-  paste0("Surv(", time_var, ",", event_var, ") ~ gene")
  if(!is.null(adjust_var)){
    fo <- paste0(fo, " + ", paste0(adjust_var, collapse = "+"))
  }
  formula <- as.formula(fo)
  
  data <- data.frame(gene = gene, pheno)
  
  # check warning
  w <- check_warning(
    formula, 
    data,
    ...
  )
  
  # fit cox regression
  cox_mod <- coxph(
    formula,
    data = data,
    ...
  )
  
  # get summary and coefficient
  cox_coef <- summary(cox_mod)$coefficients %>% data.frame()
  
  coef_df <- cox_coef[grepl("gene", rownames(cox_coef)),] %>% 
    clean_names()
  
  coef_df$warning <- w
  
  return(coef_df)
}



cox_coef <- function(exp, pheno, time_var, event_var, adjust_var = NULL, parallel = T, ...){
  
  if(parallel) doParallel::registerDoParallel(10)
  
  plyr::adply(
    exp,
    .margins = 1,
    .fun = function(ge){
      
      suppressWarnings({
        get_cox_coef(
          gene = ge, 
          pheno = pheno,
          time_var = time_var,
          event_var = event_var,
          adjust_var = adjust_var,
          ...
        )
      })
      
    }, .id = "gene", .parallel = parallel
  )
  
}

add_fdr <- function(results, method = "fdr"){
  
  results$fdr <- p.adjust(results$pr_z, method = method)
  
  return(results)
  
}
