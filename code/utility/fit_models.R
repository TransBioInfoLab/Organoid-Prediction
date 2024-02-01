# Fit ridge, RF, lm model

fit_ridge <- function(gene_select, exp, y, k = 5, a = NULL,
                      lambda = exp(seq(-2, 4, length.out = 50)), type = "regression",...){
  
  x <- t(exp[gene_select,])
  x <- scale(as.matrix(x))
  y <- y
  
  if(is.null(a)){
    a <- seq(0, 0.9, 0.1)
    # Five-fold cross-validation
    doParallel::registerDoParallel(4)
    search <- foreach(i = a, .combine = rbind) %dopar% {
      
      cv <- cv.glmnet(x = x, y = y, lambda = lambda,
                      nfolds = k, alpha = i, ...)
      data.frame(cvm = cv$cvm, lambda.1se = cv$lambda.min, alpha = i)
    }
    
    
    if(type == "survival"){
      cv <- search[which.max(search$cvm), ]
      cv_min <- max(cv$cvm)
    } else {
      cv <- search[which.min(search$cvm), ]
      cv_min <- min(cv$cvm)
    }
    l <- cv$lambda.1se
    alpha <- cv$alpha
  } else {
    
    cv <- cv.glmnet(x = x, y = y, lambda = lambda,
                    nfolds = k, alpha = a, ...)
    l <- cv$lambda.min
    alpha <- a
    if(type == "survival"){
      cv_min <- max(cv$cvm)
    } else {
      cv_min <- min(cv$cvm)
    }
  }
  
  # Fit ridge
  ridge <- glmnet(x = x, 
                  y = y, 
                  lambda = l, 
                  alpha = alpha,
                  ...)
  
  if(type == "survival"){
    coeff <- as.numeric(coef(ridge))
  } else {
    coeff <- coef(ridge)[-1]
  }
  
  if(sum(coeff) == 0){
    cv <- cv.glmnet(x = x, y = y, lambda = lambda,
                    nfolds = k, alpha = 0, ...)
    
    if(type == "survival"){
      cv_min <- max(cv$cvm)
    } else {
      cv_min <- min(cv$cvm)
    }
    l <- cv$lambda.min
    ridge <- glmnet(x = x, 
                    y = y, 
                    lambda = l, 
                    alpha = 0,
                    ...)
    
    if(type == "survival"){
      coeff <- as.numeric(coef(ridge))
    } else {
      coeff <- coef(ridge)[-1]
    }
  }
  
  names(coeff) <- colnames(x)
  
  mod <- list(
    coeff = coeff,
    cv_min = cv_min,
    lambda_min = l,
    a = alpha,
    mod = ridge
  )
  
  return(mod)
  
}

fit_rf <- function(gene_select, exp, y, mode = "VIMP", seed = 942, type = "regression",  ...){

  x <- t(exp[gene_select,])
  if(type == "regression"){
    data = data.frame(mIC50 = y,
                      scale(x))
    fo <- as.formula("mIC50 ~ .")
  } else {
    data = data.frame(y, scale(x))
    
    fo <- as.formula("Surv(time = os_month, event = death) ~ .")
  }
  
  if(mode == "perm"){
    
    rf <- rf_perm_test(fo, 
                       data = data,
                       nperm = 50,
                       thres = .01,
                       seed = seed,...)
    rf_sel <- rf$refit
  }
  
  if(mode == "VIMP"){
    
    set.seed(seed)
    rf <- rfsrc(fo, 
                data = data,
                importance = "permute",...)
    selected_genes <- names(rf$importance)[rf$importance > 0]
    
    if(type == "regression"){
      data = data.frame(mIC50 = y,
                        data.frame(scale(x))[,selected_genes])
      
    } else {
      data = data.frame(y, data.frame(scale(x))[,selected_genes])
    }
    set.seed(seed)
    rf_sel <- rfsrc(fo, 
                    data = data,
                    importance = "permute")
  }
  
  if(mode == "None"){
    set.seed(seed)
    rf_sel <- rfsrc(fo, 
                    data = data,
                    importance = "permute",...)
  }
  
  return(rf_sel)
  
}

fit_lm <- function(gene_select, exp, y, ...){
  
  x <- t(exp[gene_select,])
  
  data = data.frame(mIC50 = y,
                    scale(x))
  fo <- as.formula("mIC50 ~ .")
  
  lm_mod <- lm(fo, data = data, ...)
  coeff <- coefficients(lm_mod)[-1]
  names(coeff) <- colnames(x)
  coeff <- na.omit(coeff)
  
  return(list(coeff = coeff, mod = lm_mod))
  
}
