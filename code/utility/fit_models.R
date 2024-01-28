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

## RF variable selection

rf_perm_test <- function(formula, data, nperm = 50, thres = 0.05, refit = T, seed = 123, ...){
  
  vimp_org <- vimp_perm <- NULL
  
  response <- all.vars(as.formula(formula))[1]
  
  for(i in c(1:nperm)){
    
    mrf_mod <- rfsrc(
      as.formula(formula), 
      data = data,
      importance = "permute",
      ...
    )
    
    vimp_org <- rbind(vimp_org, mrf_mod$importance)
    
    dat_perm <- data[sample.int(nrow(data), nrow(data)),!colnames(data) %in% response]
    dat_perm <- cbind(data[[response]], dat_perm + abs(dat_perm))
    colnames(dat_perm)[1] <- response
    
    mrf_obj <- rfsrc(
      as.formula(formula), 
      data = dat_perm,
      importance = "permute",
      ...
    )
    
    vimp_perm <- rbind(vimp_perm, mrf_obj$importance)
    
  }
  
  test_df <- plyr::ldply(
    1:(ncol(data) - 1),
    .fun = function(p){
      
      w.p <- wilcox.test(vimp_org[,p], vimp_perm[,p], exact = F, paired = F, alternative = "greater")$p.value
      
      data.frame(pvalue = w.p, mean_vimp = mean(vimp_org[,p]))
    }
  )
  
  var_name <- colnames(data)[!colnames(data) %in% response]
  test_df$var <- var_name
  final_sel <- test_df %>% 
    filter(pvalue < thres & mean_vimp > 0) %>% 
    pull(var)
  
  out <- list(test = test_df, sig_genes = final_sel)
  
  if(refit){
    
    set.seed(seed)
    mrf_mod <- rfsrc(
      as.formula(formula), 
      data = data[,c(response, final_sel)],
      importance = "permute",
      ...
    )
    
    out <- c(out, list(refit = mrf_mod))
  }
  
  
  
  return(out)
  
}
