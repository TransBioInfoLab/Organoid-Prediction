# Prepare results for testing data
prepare_results <- function(model = "RF", mod, test_list, ...){
  
  if(model == "RF"){
    
    exp <- test_list$exp
    clinical <- test_list$clinical
    test <- exp[match(colnames(mod$xvar), rownames(exp)), ]
    rownames(test) <- colnames(mod$xvar)
    p <- predict(mod, 
                 data.frame(scale(t(test))), 
                 na.action = "na.impute")
    re <- get_test_results(p$predicted, 
                           time_var = "os_month", 
                           event_var = "os", 
                           pheno_mat = clinical, plot = T, ...)
    
    re <- re %>%
      mutate(num_gene_selected = ncol(mod$xvar))
  }
  
  if(model == "ridge"){
    # Get scores
    scores <- get_scores(test = test_list$exp,
                         coef = mod$coeff)
    
    # Get test
    
    re <- get_test_results(test = scores,
                           time_var = "os_month", 
                           event_var = "os", 
                           pheno_mat = test_list$clinical, 
                           plot = T,
                           ...)
    
    re <- re %>%
      mutate(num_gene_selected = length(mod$coeff[mod$coeff != 0]))
  }
  
  if(model == "LM"){
    # Get scores
    scores <- get_scores(test = test_list$exp,
                         coef = mod)
    
    # Get test
    
    re <- get_test_results(test = scores,
                           time_var = "os_month", 
                           event_var = "os", 
                           pheno_mat = test_list$clinical, 
                           plot = T,
                           ...)
    
    re <- re %>%
      mutate(num_gene_selected = length(mod))
  }
  
  
  return(re)
}


# Prepare results for testing data
prepare_KMplot <- function(model = "RF", mod, test_list, save = T, dir_save = ".", prefix2 = NULL, ...){
  
  
  exp <- test_list$exp
  clinical <- test_list$clinical
  if(model == "RF"){

    test <- exp[match(colnames(mod$xvar), rownames(exp)), ]
    rownames(test) <- colnames(mod$xvar)
    p <- predict(mod, 
                 data.frame(scale(t(test))), 
                 na.action = "na.impute")
    pl <- KM_plot(p$predicted, 
                   time_var = "os_month", 
                   event_var = "os", 
                   pheno_mat = clinical, cut = "maxstat", ...)
    
  }
  
  if(model == "ridge"){
    # Get scores
    scores <- get_scores(test = test_list$exp,
                         coef = mod$coeff)
    
    pl <- KM_plot(scores, 
            time_var = "os_month", 
            event_var = "os", 
            pheno_mat = clinical, 
            cut = "maxstat",  ...)

  }
  
  if(model == "LM"){
    # Get scores
    scores <- get_scores(test = test_list$exp,
                         coef = mod$coeff)
    
    # Get test
    
    pl <- KM_plot(scores, 
                  time_var = "os_month", 
                  event_var = "os", 
                  pheno_mat = clinical, cut = "maxstat", ...)
    
  }
  
  if(save){
    pdf(file.path(dir_save, paste0(prefix2, "_", model, "_KMplot.pdf")),
        width = 8,
        height = 6)
    print(pl, newpage = F)
    dev.off()
    
  }
  
  return(pl)
}


# Get test results
get_test_results <- function(test_var, time_var, event_var, pheno_mat, cut = "maxstat", 
                             plot = F, Y5 = F, ...){
  
  df <- data.frame(test_var = test_var, time = pheno_mat[[time_var]], death = pheno_mat[[event_var]])
  if(Y5){
    df <- df %>% mutate(death = ifelse(death == 1 & time <= 60, 1, 0),
                        time = ifelse(time > 60, 60, time) )
  }
  
  if(is.numeric(test_var)){
    
    if(cut == "median"){
      m <- median(test_var)
    }
    
    if(cut == "mean"){
      m <- mean(test_var)
    }
    
    if(cut == "maxstat"){
      
      m <- maxstat::maxstat.test(Surv(time, death) ~ test_var, data = df, smethod = "LogRank",
                                 minprop=0.15, maxprop=0.85)$estimate
      
    }
    gene_cut <- ifelse(test_var < m, "low", "high")
    
  } else {
    
    gene_cut = test_var
    
  } 
  
  df <- df %>% mutate(scores = gene_cut)
  # Create stratification
  df <- df %>% mutate(
    follow_up_Y1 = ifelse(time <= 12 & death == 1, 1, 0),
    follow_up_Y3 = ifelse(time <= 36 & death == 1, 1, 0),
    follow_up_Y5 = ifelse(time <= 60 & death == 1, 1, 0)
  )
  
  fo <- as.formula(paste0("Surv(time, death) ~ scores"))
  
  if(plot){
    KM_plot(test_var, "time", "death", df, cut = cut) %>% print()
  }
  log_rank <- survdiff(fo, data = df)
  
  p <- data.frame(log_rank_p = log_rank$p, event_low = log_rank$obs[2]/log_rank$n[2], 
                  event_high = log_rank$obs[1]/log_rank$n[1],
                  direction = ifelse(log_rank$obs[2]/log_rank$n[2] > log_rank$obs[1]/log_rank$n[1], ">", "<"))
  
  auc_ls <- plyr::llply(
    c("follow_up_Y1", "follow_up_Y3", "follow_up_Y5", "ALL"),
    .fun = function(y){
      if(y != "ALL"){
        if(length(na.omit(unique(df[[y]]))) == 1){
          auc <- NULL
        } else {
          auc <- data.frame(follow_up = paste0("AUC_", y),
                            auc = pROC::roc(factor(df[[y]], levels = c(0,1)), 
                                            df$test_var)$auc)
        }
        
      } else {
        auc <- data.frame(follow_up = paste0("AUC_", y),
                          auc = pROC::roc(factor(df[["death"]], levels = c(0,1)), 
                                          df$test_var)$auc)
      }
      
      return(auc)
      
    }
  ) %>% purrr::compact() %>%
    purrr::reduce(., rbind) %>%
    column_to_rownames("follow_up") %>%
    t() %>% 
    data.frame()
  
  cbind(auc_ls, p)
  
}

# Get scores
get_scores <- function(test, coef, sample_sel = "ALL"){
  
  if(sample_sel == "ALL"){
    test_selected <- test[rownames(test) %in% names(coef),]
  } else {
    test_selected <- test[rownames(test) %in% names(coef), sample_sel]
  }
  
  if(is.null(nrow(test_selected))){
    scores <- NULL
  } else if (nrow(test_selected) == 0){
    scores <- NULL
  } else{
    test_selected <- test_selected[rowSums(test_selected)!=0,]
    
    test_selected <- t(scale(t(test_selected)))
    coef_sel <- coef[rownames(test_selected)]
    
    scores <- (coef_sel %*% test_selected) %>%
      as.numeric()
    
    names(scores) <- colnames(test_selected)
  }
  return(scores)
}
