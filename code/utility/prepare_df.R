
# Prepare training dataset list
prepare_train_df <- function(dir.base = "../../", quantile_keep = 0.5){
  
  dir.data <- file.path(dir.base, "data/")
  # GSE171682 (exp_filtered_list)
  load(file = file.path(dir.data, 
                        "GSE171682/processed/GSE171682_colorectal_cancer_exp_data_common_genes.RDS"))
  
  # Load survival outcome for GSE171682
  GSE171682c <- read.csv(
    file.path(dir.data, "GSE171682/clinical/GSE171682_colorectal_cancer_clinical_data.csv")
  )
  
  # GSE64392 (exp_CC, pheno_CC, drug_res_by_patient)
  load(file = file.path(dir.data,
                        "GSE64392/processed/COAD_GSE64392_org_CC_with_drug_info_matched_exp_clinical_data.rda"))
  
  # Load clinical and drug response for GSE64392
  GSE64392c <- pheno_CC
  GSE64392_drug <- drug_res_by_patient
  GSE64392c <- left_join(drug_res_by_patient, GSE64392c)
  
  # Prepare drug response
  GSE64392c_5fu <- GSE64392c %>% 
    filter(`Drug Name` == "5-Fluorouracil")
  
  # Prepare survival data
  GSE171682c <- GSE171682c %>% 
    mutate(sex = factor(sex),
           tim_class = factor(tim_class),
           os_month = os_day/30,
           rfs_month = rfs_day/30)
  
  # Match common genes
  common_genes <- intersect(
    rownames(exp_filtered_list$GSE171680),
    rownames(exp_CC)
  )
  
  length(common_genes) # 15836
  
  exp_filtered <- exp_CC[common_genes,GSE64392c_5fu$geo_accession]
  exp_filtered_list <- exp_filtered_list %>%
    map(~.[common_genes,])
  
  # Add GSE39582 to list
  exp_filtered_list <- c(exp_filtered_list, list("GSE64392" = exp_filtered))
  
  # Get top variable genes
  exp_top_genes <- exp_filtered_list %>%
    purrr::map(~get_most_variable_gene(., quantile_keep))
  
  # Select common genes in top variable genes
  common_genes <- intersect(
    intersect(
      rownames(exp_top_genes$GSE171680), rownames(exp_top_genes$GSE171681)
    ), rownames(exp_top_genes$GSE64392)
  )
  
  exp_top_genes <- exp_top_genes %>%
    purrr::map(~.[common_genes,])
  
  return(list(
    exp = exp_top_genes$GSE64392,
    mIC50 = GSE64392c_5fu$mIC50,
    survival = GSE171682c,
    exps_ls = exp_top_genes
  ))
  
}

# Prepare single training dataset
prepare_single_train_df <- function(GEO, dir.base = "../../", quantile_keep = 0.5){
  
  dir.data <- file.path(dir.base, "data/")
  # GSE171682 (exp_filtered_list)
  if(GEO == "GSE64392"){
    
    # GSE64392 (exp_CC, pheno_CC, drug_res_by_patient)
    load(file = file.path(dir.data,
                          "GSE64392/processed/COAD_GSE64392_org_CC_with_drug_info_matched_exp_clinical_data.rda"))
    
    # Load clinical and drug response for GSE64392
    GSE64392c <- pheno_CC
    GSE64392_drug <- drug_res_by_patient
    GSE64392c <- left_join(drug_res_by_patient, GSE64392c)
    
    # Prepare drug response
    clinical <- GSE64392c %>% 
      filter(`Drug Name` == "5-Fluorouracil")
    
    exp <- exp_CC[,clinical$geo_accession]
    
    exp <- get_most_variable_gene(exp, quantile_keep)
    
  }
  
  if(GEO == "GSE171682"){
    
    load(file = file.path(dir.data, 
                          "GSE171682/processed/GSE171682_colorectal_cancer_exp_data_common_genes.RDS"))
    
    # Load survival outcome for GSE171682
    GSE171682c <- read.csv(
      file.path(dir.data, "GSE171682/clinical/GSE171682_colorectal_cancer_clinical_data.csv")
    )
    
    # Prepare survival data
    clinical <- GSE171682c %>% 
      mutate(sex = factor(sex),
             tim_class = factor(tim_class),
             os_month = os_day/30,
             rfs_month = rfs_day/30)
    
    exp <- exp_filtered_list %>%
      purrr::map(~get_most_variable_gene(., quantile_keep))
    
  }
  
  return(list(
    exp = exp,
    clinical = clinical
  ))
  
}

# Prepare testing dataset
prepare_test_df <- function(GEO, dir.base = "../../", drug = T){
  
  GEO0 <- GEO
  dir.data <- file.path(dir.base, "data/", GEO0)
  dir.data.raw <- file.path(dir.data, "raw")
  dir.data.processed <- file.path(dir.data, "processed")
  dir.data.clinical <- file.path(dir.data, "clinical")
  
  if(GEO == "GSE39582"){
    
    load(file.path(dir.data.processed, 
                   "GSE39582_CC_with_drug_info_matched_exp_clinical_data.rda"))
    exp <- exp2
    
    # Select drug type = 5FU
    clinical$os_month <- as.numeric(clinical$os_delay_months_ch1)
    clinical$os <- as.numeric(clinical$os_event_ch1)
    clinical$rfs_month <- as.numeric(clinical$rfs_delay_ch1)
    clinical$rfs <- as.numeric(clinical$rfs_event_ch1)
    if(drug){
      clinical <- clinical %>% 
        filter(chemotherapy_adjuvant_type_ch1 %in% c("5FU", "FOLFOX"))
    } 
    exp <- exp[,clinical$geo_accession]
    
  }
  
  if(GEO == "GSE17538"){
    load(file.path(dir.data.processed, "GSE17538_with_rfs_exp_clinical_data.rda"))
    exp <- exp2
    
    if(drug){
      clinical <- clinical %>% filter(AdjCTX == "Y")
    } 
    
    exp <- exp[,clinical$geo_accession]
    
  }
  
  if(GEO == "GSE209746"){
    load(file.path(dir.data.processed, "GSE209746_exp_cBioportal_clinical_data.rda"))
    
    if(drug){
      clinical <- clinical %>% 
        filter(Neoadj_Chemoagent == "FOLFOX") 
    }
    clinical <- clinical %>%
      mutate(os = OS_Status,
             os_month = OS_months_from_neo)
    exp <- exp[,clinical$title]
  }
  
  if(GEO == "COAD"){
    load(file.path(dir.data.processed, "TCGA_COAD_CC_matched_exp_clinical_data.rda"))
    
    clinical$os_month <- clinical$os_day/30
    clinical$OS.time <- clinical$OS.time/30
    clinical$PFI.time <- clinical$PFI.time/30
    clinical.drug <- read.csv(
      file.path(dir.data.clinical, "TCGA_COAD_colorectal_cancer_drug_type_info.csv")
    )
    if(drug){
      clinical.drug_5FU <- clinical.drug %>% 
        filter(drug_name %in% c("5 FU", "5FU", "5-FU", "Fluorouracil", 
                                "5- FU", "FOLFOX","Folfox","FolFox",
                                "5-Fluorouracil", "5-Fluoruoracil", 
                                "FLUOROURACIL", "fluorouracil"))
      clinical <- clinical %>% filter(bcr_patient_barcode %in% clinical.drug_5FU$bcr_patient_barcode)
    }
    
    exp <-   exp[,clinical$bcr_patient_barcode]
  }
  

  
  if(GEO == "GSE14333"){
    load(file.path(dir.data.processed, "GSE14333_matched_exp_clinical_data.rda"))
    exp <- exp2
    clinical$os_month <- as.numeric(clinical$DFS_Time)
    clinical$os <- as.numeric(clinical$DFS_Cens)
    
    if(drug){
      clinical <- clinical %>% filter(AdjCTX == "Y")
    }
   
    exp <- exp[,clinical$geo_accession]
  }
  
  if(GEO == "GSE106584"){
    load(file.path(dir.data.processed, "GSE106584_with_rfs_exp_clinical_data.rda"))
    if(drug){
      clinical <- clinical %>% 
        filter(!adjuvan_therapy_ch1 %in% c("none", "None", "No", "Unknown", "N/A", "None (chemoXRT for concurrent ACC)"))
    }
    
    exp <- exp2[,clinical$geo_accession]
  }
  
  if(GEO == "GSE72970"){
    load(file.path(dir.data.processed, "GSE72970_with_rfs_exp_clinical_data.rda"))
    exp <- exp2[,clinical$geo_accession]
  }
  
  if(GEO == "GSE87211"){
    load(file.path(dir.data.processed, "GSE87211_with_rfs_exp_clinical_data.rda"))
    exp <- exp2[,clinical$geo_accession]
  }
  
  return(list(
    exp = exp,
    clinical = clinical
  ))
}
