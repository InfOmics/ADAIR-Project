#Paper: Air pollution and metabolomics

#Written under R version 4.0.4


#______________________________________________________________________________________
# Open dataset for analysis ----
load("Dataset_final.RData")

ergo_final$elapse_pm25 = scale(ergo_final$elapse_pm25)
ergo_final$elapse_pm25 <- as.numeric(ergo_final$elapse_pm25)
ergo_final$elapse_bc = scale(ergo_final$elapse_bc)
ergo_final$elapse_bc <- as.numeric(ergo_final$elapse_bc)
ergo_final$elapse_no2 = scale(ergo_final$elapse_no2)
ergo_final$elapse_no2 <- as.numeric(ergo_final$elapse_no2)
ergo_final$elapse_o3a = scale(ergo_final$elapse_o3a)
ergo_final$elapse_o3a <- as.numeric(ergo_final$elapse_o3a)
ergo_final$elapse_o3w = scale(ergo_final$elapse_o3w)
ergo_final$elapse_o3w <- as.numeric(ergo_final$elapse_o3w)

# packages
library(dplyr)
library(tidyr)
library(xlsx)

list_mb <- read.xlsx("annotation_fileMetabolon.xlsx", sheetIndex = "Sheet1", header = TRUE)


#-------------------------------------------------------------------#
#                    data prep metabolomics                         #
#-------------------------------------------------------------------#

#  create vector of metabolites for analysis
metabolites <- as.character(colnames(ergo_final[,21:960]))

# create an empty results table (this needs to be changed if you use a different model function)

#-------------------------------------------------------------------#
#                        PM  -   Model 2                            #
#-------------------------------------------------------------------#

  # Loop runs through each metabolite seperately 
for(a in 1:length(ergo_final)) {results <- data.frame(Metabolite=as.character(), Beta=as.numeric(),Se=as.numeric(),Low95=as.numeric(), Up95=as.numeric(),
                                                p=as.numeric(), n=as.numeric(),stringsAsFactors=FALSE)
for (i in 1:length(metabolites)){
    ergo_final$scalemetabo <- scale(ergo_final[,metabolites[i]])
    # define your model covariates
    test<- lm(scalemetabo ~ elapse_pm25 + age + sex + rs_cohort + season + ses_UNESCO_recoded + BMI + smk + METh + Occupation + alcohol + CESD + lipidmede,
              data = ergo_final)
    # fill out the results table (this needs to be changed if you use a different model function)
    tablerow <- data.frame( Metabolite=metabolites[i], Beta=summary(test)$coefficients[2,1],Se=summary(test)$coefficients[2,2], Low95=confint(test)[2,1],
      Up95=confint(test)[2,2], p=summary(test)$coefficients[2,4], n=nobs(test), stringsAsFactors=FALSE)
    results <- rbind(results, tablerow)
  } 
  # merge full metabolite name into results table
  results$CHEMICAL_NAME <- list_mb$CHEMICAL_NAME[match(results$Metabolite, list_mb$CHEM_ID)]
  results_PM <- results[order(results$p),]
  #Generate q-value, i.e. correct for multiple testing using the false discovery rate of Benjamini-Hochberg
  results_PM$q_values <- p.adjust(results_PM$p, method = "fdr")
  save(results_PM, file="V:\\HomeDir\\051543(T_de_Crom)\\Syntax\\Paper air pollution and metabolomics\\Results\\PM_Metabolon_M1.RData")
}



#-------------------------------------------------------------------#
#                        BC  -   Model 2                            #
#-------------------------------------------------------------------#

# Loop runs through each metabolite seperately 
for(a in 1:length(ergo_final)) {results <- data.frame(Metabolite=as.character(), Beta=as.numeric(),Se=as.numeric(),Low95=as.numeric(), Up95=as.numeric(),
                                                p=as.numeric(), n=as.numeric(),stringsAsFactors=FALSE)
for (i in 1:length(metabolites)){
  ergo_final$scalemetabo <- scale(ergo_final[,metabolites[i]])
  # define your model covariates
  test<- lm(scalemetabo ~ elapse_bc + age + sex + rs_cohort + season + ses_UNESCO_recoded + BMI + smk + METh + Occupation + alcohol + CESD + lipidmede,
            data = ergo_final)
  # fill out the results table (this needs to be changed if you use a different model function)
  tablerow <- data.frame( Metabolite=metabolites[i], Beta=summary(test)$coefficients[2,1],Se=summary(test)$coefficients[2,2], Low95=confint(test)[2,1],
                          Up95=confint(test)[2,2], p=summary(test)$coefficients[2,4], n=nobs(test), stringsAsFactors=FALSE)
  results <- rbind(results, tablerow)
} 
  # merge full metabolite name into results table
  results$CHEMICAL_NAME <- list_mb$CHEMICAL_NAME[match(results$Metabolite, list_mb$CHEM_ID)]
  results_BC <- results[order(results$p),]
  #Generate q-value, i.e. correct for multiple testing using the false discovery rate of Benjamini-Hochberg
  results_BC$q_values <- p.adjust(results_BC$p, method = "fdr")
  save(results_BC, file="BC_Metabolon_M1.RData")
}



#-------------------------------------------------------------------#
#                       NO2  -   Model 2                            #
#-------------------------------------------------------------------#

# Loop runs through each metabolite seperately 
for(a in 1:length(ergo_final)) {results <- data.frame(Metabolite=as.character(), Beta=as.numeric(),Se=as.numeric(),Low95=as.numeric(), Up95=as.numeric(),
                                                p=as.numeric(), n=as.numeric(),stringsAsFactors=FALSE)
for (i in 1:length(metabolites)){
  ergo_final$scalemetabo <- scale(ergo_final[,metabolites[i]])
  # define your model covariates
  test<- lm(scalemetabo ~ elapse_no2 + age + sex + rs_cohort + season + ses_UNESCO_recoded + BMI + smk + METh + Occupation + alcohol + CESD + lipidmede,
            data = ergo_final)
  # fill out the results table (this needs to be changed if you use a different model function)
  tablerow <- data.frame( Metabolite=metabolites[i], Beta=summary(test)$coefficients[2,1],Se=summary(test)$coefficients[2,2], Low95=confint(test)[2,1],
                          Up95=confint(test)[2,2], p=summary(test)$coefficients[2,4], n=nobs(test), stringsAsFactors=FALSE)
  results <- rbind(results, tablerow)
} 
  # merge full metabolite name into results table
  results$CHEMICAL_NAME <- list_mb$CHEMICAL_NAME[match(results$Metabolite, list_mb$CHEM_ID)]
  results_NO2 <- results[order(results$p),]
  #Generate q-value, i.e. correct for multiple testing using the false discovery rate of Benjamini-Hochberg
  results_NO2$q_values <- p.adjust(results_NO2$p, method = "fdr")
  save(results_NO2, file="NO2_Metabolon_M1.RData")
}


#-------------------------------------------------------------------#
#                       O3a  -   Model 2                            #
#-------------------------------------------------------------------#

# Loop runs through each metabolite seperately 
for(a in 1:length(ergo_final)) {results <- data.frame(Metabolite=as.character(), Beta=as.numeric(),Se=as.numeric(),Low95=as.numeric(), Up95=as.numeric(),
                                                p=as.numeric(), n=as.numeric(),stringsAsFactors=FALSE)
for (i in 1:length(metabolites)){
  ergo_final$scalemetabo <- scale(ergo_final[,metabolites[i]])
  # define your model covariates
  test<- lm(scalemetabo ~ elapse_o3a + age + sex + rs_cohort + season + ses_UNESCO_recoded + BMI + smk + METh + Occupation + alcohol + CESD + lipidmede,
            data = ergo_final)
  # fill out the results table (this needs to be changed if you use a different model function)
  tablerow <- data.frame( Metabolite=metabolites[i], Beta=summary(test)$coefficients[2,1],Se=summary(test)$coefficients[2,2], Low95=confint(test)[2,1],
                          Up95=confint(test)[2,2], p=summary(test)$coefficients[2,4], n=nobs(test), stringsAsFactors=FALSE)
  results <- rbind(results, tablerow)
} 
  # merge full metabolite name into results table
  results$CHEMICAL_NAME <- list_mb$CHEMICAL_NAME[match(results$Metabolite, list_mb$CHEM_ID)]
  results_O3a <- results[order(results$p),]
  #Generate q-value, i.e. correct for multiple testing using the false discovery rate of Benjamini-Hochberg
  results_O3a$q_values <- p.adjust(results_O3a$p, method = "fdr")
  save(results_O3a, file="O3a_Metabolon_M1.RData")
}


#-------------------------------------------------------------------#
#                       O3w  -   Model 2                            #
#-------------------------------------------------------------------#

# Loop runs through each metabolite seperately 
for(a in 1:length(ergo_final)) {results <- data.frame(Metabolite=as.character(), Beta=as.numeric(),Se=as.numeric(),Low95=as.numeric(), Up95=as.numeric(),
                                                p=as.numeric(), n=as.numeric(),stringsAsFactors=FALSE)
for (i in 1:length(metabolites)){
  ergo_final$scalemetabo <- scale(ergo_final[,metabolites[i]])
  # define your model covariates
  test<- lm(scalemetabo ~ elapse_o3w + age + sex + rs_cohort + season + ses_UNESCO_recoded + BMI + smk + METh + Occupation + alcohol + CESD + lipidmede,
            data = ergo_final)
  # fill out the results table (this needs to be changed if you use a different model function)
  tablerow <- data.frame( Metabolite=metabolites[i], Beta=summary(test)$coefficients[2,1],Se=summary(test)$coefficients[2,2], Low95=confint(test)[2,1],
                          Up95=confint(test)[2,2], p=summary(test)$coefficients[2,4], n=nobs(test), stringsAsFactors=FALSE)
  results <- rbind(results, tablerow)
} 
  # merge full metabolite name into results table
  results$CHEMICAL_NAME <- list_mb$CHEMICAL_NAME[match(results$Metabolite, list_mb$CHEM_ID)]
  results_O3w <- results[order(results$p),]
  #Generate q-value, i.e. correct for multiple testing using the false discovery rate of Benjamini-Hochberg
  results_O3w$q_values <- p.adjust(results_O3w$p, method = "fdr")
  save(results_O3w, file="O3w_Metabolon_M1.RData")
}
