#Paper: Air pollution and dementia risk

#Written under R version 4.0.3

#______________________________________________________________________________________
# Open dataset for analysis ----
load("Dataset_air_pollution_and_Dementia.RData")

library(survival) 
library(survminer)
library(mice)
library(splines)
library(gtools)

# Create function to calculate hazard ratio and 95% CI from pooled summary
# Using the following formula: exp(estimate) to calculate HR and exp(estimate -+ (1.96*se)) for the 95% CI. 
present.pooledcox <- function(x) { 
  coef <- exp(summary(pool(x))[,2])
  ci_low <- exp(summary(pool(x))[,2] - (1.96 * summary(pool(x))[,3]))
  ci_upp <- exp(summary(pool(x))[,2] + (1.96 * summary(pool(x))[,3]))
  pval <- summary(pool(x))[,6]
  pooledcox <- cbind(coef, ci_low, ci_upp, pval)
  rownames(pooledcox) <- rownames(summary(pool(x)))
  colnames(pooledcox) <- c("Hazard Ratio", "Lower .95", "Upper .95", "p-value")
  round(pooledcox, 3)
}


#______________________________________________________________________________________


#Calculate dementia cases ----
ergo_imp5 <- complete(ergo_imp, 5) #5th imputated model
summary(as.factor(ergo_imp5$dementia_incident))
summary(as.factor(ergo_imp5$Alzheimer))

#______________________________________________________________________________________
#Repeat analyses with Alzheimer disease as outcome ----

#PM10
PM10_Alzheimer <- with(data = ergo_imp, coxph(Surv(Follow_up_time, Alzheimer) ~ pm10tselect + Age_baseline + sex + 
                                           ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                           Depressive_Symptoms + Monthly_income))
summary(pool(PM10_Alzheimer))
present.pooledcox(PM10_Alzheimer)

#PM2.5
PM2.5_Alzheimer <- with(data = ergo_imp, coxph(Surv(Follow_up_time, Alzheimer) ~ pm25tselect + Age_baseline + sex + 
                                            ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                            Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Alzheimer))
present.pooledcox(PM2.5_Alzheimer)

#PM2.5 ab
PM2.5ab_Alzheimer <- with(data = ergo_imp, coxph(Surv(Follow_up_time, Alzheimer) ~ absorbancet + Age_baseline + sex + 
                                              ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                              Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Alzheimer))
present.pooledcox(PM2.5ab_Alzheimer)

#NOx
NOx_Alzheimer <- with(data = ergo_imp, coxph(Surv(Follow_up_time, Alzheimer) ~ noxtselecti + Age_baseline + sex + 
                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                          Depressive_Symptoms + Monthly_income))
summary(pool(NOx_Alzheimer))
present.pooledcox(NOx_Alzheimer)

#NO2
NO2_Alzheimer <- with(data = ergo_imp, coxph(Surv(Follow_up_time, Alzheimer) ~ no2tselecti + Age_baseline + sex + 
                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                          Depressive_Symptoms + Monthly_income))
summary(pool(NO2_Alzheimer))
present.pooledcox(NO2_Alzheimer)

#PCA
PCA_Alzheimer <- with(data = ergo_imp, coxph(Surv(Follow_up_time, Alzheimer) ~ Principle_all_polutants + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PCA_Alzheimer))
present.pooledcox(PCA_Alzheimer)



#______________________________________________________________________________________
#______________________________________________________________________________________
#Repeat analyses while time variable is censor age ----

#PM10
PM10_Censor_age <- with(data = ergo_imp, coxph(Surv(Censor_age, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                           ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                           Depressive_Symptoms + Monthly_income))
summary(pool(PM10_Censor_age))
present.pooledcox(PM10_Censor_age)

#PM2.5
PM2.5_Censor_age <- with(data = ergo_imp, coxph(Surv(Censor_age, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                            ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                            Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Censor_age))
present.pooledcox(PM2.5_Censor_age)

#PM2.5 ab
PM2.5ab_Censor_age <- with(data = ergo_imp, coxph(Surv(Censor_age, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                              ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                              Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Censor_age))
present.pooledcox(PM2.5ab_Censor_age)

#NOx
NOx_Censor_age <- with(data = ergo_imp, coxph(Surv(Censor_age, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                          Depressive_Symptoms + Monthly_income))
summary(pool(NOx_Censor_age))
present.pooledcox(NOx_Censor_age)

#NO2
NO2_Censor_age <- with(data = ergo_imp, coxph(Surv(Censor_age, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                          Depressive_Symptoms + Monthly_income))
summary(pool(NO2_Censor_age))
present.pooledcox(NO2_Censor_age)

#PCA
PCA_Censor_age <- with(data = ergo_imp, coxph(Surv(Censor_age, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(PCA_Censor_age))
present.pooledcox(PCA_Censor_age)


#______________________________________________________________________________________
#______________________________________________________________________________________
#Repeat analyses for APOE4 carriers and non-carriers seperately ----
ergo_imp_APOE4 <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_imp_APOE4 <- ergo_imp_APOE4 %>% filter(apoe4>0)
summary(as.factor(ergo_imp_APOE4$dementia_incident))/6
ergo_imp_APOE4 <- as.mids(ergo_imp_APOE4)

ergo_imp_APOEno <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_imp_APOEno <- ergo_imp_APOEno %>% filter(apoe4<1)
summary(as.factor(ergo_imp_APOEno$dementia_incident))/6
ergo_imp_APOEno <- as.mids(ergo_imp_APOEno)

#APOE carriers
#PM10
PM10_APOE <- with(data = ergo_imp_APOE4, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM10_APOE))
present.pooledcox(PM10_APOE)

#PM2.5
PM2.5_APOE <- with(data = ergo_imp_APOE4, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                  ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                  Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_APOE))
present.pooledcox(PM2.5_APOE)

#PM2.5 ab
PM2.5ab_APOE <- with(data = ergo_imp_APOE4, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                    ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                    Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_APOE))
present.pooledcox(PM2.5ab_APOE)

#NOx
NOx_APOE <- with(data = ergo_imp_APOE4, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                                ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(NOx_APOE))
present.pooledcox(NOx_APOE)

#NO2
NO2_APOE <- with(data = ergo_imp_APOE4, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                                ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(NO2_APOE))
present.pooledcox(NO2_APOE)

#PCA
PCA_APOE <- with(data = ergo_imp_APOE4, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(PCA_APOE))
present.pooledcox(PCA_APOE)




#APOE non-carriers
#PM10
PM10_APOEno <- with(data = ergo_imp_APOEno, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                    ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                    Depressive_Symptoms + Monthly_income))
summary(pool(PM10_APOEno))
present.pooledcox(PM10_APOEno)

#PM2.5
PM2.5_APOEno <- with(data = ergo_imp_APOEno, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_APOEno))
present.pooledcox(PM2.5_APOEno)

#PM2.5 ab
PM2.5ab_APOEno <- with(data = ergo_imp_APOEno, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_APOEno))
present.pooledcox(PM2.5ab_APOEno)

#NOx
NOx_APOEno <- with(data = ergo_imp_APOEno, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                                   ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(NOx_APOEno))
present.pooledcox(NOx_APOEno)

#NO2
NO2_APOEno <- with(data = ergo_imp_APOEno, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                                   ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(NO2_APOEno))
present.pooledcox(NO2_APOEno)

#PCA
PCA_APOEno <- with(data = ergo_imp_APOEno, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                   ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(PCA_APOEno))
present.pooledcox(PCA_APOEno)





#______________________________________________________________________________________
#______________________________________________________________________________________
#Repeat analyses for participants <70 years and >70 years seperately ----
ergo__70y <- mice::complete(ergo_imp, "long", include = TRUE)
ergo__70y <- ergo__70y %>% filter(Age_baseline<70)
summary(as.factor(ergo__70y$dementia_incident))/6
ergo__70y <- as.mids(ergo__70y)

ergo_70y <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_70y <- ergo_70y %>% filter(Age_baseline>69.999999999999999999999)
summary(as.factor(ergo_70y$dementia_incident))/6
ergo_70y <- as.mids(ergo_70y)

#<70 years
#PM10
PM10__70y <- with(data = ergo__70y, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                            ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                            Depressive_Symptoms + Monthly_income))
summary(pool(PM10__70y))
present.pooledcox(PM10__70y)

#PM2.5
PM2.5__70y <- with(data = ergo__70y, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                             ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5__70y))
present.pooledcox(PM2.5__70y)

#PM2.5 ab
PM2.5ab__70y <- with(data = ergo__70y, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5ab__70y))
present.pooledcox(PM2.5ab__70y)

#NOx
NOx__70y <- with(data = ergo__70y, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                           ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                           Depressive_Symptoms + Monthly_income))
summary(pool(NOx__70y))
present.pooledcox(NOx__70y)

#NO2
NO2__70y <- with(data = ergo__70y, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                           ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                           Depressive_Symptoms + Monthly_income))
summary(pool(NO2__70y))
present.pooledcox(NO2__70y)

#PCA
PCA__70y <- with(data = ergo__70y, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                           ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                           Depressive_Symptoms + Monthly_income))
summary(pool(PCA__70y))
present.pooledcox(PCA__70y)




#>70 years
#PM10
PM10_70y <- with(data = ergo_70y, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                          Depressive_Symptoms + Monthly_income))
summary(pool(PM10_70y))
present.pooledcox(PM10_70y)

#PM2.5
PM2.5_70y <- with(data = ergo_70y, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                           ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                           Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_70y))
present.pooledcox(PM2.5_70y)

#PM2.5 ab
PM2.5ab_70y <- with(data = ergo_70y, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                             ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_70y))
present.pooledcox(PM2.5ab_70y)

#NOx
NOx_70y <- with(data = ergo_70y, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                         ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                         Depressive_Symptoms + Monthly_income))
summary(pool(NOx_70y))
present.pooledcox(NOx_70y)

#NO2
NO2_70y <- with(data = ergo_70y, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                         ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                         Depressive_Symptoms + Monthly_income))
summary(pool(NO2_70y))
present.pooledcox(NO2_70y)

#PCA
PCA_70y <- with(data = ergo_70y, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                         ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                         Depressive_Symptoms + Monthly_income))
summary(pool(PCA_70y))
present.pooledcox(PCA_70y)



#______________________________________________________________________________________
#______________________________________________________________________________________
#Repeat analyses for participants with lower, intermediate and higher eduction seperately ----
ergo_primary <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_primary$level_of_education <- as.numeric(ergo_primary$level_of_education)
ergo_primary <- ergo_primary %>% filter(ergo_primary$level_of_education==1)
summary(as.factor(ergo_primary$dementia_incident))/6
ergo_primary <- as.mids(ergo_primary)

ergo_lower <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_lower$level_of_education <- as.numeric(ergo_lower$level_of_education)
ergo_lower <- ergo_lower %>% filter(ergo_lower$level_of_education==2)
summary(as.factor(ergo_lower$dementia_incident))/6
ergo_lower <- as.mids(ergo_lower)

ergo_Intermediate <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_Intermediate$level_of_education <- as.numeric(ergo_Intermediate$level_of_education)
ergo_Intermediate <- ergo_Intermediate %>% filter(ergo_Intermediate$level_of_education==3)
summary(as.factor(ergo_Intermediate$dementia_incident))/6
ergo_Intermediate <- as.mids(ergo_Intermediate)

ergo_higher <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_higher$level_of_education <- as.numeric(ergo_higher$level_of_education)
ergo_higher <- ergo_higher %>% filter(ergo_higher$level_of_education==4)
summary(as.factor(ergo_higher$dementia_incident))/6
ergo_higher <- as.mids(ergo_higher)



#Lower education
#PM10
PM10_lower <- with(data = ergo_lower, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                               From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                              Depressive_Symptoms + Monthly_income))
summary(pool(PM10_lower))
present.pooledcox(PM10_lower)

#PM2.5
PM2.5_lower <- with(data = ergo_lower, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                               From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_lower))
present.pooledcox(PM2.5_lower)

#PM2.5 ab
PM2.5ab_lower <- with(data = ergo_lower, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                 From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5ab_lower))
present.pooledcox(PM2.5ab_lower)

#NOx
NOx_lower <- with(data = ergo_lower, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                             From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(NOx_lower))
present.pooledcox(NOx_lower)

#NO2
NO2_lower <- with(data = ergo_lower, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                              From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(NO2_lower))
present.pooledcox(NO2_lower)

#PCA
PCA_lower <- with(data = ergo_lower, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                             From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(PCA_lower))
present.pooledcox(PCA_lower)






#Intermediate education
#PM10
PM10_Intermediate <- with(data = ergo_Intermediate, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                            From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                            Depressive_Symptoms + Monthly_income))
summary(pool(PM10_Intermediate))
present.pooledcox(PM10_Intermediate)

#PM2.5
PM2.5_Intermediate <- with(data = ergo_Intermediate, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                             From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                             Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Intermediate))
present.pooledcox(PM2.5_Intermediate)

#PM2.5 ab
PM2.5ab_Intermediate <- with(data = ergo_Intermediate, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                               From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                               Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Intermediate))
present.pooledcox(PM2.5ab_Intermediate)

#NOx
NOx_Intermediate <- with(data = ergo_Intermediate, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                                           From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                           Depressive_Symptoms + Monthly_income))
summary(pool(NOx_Intermediate))
present.pooledcox(NOx_Intermediate)

#NO2
NO2_Intermediate <- with(data = ergo_Intermediate, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                                           From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                           Depressive_Symptoms + Monthly_income))
summary(pool(NO2_Intermediate))
present.pooledcox(NO2_Intermediate)

#PCA
PCA_Intermediate <- with(data = ergo_Intermediate, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                           From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                           Depressive_Symptoms + Monthly_income))
summary(pool(PCA_Intermediate))
present.pooledcox(PCA_Intermediate)





#higher education
#PM10
PM10_higher <- with(data = ergo_higher, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(PM10_higher))
present.pooledcox(PM10_higher)

#PM2.5
PM2.5_higher <- with(data = ergo_higher, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                 From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_higher))
present.pooledcox(PM2.5_higher)

#PM2.5 ab
PM2.5ab_higher <- with(data = ergo_higher, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                   From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_higher))
present.pooledcox(PM2.5ab_higher)

#NOx
NOx_higher <- with(data = ergo_higher, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                               From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NOx_higher))
present.pooledcox(NOx_higher)

#NO2
NO2_higher <- with(data = ergo_higher, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                               From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NO2_higher))
present.pooledcox(NO2_higher)

#PCA
PCA_higher <- with(data = ergo_higher, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                               From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PCA_higher))
present.pooledcox(PCA_higher)





#______________________________________________________________________________________
#Repeat analyses for participants who are current, former and never smokers seperately ----
ergo_Current <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_Current$Smoking2 <- as.numeric(ergo_Current$Smoking2)
ergo_Current <- ergo_Current %>% filter(ergo_Current$Smoking2==1)
summary(as.factor(ergo_Current$dementia_incident))/6
ergo_Current <- as.mids(ergo_Current)

ergo_Former <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_Former$Smoking2 <- as.numeric(ergo_Former$Smoking2)
ergo_Former <- ergo_Former %>% filter(ergo_Former$Smoking2==2)
summary(as.factor(ergo_Former$dementia_incident))/6
ergo_Former <- as.mids(ergo_Former)

ergo_Never <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_Never$Smoking2 <- as.numeric(ergo_Never$Smoking2)
ergo_Never <- ergo_Never %>% filter(ergo_Never$Smoking2==3)
summary(as.factor(ergo_Never$dementia_incident))/6
ergo_Never <- as.mids(ergo_Never)

#Current smoker
#PM10
PM10_Current <- with(data = ergo_Current, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                  ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                  Depressive_Symptoms + Monthly_income))
summary(pool(PM10_Current))
present.pooledcox(PM10_Current)

#PM2.5
PM2.5_Current <- with(data = ergo_Current, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                   ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Current))
present.pooledcox(PM2.5_Current)

#PM2.5 ab
PM2.5ab_Current <- with(data = ergo_Current, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5ab_Current))
present.pooledcox(PM2.5ab_Current)

#NOx
NOx_Current <- with(data = ergo_Current, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(NOx_Current))
present.pooledcox(NOx_Current)

#NO2
NO2_Current <- with(data = ergo_Current, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(NO2_Current))
present.pooledcox(NO2_Current)

#PCA
PCA_Current <- with(data = ergo_Current, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PCA_Current))
present.pooledcox(PCA_Current)





#Former smoker
#PM10
PM10_Former <- with(data = ergo_Former, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(PM10_Former))
present.pooledcox(PM10_Former)

#PM2.5
PM2.5_Former <- with(data = ergo_Former, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Former))
present.pooledcox(PM2.5_Former)

#PM2.5 ab
PM2.5ab_Former <- with(data = ergo_Former, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                   ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Former))
present.pooledcox(PM2.5ab_Former)

#NOx
NOx_Former <- with(data = ergo_Former, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NOx_Former))
present.pooledcox(NOx_Former)

#NO2
NO2_Former <- with(data = ergo_Former, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NO2_Former))
present.pooledcox(NO2_Former)

#PCA
PCA_Former <- with(data = ergo_Former, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PCA_Former))
present.pooledcox(PCA_Former)



#Never smoker
#PM10
PM10_Never <- with(data = ergo_Never, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                              ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                              Depressive_Symptoms + Monthly_income))
summary(pool(PM10_Never))
present.pooledcox(PM10_Never)

#PM2.5
PM2.5_Never <- with(data = ergo_Never, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Never))
present.pooledcox(PM2.5_Never)

#PM2.5 ab
PM2.5ab_Never <- with(data = ergo_Never, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Never))
present.pooledcox(PM2.5ab_Never)

#NOx
NOx_Never <- with(data = ergo_Never, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                             ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(NOx_Never))
present.pooledcox(NOx_Never)

#NO2
NO2_Never <- with(data = ergo_Never, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                             ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(NO2_Never))
present.pooledcox(NO2_Never)

#PCA
PCA_Never <- with(data = ergo_Never, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                             ses_UNESCO_recoded + From_home + Physical_activity + Alcohol + BMI + 
                                             Depressive_Symptoms + Monthly_income))
summary(pool(PCA_Never))
present.pooledcox(PCA_Never)




#______________________________________________________________________________________
#______________________________________________________________________________________
#Repeat analyses while censoring for stroke and exclude participants with prevalent stroke ----
ergo_stroke <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_stroke <- ergo_stroke %>% filter(prev_stroke_2016 == "No prevalent stroke")
summary(as.factor(ergo_stroke$dementia_incident))/6
ergo_stroke <- as.mids(ergo_stroke)


#PM10
PM10_stroke <- with(data = ergo_stroke, coxph(Surv(Follow_up_Dementia_censor_Stroke, dementia_censor_Stroke) ~ pm10tselect + Age_baseline + sex + 
                                                ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(PM10_stroke))
present.pooledcox(PM10_stroke)

#PM2.5
PM2.5_stroke <- with(data = ergo_stroke, coxph(Surv(Follow_up_Dementia_censor_Stroke, dementia_censor_Stroke) ~ pm25tselect + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_stroke))
present.pooledcox(PM2.5_stroke)

#PM2.5 ab
PM2.5ab_stroke <- with(data = ergo_stroke, coxph(Surv(Follow_up_Dementia_censor_Stroke, dementia_censor_Stroke) ~ absorbancet + Age_baseline + sex + 
                                                   ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_stroke))
present.pooledcox(PM2.5ab_stroke)

#NOx
NOx_stroke <- with(data = ergo_stroke, coxph(Surv(Follow_up_Dementia_censor_Stroke, dementia_censor_Stroke) ~ noxtselecti + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NOx_stroke))
present.pooledcox(NOx_stroke)

#NO2
NO2_stroke <- with(data = ergo_stroke, coxph(Surv(Follow_up_Dementia_censor_Stroke, dementia_censor_Stroke) ~ no2tselecti + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NO2_stroke))
present.pooledcox(NO2_stroke)

#PCA
PCA_stroke <- with(data = ergo_stroke, coxph(Surv(Follow_up_Dementia_censor_Stroke, dementia_censor_Stroke) ~ Principle_all_polutants + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PCA_stroke))
present.pooledcox(PCA_stroke)






#______________________________________________________________________________________
#______________________________________________________________________________________
#Exclude those who changed home adress during the study ----
ergo_change <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_change <- ergo_change %>% filter(change_adress==0)
summary(as.factor(ergo_change$dementia_incident))/6
ergo_change <- as.mids(ergo_change)


#PM10
PM10_change <- with(data = ergo_change, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                Depressive_Symptoms + Monthly_income))
summary(pool(PM10_change))
present.pooledcox(PM10_change)

#PM2.5
PM2.5_change <- with(data = ergo_change, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                 ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_change))
present.pooledcox(PM2.5_change)

#PM2.5 ab
PM2.5ab_change <- with(data = ergo_change, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                   ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                   Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_change))
present.pooledcox(PM2.5ab_change)

#NOx
NOx_change <- with(data = ergo_change, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NOx_change))
present.pooledcox(NOx_change)

#NO2
NO2_change <- with(data = ergo_change, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(NO2_change))
present.pooledcox(NO2_change)

#PCA
PCA_change <- with(data = ergo_change, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                               ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                               Depressive_Symptoms + Monthly_income))
summary(pool(PCA_change))
present.pooledcox(PCA_change)



#______________________________________________________________________________________
#______________________________________________________________________________________

#Exclude those who moved to Ommoord <10 years and <25 years before baseline  ----
ergo_10Ommoord <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_10Ommoord <- ergo_10Ommoord %>% filter(Residence>10)
summary(as.factor(ergo_10Ommoord$dementia_incident))/6
ergo_10Ommoord <- as.mids(ergo_10Ommoord)

ergo_25Ommoord <- mice::complete(ergo_imp, "long", include = TRUE)
ergo_25Ommoord <- ergo_25Ommoord %>% filter(Residence>25)
summary(as.factor(ergo_25Ommoord$dementia_incident))/6
ergo_25Ommoord <- as.mids(ergo_25Ommoord)

#<10 years
#PM10
PM10_10Ommoord <- with(data = ergo_10Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(PM10_10Ommoord))
present.pooledcox(PM10_10Ommoord)

#PM2.5
PM2.5_10Ommoord <- with(data = ergo_10Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_10Ommoord))
present.pooledcox(PM2.5_10Ommoord)

#PM2.5 ab
PM2.5ab_10Ommoord <- with(data = ergo_10Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                         ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                         Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_10Ommoord))
present.pooledcox(PM2.5ab_10Ommoord)

#NOx
NOx_10Ommoord <- with(data = ergo_10Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(NOx_10Ommoord))
present.pooledcox(NOx_10Ommoord)

#NO2
NO2_10Ommoord <- with(data = ergo_10Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(NO2_10Ommoord))
present.pooledcox(NO2_10Ommoord)

#PCA
PCA_10Ommoord <- with(data = ergo_10Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(PCA_10Ommoord))
present.pooledcox(PCA_10Ommoord)



#<25 years
#PM10
PM10_25Ommoord <- with(data = ergo_25Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(PM10_25Ommoord))
present.pooledcox(PM10_25Ommoord)

#PM2.5
PM2.5_10Ommoord <- with(data = ergo_25Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_10Ommoord))
present.pooledcox(PM2.5_10Ommoord)

#PM2.5 ab
PM2.5ab_25Ommoord <- with(data = ergo_25Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                         ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                         Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_25Ommoord))
present.pooledcox(PM2.5ab_25Ommoord)

#NOx
NOx_25Ommoord <- with(data = ergo_25Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(NOx_25Ommoord))
present.pooledcox(NOx_25Ommoord)

#NO2
NO2_25Ommoord <- with(data = ergo_25Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(NO2_25Ommoord))
present.pooledcox(NO2_25Ommoord)

#PCA
PCA_25Ommoord <- with(data = ergo_25Ommoord, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                     Depressive_Symptoms + Monthly_income))
summary(pool(PCA_25Ommoord))
present.pooledcox(PCA_25Ommoord)


#______________________________________________________________________________________
#______________________________________________________________________________________


