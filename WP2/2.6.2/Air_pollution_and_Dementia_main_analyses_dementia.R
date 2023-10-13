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
#Create study population by excluding all participants who are <55 years (n=662) ----
#ergo_imp <- mice::complete(ergo_imp, "long", include = TRUE)
#ergo_imp <- ergo_imp %>% filter(Age_baseline>54.999999999999)
#ergo_imp <- as.mids(ergo_imp)



#______________________________________________________________________________________


#Calculate dementia cases ----
ergo_imp5 <- complete(ergo_imp, 5) #5th imputated model
summary(as.factor(ergo_imp5$dementia_incident))

#______________________________________________________________________________________
#Table 2: PM10 and dementia risk, using cox models:----

#Model 1: 
PM10_dementia_m1 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex))
summary(pool(PM10_dementia_m1))
present.pooledcox(PM10_dementia_m1)


#_________________________________________________________________________________________________________
#Model 2: 
PM10_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                              ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                 Depressive_Symptoms + Monthly_income))
summary(pool(PM10_dementia_m2))
present.pooledcox(PM10_dementia_m2)


#_________________________________________________________________________________________________________
#Check assumptions using the 5th imputed model since we can not check assumptions using a pooled model. 
#Explanation of assumptions: (https://www.mwsug.org/proceedings/2006/stats/MWSUG-2006-SD08.pdf)
#Assumption 1: test proportional hazard assumption using Schoenfeld residuals
ergo_imp5 <- complete(ergo_imp, 5) #Test assumptions using the 5th imputated model
Fit_PM10_dementia_m2 <- (coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect + Age_baseline + sex + 
                                  ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                  Depressive_Symptoms + Monthly_income,
                               data = ergo_imp5)) #conduct cox model for 5th imputated model
Fit_PM10_dementia_m2


Fit_PM10_dementia_m2 <- cox.zph(Fit_PM10_dementia_m2, transform="km") #test proportional hazard assumption
Fit_PM10_dementia_m2
ggcoxzph(Fit_PM10_dementia_m2) #Proportional hazard is checked




#______________________________________________________________________________________
#Table 2: PM2.5 and dementia risk, using cox models:----

#Model 1: 
PM2.5_dementia_m1 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex))
summary(pool(PM2.5_dementia_m1))
present.pooledcox(PM2.5_dementia_m1)

#_________________________________________________________________________________________________________
#Model 2: 
PM2.5_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                                    ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                    Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_dementia_m2))
present.pooledcox(PM2.5_dementia_m2)


#_________________________________________________________________________________________________________

#Check assumptions using the 5th imputed model since we can not check assumptions using a pooled model. 
#Explanation of assumptions: (https://www.mwsug.org/proceedings/2006/stats/MWSUG-2006-SD08.pdf)
#Assumption 1: test proportional hazard assumption using Schoenfeld residuals
Fit_PM2.5_dementia_m2 <- (coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect + Age_baseline + sex + 
                                   ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                   Depressive_Symptoms + Monthly_income,
                                data = ergo_imp5)) #conduct cox model for 5th imputated model 
Fit_PM2.5_dementia_m2


Fit_PM2.5_dementia_m2 <- cox.zph(Fit_PM2.5_dementia_m2, transform="km") #test proportional hazard assumption
Fit_PM2.5_dementia_m2
ggcoxzph(Fit_PM2.5_dementia_m2) #Proportional hazard is checked


#______________________________________________________________________________________
#Table 2: PM2.5ab and dementia risk, using cox models:----

#Model 1: 
PM2.5ab_dementia_m1 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet  + Age_baseline + sex))
summary(pool(PM2.5ab_dementia_m1))
present.pooledcox(PM2.5ab_dementia_m1)

#_________________________________________________________________________________________________________
#Model 2: 
PM2.5ab_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5ab_dementia_m2))
present.pooledcox(PM2.5ab_dementia_m2)

#_________________________________________________________________________________________________________

#Check assumptions using the 5th imputed model since we can not check assumptions using a pooled model. 
#Explanation of assumptions: (https://www.mwsug.org/proceedings/2006/stats/MWSUG-2006-SD08.pdf)
#Assumption 1: test proportional hazard assumption using Schoenfeld residuals
Fit_PM2.5ab_dementia_m2 <- (coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet + Age_baseline + sex + 
                                     ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                     Depressive_Symptoms + Monthly_income,
                                data = ergo_imp5)) #conduct cox model for 5th imputated model



Fit_PM2.5ab_dementia_m2 <- cox.zph(Fit_PM2.5ab_dementia_m2, transform="km") #test proportional hazard assumption
Fit_PM2.5ab_dementia_m2
ggcoxzph(Fit_PM2.5ab_dementia_m2) #Proportional hazard is checked










#______________________________________________________________________________________
#Table 2: NOx and dementia risk, using cox models:----

#Model 1: 
NOx_dementia_m1 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex))
summary(pool(NOx_dementia_m1))
present.pooledcox(NOx_dementia_m1)

#_________________________________________________________________________________________________________
#Model 2: 
NOx_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                                  ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                  Depressive_Symptoms + Monthly_income))
summary(pool(NOx_dementia_m2))
present.pooledcox(NOx_dementia_m2)

#_________________________________________________________________________________________________________


#Check assumptions using the 5th imputed model since we can not check assumptions using a pooled model. 
#Explanation of assumptions: (https://www.mwsug.org/proceedings/2006/stats/MWSUG-2006-SD08.pdf)
#Assumption 1: test proportional hazard assumption using Schoenfeld residuals
Fit_NOx_dementia_m2 <- (coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti + Age_baseline + sex + 
                                 ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                 Depressive_Symptoms + Monthly_income,
                                  data = ergo_imp5)) #conduct cox model for 5th imputated model



Fit_NOx_dementia_m2 <- cox.zph(Fit_NOx_dementia_m2, transform="km") #test proportional hazard assumption
Fit_NOx_dementia_m2
ggcoxzph(Fit_NOx_dementia_m2) #Proportional hazard is checked













#______________________________________________________________________________________
#Table 2: NO2 and dementia risk, using cox models:----

#Model 1: 
NO2_dementia_m1 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex))
summary(pool(NO2_dementia_m1))
present.pooledcox(NO2_dementia_m1)

#_________________________________________________________________________________________________________
#Model 2: 
NO2_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                                  ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                  Depressive_Symptoms + Monthly_income))
summary(pool(NO2_dementia_m2))
present.pooledcox(NO2_dementia_m2)

#_________________________________________________________________________________________________________



#Check assumptions using the 5th imputed model since we can not check assumptions using a pooled model. 
#Explanation of assumptions: (https://www.mwsug.org/proceedings/2006/stats/MWSUG-2006-SD08.pdf)
#Assumption 1: test proportional hazard assumption using Schoenfeld residuals
Fit_NO2_dementia_m2 <- (coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti + Age_baseline + sex + 
                                 ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                 Depressive_Symptoms + Monthly_income,
                              data = ergo_imp5)) #conduct cox model for 5th imputated model



Fit_NO2_dementia_m2 <- cox.zph(Fit_NO2_dementia_m2, transform="km") #test proportional hazard assumption
Fit_NO2_dementia_m2
ggcoxzph(Fit_NO2_dementia_m2) #Proportional hazard is checked




#______________________________________________________________________________________
#Table 2: PCA all and dementia risk, using cox models:----

#Model 1: 
PCA_all_dementia_m1 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex))
summary(pool(PCA_all_dementia_m1))
present.pooledcox(PCA_all_dementia_m1)

#_________________________________________________________________________________________________________
#Model 2: 
PCA_all_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutants + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(PCA_all_dementia_m2))
present.pooledcox(PCA_all_dementia_m2)







#_________________________________________________________________________________________________________

#_________________________________________________________________________________________________________
#Air pollutants in quartiles and dementia risk, using cox models (model 2 only):----

#PM10
PM10abquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10quartile + Age_baseline + sex + 
                                                             ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                             Depressive_Symptoms + Monthly_income))
summary(pool(PM10abquartile_dementia_m2))
present.pooledcox(PM10abquartile_dementia_m2)

#PM2.5
PM2.5quartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm2.5quartile + Age_baseline + sex + 
                                                            ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                            Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5quartile_dementia_m2))
present.pooledcox(PM2.5quartile_dementia_m2)

#PM2.5 ab
PM2.5abquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pmabquartile + Age_baseline + sex + 
                                                              ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                              Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5abquartile_dementia_m2))
present.pooledcox(PM2.5abquartile_dementia_m2)

#NOx
NOxquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxquartile + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(NOxquartile_dementia_m2))
present.pooledcox(NOxquartile_dementia_m2)

#NO2
NO2quartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2quartile + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(NO2quartile_dementia_m2))
present.pooledcox(NO2quartile_dementia_m2)


#All air 
PCAquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutantsquartile + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(PCAquartile_dementia_m2))
present.pooledcox(PCAquartile_dementia_m2)




#_________________________________________________________________________________________________________
#Air pollutants in tertiles and dementia risk, using cox models (model 2 only):----
ergo_imp <- mice::complete(ergo_imp, "long", include = TRUE)

#Create quartile variables for air pollutants
ergo_imp$pm10tertile <- quantcut(ergo_imp$pm10tselect, q=3, na.rm=TRUE)
ergo_imp$pm2.5tertile <- quantcut(ergo_imp$pm25tselect, q=3, na.rm=TRUE)
ergo_imp$pmabtertile <- quantcut(ergo_imp$absorbancet, q=3, na.rm=TRUE)
ergo_imp$noxtertile <- quantcut(ergo_imp$noxtselecti, q=3, na.rm=TRUE)
ergo_imp$no2tertile <- quantcut(ergo_imp$no2tselecti, q=3, na.rm=TRUE)
ergo_imp$Principle_all_polutantstertile <- quantcut(ergo_imp$Principle_all_polutants, q=3, na.rm=TRUE)

ergo_imp <- as.mids(ergo_imp)

#PM10
PM10abquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tertile + Age_baseline + sex + 
                                                             ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                             Depressive_Symptoms + Monthly_income))
summary(pool(PM10abquartile_dementia_m2))
present.pooledcox(PM10abquartile_dementia_m2)

#PM2.5
PM2.5quartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm2.5tertile + Age_baseline + sex + 
                                                            ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                            Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5quartile_dementia_m2))
present.pooledcox(PM2.5quartile_dementia_m2)

#PM2.5 ab
PM2.5abquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pmabtertile + Age_baseline + sex + 
                                                              ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                              Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5abquartile_dementia_m2))
present.pooledcox(PM2.5abquartile_dementia_m2)

#NOx
NOxquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtertile + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(NOxquartile_dementia_m2))
present.pooledcox(NOxquartile_dementia_m2)

#NO2
NO2quartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tertile + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(NO2quartile_dementia_m2))
present.pooledcox(NO2quartile_dementia_m2)


#All air 
PCAquartile_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ Principle_all_polutantstertile + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(PCAquartile_dementia_m2))
present.pooledcox(PCAquartile_dementia_m2)


#______________________________________________________________________________________
#Test for interaction for sex on model 2:----

#PM10
PM10_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect*sex + Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(PM10_sex_dementia_m2))
present.pooledcox(PM10_sex_dementia_m2)

#PM2.5
PM2.5_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect*sex + Age_baseline + sex + 
                                                        ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                        Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_sex_dementia_m2))
present.pooledcox(PM2.5_sex_dementia_m2)

#PM2.5 ab
PM2.5ab_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet*sex + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_sex_dementia_m2))
present.pooledcox(PM2.5ab_sex_dementia_m2)

#NOx
NOx_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti*sex + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(NOx_sex_dementia_m2))
present.pooledcox(NOx_sex_dementia_m2)

#NO2
NO2_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti*sex + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(NO2_sex_dementia_m2))
present.pooledcox(NO2_sex_dementia_m2)


#______________________________________________________________________________________
#Test for interaction for smoking on model 2:----

#PM10
PM10_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect*Smoking +sex + Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(PM10_sex_dementia_m2))
present.pooledcox(PM10_sex_dementia_m2)

#PM2.5
PM2.5_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect*Smoking + sex + Age_baseline + sex + 
                                                        ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                        Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_sex_dementia_m2))
present.pooledcox(PM2.5_sex_dementia_m2)

#PM2.5 ab
PM2.5ab_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet*Smoking + sex + Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_sex_dementia_m2))
present.pooledcox(PM2.5ab_sex_dementia_m2)

#NOx
NOx_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti*Smoking + sex + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(NOx_sex_dementia_m2))
present.pooledcox(NOx_sex_dementia_m2)

#NO2
NO2_sex_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti*Smoking + sex + Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(NO2_sex_dementia_m2))
present.pooledcox(NO2_sex_dementia_m2)

#______________________________________________________________________________________
#Test for interaction for APOE on model 2:----

#PM10
PM10_apoe_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect*apoe4 + sex + Age_baseline + sex + 
                                                        ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                        Depressive_Symptoms + Monthly_income))
summary(pool(PM10_apoe_dementia_m2))
present.pooledcox(PM10_apoe_dementia_m2)

#PM2.5
PM2.5_apoe_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect*apoe4 + sex + Age_baseline + sex + 
                                                         ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                         Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_apoe_dementia_m2))
present.pooledcox(PM2.5_apoe_dementia_m2)

#PM2.5 ab
PM2.5ab_apoe_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet*apoe4 + sex + Age_baseline + sex + 
                                                           ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                           Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_apoe_dementia_m2))
present.pooledcox(PM2.5ab_apoe_dementia_m2)

#NOx
NOx_apoe_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti*apoe4 + sex + Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(NOx_apoe_dementia_m2))
present.pooledcox(NOx_apoe_dementia_m2)

#NO2
NO2_apoe_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti*apoe4 + sex + Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(NO2_apoe_dementia_m2))
present.pooledcox(NO2_apoe_dementia_m2)



#______________________________________________________________________________________
#Test for interaction for age on model 2:----

#PM10
PM10_Age_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm10tselect*Age_baseline + sex + 
                                                       ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                       Depressive_Symptoms + Monthly_income))
summary(pool(PM10_Age_dementia_m2))
present.pooledcox(PM10_Age_dementia_m2)

#PM2.5
PM2.5_Age_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ pm25tselect*Age_baseline + sex + 
                                                        ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                        Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Age_dementia_m2))
present.pooledcox(PM2.5_Age_dementia_m2)

#PM2.5 ab
PM2.5ab_Age_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ absorbancet*Age_baseline + sex + 
                                                          ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                          Depressive_Symptoms + Monthly_income))
summary(pool(PM2.5_Age_dementia_m2))
present.pooledcox(PM2.5ab_Age_dementia_m2)

#NOx
NOx_Age_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ noxtselecti*Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(NOx_Age_dementia_m2))
present.pooledcox(NOx_Age_dementia_m2)

#NO2
NO2_Age_dementia_m2 <- with(data = ergo_imp, coxph(Surv(Follow_up_time, dementia_incident) ~ no2tselecti*Age_baseline + sex + 
                                                      ses_UNESCO_recoded + From_home + Smoking + Physical_activity + Alcohol + BMI + 
                                                      Depressive_Symptoms + Monthly_income))
summary(pool(NO2_Age_dementia_m2))
present.pooledcox(NO2_Age_dementia_m2)

