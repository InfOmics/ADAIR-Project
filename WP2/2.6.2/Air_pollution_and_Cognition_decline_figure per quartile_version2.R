#Paper: Air pollution and dementia risk

#Written under R version 4.0.3

#______________________________________________________________________________________
#Package 
library(tidyverse) #includes packages that you're likely to use in everyday data anlsyes -> ggplot2, dplyr, tidyr, readr, purr, tibble, stringr, forcats. 
library(mice) 
library(broom) #builting fucntions
library(splitstackshape) #for long dataset 
library(nlme)
library(rio) # for the export function
library(splines)
library(nlme)
library(scales) # voor de doorzichtigheid van de confidence intervallen
library(ggh4x)


#______________________________________________________________________________________
#Open dataset ----
load("Dataset_air_pollution_and_Dementia.RData")

#______________________________________________________________________________________
#Create average dataset for analyses ----
#-- create new variable that identifies ergoid per imputation
ergo_imp <- mice::complete(ergo_imp, "long", include=TRUE)
imp_list_vol <- split(ergo_imp, ergo_imp$.imp)

imp_list_vol_2 <- lapply(imp_list_vol, function(DF)
{DF <- DF[order(DF$ergoid), ] #order by ergoid
DF$.id2 <- seq(nrow(DF)) #add new variable that counts from 1 to end of dataframe of imp == i
DF})

merge_vol <- bind_rows(imp_list_vol_2, .id="column_label")

#-- create numerical dataset to calculate average value
#- first transform all categorical factors into numerical
ergo_imp_vol <- merge_vol %>% 
  mutate_at(c("Principle_all_polutantsquartile", "no2quartile", "noxquartile", "pmabquartile", "pm2.5quartile", "pm10quartile", "Smoking2", "level_of_education", "column_label" ,"rs_cohort", "sex", "change_adress", "dementia_incident", "ses_UNESCO_recoded","Monthly_income", "Employed", "Smoking", "Alzheimer", "apoe4", 
              "Employed2", "Living_situation", "prev_stroke_2016", "incid_stroke_2016"), as.numeric)

#- then calculate average of all variables of all imputed datasets
ergo_imp_ave <- ergo_imp_vol %>% 
  group_by(.id2) %>% #indicating timepoint per ergoid per imputation
  summarise_all(mean, na.rm=T) %>% 
  mutate_at(c("Principle_all_polutantsquartile", "no2quartile", "noxquartile", "pmabquartile", "pm2.5quartile", "pm10quartile", "Smoking2", "level_of_education", "column_label" ,"rs_cohort", "sex", "change_adress", "dementia_incident", "ses_UNESCO_recoded","Monthly_income", "Employed", "Smoking", "Alzheimer", "apoe4",
              "Employed2", "Living_situation", "prev_stroke_2016", "incid_stroke_2016"), funs(round)) %>% 
  mutate_at(c("Principle_all_polutantsquartile", "no2quartile", "noxquartile", "pmabquartile", "pm2.5quartile", "pm10quartile", "Smoking2", "level_of_education", "column_label" ,"rs_cohort", "sex", "change_adress", "dementia_incident", "ses_UNESCO_recoded","Monthly_income", "Employed", "Smoking", "Alzheimer", "apoe4", 
              "Employed2", "Living_situation", "prev_stroke_2016", "incid_stroke_2016"), funs(factor)) 


#Convert ergo_imp1 to long format
# Convert to long format ----------------------------------------------------------------------
ergo_long1 <- merged.stack(ergo_imp_ave, var.stubs = c("MMSE_", "LDST_", "Stroop1_", "Stroop2_", "Stroop3_", "WFT_", "WLTimm_", "WLTdel_", "WLTrecog_", "PPB_both_", "PPB_right_", "PPB_left_", "PPB_sum_", "GFactor_", "timepoint_"), 
                           sep = "var.stubs") %>% rename(
                             MMSE = MMSE_,
                             LDST = LDST_,
                             Stroop1 = Stroop1_,
                             Stroop2 = Stroop2_,
                             Stroop3 = Stroop3_, 
                             WFT = WFT_,
                             WLTimm = WLTimm_,
                             WLTdel = WLTdel_,
                             WLTrecog = WLTrecog_,
                             PPB_both = PPB_both_, 
                             PPB_right = PPB_right_,
                             PPB_left = PPB_left_,
                             PPB_sum = PPB_sum_,
                             GFactor = GFactor_,
                             timepoint = timepoint_
                           ) # harmless warnings 


#############################
#------- Figures ----------#
############################


library(ggpubr)

effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}


setwd("V:/HomeDir/051543(T_de_Crom)/Air pollution and dementia/Air pollution and Dementia/Figure")

#-- Define names outcome variables. 
con_out_abs <- c("MMSE", "LDST", "Stroop1", "Stroop2", "Stroop3", "WFT", "WLTimm", "WLTdel", "WLTrecog", "PPB_sum", "GFactor")
con_nam_abs <- c("Mini-Mental State Examination score" ,"Letter-Digit-Substitution Test score", "Stroop Test: Reading score" , "Stroop Test: Naming score", "Stroop Test: Interference score",
                 "Word Fluency Test score", "Word Learning Test: Immediate recall score", "Word Learning Test: Delayed recall score", "Word Learning Test: Recognition score", "Purdue Pegboard Test score", "G-Factor")



#--- to create figures, you must create a new dataframe with the 'average person' regarding your covariates.
# you can do this by using the mean for continuous variables and choosing a value for categorical variables, 
# or also make categorical variables continuous and subsequently take the average (sort of the prevalence)

#create function in which we transform categorical variables into continuous variables. 
Factor_cog_function <- function (dataset, out, out_nam) { 
  listoffms <- list()
  for (i in 1:length(out)){
    y <- out[i]
    if(out[i]=="Stroop1"|out[i]=="Stroop3"|out[i]=="WLTdel"| out[i]=="WLTrecog"){
      fm <- nlme::lme(fixed = formula(paste(y, "~Principle_all_polutantsquartile*timepoint + ns(Age_ERGO5, 2) + as.numeric(sex) + 
                    as.numeric(ses_UNESCO_recoded) + From_home + as.numeric(Smoking) + Physical_activity + Alcohol + BMI + 
                    Depressive_Symptoms + as.numeric(Monthly_income)", sep = "")), 
                      data = dataset,
                      random = ~ 1 |ergoid,
                      method = "REML",
                      na.action = na.omit)}
    else{ 
      fm <- nlme::lme(fixed = formula(paste(y, "~Principle_all_polutantsquartile*timepoint + ns(Age_ERGO5, 2) + as.numeric(sex) + 
                    as.numeric(ses_UNESCO_recoded) + From_home + as.numeric(Smoking) + Physical_activity + Alcohol + BMI + 
                    Depressive_Symptoms + as.numeric(Monthly_income)", sep = "")), 
                      data = dataset,
                      random = ~ timepoint |ergoid,
                      method = "REML",
                      na.action = na.omit)
      
    }
    listoffms[[i]] <- fm  # Put all the models into a list
  }
  names(listoffms) <- out_nam # and give them the names defined by out_nam
  return(listoffms)
}


#call upon function
list_factor_cog <- Factor_cog_function(ergo_long1, out = con_out_abs, out_nam = con_nam_abs)
list_factor_cog

#---  create an average dataframe to use for predictions---------
newDF <- with(ergo_long1, expand.grid(
  timepoint = seq(0, 6.998, length.out = 10), #set follow up to maximum in your sample
  Age_ERGO5 = mean(ergo_long1$Age_ERGO5, na.rm=T),
  sex = mean(as.numeric(ergo_long1$sex)),
  ses_UNESCO_recoded =  mean(as.numeric(ergo_long1$ses_UNESCO_recoded)),
  From_home = mean(ergo_long1$From_home) ,
  Smoking = mean(as.numeric(ergo_long1$Smoking)),
  Physical_activity = mean(ergo_long1$Physical_activity),
  Alcohol = mean(ergo_long1$Alcohol),
  BMI = mean(ergo_long1$BMI), 
  Depressive_Symptoms = mean(ergo_long1$Depressive_Symptoms), #use na.rm = T if you have missings in your dataset
  Monthly_income = mean(as.numeric(ergo_long1$Monthly_income)), 
  Principle_all_polutantsquartile = levels(Principle_all_polutantsquartile)))                            


#------ create prediction dataframe ---------
#first set some variables
list_models <- list_factor_cog
model_names <- names(list_models)
trajectory_colors <- c("#be5a64", "##669b90", "#916ebd", "#d39c2a") # from https://medialab.github.io/iwanthue/
linewidth <- 4
trajectory_names <- con_nam_abs
new_df <- newDF
out <- con_out_abs

ergo_long1$ergoid <- as.numeric(ergo_long1$ergoid)
#-- rename model names, to fit in image (white matter volume iso normal appearing white matter, write in description)
y_names <- c("Mini-Mental State Examination score" ,"Letter-Digit-Substitution Test score", "Stroop Test: Reading (seconds)" , "Stroop Test: Naming (seconds)", "Stroop Test: Interference (seconds)",
             "Word Fluency Test score", "Word Learning Test: Immediate recall score", "Word Learning Test: Delayed recall score", "Word Learning Test: Recognition score", "Purdue Pegboard Test score", "G-Factor")


# run for loop to create prediction dataframe "new_df"
for(i in (1:length(model_names))){ #for each model/outcome 
  y <- out[i]
  predict <- effectPlotData(list_models[[model_names[[i]]]], newdata = newDF, ergo_long1) #predict the values for the newDF
  predict <- predict %>%
    dplyr::select(c(pred, low, upp)) #select only the newly predicted values
  # orient in same directions;
  # if (model_names[i]=='White matter hyperintensities volume'|model_names[i]=='Mean diffusivity'){ # als de lijn stijgt, doe keer -1
  #   predict$pred <- predict$pred*-1
  #   predict$upp <- predict$upp*-1
  #   predict$low <- predict$low*-1
  # }
  colnames(predict) <- paste0(c("pred", "low", "upp"), "_", out[i]) # rename variables by adding "_outcome"
  new_df <- cbind(new_df, predict) #bind these predicted variables to the prediction dataframe new_df
}
predict
pred_vector <- tidyselect::vars_select(colnames(new_df), starts_with("pred")) #select all pred_outcome variables
upr_vector <- tidyselect::vars_select(colnames(new_df), starts_with("upp")) #select all upp_outcome variables
lwr_vector <- tidyselect::vars_select(colnames(new_df), starts_with("low")) #select all low_outcome variables





new_vector <- seq(1, length(model_names), by = 1) #create a vector (1 2 3 4 5 ...) depending on the number of outcomes


myplots_factor_MMSE = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[1], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ # create a ggplot with follow-up time on x-axis, and prediction on y-axis, grouped by quartiles of air pollution value
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ #the lines have different colors based on the quartiles
    theme_light() + ggtitle("Mini-Mental State Examination") + theme(plot.title = element_text(hjust = 0.5)) +  theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    ylab(paste0("\n", y_names[1])) + #label for y axis in a vertical manner
    xlab("Follow-up time (years)") + 
    geom_ribbon( aes_string(ymin = lwr_vector[1], ymax = upr_vector[1], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1,
                 show.legend = FALSE) + #I create a ribbon for 95%CI use lower and upper values, alpha states transparency 
    scale_y_continuous(limits = c(26.65, 29), breaks = c(27, 28, 29)) + #set limits for the y-axis, so all figures look the same
    scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + #set limits for the x-axis,so all figures look the same
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(26.65, Inf), 
                                   trunc_lower=c(-Inf, 26.8))) + 
    annotate("text", x=-Inf, y=c(26.65, 26.8), label="/", angle=90) +
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) #give manual names to the legend
)
 

myplots_factor_MMSE[[1]]

#LDST
myplots_factor_LDST = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[2], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[2])) + 
    xlab("Follow-up time (years)") + ggtitle("Letter-Digit-Substitution Test") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[2], ymax = upr_vector[2], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(23.4, 30)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(23.4, Inf), 
                                   trunc_lower=c(-Inf, 23.8))) + 
    annotate("text", x=-Inf, y=c(23.4, 23.8), label="/", angle=90) + 
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")))

myplots_factor_LDST[[2]]


#Stroop1
myplots_factor_Stroop1 = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[3], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[3])) + 
    xlab("Follow-up time (years)") + ggtitle("Stroop Test: Reading") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[3], ymax = upr_vector[3], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(-20.4, -17)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(-20.4, Inf), 
                                   trunc_lower=c(-Inf, -20.2))) + 
    annotate("text", x=-Inf, y=c(-20.4, -20.2), label="/", angle=90) + 
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_Stroop1[[3]]

#Stroop2
myplots_factor_Stroop2 = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[4], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[4])) + 
    xlab("Follow-up time (years)") + ggtitle("Stroop Test: Naming") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[4], ymax = upr_vector[4], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(-27.5, -23)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(-27.5, Inf), 
                                   trunc_lower=c(-Inf, -27.2))) + 
    annotate("text", x=-Inf, y=c(-27.5, -27.2), label="/", angle=90) + 
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_Stroop2[[4]]

#Stroop3
myplots_factor_Stroop3 = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[5], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[5])) + 
    xlab("Follow-up time (years)") + ggtitle("Stroop Test: Interference") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[5], ymax = upr_vector[5], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(-71.7, -45)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(-71.7, Inf), 
                                   trunc_lower=c(-Inf, -70.5))) + 
    annotate("text", x=-Inf, y=c(-71.7, -70.5), label="/", angle=90) +
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_Stroop3[[5]]

#WFT
myplots_factor_WFT = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[6], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75) + theme_light() + ylab(paste0("\n", y_names[6])) + 
    xlab("Follow-up time (years)") + ggtitle("Word Fluency Test") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[6], ymax = upr_vector[6], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(19.5, 23)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(19.5, Inf), 
                                   trunc_lower=c(-Inf, 19.7))) + 
    annotate("text", x=-Inf, y=c(19.5, 19.7), label="/", angle=90) +
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_WFT[[6]]


#WLT: Immediate recall score
myplots_factor_WLT_imm = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[7], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ # 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[7])) +
    xlab("Follow-up time (years)") + ggtitle("Word Learning Test: Immediate recall") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[7], ymax = upr_vector[7], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(5.6, 9)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) +
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(5.6, Inf), 
                                   trunc_lower=c(-Inf, 5.8))) + 
    annotate("text", x=-Inf, y=c(5.6, 5.8), label="/", angle=90) +
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_WLT_imm[[7]]


#WLT: Delayed recall score
myplots_factor_WLT_del = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[8], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[8])) + 
    xlab("Follow-up time (years)") + ggtitle("Word Learning Test: Delayed recall") + theme(plot.title = element_text(hjust = 0.5)) +theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[8], ymax = upr_vector[8], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(3.5, 9)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(3.5, Inf), 
                                   trunc_lower=c(-Inf, 3.8))) + 
    annotate("text", x=-Inf, y=c(3.5, 3.8), label="/", angle=90) +
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_WLT_del[[8]]

#WLT: Recognition score
myplots_factor_WLT_rec = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[9], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[9])) + 
    xlab("Follow-up time (years)") + ggtitle("Word Learning Test: Recognition") + theme(plot.title = element_text(hjust = 0.5)) +theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[9], ymax = upr_vector[9], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(11.5, 15)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(11.5, Inf), 
                                   trunc_lower=c(-Inf, 11.7))) + 
    annotate("text", x=-Inf, y=c(11.5, 11.7), label="/", angle=90) +
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) #give manual names to the legend
)
myplots_factor_WLT_rec[[9]]

#Purdue Pegboard Test
myplots_factor_PPB = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[10], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[10])) + 
    xlab("Follow-up time (years)") + ggtitle("Purdue Pegboard Test") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[10], ymax = upr_vector[10], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(29.5, 36)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(29.5, Inf), 
                                   trunc_lower=c(-Inf, 29.8))) + 
    annotate("text", x=-Inf, y=c(29.5, 29.8), label="/", angle=90) +
   scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_PPB[[10]]

#G-Factor
myplots_factor_GFactor = lapply(new_vector, function(i) #lapply is somewhat similar to for loop, for each variable in the new_vector
  ggplot(new_df, aes_string(x="timepoint", y=pred_vector[11], group="Principle_all_polutantsquartile", fill="factor(Principle_all_polutantsquartile)"))+ 
    geom_line(mapping = aes(color=factor(Principle_all_polutantsquartile)), size =0.75)+ theme_light()+ ylab(paste0("\n", y_names[11])) + 
    xlab("Follow-up time (years)") + ggtitle("G-Factor") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_ribbon( aes_string(ymin = lwr_vector[11], ymax = upr_vector[11], fill = "factor(Principle_all_polutantsquartile)", color = NULL), alpha = .1, show.legend = FALSE) + 
    scale_y_continuous(limits = c(-0.88, 0.4), breaks = c(-0.8,-0.6, -0.4, -0.2, -0.0, 0.2,0.4)) + scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7)) + 
    theme(axis.text = element_text(size=8), axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10), title = element_text(size = 10),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black")) + 
    guides(y=guide_axis_truncated( trunc_upper=c(-0.88, Inf), 
                                   trunc_lower=c(-Inf, -0.83))) + 
    annotate("text", x=-Inf, y=c(-0.88, -0.83), label="/", angle=90) +
    scale_colour_discrete(name = "General marker of all air pollutants", labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4")) 
)
myplots_factor_GFactor[[11]]



