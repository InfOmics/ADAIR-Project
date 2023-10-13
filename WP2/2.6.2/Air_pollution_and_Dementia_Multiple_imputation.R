#Paper: air pollution and dementia risk


#Written under R version 4.0.3

#______________________________________________________________________________________
# Open dataset for analysis named ergo ----

library(foreign)
library(reshape2)
library(gtools)

ergo <- read.spss("Air_pollution_and_Dementia_Dataset.sav", use.value.label=TRUE, to.data.frame=TRUE)
str(ergo)


#attach dataset ergo
attach(ergo)


#create outcome variable for cox model
ergo$dementia_incident <- as.numeric(ergo$dementia_incident)
ergo$Alzheimer <- as.numeric(ergo$Alzheimer)

#Create second variable for level of education and smoking, one to adjust for and one for stratification.
ergo$level_of_education <- ergo$ses_UNESCO_recoded
ergo$Smoking2 <- ergo$Smoking




#______________________________________________________________________________________
#Imputate missing values of covariates ----
library(mice)

#Compute number of complete and incomplete cases. 
Ncc <- cbind(
   "#" = table(complete.cases(ergo)),
   "%" = round(100 * table(complete.cases(ergo))/nrow(ergo), 2)
)
rownames(Ncc) <- c("incompl.", "complete")
Ncc #Number of incomplete cases is high (100%) which is due to e7_adresidentiek

# Number and proportion of missing values per variable
cbind("# NA" = sort(colSums(is.na(ergo))),
      "% NA" = round(sort(colMeans(is.na(ergo))) * 100, 2)) #multiplied by 100 to obtain percentage and 2 decimale


# set up run of imputation
imp <- mice(ergo, maxit = 0) #1 iteractions. 
imp$loggedEvents #No collinearity between variables that I want to imputate 
meth <- imp$method
predM <- imp$predictorMatrix
pred <- imp$pred


# Setting values of variables you do not want to use as predictor to 0 in the imputation matrix.
predM[, c("ergoid", "rs_cohort", "pm10tselect", "pm25tselect", "absorbancet", "noxtselecti", "no2tselecti", "Principle_all_polutants", 
          "Principle_polutants", "change_adress", "dementia_incident", "Follow_up_time", "Censor_age", "Alzheimer", "Follow_up_time_stroke", 
          "dementia_censor_Stroke", "timepoint_1", "Smoking2", "Employed2", "Residence", "level_of_education", "incid_stroke_2016", "apoe4", 
          "Work_hours", "Age_ERGO5", "timediff", "Cognitive_1", "MMSE_1", "Cognitive_1_MMSE", "Physical_activity", "WFT_1", "Hours_home", 
          "Stroop2_1", "Stroop1_1", "LDST_1", "Stroop3_1", "PPB_left_1", "WLTrecog_1", "PPB_both_1", "PPB_right_1", "WLTimm_1", "WLTdel_1", 
          "PPB_sum_1", "GFactor_1", "Cognitive_2", "timepoint_2", "MMSE_2", "Cognitive_2_MMSE", "WFT_2", "LDST_2", "Stroop1_2", "WLTdel_2", 
          "Stroop3_2", "WLTrecog_2", "PPB_left_2", "PPB_right_2", "PPB_both_2", "PPB_sum_2", "GFactor_2", "Stroop2_2", "WLTimm_2", 
          "Monthly_income", "Living_situation")] <- 0
#Excluded the dietary scores and included the food components. 

# list any variables you do not want to be imputed.
meth[c("Smoking2", "Employed2", "Residence", "level_of_education", "incid_stroke_2016", "Living_situation", "apoe4", "Age_ERGO5","timediff", 
       "Cognitive_1", "MMSE_1", "Cognitive_1_MMSE", "WFT_1", "Stroop2_1", "Stroop1_1", "LDST_1", "Stroop3_1", "PPB_left_1", 
       "WLTrecog_1", "PPB_both_1", "PPB_right_1", "WLTimm_1", "WLTdel_1", "PPB_sum_1", "GFactor_1", "Cognitive_2", "timepoint_2", "MMSE_2", 
       "Cognitive_2_MMSE", "WFT_2", "LDST_2", "Stroop1_2", "WLTdel_2", "Stroop3_2", "WLTrecog_2", "PPB_left_2", "PPB_right_2", "PPB_both_2", 
       "PPB_sum_2", "GFactor_2", "Stroop2_2", "WLTimm_2")]=""
meth

# set imputation
ergo_imp <- mice(ergo, method=meth, predictorMatrix=predM, maxit = 10, m=5, seed=2020) 

#check logged events
#removed predictors because of multicolinearity is noted by loggedEvents
ergo_imp$loggedEvents

# check the number of missings in the first imputed dataset, to see if nothing was skipped.
cbind("# NA" = sort(colSums(is.na(complete(ergo_imp, 1)))),
      "% NA" = round(sort(colMeans(is.na(complete(ergo_imp, 1)))) * 100, 2))

#Check original database vs imputed database
#Ceck continous variables
densityplot(ergo_imp, ~ Depressive_Symptoms  + Alcohol, layout = c(1, 2)) #Check depressive symptoms and alcohol intake
densityplot(ergo_imp, ~ BMI + Physical_activity, layout  = c(1, 2)) #Check BMI and physical activity
densityplot(ergo_imp, ~ Work_hours + Hours_home, layout  = c(1, 2)) #Check work and from home hours. 


#Check catagorical variables
propplot <- function(x, formula, facet = "wrap", ...) {
   library(ggplot2)
   cd <- data.frame(mice::complete(x, "long", include = TRUE))
   cd$.imp <- factor(cd$.imp)
   r <- as.data.frame(is.na(x$data))
   impcat <- x$meth != "" & sapply(x$data, is.factor)
   vnames <- names(impcat)[impcat]
   if (missing(formula)) {
      formula <- as.formula(paste(paste(vnames, collapse = "+",
                                        sep = ""), "~1", sep = ""))
   }
   tmsx <- terms(formula[-3], data = x$data)
   xnames <- attr(tmsx, "term.labels")
   xnames <- xnames[xnames %in% vnames]
   if (paste(formula[3]) != "1") {
      wvars <- gsub("[[:space:]]*\\|[[:print:]]*", "", paste(formula)[3])
      # wvars <- all.vars(as.formula(paste("~", wvars)))
      wvars <- attr(terms(as.formula(paste("~", wvars))), "term.labels")
      if (grepl("\\|", formula[3])) {
         svars <- gsub("[[:print:]]*\\|[[:space:]]*", "", paste(formula)[3])
         svars <- all.vars(as.formula(paste("~", svars)))
      } else {
         svars <- ".imp"
      }
   } else {
      wvars <- NULL
      svars <- ".imp"
   }
   for (i in seq_along(xnames)) {
      xvar <- xnames[i]
      select <- cd$.imp != 0 & !r[, xvar]
      cd[select, xvar] <- NA
   }
   for (i in which(!wvars %in% names(cd))) {
      cd[, wvars[i]] <- with(cd, eval(parse(text = wvars[i])))
   }
   meltDF <- reshape2::melt(cd[, c(wvars, svars, xnames)], id.vars = c(wvars, svars))
   meltDF <- meltDF[!is.na(meltDF$value), ]
   wvars <- if (!is.null(wvars)) paste0("`", wvars, "`")
   a <- plyr::ddply(meltDF, c(wvars, svars, "variable", "value"), plyr::summarize,
                    count = length(value))
   b <- plyr::ddply(meltDF, c(wvars, svars, "variable"), plyr::summarize,
                    tot = length(value))
   mdf <- merge(a,b)
   mdf$prop <- mdf$count / mdf$tot
   plotDF <- merge(unique(meltDF), mdf)
   plotDF$value <- factor(plotDF$value,
                          levels = unique(unlist(lapply(x$data[, xnames], levels))),
                          ordered = T)
   p <- ggplot(plotDF, aes(x = value, fill = get(svars), y = prop)) +
      geom_bar(position = "dodge", stat = "identity") +
      theme(legend.position = "bottom", ...) +
      ylab("proportion") +
      scale_fill_manual(name = "",
                        values = c("black",
                                   colorRampPalette(
                                      RColorBrewer::brewer.pal(9, "Blues"))(x$m + 3)[1:x$m + 3])) +
      guides(fill = guide_legend(nrow = 1))
   if (facet == "wrap")
      if (length(xnames) > 1) {
         print(p + facet_wrap(c("variable", wvars), scales = "free"))
      } else {
         if (is.null(wvars)) {
            print(p)
         } else {
            print(p + facet_wrap(wvars, scales = "free"))
         }
      }
   if (facet == "grid")
      if (!is.null(wvars)) {
         print(p + facet_grid(paste(paste(wvars, collapse = "+"), "~ variable"),
                              scales = "free"))
      }
}

propplot(ergo_imp)


#Create variable 
ergo_imp <- mice::complete(ergo_imp, "long", include=TRUE)
ergo_imp$From_home <- ergo_imp$Work_hours + ergo_imp$Hours_home

#Create quartile variables for air pollutants
ergo_imp$pm10quartile <- quantcut(ergo_imp$pm10tselect, q=4, na.rm=TRUE)
ergo_imp$pm2.5quartile <- quantcut(ergo_imp$pm25tselect, q=4, na.rm=TRUE)
ergo_imp$pmabquartile <- quantcut(ergo_imp$absorbancet, q=4, na.rm=TRUE)
ergo_imp$noxquartile <- quantcut(ergo_imp$noxtselecti, q=4, na.rm=TRUE)
ergo_imp$no2quartile <- quantcut(ergo_imp$no2tselecti, q=4, na.rm=TRUE)
ergo_imp$Principle_all_polutantsquartile <- quantcut(ergo_imp$Principle_all_polutants, q=4, na.rm=TRUE)

#Determine pm ab per 0.5 increase, and nox and no2 per 5 ug/m3 increase
ergo_imp$pm10tselect = as.numeric(scale(pm10tselect))
ergo_imp$pm25tselect = as.numeric(scale(pm25tselect))
ergo_imp$absorbancet = as.numeric(scale(absorbancet))
ergo_imp$noxtselecti = as.numeric(scale(noxtselecti))
ergo_imp$no2tselecti = as.numeric(scale(no2tselecti))
ergo_imp$Monthly_income <- as.factor(ergo_imp$Monthly_income)
ergo_imp <- as.mids(ergo_imp)


save(ergo_imp,file="Dataset_air_pollution_and_Dementia.RData")




