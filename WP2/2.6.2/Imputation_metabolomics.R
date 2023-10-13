#Paper: Air pollution and metabolomics

cat("\014") 
rm(list=ls())
#Written under R version 4.0.4

#______________________________________________________________________________________
# Open dataset for analysis ----
library(foreign)
library(dplyr)
load("V:\\HomeDir\\051543(T_de_Crom)\\Syntax\\Paper air pollution and metabolomics\\Dataset.RData")

ergo1 <- ergo[,1:20]
rm(ergo)
#______________________________________________________________________________________
# summary statistics 
library(tableone)
summary(ergo1)
#table 1
colSums(is.na(ergo1))
covariates <- c("sex","age","rs_cohort","ses_UNESCO_recoded", "season","smk", "BMI", "alcohol", "Occupation","METh", "CESD", "lipidmede", "hypermed")
exposure <- c("elapse_no2", "elapse_pm25","elapse_bc", "elapse_o3w", "elapse_o3a")
table1 <- CreateTableOne(vars = c(covariates, exposure), data = ergo1)
table_1 <- print(table1, catDigits = 0, contDigits = 1, missing = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
table_1
library(xlsx)
write.xlsx(table_1, "V:/HomeDir/051543(T_de_Crom)/Syntax/Paper air pollution and metabolomics/Results/Table_1_baseline_characteristics.xlsx")

rm(table_1, table1, covariates, exposure)



#______________________________________________________________________________________
#Imputate missing values of covariates ----
library(mice)


# Number and proportion of missing values per variable
cbind("# NA" = sort(colSums(is.na(ergo1))),
      "% NA" = round(sort(colMeans(is.na(ergo1))) * 100, 2)) #multiplied by 100 to obtain percentage and 2 decimale


# set up run of imputation
imp <- mice(ergo1, maxit = 0) #0 iteractions. 
imp$loggedEvents #No collinearity 
meth <- imp$method
predM <- imp$predictorMatrix
pred <- imp$pred


# Setting values of variables you do not want to use as predictor to 0 in the imputation matrix.
predM[, c("ergoid", "rs_cohort", "season", "sampling_date", "elapse_no2", "elapse_pm25", "elapse_bc", "elapse_o3w", "elapse_o3a")] <- 0


# run the imputation 
ergo_imp <- mice(ergo1, method=meth, predictorMatrix=predM, maxit = 5, m=5, seed=2020) 


#check logged events
#removed predictors because of multicolinearity is noted by loggedEvents
ergo_imp$loggedEvents


# check the number of missings in the first imputed dataset, to see if nothing was skipped.
cbind("# NA" = sort(colSums(is.na(complete(ergo_imp, 1)))),
      "% NA" = round(sort(colMeans(is.na(complete(ergo_imp, 1)))) * 100, 2))


#Check original database vs imputed database
#Ceck continous variables
densityplot(ergo_imp, ~ alcohol + METh, layout = c(1, 2)) 
densityplot(ergo_imp, ~ BMI + CESD, layout = c(1, 2)) 

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


rm(ergo1, imp, pred, predM, meth, propplot)
#______________________________________________________________________________________
#Create average dataset ----
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
  mutate_at(c("rs_cohort", "sex", "ses_UNESCO_recoded", "smk", "Occupation", "lipidmede", "hypermed", "season"), as.numeric)

#- then calculate average of all variables of all imputed datasets
ergo_imp_ave <- ergo_imp_vol %>% 
  group_by(.id2) %>% #indicating timepoint per ergoid per imputation
  summarise_all(mean, na.rm=T) %>% 
  mutate_at(c("rs_cohort", "sex", "ses_UNESCO_recoded", "smk", "Occupation", "lipidmede", "hypermed", "season"), funs(round)) %>% 
  mutate_at(c("rs_cohort", "sex", "ses_UNESCO_recoded", "smk", "Occupation", "lipidmede", "hypermed", "season"), funs(factor)) 

#Check missings
cbind("# NA" = sort(colSums(is.na(ergo_imp_ave))),
      "% NA" = round(sort(colMeans(is.na(ergo_imp_ave))) * 100, 2)) #multiplied by 100 to obtain percentage and 2 decimale

ergo_imp_ave$.id <- NULL
ergo_imp_ave$.id2 <- NULL
ergo_imp_ave$.imp <- NULL
ergo_imp_ave$column_label <- NULL

rm(ergo_imp, ergo_imp_vol, imp_list_vol, imp_list_vol_2, merge_vol)
#______________________________________________________________________________________
#Add metabolomics ----
load("RS1_4_LDD_imputed_nlog_transformed_April13.RData")
mb_e4 <- imputeddata
load("RSIII_2_LDD_imputed_nlog_transformed_2Jun2021.RData")
mb_e5 <- imputeddata
rm(imputeddata)

# find overlapping metabolites between e4 and e5
mb_overlap <- Reduce(intersect, list(colnames(mb_e4), colnames(mb_e5)))

# subset only metabolites in both rounds
mb_e4 <- subset(mb_e4, select = mb_overlap)
mb_e5 <- subset(mb_e5, select = mb_overlap)
metabolomics <- rbind(mb_e4, mb_e5)
rm(mb_e4, mb_e5, mb_overlap)
#______________________________________________________________________________________
#Create final dataset ----
ergo_final <- merge(ergo_imp_ave, metabolomics, by.x = "ergoid", by.y = 'row.names', all = TRUE, all.y = FALSE)

save(ergo_final ,file="V:/HomeDir/051543(T_de_Crom)/Syntax/Paper air pollution and metabolomics/Dataset_final.RData")

