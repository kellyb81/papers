#title: "Korean Clonal Hematopoiesis project"
#author: "SGK, Genome Opinion Inc."
#date: '2020 6 27'
#output: html_document


install.packages("ztable")
install.packages("sjPlot")
install.packages("table1")

#Load required packages. 

library(dplyr)
library(table1)
library(data.table)
library(broom)
library(ggplot2)
library(ztable)
library(reshape2)
library(sjPlot)

#Load the helper R script to harmonize the the analysis
source("toolbox.R")

#Load the master data file
#### input data ####
covid_msk <- read.csv("input/covid_master_n112_min.csv") %>% as.data.frame(stringsAsFactors=F)


#### data processing ####
# age
agebreaks <- c(0,50,60,70,80,500)
agelabels <- c("0-50","50-60","60-70","70-80","80+")

setDT(covid_msk)[ , agegroups := cut(Age, 
                                     breaks = agebreaks, 
                                     right = FALSE, 
                                     labels = agelabels)]

covid_msk$age_scaled = as.vector(scale(covid_msk$Age))

# smoking (assign smoking status(1) to original smoking group 1 and 2, 0 to group 3)
covid_msk$smoking_msk = ifelse((covid_msk$smoking==1) | (covid_msk$smoking==2), 1, 
                               ifelse(covid_msk$smoking==3, 0, 2)) %>% as.factor()
# Assign ch category 
covid_msk$ch_all = ifelse(covid_msk$chip_all >= 1, 1, 0)
covid_msk$ch_nonpd = ifelse(covid_msk$chip_nonpd >= 1, 1, 0)
covid_msk$ch_pd = ifelse(covid_msk$chip_pd >= 1, 1, 0)
covid_msk$ch_mypd = ifelse(covid_msk$chip_mypd >= 1, 1, 0)
covid_msk$ch_silent = ifelse(covid_msk$chip_silent >= 1, 1, 0)

#assign value according to new MSK rule CH-nonPD_az 
covid_msk <- covid_msk %>% mutate(
  ch_nondriver_az = case_when(
    ch_pd == 1 ~ 0,
    ch_pd == 0 & ch_all == 1 ~ 1,
    ch_all == 0 ~ 0
  )) 

#count muntation number  for each category 
covid_msk$mutnum_all = covid_msk$chip_all
covid_msk$mutnum_pd = covid_msk$chip_pd
covid_msk$mutnum_nonpd = covid_msk$chip_nonpd
covid_msk$mutnum_nonpd_az = ifelse(covid_msk$ch_pd == 1,0,covid_msk$chip_nonpd)

#assign pts to CH_vaf_bin (variable name is kept as same as msk)
covid_msk = covid_msk %>% mutate(
  CH_vaf_bin = case_when(
    chip_all_5 >= 1 & ch_all==1 ~ 2,
    chip_all_5 == 0 & ch_all==1 ~ 1,
    ch_all == 0 ~ 0),
  CH_mutnum_bin = case_when(
    mutnum_all>1 & ch_all==1 ~ 2, 
    mutnum_all==1 & ch_all==1 ~ 1,
    ch_all == 0 ~ 0))


#------------- Printing  pts desc. table based on PD/nonPD status-----Table 1 ----
D = covid_msk %>% 
  mutate(gender_v = case_when (
    gender == 0 ~ "Male",
    gender == 1 ~ "Female"
  )) %>%
  mutate(smoke_bin = case_when(
    smoking_msk == 0 ~ "Never",
    smoking_msk == 1 ~ "Current/Former",
    smoking_msk == 2 ~ "Missing"
  )) %>%
  mutate(htn_verbose = case_when(
    HTN == 0 ~ "No",
    HTN == 1 ~ "Yes"
  )) %>%
  mutate(cad_verbose = case_when(
    CVD == 0 ~ "No",
    CVD == 1 ~ "Yes"
  )) %>%
  mutate (copd_asthma_verbose = case_when(
    COPD_asthma == 0 ~ "No",
    COPD_asthma == 1 ~ "Yes"
  )) %>%
  mutate (diabetes_verbose = case_when(
    DM == 0 ~ "No",
    DM == 1 ~ "Yes"
  )) %>%
  mutate (race_bin = case_when(
    race == 0 ~ "White",
    race == 1 ~ "Non-white"
  )) %>%
  mutate (covid_cat = case_when(
    o2.apply == 1 ~ "Severe COVID",
    o2.apply == 0 ~ "Non-Severe COVID"
  )) %>%
  mutate(CH_all_v = case_when(
    chip_all == 0 ~ "No",
    chip_all >= 1 ~ "Yes"
  )) %>%
  mutate(CH_nondriver_verbose = case_when(
    chip_nonpd == 0 ~ "No",
    chip_nonpd >= 1 ~ "Yes"
  )) %>%
  mutate(ch_pd_verbose = case_when (
    chip_pd == 0 ~ "No",
    chip_pd >= 1 ~ "Yes"
  )) %>%
  mutate(ch_mypd_verbose = case_when (
    chip_mypd == 0 ~ "No",
    chip_mypd >= 1 ~ "Yes"
  ))
  
D <- D %>% mutate(covid_cat = factor(covid_cat, levels = c("Severe COVID","Non-Severe COVID")))
D <- D %>% mutate(smoke_bin = factor(smoke_bin, levels = c("Never","Current/Former", "Missing")))

label(D$agegroups) <- "Age(y)"
label(D$gender_v) <- "Gender"
label(D$smoke_bin) <- "Smoking"
label(D$htn_verbose) <- "Hypertension"
label(D$cad_verbose) <- "Coronary Artery Disease"
label(D$copd_asthma_verbose) <- "COPD/Asthma"
label(D$diabetes_verbose) <- "Diabetes"
label(D$race_bin) <- "Race"

label(D$ch_pd_verbose) <- "CH-PD"
label(D$ch_mypd_verbose) <- "CH-MyPD"
label(D$CH_all_v)<- "Any CH"
label(D$CH_nondriver_verbose)<-"Non-Driver CH"
label(D$ch_pd_verbose)<-"Driver CH"

#--------------Plot KoCH cohort descriptive table (severe vs non-severe)  Table 1----

table1(~ agegroups + gender_v + smoke_bin + htn_verbose + cad_verbose + copd_asthma_verbose + diabetes_verbose + race_bin | covid_cat, data=D, overall=F, output="markdown", export="ktab")


#--------------Plot KoCH cohort descriptive table (CH descriptive stat) Extended Data Table 1 ----
table1(~ CH_all_v + CH_nondriver_verbose + ch_pd_verbose  | covid_cat, data=D, output="markdown", export="ktab")



#--------------Original analysis CH-PD CH-nonPD CH-silent Extended Data Fig 1 ----
# main logistic regression and plot OR  Saves the plot to /output directory: CH mutation type vs severity

glm_sev_vs_nonsev <- do_glm_and_plot_or(data = covid_msk, 
                                        response = "o2.apply", 
                                        terms_for_formula = c("Age + smoking_msk + gender + HTN + CVD + DM + COPD_asthma"), 
                                        terms_to_iterate = c("ch_all", "ch_nondriver_az", "ch_pd","ch_silent"), 
                                        max_OR = 4, 
                                        terms_to_plot = c("ch_all", "ch_nondriver_az", "ch_pd","ch_silent"), 
                                        plot_name = "output/figure1.pdf", width = 12, height = 3)


# heterogeneity test for main 
#D <- covid_msk %>% mutate(nopd_pd=case_when(ch_pd==1 ~ 0,
#                                         ch_nondriver_az==1 ~ 1),
#                       silent_pd=case_when(ch_pd==1 ~ 0,
#                                           ch_silent==1 ~ 1))


#------------- OR plot based on mutation number in CH pts - Extend Data Fig 2 CH-mutation number----
#Saves the plot to /output directory: Mutation number vs severity
covid_msk <- covid_msk %>% mutate(CH_mutnum_bin=as.factor(CH_mutnum_bin))
glm_sev_vs_non_mutnum <-do_glm_and_plot_or(data = covid_msk, 
                                           response = "o2.apply", 
                                           terms_for_formula = c("Age + smoking_msk + gender + HTN + CVD + DM + COPD_asthma"), 
                                           terms_to_iterate = c("CH_mutnum_bin"), 
                                           terms_to_plot = c("CH_mutnum_bin1", "CH_mutnum_bin2"), 
                                           term_factor = "T", 
                                           factors_for_term = c(1,2),
                                           max_OR = 3, 
                                           plot_name = "output/figure3.pdf",width = 12, height = 2)



# heterogeneity test for mutnum 
logit_gene_var = list()
ch_list <-  c("mutnum_all")

for (ch in ch_list) {
  D = covid_msk %>% filter(get(ch)>0)
  logit = glm(
    formula =o2.apply ~ Age + smoking_msk + gender + HTN + CVD + DM + COPD_asthma + get(ch),
    data = D,
    family = "binomial")
  logit_data = logit %>% sjPlot::get_model_data(type="est") %>% cbind(CH = ch)
  logit_gene_var = rbind(logit_gene_var, logit_data)
}

logit_gene_var %>% filter(term=="get(ch)")

#------------- Vaf<=5% Vaf >5% association  Extend Data figure 3 ------
#Saves the plot to /output directory: Mutation vaf vs severity
covid_msk <- covid_msk %>% mutate(CH_vaf_bin=as.factor(CH_vaf_bin))
glm_sev_vs_non_vaf <-do_glm_and_plot_or(data =covid_msk, 
                                        response = "o2.apply", 
                                        terms_for_formula = c("Age + smoking_msk + gender + HTN + CVD + DM + COPD_asthma"), 
                                        terms_to_iterate = c("CH_vaf_bin"), 
                                        terms_to_plot = c("CH_vaf_bin1", "CH_vaf_bin2"), 
                                        term_factor = TRUE, 
                                        factors_for_term = c(1,2),
                                        max_OR = 3, 
                                        plot_name = "output/figure4.pdf",width = 12, height = 2)



#----------- heterogeneity test for vaf
logit_gene_var = list()
ch_list <-  c("VAF_all")

for (ch in ch_list) {
  D = covid_msk %>% filter(get(ch)>0)
  logit = glm(
    formula = o2.apply ~ Age + smoking_msk + gender + HTN + CVD + DM + COPD_asthma + get(ch),
    data = D,
    family = "binomial")
  logit_data = logit %>% sjPlot::get_model_data(type="est") %>% cbind(CH = ch)
  logit_gene_var = rbind(logit_gene_var, logit_data)
}

logit_gene_var %>% filter(term=="get(ch)")

