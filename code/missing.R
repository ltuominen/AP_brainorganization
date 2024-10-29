library(tidyverse)

df <- read_csv('~/Documents/Research/Reettis/old_new_analyses/master/merged_data.csv')

df <- df %>% filter(df$Group %in% c('FEP', 'CHR') )

r <- df %>% group_by(Group) %>%
  summarize(NonMissingGAF = sum(!is.na(GAF)),
            NonMissingSOFAS = sum(!is.na(SOFAS)),
            NonMissingPositive_symptoms_score = sum(!is.na(Positive_symptoms_score)),
            NonMissingNegative_symptoms_score = sum(!is.na(Negative_symptoms_score)),
            NonMissingTotal_symptoms_score = sum(!is.na(Total_symptoms_score)),
            NonMissingHospital_days = sum(!is.na(Hospital_days)),
            NonMissingTimes_admitted = sum(!is.na(Times_admitted)),
            NonMissingDOI = sum(!is.na(DOI))) %>%
  t()

r

