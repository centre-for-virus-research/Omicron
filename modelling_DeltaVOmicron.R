
# Modelling script for VE analysis of Delta v Omicron for those tested between 6th December 2021 and 26th December 2021:

load("data_rect_all_2022_02_11_v2.RData") # Loads the data frame "data.rect.all" into the workspace.

# The data frame named "data.rect.all" consists of the following columns:

# * SafeHavenID : anonymised ID number, one per individual in dataset (i.e. one per CHI number).
# * date_1st_dose : date (yyyy-mm-dd) of 1st SARS-CoV-2 vaccine dose.
# * date_2nd_dose : date (yyyy-mm-dd) of 2nd SARS-CoV-2 vaccine dose.
# * date_3rd_dose : date (yyyy-mm-dd) of 3rd SARS-CoV-2 vaccine dose.
# * date_4th_dose : date (yyyy-mm-dd) of 4th SARS-CoV-2 vaccine dose.
# * date_5th_dose : date (yyyy-mm-dd) of 5th SARS-CoV-2 vaccine dose.
# * product_name_1st : product name of 1st SARS-CoV-2 vaccine dose.
# * product_name_2nd : product name of 2nd SARS-CoV-2 vaccine dose.
# * product_name_3rd : product name of 3rd SARS-CoV-2 vaccine dose.
# * product_name_4th : product name of 4th SARS-CoV-2 vaccine dose.
# * product_name_5th : product name of 5th SARS-CoV-2 vaccine dose.
# * booster_1st : indicator if 1st SARS-CoV-2 vaccine dose is a booster (1 if so, 0 if not).
# * booster_2nd : indicator if 2nd SARS-CoV-2 vaccine dose is a booster (1 if so, 0 if not).
# * booster_3rd : indicator if 3rd SARS-CoV-2 vaccine dose is a booster (1 if so, 0 if not).
# * booster_4th : indicator if 4th SARS-CoV-2 vaccine dose is a booster (1 if so, 0 if not).
# * booster_5th : indicator if 5th SARS-CoV-2 vaccine dose is a booster (1 if so, 0 if not).
# * DOB : date of birth (yyyy-mm-dd).
# * age_31Oct21 : age on 31st October 2021.
# * sex : sex (F: female; M: male).
# * SIMD_2016_VIGINTILE : vigintile of individual's postcode Scottish Index of Multiple Deprivation (SIMD) 2016 rank.
# * SIMD_2016_RANK : rank of individual's postcode Scottish Index of Multiple Deprivation (SIMD) 2016 status.
# * date_1st_pos_PCR : date (yyyy-mm-dd) of 1st positive PCR test for SARS-CoV-2.
# * date_1st_pos_PCR_after90 : date (yyyy-mm-dd) of 2nd positive PCR test for SARS-CoV-2 (at least 90 days after previous most recent positive test).
# * date_2nd_pos_PCR_after90 : date (yyyy-mm-dd) of 3rd positive PCR test for SARS-CoV-2 (at least 90 days after previous most recent positive test).
# * asp1 : pillar 1 ASP (allele-specific PCR) status for 1st positive PCR test for SARS-CoV-2.
# * asp2 : pillar 2 ASP status for 1st positive PCR test for SARS-CoV-2.
# * sg : SGTF (S-gene test failure) status for 1st positive PCR test for SARS-CoV-2.
# * M1 : indicator if Omicron (BA.1) for 1st positive PCR test for SARS-CoV-2 for small sample of data points.
# * lineage : Pango lineage for 1st positive PCR test for SARS-CoV-2.
# * asp1_1st_after90 # As for 1st positive PCR test - see above. 
# * asp2_1st_after90
# * sg_1st_after90 
# * M1_1st_after90 
# * lineage_1st_after90 
# * asp1_2nd_after90 
# * asp2_2nd_after90 
# * sg_2nd_after90 
# * M1_2nd_after90 
# * lineage_2nd_after90 
# * BA.1.any : Indicator if BA.1 or BA.1.1 (1 if so, 0 otherwise).
# * delta : Indicator if Delta (1 if so, 0 otherwise).
# * BA.1.1 : Indicator if BA.1.1 (1 if so, 0 otherwise).
# * BA.1.only : Indicator if BA.1 and not BA.1.1 (1 if so, 0 otherwise).
# * BA.1.any_1st_after90 # As for 1st positive PCR test - see above. 
# * delta_1st_after90 
# * BA.1.1_1st_after90 
# * BA.1.only_1st_after90 
# * BA.1.any_2nd_after90 
# * delta_2nd_after90 
# * BA.1.1_2nd_after90 
# * BA.1.only_2nd_after90 
# * omicron : Indicator if BA.1 or BA.1.1, from any of Pango lineage BA.1 or BA.1.1, asp1=VOC-21NOV-01, asp2=VOC-21NOV-01, M1=1 or sg=True S Gene Dropout (1 if so, 0 otherwise).
# * omicron_1st_after90 # As for 1st positive PCR test - see above. 
# * omicron_2nd_after90 
# * imm.supp.drugs : Indicator if individual is on immunosuppressent drugs (1 if so, 0 otherwise).
# * shielding.group : Indicator if individual is on shielding list due to immunosuppression (1 if so, 0 otherwise).

# Create study period:

period.study<-c(as.Date("2021-12-06"),as.Date("2021-12-26"))

# Create subset of only those tested during the study period (since they definitely still live in GGC during that time)

# (Need testing data)
SCI_STORE<-read.csv("Z:/20220210/SCI_Store_DOVE_WK4.csv")
SCI_STORE$SAMPLEDATE<-as.Date(SCI_STORE$SAMPLEDATE)

IDs.tested<-sort(unique(SCI_STORE$SafeHavenID[which(SCI_STORE$SAMPLEDATE>=period.study[1]&
                                                      SCI_STORE$SAMPLEDATE<=period.study[2])]))
length(IDs.tested)
# [1] 126418

# Before beginning, remove those with 1st doses before 8th December 2020 (as these are likely to be trial participants):
data.rect.all.2<-data.rect.all[-which(data.rect.all$date_1st_dose<=as.Date("2020-12-08")),]

data.rect.tested<-data.rect.all.2[data.rect.all.2$SafeHavenID%in%IDs.tested,]
dim(data.rect.tested)
# [1] 126279     63

# Find dates of first positive PCR during study period:

data.rect.tested$date_1st_pos_PCR_study<-ifelse(data.rect.tested$date_1st_pos_PCR>=period.study[1]&
                                                  data.rect.tested$date_1st_pos_PCR<=period.study[2],
                                                data.rect.tested$date_1st_pos_PCR,
                                                ifelse(data.rect.tested$date_1st_pos_PCR_after90>=period.study[1]&
                                                         data.rect.tested$date_1st_pos_PCR_after90<=period.study[2],
                                                       data.rect.tested$date_1st_pos_PCR_after90,
                                                       ifelse(data.rect.tested$date_2nd_pos_PCR_after90>=period.study[1]&
                                                                data.rect.tested$date_2nd_pos_PCR_after90<=period.study[2],
                                                              data.rect.tested$date_2nd_pos_PCR_after90,NA)))
data.rect.tested$date_1st_pos_PCR_study<-as.Date(data.rect.tested$date_1st_pos_PCR_study,origin="1970-01-01")

data.rect.tested$infected<-!is.na(data.rect.tested$date_1st_pos_PCR_study)

# Note: For this study, samples with "Positive S Gene" status are assumed to be Delta unless other data source classes them as Omicron.
data.rect.tested$variant.main<-ifelse(data.rect.tested$date_1st_pos_PCR>=period.study[1]&
                                        data.rect.tested$date_1st_pos_PCR<=period.study[2],
                                      ifelse(data.rect.tested$delta==1|(data.rect.tested$sg=="Positive S Gene"&!is.na(data.rect.tested$sg)&data.rect.tested$omicron==0),"Delta",
                                             ifelse(data.rect.tested$omicron==1,"Omicron",
                                                    "other_or_unknown")),
                                      ifelse(data.rect.tested$date_1st_pos_PCR_after90>=period.study[1]&
                                               data.rect.tested$date_1st_pos_PCR_after90<=period.study[2],
                                             ifelse(data.rect.tested$delta_1st_after90==1|(data.rect.tested$sg_1st_after90=="Positive S Gene"&!is.na(data.rect.tested$sg_1st_after90)&data.rect.tested$omicron_1st_after90==0),"Delta",
                                                    ifelse(data.rect.tested$omicron_1st_after90==1,"Omicron",
                                                           "other_or_unknown")),
                                             ifelse(data.rect.tested$date_2nd_pos_PCR_after90>=period.study[1]&
                                                      data.rect.tested$date_2nd_pos_PCR_after90<=period.study[2],
                                                    ifelse(data.rect.tested$delta_2nd_after90==1|(data.rect.tested$sg_2nd_after90=="Positive S Gene"&!is.na(data.rect.tested$sg_2nd_after90)&data.rect.tested$omicron_2nd_after90==0),"Delta",
                                                           ifelse(data.rect.tested$omicron_2nd_after90==1,"Omicron",
                                                                  "other_or_unknown")),NA)))

# New column with variant defined just by sequencing info:
data.rect.tested$variant.main.seqOnly<-ifelse(data.rect.tested$date_1st_pos_PCR>=period.study[1]&
                                                data.rect.tested$date_1st_pos_PCR<=period.study[2],
                                              ifelse(grepl('^AY|^B\\.1\\.617\\.2$',data.rect.tested$lineage),"Delta",
                                                     ifelse(data.rect.tested$lineage=="BA.1"|
                                                              data.rect.tested$lineage=="BA.1.1","Omicron",
                                                            "other_or_unknown")),
                                              ifelse(data.rect.tested$date_1st_pos_PCR_after90>=period.study[1]&
                                                       data.rect.tested$date_1st_pos_PCR_after90<=period.study[2],
                                                     ifelse(grepl('^AY|^B\\.1\\.617\\.2$',data.rect.tested$lineage_1st_after90),"Delta",
                                                            ifelse(data.rect.tested$lineage_1st_after90=="BA.1"|
                                                                     data.rect.tested$lineage_1st_after90=="BA.1.1","Omicron",
                                                                   "other_or_unknown")),
                                                     ifelse(data.rect.tested$date_2nd_pos_PCR_after90>=period.study[1]&
                                                              data.rect.tested$date_2nd_pos_PCR_after90<=period.study[2],
                                                            ifelse(grepl('^AY|^B\\.1\\.617\\.2$',data.rect.tested$lineage_2nd_after90),"Delta",
                                                                   ifelse(data.rect.tested$lineage_2nd_after90=="BA.1"|
                                                                            data.rect.tested$lineage_2nd_after90=="BA.1.1","Omicron",
                                                                          "other_or_unknown")),NA)))
data.rect.tested$variant.main.seqOnly[data.rect.tested$infected&is.na(data.rect.tested$variant.main.seqOnly)]<-"other_or_unknown"

data.rect.tested$previous.infection<-ifelse(data.rect.tested$date_1st_pos_PCR<=period.study[1]-90&!is.na(data.rect.tested$date_1st_pos_PCR)|
                                              data.rect.tested$date_1st_pos_PCR_after90<=period.study[1]-90&!is.na(data.rect.tested$date_1st_pos_PCR_after90)|
                                              data.rect.tested$date_2nd_pos_PCR_after90<=period.study[1]-90&!is.na(data.rect.tested$date_2nd_pos_PCR_after90),
                                            TRUE,FALSE)

# Add vaccination status 14 days before study period:

data.rect.tested$number.vacc<-ifelse(data.rect.tested$date_5th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_5th_dose),5,
                                     ifelse(data.rect.tested$date_4th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_4th_dose),4,
                                            ifelse(data.rect.tested$date_3rd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_3rd_dose),3,
                                                   ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_2nd_dose),2,
                                                          ifelse(data.rect.tested$date_1st_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_1st_dose),1,
                                                                 0)))))

data.rect.tested$status.vacc<-ifelse(data.rect.tested$date_5th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_5th_dose),"5 doses",
                                     ifelse(data.rect.tested$date_4th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_4th_dose),"4 doses",
                                            ifelse(data.rect.tested$date_3rd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_3rd_dose),paste0("3rd dose ",data.rect.tested$product_name_3rd),
                                                   ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_2nd_dose),
                                                          ifelse(data.rect.tested$product_name_1st==data.rect.tested$product_name_2nd,
                                                                 paste0("2 doses ",data.rect.tested$product_name_2nd),
                                                                 "2 doses (different)"),
                                                          ifelse(data.rect.tested$date_1st_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_1st_dose),"1 dose only",
                                                                 "unvaccinated")))))
table(data.rect.tested$status.vacc)

data.rect.tested$status.vacc.previnf.int<-ifelse(data.rect.tested$date_5th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_5th_dose),"5 doses",
                                     ifelse(data.rect.tested$date_4th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_4th_dose),"4 doses",
                                            ifelse(data.rect.tested$date_3rd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_3rd_dose),paste0("3rd dose ",data.rect.tested$product_name_3rd),
                                                   ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_2nd_dose),
                                                          ifelse(data.rect.tested$product_name_1st==data.rect.tested$product_name_2nd,
                                                                 paste0("2 doses ",data.rect.tested$product_name_2nd),
                                                                 "2 doses (different)"),
                                                          ifelse(data.rect.tested$date_1st_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_1st_dose),"1 dose only",
                                                                 "unvaccinated")))))
table(data.rect.tested$status.vacc.previnf.int)

# Interaction with previous infection added here:
data.rect.tested$status.vacc.previnf.int<-ifelse(data.rect.tested$date_5th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_5th_dose),"5 doses",
                                                 ifelse(data.rect.tested$date_4th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_4th_dose),"4 doses",
                                                        ifelse(data.rect.tested$date_3rd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_3rd_dose),
                                                               ifelse(data.rect.tested$previous.infection,
                                                                      paste0("3rd dose (prev inf) ",data.rect.tested$product_name_3rd),
                                                                      paste0("3rd dose (no prev inf) ",data.rect.tested$product_name_3rd)),
                                                               ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_2nd_dose),
                                                                      ifelse(data.rect.tested$product_name_1st==data.rect.tested$product_name_2nd,
                                                                             ifelse(data.rect.tested$previous.infection,
                                                                                    paste0("2 doses (prev inf) ",data.rect.tested$product_name_2nd),
                                                                                    paste0("2 doses (no prev inf) ",data.rect.tested$product_name_2nd)),
                                                                             "2 doses (different)"),
                                                                      ifelse(data.rect.tested$date_1st_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_1st_dose),
                                                                             ifelse(data.rect.tested$previous.infection,"1 dose only (prev inf)","1 dose only (no prev inf)"),
                                                                             ifelse(data.rect.tested$previous.infection,"unvaccinated (prev inf)","unvaccinated (no prev inf)"))))))
table(data.rect.tested$status.vacc.previnf.int)

data.rect.tested$time.since.most.recent.vacc<-ifelse(data.rect.tested$date_5th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_5th_dose),(period.study[1]-14)-data.rect.tested$date_5th_dose,
                                                     ifelse(data.rect.tested$date_4th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_4th_dose),(period.study[1]-14)-data.rect.tested$date_4th_dose,
                                                            ifelse(data.rect.tested$date_3rd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_3rd_dose),(period.study[1]-14)-data.rect.tested$date_3rd_dose,
                                                                   ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_2nd_dose),(period.study[1]-14)-data.rect.tested$date_2nd_dose,
                                                                          ifelse(data.rect.tested$date_1st_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_1st_dose),(period.study[1]-14)-data.rect.tested$date_1st_dose,
                                                                                 NA)))))

# Adding time since vaccination for 2nd doses:
data.rect.tested$status.vacc.previnf.time.int<-ifelse(data.rect.tested$date_5th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_5th_dose),"5 doses",
                                                        ifelse(data.rect.tested$date_4th_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_4th_dose),"4 doses",
                                                               ifelse(data.rect.tested$date_3rd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_3rd_dose),
                                                                      ifelse(data.rect.tested$previous.infection,
                                                                             paste0("3rd dose (prev inf) ",data.rect.tested$product_name_3rd),
                                                                             paste0("3rd dose (no prev inf) ",data.rect.tested$product_name_3rd)),
                                                                      ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_2nd_dose),
                                                                             ifelse(data.rect.tested$product_name_1st==data.rect.tested$product_name_2nd,
                                                                                    ifelse(data.rect.tested$previous.infection,
                                                                                           ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-(12*7))&!is.na(data.rect.tested$date_2nd_dose),
                                                                                                  paste0("2 doses (prev inf, >=12 weeks since dose 2) ",data.rect.tested$product_name_2nd),
                                                                                                  paste0("2 doses (prev inf, <12 weeks since dose 2) ",data.rect.tested$product_name_2nd)),
                                                                                           ifelse(data.rect.tested$date_2nd_dose<=(period.study[1]-(12*7))&!is.na(data.rect.tested$date_2nd_dose),
                                                                                                  paste0("2 doses (no prev inf, >=12 weeks since dose 2) ",data.rect.tested$product_name_2nd),
                                                                                                  paste0("2 doses (no prev inf, <12 weeks since dose 2) ",data.rect.tested$product_name_2nd))),
                                                                                    "2 doses (different)"),
                                                                             ifelse(data.rect.tested$date_1st_dose<=(period.study[1]-14)&!is.na(data.rect.tested$date_1st_dose),
                                                                                    ifelse(data.rect.tested$previous.infection,"1 dose only (prev inf)","1 dose only (no prev inf)"),
                                                                                    ifelse(data.rect.tested$previous.infection,"unvaccinated (prev inf)","unvaccinated (no prev inf)"))))))
table(data.rect.tested$status.vacc.previnf.time.int)

# Remove those with vaccine status different from those expected (to avoid including trial participants):

# (Also, need to remove those with only 2nd dose listed, since we'd need to know if they had 1st and 2nd doses of same product)

data.rect.tested.2<-data.rect.tested[which((data.rect.tested$status.vacc.previnf.int=="1 dose only (no prev inf)"|
                                              data.rect.tested$status.vacc.previnf.int=="2 doses (no prev inf) Covid-19 mRNA Vaccine Moderna"|
                                              data.rect.tested$status.vacc.previnf.int=="2 doses (no prev inf) Covid-19 mRNA Vaccine Pfizer"|
                                              data.rect.tested$status.vacc.previnf.int=="2 doses (no prev inf) Covid-19 Vaccine AstraZeneca"|
                                              data.rect.tested$status.vacc.previnf.int=="3rd dose (no prev inf) Covid-19 mRNA Vaccine Moderna"|
                                              data.rect.tested$status.vacc.previnf.int=="3rd dose (no prev inf) Covid-19 mRNA Vaccine Pfizer"|
                                              data.rect.tested$status.vacc.previnf.int=="unvaccinated (no prev inf)"|
                                              data.rect.tested$status.vacc.previnf.int=="1 dose only (prev inf)"|
                                              data.rect.tested$status.vacc.previnf.int=="2 doses (prev inf) Covid-19 mRNA Vaccine Moderna"|
                                              data.rect.tested$status.vacc.previnf.int=="2 doses (prev inf) Covid-19 mRNA Vaccine Pfizer"|
                                              data.rect.tested$status.vacc.previnf.int=="2 doses (prev inf) Covid-19 Vaccine AstraZeneca"|
                                              data.rect.tested$status.vacc.previnf.int=="3rd dose (prev inf) Covid-19 mRNA Vaccine Moderna"|
                                              data.rect.tested$status.vacc.previnf.int=="3rd dose (prev inf) Covid-19 mRNA Vaccine Pfizer"|
                                              data.rect.tested$status.vacc.previnf.int=="unvaccinated (prev inf)")&
                                             (data.rect.tested$product_name_1st==data.rect.tested$product_name_2nd|
                                                is.na(data.rect.tested$product_name_2nd))),]
dim(data.rect.tested.2)
dim(data.rect.tested)
# 124852 compared to 126279 rows (1427 removed).

# IMPORTANT! Remove those testing positive less than 90 days before the study period starts:

data.rect.tested.2<-data.rect.tested.2[-which(data.rect.tested.2$date_1st_pos_PCR>=(period.study[1]-90)&
                                         data.rect.tested.2$date_1st_pos_PCR<=period.study[1]),]
data.rect.tested.2<-data.rect.tested.2[-which(data.rect.tested.2$date_1st_pos_PCR_after90>=(period.study[1]-90)&
                                                data.rect.tested.2$date_1st_pos_PCR_after90<=period.study[1]),]
data.rect.tested.2<-data.rect.tested.2[-which(data.rect.tested.2$date_2nd_pos_PCR_after90>=(period.study[1]-90)&
                                                data.rect.tested.2$date_2nd_pos_PCR_after90<=period.study[1]),]

## Summary stats: #####################################################

table(data.rect.tested.2$variant.main)
# Delta          Omicron other_or_unknown 
#  5963            18432             2150 

table(data.rect.tested.2$status.vacc.previnf.int,data.rect.tested.2$variant.main)

table(data.rect.tested.2$status.vacc.previnf.time.int,data.rect.tested.2$variant.main)

#

table(data.rect.tested.2$status.vacc.previnf.int,data.rect.tested.2$infected)
prop.table(table(data.rect.tested.2$status.vacc.previnf.int,data.rect.tested.2$infected),1)

tapply(data.rect.tested.2$age_31Oct21,data.rect.tested.2$infected,summary)
boxplot(data.rect.tested.2$age_31Oct21~data.rect.tested.2$infected)

tapply(data.rect.tested.2$age_31Oct21,data.rect.tested.2$status.vacc.previnf.int,summary)
boxplot(data.rect.tested.2$age_31Oct21~data.rect.tested.2$status.vacc.previnf.int,las=2,cex.axis=0.5)

# plot(data.rect.tested.2$age_31Oct21,data.rect.tested.2$time.since.most.recent.vacc)

tapply(data.rect.tested.2$age_31Oct21,data.rect.tested.2$infected,summary)
boxplot(data.rect.tested.2$age_31Oct21~data.rect.tested.2$infected)

#

tapply(data.rect.tested.2$time.since.most.recent.vacc,data.rect.tested.2$infected,summary)
boxplot(data.rect.tested.2$time.since.most.recent.vacc~data.rect.tested.2$infected)

boxplot(data.rect.tested.2$time.since.most.recent.vacc~data.rect.tested.2$infected+data.rect.tested.2$number.vacc)

########################

# Prepare datasets for Delta and Omicron, for modelling:

# Set base levels of factors:
data.rect.tested.2$status.vacc.previnf.int<-as.factor(data.rect.tested.2$status.vacc.previnf.int)
levels(data.rect.tested.2$status.vacc.previnf.int)
data.rect.tested.2$status.vacc.previnf.int<-factor(data.rect.tested.2$status.vacc.previnf.int,
                                                   levels=levels(data.rect.tested.2$status.vacc.previnf.int)[
                                                     c(13:14,1:2,5,8,4,7,3,6,10,12,9,11)
                                                     ])
data.rect.tested.2$status.vacc.previnf.time.int<-as.factor(data.rect.tested.2$status.vacc.previnf.time.int)
levels(data.rect.tested.2$status.vacc.previnf.time.int)
data.rect.tested.2$status.vacc.previnf.time.int<-factor(data.rect.tested.2$status.vacc.previnf.time.int,
                                                        levels=levels(data.rect.tested.2$status.vacc.previnf.time.int)[
                                                          c(19,20,1,2,5,4,3,8,7,6,11,10,9,14,13,12,16,15,18,17)
                                                          ])

data.Delta<-data.rect.tested.2[which(data.rect.tested.2$variant.main=="Delta"|is.na(data.rect.tested.2$variant.main)),]
data.Omicron<-data.rect.tested.2[which(data.rect.tested.2$variant.main=="Omicron"|is.na(data.rect.tested.2$variant.main)),]

summary(data.Delta$infected)
summary(data.Omicron$infected)

library(mgcv)

mod1.Delta<-gam(infected~status.vacc.previnf.int+
                  sex+s(age_31Oct21)+s(SIMD_2016_VIGINTILE),
                data=data.Delta,family = binomial())
summary.mod1.Delta<-summary(mod1.Delta);summary.mod1.Delta
plot(mod1.Delta,pages=1)

mod1.Omicron<-gam(infected~status.vacc.previnf.int+
                    sex+s(age_31Oct21)+s(SIMD_2016_VIGINTILE),
                data=data.Omicron,family = binomial())
summary.mod1.Omicron<-summary(mod1.Omicron);summary.mod1.Omicron
plot(mod1.Omicron,pages=1)

# Plots and output:

names.vec<-c("Unvaccinated (prev. inf.)",paste0(rep(c("One dose only","2 x ChAdOx1","2 x BNT162b2","2 x mRNA-1273",
             "3rd dose BNT162b2","3rd dose mRNA-1273"),each=2)," ",rep(c("(no prev. inf.)","(prev. inf)"),6)))

ggdata.mod1<-cbind.data.frame(names=rep(names.vec,2),
                              variant=rep(c("Delta","Omicron"),each=length(names.vec)),
                              pred=c(100*(1-exp(summary.mod1.Delta$p.coeff))[c(2:14)],
                                     100*(1-exp(summary.mod1.Omicron$p.coeff))[c(2:14)]),
                              lwrbnd=c(100*(1-exp(summary.mod1.Delta$p.coeff+1.96*summary.mod1.Delta$se[1:length(summary.mod1.Delta$p.coeff)]))[c(2:14)],
                                       100*(1-exp(summary.mod1.Omicron$p.coeff+1.96*summary.mod1.Omicron$se[1:length(summary.mod1.Omicron$p.coeff)]))[c(2:14)]),
                              uprbnd=c(100*(1-exp(summary.mod1.Delta$p.coeff-1.96*summary.mod1.Delta$se[1:length(summary.mod1.Delta$p.coeff)]))[c(2:14)],
                                       100*(1-exp(summary.mod1.Omicron$p.coeff-1.96*summary.mod1.Omicron$se[1:length(summary.mod1.Omicron$p.coeff)]))[c(2:14)]))
ggdata.mod1$names<-as.factor(ggdata.mod1$names)
ggdata.mod1$names<-factor(ggdata.mod1$names,levels = unique(ggdata.mod1$names)[c(seq(1,13,by=2),seq(2,12,by=2))])
ggdata.mod1$variant<-as.factor(ggdata.mod1$variant)
ggdata.mod1$variant<-factor(ggdata.mod1$variant,levels = unique(ggdata.mod1$variant))

library(ggplot2)
pal<-c("#12B712","#918CF4")
hjust1<-rep(-0.3,nrow(ggdata.mod1))

ggdata2.mod1<-ggdata.mod1
ggdata2.mod1$names<-rep(c("Unvaccinated",rep(c("One dose only",
                                           "2 x ChAdOx1","2 x BNT162b2","2 x mRNA-1273",
                                           "3rd dose BNT162b2","3rd dose mRNA-1273"),each=2)),2)
ggdata2.mod1$names<-factor(ggdata2.mod1$names,levels=unique(ggdata2.mod1$names))
ggdata2.mod1$prev.inf<-rep(c("Previously confirmed infected",
                             rep(c("No previous confirmed infection","Previously confirmed infected"),6)),2)


gg2<-ggplot(ggdata2.mod1,aes(x=names,y=pred,group=variant))+
  theme_classic()+scale_color_manual(values=pal)+
  geom_errorbar(aes(ymin=lwrbnd,ymax=uprbnd,width=0.4,col=variant),
                alpha=0.5,position=position_dodge(width=0.4),size=1)+
  geom_point(size=4,position=position_dodge(width=0.4),aes(shape=variant,col=variant))+
  geom_text(aes(label=round(pred,2)),hjust=hjust1,position = position_dodge(width=0.4))+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=17,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(),
        legend.title.align = 0.5,
        strip.background = element_rect(color="white"),
        strip.placement = "inside",
        strip.text = element_text(size=15,face="bold"),
        panel.spacing = unit(2, "lines"))+
  xlab("Vaccination status")+ylab("Vaccine effectiveness (%)")+
  facet_wrap(~prev.inf)+
  # coord_cartesian(ylim=c(0,100),xlim=c(0.7,7.5))+
  geom_vline(aes(xintercept=1.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=2.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=5.5),alpha=0.8,linetype="dotted")
gg2

## After removing immunosuppressed individuals: #########

data.Delta.nonImm<-data.Delta[data.Delta$imm.supp.drugs==0&data.Delta$shielding.group==0,]
dim(data.Delta)
dim(data.Delta.nonImm)

data.Omicron.nonImm<-data.Omicron[data.Omicron$imm.supp.drugs==0&data.Omicron$shielding.group==0,]
dim(data.Omicron)
dim(data.Omicron.nonImm)

mod1.Delta.nonImm<-gam(infected~status.vacc.previnf.int+
                  sex+s(age_31Oct21)+s(SIMD_2016_VIGINTILE),
                data=data.Delta.nonImm,family = binomial())
summary.mod1.Delta.nonImm<-summary(mod1.Delta.nonImm);summary.mod1.Delta.nonImm
plot(mod1.Delta.nonImm,pages=1)

mod1.Omicron.nonImm<-gam(infected~status.vacc.previnf.int+
                    sex+s(age_31Oct21)+s(SIMD_2016_VIGINTILE),
                  data=data.Omicron.nonImm,family = binomial())
summary.mod1.Omicron.nonImm<-summary(mod1.Omicron.nonImm);summary.mod1.Omicron.nonImm
plot(mod1.Omicron.nonImm,pages=1)

# Plots and output:

names.vec<-c("Unvaccinated (prev. inf.)",paste0(rep(c("One dose only","2 x ChAdOx1","2 x BNT162b2","2 x mRNA-1273",
                                                      "3rd dose BNT162b2","3rd dose mRNA-1273"),each=2)," ",rep(c("(no prev. inf.)","(prev. inf)"),6)))

ggdata.mod1.nonImm<-cbind.data.frame(names=rep(names.vec,2),
                              variant=rep(c("Delta","Omicron"),each=length(names.vec)),
                              pred=c(100*(1-exp(summary.mod1.Delta.nonImm$p.coeff))[c(2:14)],
                                     100*(1-exp(summary.mod1.Omicron.nonImm$p.coeff))[c(2:14)]),
                              lwrbnd=c(100*(1-exp(summary.mod1.Delta.nonImm$p.coeff+1.96*summary.mod1.Delta.nonImm$se[1:length(summary.mod1.Delta.nonImm$p.coeff)]))[c(2:14)],
                                       100*(1-exp(summary.mod1.Omicron.nonImm$p.coeff+1.96*summary.mod1.Omicron.nonImm$se[1:length(summary.mod1.Omicron.nonImm$p.coeff)]))[c(2:14)]),
                              uprbnd=c(100*(1-exp(summary.mod1.Delta.nonImm$p.coeff-1.96*summary.mod1.Delta.nonImm$se[1:length(summary.mod1.Delta.nonImm$p.coeff)]))[c(2:14)],
                                       100*(1-exp(summary.mod1.Omicron.nonImm$p.coeff-1.96*summary.mod1.Omicron.nonImm$se[1:length(summary.mod1.Omicron.nonImm$p.coeff)]))[c(2:14)]))
ggdata.mod1.nonImm$names<-as.factor(ggdata.mod1.nonImm$names)
ggdata.mod1.nonImm$names<-factor(ggdata.mod1.nonImm$names,levels = unique(ggdata.mod1.nonImm$names)[c(seq(1,13,by=2),seq(2,12,by=2))])
ggdata.mod1.nonImm$variant<-as.factor(ggdata.mod1.nonImm$variant)
ggdata.mod1.nonImm$variant<-factor(ggdata.mod1.nonImm$variant,levels = unique(ggdata.mod1.nonImm$variant))

library(ggplot2)
pal<-c("#12B712","#918CF4")
hjust1<-rep(-0.3,nrow(ggdata.mod1.nonImm))

ggdata2.mod1.nonImm<-ggdata.mod1.nonImm
ggdata2.mod1.nonImm$names<-rep(c("Unvaccinated",rep(c("One dose only",
                                               "2 x ChAdOx1","2 x BNT162b2","2 x mRNA-1273",
                                               "3rd dose BNT162b2","3rd dose mRNA-1273"),each=2)),2)
ggdata2.mod1.nonImm$names<-factor(ggdata2.mod1.nonImm$names,levels=unique(ggdata2.mod1.nonImm$names))
ggdata2.mod1.nonImm$prev.inf<-rep(c("Previously confirmed infected",
                             rep(c("No previous confirmed infection","Previously confirmed infected"),6)),2)

# pred as character, to avoid .0 being removed from text on plot:
ggdata2.mod1.nonImm$pred.char<-as.character(round(ggdata2.mod1.nonImm$pred,1))
ggdata2.mod1.nonImm$pred.char[!grepl("\\.",ggdata2.mod1.nonImm$pred.char)]<-
  paste0(ggdata2.mod1.nonImm$pred.char[!grepl("\\.",ggdata2.mod1.nonImm$pred.char)],".0")

# Add rows for unvaccinated, no previous infection, to create points at zero:
ggdata2.mod1.nonImm<-rbind.data.frame(ggdata2.mod1.nonImm,
                                      data.frame(names=rep("Unvaccinated",2),
                                                 variant=c("Delta","Omicron"),
                                                 pred=rep(0,2),
                                                 lwrbnd=rep(NA,2),
                                                 uprbnd=rep(NA,2),
                                                 prev.inf=rep("No previous confirmed infection",2),
                                                 pred.char=rep("",2)))

gg2.nonImm<-ggplot(ggdata2.mod1.nonImm,aes(x=names,y=pred,group=variant))+
  theme_classic()+
  scale_color_manual(values=pal)+
  scale_shape_manual(values=c(17,19))+
  geom_errorbar(aes(ymin=lwrbnd,ymax=uprbnd,width=0.4,col=variant),
                alpha=0.5,position=position_dodge(width=0.4),size=1)+
  geom_point(size=4,position=position_dodge(width=0.4),aes(shape=variant,col=variant))+
  geom_text(aes(label=pred.char),hjust=rep(-0.3,nrow(ggdata2.mod1.nonImm)),position = position_dodge(width=0.4))+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=17,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(),
        legend.title.align = 0.5,
        strip.background = element_rect(color="white"),
        strip.placement = "inside",
        strip.text = element_text(size=15,face="bold"),
        panel.spacing = unit(2, "lines"))+
  xlab("Vaccination status")+ylab("Vaccine effectiveness (%)")+
  facet_wrap(~prev.inf)+
  # coord_cartesian(ylim=c(0,100),xlim=c(0.7,7.5))+
  geom_vline(aes(xintercept=1.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=2.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=5.5),alpha=0.8,linetype="dotted")
gg2.nonImm
pdf(file=paste0("VE_DeltaVOmicron_intPrevInf_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=14);gg2.nonImm;dev.off()

pdf(file=paste0("gamTerms_DeltaVOmicron_intPrevInf_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=8)
par(mfrow=c(2,2),mar=c(5.1,4.1,2.1,2.1))
plot(mod1.Delta.nonImm,select=1,scheme=1,bty="l")
plot(mod1.Delta.nonImm,select=2,scheme=1,bty="l")
plot(mod1.Omicron.nonImm,select=1,scheme=1,bty="l")
plot(mod1.Omicron.nonImm,select=2,scheme=1,bty="l")
par(mfrow=c(1,1),mar=c(5.1,4.1,2.1,2.1))
dev.off()

sink(paste0("mod1_Delta_DeltaVOmicron_intPrevInf_",period.study[1],"_to_",period.study[2],".txt"))
summary.mod1.Delta.nonImm
sink()

sink(paste0("mod1_Omicron_DeltaVOmicron_intPrevInf_",period.study[1],"_to_",period.study[2],".txt"))
summary.mod1.Omicron.nonImm
sink()

#

ggplot(data.rect.tested.2,aes(y=time.since.most.recent.vacc,x=status.vacc.previnf.int))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=17,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(),
        legend.title.align = 0.5,
        strip.background = element_rect(color="white"),
        strip.placement = "inside",
        strip.text = element_text(size=15,face="bold"),
        panel.spacing = unit(2, "lines"))


tapply(data.rect.tested.2$time.since.most.recent.vacc,
       data.rect.tested.2$status.vacc.previnf.int,
       summary)

boxplot(data.rect.tested.2$time.since.most.recent.vacc~
        data.rect.tested.2$status.vacc.previnf.int)

########################

# Add time since vaccination:

mod3.Delta.nonImm<-gam(infected~status.vacc.previnf.time.int+
                         sex+s(age_31Oct21)+s(SIMD_2016_VIGINTILE),
                       data=data.Delta.nonImm,family = binomial())
summary.mod3.Delta.nonImm<-summary(mod3.Delta.nonImm);summary.mod3.Delta.nonImm
plot(mod3.Delta.nonImm,pages=1)

mod3.Omicron.nonImm<-gam(infected~status.vacc.previnf.time.int+
                           sex+s(age_31Oct21)+s(SIMD_2016_VIGINTILE),
                         data=data.Omicron.nonImm,family = binomial())
summary.mod3.Omicron.nonImm<-summary(mod3.Omicron.nonImm);summary.mod3.Omicron.nonImm
plot(mod3.Omicron.nonImm,pages=1)

# Plots and output:

ggdata.mod3.nonImm<-cbind.data.frame(names=rep(c("Unvaccinated",
                                             rep("One dose only",2),
                                             rep(c("2 x ChAdOx1","2 x BNT162b2","2 x mRNA-1273"),4),
                                             rep(c("3rd dose BNT162b2","3rd dose mRNA-1273"),2)),2),
                                     prev.inf=rep(c("Previously confirmed infected",
                                                "No previous confirmed infection",
                                                "Previously confirmed infected",
                                                rep(c("No previous confirmed infection",
                                                    "Previously confirmed infected"),each=6),
                                                rep(c("No previous confirmed infection",
                                                      "Previously confirmed infected"),each=2)),2),
                                     time=rep(c(rep("",3),rep(c(rep(" <12 weeks",3),rep(" >=12 weeks",3)),2),rep("",4)),2),
                                     variant=rep(c("Delta","Omicron"),each=19),
                                     pred=c(100*(1-exp(summary.mod3.Delta.nonImm$p.coeff))[c(2:20)],
                                            100*(1-exp(summary.mod3.Omicron.nonImm$p.coeff))[c(2:20)]),
                                     lwrbnd=c(100*(1-exp(summary.mod3.Delta.nonImm$p.coeff+1.96*summary.mod3.Delta.nonImm$se[1:length(summary.mod3.Delta.nonImm$p.coeff)]))[c(2:20)],
                                              100*(1-exp(summary.mod3.Omicron.nonImm$p.coeff+1.96*summary.mod3.Omicron.nonImm$se[1:length(summary.mod3.Omicron.nonImm$p.coeff)]))[c(2:20)]),
                                     uprbnd=c(100*(1-exp(summary.mod3.Delta.nonImm$p.coeff-1.96*summary.mod3.Delta.nonImm$se[1:length(summary.mod3.Delta.nonImm$p.coeff)]))[c(2:20)],
                                              100*(1-exp(summary.mod3.Omicron.nonImm$p.coeff-1.96*summary.mod3.Omicron.nonImm$se[1:length(summary.mod3.Omicron.nonImm$p.coeff)]))[c(2:20)]))
ggdata.mod3.nonImm$names<-as.factor(ggdata.mod3.nonImm$names)
ggdata.mod3.nonImm$names<-factor(ggdata.mod3.nonImm$names,levels = unique(ggdata.mod3.nonImm$names)[c(seq(1,13,by=2),seq(2,12,by=2))])
ggdata.mod3.nonImm$variant<-as.factor(ggdata.mod3.nonImm$variant)
ggdata.mod3.nonImm$variant<-factor(ggdata.mod3.nonImm$variant,levels = unique(ggdata.mod3.nonImm$variant))

ggdata.mod3.nonImm$names.time<-paste0(ggdata.mod3.nonImm$names,ggdata.mod3.nonImm$time)

library(ggplot2)
pal<-c("#12B712","#918CF4")
hjust1<-rep(-0.3,nrow(ggdata.mod3.nonImm))

ggdata.mod3.nonImm$names.time<-as.factor(ggdata.mod3.nonImm$names.time)
levels(ggdata.mod3.nonImm$names.time)
ggdata.mod3.nonImm$names.time<-factor(ggdata.mod3.nonImm$names.time,
                                      levels=levels(ggdata.mod3.nonImm$names.time)[c(10,9,3,1,5,4,2,6,7,8)])

# pred as character, to avoid .0 being removed from text on plot:
ggdata.mod3.nonImm$pred.char<-as.character(round(ggdata.mod3.nonImm$pred,1))
ggdata.mod3.nonImm$pred.char[!grepl("\\.",ggdata.mod3.nonImm$pred.char)]<-
  paste0(ggdata.mod3.nonImm$pred.char[!grepl("\\.",ggdata.mod3.nonImm$pred.char)],".0")
# ... and set "100.0" to be "100":
ggdata.mod3.nonImm$pred.char[ggdata.mod3.nonImm$pred.char=="100.0"]<-"100"

# Add rows for unvaccinated, no previous infection, to create points at zero:
ggdata.mod3.nonImm<-rbind.data.frame(ggdata.mod3.nonImm,
                                      data.frame(names=rep("Unvaccinated",2),
                                                 variant=c("Delta","Omicron"),
                                                 pred=rep(0,2),
                                                 lwrbnd=rep(NA,2),
                                                 uprbnd=rep(NA,2),
                                                 prev.inf=rep("No previous confirmed infection",2),
                                                 pred.char=rep("",2),
                                                 time=rep("",2),
                                                 names.time=rep("Unvaccinated",2)))

gg2.nonImm<-ggplot(ggdata.mod3.nonImm,aes(x=names.time,y=pred,group=variant))+
  theme_classic()+
  scale_color_manual(values=pal)+
  scale_shape_manual(values=c(17,19))+
  geom_errorbar(aes(ymin=lwrbnd,ymax=uprbnd,width=0.4,col=variant),
                alpha=0.5,position=position_dodge(width=0.4),size=1)+
  geom_point(size=4,position=position_dodge(width=0.4),aes(shape=variant,col=variant))+
  geom_text(aes(label=pred.char),hjust=rep(-0.3,nrow(ggdata.mod3.nonImm)),position = position_dodge(width=0.4))+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=17,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(),
        legend.title.align = 0.5,
        strip.background = element_rect(color="white"),
        strip.placement = "inside",
        strip.text = element_text(size=15,face="bold"),
        panel.spacing = unit(2, "lines"))+
  xlab("Vaccination status")+ylab("Vaccine effectiveness (%)")+
  facet_wrap(~prev.inf)+
  # coord_cartesian(ylim=c(0,100),xlim=c(0.7,7.5))+
  geom_vline(aes(xintercept=1.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=2.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=5.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=8.5),alpha=0.8,linetype="dotted")
gg2.nonImm
pdf(file=paste0("VE_DeltaVOmicron_intPrevInf_time_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=18);gg2.nonImm;dev.off()

pdf(file=paste0("gamTerms_DeltaVOmicron_intPrevInf_time_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=8)
par(mfrow=c(2,2),mar=c(5.1,4.1,2.1,2.1))
plot(mod3.Delta.nonImm,select=1,scheme=1,bty="l")
plot(mod3.Delta.nonImm,select=2,scheme=1,bty="l")
plot(mod3.Omicron.nonImm,select=1,scheme=1,bty="l")
plot(mod3.Omicron.nonImm,select=2,scheme=1,bty="l")
par(mfrow=c(1,1),mar=c(5.1,4.1,2.1,2.1))
dev.off()

sink(paste0("mod3_Delta_DeltaVOmicron_intPrevInf_time_",period.study[1],"_to_",period.study[2],".txt"))
summary.mod3.Delta.nonImm
sink()

sink(paste0("mod3_Omicron_DeltaVOmicron_intPrevInf_time_",period.study[1],"_to_",period.study[2],".txt"))
summary.mod3.Omicron.nonImm
sink()

##################################################################################################################

#
# Summary tables
#

##################################################################################################################

# Age by positivity:
tapply(data.Delta.nonImm$age_31Oct21,data.Delta.nonImm$infected,summary)
tapply(data.Omicron.nonImm$age_31Oct21,data.Omicron.nonImm$infected,summary)

t.age_31Oct21<-round(cbind(Negative=tapply(data.Delta.nonImm$age_31Oct21,data.Delta.nonImm$infected,summary)[[1]],
                           Delta=tapply(data.Delta.nonImm$age_31Oct21,data.Delta.nonImm$infected,summary)[[2]],
                           Omicron=tapply(data.Omicron.nonImm$age_31Oct21,data.Omicron.nonImm$infected,summary)[[2]]),2)
write.csv(t.age_31Oct21,file=paste0("table_age_31Oct21_",period.study[1],"_to_",period.study[2],".csv"))

pdf(file=paste0("boxplot_age_31Oct21_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=5)
par(bty="l")
boxplot(data.Delta.nonImm$age_31Oct21[!data.Delta.nonImm$infected],
        data.Delta.nonImm$age_31Oct21[data.Delta.nonImm$infected],
        data.Omicron.nonImm$age_31Oct21[data.Omicron.nonImm$infected],
        xlab="Variant status",ylab="Age at 31st October 2021",names=c("Negative","Delta","Omicron"),
        col = c("grey",pal[1],pal[2]))
par(bty="o")
dev.off()

# Sex by positivity:
t.sex.Delta<-table(data.Delta.nonImm$sex,data.Delta.nonImm$infected)
t.sex.Omicron<-table(data.Omicron.nonImm$sex,data.Omicron.nonImm$infected)

pt.sex.Delta<-prop.table(t.sex.Delta,2)
pt.sex.Omicron<-prop.table(t.sex.Omicron,2)

table.sex<-cbind(matrix(paste0(t.sex.Delta," (",round(100*pt.sex.Delta,2),"%)"),ncol=2,nrow=2),
      matrix(paste0(t.sex.Omicron," (",round(100*pt.sex.Omicron,2),"%)"),ncol=2,nrow=2)[,2])
colnames(table.sex)<-c("Negative","Delta","Omicron")
rownames(table.sex)<-c("F","M")
table.sex
write.csv(table.sex,file=paste0("table_sex_",period.study[1],"_to_",period.study[2],".csv"))

# SIMD vigintile by positivity:
tapply(data.Delta.nonImm$SIMD_2016_VIGINTILE,data.Delta.nonImm$infected,summary)
tapply(data.Omicron.nonImm$SIMD_2016_VIGINTILE,data.Omicron.nonImm$infected,summary)

t.SIMD_2016_VIGINTILE<-round(cbind(Negative=tapply(data.Delta.nonImm$SIMD_2016_VIGINTILE,data.Delta.nonImm$infected,summary)[[1]],
      Delta=tapply(data.Delta.nonImm$SIMD_2016_VIGINTILE,data.Delta.nonImm$infected,summary)[[2]],
      Omicron=tapply(data.Omicron.nonImm$SIMD_2016_VIGINTILE,data.Omicron.nonImm$infected,summary)[[2]]),2)
write.csv(t.SIMD_2016_VIGINTILE,file=paste0("table_SIMD_2016_VIGINTILE_",period.study[1],"_to_",period.study[2],".csv"))

pdf(file=paste0("boxplot_SIMD_2016_VIGINTILE_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=5)
par(bty="l")
boxplot(data.Delta.nonImm$SIMD_2016_VIGINTILE[!data.Delta.nonImm$infected],
        data.Delta.nonImm$SIMD_2016_VIGINTILE[data.Delta.nonImm$infected],
        data.Omicron.nonImm$SIMD_2016_VIGINTILE[data.Omicron.nonImm$infected],
        xlab="Variant status",ylab="SIMD (2016) vigintile",names=c("Negative","Delta","Omicron"),
        col = c("grey",pal[1],pal[2]))
par(bty="o")
dev.off()

# Previous infection status:
t.previous.infection.Delta<-table(data.Delta.nonImm$previous.infection,data.Delta.nonImm$infected)
t.previous.infection.Omicron<-table(data.Omicron.nonImm$previous.infection,data.Omicron.nonImm$infected)

pt.previous.infection.Delta<-prop.table(t.previous.infection.Delta,2)
pt.previous.infection.Omicron<-prop.table(t.previous.infection.Omicron,2)

table.previous.infection<-cbind(matrix(paste0(t.previous.infection.Delta," (",round(100*pt.previous.infection.Delta,2),"%)"),ncol=2,nrow=2),
                 matrix(paste0(t.previous.infection.Omicron," (",round(100*pt.previous.infection.Omicron,2),"%)"),ncol=2,nrow=2)[,2])
colnames(table.previous.infection)<-c("Negative","Delta","Omicron")
rownames(table.previous.infection)<-c(FALSE,TRUE)
table.previous.infection
write.csv(table.previous.infection,file=paste0("table_previous_infection_",period.study[1],"_to_",period.study[2],".csv"))

# Vaccination status:
data.Delta.nonImm$status.vacc.f<-as.factor(data.Delta.nonImm$status.vacc)
levels(data.Delta.nonImm$status.vacc.f)
data.Delta.nonImm$status.vacc.f<-factor(data.Delta.nonImm$status.vacc.f,levels=levels(data.Delta.nonImm$status.vacc.f)[c(7,1,4,3,2,6,5)])

data.Omicron.nonImm$status.vacc.f<-as.factor(data.Omicron.nonImm$status.vacc)
levels(data.Omicron.nonImm$status.vacc.f)
data.Omicron.nonImm$status.vacc.f<-factor(data.Omicron.nonImm$status.vacc.f,levels=levels(data.Omicron.nonImm$status.vacc.f)[c(7,1,4,3,2,6,5)])

t.status.vacc.f.Delta<-table(data.Delta.nonImm$status.vacc.f,data.Delta.nonImm$infected)
t.status.vacc.f.Omicron<-table(data.Omicron.nonImm$status.vacc.f,data.Omicron.nonImm$infected)

pt.status.vacc.f.Delta<-prop.table(t.status.vacc.f.Delta,2)
pt.status.vacc.f.Omicron<-prop.table(t.status.vacc.f.Omicron,2)

table.status.vacc.f<-cbind(matrix(paste0(t.status.vacc.f.Delta," (",round(100*pt.status.vacc.f.Delta,2),"%)"),ncol=2,nrow=7),
                                matrix(paste0(t.status.vacc.f.Omicron," (",round(100*pt.status.vacc.f.Omicron,2),"%)"),ncol=2,nrow=7)[,2])
colnames(table.status.vacc.f)<-c("Negative","Delta","Omicron")
rownames(table.status.vacc.f)<-rownames(t.status.vacc.f.Delta)
table.status.vacc.f
write.csv(table.status.vacc.f,file=paste0("table_status_vacc_",period.study[1],"_to_",period.study[2],".csv"))

# time since most recent vaccine dose by variant and dose number:

t.time.0.Delta<-tapply((data.Delta.nonImm$time.since.most.recent.vacc+14)[data.Delta.nonImm$number.vacc==0],data.Delta.nonImm$infected[data.Delta.nonImm$number.vacc==0],summary)
t.time.1.Delta<-tapply((data.Delta.nonImm$time.since.most.recent.vacc+14)[data.Delta.nonImm$number.vacc==1],data.Delta.nonImm$infected[data.Delta.nonImm$number.vacc==1],summary)
t.time.2.Delta<-tapply((data.Delta.nonImm$time.since.most.recent.vacc+14)[data.Delta.nonImm$number.vacc==2],data.Delta.nonImm$infected[data.Delta.nonImm$number.vacc==2],summary)
t.time.3.Delta<-tapply((data.Delta.nonImm$time.since.most.recent.vacc+14)[data.Delta.nonImm$number.vacc==3],data.Delta.nonImm$infected[data.Delta.nonImm$number.vacc==3],summary)

t.time.0.Omicron<-tapply((data.Omicron.nonImm$time.since.most.recent.vacc+14)[data.Omicron.nonImm$number.vacc==0],data.Omicron.nonImm$infected[data.Omicron.nonImm$number.vacc==0],summary)
t.time.1.Omicron<-tapply((data.Omicron.nonImm$time.since.most.recent.vacc+14)[data.Omicron.nonImm$number.vacc==1],data.Omicron.nonImm$infected[data.Omicron.nonImm$number.vacc==1],summary)
t.time.2.Omicron<-tapply((data.Omicron.nonImm$time.since.most.recent.vacc+14)[data.Omicron.nonImm$number.vacc==2],data.Omicron.nonImm$infected[data.Omicron.nonImm$number.vacc==2],summary)
t.time.3.Omicron<-tapply((data.Omicron.nonImm$time.since.most.recent.vacc+14)[data.Omicron.nonImm$number.vacc==3],data.Omicron.nonImm$infected[data.Omicron.nonImm$number.vacc==3],summary)

table.time.variant<-round(cbind(Negative.1=t.time.1.Delta[[1]],Negative.2=t.time.2.Delta[[1]],Negative.3=t.time.3.Delta[[1]],
      Delta.1=t.time.1.Delta[[2]],Delta.2=t.time.2.Delta[[2]],Delta.3=t.time.3.Delta[[2]],
      Omicron.1=t.time.1.Omicron[[2]],Omicron.2=t.time.2.Omicron[[2]],Omicron.3=t.time.3.Omicron[[2]]),2)
write.csv(table.time.variant,file=paste0("table_time_variant_",period.study[1],"_to_",period.study[2],".csv"))

pdf(file=paste0("boxplot_time_since_vacc_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=5)
par(bty="l")
boxplot((data.Delta.nonImm$time.since.most.recent.vacc+14)[!data.Delta.nonImm$infected&data.Delta.nonImm$number.vacc==1],
        (data.Delta.nonImm$time.since.most.recent.vacc+14)[!data.Delta.nonImm$infected&data.Delta.nonImm$number.vacc==2],
        (data.Delta.nonImm$time.since.most.recent.vacc+14)[!data.Delta.nonImm$infected&data.Delta.nonImm$number.vacc==3],
        (data.Delta.nonImm$time.since.most.recent.vacc+14)[data.Delta.nonImm$infected&data.Delta.nonImm$number.vacc==1],
        (data.Delta.nonImm$time.since.most.recent.vacc+14)[data.Delta.nonImm$infected&data.Delta.nonImm$number.vacc==2],
        (data.Delta.nonImm$time.since.most.recent.vacc+14)[data.Delta.nonImm$infected&data.Delta.nonImm$number.vacc==3],
        (data.Omicron.nonImm$time.since.most.recent.vacc+14)[data.Omicron.nonImm$infected&data.Omicron.nonImm$number.vacc==1],
        (data.Omicron.nonImm$time.since.most.recent.vacc+14)[data.Omicron.nonImm$infected&data.Omicron.nonImm$number.vacc==2],
        (data.Omicron.nonImm$time.since.most.recent.vacc+14)[data.Omicron.nonImm$infected&data.Omicron.nonImm$number.vacc==3],
        names = rep(1:3,3),col=rep(c("grey",pal),each=3),
        ylab="Time since most recent vaccination (days)")
mtext("Negative",side=1,line=2,at=2)
mtext("Delta",side=1,line=2,at=5)
mtext("Omicron",side=1,line=2,at=8)
mtext("Variant and dose number",side=1,line=3.5)
par(bty="o")
dev.off()

## Main plot of vaccine dose by date:

data.vacc<-rbind.data.frame(cbind.data.frame(date=data.rect.all.2$date_1st_dose,product_name=data.rect.all.2$product_name_1st,dose=1,SafeHavenID=data.rect.all.2$SafeHavenID),
                            cbind.data.frame(date=data.rect.all.2$date_2nd_dose,product_name=data.rect.all.2$product_name_2nd,dose=2,SafeHavenID=data.rect.all.2$SafeHavenID),
                            cbind.data.frame(date=data.rect.all.2$date_3rd_dose,product_name=data.rect.all.2$product_name_3rd,dose=3,SafeHavenID=data.rect.all.2$SafeHavenID),
                            cbind.data.frame(date=data.rect.all.2$date_4th_dose,product_name=data.rect.all.2$product_name_4th,dose=4,SafeHavenID=data.rect.all.2$SafeHavenID),
                            cbind.data.frame(date=data.rect.all.2$date_5th_dose,product_name=data.rect.all.2$product_name_5th,dose=5,SafeHavenID=data.rect.all.2$SafeHavenID))
data.vacc<-data.vacc[which((data.vacc$product_name=="Covid-19 Vaccine AstraZeneca"|
                              data.vacc$product_name=="Covid-19 mRNA Vaccine Pfizer"|
                              data.vacc$product_name=="Covid-19 mRNA Vaccine Moderna")&
                             (data.vacc$dose==1|data.vacc$dose==2|data.vacc$dose==3)),]
data.vacc$product_name[data.vacc$product_name=="Covid-19 Vaccine AstraZeneca"]<-"ChAdOx1"
data.vacc$product_name[data.vacc$product_name=="Covid-19 mRNA Vaccine Pfizer"]<-"BNT162b2"
data.vacc$product_name[data.vacc$product_name=="Covid-19 mRNA Vaccine Moderna"]<-"mRNA-1273"
data.vacc$product_name<-as.factor(data.vacc$product_name)
data.vacc$product_name<-factor(data.vacc$product_name,levels=c("ChAdOx1","BNT162b2","mRNA-1273"))

# Remove third doses for ChAdOx1:
data.vacc<-data.vacc[-which(data.vacc$product_name=="ChAdOx1"&data.vacc$dose==3),]

# Summary: dimensions:
nrow(data.vacc) # 2330441
length(unique(data.vacc$SafeHavenID)) # 859130

pal.gg.vacc<-c("#F7FCFD","#BFD3E6","#8C96C6")

gg.vacc<-ggplot(data.vacc[data.vacc$date>=as.Date("2020-12-01"),],aes(y=date,x=factor(dose),fill=product_name))+geom_boxplot(outlier.size = 0.1)+
  theme_classic()+
  scale_fill_manual(values=pal.gg.vacc)+
  geom_point(position = position_jitterdodge(0.4),alpha=0.025,size=0.5)+
  scale_y_date(date_breaks = "months",date_labels ="%b-%y" )+
  xlab("Dose number")+ylab("Date")+labs(fill="Vaccine product")
pdf(file=paste0("boxplot_vacc_dose_by_date_",period.study[1],"_to_",period.study[2],".pdf"),height=5,width=8);gg.vacc;dev.off()
png(file=paste0("boxplot_vacc_dose_by_date_",period.study[1],"_to_",period.study[2],".png"),height=350,width=500);gg.vacc;dev.off()

# ---

# Summary stats for positivity and variant status by vaccination status plot:

length(unique(data.rect.all.2$SafeHavenID[!is.na(data.rect.all.2$date_1st_pos_PCR)|
                                          !is.na(data.rect.all.2$date_1st_pos_PCR_after90)|
                                          !is.na(data.rect.all.2$date_2nd_pos_PCR_after90)])) # 215805 (same without reinfections)
length(unique(data.rect.all.2$SafeHavenID)) # 1108260


