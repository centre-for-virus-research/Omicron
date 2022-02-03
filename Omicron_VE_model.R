
## R packages required:

library(ggplot2);library(RColorBrewer);library(mgcv)

## Datasets:

# data.rect.delta: Dataset containing data for those testing positive for Delta variant
# data.rect.omicron: Dataset containing data for those testing positive for Omicron variant

# Each dataset contains the following columns:
# . SafeHavenID:        ID of patient, generated for matching purposes (not identifiable outwith datasets).
# . infected:           TRUE if positive PCR result within study period; FALSE otherwise.
# . age_31Oct21:        Age of patient on 31st October 2021.
# . vax1_pre_study_bin: "yes" if patient received a 1st vaccine dose at least 14 days before the study period 
#                       and no further doses; "no" otherwise.
# . vax2_pre_study:     "AZ", "Pfizer" or "Moderna" if 2nd dose received at least 14 days before the study
#                       period was ChAdOx1, BNT162b2 or mRNA-1273 vaccine, respectively, with no further doses;
#                       "no" otherwise.
# . vax3_pre_study:     "Pfizer" or "Moderna" if 3rd dose received at least 14 days before the study period was 
#                       BNT162b2 or mRNA-1273 vaccine, respectively; "no" otherwise.
# . prev_infected:      TRUE if had positive PCR result at least 90 days before study period.
# . sex:                "F" for female or "M" for male.
# . SIMD_2016_QUARTILE: Scottish Index of Multiple Deprivation index (2016) quartile: 1 to 4 or NA.

## Fit models to Delta and Omicron data:

mod7=gam(infected~ s(age_31Oct21)+ vax1_pre_study_bin+vax2_pre_study+vax3_pre_study+
           sex+prev_infected+
           SIMD_2016_QUARTILE
         ,data=data.rect.delta
         ,family=binomial)

mod7.om=gam(infected~s(age_31Oct21)+vax1_pre_study_bin+vax2_pre_study+vax3_pre_study+
              sex+prev_infected+
              SIMD_2016_QUARTILE
            ,data=data.rect.omicron
            ,family=binomial)


mod=preds=list()
mod[[1]]=mod7
mod[[2]]=mod7.om

## VE estimates:

for(variant in 1:2){
  tab_rows=9
  predict_table=data.frame(sex=rep('M',tab_rows),
                               age_31Oct21=rep(60,tab_rows),
                               vax1_pre_study_bin=c('no','yes',rep('no',tab_rows-2)),
                               prev_infected=c(rep(F,tab_rows-1),T),
                               vax2_pre_study=c('none','none','AZ','Pfizer','Moderna',rep('none',tab_rows-5)),
                               vax3_pre_study=c(rep('none',5),'Moderna','Pfizer','AZ',rep('none',tab_rows-8)),
                               SIMD_2016_QUARTILE=factor(rep(2,tab_rows,levels=1:4)))
  
  # built with flexibility for vax3 terms - 3xAZ not needed
  predict_table=predict_table[c(1:7,9),]
  predict_table$plot_names=c('unvax','One dose only','2 x ChAdOx1','2 x BNT162b2',
                                 '2 x mRNA-1273','3rd dose mRNA-1273','3rd dose BNT162b2','Previously confirmed\ninfected unvaccinated')

  predict_table$plot_names=factor(predict_table$plot_names,levels=predict_table$plot_names)
  
  pred2=predict(mod[[variant]],predict_table,se.fit=T)
  ilink=family(mod[[variant]])$linkinv
  
  
  predict_table$pred=ilink(pred2$fit)
  
  predict_table$pred_upper=ilink(pred2$fit+2*pred2$se.fit)
  predict_table$pred_lower=ilink(pred2$fit-2*pred2$se.fit)

  predict_table$pred_rel=predict_table$pred/predict_table$pred[1]
  predict_table$pred_rel_upper=predict_table$pred_upper/predict_table$pred[1]
  predict_table$pred_rel_lower=predict_table$pred_lower/predict_table$pred[1]
  predict_table$variant=rep(variant,nrow(predict_table))
  preds[[variant]]=predict_table
}



predict_table_out.full=rbind.data.frame(preds[[1]],preds[[2]])
predict_table_out.full$variant=factor(var_names[predict_table_out.full$variant],levels=var_names)

# Reshuffle rows for plotting:

predict_table_out.full.2<-predict_table_out.full[-c(1,9),]
predict_table_out.full.2[1,]<-predict_table_out.full[8,]
predict_table_out.full.2[2:7,]<-predict_table_out.full[c(2:5,7,6),]
predict_table_out.full.2[8,]<-predict_table_out.full[16,]
predict_table_out.full.2[9:14,]<-predict_table_out.full[c(10:13,15,14),]

predict_table_out.full.2$plot_names<-factor(predict_table_out.full.2$plot_names,levels = predict_table_out.full.2$plot_names[1:7])

## Plot VE estimates:

pal<-c("#12B712","#918CF4")
hjust1<-rep(-0.3,14)

gg1<-ggplot(predict_table_out.full.2,aes(x=plot_names,y=100-100*pred_rel,group=variant))+
  theme_classic()+scale_color_manual(values=pal)+
  geom_errorbar(aes(ymin=100-100*pred_rel_lower,ymax=100-100*pred_rel_upper,width=0.4,col=variant),
                alpha=0.5,position=position_dodge(width=0.4),size=1)+
  geom_point(size=4,position=position_dodge(width=0.4),aes(shape=variant,col=variant))+
  geom_text(aes(label=round(100-100*pred_rel,2)),hjust=hjust1,position = position_dodge(width=0.4))+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=17,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(),
        legend.title.align = 0.5,
        strip.background = element_rect(color="white"),
        strip.placement = "inside",
        strip.text = element_text(size=15,face="bold"))+
  xlab("Vaccination status")+ylab("Vaccine effectiveness (%)")+
  coord_cartesian(ylim=c(0,100),xlim=c(0.7,7.5))+
  geom_vline(aes(xintercept=1.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=2.5),alpha=0.8,linetype="dotted")+
  geom_vline(aes(xintercept=5.5),alpha=0.8,linetype="dotted")
gg1

