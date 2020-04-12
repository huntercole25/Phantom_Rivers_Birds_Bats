## Load packages ####
library(ggplot2)
library(ggeffects)
library(DHARMa)
library(glmmTMB)
library(dplyr)
library(performance)
library(effects)

## Setup ####
PCs<-read.csv("Data/PC_Long_Removal_50m.csv")
Sites<-read.csv("Data/SiteID.csv",row.names = 1) # load in site treatment data

vars<-c("SITE","POINT")
Sites[vars]<-sapply(Sites[vars],as.character) # set vars as character for merge
PCs[vars]<-sapply(PCs[vars],as.character) # set vars as character for merge

PCs<-merge(PCs,Sites,by=vars,all.x=T) # add treatment data to PC data
PCs$Year<-as.factor(PCs$Year) # set year as factor for models

# Correlation matrix to see if model covariates are correlated:
PCc<-PCs[,c(15:18,20,22,31)]
cor(na.omit(PCc)) # nothing is correlated strongly
rm(PCc)

## global model ####

NB_global<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  # specify SPL (dB), median frequency (Hz) and the interaction between them as predictors
                      scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+ # specify spectral overlap (and its interaciton with dB) as test of masking
                      scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                      scale(YDay)+scYDaySq+Year+              # specify ordinal date, a quadratic version of that, and year as covariates
                      (1|SITE)+(1|BIRD),                      # specify site and bird as random effects (intercepts) 
                    data=PCs,
                    family=nbinom1(link="log"),               # specify negative binomial distribution
                    offset = log(p*CF_Offset)                 # specify offsets. p = removal model; CF = Correction factor (detecability experiment)
)              
# save(NB_global,file = "Model_Objects/PC_NB_global.rda")
load(file = "Model_Objects/PC_NB_global.rda")
check_collinearity(NB_global)
# plot(simulateResiduals(NB_global))
summary(NB_global)

## Percent change from global model ####


Int<-summary(NB_global)$coeff$cond[1,1]
dB_Beta<-summary(NB_global)$coeff$cond[2,1]
dB_se<-summary(NB_global)$coeff$cond[2,2]
SO_Beta<-summary(NB_global)$coeff$cond[4,1]
SO_se<-summary(NB_global)$coeff$cond[4,2]

# This is the difference in bird abundance over 12 dB (from 48 to 60 dB in this case) divided by the original estimate.. the result is the percentage decrease for 12 dB
(exp(Int+dB_Beta*((48-mean(PCs$PC_rem_L50))/sd(PCs$PC_rem_L50)))-
exp(Int+dB_Beta*((60-mean(PCs$PC_rem_L50))/sd(PCs$PC_rem_L50))) ) / # intercept + dB_Beta * X. Here we center our dB value X, by subtracting the data mean and divide by sd of original data to scale it 
  exp(Int+dB_Beta*((48-mean(PCs$PC_rem_L50))/sd(PCs$PC_rem_L50)))

-(1-exp(dB_Beta/sd(PCs$PC_rem_L50))^12) # This should be roughly the same value as above. Here the calculated value is per decibel, and multiplied by 12 at the end to put into terms of every 12 dB.
-(1-exp((dB_Beta-1.96*dB_se)/sd(PCs$PC_rem_L50))^12) # Upper 95% bound
-(1-exp((dB_Beta+1.96*dB_se)/sd(PCs$PC_rem_L50))^12) # Lower 95% bound
# -(1-exp(dB_Beta-1.96*dB_se))/sd(PCs$PC_rem_L50)*12 # Upper 95% bound # can also calculate it this way, which gives effectively the same answer
# -(1-exp(dB_Beta+1.96*dB_se))/sd(PCs$PC_rem_L50)*12 # Lower 95% bound # can also calculate it this way, which gives effectively the same answer

## spectral overlap
-(1-exp(SO_Beta/sd(PCs$Freq_Diff,na.rm=T))^2) # 
-(1-exp((SO_Beta-1.96*SO_se)/sd(PCs$Freq_Diff,na.rm=T))^2) # Upper 95% bound
-(1-exp((SO_Beta+1.96*SO_se)/sd(PCs$Freq_Diff,na.rm=T))^2) # Upper 95% bound


## Plotting global model: ####

PC_plot<-PCs
PC_plot <- PC_plot %>% group_by(SITE,SiteID,Year,Mo,Day,
                            Elev,Veg) %>% 
  dplyr::summarise(
    YDay=mean(YDay),
    scYDaySq=mean(scYDaySq),
    p=mean(p),
    PC_rem_Med=mean(PC_rem_Med),
    Freq_Diff=mean(Freq_Diff,na.rm=T),
    PC_rem_L50=mean(PC_rem_L50),
    CF_Offset=mean(CF_Offset),
    Total=sum(Number)
  )

# need to get a new intercept for the plot, since I am aggregating counts for visualization
NBint = glmmTMB(Total ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  # specify dB, Q1, and the interaction between them as predictors
                  scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+
                  scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                   as.factor(Year)+scale(YDay)+scYDaySq+   # specify year as categorical predictor and ordinal date (and date squared) as covariate
                   (1|SITE),                        # specify site and bird as random effect 
                 data=PC_plot, 
                 family=nbinom1(link="log"),               # specify negative binomial distribution
                offset = log(p*CF_Offset)                 # specify offsets. p = removal model; CF = Correction factor (detecability experiment)
) 


# Model predicitons
dBEst<-summary(NB_global)$coefficients$cond[2,1]
dBSE<-summary(NB_global)$coefficients$cond[2,2]
Int<-summary(NBint)$coefficients$cond[1,1]
IntSE<-summary(NBint)$coefficients$cond[1,2]

xP<-seq(from=25,to=82,length.out=length(PC_plot$PC_rem_L50))
yP<-exp(Int+dBEst*((xP-mean(PC_plot$PC_rem_L50))/sd(PC_plot$PC_rem_L50)))
yPup<-exp(Int+IntSE+dBEst*((xP-mean(PC_plot$PC_rem_L50))/sd(PC_plot$PC_rem_L50))+dBSE)
yPdown<-exp(Int-IntSE+dBEst*((xP-mean(PC_plot$PC_rem_L50))/sd(PC_plot$PC_rem_L50))-dBSE)

ggplot(aes(y=Total, x=PC_rem_L50), data=PC_plot)+
  theme_classic()+
  theme(axis.title = element_text(size=48),
        axis.ticks.length=unit(.3, "cm"),
        axis.ticks = element_line(size=2),
        axis.text = element_text(size=40),
        axis.line = element_line(colour = 'black', size = 3),
        axis.text.y = element_text(margin = margin(r=10)),
        axis.text.x = element_text(vjust=-0.5))+
  geom_point(alpha=0.7,size=6)+
  # stat_smooth(method = "glm",color="red",se=T)+
  coord_cartesian(ylim=c(0,29.7))+
  geom_line(aes(y=yP,x=xP),size=1.25)+
  geom_ribbon(aes(x=xP,ymin=yPdown,ymax=yPup), alpha=0.3)+
  scale_x_continuous(breaks=seq(30,80,by=10),expand = c(0, .5)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y=bquote(paste("Birds location "^-1," count day "^-1,"\n")), x="\n Sound pressure level (dB)")


FDEst<-summary(NB_global)$coefficients$cond[4,1]
FDSE<-summary(NB_global)$coefficients$cond[4,2]

xP<-seq(from=1,to=7,length.out=length(PC_plot$Freq_Diff))
yP<-exp(Int+FDEst*((xP-mean(PC_plot$Freq_Diff))/sd(PC_plot$Freq_Diff)))
yPup<-exp(Int+IntSE+FDEst*((xP-mean(PC_plot$Freq_Diff))/sd(PC_plot$Freq_Diff))+FDSE)
yPdown<-exp(Int-IntSE+FDEst*((xP-mean(PC_plot$Freq_Diff))/sd(PC_plot$Freq_Diff))-FDSE)

ggplot(aes(y=Total, x=Freq_Diff), data=PC_plot)+
  theme_classic()+
  theme(axis.title = element_text(size=48),
        axis.ticks.length=unit(.3, "cm"),
        axis.ticks = element_line(size=2),
        axis.text = element_text(size=40),
        axis.line = element_line(colour = 'black', size = 3),
        axis.text.y = element_text(margin = margin(r=10)),
        axis.text.x = element_text(vjust=-0.5))+
  geom_point(alpha=0.7,size=6)+
  coord_cartesian(ylim=c(0,29.7))+
  geom_line(aes(y=yP,x=xP),size=1.25)+
  geom_ribbon(aes(x=xP,ymin=yPdown,ymax=yPup), alpha=0.3)+
  scale_x_continuous(breaks=seq(0,10,by=1)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y=bquote(paste("Birds location "^-1," count day "^-1,"\n")), x="\n Spectral overlap (kHz)")



## Check alternative method of putting LEQ into removal model: ####
NB_global2<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  # specify dB, Q1, and the interaction between them as predictors
                      scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+
                      scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                       as.factor(Year)+scale(YDay)+scYDaySq+   # specify year as categorical predictor and ordinal date (and date squared) as covariate
                       (1|SITE)+(1|BIRD),                        # specify site and bird as random effect 
                     data=PCs, 
                     family=nbinom1(link="log"),               # specify negative binomial distribution
                     offset = log(pLEQ)                 # specify offsets. p = removal model; CF = Correction factor (detecability experiment)
                     )              

# save(NB_global2,file = "Model_Objects/PC_NB_global2.rda")
load(file = "Model_Objects/PC_NB_global2.rda")
summary(NB_global2)
# plot(simulateResiduals(NB_global2))


## Visual only counts ####
NB_Visual<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+
                     scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+
                     scale(Veg)+scale(Elev)+                 
                       scale(YDay)+scYDaySq+Year+   
                       (1|SITE)+(1|BIRD),           
                     data=PCs[PCs$Unimodal=="Y",],            # specify visual only counts
                     family=nbinom1(link="log")     
)                   # remove offsets, which don't apply since they are controlling for the effects of hearing birds - these are visual, only, counts

# save(NB_Visual,file = "Model_Objects/PC_NB_Visual.rda")
load(file = "Model_Objects/PC_NB_Visual.rda")
summary(NB_Visual)

## Control only counts ####
NB_Control<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+ 
                     scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+
                     scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                     scale(YDay)+scYDaySq+Year+   # specify year as categorical predictor and ordinal date (and date squared) as covariate
                     (1|SITE)+(1|BIRD),                     # specify site and bird as random effect 
                   data=PCs[PCs$TREAT=="Control",],  # Control only counts
                   family=nbinom1(link="log"),               # specify negative binomial distribution
                   offset = log(p*CF_Offset) 
) 
# save(NB_Control,file = "New_Order_scripting/NB_Control.rda")
load(file = "Model_Objects/PC_NB_Control.rda")
summary(NB_Control)



## Create table output for SM:####
m1<-data.frame(summary(NB_global)$coeff$cond,"Mod"="Global","Noise"="Observer Experiment")
m2<-data.frame(summary(NB_global2)$coeff$cond,"Mod"="Global","Noise"="Noise Removal Model")
m3<-data.frame(summary(NB_Visual)$coeff$cond,"Mod"="Visual Counts","Noise"="Observer Ear plugs")
m4<-data.frame(summary(NB_Control)$coeff$cond,"Mod"="Control Sites","Noise"="")

ModTable<-rbind(m1,m2,m3,m4)
ModTable$Params<-row.names(ModTable)
row.names(ModTable)<-c(1:length(ModTable$Mod))

head(ModTable)
names(ModTable)<-c("Estimate","SE","z value","p value","Model","Noise Control","Variable")
ModTable<-ModTable[,c(5,7,1:4,6)]
ModTable$Variable<-gsub("scale\\((.*)\\)","\\1",ModTable$Variable)
ModTable$Variable<-ModTable$Variable[1:11]
ModTable$Variable<-gsub("(.*)\\):scale\\((.*)","\\1:\\2",ModTable$Variable)
ModTable$Variable<-gsub("sc","",ModTable$Variable)
ModTable$Variable<-gsub("2018","",ModTable$Variable)
ModTable$Variable<-gsub("PC_rem_L50","SPL",ModTable$Variable)
ModTable$Variable<-gsub("Freq_Diff","Spectral Overlap",ModTable$Variable)
ModTable$Variable<-gsub("PC_rem_Med","Background freq.",ModTable$Variable)
ModTable$Variable<-gsub("\\(Intercept\\)","Intercept",ModTable$Variable)

ModTable$CI<-paste0("'",round(ModTable$Estimate - 1.96*ModTable$SE,3)," - ",round(ModTable$Estimate + 1.96*ModTable$SE,2))

ModTable[,c(3:6)]<-round(ModTable[,c(3:6)],3)
ModTable$`p value`<-ifelse(ModTable$`p value`<0.001,"<0.001",ModTable$`p value`)
head(ModTable)

# write.csv(ModTable, "Output/Bird_Models_Table.csv")
## Deciding on models ####

# NBINOM1 vs NBINOM2, which is a better fit?
NB_global<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  # specify SPL (dB), median frequency (Hz) and the interaction between them as predictors
                     scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+ # specify spectral overlap (and its interaciton with dB) as test of masking
                     scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                     scale(YDay)+scYDaySq+Year+              # specify ordinal date, a quadratic version of that, and year as covariates
                     (1|SITE)+(1|BIRD),                      # specify site and bird as random effects (intercepts) 
                   data=PCs,
                   family=nbinom1(link="log"),               # specify negative binomial distribution
                   offset = log(p*CF_Offset)                 # specify offsets. p = removal model; CF = Correction factor (detecability experiment)
) 

NB2_global<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  # specify SPL (dB), median frequency (Hz) and the interaction between them as predictors
                     scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+ # specify spectral overlap (and its interaciton with dB) as test of masking
                     scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                     scale(YDay)+scYDaySq+Year+              # specify ordinal date, a quadratic version of that, and year as covariates
                     (1|SITE)+(1|BIRD),                      # specify site and bird as random effects (intercepts) 
                   data=PCs,
                   family=nbinom2(link="log"),               # specify negative binomial distribution
                   offset = log(p*CF_Offset)                 # specify offsets. p = removal model; CF = Correction factor (detecability experiment)
) 

AIC(NB_global,NB2_global)
#            df      AIC
# NB_global  14 49980.12 #nbinom1 is a better fit
# NB2_global 14 51449.40


# # Are these data Zero-inflated?
ZINB = glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  # specify SPL (dB), median frequency (Hz) and the interaction between them as predictors
                 scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+ # specify spectral overlap (and its interaciton with dB) as test of masking
                 scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                 scale(YDay)+scYDaySq+Year+              # specify ordinal date, a quadratic version of that, and year as covariates
                 (1|SITE)+(1|BIRD),                      # specify site and bird as random effects (intercepts) 
               data=PCs,              
               ziformula=~1,
               family=nbinom1(link="log"), 
               offset = log(p*CF_Offset)                 # specify offsets. p = removal model; CF = Correction factor (detecability experiment)
)


AIC(NB_global,ZINB)
#           df      AIC
# NB_global 14 49980.12 # No they are not
# ZINB      15 49982.12