library(glmmTMB)
# Setup ####
PCs<-read.csv("Data/PC_Long_Removal_50m.csv")

Counts<-table(PCs$BIRD,PCs$Number)[,-1]
Counts<-Counts[-which(grepl("UN",rownames(Counts))),]
sort(rowSums(Counts)) # Numb Observations (Not Number of birds observed)

Sp<-names(which(rowSums(Counts)>=40))
# remove Sp with model convergence issues / unreliable due to limited sites: LEWO, MODO, RWBL
Sp<-Sp[-c(which(Sp=="LEWO"),which(Sp=="MODO"),which(Sp=="RWBL"))]
PCs$Year<-as.factor(PCs$Year)
# run models ####

for(i in 1:length(Sp)){
tryCatch({
  print(Sp[i])
  assign(paste0("NB",Sp[i]), glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  
                                                  scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+ # specify spectral overlap (and its interaction with dB) as test of masking
                                                  scale(Veg)+scale(Elev)+                 # specify percent riparian vegetation and elevation as covariates
                                                  scale(YDay)+scYDaySq+Year+
                                                  (1|SITE),                               # specify siteas random effect (intercept)
                                                data=PCs[PCs$Unimodal=="N"&PCs$BIRD==Sp[i],],
                                                family=nbinom1(link="log"),               # specify negative binomial distribution
                                                offset = log(p*CF_Offset)                 # specify offsets. p = removal model; CF = Correction factor (detecability experiment)
  )  
  )
  
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}, warning = function(w) {
  print(w)})
  
  }

## BRSP failed to converge as is, and dispersion parameter estimated to be 0. Change to poisson (special case of negative binomial). This will fix parameter so model isn't searching additional space.
NBBRSP<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  
          scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+
          scale(Veg)+scale(Elev)+                 
          scale(YDay)+scYDaySq+Year+
          (1|SITE),                               
        data=PCs[PCs$Unimodal=="N"&PCs$BIRD=="BRSP",],
        family=poisson(link="log"),               
        offset = log(p*CF_Offset)                 
)
NBWETA<-glmmTMB(Number ~ scale(PC_rem_L50)*scale(PC_rem_Med)+  
          scale(Freq_Diff)+ scale(PC_rem_L50):scale(Freq_Diff)+
          scale(Veg)+scale(Elev)+                 
          scale(YDay)+scYDaySq+Year+
          (1|SITE),                               
        data=PCs[PCs$Unimodal=="N"&PCs$BIRD=="WETA",],
        family=poisson(link="log"),               
        offset = log(p*CF_Offset)                 
)

# load models instead: ####

# save models:
# for(i in 1:length(Sp)){
# saveRDS(get(paste0("NB",Sp[i])),file = paste0("Model_Objects/Bird_Species_Models/PC_NB_",Sp[i],".rds"))
# }

# load models:
for(i in 1:length(Sp)){
  assign(paste0("NB",Sp[i]), readRDS(file = paste0("Model_Objects/Bird_Species_Models/PC_NB_",Sp[i],".rds")))
}

# Output for trait analysis ####

fixed<-c("Intercept","dB","kHz","FreqDiff","Veq","Elev","yDay","jDaySq","Year","dB:kHz","dB:FreqDiff")
what<-c("Est","SE","Z","P")
C_names<-"Species"
for(i in 1:length(what)){
  for(j in 1:length(fixed)){
    C_names<-c(C_names ,(paste0(fixed[j],"_",what[i])) )
  }
}


ModOut<-t(data.frame(rep(NA,45)))
colnames(ModOut)<-C_names


for(i in 1:length(Sp)){
  E<-format(round(as.numeric(summary(get(paste0("NB",Sp[i])))[["coefficients"]]$`cond`[,1]),4),scientific=FALSE)
  SE<-format(round(as.numeric(summary(get(paste0("NB",Sp[i])))[["coefficients"]]$`cond`[,2]),4),scientific=FALSE)
  Z<-format(round(as.numeric(summary(get(paste0("NB",Sp[i])))[["coefficients"]]$`cond`[,3]),4),scientific=FALSE)
  P<-format(round(as.numeric(summary(get(paste0("NB",Sp[i])))[["coefficients"]]$`cond`[,4]),4),scientific=FALSE)
  values<-c(Sp[i],E,SE,Z,P)
  ModOut<-rbind(ModOut,values)
  print(Sp[i])
}

ModOut<-ModOut[-1,]
ModOut<-as.data.frame(ModOut)
ModOut[,-1]<-sapply(ModOut[,-1],as.character)
ModOut[,-1]<-sapply(ModOut[,-1],as.numeric)

ModOut$Species<-as.character(ModOut$Species)
row.names(ModOut)<-c(1:length(ModOut$Species))

# write.csv(ModOut,"Output/Bird_Sp_ModelOutput.csv")


# building model table for supplement ####

ModOut<-read.csv("Output/Bird_Sp_ModelOutput.csv",row.names=1)


ModTable<-NULL
for(i in 1:length(Sp)){
add<-data.frame(round(summary(get(paste0("NB",Sp[i])))[["coefficients"]]$`cond`,2),Sp[i])
ModTable<-rbind(ModTable,add)
}

ModTable$Var<-row.names(ModTable)[1:11]
row.names(ModTable)<-c(1:length(ModTable$Estimate))
ModTable$Var<-gsub("scale\\((.*)\\)","\\1",ModTable$Var)
ModTable$Var<-gsub("(.*)\\):scale\\((.*)","\\1:\\2",ModTable$Var)
ModTable$Var<-gsub("PC_rem_L50","SPL",ModTable$Var)
ModTable$Var<-gsub("PC_rem_Med","Background Freq",ModTable$Var)
ModTable$Var<-gsub("Freq_Diff","Spectral Overlap",ModTable$Var)
ModTable$Var<-gsub("scYDaySq","Day2",ModTable$Var)
ModTable$Var<-gsub("YDay","Day",ModTable$Var)
ModTable$Var<-gsub("Year2018","Year",ModTable$Var)
ModTable$Var<-gsub("Veg","Vegetation",ModTable$Var)
ModTable$Var<-gsub("Elev","Elevation",ModTable$Var)
ModTable$Var<-gsub("\\((.*)\\)","\\1",ModTable$Var)

names(ModTable)<-c("Estimates","SE","Z","p","Species","Predictors")
Spnames<-read.csv("Keys/Bird_Species_Names.csv")
ModTable$Order<-c(1:length(ModTable$Estimates))
ModTable<-merge(ModTable,Spnames,by="Species",all.x=T)
ModTable<-ModTable[order(ModTable$Order),]
ModTable$p<-as.character(ModTable$p)
ModTable$p<-ifelse(ModTable$p=="0","<0.01",ModTable$p)
ModTable$CI<-paste0("'",round(ModTable$Estimates + 1.96*ModTable$SE,2)," - ",round(ModTable$Estimates - 1.96*ModTable$SE,2))
ModTable<-ModTable[,-7]

ModTable<-ModTable[,c(7,8,6,2,3,9,4,5)]
ModTable$common<-as.character(ModTable$common)
ModTable$Species.y.1<-as.character(ModTable$Species.y.1)
ModTable$common<-ifelse(duplicated(ModTable$common)==T,"", ModTable$common)
ModTable$Species.y.1<-ifelse(duplicated(ModTable$Species.y.1)==T,"", ModTable$Species.y.1)
names(ModTable)[c(1,2)]<-c("Common","Latin")

ModTable[names(ModTable)]<-sapply(ModTable[names(ModTable)],as.character)
test<-ModTable

start=1
new<-NULL

for(i in 1:length(Sp)){
new<-rbind(new,rep("",8),rep("",8),names(test), test[c(start:(start+10)),])
start=start+11
}

head(new)
# write.csv(new,"Output/Bird_Species_SM.csv", row.names=F)

# species observation table for supplement ####

Sp<-unique(PCs$BIRD)
Sp<-as.character(Sp)
Sp<-Sp[-which(grepl("UN",Sp))]
head(PCs)

SiteDay<-PCs %>% group_by(Year,Mo,Day,SiteID,BIRD) %>% 
  dplyr::summarise(Number=sum(Number))

Counts<-NULL
for(i in 1:length(Sp)){
  add<-c(Sp[i],
         sum(PCs$Number[PCs$BIRD==Sp[i]]), # total passes obs
         sum(ifelse(SiteDay$Number[SiteDay$BIRD==Sp[i]]>0,1,0  )))# unique days observed
  Counts<-rbind(Counts,add)
}

Counts<-as.data.frame(Counts)
names(Counts)<-c("Species","Observations","Unique site-days observed")

# Names<-read.csv("New_Order_scripting/Bird_Names.csv")
# Counts<-merge(Counts,Names,by="Species",all.x=T)

# write.csv(Counts,"Output/Bird_Species_Obs.csv")
