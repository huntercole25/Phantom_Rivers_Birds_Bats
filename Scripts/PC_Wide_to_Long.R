library(tidyr)
PC<-read.csv("Data/PC_Wide_Removal_50m.csv",row.names = 1)

firstBird<-which(names(PC)=="AMCR")
lastBird<-which(names(PC)=="YRWA")

Grouping<-names(PC)[1:which(names(PC)=="AMCR")-1]
start<-which(names(PC)=="AMCR")
end<-which(names(PC)=="YRWA")

PC_L<-gather(PC, key="BIRD", value="Number", firstBird:lastBird)
table(PC_L$BIRD,PC_L$Number)
PC_L<-PC_L[,-c(which(names(PC_L)=="MINUTE"))]

PCAgg<-aggregate(Number~.,data=PC_L,FUN=sum)

## Merge in peak freq from bird traits
PR_Freq<-read.csv("Data/Bird_Traits.csv",row.names = 1) # Read in species data (Peak frequency etc.)

PCAgg$BIRD<-as.character(PCAgg$BIRD) # change to character to avoid merge issues
PR_Freq$Code4<-as.character(PR_Freq$Code4) # change to character to avoid merge issues
PCs<-merge(PCAgg, PR_Freq, by.y="Code4", by.x="BIRD", all.x=T) # Merge data

PCs$Freq_Diff<-abs(PCs$PC_rem_Med/1000-PCs$Peak_Freq) # create frequency overlap metric (in kHz)

PCs<-PCs[,-which(names(PCs)=="Total")] # remove total column, since it doesn't make sense in a long format

# write.csv(PCs,"Data/PC_Long_Removal_50m.csv",row.names = F)

