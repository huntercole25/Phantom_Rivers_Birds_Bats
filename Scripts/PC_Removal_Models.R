# See: https://peter.solymos.org/code/2018/08/30/fitting-removal-models-with-the-detect-r-package.html for a detailed description of the removal model, and how to fit them. 

library(tidyr)
library(glmmTMB)
library(detect)
library(mefa4)
library(dplyr)
library(ggplot2)

PC<-read.csv("Data/PC_Wide_50m.csv", row.names=1)
PC$scYDaySq<-(scale(PC$YDay)^2) # create quadratic versions of ordinal date
PC$scTASRSq<-(scale(PC$TASR)^2) # create quadratic versions of time after sunrise

firstBird<-which(names(PC)=="AMCR")
lastBird<-which(names(PC)=="YRWA")
PC$Total<-rowSums(PC[,c(firstBird:lastBird)])
PC_R<-PC[,-c(firstBird:lastBird)]
PC_R<-PC_R[PC_R$Unimodal=="N",]

# Removal Model####
PC_R$MINUTE<-as.factor(PC_R$MINUTE)

xt <- Xtab(Total ~ SiteID + MINUTE, PC_R)
Y <- as.matrix(xt[,c("1", "2", "3")])
X <- nonDuplicated(PC_R, PC_R$SiteID, TRUE)[rownames(Y),]
h <- !is.na(X$YDay)
Y <- Y[h,]
X <- X[,c("YDay", "TASR","LEQA")]
D <- matrix(c(1 , 2, 3), nrow(Y), 3, byrow=TRUE)
dimnames(D) <- dimnames(Y)

X$scYDaySq<-(scale(X$YDay)^2)
X$scTASRSq<-(scale(X$TASR)^2)

Me0 <- cmulti(Y | D ~ 1, X, type="rem")
Me1 <- cmulti(Y | D ~ scale(YDay), X, type="rem")
Me2 <- cmulti(Y | D ~ scale(TASR), X, type="rem")
Me3 <- cmulti(Y | D ~ scale(YDay)+scale(TASR), X, type="rem")
Me4 <- cmulti(Y | D ~ scYDaySq, X, type="rem")
Me5 <- cmulti(Y | D ~ scTASRSq, X, type="rem")
Me6 <- cmulti(Y | D ~ scTASRSq+scYDaySq, X, type="rem")
Me7 <- cmulti(Y | D ~ scale(YDay)+scYDaySq, X, type="rem")
Me8 <- cmulti(Y | D ~ scTASRSq+scale(YDay), X, type="rem")
Me9 <- cmulti(Y | D ~ scYDaySq+scale(TASR), X, type="rem")
Me10 <- cmulti(Y | D ~ scTASRSq+scale(TASR), X, type="rem")
Me11 <- cmulti(Y | D ~ scTASRSq+scale(YDay)+scale(TASR), X, type="rem")
Me12 <- cmulti(Y | D ~ scYDaySq+scale(YDay)+scale(TASR), X, type="rem")
Me13 <- cmulti(Y | D ~ scTASRSq+scYDaySq+scale(YDay)+scale(TASR), X, type="rem")
Me14 <- cmulti(Y | D ~ scale(LEQA), X, type="rem")
Me15 <- cmulti(Y | D ~ scale(YDay)+ scale(LEQA), X, type="rem")
Me16 <- cmulti(Y | D ~ scale(TASR)+ scale(LEQA), X, type="rem")
Me17 <- cmulti(Y | D ~ scale(YDay)+scale(TASR)+ scale(LEQA), X, type="rem")
Me18 <- cmulti(Y | D ~ scYDaySq+ scale(LEQA), X, type="rem")
Me19 <- cmulti(Y | D ~ scTASRSq+ scale(LEQA), X, type="rem")
Me20 <- cmulti(Y | D ~ scTASRSq+scYDaySq+ scale(LEQA), X, type="rem")
Me21 <- cmulti(Y | D ~ scale(YDay)+scYDaySq+ scale(LEQA), X, type="rem")
Me22 <- cmulti(Y | D ~ scTASRSq+scale(YDay)+ scale(LEQA), X, type="rem")
Me23 <- cmulti(Y | D ~ scYDaySq+scale(TASR)+ scale(LEQA), X, type="rem")
Me24 <- cmulti(Y | D ~ scTASRSq+scale(TASR)+ scale(LEQA), X, type="rem")
Me25 <- cmulti(Y | D ~ scTASRSq+scale(YDay)+scale(TASR)+ scale(LEQA), X, type="rem")
Me26 <- cmulti(Y | D ~ scYDaySq+scale(YDay)+scale(TASR)+ scale(LEQA), X, type="rem")
Me27 <- cmulti(Y | D ~ scTASRSq+scYDaySq+scale(YDay)+scale(TASR)+ scale(LEQA), X, type="rem")

Me_AIC <- AIC(Me0, Me1, Me2,Me3,Me4,Me5,Me6,Me7,Me8,Me9,Me10,Me11,Me12,Me13,
              Me14,Me15,Me16,Me17,Me18,Me19,Me20,Me21,Me22,Me23,Me24,Me25,Me26,Me27)
Me_AIC$dAIC <- Me_AIC$AIC - min(Me_AIC$AIC)

# Creating AIC table
# Me_AIC<-round(Me_AIC,2)
# write.csv(Me_AIC,"Output/Removal_AIC_table_50m.csv")

b <- coef(Me5) # 50 m
bLEQ <- coef(Me19) # 50 m
## Mod 5; TASR**2 only
PC$p <- 1-exp(-exp(
                  b[1] +
                  b[2]*PC$scTASRSq
                     )*3)
## Mod 19; TASR**2 + TASR + LEQ
PC$pLEQ <- 1-exp(-exp(
  bLEQ[1] +
    bLEQ[2]*PC$scTASRSq+
    bLEQ[3]*scale(PC$LEQA)             
                     )*3)


# Plot detectability ~ TASR ####
ggplot(aes(x=(PC$TASR*24), y=pLEQ, color=LEQA), data=PC)+
  geom_point()+
  theme_classic()+
  labs(x="TASR")


# Add a column for offsets, due to detectability in noise 
PC$CF_Offset<-1/(1+exp(-(11.841189509-0.154101572*PC$LEQA))) # this equation is based off of our detectability experiment - see supplement

# write.csv(PC, "Data/PC_Wide_Removal_50m.csv")
