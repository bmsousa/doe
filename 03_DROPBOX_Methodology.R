
#
# Please Use this script with paper: "Enhancing Path Selection in Multihomed Nodes"
# Otherwise it can be a bit confusing
# 
#

require("reshape2")

parseRes_lm <- function(iModel){
  auxSum <- summary(iModel)
  auxAIC <- AIC(iModel)
  
  return(data.frame(Rsquare=auxSum$r.squared, RAdjsquare=auxSum$adj.r.squared, FStatistic=auxSum$fstatistic[1], aic=auxAIC))
}



# 
# Due to private reasons, the code of TOPSIS, METHODICAL and DiA approaches 
# cannot be yet made available.
#
# To use this script, please load data file:
# 
load("demoData.data")


dftempANA_METHv10_raw <- melt(dfMETHODICAL_v10, id.vars=1:7, measure.vars=c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16"), variable.name="Replication", value.name="score")
dftempANA_TOPSIS_raw <- melt(dfTOPSIS, id.vars=1:7, measure.vars=c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16"), variable.name="Replication", value.name="score")
dftempANA_DiA_raw <- melt(dfDiA, id.vars=1:7, measure.vars=c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16"), variable.name="Replication", value.name="score")


# If ANOVA assumptions are not fulfilled, transformations on data can be performed.
#
#Log transform
dftempANA_METHv10_raw$scoreLog <- log(dftempANA_METHv10_raw$score)
dftempANA_TOPSIS_raw$scoreLog <- log(dftempANA_TOPSIS_raw$score)
dftempANA_DiA_raw$scoreLog <- log(dftempANA_DiA_raw$score)

#sqrt transform
dftempANA_METHv10_raw$scoreSq <- sqrt(dftempANA_METHv10_raw$score)
dftempANA_TOPSIS_raw$scoreSq <- sqrt(dftempANA_TOPSIS_raw$score)
dftempANA_DiA_raw$scoreSq <- sqrt(dftempANA_DiA_raw$score)

dftempANA_METHv10_raw$scoreSq4 <- (dftempANA_METHv10_raw$score)^1/4
dftempANA_TOPSIS_raw$scoreSq4 <- (dftempANA_TOPSIS_raw$score)^1/4
dftempANA_DiA_raw$scoreSq4 <- (dftempANA_DiA_raw$score)^1/4


dftempANA_METHv10_raw$scoreSq3 <- (dftempANA_METHv10_raw$score)^1/3
dftempANA_TOPSIS_raw$scoreSq3 <- (dftempANA_TOPSIS_raw$score)^1/3
dftempANA_DiA_raw$scoreSq3 <- (dftempANA_DiA_raw$score)^1/3


#
#TOPSIS
#
ttest1TOP<- lm(score~BW*RTT*Jitter*Loss*Cov, data=dftempANA_TOPSIS_raw)
summary(ttest1TOP)

ttest1TOPSq<- lm(scoreSq~BW*RTT*Jitter*Loss*Cov, data=dftempANA_TOPSIS_raw)
summary(ttest1TOPSq)

ttest1TOPSq4<- lm(scoreSq4~BW*RTT*Jitter*Loss*Cov, data=dftempANA_TOPSIS_raw)
summary(ttest1TOPSq4)

ttest1TOPSq3<- lm(scoreSq3~BW*RTT*Jitter*Loss*Cov, data=dftempANA_TOPSIS_raw)
summary(ttest1TOPSq3)

TOPmodel_DROPBOX <- lm(score~BW+RTT+Jitter+Loss+Cov, data=dftempANA_TOPSIS_raw)
summary(TOPmodel_DROPBOX)
lmTOPmodel_DROPBOX <- parseRes_lm(TOPmodel_DROPBOX)

TOPmodel_MeTH_DROPBOX <- lm(score~BW+RTT+Jitter+Loss+Cov+
                          BW:Cov + 
                          BW:RTT:Cov + BW:Jitter:Cov +
                          BW:Loss:Cov +
                          BW:RTT:Jitter:Cov +
                          BW:RTT:Loss:Cov +
                          BW:Jitter:Loss:Cov
                        , data=dftempANA_TOPSIS_raw)
summary(TOPmodel_MeTH_DROPBOX)
lmTOPmodel_MeTH_DROPBOX <- parseRes_lm(TOPmodel_MeTH_DROPBOX)


#
# DiA
#
ttest1DiA<- lm(score~BW*RTT*Jitter*Loss*Cov, data=dftempANA_DiA_raw)
summary(ttest1DiA)

ttest1DiASq<- lm(scoreSq~BW*RTT*Jitter*Loss*Cov, data=dftempANA_DiA_raw)
summary(ttest1DiASq)

DiAmodel_DROPBOX <- lm(score~BW+RTT+Jitter+Loss+Cov, data=dftempANA_DiA_raw)
summary(DiAmodel_DROPBOX)
lmDiAmodel_DROPBOX <- parseRes_lm(DiAmodel_DROPBOX)
#
# the DiA model is equal to TOPSIS
DiAmodel_MeTH_DROPBOX <- lm(score~BW+RTT+Jitter+Loss+Cov+
                              BW:Cov + 
                              BW:RTT:Cov + BW:Jitter:Cov +
                              BW:Loss:Cov +
                              BW:RTT:Jitter:Cov +
                              BW:RTT:Loss:Cov +
                              BW:Jitter:Loss:Cov
                            , data=dftempANA_DiA_raw)
summary(DiAmodel_MeTH_DROPBOX)
lmDiAmodel_MeTH_DROPBOX <- parseRes_lm(DiAmodel_MeTH_DROPBOX)





#
# MetH
#
ttest1<- lm(score~BW*RTT*Jitter*Loss*Cov, data=dftempANA_METHv10_raw)
summary(ttest1)

ttest1Sq<- lm(scoreSq~BW*RTT*Jitter*Loss*Cov, data=dftempANA_METHv10_raw)
summary(ttest1Sq)


ttest1Sq4<- lm(scoreSq4~BW*RTT*Jitter*Loss*Cov, data=dftempANA_METHv10_raw)
summary(ttest1Sq4)


ttest1Sq3<- lm(scoreSq3~BW*RTT*Jitter*Loss*Cov, data=dftempANA_METHv10_raw)
summary(ttest1Sq3)


step(ttest1, direction='backward', criterion='AIC')
MeTHmodel_DROPBOX <- lm(score~BW+RTT+Jitter+Loss+Cov+
                          BW:Cov + 
                          BW:RTT:Cov + BW:Jitter:Cov +
                          BW:Loss:Cov +
                          BW:RTT:Jitter:Cov +
                          BW:RTT:Loss:Cov +
                          BW:Jitter:Loss:Cov
                          , data=dftempANA_METHv10_raw)
summary(MeTHmodel_DROPBOX)
lmMeTHmodel_DROPBOX <- parseRes_lm(MeTHmodel_DROPBOX)


MeTHmodel2_DROPBOX <- lm(score~BW+RTT+Jitter+Loss+Cov+
                          BW:RTT:Loss:Cov +
                          BW:Jitter:Loss:Cov
                        , data=dftempANA_METHv10_raw)
summary(MeTHmodel2_DROPBOX)
lmMeTHmodel2_DROPBOX <- parseRes_lm(MeTHmodel2_DROPBOX)




MeTHmodel_TOP_DROPBOX <- lm(score~BW+RTT+Jitter+Loss+Cov, data=dftempANA_METHv10_raw)
summary(MeTHmodel_TOP_DROPBOX)
lmMeTHmodel_TOP_DROPBOX <- parseRes_lm(MeTHmodel_TOP_DROPBOX)

lmAll <- NULL
lmAll <- rbind(lmAll, data.frame(method="TOPSIS", model="lmTOP", allsig="yes", interactions="no", lmTOPmodel_DROPBOX ))
lmAll <- rbind(lmAll, data.frame(method="DiA",    model="lmTOP", allsig="yes",  interactions="no", lmDiAmodel_DROPBOX ))
lmAll <- rbind(lmAll, data.frame(method="Meth",    model="lmTOP", allsig="yes",  interactions="no", lmMeTHmodel_TOP_DROPBOX ))



lmAll <- rbind(lmAll, data.frame(method="TOPSIS", model="lmMeth", allsig="no",  interactions="yes", lmTOPmodel_MeTH_DROPBOX ))
lmAll <- rbind(lmAll, data.frame(method="DiA", model="lmMeth", allsig="no",  interactions="yes", lmDiAmodel_MeTH_DROPBOX ))
lmAll <- rbind(lmAll, data.frame(method="Meth", model="lmMeth", allsig="yes",  interactions="yes", lmMeTHmodel_DROPBOX ))

print(lmAll)

