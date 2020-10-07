######################################################## Loading Libraries #################################################
library(lmerTest)
library(lme4)  
library(nlme)
library(plyr)
library(influence.ME)  
library(Matrix) 
library(rstan)
library(rstanarm)
library(brms)
library(rlang)

############################################################################################################################
########################### Hierarchical Linear Modeling (HLM) for Peng et al. 2018 ########################################
############################################################################################################################

##################################################### Loading Data #########################################################
newdata <- read.csv ('SLR_Force3_newdata_healthy.csv') 

################################################# Data Pre-processing ######################################################
newdata=newdata[1:654,] #make sure there are no extra empty raws
#convert Sex from a integer to a factor
class(newdata$Sex)
newdata$Sex=as.factor(newdata$Sex)
class(newdata$Sex)
#change Sex names from 1 and 2 to female and male
levels(newdata$Sex)
levels(newdata$Sex)=c("Female","Male")
newdata$Sex
#convert Muscle from a integer to a factor 
class(newdata$Muscle)
newdata$Muscle=as.factor(newdata$Muscle)
class(newdata$Muscle)
#change Muscle names from 1 and 2 to vmo and vm
str(newdata$Muscle)
levels(newdata$Muscle)=c("VMO","VM")
newdata$Muscle
#convert Task from a integer to a factor 
class(newdata$Task)
newdata$Task=as.factor(newdata$Task)
class(newdata$Task)
#change Task names from 2 and 3 to NR and LR
str(newdata$Task)
levels(newdata$Task)=c("Neutral Rotation","Lateral Rotation")
newdata$Task
#check data
str(newdata)
dim(newdata)

########################################## HLM for Initial Firing Rate (IFR) ###############################################
m3 <- lmer (IFR ~ Sex*Task*Muscle + Absolute_force + (1 | Subject), REML=F, data=newdata)
summary(m3) 
anova(m3)

############################################## Post-hoc Analysis for IFR ###################################################
lsmeans(m3, pairwise ~ Sex, adjust="tukey")
lsmeans(m3, pairwise ~ Task, adjust="tukey")
lsmeans(m3, pairwise ~ Task*Muscle, adjust="tukey")

############################################### Assumption Tests for IRF ###################################################
m3 <- lme(IFR ~ Sex*Task*Muscle + Absolute_force, random=~1 | Subject, data= newdata)

#Normality of residuals within each level:
hist(residuals(m3,level= 0))
hist(residuals(m3,level= 1))

#Linearity/Equal Variance:
plot(fitted(m3), residuals(m3,level=0))
scatter.smooth(fitted(m3), residuals(m3, level=0))
abline(h=0,col='red')

plot(fitted(m3), residuals(m3,level=1))
scatter.smooth(fitted(m3), residuals(m3, level=1))
abline(h=0,col='red')

############################################ Figure 2 for IFR in Peng et al. 2018 ###########################################
par(mfrow=c(1,2), mai = c(1, 1.2, 0.8, 0.1))

boxplot(newdata$IFR ~ newdata$Sex, ylab="Initail Firing Rate (Hz)")
stripchart(newdata$IFR ~ newdata$Sex, vertical=TRUE, add=TRUE, pch=21) #add data points
abline(h=9.639668, lty=2, lwd=2) #female with adjusted mean
abline(h=8.460115, lty=3, lwd=2) #male with adjusted mean
legend(0.8,25, c("Model Estimate for Female", "Model Estimate for Male"), lty=c(2,3), lwd=1.2, bty="n", cex=1)

boxplot(newdata$IFR ~ newdata$Task)
stripchart(newdata$IFR ~ newdata$Task, vertical=TRUE, add=TRUE, pch=21) #add data points
abline(h=9.285440, lty=2, lwd=2) #NR with adjusted mean
abline(h=8.814342, lty=3, lwd=2) #LR with adjusted mean
legend(0.8,25, c("Model Estimate for Neutral Rotation", "Model Estimate for Lateral Rotation"), lty=c(2,3), lwd=1.2, bty="n", cex=1)

############################################ Figure 3 for IFR in Peng et al. 2018 ###########################################
#set a new column for IFR_Inter (Muscle*Task):
for (s in 1:654){
  if(newdata$Muscle[s]=="VMO" & newdata$Task[s]=="Neutral Rotation"){
    newdata$IFR_Inter[s]="VMO-Neutral Rotation"}
  else if (newdata$Muscle[s]=="VMO" & newdata$Task[s]=="Lateral Rotation"){
    newdata$IFR_Inter[s]="VMO-Lateral Rotation"}
  else if (newdata$Muscle[s]=="VM" & newdata$Task[s]=="Neutral Rotation"){
    newdata$IFR_Inter[s]="VM-Neutral Rotation"}
  else {newdata$IFR_Inter[s]="VM-Lateral Rotation"}
}

boxplot(newdata$IFR ~ newdata$IFR_Inter, ylim = c(2,25), ylab="Initail Firing Rate (Hz)")
stripchart(newdata$IFR ~ newdata$IFR_Inter, vertical=TRUE, add=TRUE, pch=21) #add data points
abline(h=9.559731, lty=2, lwd=1.2) #VMO-NR with adjusted mean
abline(h=8.758814, lty=4, lwd=1.2) #VMO-LR with adjusted mean
abline(h=9.011150, lty=3, lwd=1.2) #VM-NR with adjusted mean
abline(h=8.869870, lty=1, lwd=1.2) #VM-LR with adjusted mean
legend(1,25, c("Model Estimate for VMO-Neutral Rotation", "Model Estimate for VMO-Lateral Rotation", "Model Estimate for VM-Neutral Rotation", "Model Estimate for VM-Lateral Rotation"), lty=c(2,4,3,1), lwd=1.2, bty="n", cex=1)

########################################## HLM for Recruitment Threshold (RT) ###############################################
newdata$Log_RT1=log(newdata$RT+1) #log-transformed the data to meet model assumption
m7 <- lmer (Log_RT1 ~ Sex*Task*Muscle + (1 | Subject), REML=F, data= newdata)
summary(m7) 
anova(m7)

############################################### Post-hoc Analysis for RT ####################################################
lsmeans(m7, pairwise ~ Muscle, adjust="tukey")
lsmeans(m7, pairwise ~ Task, adjust="tukey")
lsmeans(m7, pairwise ~ Sex*Task, adjust="tukey")

################################################ Assumption Tests for RT ####################################################
m7 <- lme (Log_RT1 ~ Sex*Task*Muscle, random=~1 | Subject, data= newdata) 

#Normality of residuals within each level:
hist(residuals(m7,level= 0))
hist(residuals(m7,level= 1))

#Linearity/Equal Variance:
plot(fitted(m7), residuals(m7,level=0))
scatter.smooth(fitted(m7), residuals(m7, level=0))
abline(h=0,col='red')

plot(fitted(m7), residuals(m7,level=1))
scatter.smooth(fitted(m7), residuals(m7, level=1))
abline(h=0,col='red')

############################################# HLM with Raw RT for Figures ###################################################
m7_raw <- lmer (RT ~ Sex*Task*Muscle + (1 | Subject), REML=F, data= newdata)
summary(m7_raw) 
anova(m7_raw)
#post-hoc Analysis for raw RT
lsmeans(m7_raw, pairwise ~ Muscle, adjust="tukey")
lsmeans(m7_raw, pairwise ~ Task, adjust="tukey")
lsmeans(m7_raw, pairwise ~ Sex*Task, adjust="tukey")

############################################ Figure 4 for IFR in Peng et al. 2018 ###########################################
newdata$RTx100=newdata$RT*100 #change into %MVC

#par(mfrow=c(2,2), mai = c(0.5, 1.2, 0.4, 0.1))
par(mfrow=c(1,2), mai = c(0.8, 1.5, 0.5, 0.2))

boxplot(newdata$RTx100 ~ newdata$Muscle, ylab="Recruitment Threshold (%MVC)", ylim=c(0,80))
stripchart(newdata$RTx100 ~ newdata$Muscle, vertical=TRUE, add=TRUE, pch=21) #add data points
abline(h=11.89351, lty=2, lwd=2) #VMO (exponentiated values)
abline(h=14.81019, lty=3, lwd=2) #VM (exponentiated values)
legend(0.8, 83, c("Model Estimate for VMO", "Model Estimate for VM"), lty=c(2,3), lwd=1.2, bty="n", cex=1)

boxplot(newdata$RTx100 ~ newdata$Task, ylim=c(0,80))
stripchart(newdata$RTx100 ~ newdata$Task, vertical=TRUE, add=TRUE, pch=21) #add data points
means=tapply(newdata$RTx100, newdata$Task, mean)
abline(h=11.48042, lty=2, lwd=2) #NR (exponentiated values)
abline(h=15.22327, lty=3, lwd=2) #LR (exponentiated values)
legend(0.7, 83, c("Model Estimate for Neutral Rotation", "Model Estimate for Lateral Rotation"), lty=c(2,3), lwd=1.2, bty="n", cex=1)

############################################ Figure 5 for IFR in Peng et al. 2018 ###########################################
#set a new column for RT_Inter (Sex*Task):
for (s in 1:654){
  if(newdata$Sex[s]=="Female" & newdata$Task[s]=="Neutral Rotation"){
    newdata$RT_Inter[s]="Female-Neutral Rotation"}
  else if (newdata$Sex[s]=="Female" & newdata$Task[s]=="Lateral Rotation"){
    newdata$RT_Inter[s]="Female-Lateral Rotation"}
  else if (newdata$Sex[s]=="Male" & newdata$Task[s]=="Neutral Rotation"){
    newdata$RT_Inter[s]="Male-Neutral Rotation"}
  else {newdata$RT_Inter[s]="Male-Lateral Rotation"}
}

boxplot(newdata$RTx100 ~ newdata$RT_Inter, ylim = c(2,80), ylab="Recruitment Threshold (%MVC)")
stripchart(newdata$RTx100 ~ newdata$RT_Inter, vertical=TRUE, add=TRUE, pch=21) #add data points
abline(h=11.37826, lty=1, lwd=1.2) #Female-NR with adjusted mean
abline(h=17.65401, lty=2, lwd=1.2) #Female-LR with adjusted mean
abline(h=11.58258, lty=3, lwd=1.2) #Male-NR with adjusted mean
abline(h=12.79254, lty=4, lwd=1.2) #Male-LR with adjusted mean
legend(1.5, 80, c("Model Estimate for Female-Neutral Rotation", "Model Estimate for Female-Lateral Rotation", "Model Estimate for Male-Neutral Rotation", "Model Estimate for Male-Lateral Rotation"), lty=c(1,2,3,4),lwd=1.2, bty="n", cex=1)

############################################################################################################################
################################# Bayesian Inference for Publication in Preparation ########################################
############################################################################################################################

##################################################### Loading Data #########################################################
newdata <- read.csv ('SLR_Force3_newdata_PFPS10.csv')

##################################################### Data Pre-processing ##################################################
newdata=newdata[1:579,] #make sure there are no extra empty raws
#convert Muscle from a integer to a factor 
class(newdata$Muscle)
newdata$Muscle=as.factor(newdata$Muscle)
class(newdata$Muscle)
#change Muscle names from 1 and 2 to vmo and vm
str(newdata$Muscle)
levels(newdata$Muscle)=c("VMO","VM")
newdata$Muscle
#convert Task from a integer to a factor 
class(newdata$Task)
newdata$Task=as.factor(newdata$Task)
class(newdata$Task)
#change Task names from 2 and 3 to NR and LR
str(newdata$Task)
levels(newdata$Task)=c("Neutral Rotation","Lateral Rotation")
newdata$Task
#convert Group from a integer to a factor 
class(newdata$Group)
newdata$Group=as.factor(newdata$Group)
class(newdata$Group)
#change Group names from 0 and 1 to healthy and PFPS
str(newdata$Group)
levels(newdata$Group)=c("Healthy","PFPS")
newdata$Group
#check data
str(newdata)
dim(newdata)

##################################### Bayesian Multilevel Modeling with brms for IFR #######################################
fit_IFR <- brm(IFR ~ Group*Task*Muscle + Absolute_force + (1 | Subject), data=newdata)
summary(fit_IFR)
posterior_summary(fit_IFR)
plot(fit_IFR)
marginal_effects(fit_IFR) #visualize effects of predictors on the expected response
pp_check(fit_IFR, nsamples = 100) #posterior probability checking

##################################### Bayesian Multilevel Modeling with brms for RT ########################################
newdata$Log_RT1=log(newdata$RT+1) #log-transformed the data to meet multilevel assumptions
fit_RT <- brm(Log_RT1 ~ Group*Task*Muscle + (1 | Subject), data= newdata)
summary(fit_RT)
posterior_summary(fit_RT)
plot(fit_RT)
marginal_effects(fit_RT) #visualize effects of predictors on the expected response
pp_check(fit_RT, nsamples = 100) #posterior probability checking

################################################## Violin Plots for RT #####################################################
newdata$RTx100=newdata$RT*100 #change into %MVC

#select data in column Task and RTx100 for "Main effect of Task (Figure 2)"
F2data <- newdata[,c(8,21)]  
colnames(F2data) <- c("HipRotation", "RecruitmentThreshold")  #change name of the heading

ggplot(F2data, aes(x = HipRotation, y = RecruitmentThreshold)) + 
  geom_violin() + stat_summary(fun.y='median' ,geom='point', shape= 3) +
  geom_point(aes(x = 1, y = 13.13), shape = 19) +
  geom_errorbar(aes(x =1, ymax = 14.32, ymin = 11.94), width=.03) +
  geom_point(aes(x = 2, y = 16.25), shape = 19) +
  geom_errorbar(aes(x =2, ymax = 17.442, ymin = 15.06), width=.03) +
  theme_bw() +
  labs(x="Hip Position", y = "Recruitment Threshold (%MVC)") +
  theme(text = element_text(size=20))

#select data in column Group, Task and RTx100 for " Interaction of Group*Task (Figure 3)"
F3data <- newdata[,c(5, 8, 21)] 

ggplot(F3data, aes(x = Task, y = RTx100, fill = Group)) + 
  geom_violin() + stat_summary(fun.y='median' ,geom='point', shape= 3, position = position_dodge(width = .9), show.legend = F) +
  geom_point(aes(x = 0.775, y = 10.3628028), shape = 19, show.legend = F) +
  geom_errorbar(aes(x =0.775, ymax = 11.9984018, ymin = 8.7272038), width=.03) +
  geom_point(aes(x = 1.225, y = 15.9624099), shape = 19, show.legend = F) +
  geom_errorbar(aes(x = 1.225, ymax = 17.6919239, ymin = 14.2328959), width=.03) +
  geom_point(aes(x = 1.775, y = 16.11865265), shape = 19, show.legend = F) +
  geom_errorbar(aes(x = 1.775, ymax = 17.74539965, ymin = 14.49190565), width=.03) +
  geom_point(aes(x = 2.225, y = 16.373059), shape = 19, show.legend = F) +
  geom_errorbar(aes(x = 2.225, ymax = 18.120377, ymin = 14.625741), width=.03) +
  scale_fill_manual(values=c("white", "grey")) +
  theme_bw() +
  labs(x="Hip Position", y = "Recruitment Threshold (%MVC)") +
  theme(text = element_text(size=20))

#select data in column Group, Muscle and RTx100 for " Interaction of Group*Muscle (Figure 4)"
F4data <- newdata[,c(5, 7, 21)] 

ggplot(F4data, aes(x = Muscle, y = RTx100, fill = Group)) + 
  geom_violin() + stat_summary(fun.y='median' ,geom='point', shape= 3, position = position_dodge(width = .9), show.legend = F) +
  geom_point(aes(x = 0.775, y = 11.80694582), shape = 19, show.legend = F) +
  geom_errorbar(aes(x =0.775, ymax = 13.38218182, ymin = 10.23170982), width=.03) +
  geom_point(aes(x = 1.225, y = 17.13712839), shape = 19, show.legend = F) +
  geom_errorbar(aes(x = 1.225, ymax = 18.86805839, ymin = 15.40619839), width=.03) +
  geom_point(aes(x = 1.775, y = 14.61881677), shape = 19, show.legend = F) +
  geom_errorbar(aes(x = 1.775, ymax = 16.31061577, ymin = 12.92701777), width=.03) +
  geom_point(aes(x = 2.225, y = 15.2059949), shape = 19, show.legend = F) +
  geom_errorbar(aes(x = 2.225, ymax = 16.9552189, ymin = 13.4567709), width=.03) +
  scale_fill_manual(values=c("white", "grey")) +
  theme_bw() +
  labs(x="Muscle", y = "Recruitment Threshold (%MVC)") +
  theme(text = element_text(size=20))


