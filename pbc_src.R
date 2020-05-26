rm(list = ls())
gc()
gc()

library(data.table)
library(survival)
library(survminer)
library(tidyverse)
library(glmnet)
library(e1071)
library(coin)

'%nin%' = Negate('%in%')

#importing 'pbc' dataset to R environment
  
data(pbc)

pbc = as.data.table(pbc)



# Descriptive Statistics

numeric_vars = c("age","bili", "chol", "albumin", "copper", "alk.phos", 
                 "ast", "platelet", "protime")

factor_vars = c("trt", "sex", "ascites", "hepato", 
                "spiders", "edema","stage")



# NA count

sapply(pbc,function(x) sum(is.na(x)))


#Analysis of Surv Obj and target variable

hist(pbc$time)
summary(pbc$time)
table(pbc$status)

Surv_obj = Surv(pbc$time,pbc$status==2)

ggsurvplot(survfit(Surv_obj~1,data = pbc),data = pbc,
           palette = "blue",linetype = 1,
           surv.median.line = "hv",risk.table = TRUE,
           cumevents = TRUE, cumcensor = TRUE, tables.height = 0.15)

survfit(Surv(time,status==2)~1,data=pbc)

# -- Num var


summary(pbc$age)

summary(pbc$bili)
summary(pbc$chol)
summary(pbc$albumin)
summary(pbc$copper)

summary(pbc$alk.phos)
summary(pbc$ast)
summary(pbc$platelet)
summary(pbc$protime)




# -- cat var

table(pbc$sex,exclude = NULL)
prop.table(table(pbc$sex,exclude = NULL))

table(pbc$stage,exclude = NULL)
prop.table(table(pbc$stage,exclude = NULL))

table(pbc$ascites,exclude = NULL)
prop.table(table(pbc$ascites,exclude = NULL))

table(pbc$hepato,exclude = NULL)
prop.table(table(pbc$hepato,exclude = NULL))

table(pbc$spiders,exclude = NULL)
prop.table(table(pbc$spiders,exclude = NULL))

table(pbc$edema,exclude = NULL)
prop.table(table(pbc$edema,exclude = NULL))


table(pbc$stage,exclude = NULL)
prop.table(table(pbc$edema,exclude = NULL))




# -- data_cleaning

pbc  = pbc[,trt:=ifelse(is.na(trt),3,trt)]

pbc = pbc[,ascites:=ifelse(is.na(ascites),-1,ascites)]

pbc = pbc[,hepato:=ifelse(is.na(hepato),-1,hepato)]

pbc = pbc[,spiders:=ifelse(is.na(spiders),-1,spiders)]

pbc[,stage:=as.character(stage)]
pbc[,stage:=ifelse(is.na(stage),'missing',stage)]


# -- num var

copper_median = median(pbc$copper,na.rm=TRUE)
pbc = pbc[,copper:=ifelse(is.na(copper),copper_median,copper)]

alkphos_median = median(pbc$alk.phos,na.rm=TRUE)
pbc = pbc[,alk.phos:=ifelse(is.na(alk.phos),alkphos_median,alk.phos)]

ast_median = median(pbc$ast,na.rm=TRUE)
pbc = pbc[,ast:=ifelse(is.na(ast),ast_median,ast)]

platelet_median = median(pbc$platelet,na.rm=TRUE)
pbc = pbc[,platelet:=ifelse(is.na(platelet),platelet_median,platelet)]

protime_median = median(pbc$protime,na.rm=TRUE)
pbc = pbc[,protime:=ifelse(is.na(protime),protime_median,protime)]

# -- rejecting trig and chol

pbc_old = pbc
pbc[,trig:=NULL]
pbc[,chol:=NULL]

# -- tranforming to log

numeric_vars = c("age","bili", "albumin", "copper", "alk.phos", 
                 "ast", "platelet", "protime")

sapply(pbc[,.SD,.SDcols=numeric_vars],function(x) skewness(x))

#pbc = pbc[,protime_log:=log(protime)]
#pbc = pbc[,ast_log:=log(ast)]

#age_mean = mean(pbc$age)
#age_sd = sd(pbc$age)
#pbc = pbc[,age_norm:=(age-age_mean)/age_sd]

pbc[,age_c:=ifelse(age<=44,"<=44",
                                  ifelse(age<=52,"<=52",
                                         ifelse(age<=59,"<=59",">59")))]

pbc[,stage:=ifelse(is.na(stage),'missing',stage)]

pbc[,stage_4_flag:=ifelse(stage==4,1,0)]



#non-Parametric

km_pbc = survfit(Surv(time,status==2)~1,data = pbc)

km_pbc = survfit(Surv_obj~1,data = pbc)

summary(survfit(Surv_obj ~ 1, data = pbc), times = 365.25)
summary(survfit(Surv_obj ~ 1, data = pbc), times = 730.5)
summary(survfit(Surv_obj ~ 1, data = pbc), times = 2000)


km_pbc_2 =  survfit(Surv_obj~age_c,data = pbc)

km_pbc_3 =  survfit(Surv_obj~stage,data = pbc)


ggsurvplot(km_pbc,data = pbc,
           palette = "blue",linetype = 1,
           surv.median.line = "hv",risk.table = TRUE,
           cumevents = TRUE, cumcensor = TRUE, tables.height = 0.15)

km_model = survreg(Surv(time,status==2)~as.factor(trt) +	as.factor(sex) +	as.factor(ascites) +	
                     as.factor(hepato) +	as.factor(spiders) +	as.factor(edema) +	as.factor(stage) +
                     age + 	log(bili) + 	log(chol) + albumin + 	log(copper) + 	log(alk.phos) + 	log(ast) + 	log(platelet) + 	
                     log(protime),data = pbc,dist = "weibull", scale = 0.5)

km_model_2 = survreg(Surv(time,status==2)~	as.factor(edema) +
                     age + 	log(bili)  + albumin + log(protime),data = pbc,dist = "weibull", scale = 0.5)

#log_rank_test

logrank_test(Surv(time,status==2) ~ sex, data = pbc)
logrank_1 <- survdiff(Surv(time, status == 2) ~ sex, data = pbc)


survdiff(Surv(time, status == 2) ~ as.factor(sex), data = pbc)
survdiff(Surv(time, status == 2) ~ as.factor(age_c), data = pbc)
survdiff(Surv(time, status == 2) ~ as.factor(stage), data = pbc)
survdiff(Surv(time, status == 2) ~ as.factor(edema), data = pbc)

plot(KM_fit, fun = "cumhaz")

ggsurvplot(KM_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups   
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

logrank_test(Surv(time,status==2) ~ as.factor(age_c), data = pbc)
KM_fit <- survfit(Surv(time, status == 2) ~ as.factor(age_c), data = pbc)
plot(KM_fit, fun = "cumhaz")

ggsurvplot(KM_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF","#009999", "#0000FF"))

logrank_test(Surv(time,status==2) ~ as.factor(edema), data = pbc)
KM_fit <- survfit(Surv(time, status == 2) ~ edema, data = pbc)
plot(KM_fit, fun = "cumhaz")

ggsurvplot(KM_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF","#009999"))



logrank_test(Surv(time,status==2) ~ as.factor(trt), data = pbc)

logrank_test(Surv(time,status==2) ~ as.factor(ascites), data = pbc)

logrank_test(Surv(time,status==2) ~ as.factor(hepato), data = pbc)

logrank_test(Surv(time,status==2) ~ as.factor(spiders), data = pbc)



#cox_ph

## The 4 candidate models

M0 <- coxph(Surv(time,status==2) ~ 1, data = pbc)


MA <- coxph(Surv(time,status==2) ~ age, data = pbc)

summary(MA)


MB <- coxph(Surv(time,status==2) ~ age+sex, data = pbc)

summary(MB)


anova(MA,MB)

MC <- coxph(Surv(time,status==2) ~ edema, data = pbc)
MC1 <- coxph(Surv(time,status==2) ~ edema+age+spiders, data = pbc)
MC2 <- coxph(Surv(time,status==2) ~ edema+age+sex, data = pbc)
MD <- coxph(Surv(time,status==2) ~ age + sex+ edema, data = pbc)
ME <- coxph(Surv(time,status==2) ~ as.factor(trt) +	as.factor(sex) +	as.factor(ascites) +	
              as.factor(hepato) +	as.factor(spiders) +	as.factor(edema) +	as.factor(stage) +
              age + 	bili + 	chol + 	albumin + 	copper + 	alk.phos + 	ast + 	platelet + 	protime,data=pbc)



fits <- list(MA = MA, MB = MB, MC = MC, MD = MD,ME = ME, MF = MF)
sapply(fits, AIC)

MF <- coxph(Surv(time,status==2) ~ as.factor(trt) +	as.factor(sex) +	as.factor(ascites) +	
              as.factor(hepato) +	as.factor(spiders) +	as.factor(edema) +	as.factor(stage) +
              age + 	log(bili) + albumin + 	log(copper) + 	log(alk.phos) + 	log(ast) + 	log(platelet) + 	
              log(protime),data=pbc)

summary(MF)

MF_1<-coxph(Surv(time,status==2) ~ 	as.factor(edema)  +
              age + 	log(bili) + albumin + log(protime),data=pbc)

summary(MF_1)

sapply(list(MF,MF_1),AIC)

par(mfrow=c(2, 3))
plot(cox.zph(MF_1))

plot(cox.zph(MC))

plot(cox.zph(MA))






# cox glm net





attach(pbc)
par(mfrow=c(1, 1))
y <- Surv(pbc$time, pbc$status==2)

x <- model.matrix(y ~ as.factor(trt) +	as.factor(sex) +	as.factor(ascites) +	
                     as.factor(hepato) +	as.factor(spiders) +	as.factor(edema) +	as.factor(stage) +
                     age + 	log(bili) + albumin + 	log(copper) + 	log(alk.phos) + 	log(ast) + 	log(platelet) + 	
                     log(protime)-1,pbc)

cv.fit <- cv.glmnet(x,y,family="cox", maxit = 1000)

plot(cv.fit)


cv.fit$lambda.min   # CV-error curve hits its minimum


cv.fit$lambda.1se   # most regularized model with CV-error within 1 standard deviation of the minimum.

coef(cv.fit, s = "lambda.min")

x_req = model.matrix(y ~ as.factor(ascites)+as.factor(edema)+as.factor(stage)+
                       age + 	log(bili) + albumin + 	log(copper) + 	log(ast) + 	
                       log(protime)-1, pbc)

glmnet_fit <- glmnet(x,y,family="cox",maxit = 1000,lambda = 0)
coef(glmnet_fit)

glmnet_fit1 <- glmnet(x,y,family="cox",maxit = 1000,lambda = cv.fit$lambda.min)
coef(glmnet_fit1)

glmnet_fit2 <- glmnet(x_req,y,family="cox",maxit = 1000,lambda = cv.fit$lambda.min)

par(mfrow=c(3, 3))
plot(glmnet_fit2)

Coefficients <- coef(glmnet_fit, s = cv.fit$lambda.min)

coef(glmnet_fit, s = cv.fit$lambda.min)

coef(glmnet_fit1, s = cv.fit$lambda.min)
