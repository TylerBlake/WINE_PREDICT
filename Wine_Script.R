library(tidyverse)
library(DataExplorer)
library(stargazer)
WINEI <- read.csv(file = "wine-training-data.csv", header = TRUE, sep = ",", na.strings = c("", " ", "NA"))
WINE <- select(WINEI, -ï..INDEX)


###DATA EXPLORATION##

str(WINE)

stargazer(WINE, type = "html", out = "WINE_SUM.doc")
plot_missing(WINE)
plot_histogram(WINE, geom_histogram_args = list(fill = "cornflowerblue"))
plot_correlation(WINE, cor_args = list("use" = "pairwise.complete.obs"))

##DATA PREP##

##Negative Values##
WINE_COMPLETE <- WINE %>%
  mutate(Alcohol = ifelse(Alcohol < 0, 0, ifelse(Alcohol > 24, 24, Alcohol)),
         LabelAppeal = as.factor(LabelAppeal),
         FixedAcidity = ifelse(FixedAcidity < 0, FixedAcidity * (-1), FixedAcidity),
         VolatileAcidity = ifelse(VolatileAcidity < 0, VolatileAcidity * (-1), VolatileAcidity),
         CitricAcid = ifelse(CitricAcid < 0, CitricAcid * (-1), CitricAcid),
         ResidualSugar = ifelse(ResidualSugar < 0, ResidualSugar * (-1), ResidualSugar),
         Chlorides = ifelse(Chlorides < 0, Chlorides * (-1), Chlorides),
         FreeSulfurDioxide = ifelse(FreeSulfurDioxide < 0, FreeSulfurDioxide * (-1), FreeSulfurDioxide),
         TotalSulfurDioxide = ifelse(TotalSulfurDioxide < 0, TotalSulfurDioxide * (-1), TotalSulfurDioxide),
         Sulphates = ifelse(Sulphates < 0, Sulphates * (-1), Sulphates))
##Missing Values##
library(mice)
WINE_COMPLETE <- WINE_COMPLETE %>%
  mutate(ALC_MISS = as.factor(ifelse(is.na(Alcohol), 1, 0)),
         TSD_MISS = as.factor(ifelse(is.na(TotalSulfurDioxide), 1, 0)),
         SUL_MISS = as.factor(ifelse(is.na(Sulphates), 1, 0)),
         FSD_MISS = as.factor(ifelse(is.na(FreeSulfurDioxide), 1, 0)),
         CHL_MISS = as.factor(ifelse(is.na(Chlorides), 1, 0)),
         RES_MISS = as.factor(ifelse(is.na(ResidualSugar), 1, 0)),
         pH_MISS = as.factor(ifelse(is.na(pH), 1, 0)),
         STARS = as.factor(ifelse(is.na(STARS), "No Rating", STARS)))
#Create predictor matrix for mice excluding flags and TARGET for prdiction     
MISSING_VARS <- c("pH", "ResidualSugar", "Chlorides", "FreeSulfurDioxide", "Alcohol", "TotalSulfurDioxide", "Sulphates")
EXCL_PRED <- c("ALC_MISS", "TSD_MISS", "SUL_MISS", "FSD_MISS", "CHL_MISS", "RES_MISS", "pH_MISS", "TARGET")
allVARS <- names(WINE_COMPLETE)
PREDMAT <- matrix(0, ncol = length(allVARS), nrow = length(allVARS))
rownames(PREDMAT) <- allVARS
colnames(PREDMAT) <- allVARS
PREDMAT[MISSING_VARS,] <- 1
PREDMAT[,"TARGET"] <- 0
diag(PREDMAT) <- 0
PREDMAT["Alcohol", "ALC_MISS"] <- 0
PREDMAT[ "TotalSulfurDioxide", "TSD_MISS"] <- 0
PREDMAT["FreeSulfurDioxide", "FSD_MISS"] <- 0
PREDMAT["Sulphates", "SUL_MISS"] <- 0
PREDMAT["Chlorides", "CHL_MISS"] <- 0
PREDMAT["ResidualSugar", "RES_MISS"] <- 0
PREDMAT["pH", "pH_MISS"] <- 0
#MICE IMPUTE
WINE_IMPUTE <- mice(WINE_COMPLETE, predictorMatrix = PREDMAT)
#REVIEW IMPUTES
densityplot(WINE_IMPUTE)
#Complete Data Set
head(WINE_IMPUTE)
WINE_IMP_CMP <- complete(WINE_IMPUTE, 5)
summary(WINE_IMP_CMP)

###MODEL BUILDING###
library(MASS)
library(pscl)

##Poisson Model 1###
MODEL_1_POISSON <- glm(TARGET ~ ., family = "poisson", data = WINE_IMP_CMP)
POISSON_BACKWARD <- stepAIC(MODEL_1_POISSON, direction = "backward")
options(scipen = 999)
summary(POISSON_BACKWARD)
anova(POISSON_BACKWARD, test = "Chisq")
WINE_POISSON <- WINE_IMP_CMP %>%
  mutate(pHunscaled = 10^(-pH),
         Impure = log(Chlorides + Sulphates),
         STARS = relevel(STARS, ref = 5))
MODEL_1a_POISSON <- glm(formula = TARGET ~ VolatileAcidity + Impure + TotalSulfurDioxide + 
                          Alcohol + LabelAppeal + AcidIndex + 
                          STARS, family = poisson(link = "log"), data = WINE_POISSON)
summary(MODEL_1a_POISSON)




Pred1 <- predict(MODEL_1a_POISSON, WINE_POISSON, type = "response")


pcount1 <- colSums(predprob(MODEL_1a_POISSON)[,1:9])
hist(WINE_POISSON$TARGET, col = "cornflowerblue", breaks=seq(-.5, 8.5, 1), xlab = "Wine Cases",main = NULL)
lines(x = seq(0, 8, 1), y = pcount1, lty=2, lwd=2, col="red")
points(x = seq(0, 8, 1), y = pcount1, cex=2, pch=16, col="red")

deviance(MODEL_1a_POISSON)/df.residual(MODEL_1a_POISSON)
anova(ssr, MODEL_1a_POISSON$df.residual, test = "Chisq")
ssr <- sum(residuals(MODEL_1a_POISSON, type = "pearson")^2)
pchisq(ssr, df = MODEL_1a_POISSON$df.residual, lower.tail = FALSE)

stargazer(MODEL_1a_POISSON, type = "html", out = "Poisson_Model.doc")
mean(var(WINE_POISSON$TARGET))
SIMRES <- simulateResiduals(MODEL_1a_POISSON, n = 1000, plot = T)
testZeroInflation(SIMRES)

##Zero-Inflated Poisson Model##
Zero_Poisson1 <- zeroinfl(formula = TARGET ~ VolatileAcidity + Impure + TotalSulfurDioxide + Alcohol + LabelAppeal + AcidIndex + 
                            STARS|VolatileAcidity + Impure + TotalSulfurDioxide + Density + Alcohol + LabelAppeal + AcidIndex + STARS, link = "logit", dist = "poisson", data = WINE_POISSON)
summary(Zero_Poisson1)
##
Zero_Poisson2 <- zeroinfl(formula = TARGET ~ Alcohol + LabelAppeal + AcidIndex + 
  STARS|VolatileAcidity + Impure + TotalSulfurDioxide + Alcohol + LabelAppeal + AcidIndex + STARS, link = "logit", dist = "poisson", data = WINE_POISSON)

WINE_ZIP <- WINE_POISSON %>%
 mutate(Rating = as.factor(ifelse(STARS == 1 | STARS == 2, "Poor", ifelse(STARS == "No Rating", "None", "High"))),
        PoorLAB = as.factor(ifelse(LabelAppeal == -1 | LabelAppeal == -2, 1, 0)),
        TYP_pH = as.factor(ifelse(pH >= 3 & pH < 3.9, 0, 1)),
        TYP_ALC = as.factor(ifelse(Alcohol >= 10 & Alcohol <= 16, 0, 1)),
        SWEET = as.factor(ifelse(ResidualSugar < 10, "Dry", ifelse(ResidualSugar >= 10 & ResidualSugar <= 30, "Off-Dry",
                                                                   ifelse(ResidualSugar > 30 & ResidualSugar < 100, "Sweet", "Dessert")))),
        LAB_Rating = as.factor(ifelse(LabelAppeal == -1 | LabelAppeal == -2, "Unappealing", ifelse(LabelAppeal == 0, "Indifferent", "Appealing"))),
        NoSTARS = as.factor(ifelse(STARS == "No Rating", 1, 0)),
        SALT_RAT = log((ResidualSugar + 1)/(Chlorides + 1)),
        TYP_SWEET = as.factor(ifelse(ResidualSugar >= .2 & ResidualSugar <= 30, 1, 0)))
Zero_Poisson3 <- zeroinfl(formula = TARGET ~ Alcohol + LabelAppeal + AcidIndex + 
                            STARS|VolatileAcidity + Impure + TotalSulfurDioxide + Rating + AcidIndex + pH + LAB_Rating + SALT_RAT, link = "logit", dist = "poisson", data = WINE_ZIP) 
summary(Zero_Poisson3)
Zero_Poisson4 <- zeroinfl(formula = TARGET ~ Alcohol + LabelAppeal + AcidIndex + 
                            STARS + TYP_ALC|VolatileAcidity + Impure + TotalSulfurDioxide + Alcohol + LabelAppeal + AcidIndex + STARS + pH, link = "logit", dist = "poisson", data = WINE_ZIP)
summary(Zero_Poisson4)
Zero_Poisson5 <- zeroinfl(formula = TARGET ~ Alcohol + LabelAppeal + AcidIndex + 
                            STARS + TYP_ALC|VolatileAcidity + Impure + TotalSulfurDioxide + Alcohol + LabelAppeal + AcidIndex + STARS + pH + SALT_RAT + ResidualSugar, link = "logit", dist = "poisson", data = WINE_ZIP)

lrtest(Zero_Poisson5, Zero_Poisson4)
vuong(MODEL_1a_POISSON, Zero_Poisson5)
pcount2 <- colSums(predprob(Zero_Poisson5)[,1:9])
hist(WINE_POISSON$TARGET, col = "cornflowerblue", breaks=seq(-.5, 8.5, 1), xlab = "Wine Cases",main = NULL)
lines(x = seq(0, 8, 1), y = pcount2, lty=2, lwd=2, col="red")
points(x = seq(0, 8, 1), y = pcount2, cex=2, pch=16, col="red")
stargazer(Zero_Poisson5, type = "html", out = "Zero-Inflated_PS_Model.doc")
stargazer(Zero_Poisson5, type = "html", out = "Zero-Inflated_PS2_Model.doc", zero.component = TRUE)
AIC(Zero_Poisson5)
AIC(NB_HURDLE)
summary(Zero_Poisson5)
ssr2 <- sum(residuals(Zero_Poisson1, type = "pearson")^2)
ssr2/df.residual(Zero_Poisson1)
Zero_Poisson5$residuals
rstd_dev <- sqrt(sum(Zero_Poisson5$residuals^2)/12794)
std_res <- Zero_Poisson5$residuals/rstd_dev
plot(WINE_ZIP$TARGET, std_res)
plot(fitted(Zero_Poisson5), std_res)
abline(0,1)
plot(WINE_ZIP$TARGET,std_res)
abline(0,0)
     ##Negative Binomial##
Wine_NB <- WINE_ZIP
residuals(Zero_Poisson5)
null_NB_Model <- glm.nb(TARGET ~ 1, Wine_NB)
full_NB_Model <- glm.nb(TARGET ~ ., Wine_NB)
summary(NB_Model1)
NB_Model2 <- stepAIC(null_NB_Model, scope = list(upper=full_NB_Model), direction = "both")
summary(NB_Model2)
NB_Model3 <- update(NB_Model2, . ~ . - Density)
options(scipen = 999)

NB_ZERO <- zeroinfl(TARGET ~ STARS + LabelAppeal + AcidIndex + VolatileAcidity + 
                      TYP_ALC + Impure + TotalSulfurDioxide + SALT_RAT + TYP_SWEET + 
                      pH|STARS + LabelAppeal + AcidIndex + VolatileAcidity + 
                      TYP_ALC + Impure + TotalSulfurDioxide + SALT_RAT + TYP_SWEET + 
                      pH, link = "logit", dist = "negbin", Wine_NB)

summary(NB_ZERO)
##Create pH and AcidIndex interaction term##
Wine_NB <- Wine_NB %>%
  mutate(Acid_Inter = pH * AcidIndex)

NB_ZERO2 <- zeroinfl(TARGET ~ STARS + LabelAppeal + AcidIndex + VolatileAcidity + 
                      TYP_ALC + Impure + TotalSulfurDioxide + SALT_RAT + TYP_SWEET + 
                      pH|STARS + LabelAppeal + AcidIndex + VolatileAcidity + 
                      TYP_ALC + Impure + TotalSulfurDioxide + SALT_RAT + TYP_SWEET + 
                      pH + Acid_Inter, link = "logit", dist = "negbin", Wine_NB)

summary(NB_ZERO2)
##Export Coefficients##
stargazer(NB_ZERO2, type = "html", out = "Zero-Inflated_NB_Model.doc")
stargazer(NB_ZERO2, type = "html", out = "Zero-Inflated_NB2_Model.doc", zero.component = TRUE)
##Plot Predictions##
pcount3 <- colSums(predprob(NB_ZERO2)[,1:9])
hist(WINE_NB$TARGET, col = "cornflowerblue", breaks=seq(-.5, 8.5, 1), xlab = "Wine Cases",main = NULL)
lines(x = seq(0, 8, 1), y = pcount3, lty=2, lwd=2, col="red")
points(x = seq(0, 8, 1), y = pcount3, cex=2, pch=16, col="red")


lmtest::lrtest(Zero_Poisson5, NB_ZERO2)
##NB2##
NB_HURDLE <- hurdle(TARGET ~ STARS + LabelAppeal + AcidIndex + TYP_ALC + SALT_RAT + pH + TYP_SWEET|STARS + LabelAppeal + AcidIndex + VolatileAcidity + 
                       TYP_ALC + Impure + TotalSulfurDioxide + SALT_RAT + TYP_SWEET + 
                       pH, dist = "negbin", Wine_NB)
stargazer(NB_HURDLE, type = "html", out = "HURDLE_NB_Model.doc")
stargazer(NB_HURDLE, type = "html", out = "HURDLE_NB2_Model.doc", zero.component = TRUE)

##Plot  Hurdle Predictions##
pcount4 <- colSums(predprob(NB_HURDLE)[,1:9])
hist(WINE_NB$TARGET, col = "cornflowerblue", breaks=seq(-.5, 8.5, 1), xlab = "Wine Cases",main = NULL)
lines(x = seq(0, 8, 1), y = pcount4, lty=2, lwd=2, col="red")
points(x = seq(0, 8, 1), y = pcount4, cex=2, pch=16, col="red")
vuong(Zero_Poisson5, NB_HURDLE)
vuong(NB_ZERO2, NB_HURDLE)
summary(NB_HURDLE)

plot(fitted(Zero_Poisson5), fitted(NB_HURDLE))
abline(0,1)
###LINEAR MODEL###
Linear_1 <- lm(TARGET ~ ., WINE_NB)
null_linear <- lm(TARGET ~ 1, WINE_NB)
Linear_2 <- stepAIC(null_linear, scope = list(upper=Linear_1), direction = "both")
summary(Linear_2)
Linear_3 <- update(Linear_2, . ~ . - RES_MISS - pHunscaled)
plot(Linear_3)
summary(Linear_3)
options(scipen = 999)
car::mmps(Linear_3)
par(mfrow=c(2,2))
plot(Linear_3)
LinFits <- fitted.values(Linear_3)
summary(LinFits)
durbinWatsonTest(Linear_3)
vif(Linear_3)
##LINEAR MODEL 2##
WINE_LIN <- WINE_NB %>%
  mutate(TARGET = log(TARGET + 1),
         logpH = log(pH),
         logcitric = log(CitricAcid + 1),
         logDensity = log(Density),
         logFSD = log(FreeSulfurDioxide + 1))
LogLinear_1 <- lm(TARGET ~ ., WINE_LIN)
null_Loglinear <- lm(TARGET ~ 1, WINE_LIN)
LogLinear_2 <- stepAIC(null_Loglinear, scope = list(upper=LogLinear_1 ), direction = "both")
summary(LogLinear_2)
plot(LogLinear_2)
car::mmps(LogLinear_2)

##MODEL COMPARISON##
hist(WINE_NB$TARGET, col = "cornflowerblue", breaks=seq(-.5, 8.5, 1), xlab = "Wine Cases",main = NULL)
lines(x = seq(0, 8, 1), y = pcount2, lty=2, lwd=2, col="red")
points(x = seq(0, 8, 1), y = pcount2, cex=2, pch=16, col="red")
lines(x = seq(0, 8, 1), y = pcount4, lty=2, lwd=2, col="green")
points(x = seq(0, 8, 1), y = pcount4, cex=2, pch=16, col="green")



####PREDICTIONS##
EVAL_WINE <- read.csv("wine-evaluation-data.csv", header = TRUE, sep = ",", na.strings = c("", " ", "NA"))
EVAL_WINE <- select(EVAL_WINE, -IN)
EVAL_COMPLETE <- EVAL_WINE %>%
  mutate(Alcohol = ifelse(Alcohol < 0, 0, ifelse(Alcohol > 24, 24, Alcohol)),
         LabelAppeal = as.factor(LabelAppeal),
         FixedAcidity = ifelse(FixedAcidity < 0, FixedAcidity * (-1), FixedAcidity),
         VolatileAcidity = ifelse(VolatileAcidity < 0, VolatileAcidity * (-1), VolatileAcidity),
         CitricAcid = ifelse(CitricAcid < 0, CitricAcid * (-1), CitricAcid),
         ResidualSugar = ifelse(ResidualSugar < 0, ResidualSugar * (-1), ResidualSugar),
         Chlorides = ifelse(Chlorides < 0, Chlorides * (-1), Chlorides),
         FreeSulfurDioxide = ifelse(FreeSulfurDioxide < 0, FreeSulfurDioxide * (-1), FreeSulfurDioxide),
         TotalSulfurDioxide = ifelse(TotalSulfurDioxide < 0, TotalSulfurDioxide * (-1), TotalSulfurDioxide),
         Sulphates = ifelse(Sulphates < 0, Sulphates * (-1), Sulphates),
         ALC_MISS = as.factor(ifelse(is.na(Alcohol), 1, 0)),
         TSD_MISS = as.factor(ifelse(is.na(TotalSulfurDioxide), 1, 0)),
         SUL_MISS = as.factor(ifelse(is.na(Sulphates), 1, 0)),
         FSD_MISS = as.factor(ifelse(is.na(FreeSulfurDioxide), 1, 0)),
         CHL_MISS = as.factor(ifelse(is.na(Chlorides), 1, 0)),
         RES_MISS = as.factor(ifelse(is.na(ResidualSugar), 1, 0)),
         pH_MISS = as.factor(ifelse(is.na(pH), 1, 0)),
         STARS = as.factor(ifelse(is.na(STARS), "No Rating", STARS)))


EVAL_IMPUTE <- mice(EVAL_COMPLETE, predictorMatrix = PREDMAT)
EVAL_IMP_CMP <- complete(EVAL_IMPUTE,5)

EVAL_IMP_CMP2 <- EVAL_IMP_CMP %>%
  mutate(Rating = as.factor(ifelse(STARS == 1 | STARS == 2, "Poor", ifelse(STARS == "No Rating", "None", "High"))),
         PoorLAB = as.factor(ifelse(LabelAppeal == -1 | LabelAppeal == -2, 1, 0)),
         TYP_pH = as.factor(ifelse(pH >= 3 & pH < 3.9, 0, 1)),
         TYP_ALC = as.factor(ifelse(Alcohol >= 10 & Alcohol <= 16, 0, 1)),
         SWEET = as.factor(ifelse(ResidualSugar < 10, "Dry", ifelse(ResidualSugar >= 10 & ResidualSugar <= 30, "Off-Dry",
                                                                    ifelse(ResidualSugar > 30 & ResidualSugar < 100, "Sweet", "Dessert")))),
         LAB_Rating = as.factor(ifelse(LabelAppeal == -1 | LabelAppeal == -2, "Unappealing", ifelse(LabelAppeal == 0, "Indifferent", "Appealing"))),
         NoSTARS = as.factor(ifelse(STARS == "No Rating", 1, 0)),
         SALT_RAT = log((ResidualSugar + 1)/(Chlorides + 1)),
         TYP_SWEET = as.factor(ifelse(ResidualSugar >= .2 & ResidualSugar <= 30, 1, 0)),
         pHunscaled = 10^(-pH),
         Impure = log(Chlorides + Sulphates),
         STARS = relevel(STARS, ref = 5))

EVAL_PRED_Zero <- predict(Zero_Poisson5, newdata = EVAL_IMP_CMP2, type = "zero")
EVAL_PRED <- predict(Zero_Poisson5, newdata = EVAL_IMP_CMP2, type = "response")
EVAL_SCORE_ZERO <- exp(EVAL_PRED_Zero)/ (1 + exp(EVAL_PRED_Zero))
EVAL_SCORE <- exp(EVAL_PRED)
EVAL_TARGET <- EVAL_SCORE * (1 - EVAL_PRED_Zero)
EVAL_SCORE_FINAL <- round(EVAL_PRED, 0)
EVAL_FINAL <- EVAL_IMP_CMP2 %>%
  mutate(TARGET = EVAL_SCORE_FINAL) %>%
  select(TARGET) 
  

write.table(EVAL_FINAL, file = "WINE_Predictions.csv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
