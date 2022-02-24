#Code for the BJGP Paper.

#Call libraries
library(glmmTMB)
library(broom.mixed)
library(data.table)
library(ggplot2)
library(sf)
library(ggsn)
library(mitml)
library(lme4)
library(miceadds)
library(abind)
library(mice)
library(micemd)
library(jomo)
library(norm2)
library(rms)

#Set working directory
setwd("")

#Source function for VPCs
devtools::source_url("https://github.com/NBakewell/CVSMedsPrescribing_BJGP/blob/main/VPCS_NB.R?raw=TRUE")

AnalysisData = read.csv("BJGP_Prescribing_CVS_Data.csv")                  

#Covariates
#65 plus -  centred at 15%
AnalysisData$CPerc65plus <- AnalysisData$Perc65plus - 15

#Female -  centred at 50%
AnalysisData$CPercFemale <- AnalysisData$PercFemale - 50

#Hypertension - centred at 15%
AnalysisData$CPercHyp <- AnalysisData$PercHyp -15

#Obese - centred at 10%
AnalysisData$CPercOB <- AnalysisData$PercOB - 10

#Diabetes - centred at 5%
AnalysisData$CPercDia<- AnalysisData$PercDia - 5

#Headcount - centred at 5 
AnalysisData$GPHeadcountCentered <- AnalysisData$GPHeadcount - 5

#Offset - re-scale to per 1000
AnalysisData$pats1000 <- AnalysisData$TotalPats/1000

#Transforming so that 2015 is time 0. This is not necessary, but 
#helpful in interpreting the interest, in case the reader is interested.
AnalysisData$yrc = AnalysisData$Year-2015

#include a quadratic
AnalysisData$yrc2 = AnalysisData$yrc^2

###############################################
############ NB Multi-level Models ############ 
###############################################

#Complete Case Analysis
colnames(AnalysisData)

#Extract only the variables we need for the analysis
GLMMTMBDATA <- AnalysisData[,c("CCGCODE", "GPPRACTICECODE", "NumItems","pats1000","yrc", "yrc2", "NHSRegion","CPercFemale", "CPerc65plus", "CPercOB", "CPercHyp", "IMD", "Urban", "GPHeadcountCentered", "CPercDia")]

#Make IDs factors for unadjusted analysis
GLMMTMBDATA$CCGCODE = as.factor(GLMMTMBDATA$CCGCODE)
GLMMTMBDATA$GPPRACTICECODE = as.factor(GLMMTMBDATA$GPPRACTICECODE)

#Take complete cases for the complete case analyses
GLMMTMBDATA_CCA <- GLMMTMBDATA[complete.cases(GLMMTMBDATA),]

#Make IDs factors
GLMMTMBDATA_CCA$CCGCODE = as.factor(GLMMTMBDATA_CCA$CCGCODE)
GLMMTMBDATA_CCA$GPPRACTICECODE = as.factor(as.character(GLMMTMBDATA_CCA$GPPRACTICECODE))

#Un-adjusted model
MLMNB_Intercept_only_Unadj <- glmmTMB(NumItems ~ 1 + offset(log(pats1000)) +
                                                (1|CCGCODE) + (1|GPPRACTICECODE), 
                                              data=GLMMTMBDATA, family=nbinom2)
summary(MLMNB_Intercept_only_Unadj)

#Summarise results
MLMNB_UnadjustedSummary <- tidy(MLMNB_Intercept_only_Unadj,exponentiate = TRUE, conf.method = "Wald", conf.int = TRUE)
write.csv(MLMNB_UnadjustedSummary, "MLMNB_UnadjustedModel_Paper_Final.csv", row.names=F)

#Note, we want to retrieve the SD of the variances of the random effects, so
#the delta method can be used to derive the SD. GLMMTMB parameterises the
#variances as log(SD), so we have exp(2*log(SD)) for the variance
#alpha is parameterised as 1/log(alpha)

#Alpha level-1
alpha_sd_unadj = sqrt( MLMNB_Intercept_only_Unadj$sdr$cov.fixed[2,2])
alpha_log_unadj  = (MLMNB_Intercept_only_Unadj$fit$par[2])
sqrt((alpha_sd_unadj^2)*(((-1)*(exp((-1)*alpha_log_unadj)))^2))
#SE = 0.0003972351

#Extract alpha as well.
1/exp(MLMNB_Intercept_only_Unadj$fit$par[2])
#Alpha = 0.04612227 

#GPP
GPP_sd_unadj = sqrt( MLMNB_Intercept_only_Unadj$sdr$cov.fixed[4,4])
GPP_logsd_unadj  = MLMNB_Intercept_only_Unadj$fit$par[4]
GPP_sd_unadj*2*(exp(2*GPP_logsd_unadj))
#SE=0.003345407 

#CCG
CCG_sd_unadj = sqrt(MLMNB_Intercept_only_Unadj$sdr$cov.fixed[3,3])
CCG_logsd_unadj  = MLMNB_Intercept_only_Unadj$fit$par[3]
CCG_sd_unadj*2*(exp(2*CCG_logsd_unadj))
#SE=0.008472335 

#Variance Partition Coefficients (VPCs)

#Code adapted from:
#http://www.bristol.ac.uk/cmm/media/leckie/articles/leckie2020.pdf

#Intercept
beta0_unadj = summary(MLMNB_Intercept_only_Unadj)$coefficients$cond[1,1]

#Level-3 (CCG) variance
sigma_CCG_unadj = summary(MLMNB_Intercept_only_Unadj)$varcor$cond$CCGCODE[1,1]

#Level-2 variance
sigma_GPP_unadj = summary(MLMNB_Intercept_only_Unadj)$varcor$cond$GPPRACTICECODE[1,1]

#Over-dispersion 
alpha_unadj = 1/(summary(MLMNB_Intercept_only_Unadj)$sigma)

#Marginal Expectation
expectation_unadj = exp(beta0_unadj + (sigma_CCG_unadj/2) + (sigma_GPP_unadj/2))

#Marginal Variance
variance_unadj = expectation_unadj + (expectation_unadj^2)*(exp(sigma_GPP_unadj + sigma_CCG_unadj)*(1+alpha_unadj)-1)

#Marginal Variance: CCG
variance_CCG_unadj = (expectation_unadj^2)*(exp(sigma_CCG_unadj)-1)

#Marginal Variance: GPP
variance_GPP_unadj = (expectation_unadj^2)*(exp(sigma_CCG_unadj)*(exp(sigma_GPP_unadj)-1))

#Marginal Variance: Level-1
variance_level1_unadj = expectation_unadj + ((expectation_unadj^2)*(exp(sigma_CCG_unadj + sigma_GPP_unadj)*alpha_unadj))

#Sum variances for denominator
vpc_denominator_unadj = variance_CCG_unadj + variance_GPP_unadj + variance_level1_unadj

#VPC CCG
vpc_CCG_unadj = variance_CCG_unadj/vpc_denominator_unadj
vpc_CCG_unadj
#VPC CCG (level 3) = 0.2250634

#VPC GPP
vpc_GPP_unadj = variance_GPP_unadj/vpc_denominator_unadj
vpc_GPP_unadj
#VPC GPP (level 2) = 0.6058095

#VPC Level-1
vpc_level1_unadj = variance_level1_unadj/vpc_denominator_unadj
vpc_level1_unadj
#VPC Level 1 (Between-GPP) = 0.1691271

#Extract random effects
RandomEffects <- ranef(MLMNB_Intercept_only_Unadj, condVar=TRUE)

#Make data frame for the CCG random effects
REs_CCGs_Unadj <- RandomEffects$cond$CCGCODE

#Make row names the first column
setDT(REs_CCGs_Unadj, keep.rownames = TRUE)[]

#Extract the variance
VarsCCG_Unadj <- attr(RandomEffects$cond$CCGCODE,"condVar")
VarsCCG_Unadj  <- t(as.data.frame(VarsCCG_Unadj))

#Covert to Standard Erros
SEsCCG_Unadj <- sqrt(VarsCCG_Unadj)

#Add to the random effects data for CCGs
REs_CCGs_Unadj$SEs <- SEsCCG_Unadj

#Rename columns
colnames(REs_CCGs_Unadj) <- c("CCG","RE_CCG","SE_CCG")

#Here, let's add NHS Region to the data.
#Read in NHS region lookup data. Could also
#join to the GLMMDATA.
NHS_Region_Lookup <- read.csv("CCG_NHSREGIONLOOKUP.csv")

#Create column for join 
NHS_Region_Lookup$CCG = NHS_Region_Lookup$CCG_CODE

#Join with CCG random effects
RandEfCCG_Unadj = merge(REs_CCGs_Unadj, NHS_Region_Lookup[, c("CCG", "NHSER19.Name")], by=c("CCG"))

#Make NHS region a factor for later
RandEfCCG_Unadj$NHSER19.Name = as.factor(RandEfCCG_Unadj$NHSER19.Name)

#Read in CCG shape file for the map
CCG <- st_read("c8aa66b8-d408-44b9-9641-ea45fb3344f02020315-1-68y7uz.trdmv.shp")

#This is a CCG Lookup table that includes codes that will link to the shape file.
CCG_Lookup <- read.csv("CCG_LU.csv")

#Add rank for map, as we want to call out the bottom and top 5 CCGs.
RandEfCCG_Unadj$rank = rank(RandEfCCG_Unadj$RE_CCG)

#Join with CCG lookup
RandEfCCG_Unadj = merge(RandEfCCG_Unadj, CCG_Lookup, by='CCG')

#Merge random effects data with the map
mergedshape <- merge(CCG, RandEfCCG_Unadj, by='ccg19cd')
mergedshape <- cbind(mergedshape, st_coordinates(st_centroid(mergedshape)))

main_mapres <- ggplot() +
                geom_sf(data = mergedshape,
                        aes(fill =RE_CCG)
                ) +
                geom_label(data = mergedshape[which(mergedshape$rank %in% c(191,190,189,188)),], aes(X, Y, label = rank), size = 3, fontface = "bold") + 
                scale_fill_gradientn(
                  colors = c("#9DBF9E", "#FCB97D", "#A84268"),
                  oob = scales::squish,
                  limits=c(-1,.6),
                  name = "Random Effects:\nCCG (Level-3)"
                ) +
                # Prevent ggplot from slightly expanding the map limits beyond the bounding box of the spatial objects
                coord_sf(expand = FALSE) +
                theme_void() +
                theme(
                  legend.justification = c(0, 1),
                  legend.position = c(0, .95)
                ) +
                north(mergedshape) +
                ggsn::scalebar(mergedshape, dist = 50, dist_unit = "km",
                               transform = TRUE, model = "WGS84", st.size=5, anchor=c(x=1.2,y=50.2))+
                ggtitle("Random Effects") +
                theme(
                  text = element_text(family = "Futura-Medium"),
                  legend.title = element_text(family = "Futura-Bold", size = 10),
                  legend.text = element_text(family = "Futura-Medium", size = 10),
                  plot.title = element_text(hjust = 0.5)
                )  

main_mapres
ggsave("MapUnadjusted_CCGREs.png", height=7, width=5)

#London inset. In this instance, it was checked that ranks 1-4 were
#in the London area, but this may not always be the case.
main_mapres +   
      geom_label(data = mergedshape[which(mergedshape$rank %in% c(1,2,3,4)),], aes(X, Y, label = rank), size = 3, fontface = "bold") + 
        coord_sf(
          xlim = c(-1, 1),
          ylim = c(51.2, 51.8),
          expand = FALSE
        ) +
        theme(legend.position = "none")

ggsave("MapUnadjusted_Inset_CCGREs.png", height=7, width=5)

#Let's calculate upper and lower bounds for the EB estimates, 
#error bars sqrt(2)*SE
RandEfCCG_Unadj$lower_ebCCG <- RandEfCCG_Unadj$RE_CCG - (sqrt(2)*RandEfCCG_Unadj$SE_CCG)
RandEfCCG_Unadj$upper_ebCCG <- RandEfCCG_Unadj$RE_CCG + (sqrt(2)*RandEfCCG_Unadj$SE_CCG)

#Create new data object for plots, although not necessary
CCGRE= RandEfCCG_Unadj
CCGRE <- CCGRE[order(CCGRE$RE_CCG),]
CCGRE$rank <- rank(CCGRE$RE_CCG)

#Change reference level so that the caterpillar plots are arranged
#North to South (West to East)
CCGRE$NHSER19.Name = factor(CCGRE$NHSER19.Name, levels=c("North East and Yorkshire", "North West", "Midlands", "East of England", "London", "South East", "South West")) 

#Save random effects (on natural, log scale)
write.csv(CCGRE, "CCGRE_UnadjustedCCGEBs_Final.csv", row.names=F)

#Make random effects (CCG) caterpillar plots
Step1_CCG <- ggplot(CCGRE, aes(x=RE_CCG, y=reorder(CCG,RE_CCG))) + 
                geom_point(colour="darkblue",size=1.75) +
                geom_vline(xintercept=0, linetype="dashed", 
                           color = "red", size=1) +
                geom_errorbarh(data=CCGRE, mapping=aes(xmin=upper_ebCCG, xmax=lower_ebCCG), height=0.2, size=1, color="blue3",alpha=0.5)  +
                facet_wrap(~NHSER19.Name, scales="free_x") +
                theme_classic() + coord_flip()

Step2_CCG <- Step1_CCG + 
             xlab("EB Estimates of Random Intercepts (Level-3, CCG)") +
             ylab("CCG Rank (within NHS Region)")  + 
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
             scale_y_discrete(labels = NULL, breaks = NULL,expand=expansion(0.02)) 
Step2_CCG

#Save image
ggsave("CaterpillarPlots_Unadjusted_CCGREs.png", height=5, width=8)

#Let's also create caterpillar plots for GP practices
REs_GPP <- RandomEffects$cond$GPPRACTICECODE

#Make row names the first column 
setDT(REs_GPP, keep.rownames = TRUE)[]

#Extract variance for the GP practice random effects. 
VarsGPP <- attr(RandomEffects$cond$GPPRACTICECODE,"condVar")
VarsGPP  <- t(as.data.frame(VarsGPP))

#Make into standard errors 
SEsGPP <- sqrt(VarsGPP)

#Add to RE data frame of GPP random effects
REs_GPP$SEs <- SEsGPP

#Rename Columns
colnames(REs_GPP) <- c("GPPRACTICECODE", "RE_GPP","SE_GPP")

#We need to get the GP practice's region from the 
#analysis data set. So, we need the CCG.
colnames(GPP_Region_Data)
GPP_Region_Data <- GLMMTMBDATA %>%
                   select(GPPRACTICECODE, NHSRegion) %>%
                   distinct()

REs_GPP_Unadj = merge(REs_GPP, GPP_Region_Data, by=c("GPPRACTICECODE"))


#Let's calculate upper and lower bounds for the EB estimates, 
#error bars sqrt(2)*SE
REs_GPP_Unadj$lower_ebGPP <- REs_GPP_Unadj$RE_GPP - (sqrt(2)*REs_GPP_Unadj$SE_GPP)
REs_GPP_Unadj$upper_ebGPP <- REs_GPP_Unadj$RE_GPP + (sqrt(2)*REs_GPP_Unadj$SE_GPP)


#Let's make a copy of the data for caterpillar plots
GPPRE= REs_GPP_Unadj
GPPRE <- GPPRE[order(GPPRE$RE_GPP),]
GPPRE$rank <- rank(GPPRE$RE_GPP)

GPPRE$NHSRegion = factor(GPPRE$NHSRegion, levels=c("North East and Yorkshire", "North West", "Midlands", "East of England", "London", "South East", "South West")) 

write.csv(GPPRE, "GPPRE_UnadjustedCCGEBs_Final.csv", row.names=F)

#Create Caterpillar plot
Step1_GPP <- ggplot(GPPRE, aes(x=RE_GPP, y=reorder(GPPRACTICECODE,RE_GPP))) + 
              geom_point(colour="darkblue",size=1.75) +
              geom_vline(xintercept=0, linetype="dashed", 
                         color = "red", size=1) +
              geom_errorbarh(data=GPPRE, mapping=aes(xmin=upper_ebGPP, xmax=lower_ebGPP), height=0.2, size=1, color="blue3",alpha=0.5)  +
              facet_wrap(~NHSRegion, scales="free_x") +
              theme_classic() + coord_flip()

Step2_GPP <- Step1_GPP + xlab("EB Estimates of Random Intercepts (Level-2, GPP)") +
              ylab("GPP Rank (within NHS Region)")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +scale_y_discrete(labels = NULL, breaks = NULL,expand=expansion(0.02)) 
Step2_GPP

ggsave("GPPCaterpillar_Unadjusted.png", height=5, width=8)


#For this analysis, we will do the complete cases (CCA) only.
GLMMTMBDATA_CCA$NHSRegion = as.factor(GLMMTMBDATA_CCA$NHSRegion)
GLMMTMBDATA_CCA$NHSRegion =relevel(GLMMTMBDATA_CCA$NHSRegion, ref="London")

#Model Adjusted (CCA) - glmmTMB
MLMNB_Intercept_only_Adjusted <- glmmTMB(NumItems ~ 1 + CPercFemale + CPerc65plus +
                                 CPercOB + CPercHyp + factor(IMD) + 
                                 factor(Urban) + GPHeadcountCentered  + NHSRegion +
                                 CPercDia + yrc + yrc2 + offset(log(pats1000)) +
                                 (1|CCGCODE) + (1|GPPRACTICECODE), 
                                 data=GLMMTMBDATA_CCA, family=nbinom2)

MLMAdjNB_AdjustedSummary <- tidy(MLMNB_Intercept_only_Adjusted, exponentiate=TRUE, conf.method = "Wald", conf.int = TRUE)
write.csv(MLMAdjNB_AdjustedSummary, "MLMAdjNB_AdjustedSummary_Final_CCA.csv", row.names=F, na="")

#As above, need to use delta method to derive SE of variance of the random effects.
#Alpha
alpha_sd_adj = sqrt( MLMNB_Intercept_only_Adjusted$sdr$cov.fixed[23,23])
alpha_log_adj  = (MLMNB_Intercept_only_Adjusted$fit$par[21])
sqrt((alpha_sd_adj^2)*(((-1)*(exp((-1)*alpha_log_adj)))^2))
#SE=8.745775e-05

#Also extract alpha
1/exp( (MLMNB_Intercept_only_Adjusted$fit$par[21]))
#alpha=0.01004833 

#GPP-level random effects' variance SE
GPP_sd_adj = sqrt( MLMNB_Intercept_only_Adjusted$sdr$cov.fixed[23,23])
GPP_logsd_adj  = MLMNB_Intercept_only_Adjusted$fit$par[23]
GPP_sd_adj*2*(exp(2*GPP_logsd_adj))
#SE=0.001315624 

#CCG-level random effects' variance SE
CCG_sd_adj = sqrt( MLMNB_Intercept_only_Adjusted$sdr$cov.fixed[22,22])
CCG_logsd_adj  = MLMNB_Intercept_only_Adjusted$fit$par[22]
CCG_sd_adj*2*(exp(2*CCG_logsd_adj))
#SE=0.001847085 

#Variance Partition Coefficients (VPCs)
xb_adj_CC <- predict(MLMNB_Intercept_only_Adjusted)
xb_adj_CC <- data.frame(xb_adj_CC)

#Level-3 (CCG) variance
sigma_CCG_adj <- summary(MLMNB_Intercept_only_Adjusted)$varcor$cond$CCGCODE[1,1]

#Level-2 (GPP) variance
sigma_GPP_adj <- summary(MLMNB_Intercept_only_Adjusted)$varcor$cond$GPPRACTICECODE[1,1]

# Over-dispersion parameter
alpha_adj <- 1/(summary(MLMNB_Intercept_only_Adjusted)$sigma)

# Marginal expectation
xb_adj_CC$expectation <- exp(xb_adj_CC$xb_adj_CC + (sigma_CCG_adj/2) + (sigma_GPP_adj/2))

# Marginal variance
xb_adj_CC$marginalvariance <- xb_adj_CC$expectation +
                             (xb_adj_CC$expectation^2)*(exp(sigma_CCG_adj + sigma_GPP_adj)*(1 + alpha_adj) - 1)

# Marginal variance: Level-3 component
xb_adj_CC$variance_CCG <- ((xb_adj_CC$expectation^2)*(exp(sigma_CCG_adj) - 1))

# Marginal variance: Level-2 component
xb_adj_CC$variance_GPP <- xb_adj_CC$expectation^2*exp(sigma_CCG_adj)*(exp(sigma_GPP_adj) - 1)

# Marginal variance: Level-1 component
xb_adj_CC$variance_level1 <- xb_adj_CC$expectation + (xb_adj_CC$expectation^2)*exp(sigma_CCG_adj + sigma_GPP_adj)*alpha_adj

#VPC denominator
xb_adj_CC$vpc_denominator_adj = xb_adj_CC$variance_CCG + xb_adj_CC$variance_GPP + xb_adj_CC$variance_level1

# Level-3 VPC
xb_adj_CC$VPC_CCG <- xb_adj_CC$variance_CCG/xb_adj_CC$vpc_denominator_adj

# Level-2 VPC
xb_adj_CC$VPC_GPP <- xb_adj_CC$variance_GPP/xb_adj_CC$vpc_denominator_adj

# Level-1 VPC
xb_adj_CC$VPC_level1 <- xb_adj_CC$variance_level1/xb_adj_CC$vpc_denominator_adj

#Get mean of VPCs
sapply(xb_adj_CC, mean)[8:10]

# VPC_CCG    VPC_GPP VPC_level1 
#0.1467169  0.7494543  0.1038288 

#Let's do the Likelihood ratio tests for the CCA
#This is the log-Likelihood for the adjusted model
#Note, only complete cases data here.
LL_MLMAdjNB = (-1)*MLMNB_Intercept_only_Adjusted$obj$fn()[[1]]


#Model with only a GPP random intercept
MLMAdjNB_AdjustedGPPIntercept_only <- glmmTMB(NumItems ~ 1 + CPercFemale + CPerc65plus +
                                    CPercOB + CPercHyp + factor(IMDScoreQuintile_LOCF) + 
                                    factor(Urban) + GPHeadcountCentered  + NHSER19.Name + 
                                    CPercDia + yrc + I(yrc^2) + offset(log(pats1000)) + 
                                    (1|GPPRACTICE), 
                                    data=GLMMTMBDATA_CCA, family=nbinom2)

LL_MLMAdjNB_GPP_RIONLY = (-1)*MLMAdjNB_AdjustedGPPIntercept_only$obj$fn()[[1]]

#Model with only a no random intercept
MLMAdjNB_AdjustedNoIntercepts <- glmmTMB(NumItems ~ 1 + CPercFemale + CPerc65plus +
                                         CPercOB + CPercHyp + factor(IMDScoreQuintile_LOCF) + 
                                         factor(Urban) + GPHeadcountCentered  + NHSER19.Name + 
                                         CPercDia + yrc + I(yrc^2) + offset(log(pats1000)), 
                                         data=GLMMTMBDATA_CCA, family=nbinom2)

LL_MLMAdjNB_NORI = (-1)*MLMAdjNB_AdjustedNoIntercepts$obj$fn()[[1]]

#Let's do a series of tests on the random intercepts

#No Random Intercepts vs. GPP
LL_RIvsNone = -2*(LL_MLMAdjNB_NORI-LL_MLMAdjNB_GPP_RIONLY)
0.5*(1-pchisq(LL_RIvsNone, df=1))

#Random Intercepts for CCG and GPP vs. only GPP random intercept model
LL_BothRIvsGPPOnly = -2*(LL_MLMAdjNB_GPP_RIONLY-LL_MLMAdjNB)
0.5*(1-pchisq(LL_BothRIvsGPPOnly, df=1))

#Imputation Data
#Read in imputed data from Blimp
imp.data = read.csv("imputations_Blimp.csv", header=F)

#Name columns
#Note, variables were centred prior to imputation (year linearly transformed so that year 0=2015),
#to improve convergence. Number of patients were rescaled to improve imputation (per 1000 (i.e. divided by 1000)), 
#but the number of items was not rescaled.
#Also, Blimp requires IDs and factor variables to be coded in numbers, hence the use of "NUM_CODE
colnames(imp.data) = c("Imputation", "CCG_NUM_CODE", "GPP_NUM_CODE", "NumItems","pats1000" ,  
                       "yrc" , "yrc2" , "CPercFemale", "C65plus",  "CPercOB" ,   
                       "CPercHyp", "IMD", "Urban","GPHeadcount", "CPercDia", "NHSRegion_NUM_CODE")


#We need to get the actual codes for the IDs/categorical variables coded with a numeric. 
Data_Codes = GLMMTMBDATA %>%
              select(GPPRACTICECODE, CCGCODE, NHSRegion) %>%
              distinct()

#Create numeric versions in the data.
Data_Codes$NHSRegion_NUM_CODE = as.numeric(as.factor(Data_Codes$NHSRegion))
Data_Codes$CCG_NUM_CODE = as.numeric(as.factor(Data_Codes$CCGCODE))
Data_Codes$GPP_NUM_CODE = as.numeric(as.factor(Data_Codes$GPPRACTICECODE))

#CCG numeric lookup table
CCG_LookupTable = Data_Codes %>%
                  select(CCG_NUM_CODE,CCGCODE) %>%
                   distinct()

#GP practice numeric lookup table
GPP_LookupTable = Data_Codes %>%
                  select(GPP_NUM_CODE,GPPRACTICECODE) %>%
                  distinct()

#NHS Region numeric lookup table
NHSRegion_LookupTable =Data_Codes %>%
                        select(NHSRegion_NUM_CODE,NHSRegion) %>%
                        distinct()


imp.data = left_join(imp.data, CCG_LookupTable, by=c("CCG_NUM_CODE"))
imp.data = left_join(imp.data, GPP_LookupTable, by=c("GPP_NUM_CODE"))
imp.data = left_join(imp.data, NHSRegion_LookupTable, by=c("NHSRegion_NUM_CODE"))

#Imputations we done in Blimp. IMD was imputed as continuous for the 25 records
#missing data, primarily for computational efficiency. 
#We need to round the values, we are treating IMD as a factor. Use simple rounding.
imp.data$IMD = if_else(imp.data$IMD <1, 1,
                             if_else(imp.data$IMD>5, 5, round(imp.data$IMD)))

imp.data$IMD = as.factor(imp.data$IMD)

#Let's make sure variables are formatted properly
imp.data$CCGCODE = as.factor(imp.data$CCGCODE)
imp.data$GPPRACTICECODE = as.factor(imp.data$GPPRACTICECODE)
imp.data$NHSRegion = as.factor(imp.data$NHSRegion)

#Make London the reference level
imp.data$NHSRegion = relevel(imp.data$NHSRegion, ref="London")

#Did not truncate or round the GP headcount variable, as per this article:
#https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-57

#Create a list of the 30 datasets
implist <-  as.mitml.list(split(imp.data, imp.data$Imputation))

#Code adapted from Applied Multiple Imputation Advantages, Pitfalls, New Developments and Applications in R
NBinModel <- est <- se <- rvarCCG <- rvarGPP <-Var_REVar_GPP <- Var_REVar_CCG <- theta <- vector(length=30, mode="list")
VPCs_MI = data.frame(VPC_CCG=NULL, VPC_GPP=NULL,VPC_level1=NULL)

#Note, not everything is in this list, such as EB estimates, but that is best done 
#on its own
for (i in 1:30) {
          cat(i,"\n")
          com = implist[[i]]
          com$CCGCODE = as.factor(com$CCGCODE)
          com$GPPRACTICECODE = as.factor(com$GPPRACTICECODE)
          #Probably do not need to do this, but just in case.
          com$NHSRegion = as.factor(com$NHSRegion)
          com$NHSRegion = relevel(com$NHSRegion, ref="London")
          NBinModel[[i]] =glmmTMB(NumItems ~ 1 + CPercFemale + C65plus +
                              CPercOB + CPercHyp + factor(IMD) + 
                              factor(Urban) + GPHeadcount  + NHSRegion + 
                              CPercDia + yrc + yrc2 + offset(log(pats1000)) +
                               (1 |CCGCODE) + (1|GPPRACTICECODE), 
                             family=nbinom2, data=com)
                          s <- summary(NBinModel[[i]])
                          #Estimates for fixed effects
                          est[[i]] <- s$coefficients$cond[,1]
                          #Standard Errors for fixed effects
                          se[[i]] <- s$coefficients$cond[,2]
                          #Random Effect Variance for GPP
                          rvarCCG[[i]] <- attr(glmmTMB::VarCorr(NBinModel[[i]])[[1]]$CCGCODE, "stddev")^2
                          #Random Effect Variance for GPP
                          rvarGPP[[i]] <- attr(glmmTMB::VarCorr(NBinModel[[i]])[[1]]$GPPRACTICECODE, "stddev")^2
                          #Variance of the Random Effect Variance for GPP
                          Var_REVar_GPP[[i]] <- (NBinModel[[i]]$sdr$cov.fixed[23,23])*(2*exp(2*NBinModel[[i]]$fit$par[23]))^2
                          #Variance of the Random Effect Variance for CCG
                          Var_REVar_CCG[[i]] <- (NBinModel[[i]]$sdr$cov.fixed[22,22])*(2*exp(2*NBinModel[[i]]$fit$par[22]))^2
                          #Over-dispersion term
                          theta[[i]] <- s$sigma
                          #VPCs
                          VPCs_MI = rbind(VPCs_MI,
                                          data.frame(Imputation=i,VPCs_NB(NBinModel[[1]])))
                          
        }


#Multiple ways for the results to be combined

#here is one way using the norm2 package.
MI_Results= miInference(est, se)

#Let's use the mice pool function
MIResults_POOL = summary(pool(NBinModel), conf.int = TRUE)

#extract necessary items
MIResults_Export = MIResults_POOL[c("term","estimate", "2.5 %", "97.5 %")]

#take exponential for RRs
MIResults_Export[,c("estimate","2.5 %", "97.5 %")] = exp(MIResults_Export[,c("estimate","2.5 %", "97.5 %")])
MIResults_Export$p.val = MIResults_POOL$p.value

#Write results to a data file
write.csv(MIResults_Export, "MIResults_Combined_Paper_Final.csv", row.names=F)

#Get VPCs
colMeans(VPCs_MI[,2:ncol(VPCs_MI)])
# VPC_CCG    VPC_GPP VPC_level1 
# 0.09736364 0.59142180 0.31121456 

#Get random variance
apply(abind(rvarCCG, along=2), 1, mean)
#Var=0.01549978 

#Get SE using Rubin's Rules for finding variance
Var_B_CCG = var(unlist(rvarCCG))
Var_W_CCG = mean(unlist(Var_REVar_CCG))
Var_T_CCG = Var_W_CCG + (1+(1/30))*Var_B_CCG
sqrt(Var_T_CCG)
#SE=0.001907356

#GPP
apply(abind(rvarGPP, along=2), 1, mean)
#Var=0.0894595

#Get SE using Rubin's Rules for finding variance
Var_B_GPP = var(unlist(rvarGPP))
Var_W_GPP = mean(unlist(Var_REVar_GPP))
Var_T_GPP = Var_W_GPP + (1+(1/30))*Var_B_GPP
sqrt(Var_T_GPP)
#SE= 0.001968267

#Alpha - overdispersion term
mean(1/(unlist(theta)))
Var_Alpha_B = var(unlist(1/(unlist(theta))))
Var_alpha_MI = NULL
for (i in (1:30)) {
  alpha_sd_adj_MI = sqrt(NBinModel[[i]]$sdr$cov.fixed[21,21])
  alpha_log_adj_MI  = (NBinModel[[i]]$fit$par[21])
  Var_alpha_MI[i] <- (alpha_sd_adj_MI^2)*(((-1)*(exp((-1)*alpha_log_adj_MI)))^2)
}

#Get SE using Rubin's Rules for finding variance
Var_Alpha_W = mean(Var_alpha_MI)
Var_T_Alpha = Var_Alpha_W + (1+(1/30))*Var_Alpha_B
sqrt(Var_T_Alpha)
#SE=0.0004036903

#Let's get the random effects
#Note, the SEs are stored in an attribute.
RandEf = RandEfCCG = RandEfGPP =vector(length=30, mode="list")

for (i in 1:30) {
                cat(i,"\n")
                RandEf[[i]] = ranef(NBinModel[[i]], condVar=TRUE)
                RandEfCCGDF =  setDT(RandEf[[i]]$cond$CCGCODE, keep.rownames = T)
                RandEfCCGDF$VarsCCG = t(as.data.frame(attr(RandEf[[i]]$cond$CCGCODE,"condVar")))
                RandEfCCG[[i]] = RandEfCCGDF
                RandEfGPPDF =  setDT(RandEf[[i]]$cond$GPPRACTICECODE, keep.rownames = T)
                RandEfGPPDF$VarsGPP = t(as.data.frame(attr(RandEf[[i]]$cond$GPPRACTICECODE,"condVar")))
                RandEfGPP[[i]] = RandEfGPPDF
}

#CCGs
RandEffects_MI_combined_CCG = bind_rows(RandEfCCG, .id = "column_label")
colnames(RandEffects_MI_combined_CCG)[3] = "RE"

#combine using Rubin's rules
RandEffects_MI_combined_CCG  = RandEffects_MI_combined_CCG %>%
                                group_by(rn) %>%
                                dplyr::summarise(RE_CCG=mean(RE), Var_B_RE_CCG=var(RE), Var_W_CCG=mean(VarsCCG)) %>%
                                mutate(Vars_T_RE_CCG=Var_W_CCG + ((1+(1/30))*Var_B_RE_CCG))

colnames(RandEffects_MI_combined_CCG)[[1]]  = "CCG_CODE"

#Need to extract region from lookup table
RandEffects_MI_combined_CCG = merge(RandEffects_MI_combined_CCG, Lookup[, c("CCG_CODE", "NHSER19.Name")], by=c("CCG_CODE"))
RandEffects_MI_combined_CCG$NHSER19.Name = as.factor(RandEffects_MI_combined_CCG$NHSER19.Name)

#rank
RandEffects_MI_combined_CCG$rank = rank(RandEffects_MI_combined_CCG$RE_CCG)

#Need to merge for mapping
RandEffects_MI_combined_CCG = left_join(RandEffects_MI_combined_CCG, CCG_Lookup, by=c('CCG_CODE'='CCG'))

#Merge for shape
mergedshape <- merge(CCG, RandEffects_MI_combined_CCG, by='ccg19cd')
mergedshape <- cbind(mergedshape, st_coordinates(st_centroid(mergedshape)))

#Make map for adjusted MI combined model - CCG level
main_mapres <-  ggplot() +
                    geom_sf(data = mergedshape,
                            aes(fill =RE_CCG)
                    ) +
                    geom_label(data = mergedshape[which(mergedshape$rank %in% c(4,190,189,188)),], aes(X, Y, label = rank), size = 3, fontface = "bold") + 
                    scale_fill_gradientn(
                      colors = c("#9DBF9E", "#FCB97D", "#A84268"),
                      oob = scales::squish,
                      limits=c(-.4,.4),
                      name = "Random Effects:\nCCG (Level-3)"
                    ) +
                    # Prevent ggplot from slightly expanding the map limits beyond the bounding box of the spatial objects
                    coord_sf(expand = FALSE) +
                    theme_void() +
                    theme(
                      legend.justification = c(0, 1),
                      legend.position = c(0, .95)
                    ) +
                    north(mergedshape) +
                    ggsn::scalebar(mergedshape, dist = 50, dist_unit = "km",
                                   transform = TRUE, model = "WGS84", st.size=5, anchor=c(x=1.2,y=50.2))+
                    ggtitle("Random Effects") +
                    theme(
                      text = element_text(family = "Futura-Medium"),
                      legend.title = element_text(family = "Futura-Bold", size = 10),
                      legend.text = element_text(family = "Futura-Medium", size = 10),
                      plot.title = element_text(hjust = 0.5)
                    )  
main_mapres
ggsave("MapAdjusted.png", height=7, width=5)

#Get London inset
main_mapres +   geom_label(data = mergedshape[which(mergedshape$rank %in% c(1,2,3,191)),], aes(X, Y, label = rank), size = 3, fontface = "bold") + 
  coord_sf(
    xlim = c(-1, 1),
    ylim = c(51.2, 51.8),
    expand = FALSE
  ) +
  theme(legend.position = "none")

ggsave("MapAdjusted_Inset.png", height=7, width=5)

#MI combined GPP random effect
RandEffects_MI_combined_GPP = bind_rows(RandEfGPP, .id = "column_label")
colnames(RandEffects_MI_combined_GPP)[3] = "RE"

#combine using Rubin's rules
RandEffects_MI_combined_GPP = RandEffects_MI_combined_GPP %>%
                              group_by(rn) %>%
                              dplyr::summarise(RE_GPP=mean(RE), Var_B_RE_GPP=var(RE), Var_W_GPP=mean(VarsGPP))  %>%
                              mutate(Vars_T_RE_GPP=Var_W_GPP + ((1+(1/30))*Var_B_RE_GPP))


colnames(RandEffects_MI_combined_GPP)[1] = "GPPRACTICE"

#Merge with the GPP region lookup table created for unadjusted analysis
RandEffects_MI_combined_GPP = left_join(RandEffects_MI_combined_GPP, GPP_Region_Data, by=c("GPPRACTICE"="GPPRACTICECODE"))

#Relevel so the order aligns with other caterpillar plots
RandEffects_MI_combined_GPP$NHSRegion = factor(RandEffects_MI_combined_GPP$NHSRegion, levels=c("North East and Yorkshire", "North West", "Midlands", "East of England", "London", "South East", "South West")) 

#Create error bars
RandEffects_MI_combined_GPP$upper_ebGPP = RandEffects_MI_combined_GPP$RE_GPP + sqrt(2)*sqrt(RandEffects_MI_combined_GPP$Vars_T_RE_GPP)
RandEffects_MI_combined_GPP$lower_ebGPP = RandEffects_MI_combined_GPP$RE_GPP - sqrt(2)*sqrt(RandEffects_MI_combined_GPP$Vars_T_RE_GPP)

#Save EB estimates
write.csv(RandEffects_MI_combined_GPP, "RandEffects_MI_combined_GPP.csv", row.names =F)

#Create GPP Caterpillar plots
Step1_GPP_MI <- ggplot(RandEffects_MI_combined_GPP, aes(x=RE_GPP, y=reorder(GPPRACTICE,RE_GPP))) + 
                geom_point(colour="darkblue",size=1.75) +
                geom_vline(xintercept=0, linetype="dashed", 
                           color = "red", size=1) +
                geom_errorbarh(data=RandEffects_MI_combined_GPP, mapping=aes(xmin=upper_ebGPP, xmax=lower_ebGPP), height=0.2, size=1, color="blue3",alpha=0.5)  +
                facet_wrap(~NHSRegion, scales="free_x") +
                theme_classic() + coord_flip()

Step2_GPP_MI <- Step1_GPP_MI + xlab("EB Estimates of Random Intercepts (Level-2, GPP)") +
                ylab("GPP Rank (within NHS Region)")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +scale_y_discrete(labels = NULL, breaks = NULL,expand=expansion(0.02)) 
Step2_GPP_MI

#Save GPP MI-combined Caterpillar plot
ggsave("GPPCaterpillarPlot_MI.png", height=5, width=8)

#Now, let's create CCG Caterpillar plots

#Relevel NHS region so the plots are in the same order as others.

RandEffects_MI_combined_CCG$NHSER19.Name = factor(RandEffects_MI_combined_CCG$NHSER19.Name, levels=c("North East and Yorkshire", "North West", "Midlands", "East of England", "London", "South East", "South West")) 

#Create error bars
RandEffects_MI_combined_CCG$upper_ebGPP = RandEffects_MI_combined_CCG$RE_CCG + sqrt(2)*sqrt(RandEffects_MI_combined_CCG$Vars_T_RE_CCG)
RandEffects_MI_combined_CCG$lower_ebGPP = RandEffects_MI_combined_CCG$RE_CCG - sqrt(2)*sqrt(RandEffects_MI_combined_CCG$Vars_T_RE_CCG)

#Save CCG EB estimates.
write.csv(RandEffects_MI_combined_CCG, "RandEffects_MI_combined_CCG.csv", row.names =F)

#Create Caterpillar plots
Step1_CCG_MI <- ggplot(RandEffects_MI_combined_CCG, aes(x=RE_CCG, y=reorder(CCG_CODE,RE_CCG))) + 
                        geom_point(colour="darkblue",size=1.75) +
                        geom_vline(xintercept=0, linetype="dashed", 
                                   color = "red", size=1) +
                        geom_errorbarh(data=RandEffects_MI_combined_CCG, mapping=aes(xmin=upper_ebGPP, xmax=lower_ebGPP), height=0.2, size=1, color="blue3",alpha=0.5)  +
                        facet_wrap(~NHSER19.Name, scales="free_x") +
                        theme_classic() + coord_flip()

Step2_CCG_MI <- Step1_CCG_MI + xlab("EB Estimates of Random Intercepts (Level-3, CCG)") +
                  ylab("CCG Rank (within NHS Region)")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +scale_y_discrete(labels = NULL, breaks = NULL,expand=expansion(0.02)) 
Step2_CCG_MI

#Save CCG MI-combined Caterpillar plot
ggsave("RandEffects_MI_combined_CCG.png", height=5, width=8)

#End of primary anlyses and figures

#######################################
########### MI-combined LRTs ##########
#######################################

#This section conducts the MI-combined LRTs

#Make empty list of length 30 for GPP intercept only model results.
NBinModel_GPPInt <- vector(length=30, mode="list")

#GPP random intercept only
for (i in 1:30) {
                cat(i,"\n")
                com = implist[[i]]
                com$GPPRACTICECODE = as.factor(com$GPPRACTICECODE)
                com$Urban = as.factor(com$Urban)
                com$IMD = as.factor(com$IMD)
                com$NHSRegion = as.factor(com$NHSRegion)
                com$NHSRegion = relevel(com$NHSRegion, ref="London")
                NBinModel_GPPInt[[i]] =glmmTMB(NumItems ~ 1 + CPercFemale + C65plus +
                                          CPercOB + CPercHyp + factor(IMD) + 
                                          factor(Urban) + GPHeadcount  + NHSRegion +
                                          CPercDia +yrc + yrc2 + offset(log(pats1000)) +
                                          (1|GPPRACTICECODE), 
                                        family=nbinom2, data=com)

}

#Make empty list of length 30 for no CCG and GPP and random intercepts model results.
NBinModel_NoRIs <- vector(length=30, mode="list")

#NB binomial model with no GPP and CCG random intercepts
for (i in 1:30) {
                  cat(i,"\n")
                  com = implist[[i]]
                  com$Urban = as.factor(com$Urban)
                  com$IMD = as.factor(com$IMD)
                  com$NHSRegion = as.factor(com$NHSRegion)
                  com$NHSRegion = relevel(com$NHSRegion, ref="London")
                  NBinModel_NoRIs[[i]] =glmmTMB(NumItems ~ 1 + CPercFemale + C65plus +
                                                   CPercOB + CPercHyp + IMD + 
                                                   Urban + GPHeadcount + + NHSRegion + CPercDia +
                                                     yrc + yrc2 + offset(log(pats1000)),
                                                 family=nbinom2, data=com)
  }

#Full Model, with random intercepts for CCG and GP practice
NBin_Full = as.mitml.result(NBinModel)

#Random GP practice intercept model
NBin_GPPInt = as.mitml.result(NBinModel_GPPInt)

#No GP practice or CCG random intercepts model
NBin_None = as.mitml.result(NBinModel_NoRIs)

#LRTs - mice package implements Li et al.'s method "D2", no need to indicate
#this, as the function will detect this is the only possible method
#using glmmTMB

anova(NBin_GPPInt, NBin_None, method="D2")
# F.value       df1       df2     P(>F)       RIV 
# 1 vs 2  15596.104         1 6.675e+05     0.000     0.007 

anova(NBin_Full, NBin_GPPInt, method="D2")
# F.value       df1       df2     P(>F)       RIV 
# 1 vs 2    695.926         1 8.076e+05     0.000     0.006 

#Can also calculate manually - for example:
#can obtain formulae nicely summarised here:
#https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-9-57/tables/2
None = unlist(lapply(NBinModel_NoRIs, logLik))
GPPInt_Only  = unlist(lapply(NBin_GPPInt, logLik))

#DF/number of parameters being tested 
#In this case, we will assume 1, rather than
#a mixture given that the variance of the random effect
#are at the boundary of the parameter space.
k=1
nimps=30

#Create vector of LRTs for all 30 models
LRT_NonevsGPPInt = (-2)*(None-GPPInt_Only)

#Take mean
Mean_LRT_NonevsGPPInt = mean(LRT_NonevsGPPInt)

#Calculate relative increase in variance
RIV2_NonevsGPPInt = ((nimps+1)/(nimps*(nimps-1)))*sum((sqrt(LRT_NonevsGPPInt)-sqrt(Mean_LRT_NonevsGPPInt))^2)

#Calculate denominator degress of freedom
V2_NonevsGPPInt = (k^(-(3/nimps)))*(nimps-1)*((1+(1/RIV2_NonevsGPPInt))^2)

TestStat_NonevsGPPInt = (1/(1 + RIV2_NonevsGPPInt))*((Mean_LRT_NonevsGPPInt/k) - (((nimps+1)/(nimps-1))*RIV2_NonevsGPPInt))
#Test Statistics=15596.1. Aligns with the anova function.

#GPP random intercepts model vs GPP AND CCG random intercepts model
#Create vector of LRTs for all 30 models
NBin_CCG_GPP  = unlist(lapply(NBin_Full, logLik))

LRT_GPPIntvsCCGGPPInts = (-2)*(GPPInt_Only-NBin_CCG_GPP)

#Take mean
Mean_LRT_GPPIntvsCCGGPPInts = mean(LRT_GPPIntvsCCGGPPInts)

#Calculate relative increase in variance
RIV2_GPPIntvsCCGGPPInt = ((nimps+1)/(nimps*(nimps-1)))*sum((sqrt(LRT_GPPIntvsCCGGPPInts)-sqrt(Mean_LRT_GPPIntvsCCGGPPInts))^2)

#Calculate denominator degress of freedom
V2_GPPIntvsCCGGPPInt = (k^(-(3/nimps)))*(nimps-1)*((1+(1/RIV2_GPPIntvsCCGGPPInt))^2)

TestStat_GPPIntvsCCGGPPInt = (1/(1 + RIV2_GPPIntvsCCGGPPInt))*((Mean_LRT_GPPIntvsCCGGPPInts/k) - (((nimps+1)/(nimps-1))*RIV2_GPPIntvsCCGGPPInt))
#Test Statistics=695.9261. Aligns with the anova function.

#End of LRTs

#######################################
######## Missing Data Analysis ########
#######################################

#We conducted logistic regression analyses for the missing data,
#both overall and for each partially observed covariate.
#We fitted logistic regression models with clustered standard errors at the
#highest level of the data structure, CCG

#Create a new dataset for missing data exploration
MAR= GLMMTMBDATA

MAR$Missing = if_else(complete.cases(MAR), 0, 1)
table(MAR$Missing)

#Look at missing overall. The outcome is the count. This will require the rms
#package

#Let's look at association between missingness and the outcome overall. 
MAR_Model_1 = lrm(Missing ~ NumItems, x=T, y=T,data=MAR)
robcov(MAR_Model_1, cluster=MAR$CCGCODE)


#Let's adjust for partially observed covariates and pat1000 (# of patients divided by 1000)
MAR_Model_2 = lrm(Missing ~ NumItems + pats1000 + Urban + CPercFemale + CPerc65plus +yrc+yrc2,x=T, y=T, data=MAR)
robcov(MAR_Model_2, cluster=MAR$CCGCODE)
#p<0.0001 for the outcome, so missingness overall appears to be associated with the outcome
#even after adjusting for all other fully observed covariates.

#Let's look at association between missingness in IMD and the outcome 
MAR$Missing_IMD = if_else(is.na(MAR$IMD), 1, 0)

#Unadjusted with outcome only
MAR_IMDModel_1 = lrm(Missing_IMD ~ NumItems, x=T, y=T,data=MAR)
robcov(MAR_IMDModel_1, cluster=MAR$CCGCODE)

#Adjusted for fully observed covariates
MAR_IMDModel_2 = lrm(Missing_IMD ~ NumItems + pats1000 + Urban + CPercFemale + CPerc65plus +yrc+yrc2,x=T, y=T, data=MAR)
robcov(MAR_IMDModel_2, cluster=MAR$CCGCODE)
#p=0.006 after adjusting for fully observed covariates, so missingness in IMD appears to be associated with the outcome
#even after adjusting for all other fully observed covariates.

#Let's look at association between missingness in % obese and the outcome 
MAR$Missing_OB= if_else(is.na(MAR$CPercOB), 1, 0)

#Unadjusted with outcome only
MAR_OBModel_1 = lrm(Missing_OB ~ NumItems, x=T, y=T,data=MAR)
robcov(MAR_OBModel_1, cluster=MAR$CCGCODE)

#Adjusted for fully observed covariates
MAR_OBModel_2 = lrm(Missing_OB ~ NumItems + pats1000 + Urban + CPercFemale + CPerc65plus +yrc+yrc2,x=T, y=T, data=MAR)
robcov(MAR_OBModel_2, cluster=MAR$CCGCODE)
#p<0.0001 after adjusting for fully observed covariates, so missingness in % obese appears to be associated with the outcome
#even after adjusting for all other fully observed covariates.

#Let's look at association between missingness in % hypertension and the outcome 
MAR$Missing_HYP= if_else(is.na(MAR$CPercHyp), 1, 0)

#Unadjusted with outcome only
MAR_HYPModel_1 = lrm(Missing_HYP ~ NumItems, x=T, y=T,data=MAR)
robcov(MAR_HYPModel_1, cluster=MAR$CCGCODE)

#Adjusted for fully observed covariates
MAR_HYPModel_2 = lrm(Missing_HYP ~ NumItems + pats1000 + Urban + CPercFemale + CPerc65plus +yrc+yrc2,x=T, y=T, data=MAR)
robcov(MAR_HYPModel_2, cluster=MAR$CCGCODE)
#p<0.0001 after adjusting for fully observed covariates, so missingness in % hypertension appears to be associated with the outcome
#even after adjusting for all other fully observed covariates.

#Let's look at association between missingness in % diabetes and the outcome 
MAR$Missing_Dia= if_else(is.na(MAR$CPercDia), 1, 0)

#Unadjusted with outcome only
MAR_DIAModel_1 = lrm(Missing_Dia ~ NumItems, x=T, y=T,data=MAR)
robcov(MAR_DIAModel_1, cluster=MAR$CCGCODE)

#Adjusted for fully observed covariates
MAR_DIAModel_2 = lrm(Missing_Dia ~ NumItems + pats1000 + Urban + CPercFemale + CPerc65plus +yrc+yrc2,x=T, y=T, data=MAR)
robcov(MAR_DIAModel_2, cluster=MAR$CCGCODE)
#p<0.0001 after adjusting for fully observed covariates, so missingness in % diabetes appears to be associated with the outcome
#even after adjusting for all other fully observed covariates.

#Let's look at association between missingness in GP headcount and the outcome 
MAR$Missing_GP= if_else(is.na(MAR$GPHeadcountCentered), 1, 0)

#Unadjusted with outcome only
MAR_GPModel_1 = lrm(Missing_GP ~ NumItems, x=T, y=T,data=MAR)
robcov(MAR_GPModel_1, cluster=MAR$CCGCODE)

#Adjusted for fully observed covariates
MAR_GPModel_2 = lrm(Missing_GP ~ NumItems + pats1000 + Urban + CPercFemale+ CPerc65plus +yrc+yrc2,x=T, y=T, data=MAR, maxit=10000)
robcov(MAR_GPModel_2, cluster=MAR$CCGCODE)
#p=0.007 after adjusting for fully observed covariates, so missingness in GP headcount appears to be associated with the outcome
#even after adjusting for all other fully observed covariates.

#END of all code for the paper.


#Alternative MI approaches for three-level, although, extremely slow to run. 
#adapted from this article:
#https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-020-01079-8#availability-of-data-and-materials

#FCS approach using MICE
IMPDATA <- AnalysisData[,c("CCGCODE", "GPPRACTICECODE", "NumItems","pats1000","yrc", "yrc2", "CPercFemale", "CPerc65plus", "CPercOB", "CPercHyp", "IMD", "Urban", "GPHeadcountCentered", "CPercDia", "NHSRegion")]

#Make factor variables.
IMPDATA$Urban=as.factor(IMPDATA$Urban)
IMPDATA$NHSRegion = as.factor(IMPDATA$Urban)

#Recommend to rescale outcome, or else the code will likely give an error
IMPDATA$NumItems = IMPDATA$NumItems/1000

#Create predictor matrix
pred=make.predictorMatrix(IMPDATA)

#Need to recode some variables to integer and factors
IMPDATA$GPPRACTICECODE = as.integer(as.factor(IMPDATA$GPPRACTICECODE))
IMPDATA$CCGCODE=as.factor(IMPDATA$CCGCODE)

#This identifies the cluster variable
pred[!rownames(pred) %in% c("GPPRACTICECODE"),"GPPRACTICECODE"] = -2

pred[!rownames(pred) %in% c("CPercDia","GPHeadcountCentered", "IMD", "CPercHyp",
                            "CPercOB"),] = 0

meth = rep("", ncol(pred))

#So for any truncated or count data, we will use 2l.pan,
#and no post-imputation round will be used for headcount even
#if below 0 (i.e., -5 for the centred), 
meth[(rownames(pred) %in% c("CPercDia", "CPercHyp",
                            "CPercOB","GPHeadcountCentered"))] <- "2l.pan"

#Since we treat IMD as a factor variable, I decided to use PMM,
#with 10 donors and a default Type 1 match. This should be sufficient for our purposes.
meth[(rownames(pred) %in% c("IMD"))] <- "2l.pmm"

#Set seed
set.seed(1234)
FCS_mice_2L_Long = mice(data=IMPDATA, m=30, method=meth, predictorMatrix = pred, maxit=10, donors=10)
plot(FCS_mice_2L_Long)
stripplot(FCS_mice_2L_Long)

#Write data
#write.mice.imputation(FCS_mice_2L_Long, "Imputations_FCS_2L_Long", include.varnames=TRUE,
#                      long=TRUE, datatype="csv")

#JOMO approach
#Create variables to be imputed
myvars = names(IMPDATA) %in% c("CPercDia", "CPercHyp",
                               "CPercOB","GPHeadcountCentered", "IMD")

dataimp = IMPDATA[myvars]

#Create CCG level indicators
CCG_DI = data.frame(model.matrix(IMPDATA$GPPRACTICECODE ~ as.factor(IMPDATA$CCGCODE)-1, IMPDATA))

names(CCG_DI)[1:ncol(CCG_DI)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep="_"), "CCG_Ind", 191))

#drop one CCG
CCG_DI = CCG_DI[,1:190]

datacomp=cbind(Intercept=rep(1, nrow(IMPDATA)),
               IMPDATA[,names(IMPDATA) %in% c("pats1000", "NumItems", "yrc", "yrc2", "CPercFemale", "CPerc65plus", "Urban","NHSRegion" )],
               CCG_DI)

datacompRE = cbind(Intercept=rep(1, nrow(IMPDATA)))

#Make sure categorical variables are factors
dataimp$IMD = as.factor(dataimp$IMD)
datacomp$Urban = as.factor(datacomp$Urban)

#run imputation
JM_JOMO_2L_Long = jomo(Y=dataimp, X=datacomp, Z=datacompRE, clus=IMPDATA$ORG_CODE, nimp=30, nburn=2000, nbetween=2000)

