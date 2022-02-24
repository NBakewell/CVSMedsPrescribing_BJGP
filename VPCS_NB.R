#This is a function to derive the Variance Partition Coefficients (VPCs) 
#for an adjusted three-level negative binomial multi-level model;
#code adpated from Leckie et al., http://www.bristol.ac.uk/cmm/media/leckie/articles/leckie2020.pdf

VPCs_NB = function(Model.Name) {
        #Variance Partition Coefficients (VPCs)
        xb_adj_CC <- predict(Model.Name)
        xb_adj_CC <- data.frame(xb_adj_CC)
        
        #Level-3 (CCG) variance
        sigma_CCG_adj <- summary(Model.Name)$varcor$cond$CCGCODE[1,1]
        
        #Level-2 (GPP) variance
        sigma_GPP_adj <- summary(Model.Name)$varcor$cond$GPPRACTICECODE[1,1]
        
        # Over-dispersion parameter
        alpha_adj <- 1/(summary(Model.Name)$sigma)
        
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
        
        t(as.data.frame(sapply(xb_adj_CC, mean)[8:10]))
}
