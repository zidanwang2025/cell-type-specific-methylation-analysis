#######################################################################################################################
## Last updated: 05/19/2022
## Author: Zidan
## Running the code should generate and save: 
## 1. the EpiDISH estimated W matrix
## 2. the pre-processed X matrix that matches the TCA results 
## 3. the covariate matrix used in the TCA model 
## 4. the TCA results 
## 5. the estimated Z matrix from tensor function 
## Note: 1. takes about 2-3 days to run on NU Quest 2. Used a little more than 50 GB (update after this run)
#######################################################################################################################

Packages <- c("TCA","ggplot2","ggpubr","pracma","matrixStats")
lapply(Packages, library, character.only = TRUE)
source("https://raw.githubusercontent.com/cozygene/TCA/master/vignettes/vignette_analysis.R")

CARDIA_Y15 <-readRDS("/projects/b1096/CARDIA/CARDIA_Y15_Combined_preprocessed_NoRCP_KNNimpute.rds")
phenotype<- readRDS("/projects/b1096/CARDIA/Basic_Phenotype_Data/CARDIA_Y15_Methy_1089_Info.rds")

library(EpiDISH)
out <- epidish(beta.m = CARDIA_Y15, ref.m = centDHSbloodDMC.m, method = "RPC") 
saveRDS(out$estF, file="EpiDish_estimated_W.RData",compress = "gzip")
cell_p = out$estF 

#Make sex into dummies
phenotype$SEX_Y15_num=0
phenotype$SEX_Y15_num[phenotype$SEX_Y15 == "Male"] <- 1
phenotype_use = phenotype[,c("AGE_Y15","SEX_Y15_num","EDUC_Y15","MLALC_Y15","BMI_Y15")]
rownames(phenotype_use) = phenotype[,"BARCODE"]
has.neg <- apply(phenotype_use, 1, function(row) any(row < 0))
phenotype_use = phenotype_use[!has.neg,]
phenotype_use = na.omit(phenotype_use)

#Make sure that X, W, and cov matrices have the same sizes and are in the same order (required by TCA)
ind = intersect(colnames(CARDIA_Y15),rownames(phenotype_use))
CARDIA_Y15 = CARDIA_Y15[,colnames(CARDIA_Y15) %in% ind] 
phenotype_use = phenotype_use[rownames(phenotype_use) %in% ind,] 
cell_p = cell_p[rownames(cell_p) %in% ind,] 
phenotype_use = phenotype_use[match(colnames(CARDIA_Y15),rownames(phenotype_use)),] #match the colnames of X and rownames of cov

Cardia<- list()
Cardia$W = cell_p 
Cardia$X = CARDIA_Y15
Cardia$cov = phenotype_use

saveRDS(Cardia$cov, file="cov.RData",compress = "gzip")
saveRDS(Cardia$X, file="X.RData",compress = "gzip")

refactor.mdl.cardia  <- refactor(X = Cardia$X,
                                 k = 7)

tca.mdl.cardia = tca(X = Cardia$X,
                     W = Cardia$W,
                     C1 = Cardia$cov[,c("AGE_Y15","SEX_Y15_num","EDUC_Y15","MLALC_Y15","BMI_Y15")],
                     C2 = refactor.mdl.cardia$scores,
                     refit_W = TRUE)

saveRDS(tca.mdl.cardia, file="TCA_results.RData",compress = "gzip")

tca.tensor =tensor(Cardia$X,
          tca.mdl.cardia,
          scale = FALSE,
          parallel = TRUE,
          num_cores = NULL,
          log_file = "tensor.log",
          debug = TRUE,
          verbose = TRUE)
saveRDS(tca.tensor, file="Tensor_Z.RData",compress = "gzip")