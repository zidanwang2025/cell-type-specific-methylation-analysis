library(qvalue)
library(ggplot2)
library(ggpubr)
library(pracma)
library(matrixStats)

tca_one_hits<- function(TCA,var_interest = colnames(TCA$C1)[1],cell_type=colnames(TCA$W)[1],correction="bonferroni",test="marginal conditional",cell_type_specific = TRUE, signif_level=0.05){
  ##############################################################################################################################
  ####Parameters of the function
  #### TCA: tca results from running tca
  #### var_interest: outcome of interest, e.g. age, sex, etc., default value is the first column name of C1 (covariate matrix)
  #### cell_type: cell type of interest, e.g. CD4, CD8, etc., default value is the first column name of W (weight matrix)
  #### cell_type (cont): only needed when test == "marginal conditional" and cell_type_specific == TRUE
  #### correction: type of multiple testing correction used, can choose between "FDR" or "bonferroni" or "hierarchical"
  #### test: type of test interested, can choose between "marginal conditional" or "joint"
  #### cell_type_specific: if TRUE, return cell-type specific results of marginal conditional test, 
  #### cell_type_specific (cont): if FALSE, return hits of all cell types under marginal conditional test
  #### cell_type_specific (cont): only needed when test == "marginal conditional"
  #### signif_level: significance level of tests, default is 0.05
  ##############################################################################################################################
  #store joint p-values for variable of interest
  if (test == "joint"){
    p_joint = TCA$gammas_hat_pvals.joint[,var_interest]
    if (correction == "bonferroni"){
      ### Exact the hits based on the joint test (Using Bonferroni Correction)
      hits.joint.bonferroni = names(which(p_joint < signif_level/nrow(TCA$mus_hat)))
      return(hits.joint.bonferroni)
    }
    else if (correction == "FDR"){
      qobj = qvalue(p=p_joint, fdr.level=signif_level)
      hits.joint.fdr = names(qobj$qvalues[qobj$qvalues<signif_level])
      return(hits.joint.fdr)
    }
  }
  else if (test == "marginal conditional"){
    pvals.marg_cond = TCA$gammas_hat_pvals[,paste(colnames(TCA$W),".",var_interest,sep="")]
    col_name = paste0(cell_type,".",var_interest)
    if (cell_type_specific == TRUE){
      if (correction == "bonferroni"){
        hits.mc.bonferroni= rownames(pvals.marg_cond[pvals.marg_cond[,col_name]<= signif_level/(nrow(TCA$mus_hat)*nrow(TCA$W)),])
        return(hits.mc.bonferroni)
      }
      else if (correction == "FDR"){
        qobj= qvalue(p=pvals.marg_cond, fdr.level=signif_level)
        hits.mc.fdr = names(which(qobj$qvalues[,col_name]<signif_level))
        return(hits.mc.fdr)
      }
      else if (correction == "hierarchical"){
        #First step: get row minimum
        min.pvals = rowMins(pvals.marg_cond)
        pvals =cbind(pvals.marg_cond,min.pvals)
        #Second step: FDR 
        qobj=qvalue(p=pvals[,'min.pvals'], fdr.level=signif_level, pi0 = 1) #since all p-values are significant, not specifying pi0 will throw an error
        hits.mc.first = names(qobj$qvalues[qobj$qvalues<signif_level])
        #Third step: get p-values of the subset of the CPG sites
        p.subset = pvals.marg_cond[hits.mc.first,]
        qobj.second= qvalue(p=p.subset, fdr.level=signif_level)
        hits.mc.fdr.second = names(qobj.second$qvalues[qobj.second$qvalues[,col_name]<signif_level])
        return(hits.mc.fdr.second)
      }
    }
    else {
      if (correction == "bonferroni"){
        hits.mc.bonferroni= rownames(pvals.marg_cond[rowSums(pvals.marg_cond<= signif_level/(nrow(TCA$mus_hat)*nrow(TCA$W)))!=0,])
        return(hits.mc.bonferroni)
      }
      else if (correction == "FDR"){
        qobj= qvalue(p=pvals.marg_cond, fdr.level=signif_level)
        hits.mc.fdr = rownames(qobj$qvalues[rowSums(qobj$qvalues<= signif_level)!=0,])
        return(hits.mc.fdr)
      }
      else if (correction == "hierarchical"){
        #First step: get row minimum
        min.pvals = rowMins(pvals.marg_cond)
        pvals =cbind(pvals.marg_cond,min.pvals)
        #Second step: FDR 
        qobj=qvalue(p=pvals[,'min.pvals'], fdr.level=signif_level, pi0 = 1) #since all p-values are significant, not specifying pi0 will throw an error
        hits.mc.first = names(qobj$qvalues[qobj$qvalues<signif_level])
        #Third step: get p-values of the subset of the CPG sites
        p.subset = pvals.marg_cond[hits.mc.first,]
        qobj.second= qvalue(p=p.subset, fdr.level=signif_level)
        hits.mc.fdr.second = names(qobj.second$qvalues[qobj.second$qvalues<signif_level])
        return(hits.mc.fdr.second)
      }
    }
  }
}


no_tca_three_hits<- function(pvals,var_interest = NULL,cell_type=NULL,correction="bonferroni",signif_level=0.05){
  ##############################################################################################################################
  ####Parameters of the function
  #### pvals: pvals from running linear regression
  #### var_interest: outcome of interest, e.g. age, sex, etc., default value is the first column name of C1 (covariate matrix)
  #### cell_type: cell type of interest, e.g. CD4, CD8, etc., default value is the first column name of W (weight matrix)
  #### cell_type (cont): only needed when test == "marginal conditional" and cell_type_specific == TRUE
  #### correction: type of multiple testing correction used, can choose between "FDR" or "bonferroni" or "hierarchical"
  #### signif_level: significance level of tests, default is 0.05
  ##############################################################################################################################
  col_name=paste0(var_interest,".",cell_type)
  if (correction=="bonferroni"){
    p_val = signif_level/(nrow(pvals)*ncol(pvals))
    return(rownames(pvals[which(pvals[,col_name]<=p_val),]))
  }
  else if (correction=="FDR"){
    qobj= qvalue(p=pvals, fdr.level=signif_level)
    return(rownames(qobj$qvalues[which(qobj$qvalues[,col_name]<=signif_level),]))
  }
}
  