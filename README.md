# Cell Type Specific Methylation Analysis

## Overview 
The DNA samples most readily available for methylation studies are derived from whole blood. However, blood consists of many functionally and developmentally distinct cell populations in varying proportions. High costs and technical limitations in techniques such as cell sorting and single-cell restrict the collection of large-scale, cell-type-specific DNA methylation data, which impedes our ability to tackle some important biological questions, such as the identification of disease-associated genes at a cell-type-specific resolution.  

The question we asked is: how do we conduct cell-type specific analyses using bulk methylation data? 

This project aims to develop and evaluate different approaches that allow us to conduct cell-type-specific analyses from bulk methylation data. We mainly considered three approaches:

***Approach 1:*** Directly using TCA (Tensor Composition Analysis). This can take very long even on HPC.
Fitting the TCA model using the ```"tca"``` function requires two arguments - an X matrix, which is a tissue-level matrix of methylation values (methylation sites by individuals), and a W matrix, which is a cell-type proportion for the individuals in X. Since the original data does not include the information on the cell-type proportions for the individuals in X, we first estimated the W matrix using EpiDish. 

***Approach 2:*** (a two-stage approach) Let $Z_{ijl}$ denote the methylation level for the $l$-th cell type at site CpG $j$ of individual $i$. Then we can first use TCA to estimate $Z_{ijl}$, Then we can use conduct a cross-sectional analysis using a fixed-effects model that models $Z_{ijl}$ in terms of covariates. 


***Approach 3:*** (a one-step approach) For this one-step approach, we did not estimate the values of $Z_ijl$. Instead, we considered the methylation level as a weighted average of cell-type-specific methylation levels, i.e. $E(M_{i j 0})=\sum_{l=1,...,k} w_{i l} Z_{i j l}$, where $w_{i l}$ is the estimated cell-type proportion for the $l$-th cell type of individual $i$. 

**Step 1:**

Use ```TCA_generate. R``` & ```run_tca.sh``` generate:

1. ```EpiDish_estimated_W.RData```

2. ```cov.RData``` (Covariates, age, gender, edu, mlalc, bmi)_

3. ```X.RData```

4. ```TCA_results.RData```

5. ```Tensor_Z.RData```

**Step 2:**
For ***Apporach 1:*** use 4. ```TCA_results.RData```, run function ```tca_one_hits()```
For  ***Apporach 2:*** use 5. ```Tensor_Z.RData```, run fixed effect model
For  ***Apporach 3:*** 1) use 1. ```EpiDish_estimated_W.RData```, 2. ```cov.RData```, and 3. ```X.RData```, run ```Approach3_Regression.py```, generate pvalues.csv. then 2) run ```tca_three_hits()```


