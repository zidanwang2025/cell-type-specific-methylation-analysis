#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyreadr
import pandas as pd
import statsmodels.api as sm


# In[2]:


X = pd.read_csv('X.csv',index_col =0) 
W = pyreadr.read_r('EpiDish_estimated_W.RData')[None]
cov = pyreadr.read_r('cov.RData')[None]


# In[85]:


cov = cov[['AGE_Y15', 'SEX_Y15_num', 'EDUC_Y15', 'MLALC_Y15', 'HEIGHT_Y15',
            'WEIGHT_Y15', 'BMI_Y15']]


# In[3]:


X = X.transpose()


# In[86]:


list_col = []
for i in cov.columns:
    for j in W.columns:
        string = i + "*" + j 
        list_col.append(string)


# In[87]:


product = pd.DataFrame(columns = list_col,index = cov.index)
product


# In[88]:


for i in cov.columns:
    for j in W.columns:
        product[i + "*" + j] = cov[i]* W[j]


# In[89]:


product


# In[90]:


predictors = pd.merge(W,product,left_index=True, right_index=True)


# In[91]:


predictors


# In[92]:


pval_df = pd.DataFrame(columns = X.columns)


# In[93]:


for cpg in X.columns:
    methylation = X.loc[:,cpg]
    model = sm.OLS(methylation, predictors).fit()
    pval_df[cpg] = model.pvalues


# In[94]:


pval_df.transpose().to_csv("pvalues.csv")

