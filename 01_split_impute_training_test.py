import sys,os
import pandas as pd
import math
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import KNNImputer
class KNNImputer_Ks(BaseEstimator, TransformerMixin):
	def __init__(self, *Ks):
		self.Ks = Ks
	def fit(self, X,Ks):
		D_imputer = {}		
		for k in [3,4,5,6,7]:
			imputer = KNNImputer(n_neighbors=k)
			D_imputer[k] = imputer.fit(X)			  
		return D_imputer
	def transform(self, X):
		Impute_train = {}
		for k in [3,4,5,6,7]:
			Impute_train[k] = pd.DataFrame(D_imputer[k].transform(X))
			Impute_train[k].index = X.index
			Impute_train[k].columns = X.columns 
			if k == 3:
				Imputed = Impute_train[k].copy(deep=True)
				Imputed.loc[:,:] = 0
			Imputed = Imputed.add(Impute_train[k],fill_value=0)
		return Imputed/5

feat = pd.read_csv('Other_feat_all_common_specific.txt',header=0,index_col=0,sep='\t')
GO = pd.read_csv('GO_OHE_all_common_specific.txt',header=0,index_col=0,sep='\t')
Domain = pd.read_csv('Domain_OHE_all_common_specific.txt',header=0,index_col=0,sep='\t')
Inter1 = pd.read_csv('Interact_additive_genetic_interaction_OHE_all_common_specific.txt',header=0,index_col=0,sep='\t')
Inter2 = pd.read_csv('Interact_physical_association_OHE_all_common_specific.txt',header=0,index_col=0,sep='\t')
Inter3 = pd.read_csv('Interact_suppressive_genetic_interaction_OHE_all_common_specific.txt',header=0,index_col=0,sep='\t')
Inter4 = pd.read_csv('Interact_synthetic_genetic_interaction_OHE_all_common_specific.txt',header=0,index_col=0,sep='\t')
GO = GO.loc[feat.index.tolist(),:]
Domain = Domain.loc[feat.index,:]
Inter1 = Inter1.loc[feat.index,:]
Inter2 = Inter2.loc[feat.index,:]
Inter3 = Inter3.loc[feat.index,:]
Inter4 = Inter4.loc[feat.index,:]
All = pd.concat([feat,GO,Domain,Inter1,Inter2,Inter3,Inter4],axis=1)
common = pd.read_csv('CommonDefectGenes_10strains_or_greater.txt',header=0,index_col=None,sep='\t')
common['Class'] = 1
common.columns = ['GeneID','Class']
df = pd.read_csv('StrainSpecifcGenes_Defect_3strainsMax.txt',header=0,index_col=None,sep='\t')
specific = pd.DataFrame({'GeneID' : [],'Class':[]})
for i in range(0,df.shape[1]):
	strain = df.columns[i]
	strain_gene = pd.DataFrame(df[strain].dropna(axis=0, how=any))
	strain_gene['Class'] = 0
	strain_gene.columns = ['GeneID','Class']
	specific = pd.concat([specific,strain_gene],axis=0)

specific = specific.drop_duplicates(keep='first')
Genes = pd.concat([common,specific],axis=0)
X = All.loc[Genes['GeneID'].tolist(),:]
y = Genes
y = y.set_index('GeneID')
df = pd.concat([y,X],axis=1)
df.to_csv('Matrix_common_vs_specific.txt', index=True, header=True,sep="\t")
# split the data to training and test
from sklearn.model_selection import train_test_split
train_set, test_set = train_test_split(df, test_size=0.2, random_state=42)
X_train = train_set.drop('Class', axis=1) 
X_test = test_set.drop('Class', axis=1)
y_train = train_set['Class']
y_test = test_set['Class']
# drop features with >50% missing data
Miss_count = X_train.count(0)
Col_to_drop = Miss_count[Miss_count <= 0.5*X_train.shape[0]].index.tolist()
X_train.drop(Col_to_drop,axis=1,inplace=True)
X_test.drop(Col_to_drop,axis=1,inplace=True)
imputer_knn = KNNImputer_Ks()
D_imputer = imputer_knn.fit(X_train, Ks="3,4,5,6,7")
X_train_KNN = imputer_knn.transform(X_train)
X_test_KNN = imputer_knn.transform(X_test)
df_train = pd.concat([y_train,X_train_KNN],axis=1)
df_test = pd.concat([y_test,X_test_KNN],axis=1)
df_train.to_csv('Matrix_common_vs_specific_training_imputed.txt', index=True, header=True,sep="\t")
df_test.to_csv('Matrix_common_vs_specific_test_imputed.txt', index=True, header=True,sep="\t")
