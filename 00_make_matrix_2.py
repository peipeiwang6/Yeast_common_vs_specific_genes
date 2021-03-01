import sys,os
import pandas as pd
import math
import numpy as np
from sklearn.preprocessing import MultiLabelBinarizer
mlb = MultiLabelBinarizer()
common = pd.read_csv('CommonDefectGenes_10strains_or_greater.txt',header=0,index_col=None,sep='\t')
common['Class'] = 1
common.columns = ['GeneID','Class']
df = pd.read_csv('StrainSpecifcGenes_Defect_3strainsMax.txt',header=0,index_col=None,sep='\t')

feat = pd.read_csv('All_features_except_OHE_features_20191028.txt',header=0,index_col=0,sep='\t')
Domain = open('domains.tab','r').readlines()
GO = open('GO_Annotations_from_Dee.txt','r').readlines()
Inter = pd.read_csv('Interaction_features.txt', sep='\t', index_col = 0, header = 0)

specific = pd.DataFrame({'GeneID' : [],'Class':[]})
for i in range(0,df.shape[1]):
	strain = df.columns[i]
	strain_gene = pd.DataFrame(df[strain].dropna(axis=0, how=any))
	strain_gene['Class'] = 0
	strain_gene.columns = ['GeneID','Class']
	specific = pd.concat([specific,strain_gene],axis=0)

specific = specific.drop_duplicates(keep='first')
Genes = pd.concat([common,specific],axis=0)
features = feat.loc[Genes['GeneID'].tolist(),:]
# for domains
D = {} # list of total domains
DOMAIN = {} # Domain[gene] = [domain1, domain2,...]
for inl in Domain:
	tem = inl.split('\t')
	if tem[0] in Genes['GeneID'].tolist():
		D[tem[4]] = 1
		if tem[0] not in DOMAIN:
			DOMAIN[tem[0]] = []
		if tem[4] not in DOMAIN[tem[0]]:
			DOMAIN[tem[0]].append(tem[4])

X = []
for gene in Genes['GeneID'].tolist():
	out = []
	if gene in DOMAIN:
		for i in range(0,len(DOMAIN[gene])):
			out.append(DOMAIN[gene][i])
	X.append(out)

X_OHE = mlb.fit_transform(X)
Domain_matrix = pd.DataFrame(X_OHE)
Domain_matrix.index = Genes['GeneID'].tolist()
Domain_matrix.columns = mlb.classes_
Domain_matrix.to_csv('Domain_OHE_all_common_specific.txt', index=True, header=True,sep="\t")

# for GO
GO_anno = {} # GO[gene] = [GO1, GO2,...]
for inl in GO:
	if inl.startswith('GO:'):
		tem = inl.split('\t')
		genes = tem[3].strip().replace('\"',"")
		genes = genes.split(',')
		for gene in genes:
			if gene in Genes['GeneID'].tolist():
				if gene not in GO_anno:
					GO_anno[gene] = []
				GO_anno[gene].append(tem[0])

X = []
for gene in Genes['GeneID'].tolist():
	out = []
	if gene in GO_anno:
		for i in range(0,len(GO_anno[gene])):
			out.append(GO_anno[gene][i])
	X.append(out)

X_OHE = mlb.fit_transform(X)
GO_matrix = pd.DataFrame(X_OHE)
GO_matrix.index = Genes['GeneID'].tolist()
GO_matrix.columns = mlb.classes_
GO_matrix.to_csv('GO_OHE_all_common_specific.txt', index=True, header=True,sep="\t")

# for interaction
I = {}
Inter = Inter.loc[Genes['GeneID'].tolist(),:]
colnames = Inter.columns.tolist()
rownames = Inter.index.tolist()
for col in colnames:
	I[col + '_numbers'] = {}
	for gene in rownames:
		if gene in Genes['GeneID'].tolist():
			try:
				if not math.isnan(Inter.loc[gene,col]):
					number = len(Inter.loc[gene,col].split(','))
					I[col + '_numbers'][gene] = number
				else:
					I[col + '_numbers'][gene] = 0
			except:
				number = len(Inter.loc[gene,col].split(','))
				I[col + '_numbers'][gene] = number

PN = pd.DataFrame.from_dict(I, orient='columns')
PN = PN.loc[Genes['GeneID'].tolist(),:]
features = pd.concat([features,PN],axis=1)
features.to_csv('Other_feat_all_common_specific.txt', index=True, header=True,sep="\t")

for col in colnames:
	X = []
	for gene in Genes['GeneID'].tolist():
		out = []
		if gene in Inter.index:
			try:
				for p in Inter.loc[gene,col].split(','):
					out.append(p)
			except:
				if not math.isnan(Inter.loc[gene,col]):
					print(gene)
		X.append(out)
	X_OHE = mlb.fit_transform(X)
	Inter_matrix = pd.DataFrame(X_OHE)
	Inter_matrix.index = Genes['GeneID'].tolist()
	Inter_matrix.columns = mlb.classes_
	Inter_matrix.to_csv('Interact_%s_OHE_all_common_specific.txt'%col, index=True, header=True,sep="\t")
