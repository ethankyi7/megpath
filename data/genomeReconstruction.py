import scanpy as sc
import numpy as np
import pandas as pd
import gzip
import anndata
from scipy.sparse import csr_matrix 

expressionOne = sc.read_10x_mtx(path = 'genePool1', prefix = 'GSE213902_MQpool1_')
expressionTwo = sc.read_10x_mtx(path = 'genePool2', prefix = 'GSE213902_MQpool2_')
###print tests

#clean data

#hto data for pool 1
hto = sc.read_mtx("sra/poolOne/umi_count/matrix.mtx.gz").T
features = []
barcodes = []
with gzip.open("sra/poolOne/umi_count/features.tsv.gz", "rt") as file:
    features = [ln.strip() for ln in file]
mask = [x.lower() != "unmapped" for x in features]
features = [x for x in features if x != "unmapped"]
hto = hto[:, mask].copy()
hto.var_names = features

with gzip.open("sra/poolOne/umi_count/barcodes.tsv.gz", "rt") as file:
    barcodes = [ln.strip() for ln in file]
hto.obs_names = barcodes

expressionOne.obs_names.astype(str)
hto.obs_names.astype(str)

strip = lambda s : s.split('-')[0]
hto.obs_names = pd.Index(hto.obs_names).map(strip)
expressionOne.obs_names = expressionOne.obs_names.map(strip)


common = expressionOne.obs_names.intersection(hto.obs_names)
if len(common) <= 0: print("no overlapping barcodes")

expressionOne = expressionOne[expressionOne.obs_names.isin(common)].copy()
hto = hto[hto.obs_names.isin(common)]
hto = hto[expressionOne.obs_names, :].copy()

hto_cols = []
for tag in hto.var_names:
    col = f"HTO_{tag}"
    expressionOne.obs[col] = np.asarray(hto[:, tag].X.todense()).ravel()
    hto_cols.append(col)

sc.external.pp.hashsolo(expressionOne, cell_hashing_columns=hto_cols)

map_htos = {
    "HTO_HTO1-GTCAACTCTTTAGCG": "W0",
    "HTO_HTO2-TGATGGCCTATTGGG": "W3",
    "HTO_HTO3-TTCCGCCTCTCTTTG": "W6",
    "HTO_HTO4-AGTAAGTTCAGCGTA": "W9",
}

is_singlet = ~expressionOne.obs["Classification"].isin(["Doublet","Negative"])
expressionOne.obs["timepoint"] = expressionOne.obs["Classification"].map(map_htos).where(is_singlet)



#do the same for pool2
htoTwo = sc.read_mtx("sra/poolTwo/umi_count/matrix.mtx.gz").T
features = []
barcodes = []

with gzip.open("sra/poolTwo/umi_count/features.tsv.gz", "rt") as file:
    features = [ln.strip() for ln in file]
mask = [x.lower() != "unmapped" for x in features]
features = [x for x in features if x != "unmapped"]
htoTwo = htoTwo[:, mask].copy()
htoTwo.var_names = features

with gzip.open("sra/poolTwo/umi_count/barcodes.tsv.gz", "rt") as file:
    barcodes = [ln.strip() for ln in file]
htoTwo.obs_names = barcodes

expressionTwo.obs_names.astype(str)
htoTwo.obs_names.astype(str)

htoTwo.obs_names = pd.Index(htoTwo.obs_names).map(strip)
expressionTwo.obs_names = expressionTwo.obs_names.map(strip)

common = expressionTwo.obs_names.intersection(htoTwo.obs_names)
if len(common) <= 0: print("expression for pool two has no overlap with hashings")

expressionTwo = expressionTwo[expressionTwo.obs_names.isin(common)].copy()
htoTwo = htoTwo[htoTwo.obs_names.isin(common)].copy()
htoTwo = htoTwo[expressionTwo.obs_names, :].copy()

hto_cols = []
for tag in hto.var_names:
    col = f"HTO_{tag}"
    expressionTwo.obs[col] = np.asarray(htoTwo[:, tag].X.todense()).ravel()
    hto_cols.append(col)

sc.external.pp.hashsolo(expressionTwo, cell_hashing_columns=hto_cols)

is_singlet = ~expressionTwo.obs["Classification"].isin(["Doublet","Negative"])
expressionTwo.obs["timepoint"] = expressionTwo.obs["Classification"].map(map_htos).where(is_singlet)

print(expressionOne.obs["Classification"].value_counts())
print(expressionOne.obs["timepoint"].value_counts())

print(expressionTwo.obs["Classification"].value_counts())
print(expressionTwo.obs["timepoint"].value_counts())

#filter pool1 for mitochondrial proportion < 12.5% 
expressionOne.var["mt"] = expressionOne.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(expressionOne, qc_vars=["mt"], inplace=True)
expressionOne = expressionOne[expressionOne.obs["pct_counts_mt"] < 12.5, :]
print(expressionOne.obs["pct_counts_mt"])
#-----------
#filter pool2 for mitochondrial proportion < 12.5%
expressionTwo.var["mt"] = expressionTwo.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(expressionTwo, qc_vars=["mt"], inplace=True)
expressionTwo = expressionTwo[expressionTwo.obs["pct_counts_mt"] < 12.5, :]
print(expressionTwo.obs["pct_counts_mt"])
#-----------
expressionTwo.obs_names = [f"{name}_2" for name in expressionTwo.obs_names]
combinedExpression = anndata.concat([expressionOne, expressionTwo], join = "outer", merge="same")
print(combinedExpression)

print(combinedExpression.obs["Classification"].value_counts())

#should be done for the combined expression
sc.pp.filter_genes(data=combinedExpression, min_cells=3)
sc.pp.filter_cells(data=combinedExpression, min_genes=350)

print(combinedExpression.obs["Classification"].value_counts())
print(combinedExpression.obs["timepoint"].value_counts())
# #convert sparse matrix, data into a dense np matrix so we can extract all the values
coefficientMatrix = combinedExpression.X.T.toarray() if not isinstance(combinedExpression, np.ndarray) else combinedExpression.X.T  
# print(type(coefficientMatrix))
# print(expressionOne.X)
# print(coefficientMatrix)

filename = 'pool_output.csv'
#np.savetxt(filename, coefficientMatrix, delimiter=',')
print(coefficientMatrix.shape)