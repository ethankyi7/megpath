import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix 

data = sc.read_10x_mtx(path = 'genePool1', prefix = 'GSE213902_MQpool1_')

###print tests

#data = data.transpose()
#print('variables: \n' + data.var)
#print('observations: \n' + data.obs)
#print('matrix: \n')
#print(data.X)
# geneCounts = (data.X > 0).sum(1)
# print(geneCounts[:10])
# print(data.var_names[:20])
#print('CX3CR1' in data.var_names)

#convert sparse matrix, data into a dense np matrix so we can extract all the values
coefficientMatrix = data.X.T.toarray() if not isinstance(data, np.ndarray) else data.X.T  
print(type(coefficientMatrix))
print(data.X)
print(coefficientMatrix)

filename = 'pool1_output.csv'
#np.savetxt(filename, coefficientMatrix, delimiter='\t')
print(coefficientMatrix.shape)