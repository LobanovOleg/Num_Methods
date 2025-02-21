import json
import numpy as np
import library.matrix as matrix
from library.QR import  qr_eigen_values

f = open("data/p1_5.json")
data = json.loads(f.read())
X = np.array(data["a"], dtype='float')
eps = 0.00001
values = qr_eigen_values(X, eps)
values_matr = matrix.MyMatrix(values)
print('Собственные значения:')
print(values_matr)
f.close()
