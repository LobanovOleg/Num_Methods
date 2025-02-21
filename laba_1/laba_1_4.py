import json
import numpy  as np
import library.matrix as matrix
from library.rotation import rotation

eps = 0.00001
f = open("data/p1_4.json")
data = json.loads(f.read())
X = np.array(data["a"], dtype='float')
eps = 0.00001
values, vectors, iterations = rotation(X, eps)
if iterations == -1:
    print('Матрица не симметрична')
else:
    values_matr = matrix.MyMatrix(values.tolist())
    vectors_matr = matrix.MyMatrix(vectors.tolist())
    print('Собственные значения:')
    print(values_matr)
    print('\nСобственные векторы:')
    print(vectors_matr)
    print('\nКоличество итераций:', iterations)
f.close()
