import json
from library.rotation import rotation
from library.iteration import MyMatrix

eps = 0.00000001
f = open("data/p1_4.json")
data = json.loads(f.read())
X = data["a"]
eps = 0.00001
values, vectors, iterations = rotation(X, eps)
if iterations == -1:
    print('Матрица не симметрична')
else:
    values_matr = MyMatrix([[val] for val in values])
    vectors_matr = MyMatrix(vectors)
    print('Собственные значения:')
    print(values_matr)
    print('\nСобственные векторы:')
    print(vectors_matr)
    print('\nКоличество итераций:', iterations)
f.close()
