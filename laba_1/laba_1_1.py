import sys
import json

import library.matrix as matrix

E = matrix.MyMatrix([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])
file = open("data/p1_1.json")
data = json.loads(file.read())
X = matrix.MyMatrix(data["a"])
L, U = X.lu_decompose()
P = matrix.MyMatrix([[1 if j == X.pivot[i] else 0 for j in range(len(X))] for i in range(len(X))])
# PA = P * X
det = X.det()
X_inverse = X.inversed()
res = X.solve_system(data["b"])
b = matrix.MyMatrix([data["b"]]).transposed()
print("LU разложение:")
print("L:", L, sep="\n")
print("U:", U, sep="\n")
print("\nLU проверка разложения: ", end='')
if P.inversed() * L * U == X:
    print("OK")
else:
    print("WRONG")
print("\nДетерминант: ", det)
print("\nОбратная матрица: ", X_inverse, sep="\n")
print("\nПроверка обратной матрицы: ", end='')
if X_inverse * X == E:
    print("OK")
else:
    print("WRONG")
print("\nРешение системы: ", res, sep='\n')
print("\nПроверка решения системы: ", end='')
if X * res == b:
    print("OK")
else:
    print("WRONG")
file.close()