import json

import library.matrix as matrix

E = matrix.MyMatrix([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])
file = open("data/p1_1.json")
data = json.loads(file.read())
X = matrix.MyMatrix(data["a"])
L, U, LU = X.lu_decompose()
det = X.det()
X_inverse = X.inversed()
res = X.solve_system(data["b"])
b = matrix.MyMatrix([data["b"]]).transposed()
print("LU разложение:", LU, sep='\n')
print('\n')
print("L:", L, sep="\n")
print('\n')
print("U:", U, sep="\n")
print('\n')
print("\nLU проверка разложения: ", end='')
if L * U == X:
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