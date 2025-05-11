import library.matrix as matrix
import json
from library.iteration import tridiagonal_to_common

f = open("data/p1_2.json")
data = json.loads(f.read())
X = matrix.MyMatrix(data["a"])
res, _ = X.solve_system_tridiagonal(data["b"])
if res is None:
    print('Невозможно использовать данный метод.')
else:
    b = matrix.MyMatrix([data["b"]]).transposed()
    print("Решение системы:\n", res)
    A = tridiagonal_to_common(X)
    print("\nТрехдиагональная матрица: ", A, sep='\n')
    print("\nПроверка решения системы: ", end='')
    if A * res == b:
        print("OK")
    else:
        print("WRONG")
f.close()