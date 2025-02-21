import library.matrix as matrix
import json
from library.progon import tridiagonal_to_common

f = open("data/p1_2.json")
data = json.loads(f.read())
X = matrix.MyMatrix(data["a"])
res, _ = X.solve_system_tridiagonal(data["b"])
if res is None:
    print('Unable to use this method')
else:
    b = matrix.MyMatrix([data["b"]]).transposed()
    print("Решение системы:\n", res)
    print("\nПроверка решения системы: ", end='')
    X = tridiagonal_to_common(X)
    if X * res == b:
        print("OK")
    else:
        print("WRONG")
f.close()