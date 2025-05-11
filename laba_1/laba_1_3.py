import json
import numpy as np
import library.matrix as matrix
from library.iteration import solve_system_iterative
from library.zeidel import solve_system_zeidel


f = open("data/p1_3.json")
data = json.loads(f.read())
X = np.array(data["a"], dtype='float')
b = np.array(data["b"], dtype='float')
eps = 0.00000001
print("Метод простых итераций:")
res, n = solve_system_iterative(X, b, eps)
if n == -1:
    print('Невозможно использовать этот метод')
else:
    for x in res:
        print(x)
    print("\nКоличество итераций: ", n, end='\n\n')
    x_matr = matrix.MyMatrix(X.tolist())
    b_matr = matrix.MyMatrix(b.tolist()).transposed()
    res_matr = matrix.MyMatrix(res.tolist()).transposed()
    print("Проверка метода простых итераций: ", end='')
    if x_matr * res_matr == b_matr:
        print("OK")
    else:
        print("WRONG")
    print("\nМетод Зейделя:")
    res, n = solve_system_zeidel(X, b, eps)
    for x in res:
        print(x)
    print("\nКоличество итераций: ", n, end='\n\n')
    x_matr = matrix.MyMatrix(X.tolist())
    b_matr = matrix.MyMatrix(b.tolist()).transposed()
    res_matr = matrix.MyMatrix(res.tolist()).transposed()
    print("Проверка метода Зейделя: ", end='')
    if x_matr * res_matr == b_matr:
        print("OK")
    else:
        print("WRONG")
f.close()