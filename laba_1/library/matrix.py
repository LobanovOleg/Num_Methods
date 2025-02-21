import copy
from numpy import *

ROUND = 6


def scam_func(i, n, l, u):
    for j in range(i, n):
        if u[i][i] != 0:
            l[j][i] = u[j][i] / u[i][i]
        else:
            raise ValueError(f"Нельзя делить на ноль! u[{i}][{i}] = {u[i][i]}")


def scam_func2(i, n, k, l, u):
    for j in range(k - 1, n):
        u[i][j] = u[i][j] - l[i][k - 1] * u[k - 1][j]


def solve_lu(l, u, b):
    """
    solves LUx=b
    """
    n = len(l)
    y = [0 for _ in range(n)]
    for i in range(n):
        s = 0
        for j in range(i):
            s += l[i][j] * y[j]
        y[i] = (b[i] - s) / l[i][i]
    x = [0 for _ in range(n)]
    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(n - 1, i - 1, -1):
            s += u[i][j] * x[j]
        x[i] = (y[i] - s) / u[i][i]
    return x


def decomposite(matrix):
    matrix = array(matrix)
    result = zeros((matrix.shape[0], matrix.shape[1]), dtype='float64')
    result[0] = matrix[0]
    for i in range(1, matrix.shape[0]):
        result[i] = matrix[i] + (- matrix[i][0] / result[0][0]) * result[0]
        result[i][0] = matrix[i][0] / result[0][0]

    for k in range(1, matrix.shape[0] - 1):
        for i in range(k + 1, matrix.shape[0]):
            result[i][k + 1:] = result[i][k + 1:] + (- result[i][k] / result[k][k]) * result[k][k + 1:]
            result[i][k] = result[i][k] / result[k][k]

    return round(result, ROUND)


def solve_lu_one_matrix(lu, b):
    size = lu.shape[0]
    solutions = array([0] * size, dtype='float64')
    solutions[-1] = b[-1] / lu[-1][-1]
    for i in range(size - 2, -1, -1):
        solutions[i] = b[i]
        for j in range(i + 1, size):
            solutions[i] -= lu[i][j] * solutions[j]
        solutions[i] = round(solutions[i] / lu[i][i], ROUND)
    solutions = round(solutions, ROUND)
    return reshape(solutions, (1, solutions.shape[0]))


class MyMatrix:
    # static field for comparing 2 matrix
    __eps = 0.001

    def __init__(self, matrix: list):
        if isinstance(matrix[0], list):
            self.__base_matrix = matrix
            if len(matrix) != 0:
                self.__size_x = len(matrix[0])
                self.__size_y = len(matrix)
            else:
                self.__size_x = 0
                self.__size_y = 0
            return
        else:
            self.__size_x = len(matrix)
            if len(matrix) == 0:
                self.__size_y = 0
                return
            else:
                self.__size_y = 1
                self.__base_matrix = [matrix]
            return

    def __len__(self):
        return self.__size_y

    def __str__(self):
        return '\n'.join([' '.join([str(elem).ljust(25) for elem in row])
                          for row in self.__base_matrix])

    def __copy__(self):
        return MyMatrix(copy.deepcopy(self.__base_matrix))

    def __getitem__(self, index):
        return self.__base_matrix[index]
    
    def __setitem__(self, index, value):
        self.__base_matrix[index] = value

    def __mul__(self, other):
        if self.__size_x != other.__size_y:
            return None
        res_matrix = MyMatrix([[0] * other.__size_x
                               for _ in range(self.__size_y)])
        for i in range(self.__size_y):
            for j in range(other.__size_x):
                for k in range(other.__size_y):
                    res_matrix[i][j] += self[i][k] * other[k][j]
        return res_matrix

    def __eq__(self, other):
        if self.__size_x != other.__size_x or self.__size_y != other.__size_y:
            return False
        else:
            is_eq = True
            for i in range(self.__size_y):
                for j in range(self.__size_x):
                    if not other[i][j] - MyMatrix.__eps <= self[i][j] \
                           <= other[i][j] + MyMatrix.__eps:
                        is_eq = False
        return is_eq

    def transposed(self):
        """
        :return: transposed copy of original matrix without of changing it
        """
        transposed = [[0 for _ in range(self.__size_y)]
                      for _ in range(self.__size_x)]
        for i in range(self.__size_y):
            for j in range(self.__size_x):
                transposed[j][i] = self[i][j]
        return MyMatrix(transposed)

    def transpose(self):
        """
        transposes original matrix
        """
        if self.__size_x == self.__size_y:
            for i in range(self.__size_x):
                for j in range(i, self.__size_y):
                    self[i][j], self[j][i] = self[j][i], self[i][j]
            return
        else:
            new_matrix = self.transposed()
            self.__base_matrix = new_matrix.__base_matrix
            self.__size_x = new_matrix.__size_x
            self.__size_y = new_matrix.__size_y
            return

    def lu_decompose(self):
        n = len(self)
        l = MyMatrix([[0] * n for _ in range(n)])
        u = MyMatrix([[0] * n for _ in range(n)])
        self.pivot = list(range(n))  # Вектор перестановок
        for i in range(n):
            for j in range(n):
                u[i][j] = self[i][j]

        for k in range(n):
            # Частичный выбор ведущего элемента
            max_row = k
            for i in range(k + 1, n):
                if abs(u[i][k]) > abs(u[max_row][k]):
                    max_row = i
            if max_row != k:
                # Перестановка строк в U и L
                u[k], u[max_row] = u[max_row], u[k]
                for i in range(k):
                    l[k][i], l[max_row][i] = l[max_row][i], l[k][i]
                # Обновление вектора перестановок
                self.pivot[k], self.pivot[max_row] = self.pivot[max_row], self.pivot[k]

            if u[k][k] == 0:
                raise ValueError("Матрица вырождена, LU-разложение невозможно.")

            for i in range(k, n):
                l[i][k] = u[i][k] / u[k][k]

            for i in range(k + 1, n):
                for j in range(k, n):
                    u[i][j] -= l[i][k] * u[k][j]

        return l, u

    def det(self):
        l, u = self.lu_decompose()
        num_perm = sum(1 for i in range(len(self)) if self.pivot[i] != i)
        det_sign = (-1) ** num_perm
        res = 1
        for i in range(len(u)):
            res *= u[i][i]
        return det_sign * res

    def inversed(self):
        n = len(self)
        E = MyMatrix([[float(i == j) for i in range(n)] for j in range(n)])
        l, u = self.lu_decompose()
        res = []
        for col in E.transposed():  # Итерация по столбцам E
            pb = [col[self.pivot[i]] for i in range(n)]
            y = solve_lu(l, u, pb)
            res.append(y)
        return MyMatrix(res).transposed()

    def inverse(self):
        """
        inverse the original matrix
        """
        new_matrix = self.inversed()
        self.__base_matrix = new_matrix.__base_matrix
        self.__size_x = new_matrix.__size_x
        self.__size_y = new_matrix.__size_y
        return

    def solve_system(self, b):
        l, u = self.lu_decompose()
        # Применение перестановок к вектору b
        pb = [b[self.pivot[i]] for i in range(len(self))]
        y = solve_lu(l, u, pb)
        return MyMatrix(y).transposed()

    def solve_system_tridiagonal(self, b):
        """
        solves self * x = b, self is tridiagonal matrix
        :param b: free members
        :return: solution in column
        """
        to_use = True
        n = len(self)
        eq = 0
        leq = 0
        for i in range(1, n-1):
            if abs(self[i][1]) < abs(self[i][0]) + abs(self[i][2]) or abs(self[i][1]) < abs(self[i-1][-1]) + abs(self[i+1][0]):
                eq += 1
                to_use = False
            else:
                leq += 1
        to_use *= leq < eq
        p, q = [0] * n, [0] * n
        ans = [0] * n
        p[0] = self[0][1] / -self[0][0]
        q[0] = b[0] / self[0][0]
        for i in range(1, n-1):
            p[i] = -self[i][2] / (self[i][1] + self[i][0]*p[i-1])
            q[i] = (b[i] - self[i][0]*q[i-1]) / (self[i][1] + self[i][0]*p[i-1])
        p[-1] = 0
        q[-1] = (b[-1] - self[-1][0]*q[-2]) / (self[-1][1] + self[-1][0]*p[-2])
        ans[-1] = q[-1]
        for i in range(n-1, 0, -1):
            ans[i-1] = p[i-1] * ans[i] + q[i-1]
        return MyMatrix(ans).transposed(), to_use