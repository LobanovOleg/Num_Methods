import math

def find_max_upper_element(X):
    """
    :return: coords of the max element above the main diagonal
    """
    n = len(X)
    i_max, j_max = 0, 1
    max_elem = abs(X[0][1])
    for i in range(n):
        for j in range(i + 1, n):
            if abs(X[i][j]) > max_elem:
                max_elem = abs(X[i][j])
                i_max = i
                j_max = j
    return i_max, j_max


def matrix_norm(X):
    """
    :return: L2 norm for elements above the main diagonal
    """
    n = len(X)
    norm = 0
    for i in range(n):
        for j in range(i + 1, n):
            norm += X[i][j]**2
    return math.sqrt(norm)


def rotation(A, eps):
    """
    :return: eigen values, eigen vectors, number of iterations
    """
    n = len(A)
    for i in range(n):
        for j in range(i, n):
            if A[i][j] != A[j][i]:
                return None, None, -1

    A_i = [row[:] for row in A]
    eigen_vectors = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    iterations = 0
    while matrix_norm(A_i) > eps:
        i_max, j_max = find_max_upper_element(A_i)
        if A_i[i_max][i_max] - A_i[j_max][j_max] == 0:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(2 * A_i[i_max][j_max] / (A_i[i_max][i_max] -
                                                           A_i[j_max][j_max]))
        
        U = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        U[i_max][j_max] = -math.sin(phi)
        U[j_max][i_max] = math.sin(phi)
        U[i_max][i_max] = math.cos(phi)
        U[j_max][j_max] = math.cos(phi)

        temp = [[0.0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    temp[i][j] += U[k][i] * A_i[k][j]
        
        A_i_new = [[0.0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    A_i_new[i][j] += temp[i][k] * U[k][j]
        A_i = A_i_new

        temp_vectors = [[0.0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    temp_vectors[i][j] += eigen_vectors[i][k] * U[k][j]
        eigen_vectors = temp_vectors

        iterations += 1
    eigen_values = [A_i[i][i] for i in range(n)]
    return eigen_values, eigen_vectors, iterations