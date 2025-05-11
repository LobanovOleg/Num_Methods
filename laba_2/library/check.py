import numpy as np


def check_phi_derivative(phi, l, r, num_points=100):
    x = np.linspace(l, r, num_points)
    h = 1e-6
    phi_prime = (phi(x + h) - phi(x)) / h
    return np.all(np.abs(phi_prime) < 1)

def check_iterations(f, phi, l, r):
    if not (f(l) * f(r) < 0):
        return False
    try:
        x = np.linspace(l, r, 10)
        y = phi(x)
        if not np.all(np.isfinite(y)):
            return False
    except:
        return False
    if not (np.all(np.diff(y) >= 0) or np.all(np.diff(y) <= 0)):
        return False
    if not (phi(l) - l) * (phi(r) - r) < 0:
        return False
    
    if not check_phi_derivative(phi, l, r):
        return False
    return True


def numerical_derivative(f, x, h=1e-6):
    return (f(x + h) - f(x - h)) / (2 * h)

def numerical_second_derivative(f, x, h=1e-6):
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h ** 2)

def check_newton_convergence_numerical(f, l, r, num_points=100):
    x = np.linspace(l, r, num_points)
    for xi in x:
        df = numerical_derivative(f, xi)
        ddf = numerical_second_derivative(f, xi)
        if abs(f(xi) * ddf) >= (df ** 2):
            return False
    return True

def check_newton(f, df, l, r):
    if not (f(l) * f(r) < 0):
        return False
    continuity = all([abs(f(x) - f(l)) < 1e-6 for x in [l, r]])
    if not continuity:
        return False
    zeros = []
    for x in np.linspace(l, r, 100):
        if abs(df(x)) < 1e-6:
            zeros.append(x)
    if zeros:
        return False
    return True


def check_system(dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2, x=None):
    if x is None:
        x = [1, 1]
    return 1 > max(abs(dphi1_dx1(x)) + abs(dphi1_dx2(x)),
                   abs(dphi2_dx1(x)) + abs(dphi2_dx2(x)))