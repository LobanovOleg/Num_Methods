import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import library.equation as eq
import library.equations_system as eqs

eps = 0.001


def lab02_01():
    l = 1
    r = 2

    def f(x):
        return math.log(x + 1) - x**3 + 1

    def phi(x):
        return (math.log(x + 1) + 1) ** (1.0/3.0)

    def df(x):
        return 1 / (x + 1) - 3 * x**2
    
    
    x_vals = np.linspace(0.5, 2.5, 400)
    y_vals = [f(x) for x in x_vals]

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, label='f(x) = ln(x+1) - x³ + 1')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(l, color='gray', linestyle='--', label=f'Интервал [{l}, {r}]')
    plt.axvline(r, color='gray', linestyle='--')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Графическое определение корня')
    plt.legend()
    plt.grid(True)
    plt.show()

    print('-' * 10, '02-01', '-' * 10, sep='', end='\n\n')
    print('log(x + 1) - x**3 + 1 = 0', end='\n\n')

    x, i = eq.iterations(f, phi, l, r, eps)
    print('Iterations method:')
    if i != -1:
        print('x =', x)
        print('f(x) =', f(x))
        print('Number of iterations:', i)
    else:
        print('Iterations limit exceeded')
    print()

    x, i = eq.newton(f, df, l, r, eps)
    print('Newton method:')
    if i != -1:
        print('x =', x)
        print('f(x) =', f(x))
        print('Number of iterations:', i)
    else:
        print('Iterations limit exceeded')
    print('-' * 25, end="\n\n\n")


def lab02_02():
    a = 4
    l = 1
    r = 1

    def f1(x):
        return a*x[0] - math.cos(x[1])

    def f2(x):
        return a*x[1] - math.exp(x[0])

    def phi1(x):
        return math.cos(x[1]) / a

    def phi2(x):
        return math.exp(x[0]) / a

    def dphi1_dx1(x):
        return 0

    def dphi1_dx2(x):
        return math.sin(x[1]) / a

    def dphi2_dx1(x):
        return math.exp(x[0]) / a

    def dphi2_dx2(x):
        return 0

    def df1_dx1(x):
        return a

    def df1_dx2(x):
        return math.sin(x[1])

    def df2_dx1(x):
        return -math.exp(x[0])

    def df2_dx2(x):
        return a
    
    x2_vals = np.linspace(0, 1.5, 400)
    x1_curve1 = np.cos(x2_vals) / 4

    x1_vals = np.linspace(0, 1, 400)
    x2_curve2 = np.exp(x1_vals) / 4

    plt.figure(figsize=(10, 6))
    plt.plot(x1_curve1, x2_vals, label='ax₁ - cos(x₂) = 0')
    plt.plot(x1_vals, x2_curve2, label='ax₂ - e^x₁ = 0')
    plt.xlabel('x₁')
    plt.ylabel('x₂')
    plt.title('Графическое решение системы')
    plt.legend()
    plt.grid(True)
    plt.show()

    print('-' * 10, '02-02', '-' * 10, sep='', end='\n\n')
    print('Equation system:')
    print('ax1 - cos(x2) = 0')
    print('ax2 - e^x1 = 0', end='\n\n')

    x, i = eqs.iterations(f1, f2, phi1, phi2, dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2, l, r, l, r, eps)
    print('Iterations method:')
    if i != -1:
        print('x1 =', x[0])
        print('x2 =', x[1])
        print('f1(x1, x2) =', f1(x))
        print('f2(x1, x2) =', f2(x))
        print('Number of iterations:', i)
    else:
        print('Iterations limit exceeded')
    print()

    x, i = eqs.newton(f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, l, r, l, r, eps)
    print('Newton method:')
    if i != -1:
        print('x1 =', x[0])
        print('x2 =', x[1])
        print('f1(x1, x2) =', f1(x))
        print('f2(x1, x2) =', f2(x))
        print('Number of iterations:', i)
    else:
        print('Iterations limit exceeded')
    print()
    print('-' * 25, end="\n\n\n")


if __name__ == '__main__':
    if '1' in sys.argv:
        lab02_01()
    if '2' in sys.argv:
        lab02_02()