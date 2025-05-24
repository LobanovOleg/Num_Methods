def rectangle(f, l, r, h):
    if l > r:
        return None
    n = int((r - l) / h)
    result = 0
    for i in range(n):
        x_left = l + i * h
        x_right = x_left + h
        result += h * f((x_left + x_right) / 2)
    return result

def trapeze(f, l, r, h):
    if l > r:
        return None
    if h <= 0:
        return None
    n = int((r - l) / h)
    if n == 0:
        return 0.0
    
    total = (f(l) + f(r)) / 2
    for i in range(1, n):
        x = l + i * h
        total += f(x)
    total *= h
    return total


def simpson(f, l, r, h):
    if l > r:
        return None
    n = int((r - l) / h)
    if n % 2 != 0:
        n += 1
    h = (r - l) / n
    result = 0
    for i in range(1, n // 2 + 1):
        x1 = l + 2 * (i - 1) * h
        x2 = x1 + h
        x3 = x1 + 2 * h
        result += (f(x1) + 4 * f(x2) + f(x3)) * h
    return result / 3.0


def runge_rombert(h1, h2, i1, i2, p):
    return i1 + (i1 - i2) / ((h2 / h1)**p - 1)