def df(x, y, x_):
    assert len(x) == len(y)
    for interval in range(len(x) - 1):
        if x[interval] <= x_ <= x[interval + 1]:
            i = interval
            break
    
    first = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
    
    a1 = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
    a2 = (((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) - a1) / (x[i + 2] - x[i])) * (2 * x_ - x[i] - x[i + 1])
    second = a1 + a2
    
    return first, second


def d2f(x, y, x_):
    assert len(x) == len(y)
    for interval in range(len(x)):
        if x[interval] <= x_ < x[interval+1]:
            i = interval
            break

    num = (y[i+2] - y[i+1]) / (x[i+2] - x[i+1]) - (y[i+1] - y[i]) / (x[i+1] - x[i])
    return 2 * num / (x[i+2] - x[i])