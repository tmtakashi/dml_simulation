import numpy as np

def mode_function(x, y, m, n, Lx, Ly):
    w = np.sin(m * np.pi * x / Lx) * np.sin(n * np.pi * y / Ly)
    return w
