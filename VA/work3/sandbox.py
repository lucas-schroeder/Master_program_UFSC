#%%
# Math Modules
import numpy as np
import math
import pandas as pd
import scipy as sp
from scipy.misc import derivative
from scipy import integrate
from scipy.sparse.linalg import eigsh

# Plot Libraries
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib import cm

# Utilities
import datetime

pi = np.pi
# %%
def f(x, y):
    return (x ** 2 * y) + (x * y ** 2)


def dbsimpson(f, limits: list, d: list):
    """Simpson's 1/3 rule for double integration

    int_{ay}^{by} int_{ax}^{bx} f(x,y) dxdy

    Args:
        f (func): two variable function, must return float or ndarray
        limits (list): limits of integration [ax, bx, ay, by]
        d (lsit): list of integral resolution [dx, dy]

    Returns:
        float: double integral of f(x,y) between the limits
    """
    ax, bx, ay, by = limits
    dx, dy = d
    nx = math.floor((bx - ax) / dx)
    ny = math.floor((by - ay) / dy)
    s = 0
    for i in range(ny + 1):  # loop of outer integral
        if i == 0 | i == ny:
            p = 1
        elif i % 2 != 0:
            p = 4
        else:
            p = 2

        for j in range(nx + 1):  # loop of inner integral
            if j == 0 | j == nx:
                q = 1
            elif j % 2 != 0:
                q = 4
            else:
                q = 2
            x = ax + j * dx
            y = ay + i * dy
            s += p * q * f(x, y)

    return dx * dy / 9 * s


dbsimpson(f, [1, 2, -1, 1], [0.01, 0.01])


# %%
def f(x, y):
    return (x ** 2 * y) + (x * y ** 2)


def dbsimpson(g: np.ndarray, dxdy: tuple = (1, 1), grid: tuple = None):
    """Simpson's 1/3 rule for double integration

    int_{ay}^{by} int_{ax}^{bx} f(x,y) dxdy
    """
    nx = g.shape[0] - 1
    ny = g.shape[1] - 1

    if grid:
        (x, y) = grid
        ax, bx = np.min(x[1]), np.max(x[1])
        ay, by = np.min(y[:, 0]), np.max(y[:, 0])
        dx = (bx - ax) / nx
        dy = (by - ay) / ny
    else:
        dx, dy = dxdy

    s = 0
    for i in range(ny + 1):  # loop of outer integral
        if i == 0 | i == ny:
            p = 1
        elif i % 2 != 0:
            p = 4
        else:
            p = 2

        for j in range(nx + 1):  # loop of inner integral
            if j == 0 | j == nx:
                q = 1
            elif j % 2 != 0:
                q = 4
            else:
                q = 2
            s += p * q * g[j, i]

    return dx * dy / 9 * s


ax, bx, ay, by = [1, 2, -1, 1]
dx, dy = [0.01, 0.01]
nx = int((bx - ax) / dx)
ny = int((by - ay) / dy)
x = np.arange(ax, bx + dx, dx)
y = np.arange(ay, by + dy, dy)
xv, yv = np.meshgrid(x, y)
g = f(xv, yv)
aa = dbsimpson(g, grid=(xv, yv))
bb = dbsimpson(g, dxdy=(dx, dy))
print(aa, bb)
