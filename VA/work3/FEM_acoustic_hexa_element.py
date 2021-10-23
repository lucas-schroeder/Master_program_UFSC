#%%
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm, markers
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import plotly.graph_objects as go
from acm_plate_model import PlateModel
from hexa_acousic_element import AcousticModel
import datetime

# %% ##########################################################################
# Problem definition: geometry, material, loads and BCs #######################
###############################################################################
geometry = {
    "W": 240e-3,  # volume width [m] (X axis)
    "D": 240e-3,  # volume depth [m] (Y axis)
    "H": 540e-3,  # volume height [m] (Z axis)
}
fluid = {  # air @ 20 degrees C
    "c": 348,  # sound speed [m/s]
    "rho": 1.7,  # density [kg/m**3]
}
model = {
    "eSize": 0.02,  # element size [m]
}

volume = AcousticModel(**{**geometry, **fluid, **model})

# %% ##########################################################################
# Generate mesh and assembly global matrices ##################################
###############################################################################
volume.generate_mesh()
start = datetime.datetime.now()
volume.get_global_matrices2()
print(f"Duration: {datetime.datetime.now()-start}")

plt.spy(volume.Hg)

# dok -> brs: 1min 49s ± 1.69 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
# lil -> brs: 43.7 s ± 274 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

# %% ##########################################################################
# Apply BC ####################################################################
###############################################################################
# Regions of zero acoutic potential delimited by Two Point Box
# as in {n: np.array([[x1,y1,z1],[x2,y2,z2]])}
# regions = {
#     1: np.array([[0.000, 0.000, 0.540], [0.240, 0.240, 0.540]]),  # top
#     2: np.array([[0.000, 0.000, 0.000], [0.240, 0.240, 0.000]]),  # bottom
#     3: np.array([[0.000, 0.000, 0.000], [0.000, 0.240, 0.540]]),  # left
#     4: np.array([[0.240, 0.000, 0.000], [0.240, 0.240, 0.540]]),  # right
#     5: np.array([[0.000, 0.240, 0.000], [0.240, 0.240, 0.540]]),  # back
#     6: np.array([[0.000, 0.000, 0.000], [0.240, 0.000, 0.540]]),  # front
# }
# volume.apply_bc(regions)

# %%
start = datetime.datetime.now()
volume.solve_eigenvalue_problem()
print(f"Duration: {datetime.datetime.now()-start}")

volume.results["fn"]
# %%
fig = plt.figure(figsize=(8, 3))
ax = fig.add_subplot(projection="3d")

x = np.linspace(0, volume.W, volume.nNodesX)
y = np.linspace(0, volume.D, volume.nNodesY)
z = np.linspace(0, volume.H, volume.nNodesZ)
X, Y, Z = np.meshgrid(x, y, z)

mode = 3
mode_shape = volume.results["P"][:, mode - 1]
mode_shape = volume.results["P"][:, mode - 1].reshape(
    (volume.nNodesX, volume.nNodesY, volume.nNodesZ)
)
mode_shape = np.moveaxis(mode_shape, 2, 0)
mode_shape = np.moveaxis(mode_shape, 2, 1)
# print(mode_shape.shape)

data_plot = ax.scatter(X, Y, Z, c=mode_shape)
ax.set_xlabel("X Label")
ax.set_ylabel("Y Label")
ax.set_zlabel("Z Label")
fig.colorbar(data_plot)
fig.show()
# %%
