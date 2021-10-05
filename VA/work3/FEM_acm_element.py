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
from acm_plate_model import A
import datetime

# %% ##########################################################################
# Problem definition: geometry, material, loads and BCs #######################
###############################################################################
geometry = {
    "h": 1.6e-3,  # thickness [m]
    "H": 500e-3,  # door height [m]
    "L": 200e-3,  # door length [m]
}
material = {
    "E": 186e9,  # modulus of elasticity [Pa]
    "nu": 0.3,  # poison coef
    "rho": 7870,  # steel density [kg/m**3]}
}
# Boundary conditions:
# Fixed points delimited by Two points rectangle
region1 = np.array([[0.000, 0.100], [0.010, 0.120]])
region2 = np.array([[0.000, 0.380], [0.010, 0.400]])
region3 = np.array([[0.180, 0.240], [0.200, 0.280]])
regions = {1: region1, 2: region2, 3: region3}

# %% ##########################################################################
# Coarse model ################################################################
###############################################################################
model = {
    "eSize": 0.01,  # element size [m]
}

start = datetime.datetime.now()
door = PlateModel(**{**geometry, **material, **model})
door.generate_mesh()
door.get_global_matrices()
door.apply_bc(regions)
door.solve_eigenvalue_problem()
time_modal = datetime.datetime.now() - start

response_node = door.get_node_nearest(x=0.150, y=0.400)

# %% ##########################################################################
# Finer model #################################################################
###############################################################################
model = {
    "eSize": 0.01 / 2,  # element size [m]
}
start = datetime.datetime.now()
door_fine = PlateModel(**{**geometry, **material, **model})
door_fine.generate_mesh()
door_fine.get_global_matrices()
door_fine.apply_bc(regions)
door_fine.solve_eigenvalue_problem()
time_modal_fine = datetime.datetime.now() - start

response_node_fine = door_fine.get_node_nearest(x=0.150, y=0.400)

# %% ##########################################################################
# Model related plots #########################################################
###############################################################################
# Coarse nodes and BC regions
door.plot_nodes(name="img/nodes_free_and_fixed.pdf")

# Finer
door_fine.plot_nodes(name="img/nodes_free_and_fixed_FINE.pdf")


def c_b(omega, E, h, rho, nu):
    return np.sqrt(omega * np.sqrt(E * h ** 3 / (12 * rho * (1 - nu ** 2))))


phase_speed = c_b(2 * np.pi * 500, 186e9, 1.6e-3, 7870 * 1.6e-3, 0.3)
print(f"Bending wave speed at 500 Hz: {phase_speed} m/s")
print(f"Wavelength at 500 Hz: {phase_speed/(2*np.pi*500)*1000} mm")

# %% ##########################################################################
# Task 1 - Find the first 10 natural frequencies and mode shapes ##############
###############################################################################
x = np.linspace(0, door.L, door.nNodesX)
y = np.linspace(0, door.H, door.nNodesY)
X, Y = np.meshgrid(x, y)
for mode in range(10):  # zero-index
    mode_shape = door.results["V"][:, mode].reshape((door.nNodesY, door.nNodesX))
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    surf = ax.plot_surface(
        X, Y, mode_shape, cmap=cm.coolwarm, linewidth=0, antialiased=True
    )
    ax.set_title(f"Mode {mode+1}")
    ax.set_xlabel("Width [m]")
    ax.set_ylabel("Height [m]")
    ax.set_zticklabels([])
    fig.savefig(f"img/mode_shape_{mode+1}.pdf")

# Print frequencies
pd.options.display.float_format = "{:.2f}".format
df = pd.DataFrame(
    data=[door.results["fn"]],
    index=[f"n={door.nElements}"],
    columns=[f"f_{j}" for j in range(1, len(door.results["fn"]) + 1)],
)
print(df.iloc[:, 0:10])
print(df.iloc[:, 0:10].to_latex())

# %% ##########################################################################
# Task 2a - Accelerance via modal method ######################################
###############################################################################
A_jk, f = A(
    door.results["V"],
    door.results["fn"],
    j=response_node,
    k=response_node,
    n_modes=40,
)
plt.plot(f, np.abs(A_jk), label="Modal Method")

# Plot options
plt.title("Accelerance via Modal Method (first 40 modes)")
plt.yscale("log")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Accelerance [$(m/s^2) / N$]")
plt.tight_layout()
plt.savefig("img/Accelerance_via_Modal_Method_first_40_modes.pdf")
plt.show()

# %% ##########################################################################
# Task 2b - Compare to Direct Method ##########################################
###############################################################################
start = datetime.datetime.now()
A_jk_direct, f = door.frf_direct_method(j=response_node, k=response_node, eta=0.03)
time_direct = datetime.datetime.now() - start

# %% different number of modes
for n in [10, 20, 30, 40]:
    A_jk, f = A(
        door.results["V"],
        door.results["fn"],
        j=response_node,
        k=response_node,
        n_modes=n,
    )
    plt.plot(f, np.abs(A_jk), label=f"{n} modes")

plt.plot(f, np.abs(A_jk_direct), label="Direct Method")

# Plot options
plt.title("Modal vs. Direct Methods")
plt.yscale("log")
plt.legend()
plt.xlabel("Frequency [Hz]")
plt.ylabel("Accelerance [$(m/s^2) / N$]")
plt.tight_layout()
plt.savefig("img/Modal_vs_Direct_Methods.pdf")
plt.show()

# %% Times
tm = {
    "Modal (coarse)": str(time_modal),
    "Modal (fine)": str(time_modal_fine),
    "Direct (coarse)": str(time_direct),
}
times = pd.DataFrame.from_dict(data=tm, columns=["Time [s]"], orient="index")
print(times)
print(times.to_latex())

# %% ##########################################################################
# Task 3a - Compare coarse and fine ###########################################
###############################################################################
A_jk_fine, f = A(
    door_fine.results["V"],
    door_fine.results["fn"],
    j=response_node_fine,
    k=response_node_fine,
    n_modes=40,
)
plt.plot(f, np.abs(A_jk), label="10 mm element")
plt.plot(f, np.abs(A_jk_fine), label="5 mm element")

# Plot options
plt.title("Accelerance, coarse vs. fine mesh")
plt.legend()
plt.yscale("log")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Accelerance [$(m/s^2) / N$]")
plt.tight_layout()
plt.savefig("img/Accelerance_coarse_vs_fine_mesh.pdf")
plt.show()

# %% ##########################################################################
# Get Ansys results on the node of interest ###################################
###############################################################################
path = "ansys_results/"

ansys_fn = np.loadtxt(path + "eigenvalues.txt", delimiter="\t", usecols=2, skiprows=1)
ansys = pd.DataFrame(columns=["Mode", "Node", "x", "y", "z"])

for file in os.listdir(path):
    if "mode_" in file:
        aux = pd.read_csv(
            path + "/" + file, sep="\t", skiprows=1, names=["Node", "x", "y", "z"]
        )
        aux.insert(0, "Mode", [int(file[5:7])])
        ansys = ansys.append(aux)
ansys.set_index("Mode", inplace=True)
ansys["fn"] = ansys_fn.T


# %% ##########################################################################
# Task 4 - Compare with Ansys results #########################################
###############################################################################
A_jk, f = A(door.results["V"], door.results["fn"], j=response_node, k=response_node)
plt.plot(f, np.abs(A_jk), label="FEM - coarse")

A_jk_ansys, f = A(
    np.array(ansys["z"]).reshape((1, -1)), np.array(ansys["fn"]), j=0, k=0
)
plt.plot(f, np.abs(A_jk_ansys), label="Ansys")

# Plot options
plt.title("Code vs. Ansys")
plt.yscale("log")
plt.legend()
plt.xlabel("Frequency [Hz]")
plt.ylabel("Accelerance [$(m/s^2) / N$]")
plt.tight_layout()
plt.savefig("img/Code_vs_Ansys.pdf")
plt.show()

# %% ##########################################################################
# Get all mode shapes from Ansys ##############################################
###############################################################################
path = "ansys_modes/"

ansys_node_coor = pd.read_csv(path + "ansys_nodes_coordinate.csv")
# Reorder acording to our convention
ansys_node_coor = ansys_node_coor.sort_values(by=["y [m]", "x [m]"])
ansys_node_coor.index = np.arange(1, door.nNodes + 1)
idx = ansys_node_coor["node"].to_list()

ansys_modes = pd.DataFrame()
for file in os.listdir(path):
    if "_mode_shape.txt" in file:
        mode = int(file.removesuffix("_mode_shape.txt"))
        aux = pd.read_csv(path + file, sep="\t")
        aux = aux.set_index("Node Number")
        aux = aux.reindex(idx)
        ansys_modes[mode] = aux

ansys_modes = ansys_modes.reindex(sorted(ansys_modes.columns), axis=1)
ansys_modes.index = np.arange(1, door.nNodes + 1)

mac = np.zeros((ansys_modes.shape[1], ansys_modes.shape[1]))
for j, phi1 in enumerate(door.results["V"].T):
    for k, phi2 in enumerate(ansys_modes.to_numpy().T):
        mac[j, k] = np.real(
            np.abs(np.dot(phi1.T, phi2)) ** 2
            / np.abs((np.dot(phi1.T, phi1) * np.dot(phi2.T, phi2)))
        )

# PLot Modal Assurance Criterion - MAC ########################################
_x = np.arange(1, len(mac) + 1)
_y = np.arange(1, len(mac) + 1)
_xx, _yy = np.meshgrid(_x, _y)
x, y = _xx.ravel(), _yy.ravel()
dx = dy = 0.5

fig = plt.figure(figsize=(8, 15))
ax = fig.add_subplot(projection="3d")
colors = cm.rainbow(mac.ravel() / mac.max())
bar = ax.bar3d(x, y, np.zeros_like(mac).ravel(), dx, dy, mac.ravel(), color=colors)
ax.set_title("Modal Assurance Criterion")
ax.set_xlabel("Mode number (code)")
ax.set_ylabel("Mode number (Ansys)")
ax.set_zlabel("MAC")
fig.savefig(f"img/MAC.pdf")


# %% ##########################################################################
# Plot modes for comparison ###################################################
###############################################################################

# Mode 10 - Ansys
fig = plt.figure(figsize=(8, 3))
ax = fig.add_subplot(projection="3d")

mode = 10
mode_shape = ansys_modes.loc[:, mode].to_numpy()
mode_shape = -mode_shape.reshape((door.nNodesY, door.nNodesX))
x = np.linspace(0, door.L, door.nNodesX)
y = np.linspace(0, door.H, door.nNodesY)
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, mode_shape, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_title("Mode 10 - Ansys")
ax.set_xlabel("Width [m]")
ax.set_ylabel("Height [m]")
ax.set_zticklabels([])
fig.savefig(f"img/Mode_10_ansys.pdf")


#  Mode 11 - Ansys
fig = plt.figure(figsize=(8, 3))
ax = fig.add_subplot(projection="3d")

mode = 11
mode_shape = ansys_modes.loc[:, mode].to_numpy()
mode_shape = -mode_shape.reshape((door.nNodesY, door.nNodesX))
x = np.linspace(0, door.L, door.nNodesX)
y = np.linspace(0, door.H, door.nNodesY)
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, -mode_shape, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_title("Mode 11 - Ansys")
ax.set_xlabel("Width [m]")
ax.set_ylabel("Height [m]")
ax.set_zticklabels([])
fig.savefig(f"img/Mode_11_ansys.pdf")

# Mode 10 and 11 from this code
door.plot_mode(10, save_img=True)
door.plot_mode(11, save_img=True)
