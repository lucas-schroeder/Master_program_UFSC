# Math Modules
import numpy as np
import pandas as pd
import scipy as sp
from scipy.misc import derivative
from scipy import integrate
from scipy.integrate import dblquad
from scipy.integrate import tplquad
from scipy.sparse import bsr_matrix
from scipy.sparse import dok_matrix
from scipy.sparse import lil_matrix
from scipy.sparse.csc import csc_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import eigs

# Plot Libraries
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib import cm

# Utilities
import datetime

pi = np.pi
# %%
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
    nx = int((bx - ax) / dx)
    ny = int((by - ay) / dy)
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


class AcousticModel:
    """Finite Element Analysis of Acoustic Vibrations of Fluids in Cavities"""

    def __init__(self, **parameters):
        # Geometry
        self.W = None  # volume width [m] (X axis)
        self.D = None  # volume depth [m] (Y axis)
        self.H = None  # volume height [m] (Z axis)

        # Fluid
        self.c = None  # sound speed [m/s]
        self.rho = None  # density [kg/m**3]

        # Model
        self.eSize = None  # element size [m]

        # Import geometry, fluid and model parameters as a dict
        self.__dict__.update(parameters)

        ## Mesh Definition
        self.nElementX = int(self.W / self.eSize)
        self.nElementY = int(self.D / self.eSize)
        self.nElementZ = int(self.H / self.eSize)

        self.nElements = self.nElementX * self.nElementY * self.nElementZ

        self.dx = self.W / self.nElementX
        self.dy = self.D / self.nElementY
        self.dz = self.H / self.nElementZ

        self.nNodesX = self.nElementX + 1
        self.nNodesY = self.nElementY + 1
        self.nNodesZ = self.nElementZ + 1

        self.nNodes = self.nNodesX * self.nNodesY * self.nNodesZ
        self.ndof = self.nNodes  # Total number of DoF
        self.node_dof = None  # ID of Dof associated with each node
        self.node_coor = None  # Node coordinates dataframe
        self.element_con = None  # Element conectivity dataframe
        self.index_table = None  # DoF index dataframe

        # Element matrices
        self.H_e = None  # Acoustic stiffness matrix
        self.Q_e = None  # Acoustic inertia matrix

        # Global matrices
        self.M: csc_matrix = None  # Global mass matrix (sparse)
        self.K: csc_matrix = None  # Global stiffness matrix (sparse)

        # Boundary conditions
        self.fixed_nodes = None  # List of fixed nodes
        self.fixed_dof = None  # List of restrained DoFs

    def get_H_e(self):
        a1 = self.dx / 2
        a2 = self.dy / 2
        a3 = self.dz / 2
        Be = np.array(
            [
                [1, -1, -1, -1, 1, 1, 1, -1],
                [1, 1, -1, -1, -1, -1, 1, 1],
                [1, 1, 1, -1, 1, -1, -1, -1],
                [1, -1, 1, -1, -1, 1, -1, 1],
                [1, -1, -1, 1, 1, -1, -1, 1],
                [1, 1, -1, 1, -1, 1, -1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1],
                [1, -1, 1, 1, -1, -1, 1, -1],
            ]
        )
        # Eq. 8.89c,d,e (Fahy), but q is a column vector
        def q1(xi1, xi2, xi3):
            return np.array([0, 1, 0, 0, xi2, xi3, 0, xi2 * xi3]).reshape((8, 1))

        def q2(xi1, xi2, xi3):
            return np.array([0, 0, 1, 0, xi1, 0, xi3, xi1 * xi3]).reshape((8, 1))

        def q3(xi1, xi2, xi3):
            return np.array([0, 0, 0, 1, 0, xi1, xi2, xi1 * xi2]).reshape((8, 1))

        # Eq. 8.95a (Fahy)
        def arg(xi3, xi2, xi1, idx):
            p1 = (1 / a1 ** 2) * (q1(xi1, xi2, xi3) @ q1(xi1, xi2, xi3).T)
            p2 = (1 / a2 ** 2) * (q2(xi1, xi2, xi3) @ q2(xi1, xi2, xi3).T)
            p3 = (1 / a3 ** 2) * (q3(xi1, xi2, xi3) @ q3(xi1, xi2, xi3).T)
            return (p1 + p2 + p3)[idx]

        I = np.zeros((8, 8))
        for i in range(8):
            for j in range(8):
                # Triple integral of arg(), with the fourth input idx=(i,j)
                I[i, j] = tplquad(arg, -1, 1, -1, 1, -1, 1, args=((i, j),))[0]

        self.H_e = (a1 * a2 * a3) * np.linalg.inv(Be.T) @ I @ np.linalg.inv(Be)
        return self.H_e

    def get_Q_e(self):
        a1 = self.dx / 2
        a2 = self.dy / 2
        a3 = self.dz / 2
        Be = np.array(
            [
                [1, -1, -1, -1, 1, 1, 1, -1],
                [1, 1, -1, -1, -1, -1, 1, 1],
                [1, 1, 1, -1, 1, -1, -1, -1],
                [1, -1, 1, -1, -1, 1, -1, 1],
                [1, -1, -1, 1, 1, -1, -1, 1],
                [1, 1, -1, 1, -1, 1, -1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1],
                [1, -1, 1, 1, -1, -1, 1, -1],
            ]
        )
        # Eq. 8.89b (Fahy), but q is a column vector
        def q(xi1, xi2, xi3):
            return np.array(
                [1, xi1, xi2, xi3, xi1 * xi2, xi1 * xi3, xi2 * xi3, xi1 * xi2 * xi3]
            ).reshape((8, 1))

        def arg(xi3, xi2, xi1, idx):
            return (q(xi1, xi2, xi3) @ q(xi1, xi2, xi3).T)[idx]

        I = np.zeros((8, 8))
        for i in range(8):
            for j in range(8):
                I[i, j] = tplquad(arg, -1, 1, -1, 1, -1, 1, args=((i, j),))[0]

        self.Q_e = (
            (a1 * a2 * a3 / self.c ** 2) * np.linalg.inv(Be.T) @ I @ np.linalg.inv(Be)
        )
        return self.Q_e

    def generate_mesh(self):
        print("Generating mesh")
        # Generate the table of DoF IDs and corresponding nodes
        node_dof = np.arange(1, self.nNodes + 1).reshape((self.nNodes, 1))
        node_dof = pd.DataFrame(
            columns=["P"],
            index=np.arange(1, self.nNodes + 1),
            data=node_dof,
            dtype=int,
        )
        node_dof.index.rename("Node", inplace=True)
        node_dof.rename_axis("DoF", axis="columns", inplace=True)
        self.node_dof = node_dof

        # Nodes coordinates
        node_coor = np.zeros((self.nNodes, 3))
        n = 0
        for i in range(self.nNodesX):
            for j in range(self.nNodesY):
                for k in range(self.nNodesZ):
                    node_coor[n, 0] = i * self.dx
                    node_coor[n, 1] = j * self.dy
                    node_coor[n, 2] = k * self.dz
                    n += 1

        node_coor = pd.DataFrame(
            columns=["x", "y", "z"], index=np.arange(1, self.nNodes + 1), data=node_coor
        )
        node_coor.index.rename("Node", inplace=True)
        node_coor.rename_axis("Coordinate", axis="columns", inplace=True)
        self.node_coor = node_coor

        # Conectivity matrix
        # Relates local nodes of elements to global nodes
        element_con = np.zeros((self.nElements, 8))
        m = n = 0
        for k in range(self.nElementZ):
            for j in range(self.nElementY):
                m = k * self.nNodesX * self.nNodesY + j * self.nElementX + 1
                for i in range(self.nElementX):
                    element_con[n, 0] = m
                    element_con[n, 1] = m + 1
                    element_con[n, 2] = m + self.nNodesX + 1
                    element_con[n, 3] = m + self.nNodesX
                    element_con[n, 4] = m + self.nNodesX * self.nNodesY
                    element_con[n, 5] = m + self.nNodesX * self.nNodesY + 1
                    element_con[n, 6] = (
                        m + self.nNodesX * self.nNodesY + self.nNodesX + 1
                    )
                    element_con[n, 7] = m + self.nNodesX * self.nNodesY + self.nNodesX
                    n += 1
                    m += 1

        element_con = pd.DataFrame(
            columns=["1", "2", "3", "4", "5", "6", "7", "8"],
            index=np.arange(1, self.nElements + 1),  # one-index
            data=element_con,
            dtype=int,
        )
        element_con.index.rename("Element", inplace=True)
        element_con.rename_axis("Local Node", axis="columns", inplace=True)
        self.element_con = element_con

        # Index table
        # Relates the gloobal DoF associated with each element
        index_table = np.zeros((self.nElements, 8))
        for n in range(self.nElements):
            index_table[n, 0] = node_dof["P"][element_con["1"][n + 1]]
            index_table[n, 1] = node_dof["P"][element_con["2"][n + 1]]
            index_table[n, 2] = node_dof["P"][element_con["3"][n + 1]]
            index_table[n, 3] = node_dof["P"][element_con["4"][n + 1]]
            index_table[n, 4] = node_dof["P"][element_con["5"][n + 1]]
            index_table[n, 5] = node_dof["P"][element_con["6"][n + 1]]
            index_table[n, 6] = node_dof["P"][element_con["7"][n + 1]]
            index_table[n, 7] = node_dof["P"][element_con["8"][n + 1]]

        col = pd.MultiIndex.from_arrays(
            [
                ["1", "2", "3", "4", "5", "6", "7", "8"],
                ["P"] * 8,
            ]
        )
        col.set_names(["Local Node", "DoF Type"], inplace=True)
        index_table = pd.DataFrame(
            columns=col,
            index=np.arange(1, self.nElements + 1),
            data=index_table,
            dtype=int,
        )
        index_table.index.rename("Element", inplace=True)
        self.index_table = index_table

    def get_global_matrices(self):
        print("Assembling global K and M matrices")
        # LIL is a convenient format for constructing sparse matrices
        K = lil_matrix((self.ndof, self.ndof), dtype=np.float64)
        M = lil_matrix((self.ndof, self.ndof), dtype=np.float64)

        self.get_Q_e()
        self.get_H_e()
        for n in range(1, self.nElements + 1):
            idx = self.index_table.loc[n].to_list()
            for i, p in enumerate(idx):
                for j, q in enumerate(idx):
                    M[p - 1, q - 1] += self.Q_e[i, j]
                    K[p - 1, q - 1] += self.H_e[i, j]

        # Compressed Sparse Column matrix
        self.M = csc_matrix(M)
        self.K = csc_matrix(K)

    def apply_bc(self, regions: dict[float, np.array], tol: list[float] = None):
        """Remove fixed DoFs from global matrices

        Args:
            regions (dict[float, np.array]): coordinates of Two Points Box
                as in {n: np.array([[x1,y1,z1],[x2,y2,z2]])}. All DoF of nodes
                inside the box will be removed from the global system of equations.
            tol (list[float]): tolerance to consider node inside the box. Defaults
                to half the element edge, in each direction.
        """
        print("Applying BC to global matrices")
        tol = tol or [self.dx / 2, self.dy / 2, self.dz / 2]
        fixed_dof = list()
        fixed_nodes = list()
        for r in regions.values():
            idx = (
                (self.node_coor["x"].between(r[0][0] - tol[0], r[1][0] + tol[0]))
                & (self.node_coor["y"].between(r[0][1] - tol[1], r[1][1] + tol[1]))
                & (self.node_coor["z"].between(r[0][2] - tol[2], r[1][2] + tol[2]))
            )
            fixed_nodes += self.node_coor[idx].index.tolist()
        self.fixed_nodes = fixed_nodes

        for node in fixed_nodes:
            fixed_dof += self.node_dof.loc[node].tolist()

        fixed_dof = np.array(fixed_dof)
        idx = fixed_dof.argsort()  # Ordering the fixed DoF to removem from [K] and [M]
        fixed_dof = fixed_dof[idx]
        self.fixed_dof = fixed_dof  # One-indexed ID's of DoF

        # Removing fixed DoF from [K] and [M]
        mask = np.ones(self.M.shape[0], dtype=bool)
        mask[fixed_dof - 1] = False  # create a boolean array
        mask = np.ix_(mask, mask)  # convert to a mesh of booleans
        self.M = self.M[mask]
        self.K = self.K[mask]

    def plot_nodes(self, save=True, name="img/nodes_free_and_fixed.pdf"):
        # Plot nodes
        plt.figure(figsize=(30 / 5, 50 / 5))
        x, y = np.meshgrid(self.node_coor["x"], self.node_coor["y"])
        plt.vlines(x[0], *y[[0, -1], 0], zorder=0, linewidth=0.3, color="tab:gray")
        plt.hlines(y[:, 0], *x[0, [0, -1]], zorder=0, linewidth=0.3, color="tab:gray")
        plt.plot(
            self.node_coor["x"],
            self.node_coor["y"],
            marker="o",
            markersize=4,
            color="tab:blue",
            linestyle="None",
            label="Free node",
        )
        plt.plot(
            self.node_coor["x"][self.fixed_nodes],
            self.node_coor["y"][self.fixed_nodes],
            marker="x",
            markersize=8,
            color="orange",
            linestyle="None",
            label="Fixed node",
        )
        plt.title("Nodes")
        plt.legend(loc="upper center", bbox_to_anchor=(1.2, 1))
        plt.xlabel("X location [m]")
        plt.ylabel("Y location [m]")
        plt.tight_layout()
        if save:
            plt.savefig(name)
        plt.show()
        return plt

    def solve_eigenvalue_problem(self):
        """Find mode shapes and natural frequencies"""
        print("Solving eigenvalue problem")

        ## Generalized eigenvalue problem. W is a 1D ndarray and Vc is a 2D ndarray with columns normalized to 1
        # Sove the sparse eig problem using using shift-invert mode
        # looking for the largest shifted eigenvalues (smalest original ones)
        # W, Vc = eigsh(A=self.K, k=40, M=self.M, which="LM", sigma=0, mode="normal")
        W, Vc = eigs(A=self.M, k=10, M=self.K, which="LM", sigma=0)

        # Ordering eigenvalues and the eigenvectors matrix
        idx = W.argsort()
        W = W[idx]
        Vc = Vc[:, idx]

        # Normalizing eigenvectors matrix by the mass matrix, such that Vc.T @ M @ Vc = I
        m_r = np.diagonal(Vc.T @ M @ Vc)
        m_r = np.reciprocal(np.sqrt(m_r))
        for a in range(Vc.shape[1]):
            Vc[:, a] *= m_r[a]  # multiply every column by the scale factor

        ## Assembling the mode shapes
        # Inserting the restricted DoF from the boundary conditions
        for c in self.fixed_dof:
            Vc = np.insert(Vc, c - 1, 0, axis=0)

        results = dict()
        results["fn"] = (W ** 0.5 / (2 * pi)).real
        results["P"] = Vc.real

        # # Make the mode shapes have the same orientation to be able to compare between different elements size used
        # for j in range(1, W.size):
        #     if np.sum(results["V"][:, j - 1]) >= 0:
        #         pass
        #     else:
        #         results["V"][:, j - 1] *= -1
        self.results = results

    def plot_mode(self, mode=1, interactive=False, save_img=False):
        mode_shape = self.results["V"][:, mode - 1].reshape(
            (self.nNodesY, self.nNodesX)
        )
        if interactive:  # Plot with Plotly
            fig = go.Figure(data=[go.Surface(z=mode_shape)])
            fig.update_layout(
                title="Mode shape",
                autosize=False,
                width=1000,
                height=1000,
                margin=dict(l=65, r=50, b=65, t=90),
            )
            fig.show()

        else:  # Plot with Matplotlib
            fig = plt.figure(figsize=(8, 3))
            ax = fig.add_subplot(projection="3d")

            x = np.linspace(0, self.L, self.nNodesX)
            y = np.linspace(0, self.H, self.nNodesY)
            X, Y = np.meshgrid(x, y)
            # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            # plt.figure(figsize=(8, 15))
            ax.plot_surface(
                X, Y, mode_shape, cmap=cm.coolwarm, linewidth=0, antialiased=True
            )

            ax.set_title(f"Mode {mode}")
            ax.set_xlabel("Width [m]")
            ax.set_ylabel("Height [m]")
            ax.set_zticklabels([])
            if save_img:
                fig.savefig(f"img/Mode_{mode}.pdf")

    def get_node_nearest(self, x, y):
        # Find nearest node within tolerance

        tol = np.sqrt(self.dx * self.dy) / 2  # tolerance
        response_node = self.node_coor[
            (self.node_coor["x"].between(x - tol, x + tol))
            & (self.node_coor["y"].between(y - tol, y + tol))
        ].index.tolist()
        response_node = response_node[0]
        return response_node

    def frf_direct_method(self, j=1, k=1, freq=np.arange(10, 502, 2), eta=0.03, p=39):
        """
        Compute acelerance A_jk

        j - response Node
        k - input Node
        """
        print("Computing FRF via direct method")
        start_time = datetime.datetime.now()

        j = self.node_dof.loc[j, "w"]  # Change from node to DoF
        k = self.node_dof.loc[k, "w"]

        # Number of DoF removed by BC
        j -= self.fixed_dof[self.fixed_dof < j].size
        k -= self.fixed_dof[self.fixed_dof < k].size

        Force = np.zeros((self.M.shape[0], 1))  # Force vector
        Force[k - 1] = 1  # Unity force at DoF k
        X_jk = np.zeros(len(freq), dtype=complex)
        for p, f in enumerate(freq):
            D = (1 + 1j * eta) * self.K - (f * 2 * pi) ** 2 * self.M
            D = sp.sparse.csc_matrix(D)
            # D = sp.sparse.csr_matrix(D)
            R = sp.sparse.linalg.spsolve(D, Force)
            X_jk[p] = R[j - 1]  # Get response at DoF j

        A_jk = -X_jk * (2 * pi * freq) ** 2

        print(f"Diret method elapsed time: {datetime.datetime.now()-start_time}")
        return A_jk, freq
