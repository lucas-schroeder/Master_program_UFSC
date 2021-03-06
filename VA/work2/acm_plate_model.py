# Math Modules
import numpy as np
import pandas as pd
import scipy as sp
from scipy.sparse.linalg import eigsh

# Plot Libraries
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib import cm

# Utilities
import datetime

pi = np.pi
# %%
class PlateModel:
    def __init__(self, **parameters):
        self.__dict__.update(parameters)
        # cross section moment of inertia per length [m**3]
        self.Iz = (self.h ** 3) / 12

        ## Mesh Definition
        self.nElementX = int(self.L / self.eSize)
        self.nElementY = int(self.H / self.eSize)
        self.nElements = self.nElementX * self.nElementY
        self.dx = self.L / self.nElementX
        self.dy = self.H / self.nElementY
        self.nNodesX = self.nElementX + 1
        self.nNodesY = self.nElementY + 1
        self.nNodes = self.nNodesX * self.nNodesY
        self.ndof = self.nNodes * 3
        self.node_dof = None  # Dof associated with each node
        self.node_coor = None  # Node coordinates dataframe
        self.element_con = None  # Element conectivity dataframe
        self.index_table = None  # DoF index dataframe

        # Global matrices
        self.M = None  # Global mass matrix
        self.K = None  # Global stiffness matrix

        # Boundary conditions
        self.fixed_nodes = None  # List of fixed nodes
        self.fixed_dof = None  # List of restrained DoFs

    def m_e(self):
        a = self.dx / 2
        b = self.dy / 2
        m11 = np.array(
            [
                [3454, 922 * b, -922 * a, 1226, 398 * b, 548 * a],
                [
                    922 * b,
                    320 * b ** 2,
                    -252 * a * b,
                    398 * b,
                    160 * b ** 2,
                    168 * a * b,
                ],
                [
                    -922 * a,
                    -252 * a * b,
                    320 * a ** 2,
                    -548 * a,
                    -168 * a * b,
                    -240 * a ** 2,
                ],
                [1226, 398 * b, -548 * a, 3454, 922 * b, 922 * a],
                [
                    398 * b,
                    160 * b ** 2,
                    -168 * a * b,
                    922 * b,
                    320 * b ** 2,
                    252 * a * b,
                ],
                [
                    548 * a,
                    168 * a * b,
                    -240 * a ** 2,
                    922 * a,
                    252 * a * b,
                    320 * a ** 2,
                ],
            ]
        )

        m21 = np.array(
            [
                [394, 323 * b, -232 * a, 1226, 548 * b, 398 * a],
                [
                    -232 * b,
                    -120 * b ** 2,
                    112 * a * b,
                    -548 * b,
                    -240 * b ** 2,
                    -168 * a * b,
                ],
                [
                    232 * a,
                    112 * a * b,
                    -120 * a ** 2,
                    398 * a,
                    168 * a * b,
                    160 * a ** 2,
                ],
                [1226, 548 * b, -398 * a, 394, 232 * b, 232 * a],
                [
                    -548 * b,
                    -240 * b ** 2,
                    168 * a * b,
                    -232 * b,
                    -120 * b ** 2,
                    -112 * a * b,
                ],
                [
                    -398 * a,
                    -168 * a * b,
                    160 * a ** 2,
                    -232 * a,
                    -112 * a * b,
                    -120 * a ** 2,
                ],
            ]
        )

        m22 = np.array(
            [
                [3454, -922 * b, 922 * a, 1226, -398 * b, -548 * a],
                [
                    -922 * b,
                    320 * b ** 2,
                    -252 * a * b,
                    -398 * b,
                    160 * b ** 2,
                    168 * a * b,
                ],
                [
                    922 * a,
                    -252 * a * b,
                    320 * a ** 2,
                    548 * a,
                    -168 * a * b,
                    -240 * a ** 2,
                ],
                [1226, -398 * b, 548 * a, 3454, -922 * b, -922 * a],
                [
                    -398 * b,
                    160 * b ** 2,
                    -168 * a * b,
                    -922 * b,
                    320 * b ** 2,
                    252 * a * b,
                ],
                [
                    -548 * a,
                    168 * a * b,
                    -240 * a ** 2,
                    -922 * a,
                    252 * a * b,
                    320 * a ** 2,
                ],
            ]
        )

        m = np.concatenate(
            (np.concatenate((m11, m21.T), axis=1), np.concatenate((m21, m22), axis=1)),
            axis=0,
        )

        return (self.rho * self.h * a * b / 6300) * m

    def k_e(self):
        a = self.dx / 2
        b = self.dy / 2
        alp = a / b
        bet = b / a
        nu = self.nu
        I1 = np.eye(3)
        I1[0, 0] = -1
        I2 = np.eye(3)
        I2[1, 1] = -1
        I3 = np.eye(3)
        I3[2, 2] = -1

        k11 = np.array(
            [
                [
                    4 * (bet ** 2 + alp ** 2) + 0.4 * (7 - 2 * nu),
                    2 * (2 * alp ** 2 + 0.2 * (1 + 4 * nu)) * b,
                    2 * (-2 * bet ** 2 - 0.2 * (1 + 4 * nu)) * a,
                ],
                [
                    2 * (2 * alp ** 2 + 0.2 * (1 + 4 * nu)) * b,
                    4 * (4 / 3 * alp ** 2 + 4 / 15 * (1 - nu)) * b ** 2,
                    -4 * nu * a * b,
                ],
                [
                    2 * (-2 * bet ** 2 - 0.2 * (1 + 4 * nu)) * a,
                    -4 * nu * a * b,
                    4 * (4 / 3 * bet ** 2 + 4 / 15 * (1 - nu)) * a ** 2,
                ],
            ]
        )

        k21 = np.array(
            [
                [
                    -(2 * (2 * bet ** 2 - alp ** 2) + 0.4 * (7 - 2 * nu)),
                    2 * (alp ** 2 - 0.2 * (1 + 4 * nu)) * b,
                    2 * (2 * bet ** 2 + 0.2 * (1 - nu)) * a,
                ],
                [
                    2 * (alp ** 2 - 0.2 * (1 + 4 * nu)) * b,
                    4 * (2 / 3 * alp ** 2 - 4 / 15 * (1 - nu)) * b ** 2,
                    0,
                ],
                [
                    -2 * (2 * bet ** 2 + 0.2 * (1 - nu)) * a,
                    0,
                    4 * (2 / 3 * bet ** 2 - 1 / 15 * (1 - nu)) * a ** 2,
                ],
            ]
        )

        k31 = np.array(
            [
                [
                    -(2 * (bet ** 2 + alp ** 2) - 0.4 * (7 - 2 * nu)),
                    2 * (-(alp ** 2) + 0.2 * (1 - nu)) * b,
                    2 * (bet ** 2 - 0.2 * (1 - nu)) * a,
                ],
                [
                    2 * (alp ** 2 - 0.2 * (1 - nu)) * b,
                    4 * (1 / 3 * alp ** 2 + 1 / 15 * (1 - nu)) * b ** 2,
                    0,
                ],
                [
                    2 * (-(bet ** 2) + 0.2 * (1 - nu)) * a,
                    0,
                    4 * (1 / 3 * bet ** 2 + 1 / 15 * (1 - nu)) * a ** 2,
                ],
            ]
        )

        k41 = np.array(
            [
                [
                    2 * (bet ** 2 - 2 * alp ** 2) - 0.4 * (7 - 2 * nu),
                    2 * (-2 * alp ** 2 - 0.2 * (1 - nu)) * b,
                    2 * (-(bet ** 2) + 0.2 * (1 + 4 * nu)) * a,
                ],
                [
                    2 * (2 * alp ** 2 + 0.2 * (1 - nu)) * b,
                    4 * (2 / 3 * alp ** 2 - 1 / 15 * (1 - nu)) * b ** 2,
                    0,
                ],
                [
                    2 * (-(bet ** 2) + 0.2 * (1 + 4 * nu)) * a,
                    0,
                    4 * (2 / 3 * bet ** 2 - 4 / 15 * (1 - nu)) * a ** 2,
                ],
            ]
        )

        k22 = I3.T @ k11 @ I3
        k32 = I3.T @ k41 @ I3
        k42 = I3.T @ k31 @ I3

        k33 = I1.T @ k11 @ I1
        k43 = I1.T @ k21 @ I1

        k44 = I2.T @ k11 @ I2

        aux1 = np.concatenate((k11, k21.T, k31.T, k41.T), axis=1)
        aux2 = np.concatenate((k21, k22, k32.T, k42.T), axis=1)
        aux3 = np.concatenate((k31, k32, k33, k43.T), axis=1)
        aux4 = np.concatenate((k41, k42, k43, k44), axis=1)
        k = np.concatenate((aux1, aux2, aux3, aux4), axis=0)

        return self.E * self.h ** 3 / (48 * (1 - nu ** 2) * a * b) * k

    def generate_mesh(self):
        print("Generating mesh")
        node_dof = np.arange(1, 3 * self.nNodes + 1).reshape((self.nNodes, 3))
        node_dof = pd.DataFrame(
            columns=["w", "thetaX", "thetaY"],
            index=np.arange(1, self.nNodes + 1),
            data=node_dof,
            dtype=int,
        )
        node_dof.index.rename("Node", inplace=True)
        node_dof.rename_axis("DoF", axis="columns", inplace=True)
        self.node_dof = node_dof

        # Nodes coordinates
        node_coor = np.zeros((self.nNodes, 2))
        for n in range(1, self.nNodes + 1):
            q, r = divmod(n - 1, self.nNodesX)  # quotient and remainder
            node_coor[n - 1, 0] = r * self.dx
            node_coor[n - 1, 1] = q * self.dy

        node_coor = pd.DataFrame(
            columns=["x", "y"], index=np.arange(1, self.nNodes + 1), data=node_coor
        )
        node_coor.index.rename("Node", inplace=True)
        node_coor.rename_axis("Coordinate", axis="columns", inplace=True)
        self.node_coor = node_coor

        # Conectivity matrix
        # Relates local nodes of elements to global nodes
        element_con = np.zeros((self.nElements, 4))
        for n in range(1, self.nElements + 1):
            q, r = divmod(n, self.nElementX)  # quotient and remainder
            if r == 0:
                q = q - 1
                r = self.nElementX
            a = q * self.nNodesX + r
            element_con[n - 1, 0] = a
            element_con[n - 1, 1] = a + 1
            element_con[n - 1, 2] = a + self.nNodesX + 1
            element_con[n - 1, 3] = a + self.nNodesX

        element_con = pd.DataFrame(
            columns=["1", "2", "3", "4"],
            index=np.arange(1, self.nElements + 1),
            data=element_con,
            dtype=int,
        )
        element_con.index.rename("Element", inplace=True)
        element_con.rename_axis("Local Node", axis="columns", inplace=True)
        self.element_con = element_con

        # Index table
        # Relates the gloobal DoF associated with each element
        index_table = np.zeros((self.nElements, 12))
        for n in range(1, self.nElements + 1):
            index_table[n - 1, 0] = node_dof["w"][element_con["1"][n]]
            index_table[n - 1, 1] = node_dof["thetaX"][element_con["1"][n]]
            index_table[n - 1, 2] = node_dof["thetaY"][element_con["1"][n]]

            index_table[n - 1, 3] = node_dof["w"][element_con["2"][n]]  # w
            index_table[n - 1, 4] = node_dof["thetaX"][element_con["2"][n]]  # thetaX
            index_table[n - 1, 5] = node_dof["thetaY"][element_con["2"][n]]  # thetaY

            index_table[n - 1, 6] = node_dof["w"][element_con["3"][n]]  # w
            index_table[n - 1, 7] = node_dof["thetaX"][element_con["3"][n]]  # thetaX
            index_table[n - 1, 8] = node_dof["thetaY"][element_con["3"][n]]  # thetaY

            index_table[n - 1, 9] = node_dof["w"][element_con["4"][n]]  # w
            index_table[n - 1, 10] = node_dof["thetaX"][element_con["4"][n]]  # thetaX
            index_table[n - 1, 11] = node_dof["thetaY"][element_con["4"][n]]  # thetaY

        col = pd.MultiIndex.from_arrays(
            [
                ["Local Node 1"] * 3
                + ["Local Node 2"] * 3
                + ["Local Node 3"] * 3
                + ["Local Node 4"] * 3,
                ["w", "thetaX", "thetaY"] * 4,
            ]
        )
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
        K = np.zeros((self.ndof, self.ndof))
        M = np.zeros((self.ndof, self.ndof))

        Ke = self.k_e()
        Me = self.m_e()
        for n in range(1, self.nElements + 1):
            idx = self.index_table.loc[n].to_list()
            for i, p in enumerate(idx):
                for j, q in enumerate(idx):
                    M[p - 1, q - 1] += Me[i, j]
                    K[p - 1, q - 1] += Ke[i, j]
        self.M = M
        self.K = K

    def apply_bc(self, regions: dict):
        print("Applying BC to global matrices")
        fixed_dof = list()
        fixed_nodes = list()
        for r in regions.values():
            fixed_nodes += self.node_coor[
                (self.node_coor["x"].between(r[0][0], r[1][0]))
                & (self.node_coor["y"].between(r[0][1], r[1][1]))
            ].index.tolist()
        self.fixed_nodes = fixed_nodes

        for node in fixed_nodes:
            fixed_dof += self.node_dof.loc[node].tolist()

        fixed_dof = np.array(fixed_dof)
        idx = fixed_dof.argsort()  # Ordering the fixed DoF to removem from [K] and [M]
        fixed_dof = fixed_dof[idx]
        self.fixed_dof = fixed_dof

        # Removing fixed DoF from [K] and [M]
        self.M = np.delete(self.M, fixed_dof - 1, axis=1)
        self.M = np.delete(self.M, fixed_dof - 1, axis=0)
        self.K = np.delete(self.K, fixed_dof - 1, axis=1)
        self.K = np.delete(self.K, fixed_dof - 1, axis=0)

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
        # Make [M] and [K] into an sparse matrix object
        M = sp.sparse.bsr_matrix(self.M)
        K = sp.sparse.bsr_matrix(self.K)

        ## Generalized eigenvalue problem. W is a 1D ndarray and Vc is a 2D ndarray with columns normalized to 1
        # Sove the sparse eig problem using using shift-invert mode
        # looking for the largest shifted eigenvalues (smalest original ones)
        W, Vc = eigsh(A=K, k=40, M=M, which="LM", sigma=0, mode="normal")

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
        results["Vc"] = Vc
        results["V"] = Vc[0::3, :].real  # Displacement shape (dof 1,4,7,10,13,...)
        results["theta"] = np.delete(
            Vc, np.arange(0, Vc.shape[1] - 1, 4), axis=1
        ).real  # Angular shapes (dof 2,3,5,6,8,9,11,12,...)

        # Make the mode shapes have the same orientation to be able to compare between different elements size used
        for j in range(1, W.size):
            if np.sum(results["V"][:, j - 1]) >= 0:
                pass
            else:
                results["V"][:, j - 1] *= -1
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


# %% Auxiliary functions
def A(U, W, j=1, k=1, freq=np.arange(10, 502, 2), eta=0.03, n_modes=40):
    """Compute the accelerance [(m/s^2)/N] given eigenvalues and eigenvectors

    U (ndarray): mode shapes matrix [m]
    W (ndarray): natural frequencies vector [Hz]
    """
    U = U[:, :n_modes]  # crop the matrices to use only n modes
    W = W[:n_modes]
    A_jk = np.zeros(len(freq), dtype=complex)
    for p, f in enumerate(freq):
        for n, fn in enumerate(W):
            A_jk[p] += (
                (-(f ** 2))
                * U[j, n]
                * U[k, n]
                / (fn ** 2 - f ** 2 + 1j * eta * fn ** 2)
            )

    return A_jk, freq
