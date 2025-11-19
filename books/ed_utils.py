from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import glob

# --- Finding GS ---
def exact_diagonalization_line(values, scan_var, opp, set, basis, L, J=1):
    '''Perform exact diagonalization across a line on the parameter space of your phase diagram'''
    if not os.path.exists("ed_results"):
        os.mkdir("ed_results")

    folder_eigenstates = f"ed_results/L{L}_{opp}{set}_{scan_var}{np.min(values)}-{np.max(values)}_states"
    folder_eigenvalues = f"ed_results/L{L}_{opp}{set}_{scan_var}{np.min(values)}-{np.max(values)}_values"

    if not os.path.exists(folder_eigenstates):
        os.mkdir(folder_eigenstates)
        os.mkdir(folder_eigenvalues)
    else:
        return 1
    
    for ind, val in tqdm(enumerate(values)):
        # determine h and k
        if scan_var == "h":
            h = val
            k_l = k_r = set
        else:
            h = set     # fixed field, change if needed
            k_l = k_r = val

        # --- Hamiltonian ---
        # polarized boundary conditions
        strength = 1e-1
        pol_bc = [[strength, 0], [-1**L * strength, L-1]]

        J_term = [[J, i, i+1] for i in range(L-1)]   
        h_term = [[h, i] for i in range(L)]  

        xy_term = [[k_r, i, i+2] for i in range(L-2)]
        yx_term = [[-k_l, i, i+2] for i in range(L-2)]
        xyz_term = [[k_l, i, i+1, i+2] for i in range(L-2)]
        zyx_term = [[-k_r, i, i+1, i+2] for i in range(L-2)]

        static = [
            ["xx", J_term],
            ["z",  h_term],
            ["x", pol_bc],
            ["xy", xy_term],
            ["yx", yx_term],
            ["xyz", xyz_term],
            ["zyx", zyx_term]
        ]

        H = hamiltonian(static, [], basis=basis, dtype=np.complex128, check_symm=False, check_herm=False)
        E, V = H.eigh()
        psi = V[:, 0]  # ground state
        print(psi.shape)
        # write it in a .txt file for easy access
        np.savetxt(os.path.join(os.getcwd(), folder_eigenstates, f"{opp}{set:.2f}_{scan_var}{val:.2f}_groundstate"), psi)
        np.savetxt(os.path.join(os.getcwd(), folder_eigenvalues, f"{opp}{set:.2f}_{scan_var}{val:.2f}_energies"), E)

# --- Magnetization ---
def magnetization_string(op, psi, basis):
    '''Magnetization across entire string'''
    op_string = []
    L = int(np.log2(len(psi)))

    for i in range(L):
        operator = hamiltonian([[op, [[1.0, i]]]], [], basis=basis, check_symm=False, check_herm=False)
        mx = np.real_if_close(psi.conj() @ (operator.dot(psi)))
        op_string.append(mx)
    return op_string

def magnetization_site(op, site, psi, basis):
    '''Magnetization on singular site'''
    op_mid = hamiltonian([[op, [[1, site]]]], [], basis=basis, check_symm=False, check_herm=False)
    expectation_value_mid = np.real_if_close(psi.conj() @ op_mid.dot(psi))

    return expectation_value_mid

# --- Correlation ---
def correlator(op1, op2, i, j, psi, basis):
    '''Connected correlator defined as g(op1, op2) = <op1, op2> - <op1><op2>'''
    operators = [ [op1, [[1, i]]],
                  [op2, [[1, j]]]
                 ]
    corr_op = hamiltonian(operators, [], basis=basis, check_symm=False, check_herm=False)
    corr = np.real_if_close(psi.conj() @ corr_op.dot(psi))
    connected_correlator = np.real(corr - (magnetization_site(op1, i, psi, basis) * magnetization_site(op2, j, psi, basis)))

    return connected_correlator

def correlator_string(op1, op2, central_site, psi, basis):
    L = int(np.log2(len(psi)))

    correlations = []
    for i in range(L):
        if i > central_site: correlations.append(correlator(op1, op2, central_site, i, psi, basis))
    return correlations

# --- Entanglement ---
def entanglement_site(psi, basis):
    L = int(np.log2(len(psi)))
    mid = list(range(L//2))
    ent_mid = basis.ent_entropy(psi, sub_sys_A=mid)["Sent_A"]
    return ent_mid

def entanglement_string(psi, basis):
    s_sites = []
    L = int(np.log2(len(psi)))

    for i in range(L):
        ent = basis.ent_entropy(psi, sub_sys_A=[i])
        s_sites.append(ent["Sent_A"])

    return s_sites

# ----- plotting ----- #
def plot_site(dmrg_path, op, values, site_measurements_ed, scan_var, opp, set):
    '''Plot the given results of a site for a range of parameters'''
    cmap = plt.cm.viridis
    cmap2 = plt.cm.inferno
    color = cmap(0.6)
    color2 = cmap2(0.6)
    L = int(len(values))
    site_measurements_dmrg = []
    loc = os.path.join(dmrg_path, f"*/OUT/*/{op}.txt")
    for i in sorted(glob.glob(loc)):
        r = np.loadtxt(i).T[1]
        site_measurements_dmrg.append(r[int(len(r)//2) - 2])

    plt.plot(values, site_measurements_ed, "-o", color=color, label="ED") # ED results
    plt.plot(values, site_measurements_dmrg, "-x", color=color2, label="DMRG") # DMRG results

    plt.title(f"Mid-chain magnetization of {op} vs {scan_var} at {opp}={set}")
    plt.xlabel(scan_var)
    plt.grid()
    plt.show()

def plot_string(dmrg_path, op, values, string_measurements_ed, scan_var, opp, set):
    '''Plot the results of the magnetization string'''
    plt.figure(figsize=(12,8))

    cmap = plt.cm.viridis  # choose any: plasma, inferno, turbo, etc.
    cmap2 = plt.cm.inferno
    colors = cmap(np.linspace(0, 1, len(values)))
    colors2 = cmap2(np.linspace(0, 1, len(values)))
    
    for ind, (v, c, c2) in enumerate(zip(values, colors, colors2)):
        if scan_var == "h":
            loc = os.path.join(dmrg_path, f"{v:.2f}_{set:.2f}_{set:.2f}/OUT/out_{v:.3f}_{set:.3f}_{set:.3f}/{op}.txt")
        else:
            loc = os.path.join(dmrg_path, f"{set:.2f}_{v:.2f}_{v:.2f}/OUT/out_{set:.3f}_{v:.3f}_{v:.3f}/{op}.txt")
        site_measurements_dmrg = np.loadtxt(loc).T[1]

        plt.plot(
            range(len(string_measurements_ed.T)),
            string_measurements_ed[ind],
            "-o",
            color=c,
            label=f"ED, {scan_var}={v:.2f}"
        )

        plt.scatter(
            np.arange(len(string_measurements_ed.T))[2:-2],
            site_measurements_dmrg, 
            marker="x",
            color="black",
            label=f"DMRG, {scan_var}={v:.2f}"
        )
    
    plt.title(f"Magnetization {op} (full string) for {opp}={set}")
    plt.xlabel("Site index")
    plt.ylabel(f"{op}")
    plt.legend()
    plt.grid()
    plt.show()

def plot_correlator(dmrg_path, op, values, measurements_ed, scan_var, opp, set, orientation="semilogy"):
    '''Plot the correlation'''
    plt.figure(figsize=(12, 8))

    cmap = plt.cm.viridis  # choose any: plasma, inferno, turbo, etc.
    cmap2 = plt.cm.inferno
    colors = cmap(np.linspace(0, 1, len(values)))
    colors2 = cmap2(np.linspace(0, 1, len(values)))



    for ind, (v, c, c2) in enumerate(zip(values, colors, colors2)):
        L = len(measurements_ed.T)
        x_ed = np.arange(1, L+1, 1)
        y_ed = np.abs(measurements_ed[ind])
        label = f"{scan_var}={v:.2f}"

        if scan_var == "h":
            loc = os.path.join(dmrg_path, f"{v:.2f}_{set:.2f}_{set:.2f}/OUT/out_{v:.3f}_{set:.3f}_{set:.3f}/{op}.txt")
        else:
            loc = os.path.join(dmrg_path, f"{set:.2f}_{v:.2f}_{v:.2f}/OUT/out_{set:.3f}_{v:.3f}_{v:.3f}/{op}.txt")


        r = np.loadtxt(loc, dtype=np.complex128)
        x_dmrg, y_dmrg = [], []

        # extract correlations we're interested in
        for i, j, corr in r:  
            if i == L//2 and j > i:
                x_dmrg.append(j-i-1)
                y_dmrg.append(corr)

        if orientation == "semilogy":
            plt.semilogy(x_ed, y_ed, "-o", color=c, label="ED, " + label)
            plt.semilogy(x_dmrg, y_dmrg, "-x", color=c2, label="DMRG, " + label)
        elif orientation == "loglog":
            plt.loglog(x_ed, y_ed, "-o", color=c, label="ED, " + label)
            plt.loglog(x_dmrg, y_dmrg, "-x", color=c2, label="DMRG, " + label)
        else:
            plt.plot(x_ed, y_ed, "-o", color=c, label="ED, " + label)
            plt.plot(x_dmrg, y_dmrg, "-x", color=c2, label="DMRG, " + label)
    
    plt.title(f"Magnetization {op} (full string) for {opp}={set}")
    plt.xlabel("Site index")
    plt.ylabel(f"{op}")
    plt.legend()
    plt.grid()
    plt.show()

def plot_entropy(dmrg_path, single_site, half_chain, scan_var, values, opp, set):
    fig, ax = plt.subplots(1, 2, figsize=(16, 10))
    fig.suptitle("Entropy")
    for i, v in enumerate(values):
        ax[0].plot(single_site[i], label=f"{scan_var}={v:.2f}")
    ax[0].set_title(f"Single-site entanglement entropy for {opp}={set}")


    ax[1].plot(values, half_chain, "-o")
    ax[1].set_title(f"Half-chain entanglement vs {scan_var}")
    ax[1].set_xlabel(f"{scan_var}")
    ax[1].set_ylabel("S_half")

    plt.legend()
    plt.grid()
    plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# import scipy.sparse as sp
# import scipy.sparse.linalg as lg

# I = sp.eye(2, format="csr", dtype=complex)
# X = sp.csr_matrix([[0, 1], [1, 0]], dtype=complex)
# Y = sp.csr_matrix([[0, -1j], [1j, 0]], dtype=complex)
# Z = sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)

# def kron_all(ops):
#     out = ops[0]
#     for op in ops[1:]:
#         out = sp.kron(out, op, format="csr")
#     return out

# def Ham(J_x, h_z, k_left, k_right, L):

#     H = sp.csr_matrix((2**L, 2**L), dtype=complex)

#     for i in range(L):
#         # XX
#         ops = []
#         for j in range(L):
#             if j == i:
#                 ops.append(X)
#             elif j == (i + 1) % L:
#                 ops.append(X)
#             else:
#                 ops.append(I)
#         H += J_x * kron_all(ops)

#         # Z
#         ops = []
#         for j in range(L):
#             if j == i:
#                 ops.append(Z)
#             else:
#                 ops.append(I)
#         H += h_z * kron_all(ops)
        
#         # XIY
#         ops = []
#         for j in range(L):
#             if j == i:
#                 ops.append(X)
#             elif j == (i + 2) % L:
#                 ops.append(Y)
#             else:
#                 ops.append(I)
#         H += k_right * kron_all(ops)

#         # YIX
#         ops = []
#         for j in range(L):
#             if j == i:
#                 ops.append(Y)
#             elif j == (i + 2) % L:
#                 ops.append(X)
#             else:
#                 ops.append(I)
#         H -= k_left * kron_all(ops)

#         # XYZ
#         ops = []
#         for j in range(L):
#             if j == i:
#                 ops.append(X)
#             elif j == (i + 1) % L:
#                 ops.append(Y)
#             elif j == (i + 2) % L:
#                 ops.append(Z)
#             else:
#                 ops.append(I)
#         H += k_left * kron_all(ops)

#         # ZYX
#         ops = []
#         for j in range(L):
#             if j == i:
#                 ops.append(Z)
#             elif j == (i + 1) % L:
#                 ops.append(Y)
#             elif j == (i + 2) % L:
#                 ops.append(X)
#             else:
#                 ops.append(I)
#         H -= k_right * kron_all(ops)

#     return H.tocsr()

# def operator(i, operator, psi0):
#     """Expectation value <psi | operator_i | psi>"""
    
#     L = int(np.log2(len(psi0)))
    
#     # Build operator list
#     op_list = []
#     for site in range(L):
#         if site == i:
#             op_list.append(operator)
#         else:
#             op_list.append(I)
    
#     O = kron_all(op_list)

#     return (psi0.conj().T @ (O @ psi0)).real # Return expectation value

# def correlator(i, j, operator1, operator2, psi0):
#     ''' 
#     i, j are the two sites you want to find the correlators for,
#     operator1 and operator2 are the two operators you are searching,
#     psi0 is your wavefunction
#     '''

#     L = np.log2(len(psi0))

#     op = []
#     for k in range(L):
#         if k == i:
#             op.append(operator1)
#         elif k == j:
#             op.append(operator2)
#         else:
#             op.append(I)
#     O = kron_all(op)
#     correlator = (psi0.conj().T @ (O @ psi0)).real
#     connected_corr = correlator - operator(i, operator1, psi0) * operator(j, operator2, psi0)

#     return connected_corr