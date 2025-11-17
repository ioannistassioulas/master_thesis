from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

# --- Finding GS ---
def exact_diagonalization_line(values, scan_var, opp, set, basis, L, J=1):
    '''Perform exact diagonalization across a line on the parameter space of your phase diagram'''

    folder = f"L{L}_{opp}{set}_{scan_var}{np.min(values)}-{np.max(values)}"
    os.mkdir(folder)
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

        # write it in a .txt file for easy access
        file_name = os.path.join(os.getcwd(), folder, f"{opp}{set:.2f}_{scan_var}{val:.2f}_groundstate")
        np.savetxt(file_name, psi)


# --- Taking measurements ---
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
def plot_site(op, values, site_measurements, scan_var, opp, set):
    '''Plot the given results of a site for a range of parameters'''
    cmap = plt.cm.viridis
    color = cmap(0.6)

    plt.plot(values, site_measurements, "-o", color=color)
    plt.plot(values, site_measurements, "-o")
    plt.title(f"Mid-chain magnetization of {op} vs {scan_var} at {opp}={set}")
    plt.xlabel(scan_var)
    plt.grid()
    plt.show()

def plot_string(op, values, string_measurements, scan_var, opp, set):
    '''Plot the results of the magnetization string'''
    plt.figure(figsize=(6,4))

    cmap = plt.cm.viridis  # choose any: plasma, inferno, turbo, etc.
    colors = cmap(np.linspace(0, 1, len(values)))

    for ind, (v, c) in enumerate(zip(values, colors)):
        plt.plot(
            range(len(string_measurements.T)),
            string_measurements[ind],
            "-o",
            color=c,
            label=f"{scan_var}={v:.2f}"
        )
    
    plt.title(f"Magnetization {op} (full string) for {opp}={set}")
    plt.xlabel("Site index")
    plt.ylabel(f"{op}")
    plt.legend()
    plt.grid()
    plt.show()

def plot_correlator(op, values, measurements, scan_var, opp, set, orientation="semilogy"):
    '''Plot the correlation'''
    plt.figure(figsize=(6, 4))

    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0, 1, len(values)))

    for ind, (v, c) in enumerate(zip(values, colors)):
        x = range(len(measurements.T))
        y = np.abs(measurements[ind])
        label = f"{scan_var}={v:.2f}"

        if orientation == "semilogy":
            plt.semilogy(x, y, "-o", color=c, label=label)
        elif orientation == "loglog":
            plt.loglog(x, y, "-o", color=c, label=label)
        else:
            plt.plot(x, y, "-o", color=c, label=label)
    
    plt.title(f"Magnetization {op} (full string) for {opp}={set}")
    plt.xlabel("Site index")
    plt.ylabel(f"{op}")
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