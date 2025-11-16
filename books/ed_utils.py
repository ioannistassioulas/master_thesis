from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
import numpy as np
import matplotlib.pyplot as plt

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
    corr = hamiltonian([[op1, [[1, i]]], [[op2, [[1, j]]]]], [], basis=basis, check_symm=False, check_herm=False)
    connected_correlator = corr - (magnetization_site(op1, i, psi, basis) * magnetization_site(op2, j, psi, basis))

    return connected_correlator

# ----- plotting ----- #
def plot_site(op, values, site_measurements, scan_var, opp, set):
    '''Plot the given results of a site for a range of parameters'''
    plt.figure(figsize=(6,4))
    plt.plot(values, site_measurements, "-o")
    plt.title(f"Mid-chain magnetization of {op} vs {scan_var} at {opp}={set}")
    plt.xlabel(scan_var)
    plt.grid()
    plt.show()

def plot_string(op, values, string_measurements, scan_var, opp, set):
    '''Plot the results of the magnetization string'''
    plt.figure(figsize=(6,4))

    for ind, v in enumerate(values):
        plt.plot(range(len(string_measurements.T)), string_measurements[ind], "-o", label=f"{scan_var}={v:.2f}")
    
    plt.title(f"Magnetization {op} (full string) for {opp}={set}")
    plt.xlabel("Site index")
    plt.ylabel(f"{op}")
    plt.legend()
    plt.grid()
    plt.show()

def plot_correlator(op1, op2, values, measurements, scan_var, opp, set, orientation="semilogy"):
    '''Plot the correlation'''
    plt.figure(figsize=(6, 4))


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