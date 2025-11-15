import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as lg

I = sp.eye(2, format="csr", dtype=complex)
X = sp.csr_matrix([[0, 1], [1, 0]], dtype=complex)
Y = sp.csr_matrix([[0, -1j], [1j, 0]], dtype=complex)
Z = sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)

def kron_all(ops):
    out = ops[0]
    for op in ops[1:]:
        out = sp.kron(out, op, format="csr")
    return out

def Ham(J_x, h_z, k_left, k_right, L):

    H = sp.csr_matrix((2**L, 2**L), dtype=complex)

    for i in range(L):
        # XX
        ops = []
        for j in range(L):
            if j == i:
                ops.append(X)
            elif j == (i + 1) % L:
                ops.append(X)
            else:
                ops.append(I)
        H += J_x * kron_all(ops)

        # Z
        ops = []
        for j in range(L):
            if j == i:
                ops.append(Z)
            else:
                ops.append(I)
        H += h_z * kron_all(ops)
        
        # XIY
        ops = []
        for j in range(L):
            if j == i:
                ops.append(X)
            elif j == (i + 2) % L:
                ops.append(Y)
            else:
                ops.append(I)
        H += k_right * kron_all(ops)

        # YIX
        ops = []
        for j in range(L):
            if j == i:
                ops.append(Y)
            elif j == (i + 2) % L:
                ops.append(X)
            else:
                ops.append(I)
        H -= k_left * kron_all(ops)

        # XYZ
        ops = []
        for j in range(L):
            if j == i:
                ops.append(X)
            elif j == (i + 1) % L:
                ops.append(Y)
            elif j == (i + 2) % L:
                ops.append(Z)
            else:
                ops.append(I)
        H += k_left * kron_all(ops)

        # ZYX
        ops = []
        for j in range(L):
            if j == i:
                ops.append(Z)
            elif j == (i + 1) % L:
                ops.append(Y)
            elif j == (i + 2) % L:
                ops.append(X)
            else:
                ops.append(I)
        H -= k_right * kron_all(ops)

    return H.tocsr()

def operator(i, operator, psi0):
    """Expectation value <psi | operator_i | psi>"""
    
    L = int(np.log2(len(psi0)))
    
    # Build operator list
    op_list = []
    for site in range(L):
        if site == i:
            op_list.append(operator)
        else:
            op_list.append(I)
    
    O = kron_all(op_list)

    return (psi0.conj().T @ (O @ psi0)).real # Return expectation value

def correlator(i, j, operator1, operator2, psi0):
    ''' 
    i, j are the two sites you want to find the correlators for,
    operator1 and operator2 are the two operators you are searching,
    psi0 is your wavefunction
    '''

    L = np.log2(len(psi0))

    op = []
    for k in range(L):
        if k == i:
            op.append(operator1)
        elif k == j:
            op.append(operator2)
        else:
            op.append(I)
    O = kron_all(op)
    correlator = (psi0.conj().T @ (O @ psi0)).real
    connected_corr = correlator - operator(i, operator1, psi0) * operator(j, operator2, psi0)

    return connected_corr