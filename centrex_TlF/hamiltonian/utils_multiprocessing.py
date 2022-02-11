import numpy as np
from centrex_TlF.hamiltonian.hamiltonian_terms_uncoupled import (
    Hff_X,
    HSx,
    HSy,
    HSz,
    HZx_X,
    HZy_X,
    HZz_X,
    Hff_X_alt,
)
from centrex_TlF.hamiltonian.hamiltonian_B_terms_coupled import (
    Hrot_B,
    H_mhf_Tl,
    H_mhf_F,
    H_LD,
    H_cp1_Tl,
    H_c_Tl,
    HZz_B,
)


def multi_transformation_matrix(i, state1, basis2):
    transform = np.zeros((len(basis2), len(basis2)))

    for j, state2 in enumerate(basis2):
        transform[i, j] = state1 @ state2

    return transform


def multi_HMatElems(H, i, a, QN):
    nstates = len(QN)
    result = np.zeros((nstates, nstates), dtype=complex)
    H = eval(H)
    for j in range(i, nstates):
        b = QN[j]
        val = (1 * a) @ H(b)
        result[i, j] = val
        if i != j:
            result[j, i] = np.conjugate(val)
    return result
