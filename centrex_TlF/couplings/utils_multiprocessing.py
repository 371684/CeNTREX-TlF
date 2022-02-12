import numpy as np
from centrex_TlF.couplings.matrix_elements import calculate_ED_ME_mixed_state


def multi_coupling_matrix(QN, ground_state, excited_states, pol_vec, reduced):
    H = np.zeros((len(QN), len(QN)), dtype=complex)
    i = QN.index(ground_state)
    for excited_state in excited_states:
        j = QN.index(excited_state)

        # calculate matrix element and add it to the Hamiltonian
        H[i, j] = calculate_ED_ME_mixed_state(
            ground_state, excited_state, pol_vec=pol_vec, reduced=reduced
        )
    return H
