import centrex_TlF
import centrex_TlF.states as states
import numpy as np
import scipy
from centrex_TlF.hamiltonian.basis_transform import generate_transform_matrix
from centrex_TlF.hamiltonian.generate_hamiltonian import (
    generate_coupled_hamiltonian_B,
    generate_uncoupled_hamiltonian_X,
)
from centrex_TlF.hamiltonian.utils import (
    generate_coupled_hamiltonian_B_function,
    generate_uncoupled_hamiltonian_X_function,
    matrix_to_states,
    reduced_basis_hamiltonian,
    reorder_evecs,
)

__all__ = [
    "generate_reduced_X_hamiltonian",
    "generate_reduced_B_hamiltonian",
    "compose_reduced_hamiltonian",
    "generate_total_reduced_hamiltonian",
    "generate_diagonalized_hamiltonian",
]


def generate_diagonalized_hamiltonian(
    hamiltonian, keep_order=True, return_V_ref=False, rtol=None
):
    D, V = np.linalg.eigh(hamiltonian)
    if keep_order:
        V_ref = np.eye(V.shape[0])
        D, V = reorder_evecs(V, D, V_ref)

    hamiltonian_diagonalized = V.conj().T @ hamiltonian @ V
    if rtol:
        hamiltonian_diagonalized[
            np.abs(hamiltonian_diagonalized)
            < np.abs(hamiltonian_diagonalized).max() * rtol
        ] = 0
    if not return_V_ref or not keep_order:
        return hamiltonian_diagonalized, V
    else:
        return hamiltonian_diagonalized, V, V_ref


def generate_reduced_X_hamiltonian(
    ground_states_approx, E=np.array([0, 0, 0]), B=np.array([0, 0, 0.001]), rtol=None
):

    QN = states.generate_uncoupled_states_ground(
        np.unique([gs.J for gs in ground_states_approx])
    )
    H_X_uc = generate_uncoupled_hamiltonian_X(QN)
    H_X_uc = generate_uncoupled_hamiltonian_X_function(H_X_uc)
    S_transform = generate_transform_matrix(QN, ground_states_approx)

    H_X = S_transform.conj().T @ H_X_uc(E, B) @ S_transform
    if rtol:
        H_X[np.abs(H_X) < np.abs(H_X).max() * rtol] = 0

    # diagonalize the Hamiltonian
    H_X_diag, V, V_ref = generate_diagonalized_hamiltonian(
        H_X, keep_order=True, return_V_ref=True, rtol=rtol
    )

    # new set of quantum numbers:
    QN_diag = matrix_to_states(V, ground_states_approx)

    ground_states = states.find_exact_states(
        [1 * gs for gs in ground_states_approx], H_X_diag, QN_diag, V_ref=V_ref
    )

    # ground_states = [gs.remove_small_components() for gs in ground_states]

    H_X_red = reduced_basis_hamiltonian(QN_diag, H_X_diag, ground_states)

    return ground_states, H_X_red


def generate_reduced_B_hamiltonian(
    excited_states_approx,
    E=np.array([0, 0, 0]),
    B=np.array([0, 0, 0.001]),
    rtol=None,
    Jmin=1,
    Jmax=3,
):
    # need to generate other states because excited states are mixed
    Ps = [-1, 1]
    I_F = 1 / 2
    I_Tl = 1 / 2
    QN_B = [
        centrex_TlF.CoupledBasisState(
            F, mF, F1, J, I_F, I_Tl, P=P, Omega=1, electronic_state="B"
        )
        for J in np.arange(Jmin, Jmax + 1)
        for F1 in np.arange(np.abs(J - I_F), J + I_F + 1)
        for F in np.arange(np.abs(F1 - I_Tl), F1 + I_Tl + 1)
        for mF in np.arange(-F, F + 1)
        for P in Ps
    ]

    for qn in excited_states_approx:
        assert qn.isCoupled, "supply list of CoupledBasisStates"
    H_B = generate_coupled_hamiltonian_B(QN_B)
    H_B = generate_coupled_hamiltonian_B_function(H_B)

    H_B_diag, V, V_ref_B = generate_diagonalized_hamiltonian(
        H_B(E, B), keep_order=True, return_V_ref=True, rtol=rtol
    )

    # new set of quantum numbers:
    QN_B_diag = matrix_to_states(V, QN_B)

    excited_states = centrex_TlF.states.find_exact_states(
        [1 * e for e in excited_states_approx], H_B_diag, QN_B_diag, V_ref=V_ref_B
    )

    H_B_red = reduced_basis_hamiltonian(QN_B_diag, H_B_diag, excited_states)
    return excited_states, H_B_red


def compose_reduced_hamiltonian(H_X_red, H_B_red, element_limit=0.1):
    H_X_red[np.abs(H_X_red) < element_limit] = 0
    H_B_red[np.abs(H_B_red) < element_limit] = 0

    H_int = scipy.linalg.block_diag(H_X_red, H_B_red)
    V_ref_int = np.eye(H_int.shape[0])

    return H_int, V_ref_int


def generate_total_reduced_hamiltonian(
    ground_states_approx, excited_states_approx, Jmin=None, Jmax=None, rtol=None
):
    ground_states, H_X_red = generate_reduced_X_hamiltonian(
        ground_states_approx, rtol=rtol
    )

    # Js to include for rotational mixing in B state
    Jexc = np.unique([s.J for s in excited_states_approx])
    if not Jmin:
        Jmin = np.min(Jexc)
        Jmin = 1 if Jmin - 1 < 1 else Jmin - 1
    if not Jmax:
        Jmax = np.max(Jexc) + 1

    excited_states, H_B_red = generate_reduced_B_hamiltonian(
        excited_states_approx, Jmin=Jmin, Jmax=Jmax, rtol=rtol
    )

    H_int, V_ref_int = compose_reduced_hamiltonian(H_X_red, H_B_red)

    QN = list(np.append(ground_states, excited_states))
    return ground_states, excited_states, QN, H_int, V_ref_int
