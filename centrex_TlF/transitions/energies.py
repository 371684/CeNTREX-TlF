import json
import scipy
import pickle
import numpy as np
from tqdm import tqdm
from pathlib import Path
from sympy import Rational
from centrex_TlF.states.states import CoupledBasisState, State
from centrex_TlF.states.utils import (
    BasisStates_from_State,
    find_state_idx_from_state,
    find_states_idxs_from_states,
    find_closest_vector_idx,
)
from centrex_TlF.hamiltonian import (
    generate_uncoupled_hamiltonian_X,
    generate_uncoupled_hamiltonian_X_function,
    generate_coupled_hamiltonian_B,
    generate_coupled_hamiltonian_B_function,
    generate_diagonalized_hamiltonian,
    matrix_to_states,
)
from functools import lru_cache

__all__ = [
    "calculate_energies",
    "calculate_transition_frequency",
    "generate_transition_frequency",
    "generate_transition_frequencies",
    "find_transition",
    "identify_transition",
]


def calculate_energies(
    ground_states,
    excited_states,
    E=np.array([0, 0, 0]),
    B=np.array([0, 0, 0]),
    offset_X=1
    * CoupledBasisState(
        J=0,
        F1=1 / 2,
        F=1,
        mF=0,
        I1=1 / 2,
        I2=1 / 2,
        electronic_state="X",
        P=+1,
        Omega=0,
    ),
    offset_B=1
    * CoupledBasisState(
        J=1,
        F1=3 / 2,
        F=2,
        mF=0,
        I1=1 / 2,
        I2=1 / 2,
        electronic_state="B",
        P=-1,
        Omega=1,
    ),
    offset=2 * np.pi * 1103407960000000.0,
    nprocs=1,
):

    # transform coupled basis to uncoupled
    QN_X_uc = [s.transform_to_uncoupled() for s in ground_states]
    # grab unique BasisStates from the State objects
    QN_X_uc = BasisStates_from_State(QN_X_uc)
    # calculate and diagonalize the X state Hamiltonian
    H_X = generate_uncoupled_hamiltonian_X(QN_X_uc, nprocs=nprocs)
    H_X = generate_uncoupled_hamiltonian_X_function(H_X)(E, B)
    H_X, V_X = generate_diagonalized_hamiltonian(H_X, keep_order=False)
    # generate the corresponding ground states in coupled basis
    QN_X = matrix_to_states(V_X, QN_X_uc)
    QN_X = [s.transform_to_coupled() for s in QN_X]

    # calculate and diagonalize the B state Hamiltonian
    H_B = generate_coupled_hamiltonian_B(excited_states, nprocs=nprocs)
    H_B = generate_coupled_hamiltonian_B_function(H_B)(E, B)
    H_B, V_B = generate_diagonalized_hamiltonian(H_B, keep_order=False)
    QN_B = matrix_to_states(V_B, excited_states)

    # add offset
    idX = find_state_idx_from_state(H_X, offset_X, QN_X)
    idB = find_state_idx_from_state(H_B, offset_B, QN_B)
    Eshift = (H_X[idX, idX] - H_B[idB, idB]) + offset
    H_B += Eshift * np.eye(len(excited_states))

    H_tot = scipy.linalg.block_diag(H_X, H_B)
    QN = list(QN_X) + list(QN_B)

    return QN, H_tot


def calculate_transition_frequency(state1, state2, H, QN):
    """
    Function that outputs the frequency of the transition between state1 and state2
    which are assumed to be eigenstates of the Hamiltonian H whose basis is defined by
    QN. The function does not check whether or not the transition is allowed. If state1
    state2 are not eigenstates of H, the function will find the eigenstates with
    maximum overlap for state1 and state2 and calculate frequency between those.

    inputs:
    state1      : State object or state vector (np.ndarray) representing the first state of interest
    state2      : State object or state vector (np.ndarray) representing the second state of interest
    H           : Hamiltonian that is used to calculate the energies of states 1 and 2
                  (assumed to be in angular frequency units - 2pi*Hz)
    QN          : List of State objects that defines the basis for the Hamiltonian

    returns:
    freq        : Transition frequency between states 1 and 2 in Hz
    """
    # If provided with State objects, convert them to state vectors (np.ndarray)
    if not isinstance(state1, np.ndarray):
        state1 = state1.state_vector(QN)

    if not isinstance(state2, np.ndarray):
        state2 = state2.state_vector(QN)

    # Diagonalize hamiltonian
    D, V = np.linalg.eigh(H)

    # Find the indices that correspond to each state
    i1 = find_closest_vector_idx(state1, V)
    i2 = find_closest_vector_idx(state2, V)

    # Find energies
    E1 = D[i1]
    E2 = D[i2]

    # Return transitions frequency in Hz
    return np.abs(E1 - E2) / (2 * np.pi)


def _check_precached(state, config=None):
    if not config:
        path = Path(__file__).parent.parent / "pre_calculated"
        js = path / "precalculated.json"
        with open(js) as json_file:
            config = json.load(json_file)

    Js = np.unique([s.J for _, s in state])
    es = state.find_largest_component().electronic_state
    for J in Js:
        # if B state check if J + 1 is included in pre-cached values to account
        # for mixing
        if es == "B":
            J += 1
        assert (
            J in config["transitions"][es]
        ), f"Hamiltonian term for |{es}, {J}> not pre-cached"


@lru_cache(maxsize=2 ** 16)
def generate_transition_frequency(state1, state2):
    """Get the field free transition frequency between state1 and state2 in 2π⋅Hz.
    Grabs pre-cached version, otherwise throws exception.

    Args:
        state1 (State): TlF State
        state2 (State): TlF State

    Returns:
        frequency: transition frequency in 2π⋅Hz
    """
    path = Path(__file__).parent.parent / "pre_calculated"
    js = path / "precalculated.json"

    # check if state1 and state2 are included in the pre-cached transitions
    _check_precached(state1)
    _check_precached(state2)

    with open(path / "transitions.pickle", "rb") as f:
        _ = pickle.load(f)
        QN, H = _["QN"], _["H"]

    return calculate_transition_frequency(state1, state2, H, QN)


def generate_transition_frequencies(states1, states2):
    """Get the field free transition frequencies between states1 and states2 in 
    2π⋅Hz.
    Grabs pre-cached version, otherwise throws exception.

    Args:
        state1 (State): TlF State
        state2 (State): TlF State

    Returns:
        frequency: transition frequency in 2π⋅Hz
    """
    path = Path(__file__).parent.parent / "pre_calculated"
    js = path / "precalculated.json"

    # check if states1 and states2 are included in the pre-cached transitions
    for state1, state2 in zip(states1, states2):
        _check_precached(state1)
        _check_precached(state2)

    with open(path / "transitions.pickle", "rb") as f:
        _ = pickle.load(f)
        QN, H = _["QN"], _["H"]

    indices1 = find_states_idxs_from_states(H, states1, QN)
    indices2 = find_states_idxs_from_states(H, states2, QN)
    E = np.diag(H)
    E1 = E[indices1]
    E2 = E[indices2]
    return E2 - E1


def find_transition(transition, F1, F, eg="X", ee="B", return_states=False):
    transitions = {"R": +1, "P": -1, "Q": 0, "S": +2, "O": -2, "T": +3}

    Je, Jg = transition
    Jg = int(Jg)
    Je = Jg + transitions[Je]

    ground_state = 1 * CoupledBasisState(
        J=Jg,
        F1=Jg + 1 / 2,
        F=Jg,
        mF=0,
        I1=1 / 2,
        I2=1 / 2,
        P=(-1) ** Jg,
        Omega=0,
        electronic_state=eg,
    )
    excited_state = 1 * CoupledBasisState(
        J=Je,
        F1=F1,
        F=F,
        mF=0,
        I1=1 / 2,
        I2=1 / 2,
        P=(-1) ** (Jg + 1),
        Omega=1,
        electronic_state=ee,
    )
    frequency = generate_transition_frequency(ground_state, excited_state)
    if return_states:
        return ground_state, excited_state, frequency
    else:
        return frequency


def identify_transition(state1, state2):
    assert state1 != state2, "no transition between same states"
    transitions = {0: "Q", +1: "R", -1: "P", +2: "S", -2: "O", +3: "T"}
    state1 = state1.find_largest_component()
    state2 = state2.find_largest_component()
    assert state1.electronic_state != "B", "state1 required to be a ground state"
    assert state2.electronic_state != "X", "state2 required to be an excited state"
    assert state1.isCoupled, "supply state1 in coupled basis"
    assert state2.isCoupled, "supply state2 in coupled basis"
    Jg, Je = state1.J, state2.J
    transition = transitions[Je - Jg]
    string = f"{transition}({Jg}) F1'={Rational(state2.F1)}, F'={Rational(state2.F)}"
    return string


# class Transition:
#     transitions_nom = {0: 'Q', +1: 'R', -1: 'P', +2: 'S', -2: 'O', +3: 'T'}
#     transitions_ΔJ = {'R': +1, 'P': -1, 'Q': 0, 'S': +2, 'O': -2, 'T': +3}

#     def __init__(self, transition, F1, F, eg = 'X', ee = 'B'):
#         Je,Jg = transition
#         self.transition = Je
#         self.Jg = int(Jg)
#         self.Je = self.Jg+self.transitions_ΔJ[Je]
#         self.F1 = F1
#         self.F = F


#         self.ground_state = 1*CoupledBasisState(
#                                 J=self.Jg, F1 = self.Jg+1/2, F = self.Jg, mF = 0,
#                                 I1 = 1/2, I2 = 1/2, P = (-1)**self.Jg, Omega = 0,
#                                 electronic_state = eg
#                                 )
#         self.excited_state = 1*CoupledBasisState(
#                                 J=self.Je, F1 = F1, F = F, mF = 0, I1 = 1/2,
#                                 I2 = 1/2, P = (-1)**(self.Jg+1), Omega = 1,
#                                 electronic_state = ee
#                                 )
#         self.frequency = generate_transition_frequency(
#                             self.ground_state, self.excited_state
#                         )

#     def __repr__(self):
#         string = \
#             f"{self.transition}({self.Jg}) F1'={Rational(self.F1)}, F'={Rational(self.F)}"
#         # string += f' -> {self.frequency/(2*np.pi*1e9):.2f} GHz'
#         return f"Transition({string})"
