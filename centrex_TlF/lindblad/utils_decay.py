import copy
from dataclasses import dataclass

import numpy as np
import sympy as smp
from centrex_TlF.states import CoupledBasisState

__all__ = ["DecayChannel"]


@dataclass
class DecayChannel:
    ground: CoupledBasisState
    excited: np.ndarray
    branching: float
    description: str = ""


def add_levels_symbolic_hamiltonian(hamiltonian, decay_channels, QN, excited_states):
    arr = hamiltonian.copy()
    indices = [i + len(QN) - len(excited_states) for i in range(len(decay_channels))]
    for idx in indices:
        arr = add_level_symbolic_hamiltonian(arr, idx)
    return indices, arr


def add_level_symbolic_hamiltonian(hamiltonian, idx):
    arr = hamiltonian.copy()
    arr = arr.row_insert(idx, smp.zeros(1, arr.shape[1]))
    arr = arr.col_insert(idx, smp.zeros(arr.shape[0], 1))
    return arr


def add_states_QN(decay_channels, QN, indices):
    states = copy.copy(QN)
    for idx, decay_channel in zip(indices, decay_channels):
        states.insert(idx, decay_channel.ground)
    return states


def add_levels_C_array(C_array, indices):
    arr = C_array.copy()
    # inserting rows and columns of zeros to account for the new decay levels
    for idx in indices:
        arr = np.insert(arr, idx, np.zeros(arr.shape[2]), 1)
        arr = np.insert(arr, idx, np.zeros(arr.shape[1]), 2)
    return arr


def add_decays_C_arrays(decay_channels, indices, QN, C_array, Γ):
    # converting the C arrays to branching ratio arrays and adding the new
    # levels
    BR = add_levels_C_array(C_array, indices)
    BR = BR ** 2 / Γ

    # getting the excited state indices
    indices_excited = [
        decay_channel.excited.get_indices(QN) for decay_channel in decay_channels
    ]
    # getting the total added branching ratios for each excited state
    BR_added = {}
    for ides, decay_channel in zip(indices_excited, decay_channels):
        for ide in ides:
            if BR_added.get(ide) is None:
                BR_added[ide] = decay_channel.branching
            else:
                BR_added[ide] += decay_channel.branching
    # renormalizing the old branching ratios to ensure the sum is 1 when adding
    # the new branching ratios
    for ide, BR_add in BR_added.items():
        BR[:, :, ide] *= 1 - BR_add

    # adding the new branching ratios
    for idg, ides, decay_channel in zip(indices, indices_excited, decay_channels):
        for ide in ides:
            BR_new = np.zeros([1, *BR[0, :, :].shape], dtype=complex)
            BR_new[:, idg, ide] = decay_channel.branching
            BR = np.append(BR, BR_new, axis=0)
    # converting the branching ratios to C arrays
    return np.sqrt(BR * Γ)
