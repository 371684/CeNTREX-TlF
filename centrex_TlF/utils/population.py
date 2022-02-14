import numpy as np
import scipy.constants as cst
from centrex_TlF.couplings.utils_compact import delete_row_column
from centrex_TlF.states.utils import QuantumSelector

__all__ = [
    "thermal_population",
    "J_levels",
    "J_slice",
    "generate_thermal_J",
    "generate_population_states",
    "generate_thermal_population_states",
]


def thermal_population(J, T, B=6.66733e9, n=100):
    """calculate the thermal population of a given J sublevel

    Args:
        J (int): rotational level
        T (float): temperature [Kelvin]
        B (float, optional): rotational constant. Defaults to 6.66733e9.
        n (int, optional): number of rotational levels to normalize with.
                            Defaults to 100.

    Returns:
        float: relative population in a rotational sublevel
    """
    c = 2 * np.pi * cst.hbar * B / (cst.k * T)

    def a(J):
        return -c * J * (J + 1)

    Z = np.sum([J_levels(i) * np.exp(a(i)) for i in range(n)])
    return J_levels(J) * np.exp(a(J)) / Z


def J_levels(J):
    """calculate the number of hyperfine sublevels per J rotational level

    Args:
        J (int): rotational level

    Returns:
        int: number of levels
    """
    return 4 * (2 * J + 1)


def J_slice(J):
    """generate a slice object for a rotational sublevel

    Args:
        J (int): rotational level

    Returns:
        numpy slice: numpy slice object
    """
    if J == 0:
        return np.s_[0 : J_levels(0)]
    else:
        levels = J_levels(np.arange(J + 1))
        return np.s_[np.sum(levels[:-1]) : np.sum(levels)]


def generate_thermal_J(Js, n_excited, T, normalized=True, slice_compact=None):
    """generate a thermal distribution over the rotational states

    Args:
        Js (list,array): included J levels
        T (int, float): Temperature in Kelvin

    Returns:
        np.ndarray: density matrix
    """
    # calculate the total number of levels in the system
    levels = np.sum([J_levels(J) for J in Js])

    # initialize empty density matrix
    ρ = np.zeros([levels + n_excited, levels + n_excited], dtype=complex)

    index = 0
    for J in Js:
        p = thermal_population(J, T)
        levels = J_levels(J)
        sl = np.s_[index : index + levels]
        np.fill_diagonal(ρ[sl, sl], p / levels)
        index += levels

    if normalized:
        # normalize the density matrix trace to 1
        ρ /= np.sum(np.diag(ρ))

    if slice_compact:
        ρ_compact = delete_row_column(ρ.copy(), slice_compact)
        range_compact = range(slice_compact.start, slice_compact.stop - 1)
        for idx in range_compact:
            ρ_compact[slice_compact.start, slice_compact.start] += ρ[idx, idx]
        ρ = ρ_compact

    return ρ


def generate_thermal_population_states(states_to_fill, states, T):
    """Generate a thermal distrubtion over the states specified in
    states_to_fill, a QuantumSelector or list of Quantumselectors

    Args:
        states_to_fill (QuantumSelector): Quantumselector specifying states to
        fill
        states (list, np.ndarray): all states used in simulation
        T (float): temperature in Kelvin

    Returns:
        np.ndarray: density matrix with trace normalized to 1
    """
    # branch for single QuantumSelector use
    if isinstance(states_to_fill, QuantumSelector):
        # get all involved Js
        Js = states_to_fill.J
        # check if J was a list
        if not isinstance(Js, (np.ndarray, list, tuple)):
            Js = [Js]
        # get indices of states to fill
        indices_to_fill = states_to_fill.get_indices(states)
    # branch for multiple QuantumSelectors use
    elif isinstance(states_to_fill, (list, np.ndarray, tuple)):
        # get all involved Js
        Js = []
        for stf in states_to_fill:
            J = stf.J
            # check if J was a list
            if not isinstance(J, (np.ndarray, list, tuple)):
                J = [J]
            Js.extend(J)
        # get indices of states to fill
        indices_to_fill = []
        for stf in states_to_fill:
            indices_to_fill.extend(stf.get_indices(states))

    # remove duplicates from Js and indices_to_fill
    Js = np.unique(Js)
    indices_to_fill = np.unique(indices_to_fill)

    # thermal population per hyperfine level for each involved J
    thermal_populations = dict(
        [(Ji, thermal_population(Ji, T) / J_levels(Ji)) for Ji in Js]
    )
    # generate an empty density matrix
    ρ = np.zeros([len(states), len(states)], dtype=complex)
    # fill the density matrix
    for idρ in indices_to_fill:
        state = states[idρ].find_largest_component()
        thermal_pop = thermal_populations[state.J]
        ρ[idρ, idρ] = thermal_pop
    # normalize the trace to 1 and return the density matrix
    return ρ / np.trace(ρ)


def generate_population_states(states, levels):
    """generate a uniform population distribution with population in the
    specified states

    Args:
        states (list, np.ndarray): indices to put population into
        levels (int): total number of levels

    Returns:
        np.ndarray: density matrix
    """
    ρ = np.zeros([levels, levels], dtype=complex)
    for state in states:
        ρ[state, state] = 1
    return ρ / np.trace(ρ)
