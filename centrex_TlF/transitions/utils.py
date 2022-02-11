import numpy as np
from centrex_TlF.states.utils import QuantumSelector
from centrex_TlF.states.generate_states import (
    generate_coupled_states_ground_X,
    generate_coupled_states_excited_B,
)

__all__ = [
    "check_transition_coupled_allowed",
    "assert_transition_coupled_allowed",
    "construct_ground_states_allowed",
]


def construct_ground_states_allowed(Jg, Je, F1, F, mF=None, P=None, ΔmF=None):
    """Construct the coupled ground states given a ground state J and excited
    state quantum numbers J', F1', F' and optionally mF'

    Args:
        Jg (int): rotational ground state
        Je (int): rotational excited state
        F1 (float): excited F1 state
        F (int): excited F state
        mF (int, optional): excited mF state. Defaults to None, which includes 
                            all excited state mF.
        P (int, optional): excited state parity. Defaults to None, which picks
                            the allowed excited state parity given the ground
                            state.
        ΔmF (int, optional): ΔmF of transition. Defaults to None, which checks
                            for all ΔmF.

    Returns:
        np.ndarray: np.ndarray of allowed ground states as CoupledBasisStates
    """
    if P is None:
        P = (-1) ** (Jg + 1)
    if ΔmF is None:
        ΔmF = [-1, 0, 1]
    if not isinstance(ΔmF, (list, tuple, np.ndarray)):
        ΔmF = [ΔmF]
    # generate states
    ground = QuantumSelector(J=Jg, electronic="X", P=(-1) ** Jg)
    excited = QuantumSelector(J=Je, F1=F1, F=F, mF=mF, electronic="B", P=P)
    ground = generate_coupled_states_ground_X(ground)
    excited = generate_coupled_states_excited_B(excited)

    # brute force loop over states to find allowed ground states
    ground_allowed = []
    for idg, g in enumerate(ground):
        for e in excited:
            for Δ in ΔmF:
                if check_transition_coupled_allowed(
                    g, e, ΔmF_allowed=Δ, return_err=False
                ):
                    ground_allowed.append(idg)
    ground_allowed = np.unique(ground_allowed)

    if len(ground_allowed) == 0:
        return []
    else:
        return ground[ground_allowed]


def check_transition_coupled_allowed(state1, state2, ΔmF_allowed, return_err=True):
    """Check whether the transition is allowed based on the quantum numbers

    Args:
        state1 (CoupledBasisState): ground CoupledBasisState
        state2 (CoupledBasisState): excited CoupledBasisState
        return_err (boolean): boolean flag for returning the error message

    Returns:
        tuple: (allowed boolean, error message)
    """
    ΔF = int(state2.F - state1.F)
    ΔmF = int(state2.mF - state1.mF)
    ΔP = int(state2.P - state1.P)

    flag_ΔP = np.abs(ΔP) != 2
    flag_ΔF = np.abs(ΔF) > 1
    flag_ΔmF = ΔmF != ΔmF_allowed
    flag_ΔFΔmF = (not flag_ΔmF) & ((ΔF == 0) & (ΔmF == 0) & (state1.mF == 0))

    errors = ""
    if flag_ΔP:
        errors += f"parity invalid"
    if flag_ΔF:
        if len(errors) != 0:
            errors += ", "
        errors += f"ΔF invalid -> ΔF = {ΔF}"
    if flag_ΔmF:
        if len(errors) != 0:
            errors += ", "
        errors += f"ΔmF invalid -> ΔmF = {ΔmF}"

    if flag_ΔFΔmF:
        if len(errors) != 0:
            errors += ", "
        errors += f"ΔF = 0 & ΔmF = 0 invalid"

    if len(errors) != 0:
        errors = f"transition not allowed; {errors}"

    if return_err:
        return not (flag_ΔP | flag_ΔF | flag_ΔmF | flag_ΔFΔmF), errors
    else:
        return not (flag_ΔP | flag_ΔF | flag_ΔmF | flag_ΔFΔmF)


def assert_transition_coupled_allowed(state1, state2, ΔmF_allowed):
    """Check whether the transition is allowed based on the quantum numbers.
    Raises an AssertionError if the transition is not allowed.

    Args:
        state1 (CoupledBasisState): ground CoupledBasisState
        state2 (CoupledBasisState): excited CoupledBasisState

    Returns:
        tuple: (allowed boolean, error message)
    """
    allowed, errors = check_transition_coupled_allowed(state1, state2, ΔmF_allowed)
    assert allowed, errors
    return allowed
