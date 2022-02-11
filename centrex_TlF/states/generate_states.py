import copy
import numpy as np
from typing import Union
from centrex_TlF.constants.constants import I_Tl, I_F
from centrex_TlF.states.utils import QuantumSelector, get_unique_basisstates, parity_X
from centrex_TlF.states.states import UncoupledBasisState, CoupledBasisState

__all__ = [
    "generate_uncoupled_states_ground",
    "generate_uncoupled_states_excited",
    "generate_coupled_states_ground",
    "generate_coupled_states_excited",
    "generate_coupled_states_base",
    "generate_coupled_states_ground_X",
    "generate_coupled_states_excited_B",
]


def generate_uncoupled_states_ground(Js):
    # convert J to int(J); np.int with (-1)**J throws an exception for negative J
    QN = np.array(
        [
            UncoupledBasisState(
                int(J),
                mJ,
                I_Tl,
                m1,
                I_F,
                m2,
                Omega=0,
                P=parity_X(J),
                electronic_state="X",
            )
            for J in Js
            for mJ in range(-J, J + 1)
            for m1 in np.arange(-I_Tl, I_Tl + 1)
            for m2 in np.arange(-I_F, I_F + 1)
        ]
    )
    return QN


def generate_uncoupled_states_excited(Js, Ωs=[-1, 1]):
    QN = np.array(
        [
            UncoupledBasisState(J, mJ, I_Tl, m1, I_F, m2, Omega=Ω, electronic_state="B")
            for Ω in Ωs
            for J in Js
            for mJ in range(-J, J + 1)
            for m1 in np.arange(-I_Tl, I_Tl + 1)
            for m2 in np.arange(-I_F, I_F + 1)
        ]
    )
    return QN


def generate_coupled_states_ground(Js):
    QN = np.array(
        [
            CoupledBasisState(
                F, mF, F1, J, I_F, I_Tl, electronic_state="X", P=parity_X(J), Omega=0
            )
            for J in Js
            for F1 in np.arange(np.abs(J - I_F), J + I_F + 1)
            for F in np.arange(np.abs(F1 - I_Tl), F1 + I_Tl + 1)
            for mF in np.arange(-F, F + 1)
        ]
    )
    return QN


def generate_coupled_states_base(qn_selector: QuantumSelector) -> np.ndarray:
    """generate CoupledBasisStates for the quantum numbers given by qn_selector

    Args:
        qn_selector (QuantumSelector): quantum numbers to use to generate the 
                                        CoupledBasisStates

    Returns:
        np.ndarray: array of CoupledBasisStates for the excited state
    """
    assert qn_selector.P is not None, "function requires a parity to be set"
    assert (
        qn_selector.J is not None
    ), "function requires a rotational quantum number to be set"
    assert (
        qn_selector.electronic is not None
    ), "function requires electronic state to be set"
    assert qn_selector.Ω is not None, "function requires Ω to be set"

    # generate all combinations
    quantum_numbers = []
    for par in ["J", "F1", "F", "mF", "electronic", "P", "Ω"]:
        par = getattr(qn_selector, par)
        quantum_numbers.append(
            [par] if not isinstance(par, (list, tuple, np.ndarray)) else par
        )

    QN = []
    # the worst nested loops I've ever created
    Js, F1s, Fs, mFs, estates, Ps, Ωs = quantum_numbers
    for estate in estates:
        for J in Js:
            F1_allowed = np.arange(np.abs(J - I_F), J + I_F + 1)
            F1sl = F1s if F1s[0] is not None else F1_allowed
            for F1 in F1sl:
                if F1 not in F1_allowed:
                    continue
                Fs_allowed = np.arange(np.abs(F1 - I_Tl), F1 + I_Tl + 1)
                Fsl = Fs if Fs[0] is not None else Fs_allowed
                for F in Fsl:
                    if F not in Fs_allowed:
                        continue
                    mF_allowed = np.arange(-F, F + 1)
                    mFsl = mFs if mFs[0] is not None else mF_allowed
                    for mF in mFsl:
                        if mF not in mF_allowed:
                            continue
                        for P in Ps:
                            P = P if not callable(P) else P(J)
                            for Ω in Ωs:
                                QN.append(
                                    CoupledBasisState(
                                        F,
                                        mF,
                                        F1,
                                        J,
                                        I_F,
                                        I_Tl,
                                        electronic_state=estate,
                                        P=P,
                                        Ω=Ω,
                                    )
                                )
    return np.asarray(QN)


def generate_coupled_states_ground_X(
    qn_selector: Union[QuantumSelector, list, np.ndarray]
) -> np.ndarray:
    """generate ground X state CoupledBasisStates for the quantum numbers given 
    by qn_selector

    Args:
        qn_selector (Union[QuantumSelector, list, np.ndarray]): QuantumSelector
            or list/array of QuantumSelectors to use for generating the 
            CoupledBasisStates

    Returns:
        np.ndarray: array of CoupledBasisStates for the excited state
    """
    if isinstance(qn_selector, QuantumSelector):
        qns = copy.copy(qn_selector)
        qns.Ω = 0
        qns.P = parity_X
        qns.electronic = "X"
        return generate_coupled_states_base(qns)
    elif isinstance(qn_selector, (list, np.ndarray)):
        coupled_states = []
        for qns in qn_selector:
            qns = copy.copy(qns)
            qns.Ω = 0
            qns.P = parity_X
            qns.electronic = "X"
            coupled_states.append(generate_coupled_states_base(qns))
        return get_unique_basisstates(np.concatenate(coupled_states))
    else:
        raise AssertionError(
            "qn_selector required to be of type QuantumSelector, list or np.ndarray"
        )


def generate_coupled_states_excited_B(
    qn_selector: Union[QuantumSelector, list, np.ndarray]
) -> np.ndarray:
    """generate excited B state CoupledBasisStates for the quantum numbers given 
    by qn_selector

    Args:
        qn_selector (Union[QuantumSelector, list, np.ndarray]): QuantumSelector
            or list/array of QuantumSelectors to use for generating the 
            CoupledBasisStates

    Returns:
        np.ndarray: array of CoupledBasisStates for the excited state
    """
    if isinstance(qn_selector, QuantumSelector):
        qns = copy.copy(qn_selector)
        qns.Ω = 1
        qns.electronic = "B"
        return generate_coupled_states_base(qns)
    elif isinstance(qn_selector, (list, np.ndarray)):
        coupled_states = []
        for qns in qn_selector:
            qns = copy.copy(qns)
            qns.Ω = 1
            qns.electronic = "B"
            coupled_states.append(generate_coupled_states_base(qns))
        return get_unique_basisstates(np.concatenate(coupled_states))
    else:
        raise AssertionError(
            "qn_selector required to be of type QuantumSelector, list or np.ndarray"
        )


def generate_coupled_states_excited(Js, Fs=None, F1s=None, Ps=None):
    if not Fs:
        if not Ps:
            Ps = [+1]
        QN = np.array(
            [
                CoupledBasisState(
                    F, mF, F1, J, I_F, I_Tl, electronic_state="B", P=P, Omega=1
                )
                for J in Js
                for F1 in np.arange(np.abs(J - I_F), J + I_F + 1)
                for F in np.arange(np.abs(F1 - I_Tl), F1 + I_Tl + 1)
                for mF in np.arange(-F, F + 1)
                for P in Ps
            ]
        )
    else:
        assert None not in [
            Fs,
            F1s,
            Ps,
        ], "need to supply lists of F, F1 or P if one of them is used as an input parameter"
        QN = np.array(
            [
                CoupledBasisState(
                    F, mF, F1, J, I_F, I_Tl, electronic_state="B", P=P, Omega=1
                )
                for J, F1, F in zip(Js, F1s, Fs)
                for mF in np.arange(-F, F + 1)
                for P in Ps
            ]
        )
    return QN
