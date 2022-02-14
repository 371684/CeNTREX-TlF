import multiprocessing
import sqlite3
from pathlib import Path

import numpy as np
from centrex_TlF.couplings.matrix_elements import (
    calculate_ED_ME_mixed_state,
    generate_ED_ME_mixed_state,
)
from centrex_TlF.couplings.utils import generate_D, select_main_states
from centrex_TlF.couplings.utils_multiprocessing import multi_coupling_matrix
from centrex_TlF.couplings.utils_sqlite import check_states_in_ED_ME_coupled
from centrex_TlF.states.utils import check_approx_state_exact_state, find_exact_states
from centrex_TlF.transitions.utils import assert_transition_coupled_allowed

__all__ = [
    "calculate_coupling_matrix",
    "generate_coupling_field",
    "generate_coupling_matrix",
    "calculate_coupling_field",
    "generate_coupling_field_automatic",
    "calculate_coupling_field_automatic",
]


def generate_coupling_matrix(
    QN,
    ground_states,
    excited_states,
    pol_vec=np.array([0, 0, 1]),
    reduced=False,
    normalize_pol=True,
    nprocs=1,
):
    """generate optical coupling matrix for given ground and excited states
    Checks if couplings are already pre-cached, otherwise falls back to
    calculate_coupling_matrix.

    Args:
        QN (list): list of basis states
        ground_states (list): list of ground states coupling to excited states
        excited_states (list): list of excited states
        pol_vec (np.ndarray, optional): polarization vector. Defaults to
                                        np.array([0,0,1]).
        reduced (bool, optional): [description]. Defaults to False.
        normalize_pol (bool, optional): Normalize the polarization vector.
                                        Defaults to True.
        nrpocs (int): # processes to use for multiprocessing

    Returns:
        np.ndarray: optical coupling matrix
    """
    assert isinstance(QN, list), "QN required to be of type list"

    # check if states are pre-cached
    Jg = np.unique(
        [s[1].J for state in QN for s in state if s[1].electronic_state == "X"]
    )
    Je = np.unique(
        [s[1].J for state in QN for s in state if s[1].electronic_state == "B"]
    )
    pre_cached = check_states_in_ED_ME_coupled(Jg, Je, pol_vec)

    if pre_cached:
        # connect to sqlite3 database on file, not used when multiprocessing
        path = Path(__file__).parent.parent / "pre_calculated"
        db = path / "matrix_elements.db"

        con = sqlite3.connect(db)
        cur = con.cursor()
        cur.execute("PRAGMA synchronous = OFF")
        cur.execute("PRAGMA journal_mode = MEMORY")
        # initialize the coupling matrix
        H = np.zeros((len(QN), len(QN)), dtype=complex)

        # start looping over ground and excited states
        for ground_state in ground_states:
            i = QN.index(ground_state)
            for excited_state in excited_states:
                j = QN.index(excited_state)

                # calculate matrix element and add it to the Hamiltonian
                H[i, j] = generate_ED_ME_mixed_state(
                    ground_state,
                    excited_state,
                    pol_vec=pol_vec,
                    reduced=reduced,
                    normalize_pol=normalize_pol,
                    con=con,
                )
        con.close()
        H = H + H.conj().T
        return H
    else:
        return calculate_coupling_matrix(
            QN, ground_states, excited_states, pol_vec, reduced, nprocs
        )


def calculate_coupling_matrix(
    QN,
    ground_states,
    excited_states,
    pol_vec=np.array([0, 0, 1]),
    reduced=False,
    normalize_pol=True,
    nprocs=1,
):
    """calculate optical coupling matrix for given ground and excited states

    Args:
        QN (list): list of basis states
        ground_states (list): list of ground states coupling to excited states
        excited_states (list): list of excited states
        pol_vec (np.ndarray, optional): polarization vector. Defaults to
                                        np.array([0,0,1]).
        reduced (bool, optional): [description]. Defaults to False.
        normalize_pol (bool, optional): Normalize the polarization vector.
                                        Defaults to True.
        nrpocs (int): # processes to use for multiprocessing

    Returns:
        np.ndarray: optical coupling matrix
    """
    assert isinstance(QN, list), "QN required to be of type list"

    if nprocs > 1:
        with multiprocessing.Pool(nprocs) as pool:
            result = pool.starmap(
                multi_coupling_matrix,
                [(QN, gs, excited_states, pol_vec, reduced) for gs in ground_states],
            )
        H = np.sum(result, axis=0)
    else:
        # initialize the coupling matrix
        H = np.zeros((len(QN), len(QN)), dtype=complex)

        # start looping over ground and excited states
        for ground_state in ground_states:
            i = QN.index(ground_state)
            for excited_state in excited_states:
                j = QN.index(excited_state)

                # calculate matrix element and add it to the Hamiltonian
                H[i, j] = calculate_ED_ME_mixed_state(
                    ground_state,
                    excited_state,
                    pol_vec=pol_vec,
                    reduced=reduced,
                    normalize_pol=normalize_pol,
                )

    # make H hermitian
    H = H + H.conj().T

    return H


def generate_coupling_field(
    ground_main_approx,
    excited_main_approx,
    ground_states_approx,
    excited_states_approx,
    H_rot,
    QN,
    V_ref,
    pol_main=np.array([0, 0, 1]),
    pol_vec=[],
    relative_coupling=1e-3,
    absolute_coupling=1e-6,
    nprocs=2,
):
    ground_states = find_exact_states(ground_states_approx, H_rot, QN, V_ref=V_ref)
    excited_states = find_exact_states(excited_states_approx, H_rot, QN, V_ref=V_ref)
    ground_main = find_exact_states([ground_main_approx], H_rot, QN, V_ref=V_ref)[0]
    excited_main = find_exact_states([excited_main_approx], H_rot, QN, V_ref=V_ref)[0]

    check_approx_state_exact_state(ground_main_approx, ground_main)
    check_approx_state_exact_state(excited_main_approx, excited_main)
    ME_main = generate_ED_ME_mixed_state(excited_main, ground_main, pol_vec=pol_main)

    assert_transition_coupled_allowed(
        ground_main.find_largest_component(),
        excited_main.find_largest_component(),
        ΔmF_allowed=0 if pol_main[2] != 0 else 1,
    )
    # assert_msg = f"main coupling element small, {ME_main:.2e}" + \
    #               ", check states and/or polarization"
    # assert np.abs(ME_main) > 1e-2, assert_msg

    couplings = {
        "ground main": ground_main,
        "excited main": excited_main,
        "main coupling": ME_main,
        "ground states": ground_states,
        "excited states": excited_states,
        "D": generate_D(H_rot, QN, ground_main, excited_main, excited_states),
        "fields": [],
    }

    for pol in pol_vec:
        coupling = generate_coupling_matrix(
            QN, ground_states, excited_states, pol_vec=pol, reduced=False, nprocs=nprocs
        )

        coupling[np.abs(coupling) < relative_coupling * np.max(np.abs(coupling))] = 0
        coupling[np.abs(coupling) < absolute_coupling] = 0
        d = {"pol": pol, "field": coupling}
        couplings["fields"].append(d)
    return couplings


def calculate_coupling_field(
    ground_main_approx,
    excited_main_approx,
    ground_states_approx,
    excited_states_approx,
    H_rot,
    QN,
    V_ref,
    pol_main=np.array([0, 0, 1]),
    pol_vec=[],
    relative_coupling=1e-3,
    absolute_coupling=1e-6,
    nprocs=2,
):
    ground_states = find_exact_states(ground_states_approx, H_rot, QN, V_ref=V_ref)
    excited_states = find_exact_states(excited_states_approx, H_rot, QN, V_ref=V_ref)
    ground_main = find_exact_states([ground_main_approx], H_rot, QN, V_ref=V_ref)[0]
    excited_main = find_exact_states([excited_main_approx], H_rot, QN, V_ref=V_ref)[0]

    check_approx_state_exact_state(ground_main_approx, ground_main)
    check_approx_state_exact_state(excited_main_approx, excited_main)
    ME_main = calculate_ED_ME_mixed_state(excited_main, ground_main, pol_vec=pol_main)
    assert_transition_coupled_allowed(
        ground_main.find_largest_component(),
        excited_main.find_largest_component(),
        ΔmF_allowed=0 if pol_main[2] != 0 else 1,
    )
    # assert_msg = f"main coupling element small, {ME_main:.2e}" + \
    #               ", check states and/or polarization"
    # assert np.abs(ME_main) > 1e-2, assert_msg

    couplings = {
        "ground main": ground_main,
        "excited main": excited_main,
        "main coupling": ME_main,
        "ground states": ground_states,
        "excited states": excited_states,
        "D": generate_D(H_rot, QN, ground_main, excited_main, excited_states),
        "fields": [],
    }

    for pol in pol_vec:
        coupling = calculate_coupling_matrix(
            QN, ground_states, excited_states, pol_vec=pol, reduced=False, nprocs=nprocs
        )

        coupling[np.abs(coupling) < relative_coupling * np.max(np.abs(coupling))] = 0
        coupling[np.abs(coupling) < absolute_coupling] = 0
        d = {"pol": pol, "field": coupling}
        couplings["fields"].append(d)

    return couplings


# automatically generate main coupling


def generate_coupling_field_automatic(
    ground_states_approx,
    excited_states_approx,
    H_rot,
    QN,
    V_ref,
    pol_vec=[],
    relative_coupling=1e-3,
    absolute_coupling=1e-6,
    nprocs=2,
):
    """Generate the coupling fields for a transition for one or multiple
    polarizations. Uses pre-cached values where possible.

    Args:
        ground_states_approx (list): list of approximate ground states
        excited_states_approx (list): list of approximate excited states
        H_rot (np.ndarray): System hamiltonian in the rotational frame
        QN (list): list of states in the system
        V_ref ([type]): [description]
        pol_vec (list, optional): list of polarizations. Defaults to [].
        relative_coupling (float, optional): minimum relative coupling, set
                                            smaller coupling to zero.
                                            Defaults to 1e-3.
        absolute_coupling (float, optional): minimum absolute coupling, set
                                            smaller couplings to zero.
                                            Defaults to 1e-6.
        nprocs (int, optional): number of processes to employ. Defaults to 2.

    Returns:
        dictionary: dictionary with the coupling field information.
                    Dict entries:
                        ground main: main ground state
                        excited main: main excited state
                        main coupling: coupling strenght between main_ground
                                        and main_excited
                        ground_states: ground states in coupling
                        excited states: excited_states in coupling
                        fields: list of dictionaries, one for each polarization,
                                containing the polarization and coupling field
    """
    assert len(pol_vec) != 0, "define polarization vectors for transitions"
    pol_main = pol_vec[0]
    ground_main_approx, excited_main_approx = select_main_states(
        ground_states_approx, excited_states_approx, pol_main
    )

    ground_states = find_exact_states(ground_states_approx, H_rot, QN, V_ref=V_ref)
    excited_states = find_exact_states(excited_states_approx, H_rot, QN, V_ref=V_ref)
    ground_main = find_exact_states([ground_main_approx], H_rot, QN, V_ref=V_ref)[0]
    excited_main = find_exact_states([excited_main_approx], H_rot, QN, V_ref=V_ref)[0]

    check_approx_state_exact_state(ground_main_approx, ground_main)
    check_approx_state_exact_state(excited_main_approx, excited_main)
    ME_main = generate_ED_ME_mixed_state(excited_main, ground_main, pol_vec=pol_main)
    assert_transition_coupled_allowed(
        ground_main.find_largest_component(),
        excited_main.find_largest_component(),
        ΔmF_allowed=0 if pol_main[2] != 0 else 1,
    )
    # assert_msg = f"main coupling element small, {ME_main:.2e}" + \
    #               ", check states and/or polarization"
    # assert np.abs(ME_main) > 1e-2, assert_msg

    couplings = {
        "ground main": ground_main,
        "excited main": excited_main,
        "main coupling": ME_main,
        "ground states": ground_states,
        "excited states": excited_states,
        "fields": [],
    }

    for pol in pol_vec:
        coupling = generate_coupling_matrix(
            QN,
            ground_states,
            excited_states,
            pol_vec=pol,
            reduced=False,
            normalize_pol=True,
            nprocs=nprocs,
        )

        coupling[np.abs(coupling) < relative_coupling * np.max(np.abs(coupling))] = 0
        coupling[np.abs(coupling) < absolute_coupling] = 0
        d = {"pol": pol, "field": coupling}
        couplings["fields"].append(d)
    return couplings


def calculate_coupling_field_automatic(
    ground_states_approx,
    excited_states_approx,
    H_rot,
    QN,
    V_ref,
    pol_vec=[],
    relative_coupling=1e-3,
    absolute_coupling=1e-6,
    nprocs=2,
):
    """Calculate the coupling fields for a transition for one or multiple
    polarizations.

    Args:
        ground_states_approx (list): list of approximate ground states
        excited_states_approx (list): list of approximate excited states
        H_rot (np.ndarray): System hamiltonian in the rotational frame
        QN (list): list of states in the system
        V_ref ([type]): [description]
        pol_vec (list, optional): list of polarizations. Defaults to [].
        relative_coupling (float, optional): minimum relative coupling, set
                                            smaller coupling to zero.
                                            Defaults to 1e-3.
        absolute_coupling (float, optional): minimum absolute coupling, set
                                            smaller couplings to zero.
                                            Defaults to 1e-6.
        nprocs (int, optional): number of processes to employ. Defaults to 2.

    Returns:
        dictionary: dictionary with the coupling field information.
                    Dict entries:
                        ground main: main ground state
                        excited main: main excited state
                        main coupling: coupling strenght between main_ground
                                        and main_excited
                        ground_states: ground states in coupling
                        excited states: excited_states in coupling
                        fields: list of dictionaries, one for each polarization,
                                containing the polarization and coupling field
    """
    assert len(pol_vec) != 0, "define polarization vectors for transitions"
    pol_main = pol_vec[0]
    ground_main_approx, excited_main_approx = select_main_states(
        ground_states_approx, excited_states_approx, pol_main
    )
    ground_states = find_exact_states(ground_states_approx, H_rot, QN, V_ref=V_ref)
    excited_states = find_exact_states(excited_states_approx, H_rot, QN, V_ref=V_ref)
    ground_main = find_exact_states([ground_main_approx], H_rot, QN, V_ref=V_ref)[0]
    excited_main = find_exact_states([excited_main_approx], H_rot, QN, V_ref=V_ref)[0]

    check_approx_state_exact_state(ground_main_approx, ground_main)
    check_approx_state_exact_state(excited_main_approx, excited_main)
    assert_transition_coupled_allowed(
        ground_main.find_largest_component(),
        excited_main.find_largest_component(),
        ΔmF_allowed=0 if pol_main[2] != 0 else 1,
    )
    ME_main = calculate_ED_ME_mixed_state(excited_main, ground_main, pol_vec=pol_main)

    # assert_msg = f"main coupling element small, {ME_main:.2e}" + \
    #               ", check states and/or polarization"
    # assert np.abs(ME_main) > 1e-2, assert_msg

    couplings = {
        "ground main": ground_main,
        "excited main": excited_main,
        "main coupling": ME_main,
        "ground states": ground_states,
        "excited states": excited_states,
        "fields": [],
    }

    for pol in pol_vec:
        coupling = calculate_coupling_matrix(
            QN,
            ground_states,
            excited_states,
            pol_vec=pol,
            reduced=False,
            normalize_pol=True,
            nprocs=nprocs,
        )

        coupling[np.abs(coupling) < relative_coupling * np.max(np.abs(coupling))] = 0
        coupling[np.abs(coupling) < absolute_coupling] = 0
        d = {"pol": pol, "field": coupling}
        couplings["fields"].append(d)

    return couplings
