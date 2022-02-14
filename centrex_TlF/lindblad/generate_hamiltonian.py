import copy

import numpy as np
from centrex_TlF.couplings.utils import TransitionSelector
from centrex_TlF.lindblad.utils_compact import (
    compact_symbolic_hamiltonian_indices,
    delete_row_column_symbolic,
)
from centrex_TlF.states.utils import QuantumSelector, get_indices_quantumnumbers
from centrex_TlF.states.utils_compact import compact_QN_coupled_indices
from sympy import Matrix, Symbol, diff
from sympy import exp as symb_exp
from sympy import eye, simplify, solve, symbols, zeros

__all__ = [
    "generate_symbolic_hamiltonian",
    "generate_symbolic_detunings",
    "generate_total_symbolic_hamiltonian",
]


def symbolic_hamiltonian_to_rotating_frame(hamiltonian, QN, H_int, couplings, δs):
    """Transform a symbolic hamiltonian to the rotating frame. Exponential terms
    with the transition frequencies are required to be present in the
    hamiltonian matrix, as well as symbolic energies on the diagonal.

    Args:
        hamiltonian (sympy.Matrix): symbolic hamiltonian
        QN (list/array): list/array of states in the system
        H_int (np.ndarray): numerical hamiltonian, energies only
        couplings (list): list of couplings in system

    Returns:
        sympy.Matrix: symbolic hamiltonian in the rotating frame
    """
    n_states = H_int.shape[0]
    energies = np.diag(hamiltonian)

    # generate t symbol for non-rotating frame
    t = Symbol("t", real=True)

    coupled_states = []
    for i, j in zip(*np.nonzero(hamiltonian)):
        if i < j:
            syms = hamiltonian[i, j].free_symbols
            syms = [s for s in syms if str(s)[0] == "ω"]
            assert len(syms) == 1, f"Too many/few couplings, syms = {syms}"
            coupled_states.append((i, j, syms[0]))

    # solve equations to generate unitary transformation to rotating frame
    A = symbols(f"a:{n_states}")
    Eqns = []
    # generate equations
    for i, j, ω in coupled_states:
        Eqns.append(ω - (A[i] - A[j]))
    # solve system of equations
    sol = solve(Eqns, A)
    # set free parameters to zero in the solution
    free_params = [value for value in A if value not in list(sol.keys())]
    for free_param in free_params:
        for key, val in sol.items():
            sol[key] = val.subs(free_param, 0)

    # generate unitary transformation matrix
    T = eye(*H_int.shape)
    for var in sol.keys():
        ida = int(str(var)[1:])
        T[ida, ida] = symb_exp(1j * sol[var] * t)

    # use unitary matrix to transform to rotating frame
    transformed = T.adjoint() @ hamiltonian @ T - 1j * T.adjoint() @ diff(T, t)
    transformed = simplify(transformed)

    transformed = Matrix(transformed)

    for idc, (δ, coupling) in enumerate(zip(δs, couplings)):
        # generate transition frequency symbol
        ω = Symbol(f"ω{idc}", real=True)
        # get indices of ground and excited states
        idg = QN.index(coupling["ground main"])
        ide = QN.index(coupling["excited main"])
        # transform to δ instead of ω and E
        transformed = transformed.subs(ω, energies[ide] - energies[idg] + δ)

    # substitute level energies for symbolic values
    transformed = transformed.subs(
        [(E, val) for E, val in zip(energies, np.diag(H_int))]
    )

    # set energie difference between excited and ground states to zero
    # should be done automatically when solving for the unitary matrix, not sure
    # why this is not happening currently
    for coupling in couplings:
        idg = QN.index(coupling["ground main"])
        ide = QN.index(coupling["excited main"])
        indices_ground = [QN.index(s) for s in coupling["ground states"]]
        indices_excited = [QN.index(s) for s in coupling["excited states"]]
        g = transformed[idg, idg].subs(
            [(s, 0) for s in transformed[idg, idg].free_symbols]
        )
        e = transformed[ide, ide].subs(
            [(s, 0) for s in transformed[ide, ide].free_symbols]
        )
        for idg in indices_ground:
            transformed[idg, idg] -= g
        for ide in indices_excited:
            transformed[ide, ide] -= e

    return transformed


def generate_symbolic_hamiltonian(QN, H_int, couplings, Ωs, δs, pols):
    n_states = H_int.shape[0]
    # initialize empty hamiltonian
    hamiltonian = zeros(*H_int.shape)
    energies = symbols(f"E:{n_states}")
    hamiltonian += eye(n_states) * np.asarray(energies)

    # generate t symbol for non-rotating frame
    t = Symbol("t", real=True)

    # iterate over couplings
    for idc, (Ω, coupling) in enumerate(zip(Ωs, couplings)):
        # generate transition frequency symbol
        ω = Symbol(f"ω{idc}", real=True)
        # main coupling matrix element
        main_coupling = coupling["main coupling"]
        # iterate over fields (polarizations) in the coupling
        for idf, field in enumerate(coupling["fields"]):
            if pols:
                P = pols[idc]
                if P:
                    P = P[idf]
                    val = (P * Ω / main_coupling) / 2
                    for i, j in zip(*np.nonzero(field["field"])):
                        if i < j:
                            hamiltonian[i, j] += (
                                val * field["field"][i, j] * symb_exp(1j * ω * t)
                            )
                            hamiltonian[j, i] += (
                                val * field["field"][j, i] * symb_exp(-1j * ω * t)
                            )
                else:
                    val = (Ω / main_coupling) / 2
                    for i, j in zip(*np.nonzero(field["field"])):
                        if i < j:
                            hamiltonian[i, j] += (
                                val * field["field"][i, j] * symb_exp(1j * ω * t)
                            )
                            hamiltonian[j, i] += (
                                val * field["field"][j, i] * symb_exp(-1j * ω * t)
                            )
            else:
                val = (Ω / main_coupling) / 2
                for i, j in zip(*np.nonzero(field["field"])):
                    if i < j:
                        hamiltonian[i, j] += (
                            val * field["field"][i, j] * symb_exp(1j * ω * t)
                        )
                        hamiltonian[j, i] += (
                            val * field["field"][j, i] * symb_exp(-1j * ω * t)
                        )

    hamiltonian = simplify(hamiltonian)

    transformed = symbolic_hamiltonian_to_rotating_frame(
        hamiltonian, QN, H_int, couplings, δs
    )
    transformed = Matrix(transformed)

    Ωsᶜ = [Symbol(str(Ω) + "ᶜ", complex=True) for Ω in Ωs]
    for idx in range(n_states):
        for idy in range(0, idx):
            for Ω, Ωᶜ in zip(Ωs, Ωsᶜ):
                transformed[idx, idy] = transformed[idx, idy].subs(Ω, Ωᶜ)

    return transformed


# def generate_symbolic_hamiltonian(QN, H_int, couplings, Ωs = None,  δs = None,
#                                     pols = None):

#     n_states = H_int.shape[0]

#     if not Ωs:
#         Ωs = [Symbol(f'Ω{idx}', complex = True) for idx in range(len(couplings))]
#     elif len(Ωs) != len(couplings):
#         Ωs = [Symbol(f'Ω{idx}', complex = True) for idx in range(len(couplings))]
#         logging.warning("Warning in generate_symbolic_hamiltonian: supplied " +
#                     f"Ωs length does not match # couplings ({len(Ωs)} != {len(couplings)})"
#             )
#     if not δs:
#         δs = [Symbol(f'δ{idx}') for idx in range(len(couplings))]
#     elif len(δs) != len(couplings):
#         δs = [Symbol(f'δ{idx}') for idx in range(len(couplings))]
#         logging.warning("Warning in generate_symbolic_hamiltonian: supplied " +
#                     f"δs length does not match # couplings ({len(δs)} != {len(couplings)})"
#             )

#     # initialize empty Hamiltonian
#     hamiltonian = zeros(*H_int.shape)

#     # add the couplings to the fields
#     for idc, (Ω, coupling) in enumerate(zip(Ωs, couplings)):
#         # check if Ω symbol exists, else create
#         if not Ω:
#             _ = idc
#             while True:
#                 Ω = Symbol(f'Ω{_}', complex = True)
#                 _ += 1
#                 if Ω not in Ωs:
#                     break
#             Ωs[idc] = Ω
#         main_coupling = coupling['main coupling']
#         for idf, field in enumerate(coupling['fields']):
#             if pols:
#                 P = pols[idc]
#                 if P:
#                     P = P[idf]
#                     hamiltonian += (P*Ω/main_coupling)/2 * field['field']
#                 else:
#                     hamiltonian += (Ω/main_coupling)/2 * field['field']
#             else:
#                 hamiltonian += (Ω/main_coupling)/2 * field['field']

#     # add energies
#     hamiltonian += H_int

#     # add detunings to the hamiltonian
#     for idc, (δ, coupling) in enumerate(zip(δs, couplings)):
#         # check if Δ symbol exists, else create
#         if not δ:
#             _ = idc
#             while True:
#                 δ = Symbol(f'δ{_}')
#                 _ += 1
#                 if δ not in δs:
#                     break
#             δs[idc] = δ
#         indices_ground = [QN.index(s) for s in coupling['ground states']]
#         indices_excited = [QN.index(s) for s in coupling['excited states']]
#         idg = QN.index(coupling['ground main'])
#         ide = QN.index(coupling['excited main'])
#         # subtract excited state energy over diagonal for first entry:
#         if idc == 0:
#             hamiltonian -= eye(hamiltonian.shape[0])*hamiltonian[ide,ide]
#         Δ = hamiltonian[ide,ide] - hamiltonian[idg,idg]
#         for idx in indices_ground:
#             hamiltonian[idx, idx] += Δ
#         for idx in indices_ground:
#             hamiltonian[idx,idx] += -δ

#     # ensure hermitian Hamiltonian for complex Ω
#     # complex conjugate Rabi rates
#     Ωsᶜ = [Symbol(str(Ω)+"ᶜ", complex = True) for Ω in Ωs]
#     for idx in range(n_states):
#         for idy in range(0,idx):
#             for Ω,Ωᶜ in zip(Ωs, Ωsᶜ):
#                 hamiltonian[idx,idy] = hamiltonian[idx,idy].subs(Ω, Ωᶜ)


#     return hamiltonian#, symbols


def generate_symbolic_detunings(n_states, detunings):
    detuning = zeros(n_states, n_states)

    if len(detunings) == 1:
        for idd, indices in enumerate(detunings):
            Δ = Symbol("Δ", real=True)
            for idx in indices:
                detuning[idx, idx] += Δ

        symbols = [Symbol("Δ", complex=True)]
    else:
        for idd, indices in enumerate(detunings):
            Δ = Symbol(f"Δ{idd+1}", real=True)
            for idx in indices:
                detuning[idx, idx] += Δ

        symbols = [Symbol(f"Δ{idd+1}", complex=True) for idd in range(len(detunings))]
    return detuning, symbols


def generate_total_symbolic_hamiltonian(
    QN, H_int, couplings, transitions, slice_compact=None, qn_compact=None
):
    if isinstance(transitions[0], TransitionSelector):
        return generate_total_symbolic_hamiltonian_TransitionSelector(
            QN,
            H_int,
            couplings,
            transitions,
            slice_compact=slice_compact,
            qn_compact=qn_compact,
        )
    elif isinstance(transitions[0], dict):
        return generate_total_symbolic_hamiltonian_transitiondict(
            QN,
            H_int,
            couplings,
            transitions,
            slice_compact=slice_compact,
            qn_compact=qn_compact,
        )
    else:
        raise AssertionError(
            "transitions required to be a list of TransitionSelectors or a list of "
            "dicts"
        )


def generate_total_symbolic_hamiltonian_transitiondict(
    QN, H_int, couplings, transitions, slice_compact=None, qn_compact=None
):
    """Generate the total symbolic hamiltonian for the given system

    Args:
        QN (list): states
        H_int (array): internal hamiltonian
        couplings (list): list of dictionaries with all couplings of the system
        transitions (list): list of dictionaries with all transitions of the
                            system
        slice_compact (slice operator, optional): numpy slice operator for which
                                                    part of the system to compact
                                                    to a single state.
                                                    Defaults to None.
        qn_compact (list, optional): list of QuantumSelectors or lists of
                                    QuantumSelectors with each
                                    QuantumSelector containing the quantum
                                    numbers to compact into a single state.
                                    Defaults to None.

    Returns:
        sympy matrix: symbolic hamiltonian
        if qn_compact is provided, also returns the states corresponding to the
        compacted hamiltonian, i.e. ham, QN_compact
    """
    Ωs = [t.get("Ω symbol") for t in transitions]
    Δs = [t.get("Δ symbol") for t in transitions]
    pols = []
    for transition in transitions:
        if not transition.get("polarization symbols"):
            pols.append(None)
        else:
            pols.append(transition["polarization symbols"])

    H_symbolic = generate_symbolic_hamiltonian(QN, H_int, couplings, Ωs, Δs, pols)

    if slice_compact:
        H_symbolic = delete_row_column_symbolic(H_symbolic, slice_compact)
    elif qn_compact:
        if isinstance(qn_compact, QuantumSelector):
            qn_compact = [qn_compact]
        QN_compact = copy.deepcopy(QN)
        for qnc in qn_compact:
            indices_compact = get_indices_quantumnumbers(qnc, QN_compact)
            QN_compact = compact_QN_coupled_indices(QN_compact, indices_compact)
            H_symbolic = compact_symbolic_hamiltonian_indices(
                H_symbolic, indices_compact
            )
        return H_symbolic, QN_compact

    return H_symbolic


def generate_total_symbolic_hamiltonian_TransitionSelector(
    QN, H_int, couplings, transitions, slice_compact=None, qn_compact=None
):
    """Generate the total symbolic hamiltonian for the given system

    Args:
        QN (list): states
        H_int (array): internal hamiltonian
        couplings (list): list of dictionaries with all couplings of the system
        transitions (list): list of dictionaries with all transitions of the
                            system
        slice_compact (slice operator, optional): numpy slice operator for which
                                                    part of the system to compact
                                                    to a single state.
                                                    Defaults to None.
        qn_compact (list, optional): list of QuantumSelectors or lists of
                                    QuantumSelectors with each
                                    QuantumSelector containing the quantum
                                    numbers to compact into a single state.
                                    Defaults to None.

    Returns:
        sympy matrix: symbolic hamiltonian
        if qn_compact is provided, also returns the states corresponding to the
        compacted hamiltonian, i.e. ham, QN_compact
    """
    Ωs = [t.Ω for t in transitions]
    Δs = [t.δ for t in transitions]
    pols = []
    for transition in transitions:
        if not transition.polarization_symbols:
            pols.append(None)
        else:
            pols.append(transition.polarization_symbols)

    H_symbolic = generate_symbolic_hamiltonian(QN, H_int, couplings, Ωs, Δs, pols)
    if slice_compact:
        H_symbolic = delete_row_column_symbolic(H_symbolic, slice_compact)
    elif qn_compact is not None:
        if isinstance(qn_compact, QuantumSelector):
            qn_compact = [qn_compact]
        QN_compact = copy.deepcopy(QN)
        for qnc in qn_compact:
            indices_compact = get_indices_quantumnumbers(qnc, QN_compact)
            QN_compact = compact_QN_coupled_indices(QN_compact, indices_compact)
            H_symbolic = compact_symbolic_hamiltonian_indices(
                H_symbolic, indices_compact
            )
        return H_symbolic, QN_compact

    return H_symbolic
