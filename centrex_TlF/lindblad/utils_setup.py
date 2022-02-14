import logging
from dataclasses import dataclass

import numpy as np
from centrex_TlF.couplings.collapse import collapse_matrices
from centrex_TlF.couplings.coupling_matrix import (
    generate_coupling_field,
    generate_coupling_field_automatic,
)
from centrex_TlF.hamiltonian.generate_reduced_hamiltonian import (
    generate_total_reduced_hamiltonian,
)
from centrex_TlF.lindblad.generate_hamiltonian import (
    generate_total_symbolic_hamiltonian,
)
from centrex_TlF.lindblad.generate_julia_code import (
    generate_preamble,
    system_of_equations_to_lines,
)
from centrex_TlF.lindblad.generate_system_of_equations import (
    generate_system_of_equations_symbolic,
)
from centrex_TlF.lindblad.utils_decay import (
    add_decays_C_arrays,
    add_levels_symbolic_hamiltonian,
    add_states_QN,
)
from centrex_TlF.lindblad.utils_julia import generate_ode_fun_julia, initialize_julia
from centrex_TlF.states.generate_states import (
    generate_coupled_states_excited_B,
    generate_coupled_states_ground_X,
)
from julia import Main

__all__ = ["generate_OBE_system", "setup_OBE_system_julia"]


@dataclass
class OBESystem:
    ground: np.ndarray
    excited: np.ndarray
    QN: np.ndarray
    H_int: np.ndarray
    V_ref_int: np.ndarray
    couplings: list
    H_symbolic: np.ndarray
    C_array: np.ndarray
    system: np.ndarray
    code_lines: list
    full_output: bool = False
    preamble: str = ""
    QN_original: np.ndarray = None
    decay_channels: np.ndarray = None


def generate_OBE_system(
    system_parameters, transitions, qn_compact=None, decay_channels=None, verbose=False
):
    """Convenience function for generating the symbolic OBE system of equations
    and Julia code.

    Args:
        system_parameters (SystemParameters): dataclass holding system parameters

        transitions (list): list of TransitionSelectors defining the transitions
                            used in the OBE system.
        qn_compact (QuantumSelector): dataclass specifying a subset of states to
                                        select based on the quantum numbers
        decay_channels (DecayChannel): dataclass specifying the decay channel to
                                        add
        verbose (bool, optional): Log progress to INFO. Defaults to False.

    Returns:
        OBESystem: dataclass designed to hold the generated values
                    ground, exxcited, QN, H_int, V_ref_int, couplings, H_symbolic,
                    C_array, system, code_lines
    """
    # values above and below of excited J states to include in hamiltonian
    # to take into account excited state mixing
    # default input None does min(J)-1 and max(J)+1
    Jmin = None
    Jmax = None
    rtol = None
    if verbose:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        logger.info("generate_OBE_system: 1/6 -> Generating the reduced Hamiltonian")
    (
        ground_states,
        excited_states,
        QN,
        H_int,
        V_ref_int,
    ) = generate_total_reduced_hamiltonian(
        ground_states_approx=generate_coupled_states_ground_X(system_parameters.ground),
        excited_states_approx=generate_coupled_states_excited_B(
            system_parameters.excited
        ),
        Jmin=Jmin,
        Jmax=Jmax,
        rtol=rtol,
    )
    if verbose:
        logger.info(
            "generate_OBE_system: 2/6 -> "
            "Generating the couplings corresponding to the transitions"
        )
    couplings = []
    for transition in transitions:
        if transition.ground_main is not None and transition.excited_main is not None:
            couplings.append(
                generate_coupling_field(
                    transition.ground_main,
                    transition.excited_main,
                    transition.ground,
                    transition.excited,
                    H_int,
                    QN,
                    V_ref_int,
                    pol_vec=transition.polarizations,
                    pol_main=transition.polarizations[0],
                    nprocs=system_parameters.nprocs,
                )
            )
        else:
            couplings.append(
                generate_coupling_field_automatic(
                    transition.ground,
                    transition.excited,
                    H_int,
                    QN,
                    V_ref_int,
                    pol_vec=transition.polarizations,
                    nprocs=system_parameters.nprocs,
                )
            )

    if verbose:
        logger.info("generate_OBE_system: 3/6 -> Generating the symbolic Hamiltonian")
    if qn_compact is not None:
        H_symbolic, QN_compact = generate_total_symbolic_hamiltonian(
            QN, H_int, couplings, transitions, qn_compact=qn_compact
        )
    else:
        H_symbolic = generate_total_symbolic_hamiltonian(
            QN, H_int, couplings, transitions
        )

    if verbose:
        logger.info("generate_OBE_system: 4/6 -> Generating the collapse matrices")
    C_array = collapse_matrices(
        QN,
        ground_states,
        excited_states,
        gamma=system_parameters.Γ,
        qn_compact=qn_compact,
    )

    if decay_channels is not None:
        if not isinstance(decay_channels, (list, np.ndarray)):
            decay_channels = np.ndarray(decay_channels)
        if qn_compact is not None:
            indices, H_symbolic = add_levels_symbolic_hamiltonian(
                H_symbolic, decay_channels, QN_compact, excited_states
            )
            QN_compact = add_states_QN(decay_channels, QN_compact, indices)
            C_array = add_decays_C_arrays(
                decay_channels, indices, QN_compact, C_array, system_parameters.Γ
            )
        else:
            indices, H_symbolic = add_levels_symbolic_hamiltonian(
                H_symbolic, decay_channels, QN, excited_states
            )
            QN = add_states_QN(decay_channels, QN_compact, indices)
            C_array = add_decays_C_arrays(
                decay_channels, indices, QN, C_array, system_parameters.Γ
            )
    if verbose:
        logger.info(
            "generate_OBE_system: 5/6 -> Transforming the Hamiltonian and collapse "
            "matrices into a symbolic system of equations"
        )
    system = generate_system_of_equations_symbolic(
        H_symbolic, C_array, progress=False, fast=True
    )
    if verbose:
        logger.info(
            "generate_OBE_system: 6/6 -> Generating Julia code representing the system "
            "of equations"
        )
        logging.basicConfig(level=logging.WARNING)
    code_lines = system_of_equations_to_lines(system, nprocs=system_parameters.nprocs)
    if qn_compact is not None:
        QN_original = QN
        QN = QN_compact
    else:
        QN_original = None
    obe_system = OBESystem(
        QN=QN,
        ground=ground_states,
        excited=excited_states,
        couplings=couplings,
        H_symbolic=H_symbolic,
        H_int=H_int,
        V_ref_int=V_ref_int,
        C_array=C_array,
        system=system,
        code_lines=code_lines,
        QN_original=QN_original,
        decay_channels=decay_channels,
    )
    return obe_system


def setup_OBE_system_julia(
    system_parameters,
    ode_parameters,
    transitions,
    qn_compact=None,
    full_output=False,
    decay_channels=None,
    verbose=False,
):
    """Convenience function for generating the OBE system and initializing it in
    Julia

    Args:
        system_parameters (SystemParameters): dataclass holding the system
                                                parameters, e.g. Γ,
                                                (laser) ground states,
                                                (laser) excited states
        ode_parameters (odeParameters): dataclass containing the ode parameters.
                                        e.g. Ω, δ, vz, ..., etc.
        transitions (TransitionSelector): object containing all information
                                            required to generate the coupling
                                            matrices and symbolic matrix for
                                            each transition
        qn_compact (QuantumSelector): dataclass specifying a subset of states to
                                        select based on the quantum numbers
        full_output (bool, optional): Returns all matrices, states etc. if True,
                                        Returns only QN if False.
                                        Defaults to False.
        verbose (bool, optional): Log progress to INFO. Defaults to False.

    Returns:
        full_output == True:
            list: list of states in system
        full_output == False:
            OBESystem: dataclass designed to hold the generated values
                        ground, exxcited, QN, H_int, V_ref_int, couplings,
                        H_symbolic, C_array, system, code_lines, preamble
    """
    obe_system = generate_OBE_system(
        system_parameters,
        transitions,
        qn_compact=qn_compact,
        decay_channels=decay_channels,
        verbose=verbose,
    )
    obe_system.full_output = full_output
    if verbose:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        logger.info("setup_OBE_system_julia: 1/3 -> Generating the preamble")
    obe_system.preamble = generate_preamble(ode_parameters, transitions)

    if verbose:
        logger.info(
            "setup_OBE_system_julia: 2/3 -> Initializing Julia on "
            f"{system_parameters.nprocs} cores"
        )
    initialize_julia(nprocs=system_parameters.nprocs)

    if verbose:
        logger.info(
            "setup_OBE_system_julia: 3/3 -> Defining the ODE equation and parameters in"
            " Julia"
        )
        logging.basicConfig(level=logging.WARNING)
    generate_ode_fun_julia(obe_system.preamble, obe_system.code_lines)
    Main.eval(f"@everywhere Γ = {system_parameters.Γ}")
    ode_parameters.generate_p_julia()
    if not full_output:
        return obe_system.QN
    else:
        return obe_system
