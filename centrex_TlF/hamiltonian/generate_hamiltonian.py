import centrex_TlF
import logging
import numpy as np
from tqdm import tqdm
import multiprocessing
from pathlib import Path
from centrex_TlF.hamiltonian.utils_multiprocessing import multi_HMatElems
from centrex_TlF.hamiltonian.hamiltonian_terms_uncoupled import (
    Hff_X,
    HSx,
    HSy,
    HSz,
    HZx_X,
    HZy_X,
    HZz_X,
    Hff_X_alt,
)
from centrex_TlF.hamiltonian.hamiltonian_B_terms_coupled import (
    Hrot_B,
    H_mhf_Tl,
    H_mhf_F,
    H_LD,
    H_cp1_Tl,
    H_c_Tl,
    HZz_B,
)
from centrex_TlF.hamiltonian.utils_sqlite import (
    retrieve_uncoupled_hamiltonian_X_sqlite,
    retrieve_coupled_hamiltonian_B_sqlite,
    check_states_uncoupled_hamiltonian_X,
    check_states_coupled_hamiltonian_B,
)

__all__ = [
    "generate_uncoupled_hamiltonian_X",
    "generate_coupled_hamiltonian_B",
    "calculate_uncoupled_hamiltonian_X",
    "calculate_coupled_hamiltonian_B",
]


def HMatElems(H, QN, progress=False, nprocs=1):
    if nprocs > 1:
        with multiprocessing.Pool(nprocs) as pool:
            result = pool.starmap(
                multi_HMatElems, [(H.__name__, i, a, QN) for i, a in enumerate(QN)]
            )
            result = np.sum(result, axis=0)
    else:
        result = np.zeros((len(QN), len(QN)), dtype=complex)
        for i, a in tqdm(enumerate(QN), total=len(QN), disable=not progress):
            for j in range(i, len(QN)):
                b = QN[j]
                val = (1 * a) @ H(b)
                result[i, j] = val
                if i != j:
                    result[j, i] = np.conjugate(val)
    return result


def generate_uncoupled_hamiltonian_X(QN, nprocs=1):
    """
    Generate the uncoupled X state hamiltonian for the supplied set of 
    basis states.
    Retrieved from a pre-calculated sqlite3 database

    Args:
        QN (array): array of UncoupledBasisStates

    Returns:
        dict: dictionary with all X state hamiltonian terms
    """
    for qn in QN:
        assert qn.isUncoupled, "supply list with UncoupledBasisStates"

    pre_cached = check_states_uncoupled_hamiltonian_X(QN)

    if pre_cached:
        path = Path(__file__).parent.parent / "pre_calculated"
        db = path / "uncoupled_hamiltonian_X.db"

        return retrieve_uncoupled_hamiltonian_X_sqlite(QN, db)
    else:
        logging.warning(
            "X state Hamiltonian not pre-cached for supplied states, calculating"
        )
        return calculate_uncoupled_hamiltonian_X(QN, nprocs)


def calculate_uncoupled_hamiltonian_X(QN, nprocs=1):
    """Calculate the uncoupled X state hamiltonian for the supplies set of 
    basis states.
    Calculated directly from supplied basis states

    Args:
        QN (array): array of UncoupledBasisStates

    Returns:
        dict: dictionary with all X state hamiltonian terms
    """
    for qn in QN:
        assert qn.isUncoupled, "supply list with UncoupledBasisStates"

    return {
        "Hff": HMatElems(Hff_X_alt, QN, nprocs=nprocs),
        "HSx": HMatElems(HSx, QN, nprocs=nprocs),
        "HSy": HMatElems(HSy, QN, nprocs=nprocs),
        "HSz": HMatElems(HSz, QN, nprocs=nprocs),
        "HZx": HMatElems(HZx_X, QN, nprocs=nprocs),
        "HZy": HMatElems(HZy_X, QN, nprocs=nprocs),
        "HZz": HMatElems(HZz_X, QN, nprocs=nprocs),
    }


def generate_coupled_hamiltonian_B(QN, nprocs=1):
    """Calculate the coupled B state hamiltonian for the supplied set of 
    basis states.
    Retrieved from a pre-calculated sqlite3 database

    Args:
        QN (array): array of UncoupledBasisStates

    Returns:
        dict: dictionary with all B state hamiltonian terms
    """
    for qn in QN:
        assert qn.isCoupled, "supply list withCoupledBasisStates"

    pre_cached = check_states_coupled_hamiltonian_B(QN)

    if pre_cached:
        path = Path(__file__).parent.parent / "pre_calculated"
        db = path / "coupled_hamiltonian_B.db"

        return retrieve_coupled_hamiltonian_B_sqlite(QN, db)
    else:
        logging.warning(
            "B state Hamiltonian not pre-cached for supplied states, calculating"
        )
        return calculate_coupled_hamiltonian_B(QN, nprocs)


def calculate_coupled_hamiltonian_B(QN, nprocs=1):
    for qn in QN:
        assert qn.isCoupled, "supply list withCoupledBasisStates"
    return {
        "Hrot": HMatElems(Hrot_B, QN, nprocs=nprocs),
        "H_mhf_Tl": HMatElems(H_mhf_Tl, QN, nprocs=nprocs),
        "H_mhf_F": HMatElems(H_mhf_F, QN, nprocs=nprocs),
        "H_LD": HMatElems(H_LD, QN, nprocs=nprocs),
        "H_cp1_Tl": HMatElems(H_cp1_Tl, QN, nprocs=nprocs),
        "H_c_Tl": HMatElems(H_c_Tl, QN, nprocs=nprocs),
        "HZz": HMatElems(HZz_B, QN, nprocs=nprocs),
    }
