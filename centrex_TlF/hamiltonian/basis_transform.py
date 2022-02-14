import multiprocessing
from pathlib import Path

import numpy as np
from centrex_TlF.hamiltonian.utils_multiprocessing import multi_transformation_matrix
from centrex_TlF.hamiltonian.utils_sqlite import (
    retrieve_S_transform_uncoupled_to_coupled_sqlite,
)
from tqdm import tqdm

__all__ = ["generate_transform_matrix", "calculate_transform_matrix"]


def generate_transform_matrix(basis1, basis2, progress=False):
    """
    Function that generates a transform matrix that takes Hamiltonian expressed
    in basis1 to basis2: H_2 = S.conj().T @ H_1 @ S
    Retrieved from a pre-calculated sqlite3 database.

    inputs:
    basis1 = list of basis states that defines basis1
    basis2 = list of basis states that defines basis2
    progress = boolean to display tqdm progress bar

    returns:
    S = transformation matrix that takes Hamiltonian (or any operator) from
    basis1 to basis2
    """
    path = Path(__file__).parent.parent / "pre_calculated"
    db = path / "transformation.db"

    # Check that the two bases have the same dimension
    assert len(basis1) == len(basis2), "Bases don't have the same dimension"

    return retrieve_S_transform_uncoupled_to_coupled_sqlite(basis1, basis2, db)


def calculate_transform_matrix(basis1, basis2, progress=False, nprocs=2):
    """
    Function that generates a transform matrix that takes Hamiltonian expressed
    in basis1 to basis2: H_2 = S.conj().T @ H_1 @ S
    Calculated directly from supplies basis sets.

    inputs:
    basis1 = list of basis states that defines basis1
    basis2 = list of basis states that defines basis2
    progress = boolean to display tqdm progress bar
    nprocs = number of processes to utilize with multiprocessing

    returns:
    S = transformation matrix that takes Hamiltonian (or any operator) from
    basis1 to basis2
    """

    # Check that the two bases have the same dimension
    assert len(basis1) == len(basis2), "Bases don't have the same dimension"

    # multiprocessing
    if nprocs > 1:
        with multiprocessing.Pool(nprocs) as pool:
            result = pool.starmap(
                multi_transformation_matrix,
                [(i, state1, basis2) for i, state1 in enumerate(basis1)],
            )
        S = np.sum(result, axis=0)

    else:
        # Initialize S
        S = np.zeros((len(basis1), len(basis1)), dtype=complex)

        # Loop over the bases and calculate inner products
        for i, state1 in enumerate(tqdm(basis1, disable=not progress)):
            for j, state2 in enumerate(basis2):
                S[i, j] = state1 @ state2

    # Return the transform matrix
    return S
