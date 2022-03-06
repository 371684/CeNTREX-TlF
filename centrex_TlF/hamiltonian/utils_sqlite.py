import json
import sqlite3
import numpy as np
from pathlib import Path


def retrieve_uncoupled_hamiltonian_X_sqlite(QN, db):
    con = sqlite3.connect(db)
    H = {}
    with con:
        cur = con.cursor()
        for term in ["Hff", "HSx", "HSy", "HSz", "HZx", "HZy", "HZz"]:
            result = np.zeros((len(QN), len(QN)), complex)
            for i, a in enumerate(QN):
                for j in range(i, len(QN)):
                    b = QN[j]
                    string = (
                        f"J₁ = {a.J} AND mJ₁ = {a.mJ} AND I1₁ = {a.I1} AND m1₁ = {a.m1}"
                        f" AND I2₁ = {a.I2} AND m2₁ = {a.m2} AND J₂ = {b.J} "
                        f"AND mJ₂ = {b.mJ} AND I1₂ = {b.I1} AND m1₂ = {b.m1} "
                        f"AND I2₂ = {b.I2} AND m2₂ = {b.m2}"
                    )
                    cur.execute(
                        f"select value_real, value_imag from {term} WHERE {string}"
                    )
                    values = cur.fetchall()
                    if values:
                        values = values[0]
                        result[i, j] = values[0] + 1j * values[1]
                        if i != j:
                            result[j, i] = np.conjugate(values[0] + 1j * values[1])
            H[term] = result.copy()
    con.close()
    return H


def retrieve_coupled_hamiltonian_B_sqlite(QN, db):
    con = sqlite3.connect(db)
    H = {}
    with con:
        cur = con.cursor()
        for term in [
            "Hrot",
            "H_mhf_Tl",
            "H_mhf_F",
            "H_LD",
            "H_cp1_Tl",
            "H_c_Tl",
            "HZz",
        ]:
            result = np.zeros((len(QN), len(QN)), complex)
            for i, a in enumerate(QN):
                for j in range(i, len(QN)):
                    b = QN[j]
                    string = (
                        f"J₁ = {a.J} AND F1₁ = {a.F1} AND F₁ = {a.F} AND mF₁ = {a.mF} "
                        f"AND I1₁ = {a.I1} AND I2₁ = {a.I2} AND P₁ = {a.P} "
                        f"AND J₂ = {b.J} AND F1₂ = {b.F1} AND F₂ = {b.F} "
                        f"AND mF₂ = {b.mF} AND I1₂ = {b.I1} AND I2₂ = {b.I2} "
                        f"AND P₂ = {b.P}"
                    )
                    cur.execute(
                        f"select value_real, value_imag from {term} WHERE {string}"
                    )
                    values = cur.fetchall()
                    if values:
                        values = values[0]
                        result[i, j] = values[0] + 1j * values[1]
                        if i != j:
                            result[j, i] = np.conjugate(values[0] + 1j * values[1])
            H[term] = result.copy()
    con.close()
    return H


def retrieve_S_transform_uncoupled_to_coupled_sqlite(basis1, basis2, db):
    con = sqlite3.connect(db)
    cur = con.cursor()

    S_transform = np.zeros((len(basis1), len(basis2)), dtype=complex)
    with con:
        for i, a in enumerate(basis1):
            for j, b in enumerate(basis2):
                string = (
                    f"J = {a.J} AND mJ = {a.mJ} AND I1 = {a.I1} AND m1 = {a.m1} "
                    f"AND I2 = {a.I2} AND m2 = {a.m2} AND Jc = {b.J} AND F1 = {b.F1} "
                    f"AND F = {b.F} AND mF = {b.mF} AND I1 = {b.I1} AND I2 = {b.I2}"
                )
                cur.execute(
                    "select value_real, value_imag from uncoupled_to_coupled WHERE "
                    f"{string}"
                )
                values = cur.fetchall()
                if values:
                    values = values[0]
                    S_transform[i, j] = values[0] + 1j * values[1]
    con.close()
    return S_transform


def check_states_hamiltonian(QN, ham):
    # load json
    path = Path(__file__).parent.parent / "pre_calculated"
    js = path / "precalculated.json"
    with open(js) as json_file:
        f = json.load(json_file)

    # check if Js are pre-cached
    Js = np.unique([s.J for s in QN])
    if not np.all([J in f[ham] for J in Js]):
        return False
    else:
        return True


def check_states_coupled_hamiltonian_B(QN):
    return check_states_hamiltonian(QN, "coupled_hamiltonian_B")


def check_states_uncoupled_hamiltonian_X(QN):
    return check_states_hamiltonian(QN, "uncoupled_hamiltonian_X")
