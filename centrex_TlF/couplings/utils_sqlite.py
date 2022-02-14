import json
from pathlib import Path

import numpy as np


def retrieve_ED_ME_coupled_sqlite_single_rme(a, b, pol_vec, con):
    table = "ED_ME_coupled_rme"
    with con:
        cur = con.cursor()
        string = (
            f"J₁ = {a.J} AND F1₁ = {a.F1} AND F₁ = {int(a.F)} AND mF₁ = {int(a.mF)} "
            f"AND I1₁ = {a.I1} AND I2₁ = {a.I2} AND Ω₁ = {int(a.Omega)} "
            f"AND state₁ = '{a.electronic_state}' "
        )
        string += (
            f"AND J₂ = {b.J} AND F1₂ = {b.F1} AND F₂ = {int(b.F)} AND mF₂ = {int(b.mF)}"
            f" AND I1₂ = {b.I1} AND I2₂ = {b.I2} AND Ω₂ = {int(b.Omega)} "
            f"AND state₂ = '{b.electronic_state}' "
        )
        cur.execute(f"select value_real, value_imag from {table} WHERE {string}")
        values = cur.fetchall()
        if values:
            values = values[0]
            values = values[0] + 1j * values[1]
        else:
            values = complex(0)
    return values


def retrieve_ED_ME_coupled_sqlite_single(a, b, pol_vec, con):
    table = "ED_ME_coupled"
    with con:
        cur = con.cursor()
        string = (
            f"J₁ = {a.J} AND F1₁ = {a.F1} AND F₁ = {int(a.F)} AND mF₁ = {int(a.mF)} "
            f"AND I1₁ = {a.I1} AND I2₁ = {a.I2} AND Ω₁ = {int(a.Omega)} "
            f"AND state₁ = '{a.electronic_state}' "
        )
        string += (
            f"AND J₂ = {b.J} AND F1₂ = {b.F1} AND F₂ = {int(b.F)} AND mF₂ = {int(b.mF)}"
            f" AND I1₂ = {b.I1} AND I2₂ = {b.I2} AND Ω₂ = {int(b.Omega)} "
            f"AND state₂ = '{b.electronic_state}' "
        )
        string += (
            f"AND Px = {int(pol_vec[0])} AND Py = {int(pol_vec[1])} "
            f"AND Pz = {int(pol_vec[2])}"
        )
        cur.execute(f"select value_real, value_imag from {table} WHERE {string}")
        values = cur.fetchall()
        if values:
            values = values[0]
            values = values[0] + 1j * values[1]
        else:
            values = complex(0)
    return values


def check_states_in_ED_ME_coupled(Jg, Je, pol_vec):
    # load json
    path = Path(__file__).parent.parent / "pre_calculated"
    js = path / "precalculated.json"
    with open(js) as json_file:
        f = json.load(json_file)

    # check if ground state J is pre-cached
    if not np.all(J in f["matrix_elements"]["X"] for J in Jg):
        return False
    # check if excited state J is pre-cached
    if not np.all(J in f["matrix_elements"]["B"] for J in Je):
        return False
    # check if the pol-vec is pre-cached
    if not list(pol_vec) in f["matrix_elements"]["pol_vec"]:
        return False
    else:
        return True
