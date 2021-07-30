import sqlite3
import numpy as np

def retrieve_uncoupled_hamiltonian_X_sqlite(QN, db):
    con = sqlite3.connect(db)
    H = {}
    with con:
        cur = con.cursor()
        for term in ['Hff', 'HSx', 'HSy', 'HSz', 'HZx', 'HZy', 'HZz']:
            result = np.zeros((len(QN),len(QN)), complex)
            for i,a in enumerate(QN):
                for j in range(i,len(QN)):
                    b = QN[j]
                    string = f"J₁ = {a.J} AND mJ₁ = {a.mJ} AND I1₁ = {a.I1} AND m1₁ = {a.m1} AND I2₁ = {a.I2} AND m2₁ = {a.m2} AND J₂ = {b.J} AND mJ₂ = {b.mJ} AND I1₂ = {b.I1} AND m1₂ = {b.m1} AND I2₂ = {b.I2} AND m2₂ = {b.m2}"
                    cur.execute(f"select value_real, value_imag from {term} WHERE {string}")
                    values = cur.fetchall()
                    if values:
                        values = values[0]
                        result[i,j] = values[0] + 1j*values[1]
                        if i != j:
                            result[j,i] = np.conjugate(values[0] + 1j*values[1])
            H[term] = result.copy()
    con.close()
    return H

def retrieve_coupled_hamiltonian_B_sqlite(QN, db):
    con = sqlite3.connect(db)
    H = {}
    with con:
        cur = con.cursor()
        for term in ['Hrot', 'H_mhf_Tl', 'H_mhf_F', 'H_LD', 'H_cp1_Tl', 'H_c_Tl', 'HZz']:
            result = np.zeros((len(QN),len(QN)), complex)
            for i,a in enumerate(QN):
                for j in range(i,len(QN)):
                    b = QN[j]
                    string = f"J₁ = {a.J} AND F1₁ = {a.F1} AND F₁ = {a.F} AND mF₁ = {a.mF} AND I1₁ = {a.I1} AND I2₁ = {a.I2} AND P₁ = {a.P} AND J₂ = {b.J} AND F1₂ = {b.F1} AND F₂ = {b.F} AND mF₂ = {b.mF} AND I1₂ = {b.I1} AND I2₂ = {b.I2} AND P₂ = {b.P}"
                    cur.execute(f"select value_real, value_imag from {term} WHERE {string}")
                    values = cur.fetchall()
                    if values:
                        values = values[0]
                        result[i,j] = values[0] + 1j*values[1]
                        if i != j:
                            result[j,i] = np.conjugate(values[0] + 1j*values[1])
            H[term] = result.copy()
    con.close()
    return H

def retrieve_S_transform_uncoupled_to_coupled_sqlite(basis1, basis2, db):
    con = sqlite3.connect(db)
    cur = con.cursor()

    S_transform = np.zeros((len(basis1), len(basis2)), dtype = complex)
    with con:
        for i,a in enumerate(basis1):
            for j,b in enumerate(basis2):
                string = f"J = {a.J} AND mJ = {a.mJ} AND I1 = {a.I1} AND m1 = {a.m1} AND I2 = {a.I2} AND m2 = {a.m2} AND Jc = {b.J} AND F1 = {b.F1} AND F = {b.F} AND mF = {b.mF} AND I1 = {b.I1} AND I2 = {b.I2}"
                cur.execute(f"select value_real, value_imag from uncoupled_to_coupled WHERE {string}")
                values = cur.fetchall()
                if values:
                    values = values[0]
                    S_transform[i,j] = values[0] + 1j*values[1]
    con.close()
    return S_transform