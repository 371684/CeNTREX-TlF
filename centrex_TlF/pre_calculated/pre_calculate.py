import json
import copy
import pickle
import sqlite3
import numpy as np
from tqdm import tqdm
from pathlib import Path
import centrex_TlF as centrex

def create_coupled_hamiltonian_B_sqlite(con):
    cur = con.cursor()
    try:
        for table in ['Hrot', 'H_mhf_Tl', 'H_mhf_F', 'H_LD', 'H_cp1_Tl', 'H_c_Tl', 'HZz']:
            table = table.strip()
            string = f'''CREATE TABLE {table} (J₁ int, F1₁ real, F₁ int, mF₁ int, I1₁ real, I2₁ real, P₁ int, J₂ int, F1₂ real, F₂ int, mF₂ int, I1₂ real, I2₂ real, P₂ int, value_real real, value_imag real, 
                            unique (J₁, F1₁, F₁, mF₁, I1₁, I2₁, P₁, J₂, F1₂, F₂, mF₂, I1₂, I2₂, P₂))'''
            cur.execute(string)
            con.commit()
    except Exception as e:
        pass
    return

def generate_coupled_hamiltonian_B_sqlite(QN, con, progress = False):
    with con:
        cur = con.cursor()
        for table in ['Hrot', 'H_mhf_Tl', 'H_mhf_F', 'H_LD', 'H_cp1_Tl', 'H_c_Tl', 'HZz']:
            tab = copy.copy(table)
            if tab in ['Hrot', 'HZz']:
                tab += '_B'
            desc = f"pre-calculating coupled B state Hamiltonian; {table}"
            for i,a in tqdm(enumerate(QN), total = len(QN), disable = not progress, desc = desc):
                for j in range(i,len(QN)):
                    b = QN[j]
                    val = (1*a)@(eval(f"centrex.hamiltonian.hamiltonian_B_terms_coupled.{tab}")(b))
                    if val != 0:
                        try:
                            string = f"{a.J}, {a.F1}, {a.F}, {a.mF}, {a.I1}, {a.I2}, {a.P}, {b.J}, {b.F1}, {b.F}, {b.mF}, {b.I1}, {b.I2}, {b.P}, {val.real}, {val.imag}"
                            cmd = f"INSERT INTO {table} VALUES ({string})"
                            cur.execute(cmd)
                        except Exception as e:
                            raise e
            con.commit()
    return

def create_uncoupled_hamiltonian_X_sqlite(con):
    cur = con.cursor()
    try:
        for table in 'Hff, HSx, HSy, HSz, HZx, HZy, HZz'.split(","):
            table = table.strip()
            cur.execute(f'''CREATE TABLE {table} (J₁ int, mJ₁ int, I1₁ real, m1₁ real, I2₁ real, m2₁ real, J₂ int, mJ₂ int, I1₂ real, m1₂ real, I2₂ real, m2₂ real, value_real real, value_imag real, 
                            unique (J₁, mJ₁, I1₁, m1₁, I2₁, m2₁, J₂, mJ₂, I1₂, m1₂, I2₂, m2₂))''')
            con.commit()
    except:
        pass
    return

def generate_uncoupled_hamiltonian_X_sqlite(QN, con, progress = False):
    cur = con.cursor()
    for table in 'Hff_alt, HSx, HSy, HSz, HZx, HZy, HZz'.split(","):
        table = table.strip()
        if table =='Hff_alt':
            tab = 'Hff_X_alt'
            table = 'Hff'
        elif 'HZ' in table:
            tab = table + '_X'
        else:
            tab = table
        desc = f"pre-calculating uncoupled X state Hamiltonian; {table}"
        for a in tqdm(QN, desc = desc):
            for b in QN:
                val = (1*a)@(eval(f"centrex.hamiltonian.hamiltonian_terms_uncoupled.{tab}")(b))
                if val != 0:
                    try:
                        string = f"{a.J}, {a.mJ}, {a.I1}, {a.m1}, {a.I2}, {a.m2}, {b.J}, {b.mJ}, {b.I1}, {b.m1}, {b.I2}, {b.m2}, {val.real}, {val.imag}"
                        cmd = f"INSERT INTO {table} VALUES ({string})"
                        cur.execute(cmd)
                    except Exception as e:
                        pass
        con.commit()

def create_transformation_uncoupled_to_coupled(con):
    cur = con.cursor()
    try:
        table = 'uncoupled_to_coupled'
        cur.execute(f'''CREATE TABLE {table} (J int, mJ int, I1 real, m1 real, I2 real, m2 real, Jc int, F1 real, F int, mF int, I1c real, I2c real, value_real real, value_imag real,
                        unique (J, mJ, I1, m1, I2, m2, Jc, F1, F, mF, I1c, I2c))''')
    except Exception as e:
        raise e
    con.commit()
    return

def generate_transformation_uncoupled_to_coupled(QN, QNc, con, progress = False):
    cur = con.cursor()
    desc = "pre-calculating transformation uncoupled to coupled"
    with con:
        for a in tqdm(QN, disable = not progress, desc = desc):
            for b in QNc:
                val = a@b
                if val != 0:
                    try:
                        string = f"{a.J}, {a.mJ}, {a.I1}, {a.m1}, {a.I2}, {a.m2}, {b.J}, {b.F1}, {b.F}, {b.mF}, {b.I1}, {b.I2}, {val.real}, {val.imag}"
                        cmd = f"INSERT INTO uncoupled_to_coupled VALUES ({string})"
                        cur.execute(cmd)
                    except Exception as e:
#                         raise e
                        pass
        con.commit()
    return

def create_ED_ME(con):
    cur = con.cursor()
    try:
        for table in ['ED_ME_coupled']:
            table = table.strip()
            string = f'''CREATE TABLE {table} (J₁ int, F1₁ real, F₁ int, mF₁ int, I1₁ real, I2₁ real, Ω₁ int, state₁ text, J₂ int, F1₂ real, F₂ int, mF₂ int, I1₂ real, I2₂ real, Ω₂ int, state₂ text, Px int, Py int, Pz int, value_real real, value_imag real, 
                            unique (J₁, F1₁, F₁, mF₁, I1₁, I2₁, Ω₁, state₁, J₂, F1₂, F₂, mF₂, I1₂, I2₂, Ω₂, state₂, Px, Py, Pz))'''
            cur.execute(string)
            con.commit()
    except Exception as e:
        pass

    try:
        for table in ['ED_ME_coupled_rme']:
            table = table.strip()
            string = f'''CREATE TABLE {table} (J₁ int, F1₁ real, F₁ int, mF₁ int, I1₁ real, I2₁ real, Ω₁ int, state₁ text, J₂ int, F1₂ real, F₂ int, mF₂ int, I1₂ real, I2₂ real, Ω₂ int, state₂ text, value_real real, value_imag real, 
                            unique (J₁, F1₁, F₁, mF₁, I1₁, I2₁, Ω₁, state₁, J₂, F1₂, F₂, mF₂, I1₂, I2₂, Ω₂, state₂))'''
            cur.execute(string)
            con.commit()
    except Exception as e:
        pass
    return

def generate_ED_ME(QN, pol_vecs, con, progress = False):
    with con:
        for pol_vec in pol_vecs:
            cur = con.cursor()
            for table in ['ED_ME_coupled']:
                desc = f"coupling elements: pol = {pol_vec}, reduced = False"
                for i,a in tqdm(enumerate(QN), total = len(QN), disable = not progress, desc = desc):
                    for j,b in enumerate(QN):
                        val = centrex.couplings.matrix_elements.ED_ME_coupled(a,b,pol_vec, rme_only = False)
                        if val != 0:
                            try:
                                string = f"{a.J}, {a.F1}, {a.F}, {a.mF}, {a.I1}, {a.I2}, {a.Omega}, '{a.electronic_state}', {b.J}, {b.F1}, {b.F}, {b.mF}, {b.I1}, {b.I2}, {b.Omega}, '{b.electronic_state}', {int(pol_vec[0])}, {int(pol_vec[1])}, {int(pol_vec[2])}, {val.real}, {val.imag}"
                                cmd = f"INSERT INTO {table} VALUES ({string})"
                                cur.execute(cmd)
                            except Exception as e:
                                pass
                con.commit()
    with con:
        cur = con.cursor()
        desc = 'coupling elements: reduced = True'
        for table in ['ED_ME_coupled_rme']:
            for i,a in tqdm(enumerate(QN), total = len(QN), desc = desc, disable = not progress):
                for j,b in enumerate(QN):
                    val = centrex.couplings.matrix_elements.ED_ME_coupled(a,b, rme_only=True)
                    if val != 0:
                        try:
                            string = f"{a.J}, {a.F1}, {a.F}, {a.mF}, {a.I1}, {a.I2}, {a.Omega}, '{a.electronic_state}', {b.J}, {b.F1}, {b.F}, {b.mF}, {b.I1}, {b.I2}, {b.Omega}, '{b.electronic_state}', {val.real}, {val.imag}"
                            cmd = f"INSERT INTO {table} VALUES ({string})"
                            cur.execute(cmd)
                        except Exception as e:
                            pass
            con.commit()
    return

def generate_transitions(fname, QN_X, QN_B, E, B, nprocs):
    QN, H_tot = centrex.transitions.calculate_energies(QN_X, QN_B, E, B, 
                                nprocs = nprocs)

    QN = [s.remove_small_components(1e-3) for s in QN]

    with open(fname, 'wb') as f:
        pickle.dump({'QN': QN, 'H': H_tot}, f)

    return


def generate_pre_calculated(nprocs, db_path = None):
    path = Path(__file__).parent.absolute()
    js = path / "precalculated.json"
    with open(js) as json_file:
        config = json.load(json_file)
    
    if not db_path:
        db_path = path

    # db = 'uncoupled_hamiltonian_X'
    # Js = config[db]
    # print(f"pre-calculating {db} for J = {Js}")
    # try:
    #     (db_path / (db + '.db')).unlink()
    # except:
    #     pass
    # con = sqlite3.connect(db_path / (db + '.db'))
    # create_uncoupled_hamiltonian_X_sqlite(con)
    # QN = centrex.states.generate_uncoupled_states_ground(Js)
    # generate_uncoupled_hamiltonian_X_sqlite(QN, con, progress = True)

    # db = 'coupled_hamiltonian_B'
    # Js = config[db]
    # print(f"pre-calculating {db} for J = {Js}")
    # try:
    #     (db_path / (db + '.db')).unlink()
    # except:
    #     pass
    # con = sqlite3.connect(db_path / (db + '.db'))
    # create_coupled_hamiltonian_B_sqlite(con)
    # QN = centrex.states.generate_coupled_states_excited(Js, Ps = [-1,1])
    # generate_coupled_hamiltonian_B_sqlite(QN, con, progress = True)

    # db = 'matrix_elements'
    # Jg = config[db]['X']
    # Je = config[db]['B']
    # QN = list(np.append(
    #     centrex.states.generate_coupled_states_ground(Jg),
    #     centrex.states.generate_coupled_states_excited(Je, Ps = [+1])
    # ))
    # QN_B = centrex.states.generate_coupled_states_excited(Je, Ps = [+1])
    # for idx in range(len(QN_B)):
    #     QN_B[idx].Omega *= -1
    # QN.extend(QN_B)
    # pol_vecs = config[db]['pol_vec']
    # print(f"pre-calculating {db} for Jg = {Jg} and Je = {Je}")
    # try:
    #     (db_path / (db + '.db')).unlink()
    # except:
    #     pass
    # con = sqlite3.connect(db_path / (db + '.db'))
    # create_ED_ME(con)
    # generate_ED_ME(QN, pol_vecs, con, progress = True)

    # db = 'transformation'
    # Js = config[db]
    # print(f"pre-calculating {db} for J = {Js}")
    # try:
    #     (db_path / (db + '.db')).unlink()
    # except:
    #     pass
    # con = sqlite3.connect(db_path / (db + '.db'))
    # create_transformation_uncoupled_to_coupled(con)
    # QN = centrex.states.generate_uncoupled_states_ground(Js)
    # QNc = centrex.states.generate_coupled_states_ground(Js)
    # generate_transformation_uncoupled_to_coupled(QN, QNc, con, progress = True)

    db = 'transitions'
    JX = config[db]['X']
    JB = config[db]['B']
    E,B = config[db]['field'][0]
    QN_X = centrex.states.generate_coupled_states_ground(JX)
    QN_B = centrex.states.generate_coupled_states_excited(JB, Ps = [-1,+1])
    generate_transitions(db_path / (db + '.pickle'), QN_X, QN_B, E, B, 
                        nprocs = nprocs)
    
    # con.close()


if __name__ == "__main__":
    generate_pre_calculated(nprocs = 1)