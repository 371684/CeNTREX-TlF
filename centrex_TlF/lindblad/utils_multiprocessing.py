def multi_C_ρ_Cconj(C, Cᶜ, ρ):
    return C @ ρ @ Cᶜ


def multi_system_of_equations_to_lines(system, ρ, idx):
    n_states = system.shape[0]
    code_lines = []
    for idy in range(n_states):
        if system[idx, idy] != 0:
            if idy >= idx:
                cline = str(system[idx, idy])
                cline = f"du[{idx+1},{idy+1}] = " + cline
                cline = cline.replace("(t)", "")
                cline = cline.replace("(t)", "")
                cline = cline.replace("I", "1im")
                cline += "\n"
                for i in range(system.shape[0]):
                    for j in range(system.shape[1]):
                        _ = str(ρ[i, j])
                        cline = cline.replace(_ + "*", f"ρ[{i+1},{j+1}]*")
                        cline = cline.replace(_ + " ", f"ρ[{i+1},{j+1}] ")
                        cline = cline.replace(_ + "\n", f"ρ[{i+1},{j+1}]")
                        cline = cline.replace(_ + ")", f"ρ[{i+1},{j+1}])")
                cline = cline.strip()
                code_lines.append(cline)
            else:
                if system[idx, idy] != 0:
                    cline = f"du[{idx+1},{idy+1}] = conj(du[{idy+1},{idx+1}])"
                    code_lines.append(cline)
    return code_lines
