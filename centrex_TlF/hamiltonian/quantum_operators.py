import numpy as np
from centrex_TlF.states.states import State, UncoupledBasisState

########################################################
### Diagonal operators multiple state by eigenvalue
########################################################


def J2(psi):
    return State([(psi.J * (psi.J + 1), psi)])


def J4(psi):
    return State([((psi.J * (psi.J + 1)) ** 2, psi)])


def J6(psi):
    return State([((psi.J * (psi.J + 1)) ** 3, psi)])


def Jz(psi):
    return State([(psi.mJ, psi)])


def I1z(psi):
    return State([(psi.m1, psi)])


def I2z(psi):
    return State([(psi.m2, psi)])


########################################################
###
########################################################


def Jp(psi):
    amp = np.sqrt((psi.J - psi.mJ) * (psi.J + psi.mJ + 1))
    ket = UncoupledBasisState(
        psi.J,
        psi.mJ + 1,
        psi.I1,
        psi.m1,
        psi.I2,
        psi.m2,
        Omega=psi.Omega,
        P=psi.P,
        electronic_state=psi.electronic_state,
    )
    return State([(amp, ket)])


def Jm(psi):
    amp = np.sqrt((psi.J + psi.mJ) * (psi.J - psi.mJ + 1))
    ket = UncoupledBasisState(
        psi.J,
        psi.mJ - 1,
        psi.I1,
        psi.m1,
        psi.I2,
        psi.m2,
        Omega=psi.Omega,
        P=psi.P,
        electronic_state=psi.electronic_state,
    )
    return State([(amp, ket)])


def I1p(psi):
    amp = np.sqrt((psi.I1 - psi.m1) * (psi.I1 + psi.m1 + 1))
    ket = UncoupledBasisState(
        psi.J,
        psi.mJ,
        psi.I1,
        psi.m1 + 1,
        psi.I2,
        psi.m2,
        Omega=psi.Omega,
        P=psi.P,
        electronic_state=psi.electronic_state,
    )
    return State([(amp, ket)])


def I1m(psi):
    amp = np.sqrt((psi.I1 + psi.m1) * (psi.I1 - psi.m1 + 1))
    ket = UncoupledBasisState(
        psi.J,
        psi.mJ,
        psi.I1,
        psi.m1 - 1,
        psi.I2,
        psi.m2,
        Omega=psi.Omega,
        P=psi.P,
        electronic_state=psi.electronic_state,
    )
    return State([(amp, ket)])


def I2p(psi):
    amp = np.sqrt((psi.I2 - psi.m2) * (psi.I2 + psi.m2 + 1))
    ket = UncoupledBasisState(
        psi.J,
        psi.mJ,
        psi.I1,
        psi.m1,
        psi.I2,
        psi.m2 + 1,
        Omega=psi.Omega,
        P=psi.P,
        electronic_state=psi.electronic_state,
    )
    return State([(amp, ket)])


def I2m(psi):
    amp = np.sqrt((psi.I2 + psi.m2) * (psi.I2 - psi.m2 + 1))
    ket = UncoupledBasisState(
        psi.J,
        psi.mJ,
        psi.I1,
        psi.m1,
        psi.I2,
        psi.m2 - 1,
        Omega=psi.Omega,
        P=psi.P,
        electronic_state=psi.electronic_state,
    )
    return State([(amp, ket)])


########################################################
###
########################################################


def Jx(psi):
    return 0.5 * (Jp(psi) + Jm(psi))


def Jy(psi):
    return -0.5j * (Jp(psi) - Jm(psi))


def I1x(psi):
    return 0.5 * (I1p(psi) + I1m(psi))


def I1y(psi):
    return -0.5j * (I1p(psi) - I1m(psi))


def I2x(psi):
    return 0.5 * (I2p(psi) + I2m(psi))


def I2y(psi):
    return -0.5j * (I2p(psi) - I2m(psi))


########################################################
### Composition of operators
########################################################


def com(A, B, psi):
    ABpsi = State()
    # operate with A on all components in B|psi>
    for amp, cpt in B(psi):
        ABpsi += amp * A(cpt)
    return ABpsi
