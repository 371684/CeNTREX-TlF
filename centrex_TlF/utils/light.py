from typing import Optional, List

import numpy as np
import scipy.constants as cst

__all__ = [
    "generate_1D_multipass",
    "generate_2D_multipass",
    "calculate_power_from_rabi_gaussian_beam",
    "calculate_rabi_from_power_gaussian_beam",
    "calculate_power_from_rabi_gaussian_beam_microwave",
    "calculate_rabi_from_power_gaussian_beam_microwave",
    "calculate_intensity_from_power_gaussian_beam",
    "intensity_to_electric_field",
]


def gaussian(x: float, μ: float, σ: float) -> float:
    """
    Non-normalized gaussian function

    Args:
        x (float): x values to evaluate gaussian
        μ (float): mean of gaussian
        σ (float): standard deviation of gaussian

    Returns:
        (float): gaussian evaluated at x
    """
    return np.exp(-((x - μ) ** 2) / (2 * σ ** 2))


def gaussian_2d(
    y: float, z: float, μy: float, μz: float, σy: float, σz: float
) -> float:
    """
    Non-normalized 2D gaussian function

    Args:
        y (float): y values to evaluate gaussian
        z (float): z values to evaluate gaussian
        μy (float): mean of gaussian in y
        μz (float): mean of gaussian in z
        σy (float): standard deviation of gaussian in y
        σz (float): standard deviation of gaussian in y

    Returns:
        (float): gaussian evaluated at y and z
    """
    a = (y - μy) ** 2 / (2 * σy ** 2)
    b = (z - μz) ** 2 / (2 * σz ** 2)
    return np.exp(-(a + b))


def gaussian_amp(x: float, a: float, μ: float, σ: float) -> float:
    """
    Non-normalized gaussian function with amplitude a

    Args:
        x (float): x values to evaluate gaussian
        a (float): amplitude
        μ (float): mean of gaussian
        σ (float): standard deviation of gaussian

    Returns:
        (float): gaussian evaluated at x
    """
    return a * np.exp(-((x - μ) ** 2) / (2 * σ ** 2))


def gaussian_2d_amp(
    a: float, y: float, z: float, μy: float, μz: float, σy: float, σz: float
) -> float:
    """
    Non-normalized 2D gaussian function with amplitude a

    Args:
        a (float): amplitude of gaussian
        y (float): y values to evaluate gaussian
        z (float): z values to evaluate gaussian
        μy (float): mean of gaussian in y
        μz (float): mean of gaussian in z
        σy (float): standard deviation of gaussian in y
        σz (float): standard deviation of gaussian in y

    Returns:
        (float): gaussian evaluated at y and z
    """
    return a * gaussian_2d(y, z, μy, μz, σy, σz)


def multipass_prism_order(passes: int) -> List:
    """
    Generate the multipass prism pass order.
    Passes in list are ordered by geometry (left-right or right-left),
    number indicate the pass number.

    Args:
        passes (int): total number of passes

    Returns:
        (list): entries indicate pass number
    """
    npass = [1]
    for p in range(1, passes):
        if p % 2 == 0:
            npass.append(p + 1)
        else:
            npass.append(passes - p)
    return npass


def generate_1D_multipass(
    x: np.ndarray, npasses: int, loss: float, σ: float, spacing: float
) -> np.ndarray:
    """
    Generate a 1D multipass for coordinates x

    Args:
        x (np.ndarray): x coordinates
        npasses (int): number of passes
        loss (float): loss per pass
        σ (float): laser 1-sigma width
        spacing (float): spacing between passes

    Returns:
        (np.ndarray): 1D multipass for coordinates x
    """
    reflections = np.array(multipass_prism_order(npasses)) - 1
    amplitudes = np.array([(1 - loss) ** r for r in reflections])
    beam_locs = np.array([i * spacing for i in range(npasses)])
    multipass = gaussian_amp(
        x[:, np.newaxis], amplitudes[np.newaxis, :], beam_locs[np.newaxis, :], σ
    ).sum(axis=1)
    return multipass


def generate_2D_multipass(
    X: np.ndarray,
    Y: np.ndarray,
    npasses: int,
    loss: float,
    σx: float,
    σy: float,
    spacing: float,
) -> np.ndarray:
    """
    Generate a 2D multipass for coordinates X and Y.

    Args:
        X (np.ndarray): X coordinates
        Y (np.ndarray): Y coordinates
        npasses (int): number of passes
        loss (float): loss per pass
        σx (float): 1-sigma laser width in x
        σy (float): 1-sigma laser width in y
        spacing (float): spacing between passes
    """
    reflections = np.array(multipass_prism_order(npasses)) - 1
    amplitudes = np.array([(1 - loss) ** r for r in reflections])
    beam_locs = np.array([i * spacing for i in range(npasses)])
    multipass = gaussian_2d_amp(
        amplitudes[np.newaxis, :],
        X[:, :, np.newaxis],
        Y[:, :, np.newaxis],
        0,
        beam_locs[np.newaxis, :],
        σx,
        σy,
    ).sum(axis=2)
    return multipass


def calculate_intensity_from_power_gaussian_beam(power: float, σx: float, σy: float):
    """Calculate the maximum laser intensity of a gaussian beam from the total
    laser power given the beam parameters σx and σy

    Args:
        power (float): beam power [W]
        σx (float): beam standard deviation in x [m]
        σy (float): beam standard deviation in y [m]

    Returns:
        float: intensity in W/m^2
    """
    return power / (2 * np.pi * σx * σy)


def calculate_power_from_rabi_gaussian_beam(
    Ω: float, main_coupling: float, σx: float, σy: float, D: float = 2.6675506e-30
) -> float:
    """calculate the required power for a given Ω, given a main transition
    matrix element

    Args:
        Ω (float): Rabi rate [2π⋅Hz]
        main_coupling (complex): main transition matrix element
        σx (float): σx of beam
        σy (float): σy of beam
        D (float, optional): effective dipole moment for transitions.
                            Defaults to 2.6675506e-30 for optical TlF.

    Returns:
        float: power [W]
    """
    # Electric field
    E = Ω * cst.hbar / (main_coupling * D)

    # convert to peak intensity
    intensity = 1 / 2 * cst.c * cst.epsilon_0 * E ** 2

    # convert power to amplitude of the gaussian
    P = intensity * (2 * np.pi * σx * σy)

    return P


def intensity_to_electric_field(intensity: float) -> float:
    """Convert intensity in W/m^2 to the electric field

    Args:
        intensity (float): intensity in W/m^2

    Returns:
        float: electric field E
    """
    return np.sqrt((2 / (cst.c * cst.epsilon_0)) * intensity)


def calculate_rabi_from_power_gaussian_beam(
    P: float,
    main_coupling: float,
    σx: float,
    σy: float,
    D: float = 2.6675506e-30,
    Γ: Optional[float] = None,
) -> float:
    """calculate Ω for a given power, given a main transition matrix element

    Args:
        P (float): optical power [W]
        main_coupling (complex): main transition matrix element
        σx (float): σx of beam
        σy (float): σy of beam
        D ([type], optional): effective dipole moment for transitions.
                                Defaults to 2.6675506e-30 for optical TlF.
        Γ (float, optional): Γ to normalize Rabi rate with
    Returns:
        float: power [W]
    """
    # intensity from power
    intensity = P / (2 * np.pi * σx * σy)

    # electric field from intensity
    E = np.sqrt(intensity * 2 / (cst.c * cst.epsilon_0))

    # rabi rate from electric field
    Ω = (E * main_coupling * D) / cst.hbar

    # normalize with Γ if provided
    if Γ:
        Ω /= Γ

    return Ω


def calculate_power_from_rabi_gaussian_beam_microwave(
    Ω: float, main_coupling: float, σx: float, σy: float, D: float = 1.4103753e-29
) -> float:
    """Calculate the microwave Ω for a given power, given a main transition matrix element

    Args:
        P (float): optical power [W]
        main_coupling (complex): main transition matrix element
        σx (float): σx of beam
        σy (float): σy of beam
        D ([type], optional): effective dipole moment for transitions.
                                Defaults to 1.4103753e-29 for microwave TlF.
        Γ (float, optional): Γ to normalize Rabi rate with
    Returns:
        float: power [W]
    """
    return calculate_power_from_rabi_gaussian_beam(Ω, main_coupling, σx, σy, D)


def calculate_rabi_from_power_gaussian_beam_microwave(
    P: float,
    main_coupling: float,
    σx: float,
    σy: float,
    D: float = 1.4103753e-29,
    Γ: Optional[float] = None,
) -> float:
    """Calculate the microwave Ω for a given power, given a main transition matrix element

    Args:
        P (float): optical power [W]
        main_coupling (complex): main transition matrix element
        σx (float): σx of beam
        σy (float): σy of beam
        D ([type], optional): effective dipole moment for transitions.
                                Defaults to 1.4103753e-29 for microwave TlF.
        Γ (float, optional): Γ to normalize Rabi rate with
    Returns:
        float: power [W]
    """
    return calculate_rabi_from_power_gaussian_beam(P, main_coupling, σx, σy, D, Γ)
