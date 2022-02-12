import numpy as np
import scipy.constants as cst

__all__ = ["doppler_shift", "velocity_to_detuning"]


def doppler_shift(velocity, frequency):
    """calculate the doppler shifted freuqency

    Args:
        velocity (float): velocity [m/s]
        frequency (float): frequency [Hz]

    Returns:
        float: doppler shifted frequency in Hz
    """
    return frequency * (1 + velocity / cst.c)


def velocity_to_detuning(velocity, frequency, Γ):
    """convert velocity to detuning in units of Γ

    Args:
        velocity (float): velocity [m/s]
        frequency (float): frequency [Hz]
        Γ (float): Γ [2π⋅Hz]

    Returns:
        float: detuning in units of Γ
    """
    return (doppler_shift(velocity, frequency) - frequency) * 2 * np.pi / Γ
