import numpy as np

__all__ = [
    "generate_random_coords_circle_2D",
    "generate_random_longitudinal_velocities",
    "generate_random_transverse_velocities",
]


def generate_random_coords_circle_2D(n, σ):
    """generate random xy coordinates in a circle

    Args:
        n (int): number of coordinates to generate
        σ (float): σ width of circle to generate coordinates in

    Returns:
        np.ndarray: 2D array with xy coordinates
    """

    def generate_random_r(n, σ):
        return np.random.normal(0, σ, n)

    xy = np.empty([2, n])
    xy[0, :] = generate_random_r(n, σ)
    xy[1, :] = generate_random_r(n, σ)
    return xy


def generate_random_longitudinal_velocities(n, vμ, vσ):
    """generate random longitudinal velocities

    Args:
        n (int): number of velocities to genrate
        vμ (float): mean velocity
        vσ (float): σ velocity width

    Returns:
        np.ndarray: 1D array with longitudinal velocities
    """
    return np.random.normal(vμ, vσ, n)


def generate_random_transverse_velocities(n, vσ):
    """generate random transverse velocities

    Args:
        n (int): number of velocities to generate
        vσ (float): σ width of transverse velocity distribution

    Returns:
        np.ndarray: 2D array with transverse velocities
    """
    velocities = np.random.multivariate_normal(
        [0, 0], np.array([[vσ ** 2, 0], [0, vσ ** 2]]), n
    ).T
    return velocities
