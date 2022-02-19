import numpy as np
from centrex_TlF.utils.light import (
    calculate_intensity_from_power_gaussian_beam,
    calculate_power_from_rabi_gaussian_beam,
    calculate_power_from_rabi_gaussian_beam_microwave,
    calculate_rabi_from_power_gaussian_beam,
    calculate_rabi_from_power_gaussian_beam_microwave,
    gaussian,
    gaussian_2d,
    gaussian_2d_amp,
    gaussian_amp,
    generate_1D_multipass,
    generate_2D_multipass,
    intensity_to_electric_field,
    multipass_prism_order,
)


def test_gaussian():
    assert gaussian(0.38, 2.1, 0.81) == 0.10492238740021324


def test_gaussian_2d():
    assert gaussian_2d(0.38, -0.38, 1.2, 2.1, 0.81, 0.78) == 0.0038219624834126514


def test_gaussian_amp():
    assert gaussian_amp(0.38, 2.4, 2.1, 0.81) == 0.25181372976051175


def test_gaussian_2d_amp():
    assert (
        gaussian_2d_amp(2.4, 0.38, -0.38, 1.2, 2.1, 0.81, 0.78) == 0.009172709960190362
    )


def test_multipass_prism_order():
    assert multipass_prism_order(11) == [1, 10, 3, 8, 5, 6, 7, 4, 9, 2, 11]


def test_generate_1D_multipass():
    np.testing.assert_allclose(
        generate_1D_multipass(np.linspace(0, 1, 11), 11, 5e-2, 0.1, 0.3),
        np.array(
            [
                1.00700145,
                0.69182901,
                0.51790363,
                0.65138429,
                0.50474375,
                0.6329269,
                0.91725929,
                0.64211806,
                0.54597866,
                0.71741153,
                0.53410005,
            ]
        ),
    )


def test_generate_2D_multipass():
    X, Y = np.linspace(-1, 1, 3), np.linspace(-1, 1, 3)
    X, Y = np.meshgrid(X, Y)

    np.testing.assert_allclose(
        generate_2D_multipass(X, Y, 11, 5e-2, 0.81, 0.94, 1.1),
        np.array(
            [
                [0.29056624, 0.62260391, 0.29056624],
                [0.64293287, 1.37762912, 0.64293287],
                [0.76084547, 1.63028353, 0.76084547],
            ]
        ),
    )


def test_calculate_intensity_from_power_gaussian_beam():
    assert (
        calculate_intensity_from_power_gaussian_beam(1.0, 1e-3, 1e-3)
        == 159154.94309189534
    )


def test_calculate_power_from_rabi_gaussian_beam():
    assert (
        calculate_power_from_rabi_gaussian_beam(1e6, 0.5, 1e-3, 1e-3)
        == 5.213222807755201e-05
    )


def test_intensity_to_electric_field():
    assert intensity_to_electric_field(1.5) == 33.618312584075504


def test_calculate_rabi_from_power_gaussian_beam():
    assert (
        calculate_rabi_from_power_gaussian_beam(1, 0.5, 1e-3, 1e-3) == 138499069.8405574
    )


def test_calculate_power_from_rabi_gaussian_beam_microwave():
    assert (
        calculate_power_from_rabi_gaussian_beam_microwave(1e6, 0.5, 1e-3, 1e-3)
        == 1.8649286975200523e-06
    )


def test_calculate_rabi_from_power_gaussian_beam_microwave():
    assert (
        calculate_rabi_from_power_gaussian_beam_microwave(1, 0.5, 1e-3, 1e-3)
        == 732266024.0300486
    )


def test_conversion_power_rabi():
    Ω = 2 * 3.1415 * 1e6
    main_coupling = 0.5
    σ = 1e-3
    P = calculate_power_from_rabi_gaussian_beam(Ω, main_coupling, σ, σ)
    assert Ω == calculate_rabi_from_power_gaussian_beam(P, main_coupling, σ, σ)
