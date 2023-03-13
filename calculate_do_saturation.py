# Author: A.Wachholz
# Date: 19.10.22
"""
This script calculates the dissolved oxygen (DO) saturation in percent based on
the temperature, atmospheric pressure and the DO concentration in mg/l.
 The approach used is based on Weiss's method as described in the paper
"https://water.usgs.gov/water-resources/memos/documents/WQ.2011.03.pdf"
"""

# imports
import numpy as np
from matplotlib import pyplot as plt

# define constants
ABSOLUTE_ZERO = - 273.15

# define functions


def deg2kelvin(x):
    """
    This function converts temperature in degree celsius to kelvin
    :param x: temperature in degree celsius
    :return: temperature in kelvin
    """
    return x - ABSOLUTE_ZERO


def kelvin2deg(x):
    """
    This function converts temperature in kelvin to degree celsius
    :param x: temperature in kelvin
    :return: temperature in degree
    """
    return x + ABSOLUTE_ZERO


def vapor_pressure_of_water(t):
    """
    This function calculates the vapor pressure of water at a given temperature
    :param t: temperature in degree celsius
    :return: vapor pressure of water in mgHg (or Torr)
    """
    u = 10 ** (
        8.10765 - (1750.286/(235+t))
    )
    return u


def pressure_correction(p, t):
    """
    This function calculates the correction factor for atmospheric pressure
    :param p: atmospheric pressure in mgHg (or Torr)
    :param t: water temperature t in degree Celsius
    :return: factor F_p
    """
    f_p = (p - vapor_pressure_of_water(t)) / (760 - vapor_pressure_of_water(t))
    return f_p


def do_at_saturation(t, p):
    """
    This function calculates the dissolved oxygen (DO) conc in mg/l at saturation
    :param t: temperature in degree celsius
    :param p: atmospheric pressure in Torr or mmHg
    :return: DO conc in mg/l at saturation
    """
    # convert degree Celsius to Kelvin
    t_Kelvin = deg2kelvin(t)

    do_sat = 1.42905 * np.exp(
        - 173.4292
        + 249.6339 * (100/t_Kelvin)
        + 143.3483 * np.log((t_Kelvin/100))
        - 21.8493 * (t_Kelvin/100)
    )
    return do_sat * pressure_correction(p=p, t=t)


def do_saturation(do, t, p):
    """
    This function calculates DO saturation in percent based on water temperature
    and air pressure
    :param do: DO concentration in mg/l
    :param t: temperature in degree Celsius
    :param p: atmospheric pressure in mmHg or Torr
    :return: DO saturation in percent
    """

    do_max = do_at_saturation(t=t, p=p)
    return (do / do_max) * 100


def unit_test(
        func: callable,
        test_input: list[tuple],
        test_outputs: list[float]
) -> None:
    """
    Perform unit testing on a given function.

    Args:
    - func: callable function to be tested
    - test_input: list of input values in tuple format to be passed to the function
    - test_outputs: list of expected output values for the function given the provided inputs

    Returns:
    None

    Raises:
    ValueError: if the output of the function does not match the expected output.
    The error message includes the name of the function that failed the test.
    """

    test_returns = [np.round(func(*[val]), 1) for val in test_input]

    if not test_outputs == test_returns:
        raise ValueError('function', func.__name__, 'failed test')
    else:
        pass


if __name__ == '__main__':

    # unit testing
    # * vapor pressure of water
    #   (see https://en.wikipedia.org/wiki/Vapour_pressure_of_water)
    unit_test(
        func=vapor_pressure_of_water,
        test_input=[0, 5, 50],
        test_outputs=[4.6, 6.5, 92.5]
    )
    # * pressure correction
    #   should be 1 at one atmosphere (760 mm Hg) and 22 deg Celsius
    if not np.round(pressure_correction(p=760, t=22), 1) == np.round(1, 1):
        raise ValueError('function', 'pressure_correction', 'failed')

    # * do_at_saturation
    #   should be ca 9 at 20 deg celsius and 760 mmHg
    if not np.round(do_at_saturation(t=20, p=760), 1) == 9.1:
        raise ValueError('function', 'do_at_saturation', 'failed')

    # test plot
    # define input
    do_mgl_space = np.linspace(5, 12, 100)
    temperature_DegCelsius_space = np.linspace(0, 30, 100)
    pressure_mmHg = 760

    sats = []
    dos = []
    ts = []
    for do in do_mgl_space:
        for t in temperature_DegCelsius_space:
            sat = do_saturation(do=do, t=t, p=pressure_mmHg)
            sats.append(sat)
            dos.append(do)
            ts.append(t)

    sc = plt.scatter(
        dos,
        ts,
        c=sats
    )
    plt.xlabel('DO [mg/l]')
    plt.ylabel('water temperature [$\degree$ C]')
    plt.colorbar(mappable=sc, label='DO saturation [%]')
    plt.show()



