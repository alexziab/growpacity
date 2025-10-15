import pytest
import numpy as np
import growpacity as op
import astropy.units as u
from scipy.interpolate import interpn


@pytest.fixture(scope="module")
def opacity_data():
    OC = op.OpacityCalculator(
        amax_min=0.1*u.um, amax_max=1*u.um, Namax=3,
        T_min=1*u.K, T_max=3000*u.K, NT=101,
        q_min=-4.5, q_max=-2.5, Nq=9,
        optool_args='', dirc='test-data/'
    )
    OC.execute_optool(quiet=True)
    OC.build_mean_opacities()
    OC.compute_and_store_master_arrays()
    kR_arr, kP_arr = OC.load_master_arrays()
    q    = OC.q
    amax = OC.amax.to_value('um')
    T    = OC.T.to_value('K')
    arrs = (q, amax, T)
    return arrs, kR_arr, kP_arr


def test_rosseland_opacity(opacity_data):
    arrs, kR_arr, kP_arr = opacity_data
    qtest, atest, Ttest = -2.65, 0.5, 120

    ev1 = op.evaluate_mean_opacity(*arrs, kR_arr, qtest, atest/1e4, Ttest)
    ev2 = 10 ** interpn(
        (arrs[0], np.log10(arrs[1]), np.log10(arrs[2])),
        np.log10(kR_arr),
        [(qtest, np.log10(atest), np.log10(Ttest))],
        bounds_error=True)[0]
    assert np.isclose(ev1, ev2, rtol=1e-6)


def test_planck_opacity(opacity_data):
    arrs, kR_arr, kP_arr = opacity_data
    qtest, atest, Ttest = -2.65, 0.5, 120

    ev1 = op.evaluate_mean_opacity(*arrs, kP_arr, qtest, atest/1e4, Ttest)
    ev2 = 10 ** interpn(
        (arrs[0], np.log10(arrs[1]), np.log10(arrs[2])),
        np.log10(kP_arr),
        [(qtest, np.log10(atest), np.log10(Ttest))],
        bounds_error=True)[0]
    assert np.isclose(ev1, ev2, rtol=1e-6)
