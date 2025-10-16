import pytest
import shutil
import os
import zipfile
import numpy as np
import growpacity as op
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import interpn
from pathlib import Path
from importlib.resources import files


@pytest.fixture(scope="module")
def opacity_data():
    OC = op.OpacityCalculator(
        amax_min=0.1*u.um, amax_max=1*u.um, Namax=3,
        T_min=1*u.K, T_max=3000*u.K, NT=101,
        q_min=-4.5, q_max=-2.5, Nq=9,
        optool_args='', dirc='test-data/'
    )
    OC.execute_optool(quiet=True, overwrite=True)
    OC.build_mean_opacities(overwrite=True)
    OC.compute_and_store_master_arrays(overwrite=True)
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

def test_rosseland_opacities(opacity_data):
    arrs, kR_arr, kP_arr = opacity_data
    qtest, atest, Ttest = np.array([-2.65, -3.1]), 0.5, np.array([120, 150, 2000])

    ev1 = op.evaluate_mean_opacities(*arrs, kR_arr, qtest, atest/1e4, Ttest)

    x, y, z = np.meshgrid(qtest, np.log10(atest), np.log10(Ttest), indexing='ij')
    points = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T
    ev2 = 10 ** interpn(
        (arrs[0], np.log10(arrs[1]), np.log10(arrs[2])),
        np.log10(kR_arr),
        points,
        bounds_error=True).reshape(x.shape)
    assert np.allclose(ev1, ev2, rtol=1e-6)

def test_black_body_flux():
    T = 1000 * u.K
    
    wl = op.toQuantity(np.geomspace(0.1, 10000), u.um)
    prefactor = 2 * const.h * const.c**2 / wl ** 5
    arg = (const.h * const.c / (const.k_B * wl * T)).to_value('')
    expr = np.expm1(arg)
    BT = (prefactor / expr).to(u.erg/u.s/u.cm**3)

    BT_op, uT_op = op.BBflux(T, wl=wl)
    BT_op = BT_op.to(u.erg/u.s/u.cm**3)

    assert np.allclose(BT, BT_op)

def test_toQuantity():
    arr = np.array([1, 10, 100])
    qarr = op.toQuantity(arr, u.um)
    assert np.all(qarr == arr * u.um)


@pytest.fixture(scope="module")
def tmp_opacity_dir(tmp_path_factory):
    """
    Create a temporary directory 'test-data', extract the example files, and yield its path as str.
    pytest removes tmp_path after the test run automatically.
    """
    d = tmp_path_factory.mktemp("tmp_opacity")

    # copy and extract the example files from the data folder using importlib.resources

    src = files("growpacity").joinpath("data", "default-data.zip")
    _safe_extract_zip(src, extract_dir=d)

    yield str(d / 'default-data')
    # no cleanup needed, pytest handles tmp_path


def _safe_extract_zip(zip_path, extract_dir: Path):
    """
    Safely extract a zip file to a directory, preventing path traversal attacks.
    """

    extract_dir = extract_dir.resolve()
    with zipfile.ZipFile(zip_path, "r") as zf:
        for info in zf.infolist():
            member_path = (extract_dir / info.filename).resolve()
            if not str(member_path).startswith(str(extract_dir) + os.sep) and member_path != extract_dir:
                raise Exception("Path traversal detected in zip archive")
            if info.is_dir():
                member_path.mkdir(parents=True, exist_ok=True)
            else:
                member_path.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info) as src_f, open(member_path, "wb") as dst_f:
                    shutil.copyfileobj(src_f, dst_f)


def test_C(tmp_opacity_dir, opacity_data):

    # read in the data using the C library
    oldcwd = os.getcwd()
    try:
        os.chdir(tmp_opacity_dir)
        op.growpacity_c.read_opacity_data()
    finally:
        os.chdir(oldcwd)

    # read in the data using python
    OC = op.OpacityCalculator(dirc=tmp_opacity_dir, read_in=True)
    kR_arr, kP_arr = OC.load_master_arrays()
    q = OC.q
    amax = OC.amax.to_value('um')
    T = OC.T.to_value('K')
    arrs = (q, amax, T)

    qtest, atest, Ttest = -2.65, 0.5, 120

    # call the python function
    ev1 = op.evaluate_mean_opacity(*arrs, kR_arr, qtest, atest/1e4, Ttest)

    # call the C function
    ev2 = op.growpacity_c.rosseland_opacity(qtest, atest/1e4, Ttest)

    assert np.allclose(ev1, ev2, rtol=1e-6)
