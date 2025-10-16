import numpy as np
from astropy import constants as const, units as u
import os
from numba import njit

def toQuantity(array, unit=None):
    """
    Checks if an array is an astropy Quantity. Converts it to one if not.

    In that case, `unit` must be provided in the form of either plain text
    (compatible with astropy.units.Unit), or as an astropy.units.Unit.
    
    Arguments
    ----------
    array : array-like or astropy.units.Quantity
        Input array or quantity.
    unit : str or astropy.units.Unit, optional
        The unit to convert to if `array` is not a Quantity.

    Returns
    -------
    astropy.units.Quantity
        The input array as a Quantity.
    """
    if isinstance(array, u.Quantity): return array.to(unit)
    else: return np.array(array) * u.Unit(unit)

def BBflux(T=5780, wl=None):
    """
    Returns the flux spectrum of a black body at wavelength λ and temperature T.
    The flux is calculated assuming::

        f = exp(h*c/(λkT))
        B(λ) = 2hc^2/λ^5 / (f - 1)
        u(λ) = 2h^2c^3 / (kλ^6T^2) * exp / (exp-1)

    Arguments
    ----------
    T : array-like or astropy.units.Quantity
        the temperature of the black body. Defaults to K.
    wl : array-like or astropy.units.Quantity
        the wavelengths to sample on. Defaults to cm.
    
    Returns
    -------
    B : astropy.units.Quantity [erg/(cm^2 s cm)]
        the black body flux at the requested wavelengths
    dB/dT : astropy.units.Quantity [erg/(cm^2 s K cm)]
        the derivative of the black body flux with respect to temperature
    """

    tmp = toQuantity(T, u.K)
    lam = toQuantity(wl, u.cm)

    h = const.h
    c = const.c
    kB= const.k_B

    max_limit   = np.log(1e308)
    expr = np.minimum(h*c/(lam*kB*tmp), max_limit)

    constant = 2*h*c*c/lam**5
    B = constant / np.expm1(expr)
    B = B.to('erg/(cm2*s *cm)')

    constant = 2*h**2*c**3/(kB*lam**6*tmp**2)
    with np.errstate(over='ignore'):
        frac = np.exp(expr)/np.expm1(expr)**2
    dBdT = constant * frac
    dBdT = dBdT.to('erg/(cm2*s*K *cm)')

    return B, dBdT
    
                
class OpacityCalculator:
    """
    Master class that handles the opacity calculations.
    This is the meat of the code: basically call OpTool repeatedly 
    over the range of grain size distributions requested, and then
    compute Rosseland and Planck means and store them into files.
    I'm defining it as a class, in case you would like to use it in
    python format in another project where you want access to metadata.
    
    Input (when marked with [...], can be an astropy Quantity that defaults to unit "...")

    Arguments
    ----------
    amin : float or astropy.units.Quantity [μm]
        minimum grain size for all distributions
    q_min : float
        smallest powerlaw exponent within the range of distributions
    q_max : float
        largest  powerlaw exponent within the range of distributions
    Nq : int
        number of samples between q_min and q_max. Linear spacing.
    amax_min : float or astropy.units.Quantity [μm]
        smallest maximum grain size within the range of distributions
    amax_max : float or astropy.units.Quantity [μm]
        largest  maximum grain size within the range of distributions
    Namax : int
        number of samples between amax_min and amax_max. Logarithmic spacing.
    T_min : float or astropy.units.Quantity [K]
        smallest temperature to sample when computing means
    T_max : float or astropy.units.Quantity [K]
        largest  temperature to sample when computing means
    NT : int
        number of samples between T_min and T_max. Logarithmic spacing.
    optool_args : str
        additional arguments to be passed to optool,
        beyond grain size distributions (e.g., composition).
    dirc : str
        directory to save the data
        (about 20KB/file for absorption coefficients, 230B*NT/file for mean opacities).
    name : str
        identifier to be appended to the file name.

    The complete filenames will be:
        {dirc}/{name}_a{amax}_q{q}.inp for absorption coefficients
        {dirc}/{name}_a{amax}_q{q}.dat for mean opacities
    """

    def __init__(self, amin=0.1*u.um,
                 q_min=-4.5, q_max=-2.5, Nq=9,
                 amax_min=0.1*u.um, amax_max=1e4*u.um, Namax=14,
                 T_min=1*u.K, T_max=3000*u.K, NT=100,
                 optool_args='', dirc='data/', name=''):
        
        self.amin = toQuantity(amin, u.um)
        self.q = np.linspace(q_min, q_max, Nq)
        self.amax = np.geomspace(toQuantity(amax_min, u.um),
                                toQuantity(amax_max, u.um), Namax)
        self.T = np.geomspace(toQuantity(T_min, u.K),
                               toQuantity(T_max, u.K), NT)
        
        self.optool_args = optool_args

        self.dirc = dirc
        self.name = name

    def get_filename(self, amax, q):
        """
        Helper function that returns the file name
        (excluding prepended 'dustkappa/kappaRP')
        for a given amax and q.
        """

        sep = '_' if len(self.name) else ''
        full_name = f"{self.name}{sep}a{amax:.2e}_q{q:.2f}"
        return full_name

    def execute_optool(self, quiet=True, overwrite=False):
        """
        Repeatedly calls OpTool over the range of amax and q requested.
        If quiet, suppressed OpTool output. Useful to hide all the dots.
        """

        dirc = self.dirc
        name = self.name
        amin = self.amin.to_value('um')

        if not os.path.exists(dirc): os.makedirs(dirc)

        for amax in self.amax.to_value('um'):
            for q in self.q:
                full_name = self.get_filename(amax, q)
                filename  = f"{dirc}/dustkappa_{full_name}.inp"
                if os.path.exists(filename) and not overwrite: continue

                cmd = f'optool -a {amin} {amax} {-q} ' + self.optool_args
                cmd += f' -o {dirc} -radmc {full_name}'
                if quiet: cmd += ' >/dev/null 2>&1' # ignore default output

                print(f'Running: {cmd}')
                os.system(cmd)
    
    def __compute_mean_opacities(self, filename):
        """
        Internal function that computes Rosseland and Planck means
        from a given file with absorption coefficients.
        """

        # first skip to start of the file
        idx = 0
        with open(filename, 'r') as f:
            while idx < 100:
                l = f.readline()
                if not l.startswith('#'): break
                idx += 1
            idx += 2 # radmc format lines

        # read data and apply units
        wl, kabs, ksca, g = np.genfromtxt(filename, skip_header=idx, unpack=True)
        wl, kabs, ksca = wl * u.um, kabs * u.cm**2/u.g, ksca * u.cm**2/u.g

        # compute Rosseland and Planck means
        B, dB = BBflux(self.T[:,None], wl=wl)
        def trapz(y, x): return np.trapezoid(y, x=x, axis=1)
        kR = (trapz(dB, wl) / trapz(dB/(kabs+ksca*(1-g)), wl)).to('cm2/g')
        kP = (trapz(B*kabs, wl) / trapz(B, wl)).to('cm2/g')
        return kR, kP

    def build_mean_opacities(self, overwrite=False):
        """
        Loop over all files of absorption coefficients and print ASCII
        files with the respective temperature-dependent mean opacities.
        If not overwrite, it skips existing files with that name.
        """

        dirc = self.dirc
        name = self.name
        T    = self.T.to('K')

        for amax in self.amax.to_value('um'):
            for q in self.q:
                full_name = self.get_filename(amax, q)
                opac_name = f"{dirc}/kappaRP_{full_name}.dat"
                if os.path.exists(opac_name) and not overwrite: continue

                # read absorption opacities from file
                filename  = f"{dirc}/dustkappa_{full_name}.inp"

                kR, kP = self.__compute_mean_opacities(filename)
                data = np.array([T.to_value('K'), kR.to_value('cm2/g'), kP.to_value('cm2/g')])
                np.savetxt(opac_name, data.T)

    def compute_master_arrays(self):
        """
        Loads all files of mean opacities and combines them into master arrays
        of κ(q, amax, T).
        """

        dirc = self.dirc
        name = self.name
        T    = self.T.to('K')

        shape = (self.q.size, self.amax.size, self.T.size)
        kR_arr = np.zeros(shape, dtype=float)
        kP_arr = np.zeros_like(kR_arr)

        for q_idx, q in enumerate(self.q):
            for a_idx, amax in enumerate(self.amax.to_value('um')):
                full_name = self.get_filename(amax, q)
                opac_name = f"{dirc}/kappaRP_{full_name}.dat"

                T_, kR, kP = np.genfromtxt(opac_name, unpack=True)
                kR_arr[q_idx, a_idx] = kR
                kP_arr[q_idx, a_idx] = kP

        return kR_arr, kP_arr

    def compute_and_store_master_arrays(self, overwrite=False):
        """
        Stores 1D arrays of q, amax T, and the 3D arrays of 
        kappa_Rosseland(q, amax, T) and kappa_Planck(q, amax, T).
        If overwrite and a file exists, skip it.
        """
        
        dirc = self.dirc
        have = os.path.exists

        def write_array(flnm, arr):
            if have(flnm) and not overwrite: return
            with open(flnm, 'w') as f:
                f.write(f'{arr.size}\n')
                np.savetxt(f, arr)

        # these are in ASCII
        write_array(f'{dirc}/q.dat', self.q)
        write_array(f'{dirc}/amax_um.dat', self.amax.to_value('um'))
        write_array(f'{dirc}/T_K.dat', self.T.to_value('K'))

        # these are binary

        flnm_kR, flnm_kP = f'{dirc}/kR_cm2g.dbl', f'{dirc}/kP_cm2g.dbl'
        if have(flnm_kR) and have(flnm_kP) and (not overwrite): return

        kR_arr, kP_arr = self.compute_master_arrays()
        if not have(flnm_kR) or overwrite: kR_arr.tofile(flnm_kR)
        if not have(flnm_kP) or overwrite: kP_arr.tofile(flnm_kP)

    def load_master_arrays(self):
        """If master arrays have already been computed, load them directly."""
        dirc = self.dirc

        shape = (self.q.size, self.amax.size, self.T.size)
        kR_arr = np.fromfile(f'{dirc}/kR_cm2g.dbl').reshape(shape)
        kP_arr = np.fromfile(f'{dirc}/kP_cm2g.dbl').reshape(shape)

        return kR_arr, kP_arr

@njit
def evaluate_mean_opacity(q_arr, amax_arr, T_arr, kappa_arr, q, amax, T):
    """
    Given the arrays of q, amax, T, and kappa(q, amax, T),
    evaluates kappa at the requested values using
    trilinear interpolation.

    Note that amax_arr is expected in μm, and T_arr in K.
    This is done to match the output of OpacityCalculator.
    However, amax (the requested value) is expected in cm,
    and will be converted to μm internally.

    This function is JIT-compiled with numba for speed.

    Arguments
    ----------
    q_arr
        1D array of powerlaw exponents
    amax_arr
        1D array of maximum grain sizes [um]
    T_arr
        1D array of temperatures [K]
    kappa_arr
        3D array of opacities [cm^2/g]
    q
        requested powerlaw exponent
    amax
        requested maximum grain size [cm]
    T
        requested temperature [K]

    Returns
    -------
    kappa
        the interpolated opacity [cm^2/g]
    """

    def get_idx(arr, v):
        """
        Assuming regular spacing, compute the index
        to the left of the requested value using
        a simple linear interpolation.
        """
        Amin, Amax, NA = arr.min(), arr.max(), arr.size
        return int((NA-1) * (v-Amin) / (Amax-Amin))

    def get_edges(N, idx):
        """
        Returns the neighboring indices around the target value.
        Also clamps the arrays to avoid over/undershooting.
        """
        if idx < 0: return 0, 0 # undershooting
        if idx > N - 1: return N-1, N-1 # overshooting
        Lidx = min(max(idx, 0), N-1) # normal
        return Lidx, Lidx+1

    def get_weights(vL, vR, dv, v):
        """linear interpolation weights"""
        w = (v - vL) / dv
        return 1-w, w

    lamax_arr = np.log10(amax_arr)
    lT_arr    = np.log10(T_arr)

    lamax, lT = np.log10(amax*1e4), np.log10(T)

    qLidx, qRidx = get_edges(q_arr.size,     get_idx(q_arr, q))
    aLidx, aRidx = get_edges(lamax_arr.size, get_idx(lamax_arr, lamax))
    TLidx, TRidx = get_edges(lT_arr.size,    get_idx(lT_arr, lT))

    def get_dv(arr): return (arr[-1] - arr[0]) / (arr.size - 1)
    wqL, wqR = get_weights(q_arr[qLidx], q_arr[qRidx], get_dv(q_arr), q)
    waL, waR = get_weights(lamax_arr[aLidx], lamax_arr[aRidx], get_dv(lamax_arr), lamax)
    wTL, wTR = get_weights(lT_arr[TLidx], lT_arr[TRidx], get_dv(lT_arr), lT)

    # we're interpolating in q, logamax, logT, logkappa space because
    # kappa is roughly a power law in amax and T
    # indexing is q, a, T
    kappaLLL, kappaLLR = kappa_arr[qLidx, aLidx, TLidx], kappa_arr[qLidx, aLidx, TRidx]
    kappaLRL, kappaLRR = kappa_arr[qLidx, aRidx, TLidx], kappa_arr[qLidx, aRidx, TRidx]
    kappaRLL, kappaRLR = kappa_arr[qRidx, aLidx, TLidx], kappa_arr[qRidx, aLidx, TRidx]
    kappaRRL, kappaRRR = kappa_arr[qRidx, aRidx, TLidx], kappa_arr[qRidx, aRidx, TRidx]

    # collapse to index a, T using q weights
    lkappaLL = np.log10(kappaLLL) * wqL + np.log10(kappaRLL) * wqR
    lkappaLR = np.log10(kappaLLR) * wqL + np.log10(kappaRLR) * wqR
    lkappaRL = np.log10(kappaLRL) * wqL + np.log10(kappaRRL) * wqR
    lkappaRR = np.log10(kappaLRR) * wqL + np.log10(kappaRRR) * wqR

    # collapse to index T using amax weights
    lkappaL  = lkappaLL * waL + lkappaRL * waR
    lkappaR  = lkappaLR * waL + lkappaRR * waR

    # collapse to final value using T weights
    lkappa = lkappaL * wTL + lkappaR * wTR

    kappa = 10.0 ** lkappa
    return kappa

def evaluate_mean_opacities(q_arr, amax_arr, T_arr, kappa_arr, q, amax, T):
    """
    Wrapper around the numba-jitted function `evaluate_mean_opacity`
    to allow for array-like inputs for q, amax, and T.

    Arguments
    ----------
    q_arr : 1D array
        1D array of powerlaw exponents
    amax_arr : 1D array
        1D array of maximum grain sizes [um]
    T_arr : 1D array
        1D array of temperatures [K]
    kappa_arr : 3D array
        3D array of opacities [cm^2/g]
    q : float or array-like
        requested powerlaw exponent
    amax : float or array-like
        requested maximum grain size [cm]
    T : float or array-like
        requested temperature [K]

    Returns
    -------
    kappa
        the interpolated opacity [cm^2/g]
    """

    q_ = np.atleast_1d(q)
    amax_ = np.atleast_1d(amax)
    T_ = np.atleast_1d(T)
    kappa = __evaluate_mean_opacities(q_arr, amax_arr, T_arr, kappa_arr, q_, amax_, T_)
    if kappa.size == 1: return kappa[0,0,0]
    return kappa

@njit
def __evaluate_mean_opacities(q_arr, amax_arr, T_arr, kappa_arr, q_, amax_, T_):
    """
    Internal function that loops over all requested values (as arrays)
    and calls the jitted `evaluate_mean_opacity` function.
    """
    shape = (q_.size, amax_.size, T_.size)
    kappa = np.zeros(shape, dtype=float)

    for k in range(q_.size):
        for j in range(amax_.size):
            for i in range(T_.size):
                kappa[k,j,i] = evaluate_mean_opacity(q_arr, amax_arr, T_arr, kappa_arr,
                                                    q_[k], amax_[j], T_[i])
    return kappa