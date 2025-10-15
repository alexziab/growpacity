import numpy as np
import growpacity as op
import astropy.units as u
from scipy.interpolate import interpn

if __name__ == "__main__":

    OC = op.OpacityCalculator(amax_min=0.1*u.um, amax_max=1*u.um, Namax=3, # 2 bins per decade
                               T_min=1*u.K, T_max=3000*u.K, NT=101, # good enough
                               q_min=-4.5, q_max=-2.5, Nq=9,                           optool_args='', dirc='test-data/')

    OC.execute_optool(quiet=True) # run optool if necessary
    OC.build_mean_opacities() # build mean opacities if necessary
    OC.compute_and_store_master_arrays() # build master arrays if necessary
    kR_arr, kP_arr = OC.load_master_arrays() # load the master arrays

    q    = OC.q
    amax = OC.amax.to_value('um')
    T    = OC.T.to_value('K')
    arrs = (q, amax, T)

    # check interpolator against scipy: pick some random point
    qtest, atest, Ttest = -2.65, 0.5, 120 # amax in um, T in K

    # note that evaluate_mean_opacity expects amax in cm
    ev1 = op.evaluate_mean_opacity(*arrs, kR_arr, qtest, atest/1e4, Ttest)

    ev2 = 10 ** interpn((q, np.log10(amax), np.log10(T)), np.log10(kR_arr),
                (qtest, np.log10(atest), np.log10(Ttest)), bounds_error=True)[0]

    print("Rosseland:")
    print(f'Computed:   {ev1:.16e}')
    print(f'Scipy:      {ev2:.16e}')
    print(f'Difference: {ev1/ev2-1:.16e}')

    # note that evaluate_mean_opacity expects amax in cm
    ev1 = op.evaluate_mean_opacity(*arrs, kP_arr, qtest, atest/1e4, Ttest)

    ev2 = 10 ** interpn((q, np.log10(amax), np.log10(T)), np.log10(kP_arr),
                (qtest, np.log10(atest), np.log10(Ttest)), bounds_error=True)[0]

    print("Planck:")
    print(f'Computed:   {ev1:.16e}')
    print(f'Scipy:      {ev2:.16e}')
    print(f'Difference: {ev1/ev2-1:.16e}')

