import numpy as np

from openmc.mgxs.groups import EnergyGroups
from openmc.mgxs.library import Library
from openmc.mgxs.mgxs import *
from openmc.mgxs.mdgxs import *

GROUP_STRUCTURES = {}
"""Dictionary of commonly used energy group structures:

- "CASMO-X" (where X is 2, 4, 8, 16, 25, 40 or 70) from the CASMO_ lattice
  physics code
- "XMAS-172_" designed for LWR analysis ([SAR1990]_, [SAN2004]_)
- "SHEM-361_" designed for LWR analysis to eliminate self-shielding calculations
  of thermal resonances ([HFA2005]_, [SAN2007]_, [HEB2008]_)
- activation_ energy group structures "VITAMIN-J-175", "TRIPOLI-315",
  "CCFE-709_" and "UKAEA-1102_"

.. _CASMO: https://www.studsvik.com/SharepointFiles/CASMO-5%20Development%20and%20Applications.pdf
.. _XMAS-172: https://www-nds.iaea.org/wimsd/energy.htm
.. _SHEM-361: https://www.polymtl.ca/merlin/downloads/FP214.pdf
.. _activation: https://fispact.ukaea.uk/wiki/Keyword:GETXS
.. _CCFE-709: https://fispact.ukaea.uk/wiki/CCFE-709_group_structure
.. _UKAEA-1102: https://fispact.ukaea.uk/wiki/UKAEA-1102_group_structure
.. [SAR1990] Sartori, E., OECD/NEA Data Bank: Standard Energy Group Structures
   of Cross Section Libraries for Reactor Shielding, Reactor Cell and Fusion
   Neutronics Applications: VITAMIN-J, ECCO-33, ECCO-2000 and XMAS JEF/DOC-315
   Revision 3 - DRAFT (December 11, 1990).
.. [SAN2004] Santamarina, A., Collignon, C., & Garat, C. (2004). French
   calculation schemes for light water reactor analysis. United States:
   American Nuclear Society - ANS.
.. [HFA2005] Hfaiedh, N. & Santamarina, A., "Determination of the Optimized
   SHEM Mesh for Neutron Transport Calculations," Proc. Top. Mtg. in
   Mathematics & Computations, Supercomputing, Reactor Physics and Nuclear and
   Biological Applications, September 12-15, Avignon, France, 2005.
.. [SAN2007] Santamarina, A. & Hfaiedh, N. (2007). The SHEM energy mesh for
   accurate fuel depletion and BUC calculations. Proceedings of the International
   Conference on Safety Criticality ICNC 2007, St Peterburg (Russia), Vol. I pp.
   446-452.
.. [HEB2008] Hébert, Alain & Santamarina, Alain. (2008). Refinement of the
   Santamarina-Hfaiedh energy mesh between 22.5 eV and 11.4 keV. International
   Conference on the Physics of Reactors 2008, PHYSOR 08. 2. 929-938.
"""

GROUP_STRUCTURES['CASMO-2'] = np.array([
    0., 6.25e-1, 2.e7])
GROUP_STRUCTURES['CASMO-4'] = np.array([
    0., 6.25e-1, 5.53e3, 8.21e5, 2.e7])
GROUP_STRUCTURES['CASMO-8'] = np.array([
    0., 5.8e-2, 1.4e-1, 2.8e-1, 6.25e-1, 4., 5.53e3, 8.21e5, 2.e7])
GROUP_STRUCTURES['CASMO-16'] = np.array([
    0., 3.e-2, 5.8e-2, 1.4e-1, 2.8e-1, 3.5e-1, 6.25e-1, 8.5e-1,
    9.72e-1, 1.02, 1.097, 1.15, 1.3,  4., 5.53e3, 8.21e5, 2.e7])
GROUP_STRUCTURES['CASMO-25']  = np.array([
    0., 3.e-2, 5.8e-2, 1.4e-1, 2.8e-1, 3.5e-1, 6.25e-1, 9.72e-1, 1.02, 1.097,
    1.15, 1.855, 4., 9.877, 1.5968e1, 1.4873e2, 5.53e3, 9.118e3, 1.11e5, 5.e5,
    8.21e5, 1.353e6, 2.231e6, 3.679e6, 6.0655e6, 2.e7])
GROUP_STRUCTURES['CASMO-40'] = np.array([
    0., 1.5e-2, 3.e-2, 4.2e-2, 5.8e-2, 8.e-2, 1.e-1, 1.4e-1,
    1.8e-1, 2.2e-1, 2.8e-1, 3.5e-1, 6.25e-1, 8.5e-1, 9.5e-1,
    9.72e-1, 1.02, 1.097, 1.15, 1.3, 1.5, 1.855, 2.1, 2.6, 3.3, 4.,
    9.877, 1.5968e1, 2.77e1, 4.8052e1, 1.4873e2, 5.53e3, 9.118e3,
    1.11e5, 5.e5, 8.21e5, 1.353e6, 2.231e6, 3.679e6, 6.0655e6, 2.e7])
GROUP_STRUCTURES['CASMO-70'] = np.array([
    0., 5.e-3, 1.e-2, 1.5e-2, 2.e-2, 2.5e-2, 3.e-2, 3.5e-2, 4.2e-2,
    5.e-2, 5.8e-2, 6.7e-2, 8.e-2, 1.e-1, 1.4e-1, 1.8e-1, 2.2e-1,
    2.5e-1, 2.8e-1, 3.e-1, 3.2e-1, 3.5e-1, 4.e-1, 5.e-1, 6.25e-1,
    7.8e-1, 8.5e-1, 9.1e-1, 9.5e-1, 9.72e-1, 9.96e-1, 1.02, 1.045,
    1.071, 1.097, 1.123, 1.15, 1.3, 1.5, 1.855, 2.1, 2.6, 3.3, 4.,
    9.877, 1.5968e1, 2.77e1, 4.8052e1, 7.5501e1, 1.4873e2,
    3.6726e2, 9.069e2, 1.4251e3, 2.2395e3, 3.5191e3, 5.53e3,
    9.118e3, 1.503e4, 2.478e4, 4.085e4, 6.734e4, 1.11e5, 1.83e5,
    3.025e5, 5.e5, 8.21e5, 1.353e6, 2.231e6, 3.679e6, 6.0655e6, 2.e7])
GROUP_STRUCTURES['XMAS-172'] = np.array([
    1.00001e-05, 3.00000e-03, 5.00000e-03, 6.90000e-03, 1.00000e-02,
    1.50000e-02, 2.00000e-02, 2.50000e-02, 3.00000e-02, 3.50000e-02,
    4.20000e-02, 5.00000e-02, 5.80000e-02, 6.70000e-02, 7.70000e-02,
    8.00000e-02, 9.50000e-02, 1.00001e-01, 1.15000e-01, 1.34000e-01,
    1.40000e-01, 1.60000e-01, 1.80000e-01, 1.89000e-01, 2.20000e-01,
    2.48000e-01, 2.80000e-01, 3.00000e-01, 3.14500e-01, 3.20000e-01,
    3.50000e-01, 3.91000e-01, 4.00000e-01, 4.33000e-01, 4.85000e-01,
    5.00000e-01, 5.40000e-01, 6.25000e-01, 7.05000e-01, 7.80000e-01,
    7.90000e-01, 8.50000e-01, 8.60000e-01, 9.10000e-01, 9.30000e-01,
    9.50000e-01, 9.72000e-01, 9.86000e-01, 9.96000e-01, 1.02000e+00,
    1.03500e+00, 1.04500e+00, 1.07100e+00, 1.09700e+00, 1.11000e+00,
    1.12535e+00, 1.15000e+00, 1.17000e+00, 1.23500e+00, 1.30000e+00,
    1.33750e+00, 1.37000e+00, 1.44498e+00, 1.47500e+00, 1.50000e+00,
    1.59000e+00, 1.67000e+00, 1.75500e+00, 1.84000e+00, 1.93000e+00,
    2.02000e+00, 2.10000e+00, 2.13000e+00, 2.36000e+00, 2.55000e+00,
    2.60000e+00, 2.72000e+00, 2.76792e+00, 3.30000e+00, 3.38075e+00,
    4.00000e+00, 4.12925e+00, 5.04348e+00, 5.34643e+00, 6.16012e+00,
    7.52398e+00, 8.31529e+00, 9.18981e+00, 9.90555e+00, 1.12245e+01,
    1.37096e+01, 1.59283e+01, 1.94548e+01, 2.26033e+01, 2.49805e+01,
    2.76077e+01, 3.05113e+01, 3.37201e+01, 3.72665e+01, 4.01690e+01,
    4.55174e+01, 4.82516e+01, 5.15780e+01, 5.55951e+01, 6.79041e+01,
    7.56736e+01, 9.16609e+01, 1.36742e+02, 1.48625e+02, 2.03995e+02,
    3.04325e+02, 3.71703e+02, 4.53999e+02, 6.77287e+02, 7.48518e+02,
    9.14242e+02, 1.01039e+03, 1.23410e+03, 1.43382e+03, 1.50733e+03,
    2.03468e+03, 2.24867e+03, 3.35463e+03, 3.52662e+03, 5.00451e+03,
    5.53084e+03, 7.46586e+03, 9.11882e+03, 1.11378e+04, 1.50344e+04,
    1.66156e+04, 2.47875e+04, 2.73944e+04, 2.92830e+04, 3.69786e+04,
    4.08677e+04, 5.51656e+04, 6.73795e+04, 8.22975e+04, 1.11090e+05,
    1.22773e+05, 1.83156e+05, 2.47235e+05, 2.73237e+05, 3.01974e+05,
    4.07622e+05, 4.50492e+05, 4.97871e+05, 5.50232e+05, 6.08101e+05,
    8.20850e+05, 9.07180e+05, 1.00259e+06, 1.10803e+06, 1.22456e+06,
    1.35335e+06, 1.65299e+06, 2.01897e+06, 2.23130e+06, 2.46597e+06,
    3.01194e+06, 3.67879e+06, 4.49329e+06, 5.48812e+06, 6.06531e+06,
    6.70320e+06, 8.18731e+06, 1.00000e+07, 1.16183e+07, 1.38403e+07,
    1.49182e+07, 1.73325e+07, 1.96403e+07])
GROUP_STRUCTURES['VITAMIN-J-175'] = np.array([
    1.0000e-5, 1.0000e-1, 4.1399e-1, 5.3158e-1, 6.8256e-1,
    8.7643e-1, 1.1253, 1.4450, 1.8554, 2.3824, 3.0590,
    3.9279, 5.0435, 6.4759, 8.3153, 1.0677e1, 1.3710e1,
    1.7604e1, 2.2603e1, 2.9023e1, 3.7266e1, 4.7851e1, 6.1442e1,
    7.8893e1, 1.0130e2, 1.3007e2, 1.6702e2, 2.1445e2, 2.7536e2,
    3.5358e2, 4.5400e2, 5.8295e2, 7.4852e2, 9.6112e2, 1.2341e3,
    1.5846e3, 2.0347e3, 2.2487e3, 2.4852e3, 2.6126e3, 2.7465e3,
    3.0354e3, 3.3546e3, 3.7074e3, 4.3074e3, 5.5308e3, 7.1017e3,
    9.1188e3, 1.0595e4, 1.1709e4, 1.5034e4, 1.9304e4, 2.1875e4,
    2.3579e4, 2.4176e4, 2.4788e4, 2.6058e4, 2.7000e4, 2.8501e4,
    3.1828e4, 3.4307e4, 4.0868e4, 4.6309e4, 5.2475e4, 5.6562e4,
    6.7380e4, 7.2024e4, 7.9499e4, 8.2503e4, 8.6517e4, 9.8036e4,
    1.1109e5, 1.1679e5, 1.2277e5, 1.2907e5, 1.3569e5, 1.4264e5,
    1.4996e5, 1.5764e5, 1.6573e5, 1.7422e5, 1.8316e5, 1.9255e5,
    2.0242e5, 2.1280e5, 2.2371e5, 2.3518e5, 2.4724e5, 2.7324e5,
    2.8725e5, 2.9452e5, 2.9721e5, 2.9849e5, 3.0197e5, 3.3373e5,
    3.6883e5, 3.8774e5, 4.0762e5, 4.5049e5, 4.9787e5, 5.2340e5,
    5.5023e5, 5.7844e5, 6.0810e5, 6.3928e5, 6.7206e5, 7.0651e5,
    7.4274e5, 7.8082e5, 8.2085e5, 8.6294e5, 9.0718e5, 9.6167e5,
    1.0026e6, 1.1080e6, 1.1648e6, 1.2246e6, 1.2874e6, 1.3534e6,
    1.4227e6, 1.4957e6, 1.5724e6, 1.6530e6, 1.7377e6, 1.8268e6,
    1.9205e6, 2.0190e6, 2.1225e6, 2.2313e6, 2.3069e6, 2.3457e6,
    2.3653e6, 2.3851e6, 2.4660e6, 2.5924e6, 2.7253e6, 2.8650e6,
    3.0119e6, 3.1664e6, 3.3287e6, 3.6788e6, 4.0657e6, 4.4933e6,
    4.7237e6, 4.9658e6, 5.2205e6, 5.4881e6, 5.7695e6, 6.0653e6,
    6.3763e6, 6.5924e6, 6.7032e6, 7.0469e6, 7.4082e6, 7.7880e6,
    8.1873e6, 8.6071e6, 9.0484e6, 9.5123e6, 1.0000e7, 1.0513e7,
    1.1052e7, 1.1618e7, 1.2214e7, 1.2523e7, 1.2840e7, 1.3499e7,
    1.3840e7, 1.4191e7, 1.4550e7, 1.4918e7, 1.5683e7, 1.6487e7,
    1.6905e7, 1.7332e7, 1.9640e7])
GROUP_STRUCTURES['TRIPOLI-315,'] = np.array([
    1.0e-5, 1.1e-4, 3.000e-3, 5.500e-3, 1.000e-2, 1.500e-2, 2.000e-2, 3.000e-2,
    3.200e-2, 3.238e-2, 4.300e-2, 5.900e-2, 7.700e-2, 9.500e-2, 1.000e-1,
    1.150e-1, 1.340e-1, 1.600e-1, 1.890e-1, 2.200e-1, 2.480e-1, 2.825e-1,
    3.145e-1, 3.520e-1, 3.910e-1, 4.140e-1, 4.330e-1, 4.850e-1, 5.316e-1,
    5.400e-1, 6.250e-1, 6.826e-1, 7.050e-1, 7.900e-1, 8.600e-1, 8.764e-1,
    9.300e-1, 9.860e-1, 1.010, 1.035, 1.070, 1.080, 1.090,
    1.110, 1.125, 1.170, 1.235, 1.305, 1.370, 1.440,
    1.445, 1.510, 1.590, 1.670, 1.755, 1.840, 1.855,
    1.930, 2.020, 2.130, 2.360, 2.372, 2.768, 3.059,
    3.381, 3.928, 4.129, 4.470, 4.670, 5.043, 5.623,
    6.160, 6.476, 7.079, 7.524, 7.943, 8.315, 8.913,
    9.190, 1.000e1, 1.068e1, 1.122e1, 1.259e1, 1.371e1, 1.523e1,
    1.674e1, 1.760e1, 1.903e1, 2.045e1, 2.260e1, 2.498e1, 2.792e1,
    2.920e1, 3.051e1, 3.389e1, 3.727e1, 3.981e1, 4.552e1, 4.785e1,
    5.012e1, 5.559e1, 6.144e1, 6.310e1, 6.790e1, 7.079e1, 7.889e1,
    8.528e1, 9.166e1, 1.013e2, 1.122e2, 1.301e2, 1.367e2, 1.585e2,
    1.670e2, 1.778e2, 2.040e2, 2.145e2, 2.430e2, 2.754e2, 3.043e2,
    3.536e2, 3.981e2, 4.540e2, 5.145e2, 5.830e2, 6.310e2, 6.773e2,
    7.079e2, 7.485e2, 8.482e2, 9.611e2, 1.010e3, 1.117e3, 1.234e3,
    1.364e3, 1.507e3, 1.585e3, 1.796e3, 2.035e3, 2.113e3, 2.249e3,
    2.371e3, 2.485e3, 2.613e3, 2.661e3, 2.747e3, 2.818e3, 3.035e3,
    3.162e3, 3.355e3, 3.548e3, 3.707e3, 3.981e3, 4.307e3, 4.643e3,
    5.004e3, 5.531e3, 6.267e3, 7.102e3, 7.466e3, 8.251e3, 9.119e3,
    1.008e4, 1.114e4, 1.171e4, 1.273e4, 1.383e4, 1.503e4, 1.585e4,
    1.662e4, 1.778e4, 1.931e4, 1.995e4, 2.054e4, 2.113e4, 2.187e4,
    2.239e4, 2.304e4, 2.358e4, 2.418e4, 2.441e4, 2.479e4, 2.512e4,
    2.585e4, 2.606e4, 2.661e4, 2.700e4, 2.738e4, 2.818e4, 2.850e4,
    2.901e4, 2.985e4, 3.073e4, 3.162e4, 3.183e4, 3.431e4, 3.698e4,
    4.087e4, 4.359e4, 4.631e4, 4.939e4, 5.248e4, 5.517e4, 5.656e4,
    6.173e4, 6.738e4, 7.200e4, 7.499e4, 7.950e4, 8.230e4, 8.250e4,
    8.652e4, 9.804e4, 1.111e5, 1.168e5, 1.228e5, 1.291e5, 1.357e5,
    1.426e5, 1.500e5, 1.576e5, 1.657e5, 1.742e5, 1.832e5, 1.925e5,
    2.024e5, 2.128e5, 2.237e5, 2.352e5, 2.472e5, 2.732e5, 2.873e5,
    2.945e5, 2.972e5, 2.985e5, 3.020e5, 3.337e5, 3.688e5, 3.877e5,
    4.076e5, 4.505e5, 5.234e5, 5.502e5, 5.784e5, 6.081e5, 6.393e5,
    6.721e5, 7.065e5, 7.427e5, 7.808e5, 8.209e5, 8.629e5, 9.072e5,
    9.616e5, 1.003e6, 1.108e6, 1.165e6, 1.225e6, 1.287e6, 1.353e6,
    1.423e6, 1.496e6, 1.572e6, 1.653e6, 1.738e6, 1.827e6, 1.921e6,
    2.019e6, 2.122e6, 2.231e6, 2.307e6, 2.346e6, 2.365e6, 2.385e6,
    2.466e6, 2.592e6, 2.725e6, 2.865e6, 3.012e6, 3.166e6, 3.329e6,
    3.679e6, 4.066e6, 4.493e6, 4.724e6, 4.966e6, 5.220e6, 5.488e6,
    5.769e6, 6.065e6, 6.376e6, 6.592e6, 6.703e6, 7.047e6, 7.408e6,
    7.788e6, 8.187e6, 8.607e6, 9.048e6, 9.512e6, 1.000e7, 1.051e7,
    1.105e7, 1.162e7, 1.221e7, 1.284e7, 1.350e7, 1.384e7, 1.419e7,
    1.455e7, 1.492e7, 1.568e7, 1.649e7, 1.691e7, 1.733e7, 1.964e7])
GROUP_STRUCTURES['SHEM-361'] = np.array([
    0.00000e+00, 2.49990e-03, 4.55602e-03, 7.14526e-03, 1.04505e-02,
    1.48300e-02, 2.00104e-02, 2.49394e-02, 2.92989e-02, 3.43998e-02,
    4.02999e-02, 4.73019e-02, 5.54982e-02, 6.51999e-02, 7.64969e-02,
    8.97968e-02, 1.04298e-01, 1.19995e-01, 1.37999e-01, 1.61895e-01,
    1.90005e-01, 2.09610e-01, 2.31192e-01, 2.54997e-01, 2.79989e-01,
    3.05012e-01, 3.25008e-01, 3.52994e-01, 3.90001e-01, 4.31579e-01,
    4.75017e-01, 5.20011e-01, 5.54990e-01, 5.94993e-01, 6.24999e-01,
    7.19999e-01, 8.00371e-01, 8.80024e-01, 9.19978e-01, 9.44022e-01,
    9.63960e-01, 9.81959e-01, 9.96501e-01, 1.00904e+00, 1.02101e+00,
    1.03499e+00, 1.07799e+00, 1.09198e+00, 1.10395e+00, 1.11605e+00,
    1.12997e+00, 1.14797e+00, 1.16999e+00, 1.21397e+00, 1.25094e+00,
    1.29304e+00, 1.33095e+00, 1.38098e+00, 1.41001e+00, 1.44397e+00,
    1.51998e+00, 1.58803e+00, 1.66895e+00, 1.77997e+00, 1.90008e+00,
    1.98992e+00, 2.07010e+00, 2.15695e+00, 2.21709e+00, 2.27299e+00,
    2.33006e+00, 2.46994e+00, 2.55000e+00, 2.59009e+00, 2.62005e+00,
    2.64004e+00, 2.70012e+00, 2.71990e+00, 2.74092e+00, 2.77512e+00,
    2.88405e+00, 3.14211e+00, 3.54307e+00, 3.71209e+00, 3.88217e+00,
    4.00000e+00, 4.21983e+00, 4.30981e+00, 4.41980e+00, 4.76785e+00,
    4.93323e+00, 5.10997e+00, 5.21008e+00, 5.32011e+00, 5.38003e+00,
    5.41025e+00, 5.48817e+00, 5.53004e+00, 5.61979e+00, 5.72015e+00,
    5.80021e+00, 5.96014e+00, 6.05991e+00, 6.16011e+00, 6.28016e+00,
    6.35978e+00, 6.43206e+00, 6.48178e+00, 6.51492e+00, 6.53907e+00,
    6.55609e+00, 6.57184e+00, 6.58829e+00, 6.60611e+00, 6.63126e+00,
    6.71668e+00, 6.74225e+00, 6.75981e+00, 6.77605e+00, 6.79165e+00,
    6.81070e+00, 6.83526e+00, 6.87021e+00, 6.91778e+00, 6.99429e+00,
    7.13987e+00, 7.38015e+00, 7.60035e+00, 7.73994e+00, 7.83965e+00,
    7.97008e+00, 8.13027e+00, 8.30032e+00, 8.52407e+00, 8.67369e+00,
    8.80038e+00, 8.97995e+00, 9.14031e+00, 9.50002e+00, 1.05793e+01,
    1.08038e+01, 1.10529e+01, 1.12694e+01, 1.15894e+01, 1.17094e+01,
    1.18153e+01, 1.19795e+01, 1.21302e+01, 1.23086e+01, 1.24721e+01,
    1.26000e+01, 1.33297e+01, 1.35460e+01, 1.40496e+01, 1.42505e+01,
    1.44702e+01, 1.45952e+01, 1.47301e+01, 1.48662e+01, 1.57792e+01,
    1.60498e+01, 1.65501e+01, 1.68305e+01, 1.74457e+01, 1.75648e+01,
    1.77590e+01, 1.79591e+01, 1.90848e+01, 1.91997e+01, 1.93927e+01,
    1.95974e+01, 2.00734e+01, 2.02751e+01, 2.04175e+01, 2.05199e+01,
    2.06021e+01, 2.06847e+01, 2.07676e+01, 2.09763e+01, 2.10604e+01,
    2.11448e+01, 2.12296e+01, 2.13360e+01, 2.14859e+01, 2.17018e+01,
    2.20011e+01, 2.21557e+01, 2.23788e+01, 2.25356e+01, 2.46578e+01,
    2.78852e+01, 3.16930e+01, 3.30855e+01, 3.45392e+01, 3.56980e+01,
    3.60568e+01, 3.64191e+01, 3.68588e+01, 3.73038e+01, 3.77919e+01,
    3.87874e+01, 3.97295e+01, 4.12270e+01, 4.21441e+01, 4.31246e+01,
    4.41721e+01, 4.52904e+01, 4.62053e+01, 4.75173e+01, 4.92591e+01,
    5.17847e+01, 5.29895e+01, 5.40600e+01, 5.70595e+01, 5.99250e+01,
    6.23083e+01, 6.36306e+01, 6.45923e+01, 6.50460e+01, 6.55029e+01,
    6.58312e+01, 6.61612e+01, 6.64929e+01, 6.68261e+01, 6.90682e+01,
    7.18869e+01, 7.35595e+01, 7.63322e+01, 7.93679e+01, 8.39393e+01,
    8.87741e+01, 9.33256e+01, 9.73287e+01, 1.00594e+02, 1.01098e+02,
    1.01605e+02, 1.02115e+02, 1.03038e+02, 1.05646e+02, 1.10288e+02,
    1.12854e+02, 1.15480e+02, 1.16524e+02, 1.17577e+02, 1.20554e+02,
    1.26229e+02, 1.32701e+02, 1.39504e+02, 1.46657e+02, 1.54176e+02,
    1.63056e+02, 1.67519e+02, 1.75229e+02, 1.83295e+02, 1.84952e+02,
    1.86251e+02, 1.87559e+02, 1.88877e+02, 1.90204e+02, 1.93078e+02,
    1.95996e+02, 2.00958e+02, 2.12108e+02, 2.24325e+02, 2.35590e+02,
    2.41796e+02, 2.56748e+02, 2.68297e+02, 2.76468e+02, 2.84888e+02,
    2.88327e+02, 2.95922e+02, 3.19928e+02, 3.35323e+02, 3.53575e+02,
    3.71703e+02, 3.90760e+02, 4.19094e+02, 4.53999e+02, 5.01746e+02,
    5.39204e+02, 5.77146e+02, 5.92941e+02, 6.00099e+02, 6.12834e+02,
    6.46837e+02, 6.77287e+02, 7.48517e+02, 8.32218e+02, 9.09681e+02,
    9.82494e+02, 1.06432e+03, 1.13467e+03, 1.34358e+03, 1.58620e+03,
    1.81183e+03, 2.08410e+03, 2.39729e+03, 2.70024e+03, 2.99618e+03,
    3.48107e+03, 4.09735e+03, 5.00451e+03, 6.11252e+03, 7.46585e+03,
    9.11881e+03, 1.11377e+04, 1.36037e+04, 1.48997e+04, 1.62005e+04,
    1.85847e+04, 2.26994e+04, 2.49991e+04, 2.61001e+04, 2.73944e+04,
    2.92810e+04, 3.34596e+04, 3.69786e+04, 4.08677e+04, 4.99159e+04,
    5.51656e+04, 6.73794e+04, 8.22974e+04, 9.46645e+04, 1.15624e+05,
    1.22773e+05, 1.40000e+05, 1.64999e+05, 1.95008e+05, 2.30014e+05,
    2.67826e+05, 3.20646e+05, 3.83884e+05, 4.12501e+05, 4.56021e+05,
    4.94002e+05, 5.78443e+05, 7.06511e+05, 8.60006e+05, 9.51119e+05,
    1.05115e+06, 1.16205e+06, 1.28696e+06, 1.33694e+06, 1.40577e+06,
    1.63654e+06, 1.90139e+06, 2.23130e+06, 2.72531e+06, 3.32871e+06,
    4.06569e+06, 4.96585e+06, 6.06530e+06, 6.70319e+06, 7.40817e+06,
    8.18730e+06, 9.04836e+06, 9.99999e+06, 1.16183e+07, 1.38403e+07,
    1.49182e+07, 1.96403e+07])
GROUP_STRUCTURES['CCFE-709'] = np.array([
    1.e-5, 1.0471e-5, 1.0965e-5, 1.1482e-5, 1.2023e-5,
    1.2589e-5, 1.3183e-5, 1.3804e-5, 1.4454e-5, 1.5136e-5,
    1.5849e-5, 1.6596e-5, 1.7378e-5, 1.8197e-5, 1.9055e-5,
    1.9953e-5, 2.0893e-5, 2.1878e-5, 2.2909e-5, 2.3988e-5,
    2.5119e-5, 2.6303e-5, 2.7542e-5, 2.8840e-5, 3.0200e-5,
    3.1623e-5, 3.3113e-5, 3.4674e-5, 3.6308e-5, 3.8019e-5,
    3.9811e-5, 4.1687e-5, 4.3652e-5, 4.5709e-5, 4.7863e-5,
    5.0119e-5, 5.2481e-5, 5.4954e-5, 5.7544e-5, 6.0256e-5,
    6.3096e-5, 6.6069e-5, 6.9183e-5, 7.2444e-5, 7.5858e-5,
    7.9433e-5, 8.3176e-5, 8.7096e-5, 9.1201e-5, 9.5499e-5,
    1.0000e-4, 1.0471e-4, 1.0965e-4, 1.1482e-4, 1.2023e-4,
    1.2589e-4, 1.3183e-4, 1.3804e-4, 1.4454e-4, 1.5136e-4,
    1.5849e-4, 1.6596e-4, 1.7378e-4, 1.8197e-4, 1.9055e-4,
    1.9953e-4, 2.0893e-4, 2.1878e-4, 2.2909e-4, 2.3988e-4,
    2.5119e-4, 2.6303e-4, 2.7542e-4, 2.8840e-4, 3.0200e-4,
    3.1623e-4, 3.3113e-4, 3.4674e-4, 3.6308e-4, 3.8019e-4,
    3.9811e-4, 4.1687e-4, 4.3652e-4, 4.5709e-4, 4.7863e-4,
    5.0119e-4, 5.2481e-4, 5.4954e-4, 5.7544e-4, 6.0256e-4,
    6.3096e-4, 6.6069e-4, 6.9183e-4, 7.2444e-4, 7.5858e-4,
    7.9433e-4, 8.3176e-4, 8.7096e-4, 9.1201e-4, 9.5499e-4,
    1.0000e-3, 1.0471e-3, 1.0965e-3, 1.1482e-3, 1.2023e-3,
    1.2589e-3, 1.3183e-3, 1.3804e-3, 1.4454e-3, 1.5136e-3,
    1.5849e-3, 1.6596e-3, 1.7378e-3, 1.8197e-3, 1.9055e-3,
    1.9953e-3, 2.0893e-3, 2.1878e-3, 2.2909e-3, 2.3988e-3,
    2.5119e-3, 2.6303e-3, 2.7542e-3, 2.8840e-3, 3.0200e-3,
    3.1623e-3, 3.3113e-3, 3.4674e-3, 3.6308e-3, 3.8019e-3,
    3.9811e-3, 4.1687e-3, 4.3652e-3, 4.5709e-3, 4.7863e-3,
    5.0119e-3, 5.2481e-3, 5.4954e-3, 5.7544e-3, 6.0256e-3,
    6.3096e-3, 6.6069e-3, 6.9183e-3, 7.2444e-3, 7.5858e-3,
    7.9433e-3, 8.3176e-3, 8.7096e-3, 9.1201e-3, 9.5499e-3,
    1.0000e-2, 1.0471e-2, 1.0965e-2, 1.1482e-2, 1.2023e-2,
    1.2589e-2, 1.3183e-2, 1.3804e-2, 1.4454e-2, 1.5136e-2,
    1.5849e-2, 1.6596e-2, 1.7378e-2, 1.8197e-2, 1.9055e-2,
    1.9953e-2, 2.0893e-2, 2.1878e-2, 2.2909e-2, 2.3988e-2,
    2.5119e-2, 2.6303e-2, 2.7542e-2, 2.8840e-2, 3.0200e-2,
    3.1623e-2, 3.3113e-2, 3.4674e-2, 3.6308e-2, 3.8019e-2,
    3.9811e-2, 4.1687e-2, 4.3652e-2, 4.5709e-2, 4.7863e-2,
    5.0119e-2, 5.2481e-2, 5.4954e-2, 5.7544e-2, 6.0256e-2,
    6.3096e-2, 6.6069e-2, 6.9183e-2, 7.2444e-2, 7.5858e-2,
    7.9433e-2, 8.3176e-2, 8.7096e-2, 9.1201e-2, 9.5499e-2,
    1.0000e-1, 1.0471e-1, 1.0965e-1, 1.1482e-1, 1.2023e-1,
    1.2589e-1, 1.3183e-1, 1.3804e-1, 1.4454e-1, 1.5136e-1,
    1.5849e-1, 1.6596e-1, 1.7378e-1, 1.8197e-1, 1.9055e-1,
    1.9953e-1, 2.0893e-1, 2.1878e-1, 2.2909e-1, 2.3988e-1,
    2.5119e-1, 2.6303e-1, 2.7542e-1, 2.8840e-1, 3.0200e-1,
    3.1623e-1, 3.3113e-1, 3.4674e-1, 3.6308e-1, 3.8019e-1,
    3.9811e-1, 4.1687e-1, 4.3652e-1, 4.5709e-1, 4.7863e-1,
    5.0119e-1, 5.2481e-1, 5.4954e-1, 5.7544e-1, 6.0256e-1,
    6.3096e-1, 6.6069e-1, 6.9183e-1, 7.2444e-1, 7.5858e-1,
    7.9433e-1, 8.3176e-1, 8.7096e-1, 9.1201e-1, 9.5499e-1,
    1.0000e0, 1.0471e0, 1.0965e0, 1.1482e0, 1.2023e0,
    1.2589e0, 1.3183e0, 1.3804e0, 1.4454e0, 1.5136e0,
    1.5849e0, 1.6596e0, 1.7378e0, 1.8197e0, 1.9055e0,
    1.9953e0, 2.0893e0, 2.1878e0, 2.2909e0, 2.3988e0,
    2.5119e0, 2.6303e0, 2.7542e0, 2.8840e0, 3.0200e0,
    3.1623e0, 3.3113e0, 3.4674e0, 3.6308e0, 3.8019e0,
    3.9811e0, 4.1687e0, 4.3652e0, 4.5709e0, 4.7863e0,
    5.0119e0, 5.2481e0, 5.4954e0, 5.7544e0, 6.0256e0,
    6.3096e0, 6.6069e0, 6.9183e0, 7.2444e0, 7.5858e0,
    7.9433e0, 8.3176e0, 8.7096e0, 9.1201e0, 9.5499e0,
    1.0000e1, 1.0471e1, 1.0965e1, 1.1482e1, 1.2023e1,
    1.2589e1, 1.3183e1, 1.3804e1, 1.4454e1, 1.5136e1,
    1.5849e1, 1.6596e1, 1.7378e1, 1.8197e1, 1.9055e1,
    1.9953e1, 2.0893e1, 2.1878e1, 2.2909e1, 2.3988e1,
    2.5119e1, 2.6303e1, 2.7542e1, 2.8840e1, 3.0200e1,
    3.1623e1, 3.3113e1, 3.4674e1, 3.6308e1, 3.8019e1,
    3.9811e1, 4.1687e1, 4.3652e1, 4.5709e1, 4.7863e1,
    5.0119e1, 5.2481e1, 5.4954e1, 5.7544e1, 6.0256e1,
    6.3096e1, 6.6069e1, 6.9183e1, 7.2444e1, 7.5858e1,
    7.9433e1, 8.3176e1, 8.7096e1, 9.1201e1, 9.5499e1,
    1.0000e2, 1.0471e2, 1.0965e2, 1.1482e2, 1.2023e2,
    1.2589e2, 1.3183e2, 1.3804e2, 1.4454e2, 1.5136e2,
    1.5849e2, 1.6596e2, 1.7378e2, 1.8197e2, 1.9055e2,
    1.9953e2, 2.0893e2, 2.1878e2, 2.2909e2, 2.3988e2,
    2.5119e2, 2.6303e2, 2.7542e2, 2.8840e2, 3.0200e2,
    3.1623e2, 3.3113e2, 3.4674e2, 3.6308e2, 3.8019e2,
    3.9811e2, 4.1687e2, 4.3652e2, 4.5709e2, 4.7863e2,
    5.0119e2, 5.2481e2, 5.4954e2, 5.7544e2, 6.0256e2,
    6.3096e2, 6.6069e2, 6.9183e2, 7.2444e2, 7.5858e2,
    7.9433e2, 8.3176e2, 8.7096e2, 9.1201e2, 9.5499e2,
    1.0000e3, 1.0471e3, 1.0965e3, 1.1482e3, 1.2023e3,
    1.2589e3, 1.3183e3, 1.3804e3, 1.4454e3, 1.5136e3,
    1.5849e3, 1.6596e3, 1.7378e3, 1.8197e3, 1.9055e3,
    1.9953e3, 2.0893e3, 2.1878e3, 2.2909e3, 2.3988e3,
    2.5119e3, 2.6303e3, 2.7542e3, 2.8840e3, 3.0200e3,
    3.1623e3, 3.3113e3, 3.4674e3, 3.6308e3, 3.8019e3,
    3.9811e3, 4.1687e3, 4.3652e3, 4.5709e3, 4.7863e3,
    5.0119e3, 5.2481e3, 5.4954e3, 5.7544e3, 6.0256e3,
    6.3096e3, 6.6069e3, 6.9183e3, 7.2444e3, 7.5858e3,
    7.9433e3, 8.3176e3, 8.7096e3, 9.1201e3, 9.5499e3,
    1.0000e4, 1.0471e4, 1.0965e4, 1.1482e4, 1.2023e4,
    1.2589e4, 1.3183e4, 1.3804e4, 1.4454e4, 1.5136e4,
    1.5849e4, 1.6596e4, 1.7378e4, 1.8197e4, 1.9055e4,
    1.9953e4, 2.0893e4, 2.1878e4, 2.2909e4, 2.3988e4,
    2.5119e4, 2.6303e4, 2.7542e4, 2.8840e4, 3.0200e4,
    3.1623e4, 3.3113e4, 3.4674e4, 3.6308e4, 3.8019e4,
    3.9811e4, 4.1687e4, 4.3652e4, 4.5709e4, 4.7863e4,
    5.0119e4, 5.2481e4, 5.4954e4, 5.7544e4, 6.0256e4,
    6.3096e4, 6.6069e4, 6.9183e4, 7.2444e4, 7.5858e4,
    7.9433e4, 8.3176e4, 8.7096e4, 9.1201e4, 9.5499e4,
    1.0000e5, 1.0471e5, 1.0965e5, 1.1482e5, 1.2023e5,
    1.2589e5, 1.3183e5, 1.3804e5, 1.4454e5, 1.5136e5,
    1.5849e5, 1.6596e5, 1.7378e5, 1.8197e5, 1.9055e5,
    1.9953e5, 2.0893e5, 2.1878e5, 2.2909e5, 2.3988e5,
    2.5119e5, 2.6303e5, 2.7542e5, 2.8840e5, 3.0200e5,
    3.1623e5, 3.3113e5, 3.4674e5, 3.6308e5, 3.8019e5,
    3.9811e5, 4.1687e5, 4.3652e5, 4.5709e5, 4.7863e5,
    5.0119e5, 5.2481e5, 5.4954e5, 5.7544e5, 6.0256e5,
    6.3096e5, 6.6069e5, 6.9183e5, 7.2444e5, 7.5858e5,
    7.9433e5, 8.3176e5, 8.7096e5, 9.1201e5, 9.5499e5,
    1.0000e6, 1.0471e6, 1.0965e6, 1.1482e6, 1.2023e6,
    1.2589e6, 1.3183e6, 1.3804e6, 1.4454e6, 1.5136e6,
    1.5849e6, 1.6596e6, 1.7378e6, 1.8197e6, 1.9055e6,
    1.9953e6, 2.0893e6, 2.1878e6, 2.2909e6, 2.3988e6,
    2.5119e6, 2.6303e6, 2.7542e6, 2.8840e6, 3.0200e6,
    3.1623e6, 3.3113e6, 3.4674e6, 3.6308e6, 3.8019e6,
    3.9811e6, 4.1687e6, 4.3652e6, 4.5709e6, 4.7863e6,
    5.0119e6, 5.2481e6, 5.4954e6, 5.7544e6, 6.0256e6,
    6.3096e6, 6.6069e6, 6.9183e6, 7.2444e6, 7.5858e6,
    7.9433e6, 8.3176e6, 8.7096e6, 9.1201e6, 9.5499e6,
    1.0000e7, 1.0200e7, 1.0400e7, 1.0600e7, 1.0800e7,
    1.1000e7, 1.1200e7, 1.1400e7, 1.1600e7, 1.1800e7,
    1.2000e7, 1.2200e7, 1.2400e7, 1.2600e7, 1.2800e7,
    1.3000e7, 1.3200e7, 1.3400e7, 1.3600e7, 1.3800e7,
    1.4000e7, 1.4200e7, 1.4400e7, 1.4600e7, 1.4800e7,
    1.5000e7, 1.5200e7, 1.5400e7, 1.5600e7, 1.5800e7,
    1.6000e7, 1.6200e7, 1.6400e7, 1.6600e7, 1.6800e7,
    1.7000e7, 1.7200e7, 1.7400e7, 1.7600e7, 1.7800e7,
    1.8000e7, 1.8200e7, 1.8400e7, 1.8600e7, 1.8800e7,
    1.9000e7, 1.9200e7, 1.9400e7, 1.9600e7, 1.9800e7,
    2.0000e7, 2.1000e7, 2.2000e7, 2.3000e7, 2.4000e7,
    2.5000e7, 2.6000e7, 2.7000e7, 2.8000e7, 2.9000e7,
    3.0000e7, 3.2000e7, 3.4000e7, 3.6000e7, 3.8000e7,
    4.0000e7, 4.2000e7, 4.4000e7, 4.6000e7, 4.8000e7,
    5.0000e7, 5.2000e7, 5.4000e7, 5.6000e7, 5.8000e7,
    6.0000e7, 6.5000e7, 7.0000e7, 7.5000e7, 8.0000e7,
    9.0000e7, 1.0000e8, 1.1000e8, 1.2000e8, 1.3000e8,
    1.4000e8, 1.5000e8, 1.6000e8, 1.8000e8, 2.0000e8,
    2.4000e8, 2.8000e8, 3.2000e8, 3.6000e8, 4.0000e8,
    4.4000e8, 4.8000e8, 5.2000e8, 5.6000e8, 6.0000e8,
    6.4000e8, 6.8000e8, 7.2000e8, 7.6000e8, 8.0000e8,
    8.4000e8, 8.8000e8, 9.2000e8, 9.6000e8, 1.0000e9,])
GROUP_STRUCTURES['UKAEA-1102'] = np.array([
    1.0000e-5, 1.0471e-5, 1.0965e-5, 1.1482e-5, 1.2023e-5,
    1.2589e-5, 1.3183e-5, 1.3804e-5, 1.4454e-5, 1.5136e-5,
    1.5849e-5, 1.6596e-5, 1.7378e-5, 1.8197e-5, 1.9055e-5,
    1.9953e-5, 2.0893e-5, 2.1878e-5, 2.2909e-5, 2.3988e-5,
    2.5119e-5, 2.6303e-5, 2.7542e-5, 2.8840e-5, 3.0200e-5,
    3.1623e-5, 3.3113e-5, 3.4674e-5, 3.6308e-5, 3.8019e-5,
    3.9811e-5, 4.1687e-5, 4.3652e-5, 4.5709e-5, 4.7863e-5,
    5.0119e-5, 5.2481e-5, 5.4954e-5, 5.7544e-5, 6.0256e-5,
    6.3096e-5, 6.6069e-5, 6.9183e-5, 7.2444e-5, 7.5858e-5,
    7.9433e-5, 8.3176e-5, 8.7096e-5, 9.1201e-5, 9.5499e-5,
    1.0000e-4, 1.0471e-4, 1.0965e-4, 1.1482e-4, 1.2023e-4,
    1.2589e-4, 1.3183e-4, 1.3804e-4, 1.4454e-4, 1.5136e-4,
    1.5849e-4, 1.6596e-4, 1.7378e-4, 1.8197e-4, 1.9055e-4,
    1.9953e-4, 2.0893e-4, 2.1878e-4, 2.2909e-4, 2.3988e-4,
    2.5119e-4, 2.6303e-4, 2.7542e-4, 2.8840e-4, 3.0200e-4,
    3.1623e-4, 3.3113e-4, 3.4674e-4, 3.6308e-4, 3.8019e-4,
    3.9811e-4, 4.1687e-4, 4.3652e-4, 4.5709e-4, 4.7863e-4,
    5.0119e-4, 5.2481e-4, 5.4954e-4, 5.7544e-4, 6.0256e-4,
    6.3096e-4, 6.6069e-4, 6.9183e-4, 7.2444e-4, 7.5858e-4,
    7.9433e-4, 8.3176e-4, 8.7096e-4, 9.1201e-4, 9.5499e-4,
    1.0000e-3, 1.0471e-3, 1.0965e-3, 1.1482e-3, 1.2023e-3,
    1.2589e-3, 1.3183e-3, 1.3804e-3, 1.4454e-3, 1.5136e-3,
    1.5849e-3, 1.6596e-3, 1.7378e-3, 1.8197e-3, 1.9055e-3,
    1.9953e-3, 2.0893e-3, 2.1878e-3, 2.2909e-3, 2.3988e-3,
    2.5119e-3, 2.6303e-3, 2.7542e-3, 2.8840e-3, 3.0200e-3,
    3.1623e-3, 3.3113e-3, 3.4674e-3, 3.6308e-3, 3.8019e-3,
    3.9811e-3, 4.1687e-3, 4.3652e-3, 4.5709e-3, 4.7863e-3,
    5.0119e-3, 5.2481e-3, 5.4954e-3, 5.7544e-3, 6.0256e-3,
    6.3096e-3, 6.6069e-3, 6.9183e-3, 7.2444e-3, 7.5858e-3,
    7.9433e-3, 8.3176e-3, 8.7096e-3, 9.1201e-3, 9.5499e-3,
    1.0000e-2, 1.0471e-2, 1.0965e-2, 1.1482e-2, 1.2023e-2,
    1.2589e-2, 1.3183e-2, 1.3804e-2, 1.4454e-2, 1.5136e-2,
    1.5849e-2, 1.6596e-2, 1.7378e-2, 1.8197e-2, 1.9055e-2,
    1.9953e-2, 2.0893e-2, 2.1878e-2, 2.2909e-2, 2.3988e-2,
    2.5119e-2, 2.6303e-2, 2.7542e-2, 2.8840e-2, 3.0200e-2,
    3.1623e-2, 3.3113e-2, 3.4674e-2, 3.6308e-2, 3.8019e-2,
    3.9811e-2, 4.1687e-2, 4.3652e-2, 4.5709e-2, 4.7863e-2,
    5.0119e-2, 5.2481e-2, 5.4954e-2, 5.7544e-2, 6.0256e-2,
    6.3096e-2, 6.6069e-2, 6.9183e-2, 7.2444e-2, 7.5858e-2,
    7.9433e-2, 8.3176e-2, 8.7096e-2, 9.1201e-2, 9.5499e-2,
    1.0000e-1, 1.0471e-1, 1.0965e-1, 1.1482e-1, 1.2023e-1,
    1.2589e-1, 1.3183e-1, 1.3804e-1, 1.4454e-1, 1.5136e-1,
    1.5849e-1, 1.6596e-1, 1.7378e-1, 1.8197e-1, 1.9055e-1,
    1.9953e-1, 2.0893e-1, 2.1878e-1, 2.2909e-1, 2.3988e-1,
    2.5119e-1, 2.6303e-1, 2.7542e-1, 2.8840e-1, 3.0200e-1,
    3.1623e-1, 3.3113e-1, 3.4674e-1, 3.6308e-1, 3.8019e-1,
    3.9811e-1, 4.1687e-1, 4.3652e-1, 4.5709e-1, 4.7863e-1,
    5.0119e-1, 5.2481e-1, 5.5000e-1, 5.7500e-1, 6.0000e-1,
    6.2500e-1, 6.5000e-1, 6.7500e-1, 7.0000e-1, 7.2500e-1,
    7.5000e-1, 7.7500e-1, 8.0000e-1, 8.2500e-1, 8.5000e-1,
    8.7500e-1, 9.0000e-1, 9.2500e-1, 9.5000e-1, 9.7500e-1,
    1.0000, 1.0250, 1.0500, 1.0750, 1.1000,
    1.1250, 1.1500, 1.1750, 1.2000, 1.2250,
    1.2500, 1.2750, 1.3000, 1.3250, 1.3500,
    1.3750, 1.4000, 1.4250, 1.4500, 1.4750,
    1.5000, 1.5250, 1.5500, 1.5750, 1.6000,
    1.6250, 1.6500, 1.6750, 1.7000, 1.7250,
    1.7500, 1.7750, 1.8000, 1.8250, 1.8500,
    1.8750, 1.9000, 1.9250, 1.9500, 1.9750,
    2.0000, 2.0250, 2.0500, 2.0750, 2.1000,
    2.1250, 2.1500, 2.1750, 2.2000, 2.2250,
    2.2500, 2.2750, 2.3000, 2.3250, 2.3500,
    2.3750, 2.4000, 2.4250, 2.4500, 2.4750,
    2.5000, 2.5250, 2.5500, 2.5750, 2.6000,
    2.6250, 2.6500, 2.6750, 2.7000, 2.7250,
    2.7500, 2.7750, 2.8000, 2.8250, 2.8500,
    2.8750, 2.9000, 2.9250, 2.9500, 2.9750,
    3.0000, 3.0250, 3.0500, 3.0750, 3.1000,
    3.1250, 3.1500, 3.1750, 3.2000, 3.2250,
    3.2500, 3.2750, 3.3000, 3.3250, 3.3500,
    3.3750, 3.4000, 3.4250, 3.4500, 3.4750,
    3.5000, 3.5250, 3.5500, 3.5750, 3.6000,
    3.6250, 3.6500, 3.6750, 3.7000, 3.7250,
    3.7500, 3.7750, 3.8000, 3.8250, 3.8500,
    3.8750, 3.9000, 3.9250, 3.9500, 3.9750,
    4.0000, 4.0250, 4.0500, 4.0750, 4.1000,
    4.1250, 4.1500, 4.1750, 4.2000, 4.2250,
    4.2500, 4.2750, 4.3000, 4.3250, 4.3500,
    4.3750, 4.4000, 4.4250, 4.4500, 4.4750,
    4.5000, 4.5250, 4.5500, 4.5750, 4.6000,
    4.6250, 4.6500, 4.6750, 4.7000, 4.7250,
    4.7500, 4.7750, 4.8000, 4.8250, 4.8500,
    4.8750, 4.9000, 4.9250, 4.9500, 4.9750,
    5.0000, 5.0250, 5.0500, 5.0750, 5.1000,
    5.1250, 5.1500, 5.1750, 5.2000, 5.2250,
    5.2500, 5.2750, 5.3000, 5.3250, 5.3500,
    5.3750, 5.4000, 5.4250, 5.4500, 5.4750,
    5.5000, 5.5250, 5.5500, 5.5750, 5.6000,
    5.6250, 5.6500, 5.6750, 5.7000, 5.7250,
    5.7500, 5.7750, 5.8000, 5.8250, 5.8500,
    5.8750, 5.9000, 5.9250, 5.9500, 5.9750,
    6.0000, 6.0250, 6.0500, 6.0750, 6.1000,
    6.1250, 6.1500, 6.1750, 6.2000, 6.2250,
    6.2500, 6.2750, 6.3000, 6.3250, 6.3500,
    6.3750, 6.4000, 6.4250, 6.4500, 6.4750,
    6.5000, 6.5250, 6.5500, 6.5750, 6.6000,
    6.6250, 6.6500, 6.6750, 6.7000, 6.7250,
    6.7500, 6.7750, 6.8000, 6.8250, 6.8500,
    6.8750, 6.9000, 6.9250, 6.9500, 6.9750,
    7.0000, 7.0250, 7.0500, 7.0750, 7.1000,
    7.1250, 7.1500, 7.1750, 7.2000, 7.2250,
    7.2500, 7.2750, 7.3000, 7.3250, 7.3500,
    7.3750, 7.4000, 7.4250, 7.4500, 7.4750,
    7.5000, 7.5250, 7.5500, 7.5750, 7.6000,
    7.6250, 7.6500, 7.6750, 7.7000, 7.7250,
    7.7500, 7.7750, 7.8000, 7.8250, 7.8500,
    7.8750, 7.9000, 7.9250, 7.9500, 7.9750,
    8.0000, 8.0250, 8.0500, 8.0750, 8.1000,
    8.1250, 8.1500, 8.1750, 8.2000, 8.2250,
    8.2500, 8.2750, 8.3000, 8.3250, 8.3500,
    8.3750, 8.4000, 8.4250, 8.4500, 8.4750,
    8.5000, 8.5250, 8.5500, 8.5750, 8.6000,
    8.6250, 8.6500, 8.6750, 8.7000, 8.7250,
    8.7500, 8.7750, 8.8000, 8.8250, 8.8500,
    8.8750, 8.9000, 8.9250, 8.9500, 8.9750,
    9.0000, 9.0250, 9.0500, 9.0750, 9.1000,
    9.1250, 9.1500, 9.1750, 9.2000, 9.2250,
    9.2500, 9.2750, 9.3000, 9.3250, 9.3500,
    9.3750, 9.4000, 9.4250, 9.4500, 9.4750,
    9.5000, 9.5250, 9.5500, 9.5750, 9.6000,
    9.6250, 9.6500, 9.6750, 9.7000, 9.7250,
    9.7500, 9.7750, 9.8000, 9.8250, 9.8500,
    9.8750, 9.9000, 9.9250, 9.9500, 9.9750,
    1.0000e1, 1.0471e1, 1.0965e1, 1.1482e1, 1.2023e1,
    1.2589e1, 1.3183e1, 1.3804e1, 1.4454e1, 1.5136e1,
    1.5849e1, 1.6596e1, 1.7378e1, 1.8197e1, 1.9055e1,
    1.9953e1, 2.0893e1, 2.1878e1, 2.2909e1, 2.3988e1,
    2.5119e1, 2.6303e1, 2.7542e1, 2.8840e1, 3.0200e1,
    3.1623e1, 3.3113e1, 3.4674e1, 3.6308e1, 3.8019e1,
    3.9811e1, 4.1687e1, 4.3652e1, 4.5709e1, 4.7863e1,
    5.0119e1, 5.2481e1, 5.4954e1, 5.7544e1, 6.0256e1,
    6.3096e1, 6.6069e1, 6.9183e1, 7.2444e1, 7.5858e1,
    7.9433e1, 8.3176e1, 8.7096e1, 9.1201e1, 9.5499e1,
    1.0000e2, 1.0471e2, 1.0965e2, 1.1482e2, 1.2023e2,
    1.2589e2, 1.3183e2, 1.3804e2, 1.4454e2, 1.5136e2,
    1.5849e2, 1.6596e2, 1.7378e2, 1.8197e2, 1.9055e2,
    1.9953e2, 2.0893e2, 2.1878e2, 2.2909e2, 2.3988e2,
    2.5119e2, 2.6303e2, 2.7542e2, 2.8840e2, 3.0200e2,
    3.1623e2, 3.3113e2, 3.4674e2, 3.6308e2, 3.8019e2,
    3.9811e2, 4.1687e2, 4.3652e2, 4.5709e2, 4.7863e2,
    5.0119e2, 5.2481e2, 5.4954e2, 5.7544e2, 6.0256e2,
    6.3096e2, 6.6069e2, 6.9183e2, 7.2444e2, 7.5858e2,
    7.9433e2, 8.3176e2, 8.7096e2, 9.1201e2, 9.5499e2,
    1.0000e3, 1.0471e3, 1.0965e3, 1.1482e3, 1.2023e3,
    1.2589e3, 1.3183e3, 1.3804e3, 1.4454e3, 1.5136e3,
    1.5849e3, 1.6596e3, 1.7378e3, 1.8197e3, 1.9055e3,
    1.9953e3, 2.0893e3, 2.1878e3, 2.2909e3, 2.3988e3,
    2.5119e3, 2.6303e3, 2.7542e3, 2.8840e3, 3.0200e3,
    3.1623e3, 3.3113e3, 3.4674e3, 3.6308e3, 3.8019e3,
    3.9811e3, 4.1687e3, 4.3652e3, 4.5709e3, 4.7863e3,
    5.0119e3, 5.2481e3, 5.4954e3, 5.7544e3, 6.0256e3,
    6.3096e3, 6.6069e3, 6.9183e3, 7.2444e3, 7.5858e3,
    7.9433e3, 8.3176e3, 8.7096e3, 9.1201e3, 9.5499e3,
    1.0000e4, 1.0471e4, 1.0965e4, 1.1482e4, 1.2023e4,
    1.2589e4, 1.3183e4, 1.3804e4, 1.4454e4, 1.5136e4,
    1.5849e4, 1.6596e4, 1.7378e4, 1.8197e4, 1.9055e4,
    1.9953e4, 2.0893e4, 2.1878e4, 2.2909e4, 2.3988e4,
    2.5119e4, 2.6303e4, 2.7542e4, 2.8840e4, 3.0200e4,
    3.1623e4, 3.3113e4, 3.4674e4, 3.6308e4, 3.8019e4,
    3.9811e4, 4.1687e4, 4.3652e4, 4.5709e4, 4.7863e4,
    5.0119e4, 5.2481e4, 5.4954e4, 5.7544e4, 6.0256e4,
    6.3096e4, 6.6069e4, 6.9183e4, 7.2444e4, 7.5858e4,
    7.9433e4, 8.3176e4, 8.7096e4, 9.1201e4, 9.5499e4,
    1.0000e5, 1.0471e5, 1.0965e5, 1.1482e5, 1.2023e5,
    1.2589e5, 1.3183e5, 1.3804e5, 1.4454e5, 1.5136e5,
    1.5849e5, 1.6596e5, 1.7378e5, 1.8197e5, 1.9055e5,
    1.9953e5, 2.0893e5, 2.1878e5, 2.2909e5, 2.3988e5,
    2.5119e5, 2.6303e5, 2.7542e5, 2.8840e5, 3.0200e5,
    3.1623e5, 3.3113e5, 3.4674e5, 3.6308e5, 3.8019e5,
    3.9811e5, 4.1687e5, 4.3652e5, 4.5709e5, 4.7863e5,
    5.0119e5, 5.2481e5, 5.4954e5, 5.7544e5, 6.0256e5,
    6.3096e5, 6.6069e5, 6.9183e5, 7.2444e5, 7.5858e5,
    7.9433e5, 8.3176e5, 8.7096e5, 9.1201e5, 9.5499e5,
    1.0000e6, 1.0471e6, 1.0965e6, 1.1482e6, 1.2023e6,
    1.2589e6, 1.3183e6, 1.3804e6, 1.4454e6, 1.5136e6,
    1.5849e6, 1.6596e6, 1.7378e6, 1.8197e6, 1.9055e6,
    1.9953e6, 2.0893e6, 2.1878e6, 2.2909e6, 2.3988e6,
    2.5119e6, 2.6303e6, 2.7542e6, 2.8840e6, 3.0200e6,
    3.1623e6, 3.3113e6, 3.4674e6, 3.6308e6, 3.8019e6,
    3.9811e6, 4.1687e6, 4.3652e6, 4.5709e6, 4.7863e6,
    5.0000e6, 5.2000e6, 5.4000e6, 5.6000e6, 5.8000e6,
    6.0000e6, 6.2000e6, 6.4000e6, 6.6000e6, 6.8000e6,
    7.0000e6, 7.2000e6, 7.4000e6, 7.6000e6, 7.8000e6,
    8.0000e6, 8.2000e6, 8.4000e6, 8.6000e6, 8.8000e6,
    9.0000e6, 9.2000e6, 9.4000e6, 9.6000e6, 9.8000e6,
    1.0000e7, 1.0200e7, 1.0400e7, 1.0600e7, 1.0800e7,
    1.1000e7, 1.1200e7, 1.1400e7, 1.1600e7, 1.1800e7,
    1.2000e7, 1.2200e7, 1.2400e7, 1.2600e7, 1.2800e7,
    1.3000e7, 1.3200e7, 1.3400e7, 1.3600e7, 1.3800e7,
    1.4000e7, 1.4200e7, 1.4400e7, 1.4600e7, 1.4800e7,
    1.5000e7, 1.5200e7, 1.5400e7, 1.5600e7, 1.5800e7,
    1.6000e7, 1.6200e7, 1.6400e7, 1.6600e7, 1.6800e7,
    1.7000e7, 1.7200e7, 1.7400e7, 1.7600e7, 1.7800e7,
    1.8000e7, 1.8200e7, 1.8400e7, 1.8600e7, 1.8800e7,
    1.9000e7, 1.9200e7, 1.9400e7, 1.9600e7, 1.9800e7,
    2.0000e7, 2.0200e7, 2.0400e7, 2.0600e7, 2.0800e7,
    2.1000e7, 2.1200e7, 2.1400e7, 2.1600e7, 2.1800e7,
    2.2000e7, 2.2200e7, 2.2400e7, 2.2600e7, 2.2800e7,
    2.3000e7, 2.3200e7, 2.3400e7, 2.3600e7, 2.3800e7,
    2.4000e7, 2.4200e7, 2.4400e7, 2.4600e7, 2.4800e7,
    2.5000e7, 2.5200e7, 2.5400e7, 2.5600e7, 2.5800e7,
    2.6000e7, 2.6200e7, 2.6400e7, 2.6600e7, 2.6800e7,
    2.7000e7, 2.7200e7, 2.7400e7, 2.7600e7, 2.7800e7,
    2.8000e7, 2.8200e7, 2.8400e7, 2.8600e7, 2.8800e7,
    2.9000e7, 2.9200e7, 2.9400e7, 2.9600e7, 2.9800e7,
    3.0000e7, 3.0200e7, 3.1623e7, 3.3113e7, 3.4674e7,
    3.6308e7, 3.8019e7, 3.9811e7, 4.1687e7, 4.3652e7,
    4.5709e7, 4.7863e7, 5.0119e7, 5.2481e7, 5.4954e7,
    5.7544e7, 6.0256e7, 6.3096e7, 6.6069e7, 6.9183e7,
    7.2444e7, 7.5858e7, 7.9433e7, 8.3176e7, 8.7096e7,
    9.1201e7, 9.5499e7, 1.0000e8, 1.0471e8, 1.0965e8,
    1.1482e8, 1.2023e8, 1.2589e8, 1.3183e8, 1.3804e8,
    1.4454e8, 1.5136e8, 1.5849e8, 1.6596e8, 1.7378e8,
    1.8197e8, 1.9055e8, 1.9953e8, 2.0893e8, 2.1878e8,
    2.2909e8, 2.3988e8, 2.5119e8, 2.6303e8, 2.7542e8,
    2.8840e8, 3.0200e8, 3.1623e8, 3.3113e8, 3.4674e8,
    3.6308e8, 3.8019e8, 3.9811e8, 4.1687e8, 4.3652e8,
    4.5709e8, 4.7863e8, 5.0119e8, 5.2481e8, 5.4954e8,
    5.7544e8, 6.0256e8, 6.3096e8, 6.6069e8, 6.9183e8,
    7.2444e8, 7.5858e8, 7.9433e8, 8.3176e8, 8.7096e8,
    9.1201e8, 9.5499e8, 1.e9])