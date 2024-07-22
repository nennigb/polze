# -*- coding: utf-8 -*-
# This file is part of polze, a package to locate poles and zeros of a
# meromorphic function with their multiplicities.
# Copyright (C) 2019  Benoit Nennig, benoit.nennig@isae-supmeca.fr

# polze is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# polze is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with polze. If not, see <https://www.gnu.org/licenses/>.
""" This module contains more involved test cases.
"""
import unittest
import doctest
import polze
import numpy as np
import cmath
import sys

# Numpy 2.0 change default printing options making doctest failing.
# https://numpy.org/neps/nep-0051-scalar-representation.html
# Use legacy mode for testing
if np.lib.NumpyVersion(np.__version__) >= '2.0.0b1':
    np.set_printoptions(legacy="1.25")

# Define a non vectorized function
# For parallel processing we need a function defined outside the class.
def f(z):
    return cmath.tan(z)           # define the function


def df(z):
    return cmath.tan(z)**2 + 1    # and its derivative [optional]

# Define a vectorized function
# For parallel processing we need a function defined outside the class.
def f_np(z):
    return np.tan(z)

def df_np(z):
    return np.tan(z)**2 + 1

class TestBasic(unittest.TestCase):
    """ Test suite.
    """

    def test_non_vectorized(self):
        """ Test with non vectorized function.
        """
        # Reuse a non vectorized function
        pz = polze.PZ((f, df), Rmax=5, Npz=10, Ni=1024,
                      options={'_vectorized': False})
        _ = pz.solve()
        p, z = pz.dispatch()

        p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
        z_ref = np.array([-1, 0, 1])*(np.pi)
        p.sort()
        self.assertTrue(np.linalg.norm(p - p_ref) < 1e-12)

        z.sort()
        self.assertTrue(np.linalg.norm(z - z_ref) < 1e-12)

    def test_iterative_ref(self):
        """ Test iterative refinement for simple roots.
        """
        tol = 1e-12
        # Define a vectorized function
        f = lambda z: np.tan(z)            # define the function
        df = lambda z: np.tan(z)**2 + 1    # and its derivative [optional]
        # Low quality initial solution (small Ni)
        pz = polze.PZ((f, df), Rmax=5, Npz=10, Ni=100,
                      options={'_vectorized': True,
                               '_tol': 1e-3,
                               '_NR_tol': tol})
        # Get poles and zeros from initial solution
        _ = pz.solve()
        p0, z0 = pz.dispatch()
        p0.sort()
        z0.sort()
        # Refine step
        pz_ref, mu_ref, info = pz.iterative_ref()
        p, z = pz._dispatch(pz_ref, mu_ref)
        z.sort()
        p.sort()
        # Check reference solutions
        p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
        z_ref = np.array([-1, 0, 1])*(np.pi)
        self.assertTrue(np.linalg.norm(p - p_ref) < tol)
        self.assertTrue(np.linalg.norm(z - z_ref) < tol)
        # check relative convergence
        self.assertTrue(info['ErrTot'] < tol)
        # initial sols are not accurate enougth
        self.assertFalse(np.linalg.norm(z0 - z_ref) < tol)

    def test_iterative_ref_no_df(self):
        """ Test iterative refinement for simple roots without df.
        """
        tol = 1e-10
        # Define a vectorized function
        f = lambda z: np.tan(z)            # define the function
        # Low quality initial solution (small Ni)
        # Here we prescribe a big '_clean_tol' to avoid spurious root that
        # break the number of poles in the test
        pz = polze.PZ(f, Rmax=5, Npz=10, Ni=100,
                      options={'_vectorized': True,
                               '_tol': 1e-3,
                               '_NR_tol': tol,
                               '_clean_tol': 1e-2})
        # Get poles and zeros from initial solution
        _ = pz.solve()
        p0, z0 = pz.dispatch()
        p0.sort()
        z0.sort()
        # Refine step
        pz_ref, mu_ref, info = pz.iterative_ref()
        p, z = pz._dispatch(pz_ref, mu_ref)
        z.sort()
        p.sort()
        # Check reference solutions
        p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
        z_ref = np.array([-1, 0, 1])*(np.pi)
        self.assertTrue(np.linalg.norm(p - p_ref) < tol)
        self.assertTrue(np.linalg.norm(z - z_ref) < tol)
        # check relative convergence
        self.assertTrue(info['ErrTot'] < tol)
        # initial sols are not accurate enougth
        self.assertFalse(np.linalg.norm(z0 - z_ref) < tol)
        
    def test_contour_ref(self):
        """ Test contour refinement for simple roots.
        """
        tol = 1e-12
        # Define a vectorized function
        f = lambda z: np.tan(z)            # define the function
        df = lambda z: np.tan(z)**2 + 1    # and its derivative [optional]
        # Low quality initial solution (small Ni)
        pz = polze.PZ((f, df), Rmax=5, Npz=10, Ni=100,
                      options={'_vectorized': True,
                               '_tol': 1e-3,
                               '_NR_tol': tol})
        # Get poles and zeros from initial solution
        _ = pz.solve()
        p0, z0 = pz.dispatch()
        p0.sort()
        z0.sort()
        # Refine step
        pz_ref, mu_ref, info = pz.contour_ref(Ni=1000)
        p, z = pz._dispatch(pz_ref, mu_ref)
        z.sort()
        p.sort()
        # Check reference solutions
        p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
        z_ref = np.array([-1, 0, 1])*(np.pi)
        self.assertTrue(np.linalg.norm(p - p_ref) < tol)
        self.assertTrue(np.linalg.norm(z - z_ref) < tol)
        # initial sols are not accurate enougth
        self.assertFalse(np.linalg.norm(z0 - z_ref) < tol)

    def test_with_finite_difference_on_circle(self):
        """ Test using the finite difference on the circle for a
        non-vectorized function.
        """
        # Here tol bigger than when the derivative are given.
        test_tol = 1e-9
        # define a non vectorized function
        f = lambda z: cmath.tan(z)

        pz = polze.PZ(f, Rmax=5, Npz=10, Ni=1024,
                      options={'_vectorized': False})
        _ = pz.solve()
        p, z = pz.dispatch()

        p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
        z_ref = np.array([-1, 0, 1])*(np.pi)
        p = np.sort_complex(p)
        self.assertTrue(np.linalg.norm(p - p_ref) < test_tol)

        z = np.sort_complex(z)
        self.assertTrue(np.linalg.norm(z - z_ref) < test_tol)

    def test_zeros_only(self):
        """ Test with function without poles.
        """
        # define a non vectorized function
        f = lambda z: np.sin(z)            # define the function
        df = lambda z: np.cos(z)          # and its derivative [optional]
        # Npz is too small and will be corrected by `estimate_Npz` method
        pz = polze.PZ((f, df), Rmax=5, Npz=1, Ni=2048,
                      options={'_zeros_only': True})
        _ = pz.solve()
        p, z = pz.dispatch()

        z_ref = np.array([-1, 0, 1])*(np.pi)
        z.sort()
        self.assertTrue(np.linalg.norm(z - z_ref) < 1e-12)

        z.sort()
        self.assertTrue(np.linalg.norm(z - z_ref) < 1e-12)

    def test_rational_fraction(self):
        """ Test with a rational fraction.
        """
        tol = 1e-10
        # Define the rational fraction
        pol = np.polynomial.polynomial
        z_ref = np.array([1. + 0.2j, 2. - 0.3j, -1 + 0.4j, 7.])
        p_ref = np.array([.5 + 0.2j, 3 + 0.3j, -7, 5.])

        def f(z):
            """ Define the rational fraction from roots.
            """
            n = pol.polyval(z, pol.polyfromroots(z_ref))
            d = pol.polyval(z, pol.polyfromroots(p_ref))
            return n/d

        pz = polze.PZ(f, Rmax=8, Npz=10, Ni=2048)
        _ = pz.solve()
        p, z = pz.dispatch()
        z.sort()
        p.sort()
        self.assertTrue(np.linalg.norm(z - np.sort(z_ref)) < tol)
        self.assertTrue(np.linalg.norm(p - np.sort(p_ref)) < tol)

    def test_rational_fraction_multiple_roots(self):
        """ Test with a rational fraction with multiple roots (poles and zeros).
        """
        tol = 1e-10
        # Define the rational fraction
        pol = np.polynomial.polynomial
        z_ref = np.array([1. + 0.2j, 1. + 0.2j, -1 + 0.4j, 7.])
        p_ref = np.array([.5 + 0.2j, -7, -7, -7])

        def f(z):
            """ Define the rational fraction from roots.
            """
            n = pol.polyval(z, pol.polyfromroots(z_ref))
            d = pol.polyval(z, pol.polyfromroots(p_ref))
            return n/d

        pz = polze.PZ(f, Rmax=8, Npz=10, Ni=2048)
        _ = pz.solve()
        (p, pm), (z, zm) = pz.dispatch(multiplicities=True)
        z.sort()
        p.sort()
        self.assertTrue(np.linalg.norm(z - np.sort(np.unique(z_ref))) < tol)
        self.assertTrue(np.linalg.norm(p - np.sort(np.unique(p_ref))) < tol)
        # check their is pole-3
        self.assertAlmostEqual(pm.min(), -3, places=4)
        # check their is zeros-2
        self.assertAlmostEqual(zm.max(), 2, places=4)

    def test_parallel(self):
        """ Test accuracy with parallel execution.

        In this example, each task is very short, thus no speed up
        is expected here. Just check if the results are good.
        """
        # Reuse functions defined at top level
        pz = polze.PZ((f, df), Rmax=5, Npz=10, Ni=1000,
                      options={'_vectorized': False,
                               '_parallel': True, '_max_workers': 4})
        _ = pz.solve()
        p, z = pz.dispatch()

        p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
        z_ref = np.array([-1, 0, 1])*(np.pi)
        p.sort()
        self.assertTrue(np.linalg.norm(p - p_ref) < 1e-12)

        z.sort()
        self.assertTrue(np.linalg.norm(z - z_ref) < 1e-12)


    def test_parallel_vectorized(self):
        """ Test accuracy with parallel execution with vectorized (numpy) func.

        In this example, each task is very short, thus no speed up
        is expected here. Just check if the results are good.
        """
        # Reuse functions defined at top level
        pz = polze.PZ((f_np, df_np), Rmax=5, Npz=10, Ni=1000,
                      options={'_vectorized': True,
                               '_parallel': True, '_max_workers': 4})
        _ = pz.solve()
        p, z = pz.dispatch()

        p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
        z_ref = np.array([-1, 0, 1])*(np.pi)
        p.sort()
        self.assertTrue(np.linalg.norm(p - p_ref) < 1e-12)

        z.sort()
        self.assertTrue(np.linalg.norm(z - z_ref) < 1e-12)


if __name__ == '__main__':
    # Run unittest test suite
    print('> Running tests...')
    loadTestsFromTestCase = unittest.defaultTestLoader.loadTestsFromTestCase
    suite = loadTestsFromTestCase(TestBasic)
    # Add doctest from _polze module
    suite.addTest(doctest.DocTestSuite(polze._polze,
                                       optionflags=(doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)))
    # define the runner
    runner = unittest.TextTestRunner(verbosity=3)
    # run all the suite
    result = runner.run(suite)
    # runner doesn't change exit status
    if result.wasSuccessful():
        sys.exit(0)
    else:
        sys.exit(1)
