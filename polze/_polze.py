#!/usr/bin/env python3
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
# along with polze.  If not, see <https://www.gnu.org/licenses/>.

r"""
# Tutorials

_This tutorial is part of `polze._polze.py` module doctest suite. The following
examples illustrate the main functionnalities of the API._

Let us consider the `tan` function that contains *poles* and *zeros*
>>> from polze import PZ
>>> import numpy as np
>>> f = lambda z: np.tan(z)            # define the function
>>> df = lambda z: np.tan(z)**2 + 1    # and its derivative [optional]

Initiatlize the solver
>>> pz = PZ((f, df), Rmax=5, Npz=10, Ni=1024)

Compute a 1st solution
>>> _ = pz.solve()

Separate zeros and poles using multiplicities' sign
>>> p, z = pz.dispatch()

If `df` is not provided, the derivative is computed numerically using high order
finite difference scheme.

Check error with analytic solution
>>> p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
>>> z_ref = np.array([-1, 0, 1])*(np.pi)
>>> p.sort(); np.linalg.norm(p - p_ref) < 1e-12
True
>>> z.sort(); np.linalg.norm(z- z_ref) < 1e-12
True

To improve accuracy, optional iterative refinement step can be done with
`pz.iterative_ref()`. This method use a modified Newton-Raphson scheme
to take into account the roots multiplicity. See also `pz.contour_ref()`.

To be sure that all the roots are _genuine_ poles or zeros, the `check` method
can be used to eval the function.
The default ouput of `dispatch` remove root with multiplicity
close to 0. The `check` method allows also to inspect all roots.
>>> status = pz.check()
>>> status.all_roots
True

Another possibility is to solve in a circle centered on R0 to solve in a
specific area of the complex plane :
>>> pz = PZ((f, df), Rmax=5, Npz=10, Ni=1500, R0=10)
>>> _ = pz.solve()
>>> p, z = pz.dispatch()
>>> p_ref = np.array([2.5, 3.5, 4.5])*np.pi
>>> z_ref = np.array([2, 3, 4])*np.pi
>>> p.sort(); np.linalg.norm(p - p_ref) < 1e-10
True
>>> z.sort(); np.linalg.norm(z - z_ref) < 1e-10
True

The contour may be split into several annular regions if the upper bound for
the number of roots `Npz` becomes bigger than the `_Npz_limit` default value.
It could be changed as all preset options thank to the `options` dict
>>> pz = PZ((f, df), Rmax=40.8, Npz=60, Ni=3000, R0=10, options={'_Npz_limit': 8})

Here, the circle is split into 8 (`Npz/_Npz_limit`) annular regions. Each
radius can be found in `pz._Ri` attribut
>>> _ = pz.solve()  # doctest: +ELLIPSIS
>>> p, z = pz.dispatch()
>>> p_ref = np.arange(-9.5, 16.5)*np.pi
>>> p.sort(); np.linalg.norm(p - p_ref) < 1e-3
True

The other way is to used `split` argument to split into a prescribed number of
annular regions
>>> pz = PZ((f, df), Rmax=40, Npz=60, Ni=3000, R0=10, split=16)
>>> len(pz._Ri)
16

If poles or zeros lies in the contour or the subcontours, `polze` can detect it
and shift the radius Ri by multiply it by `_RiShift`. For instance
>>> pz = PZ((f, df), Rmax=2*np.pi, Npz=20, Ni=1000, R0=0, options={'_Npz_limit': 10}) # doctest: +ELLIPSIS
>>> _ = pz.solve()

```
INFO    > Warning : The first moment is not an integer (Ri=3.141...
INFO    >   Try new radius Ri=3.20...
INFO    > Warning : The first moment is not an integer (Ri=6.283...
INFO    >   Try new radius Ri=6.40...
```
The radius is changed and accurate solution is reached.

To speed up the computation the evaluation of `f` and `df` can be vectorized.
For instance, here `f` can be evaluate directly for the whole `z` vector
since it involved numpy function. This behavior is driven by the boolean
options `_vectorized` (défault is `False`\). Like in
>>> pz = PZ((f, df), Rmax=5, Npz=10, Ni=1024, options={'_vectorized': True})
>>> _ = pz.solve()
>>> p, z = pz.dispatch()
>>> p_ref = np.array([-3, -1, 1, 3])*(np.pi/2)
>>> z_ref = np.array([-1, 0, 1])*(np.pi)
>>> p.sort(); np.linalg.norm(p - p_ref) < 1e-12
True
>>> z.sort(); np.linalg.norm(z- z_ref) < 1e-12
True

When it is not possible to vectorized the function, parallel evaluations can be
achieve with the `_vectorized` and `_max_workers` options. The parallelism is
done with `concurrent.futures` from the standard library. When used, it may be
better to desactivate openMP threading.

## Troubleshootings

The standard problem you may encounter using `polze` are :

  - **The first moment is not an integer**. It generally means that some roots
    may lie on the integration path; that the number of integration point is
    too small (mostly for hudge circle). `polze` will change the radius and
    try it again.

  - **Found inf in eig at Ri=...**. As the number of roots is assumed to be
    constant in all annular regions, inf may occured sometimes in the
    generalized eigenvalue problem. Those values are filtered out and the
    computation continued.

  - **Some multiplicities are far from integer**. It generally means that
    some roots are very close (the sum of their multiplicity is an int).
    Need to limit the number of root in this annular region.

`polze` outputs are manage with the `logging` module. The default verbosity
show all the trouble encountered by the algorithm. The verbosity can reduce
incressing the logger level
>>> import logging
>>> import polze
>>> logging.basicConfig(level=logging.ERROR)
>>> logging.getLogger('polze').level = logging.DEBUG

You have now complete this tuto :
```
                    ,|
                    ~i_
                     sc.
                    `Xw~
                   .7%V_
                  .rnT>.
                 ,JT>".
                ;{v".
              :+c(|.  '`
             ;*</+|.. =Ww(r="~,
           '==-`:|' .`v@aS1v+"|_:
         :\+;. .`. . Xs}JT/"-_/}v-
       :!J<"-~,  .. oRyGSsw]{))}=:
    .~!Foc<*+\~.`. id$NK4x7T{<>>:
,~"={FoF[c{v>'.'`.sHHg#iVIv>;-"!*.
*?vJ[Fo]I}!*_``-\tHhD8OjoFTI)+|_'
}IcrccJ(*(<"-'~)j@0kyamJ(?+;'-/+.
cJ]I?})??)*~,.'(TFS6yWW1r}?=;|)>.
cc]c(!>!)*\,   `|*r1stSWV{)/;_".
IIvv))?cr})^:`::,_>((!+/==/\^/\`
?!?*?(!))<*|_',`.`.
|";-:.
```
"""
import scipy.linalg as spl
from scipy.spatial import KDTree
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor

import logging
import sys

# Create console logger handler
_logging_level_default = logging.INFO
logger = logging.getLogger('polze')
# logger.addHandler(logging.NullHandler())
logger.setLevel(_logging_level_default)
# The logger is created once, but multiple handlers are created.
if not logger.hasHandlers():
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(_logging_level_default)
    formatter = logging.Formatter('%(levelname)-8s> %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


class PZ(object):
    """Poles-zeros solver class.

    Attributes
    -----------
    fdf : function or tuple of function
        the function and optionally its derivative from which the roots are looking
    Rmax : Float
        the max radius where the roots are tracked
    Npz : int
        An upper bound for the number of roots. If `_zeros_only` option is
        `True`. This parameter is ignored and computed with 0-order moment but
        still used to determined the number of required contours.
        The `_zeros_only` options is **only allowed for function without poles**.
    Ni : int
        Number of integration points. default 1024
    R0 : complex
        Origin of the countour in the complex plane. default 0
    split : int, optional
        Force to split the integration path in `split` paths.
    options : dict
        Overwrite defaut parameters. The keys should be taken in the defaut
        parameters (see below).
    """

    _Npz_limit = 5;       """Number of Pole and Zeros max befor splitting path"""
    _tol = 1e-5;          """General tolerance"""
    _clean_tol = 1e-5;    """Spurious roots multiplicity tolerance"""
    _NR_tol = 1e-12;      """Newton Raphson tolerance"""
    _NiterMax  = 15;       """Max number of Newton Raphson iterations"""
    _vectorized = False;  """Vectorized function evaluation"""
    _RiShift = 1.02;      """Scaling factor of the radius when s0 is not an integer"""
    _zeros_only = False;  """Use 0-moment for Npz estimate"""
    _Niter_for_int = 5;   """Max Number of iterations to estimate the moments"""
    _parallel = False;    """Split evaluation loop in a process pool executor"""
    _max_workers = 1;     """Number of workers in the process pool executor"""

    def __init__(self, fdf, Rmax, Npz, Ni=1024, R0=0., split=0, options=None):

        # First update defaut options
        if options:
            vars(self).update(options)

        # test if the derivitive is given
        try:
            self.f, self.df = fdf
        except:
            self.f = fdf
            self.df = None

        # dict for caching some data
        self._h_cached = dict()
        self.Rmax = Rmax

        self.R0 = R0
        self.Ni = Ni

        self._rho = Rmax
        # Initialize solution attribut
        self._reset_sol()
        # TODO change Ni if several Ri
        if (Npz > self._Npz_limit) or (split > 0):
            Nsplit = max(round(Npz / self._Npz_limit), 1)
            Nsplit = max(split, Nsplit)
            self._Ri = np.linspace(Rmax/Nsplit, Rmax, Nsplit).tolist()
            self.Npz = round(Npz/Nsplit)
        else:
            self._Ri = [Rmax]
            self.Npz = Npz

    def _reset_sol(self):
        """ Clean previous solution storage.

        The cached values for function evaluation are concerved. Maybe usefull
        if only some radius are changed.
        """
        self.pz = []
        self.m = []
        self._s = []

    def _eval(self, f, z):
        """ Eval the function `f` depending on the context (vectorized,
        parallel or not).

        Parameters
        ----------
        f : function
            The function to evaluate.
        z : array-like
            The vector where f should be evaluated.

        Returns
        ---------
        fz : array-like
            f(z) evaluations.

        Remarks
        -------
        Possible to use `_parallel` and `_vectorized` for function involving
        numpy array.
        Process parallelism is suitable for **only for heavy computation**
        since the overhead may be significant.
        """
        max_workers = self._max_workers

        # Compute f with loop or vectorized computation
        if self._vectorized:
            if self._parallel:
                # Split the array for each worker
                z_split = np.array_split(z, max_workers)
                logger.debug('Parallel run with {} workers.'.format(max_workers))
                # array for storing each process output
                fz_split = np.empty((max_workers,), dtype=object)
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    for i, fzi in enumerate(executor.map(f, z_split, chunksize=1)):
                        fz_split[i] = fzi
                # merge all sub process array
                fz = np.concatenate(fz_split)
            else:
                fz = f(z)
        else:
            # Non vectorized case
            fz = np.zeros(z.shape, dtype=complex)
            if self._parallel:
                chunksize = z.size // max_workers
                logger.debug('Parallel run with {} workers.'.format(max_workers))
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    for i, fzi in enumerate(executor.map(f, z, chunksize=chunksize)):
                        fz[i] = fzi
            else:
                for i, zi in enumerate(z):
                    fz[i] = f(zi)
        return fz

    def _eval_or_cache(self, Ri, z):
        """ Get the ratio of dfz/fz for the circle of radius Ri.

        The value are computed or reloaded from cache especially if
        `estimate_Npz` has been used.

        Parameters
        ----------
        Ri : float
            The contour radius
        z : np.array
            The coutour value.

        Returns
        -------
        h : np.array
            The dfz/fz ratio at z.

        Attributes
        -----------
        _h_cached : dict
            Update this dict with new cached values. The Keys are `Ri`
            floatting value.
        """

        # Look for cached value
        if Ri in self._h_cached.keys():
            h = self._h_cached[Ri]
            # Check only if both have the same number of integration points
            if h.size == z.size:
                logger.debug('Used cached value for Ri ={}'.format(Ri))
                return h
        else:
            # Compute f with loop or vectorized computation
            fz = self._eval(self.f, z)
            # Select analytic or finite difference derivative
            if self.df is None:
                dfz = self.diff_circle(fz, z)
            else:
                dfz = self._eval(self.df, z)
            h = dfz / fz
            # Cache new contour
            self._h_cached[Ri] = h.copy()
            return h

    def moment(self, K, Ri):
        """Compute the K-first scaled moments at fixed Ri.

        Parameters
        ----------
        K : int
            The maximal required moment order.
        Ri : float
            The radius of the contour.

        Returns
        -------
        s : np.array
            The computed moments on the circle.

        Attribut
        --------
        _h_cached : dict
            Reuse cached value for h=dfz/fz evaluation. To acces to h at Ri, the
            data are stored as a dict _h_cached[Ri]

        """
        # FIXME need to be modify, need to substract R0 for statibility ?
        R0 = self.R0
        pi = np.pi
        # integration path
        theta = 2*pi*np.arange(0, self.Ni+1)/self.Ni
        z = Ri*np.exp(1j*theta) + R0

        # Computed h=dfz/fz
        h = self._eval_or_cache(Ri, z)

        s = np.zeros((K,), dtype=complex)
        alpha = 1/(2j*pi)
        z_rho = (z)/self._rho
        for k in range(0, K):
            # TODO try to limit numerical roundoff
            s[k] = alpha * np.trapz(h * z_rho**k, z)

        # Need modify the contour in this case.
        if abs(s[0]-np.round(s[0])) > self._tol:
            logger.info('The first moment is not an integer (Ri={}, s0={}, tol={}).'.format(Ri, s[0], self._tol)
                  + ' Some roots may lie on the integration path or increase Ni ({}).'.format(self.Ni))
            return None
        else:
            return s.copy()

    def diff_circle(self, fz, z):
        """Estimate the derivatives by 9th order finite difference on a circle

        Remarks
        -------
        On a circle, the derivative with respect to z is deduce from to the
        theta derivative df/dz = df/dtheta * 1/(iz).
        """
        # TODO test vs fft derivative
        n = len(fz)
        z = z - self.R0

        h = np.abs(np.angle(z[0]/z[1]))
        #  schema centré d'ordre 9
        Coef = np.array([0., 224/280., -56/280., 32/840., -1/280.],
                        dtype=complex)
        nc = len(Coef) - 1
        # f[-1]= f[0], we must remove this redondancy
        f_ = np.concatenate((fz[-(nc+1):-1], fz, fz[1:nc+1]))

        dfdz = np.zeros(fz.shape, dtype=complex)
        for i, c in enumerate(Coef):
            # first term is 0
            if i > 0:
                dfdz += c*(f_[(nc+i):(n+nc+i)] - f_[(nc-i):(n+nc-i)])

        dfdz = dfdz / (h * 1j * z)
        return dfdz

    def estimate_Npz(self):
        """Estimate the number of zeros if there is no poles with the 0-moment.

        Returns
        -------
        Nz_tot : int
            The total number of zeros in the countour.
        Npzi_max : int
            The max number of zeros in all the annular region.
        Nzeros : array
            The number of zeros in each annular region.

        Remarks
        -------
        The multiplicity is included.
        """
        # Max number of iteration with RiShift to ensure that s0 is an int
        Niter_for_int = self._Niter_for_int
        Nzeros = np.zeros_like(self._Ri, dtype=int)
        S = 0
        for i, Ri in enumerate(self._Ri):
            s0_is_int = False
            Niter = 0
            # Compute moment and shift radius if needed
            while (not s0_is_int) and (Niter <= Niter_for_int):
                # return the moment or None
                s = self.moment(1, Ri)
                if s is not None:
                    s0_is_int = True
                else:
                    # shift the radius to avoid roots.
                    Ri = self._RiShift*Ri
                    self._Ri[i] = Ri
                    Niter += 1
                    logger.info('  Try new radius Ri={} (Niter={})'.format(Ri, Niter))
                    if Niter > Niter_for_int:
                        raise RuntimeError(('Max number of iteration is '
                                            'reached. The moment s0 is still '
                                            'not an integer a the given `_tol`.'))
            # Continue with the obtained moment s
            si = int(s.real.round()[0])
            Nzeros[i] = si - S
            # Update
            S = si

        # Output
        Nz_tot = Nzeros.sum()
        Npzi_max = Nzeros.max()
        return Nz_tot, Npzi_max, Nzeros

    def solve(self):
        """Invokes the Poles and zeros solver.

        Remarks
        -------
        Since the contour integral is splitted into several annular regions,
        the same moments orders have to be computed in all regions. `Npz`
        should be the max number of Npzi in all of these regions.
        """
        # Max number of iteration with RiShift to ensure that s0 is an int
        Niter_for_int = 5
        # Reset solution from previous run
        self._reset_sol()
        # For function without poles, Npz can be obtained numerically
        if self._zeros_only:
            _, Npz, _ = self.estimate_Npz()
            self.Npz = Npz
        else:
            Npz = self.Npz
        # Loop over all annular regions
        for i, Ri in enumerate(self._Ri):
            s0_is_int = False
            Niter = 0
            # Compute moment and shift radius if needed
            while (not s0_is_int) and (Niter <= Niter_for_int):
                # return the moment or None
                s = self.moment(2*Npz+1, Ri)
                if s is not None:
                    # store moments of Ri (for annular path) and continue
                    self._s.append(s.copy())
                    s0_is_int = True
                else:
                    # shift the radius to avoid roots.
                    Ri = self._RiShift*Ri
                    self._Ri[i] = Ri
                    Niter += 1
                    logger.info('  Try new radius Ri={} (Niter={})'.format(Ri, Niter))
                    if Niter > Niter_for_int:
                        raise RuntimeError(('  Max number of iteration is '
                                            'reached. The moment s0 is still '
                                            'not an integer a the given `_tol`.'))
            # Continue with the obtained moment s
            if i > 0:
                # It requires that all moment vector have the same size!
                s -= self._s[i-1]
            if np.round(s[0]) > Npz:
                # necessary but not sufficient condition, ex Np ~ Nz -> s0~0
                raise ValueError('The Npz upperbound ({}) if not sufficient (s0={})'.format(Npz, s[0]))
            # build Hankel matrix
            H = spl.hankel(s[0:Npz], s[Npz-1:2*Npz-1])
            H2 = spl.hankel(s[1:Npz+1], s[Npz:2*Npz])
            # solve eigenvalue problem (eig used lapack zggev that use QZ)
            D, V = spl.eig(H2, H)
            # filter if inf after since Npz is an overbound for the rank
            # TODO need to check it theoritically.
            isinf = np.isinf(D)
            if isinf.any():
                D = D[~isinf]
                NpzFinite = len(D)
                logger.warning('Found inf in eig at ' +
                               'Ri={}. Keep {} finite eigenvalue on {}.'.format(Ri,
                                                                                NpzFinite,
                                                                                Npz))
            else:
                NpzFinite = Npz

            # Compute the multiplicity by solving the Vandermonde systeme.
            Vander = np.vander(D, increasing=True)
            m = spl.solve(Vander.T, s[0:NpzFinite])
            # rescale the roots
            pz = D*self._rho  # + self.R0
            # store it
            self.m.append(m.copy())
            self.pz.append(pz.copy())

        # and return
        return self.pz, self.m

    def iterative_ref(self):
        r"""Use Newton-Raphson method to refine the roots.

        Suitable when f and df are not poluted by roundoff close to the roots.
        The method requires that df is given. If not, an approximation using
        finite difference is computed.

        As the root multiplicity \(\mu\) is known, this modified NR method
        preserves the quadratic convergence. Following [1], for zeros
        $$
        z_{n+1}=z_{n}-\mu{\frac {f(z_{n})}{f'(z_{n})}}.
        $$
        and for the Poles
        $$
        z_{n+1}=z_{n}+\mu{\frac {f(z_{n})}{f'(z_{n})}}.
        $$
        Since the stored multiplicities are \(-\mu_n\) for pole, the formula
        are finally identical.

        Returns
        -------
        pz_raf : array
            The refined roots.
        mus : array
            The refine multiplicities.
        info : dict
            Information about the refinement step.

        Remarks
        -------
        The orginal poles and zeros are not updated.
        """
        f = self.f
        if self.df is None:
            logger.info('NR requires the derivative. Use FD approximation...')
            def df(z):
                """ Five order finite difference approximation of df.
                """
                # eps = np.finfo(float).eps
                h = 1e-4  # good compromize with cancellation error for f~1 with O(h**4)
                # if f has zeros at z, Taylor series converge
                if mu_ > 0:
                    fm2, f2 = f(z-2*h), f(z+2*h)
                    f1, fm1 = f(z+h), f(z-h)
                    dfz = ((fm2 - f2)/12. + (f1 - fm1)*(2./3.))/h

                else:
                    # if f has pole at z, Taylor series DO NOT converge
                    # but 1/f Taylor series does...
                    g = lambda x: 1/f(x)
                    gm2, g2 = g(z-2*h), g(z+2*h)
                    g1, gm1 = g(z+h), g(z-h)
                    dgz = ((gm2 - g2)/12. + (g1 - gm1)*(2./3.))/h
                    g0 = 1/f(z)
                    dfz = - dgz/g0**2
                return dfz
        else:
            df = self.df

        pzs, mus = self.clean()
        pz_raf = np.zeros(pzs.shape, dtype=complex)
        Err = -np.ones(pzs.shape, dtype=complex)
        Niter = np.zeros(pzs.shape, dtype=np.int32)
        # loop over roots
        for n, (pz, mu) in enumerate(zip(pzs, mus)):
            # init NR method
            pz_, mu_ = pz, np.round(mu)
            niter = 0
            err = 1/self._NR_tol
            # apply NR method
            while (err > self._NR_tol) and (niter < self._NiterMax):
                fz = f(pz_)
                dfz = df(pz_)
                pz_old = pz_
                pz_ = pz_ - mu_*fz/dfz
                niter += 1
                # compute relative error for non 0 roots
                if abs(pz) > self._tol:
                    err = abs(pz_ - pz_old)/abs(pz_)
                else:
                    err = abs(pz_ - pz_old)
            # store solution and Error for check
            pz_raf[n] = pz_
            Err[n] = err
            Niter[n] = niter

        ErrTot = np.sum(np.abs(Err))
        if (np.abs(Err) > self._NR_tol).any():
            logger.info('Refine step finalized. Some roots have not converged. Please check `info`.')

        # package output
        info = {'Err': Err, 'Niter': Niter, 'ErrTot': ErrTot}
        return pz_raf, mus, info

    def contour_ref(self, Ni=500, scaling=0.01):
        r"""Use contour solver to refine the roots.

        This approach creates a new PZ object around already found roots
        to get a better accuracy. The contour radius is obtained by
        `scaling` the distance to the closest roots.

        This approach is sometimes more robust than NR refinement for functions
        subject to rounding error close to the roots (ex: det M(z) = 0) or
        close to poles, but is significantly more longer.

        Parameters
        ----------
        Ni : int
            The number of integration points.
        scaling : float
            The scaling used to compute the contour radius. Default 0.01.

        Returns
        -------
        pz_raf : array
            The refined roots.
        mus : array
            The refine multiplicities.
        info : dict
            Information about the refinement step.

        Remarks
        -------
        The orginal poles and zeros are not updated.
        """
        pzs, mus = self.clean()
        pz_raf = np.zeros(pzs.shape, dtype=complex)
        mus_raf = np.zeros(mus.shape, dtype=complex)
        Err = np.zeros(pzs.shape)
        radii = np.zeros(pzs.shape)
        # Estimate minimal distance between roots using Kdtree computation
        # perhaps overkilling...
        points = np.array([pzs.real, pzs.imag]).T
        tree = KDTree(points)
        for i, p in enumerate(points):
            d, _ = tree.query(p, k=2)
            # d[0] is always 0 (distance to itself)
            radii[i] = d[1] * scaling

        # define the function to solve
        if self.df is None:
            fdf = self.f
        else:
            fdf = (self.f, self.df)

        # loop over roots
        for n, (pz, mu) in enumerate(zip(pzs, mus)):
            pz_ref = PZ(fdf, Rmax=radii[n], Npz=1, Ni=Ni, R0=pz,
                        options={'_vectorized': self._vectorized,
                                 '_parallel': self._parallel,
                                 '_max_workers': self._max_workers})
            pz_ref.solve()
            pz_, mus_ = pz_ref.clean()
            # Expect only one value
            if pz_.size > 1:
                logger.warning('Several roots have been found during contour refinement. Expect one.')
            # store solution and Error for check
            pz_raf[n] = pz_
            mus_raf[n] = mus_
            Err[n] = abs(pz_ - pz)

        # Package output
        info = {'deviation': Err, 'radii': radii}
        return pz_raf, mus_raf, info

    def clean(self, tol=None, split=False):
        """Clean spurious roots by removing 0 multiplicity roots.

        Parameters
        ----------
        tol : float, optional
            The tolerance use to filter spurious roots. The default is store in
            `_clean_tol` attribut.
        split : bool, optional
            Return clean roots for each annular contour.

        Returns
        -------
        pz : np.array
            All found poles and zeros.
        m : np.array
            All the associated multiplicities.
        """
        if tol is None:
            tol = self._clean_tol
        pz = []
        mu = []
        for i, m in enumerate(self.m):
            ind = np.abs(m) > tol
            # Warn if multiplicity is too far from an integer
            dist_mult_to_int = np.abs(m[ind] - np.round(m[ind].real))
            if (dist_mult_to_int > tol).any():
                logger.warning('Some multiplicities are far from integer ' +
                               'at Ri={}. (Max dist = {}, tol = {}).'.format(self._Ri[i],
                                                                             dist_mult_to_int.max(),
                                                                             tol))
            # TODO add a sort option
            pz.append(self.pz[i][ind])
            mu.append(m[ind])

        if split:
            # return the filtered lists
            return pz, mu
        else:
            # return an single array
            return np.concatenate(pz), np.concatenate(mu)

    def display(self, clean=True):
        """Pretty print of poles and zeros."""

        if clean:
            pz, m = self.clean(split=False)
        else:
            pz, m = np.concatenate(self.pz), np.concatenate(self.m)
        print('\n')
        print('-----------------------------------------------------------------------------------------------------------------------------')
        print('|                  poles                 |                   zeros                    |                    mu                ')
        print('-----------------------------------------------------------------------------------------------------------------------------')
        for i, (p, mu) in enumerate(zip(pz, m)):
            if np.round(mu) > 0:
                print('                                             {}    {}'.format(p, mu))
            else:
                print('{}  |                                            {}'.format(p, mu))
        print('-----------------------------------------------------------------------------------------------------------------------------')

    def plot(self, ax=None, variable='\\nu'):
        """ Plot the zeros (o) and the poles (+).

        Parameters
        -----------
        ax : axis
            An axis can be provided to overplot roots.

        Returns
        ---------
        ax : axis
            The axis object.
        """
        split = True
        # list of colors
        col = [None, 'k', 'b']
        # get geniuine  as list of array
        pz, m = self.clean(split=split)
        col_index = -1
        # plt.xkcd()
        if ax is None:
            fig = plt.figure()
            ax = fig.gca()

        for pzi, mi in zip(pz, m):
            ind = mi.real > 0
            ax.plot(pzi[ind].real, pzi[ind].imag, '.', color=col[col_index])  # zeros
            ax.plot(pzi[~ind].real, pzi[~ind].imag, '+', color=col[col_index])  # poles
            # color alternance
            col_index *= -1

        ax.set_xlabel(r'$\mathrm{Re}\,' + variable + r' $')
        ax.set_ylabel(r'$\mathrm{Im}\,' + variable + r' $')
        ax.axis('scaled')
        for Ri in self._Ri:
            ax.add_artist(plt.Circle((self.R0.real, self.R0.imag), radius=Ri,
                                     edgecolor='gray', facecolor=None,
                                     fill=False))
        return ax

    def mapz(self, NgridPoint=51, variable='z', bounds=None, ax=None,
             addsols=False, part='modulus'):
        """ Create a map to explore the function in the complex plane.

        Parameters
        ----------
        NgridPoint : int
            number of grid point in each direction.
        variable: str, default: '\\nu'
            Variable name for axis labels
        bounds : tuple
            contains the z-plane corners ex : (-1-1j,1+1j). If None, `Rmax`
            attribute is use.
        ax : axis
            The plotting axis.
        addsols : bool
            If true overlay the found poles and zeros.
        part : string
            The quantity to plot when dealing with complex function. Default
            is 'modulus'. Allowed values are {'real', 'imag', 'modulus'}.

        Returns
        ---------
        ax : axis
            the axis object
        """

        if bounds is None:
            bounds = (-self.Rmax*(1+1j) + self.R0, self.Rmax * (1+1j) + self.R0)

        z_real, z_imag = np.meshgrid(np.linspace(bounds[0].real, bounds[1].real, NgridPoint),
                                     np.linspace(bounds[0].imag, bounds[1].imag, NgridPoint))
        Z = 1.*(z_real+1j*z_imag)
        # Use eval common interface
        FZ = self._eval(self.f, Z.ravel()).reshape(Z.shape)

        # plot
        if not(ax):
            fig = plt.figure()
            ax = plt.gca()
        if part == 'modulus':
            contour_F = ax.pcolormesh(Z.real, Z.imag, np.log10(np.abs(FZ)),
                                      zorder=-1, shading='auto')
        elif part == 'real':
            contour_F = ax.pcolormesh(Z.real, Z.imag, FZ.real,
                                      zorder=-1, shading='auto')
        elif part == 'imag':
            contour_F = ax.pcolormesh(Z.real, Z.imag, FZ.imag,
                                      zorder=-1, shading='auto')
        else:
            raise ValueError("Input argument `part` should belong "
                             "to (`modulus`, `real`, `imag`).")
        ax.set_xlabel(r'$\mathrm{Re}\,' + variable + r' $')
        ax.set_ylabel(r'$\mathrm{Im}\,' + variable + r' $')
        plt.colorbar(contour_F, ax=ax)

        if addsols:
            self.plot(ax=ax, variable=variable)

        return ax, FZ

    @staticmethod
    def _dispatch(pz, m, multiplicities=False):
        """Separate Zeros and poles from `pz` using multiplicity `m`.

        Parameters
        -----------
        pz : np.array
            The array with the poles and the zeros
        m : np.array
            The array containing the roots multiplicities.
        multiplicities : bool
            If multiplicities is `True`, return 2 tuples with (pole, mult.) and
            (zeros, mult.). Else only poles and zeros are returned.
        """
        ind = m.real > 0
        p = pz[~ind]
        z = pz[ind]
        if multiplicities:
            return (p, m[~ind]), (z, m[ind])
        else:
            return p, z

    def dispatch(self, refine=False, multiplicities=False):
        """Separate Zeros and pole using multiplicity.

        This method is the short hand for the method `_distpach` (see its doc).

        Parameters
        ----------
        refine : bool
            If `True` performed an iterative refinement before returning the
            roots. See also `PZ.iterative_ref` method.
        multiplicities : bool
            See `_distpach`.

        If `multiplicties` is `False` (default)

        Returns
        --------
        p : array
            poles locations.
        z : array
            zeros locations.

        Else returns 2 tuples for p and z respectivelly and their associated
        multiplities.
        """
        if refine is False:
            pz, m = self.clean(split=False)
        else:
            pz, m, info = self.iterative_ref()

        return self._dispatch(pz, m, multiplicities=multiplicities)

    def check(self, clean=True, tol=1e-5, shift=1.0095+0.00295j):
        """ Check if the roots are accurate using function evaluation.

        Parameters
        ----------
        clean : bool
            If the root must be _cleaned_ using the multiplicities before the
            check. If `True` it will excluded _spurious_ root from the check.
        tol : float
            The tolerance use to validate the check.
        shift : complex
            Small shift to compare the function value at `root*shift`.

        Returns
        -------
        status : Status
            Object containing `errz`,  `errp` and `all_roots`. This last
            boolean parameter says that all the roots satisfy the `tol`. The
            other contains the function values at the roots and the ratio with
            the shifted position.

        """
        # Define a namedtupled for summurize the output results
        Status = namedtuple('Status', 'errz errp all_roots')
        # check after filter out spurious from multiplicity
        if clean:
            pz, m = self.clean(split=False)
        else:
            pz, m = self.pz, self.m
        # if list -> np.array
        if isinstance(pz, list):
            pz = np.concatenate(pz)
        if isinstance(m, list):
            m = np.concatenate(m)
        # treat separatly poles and zeros
        ps, zs = self._dispatch(pz, m)
        errz = np.zeros((3, zs.size), dtype=complex)
        errp = np.zeros((3, ps.size), dtype=complex)
        # zeros
        # TODO see if it robust if 1/0 or other exceptions...
        for i, z in enumerate(zs):
            errz0 = abs(self.f(z))
            errz_shift = abs(self.f(z*shift))
            errz[0, i] = errz0
            errz[1, i] = abs(errz_shift) / abs(errz0)
            errz[2, i] = z
        # poles
        for i, p in enumerate(ps):
            errp0 = abs(self.f(p))
            errp_shift = abs(self.f(z*shift))
            errp[0, i] = errp0
            errp[1, i] = abs(errp_shift) / abs(errp0)
            errp[2, i] = p

        # True if all the roots are real roots at the fixed tol.
        all_roots = (errz[0, :] < tol).all() and (1/errp[0, :] < tol).all()
        # Encapsulate into Status named-tuple
        status = Status(all_roots=all_roots, errp=errp, errz=errz)
        return status


if __name__ == '__main__':
    """ Run doctests.
    """
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=True)
