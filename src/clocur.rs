//! Closed Periodic Curve (K-Degree) Spline-Fit for Multi-Dimensional (N) Data
//!
//! Rust wrapper of Dierckx' `clocur` Fortran subroutine for fitting closed/periodic splines:
//!
//! ```fortran
//!  subroutine clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)
//! ```
//! Wrapper `ClosedParameterSplineCurveFit<K,N>` encapsulates most of this data, and implements
//! methods to build closed spline representations of a set of parameterized input coordinates.
//! To build periodic splines where the curve forms a closed loop, the first and last data points must
//! coincide.

use super::FitError;
use crate::Result;
use dierckx_sys::clocur_;
use spliny::SplineCurve;

#[derive(Clone)]
pub struct ClosedParameterSplineCurveFit<const K: usize, const N: usize> {
    // input values
    xn: Vec<f64>, // data (x,y,..) coordinates
    u: Vec<f64>,
    w: Vec<f64>, // weight factors
    ipar: i32,

    t: Vec<f64>,
    c: Vec<f64>,
    e_rms: Option<f64>,
    n: i32,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries,
    m: i32,
    mx: i32,
    nest: i32,
    k: i32,
    idim: i32,
}

/**

Fit closed parametric B-Spline curve to a set of coordinates

Wrapper for Dierckx' `clocur` subroutine.

*/
impl<const K: usize, const N: usize> ClosedParameterSplineCurveFit<K, N> {
    /// Constructor, taking as input cuve parameter, curve coordinates, and end point constraints, and setting up
    /// remaining datastructures and values for `clocur`.
    /// The curve parameter is `u`.
    ///
    /// Coordinates are represented by the vector `xn`, starting with the coordinates of the first point;
    /// for example, if N=3, a three dimensional space, with coordinates given as (x,y,z), the array can be
    /// constructed as [x0, y0, z0, x1, y1, z1, x2 ...]. Its the number of coordinates is m, its size is
    /// N * m. The first and last data points must coincide.

    pub fn new(u: Vec<f64>, xn: Vec<f64>) -> Result<Self> {
        let k = K as i32;
        if ![1, 3, 5].contains(&k) {
            return Err(FitError(208).into());
        };
        let idim = if (1..=10).contains(&N) {
            N as i32
        } else {
            return Err(FitError(200).into());
        };
        let m = u.len() as i32;
        if m < 2 {
            return Err(FitError(201).into());
        };
        let mx = m * idim;
        if xn.len() as i32 != mx {
            return Err(FitError::new(202).into());
        }

        // Validate that first and last points coincide
        let n_dim = N;
        let m_usize = m as usize;
        for d in 0..n_dim {
            let first = xn[d];
            let last = xn[(m_usize - 1) * n_dim + d];
            if (first - last).abs() > 1e-10 {
                return Err(
                    format!("clocur requires first and last points to coincide, but dimension {d} differs: {first} vs {last}").into()
                );
            }
        }

        let ipar = 1; // user-supplied parameter values
        let w_vec = vec![1.0; m as usize];

        let nest = m + 2 * k;
        let n = 0;
        let t_vec = vec![0.0; nest as usize];
        let c_vec = vec![0.0; (nest * idim) as usize];

        let iwrk_vec = vec![0i32; nest as usize];

        let lwrk = m * (k + 1) + nest * (7 + idim + 5 * k);
        let wrk_vec = vec![0f64; lwrk as usize];

        Ok(Self {
            u,
            xn,
            w: w_vec,
            ipar,
            t: t_vec,
            c: c_vec,
            wrk: wrk_vec,
            iwrk: iwrk_vec,
            m,
            mx,
            nest,
            k,
            idim,
            n,
            e_rms: None,
        })
    }

    pub fn weights(mut self, weights: Vec<f64>) -> Result<Self> {
        if weights.len() == self.u.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err(FitError(203).into())
        }
    }

    fn clocur(&mut self, iopt: i32, e_rms: Option<f64>, knots: Option<Vec<f64>>) -> i32 {
        let mut fp = 0.0;
        let s = if let Some(e) = e_rms {
            self.m as f64 * e.powi(2)
        } else {
            0.0
        };

        let nc = self.c.len() as i32;
        let lwrk = self.wrk.len() as i32;

        if let Some(knots) = knots {
            self.n = knots.len() as i32;
            self.t = knots;
        }
        let mut ierr = 0;
        unsafe {
            clocur_(
                &iopt,
                &self.ipar,
                &self.idim,
                &self.m,
                self.u.as_ptr(),
                &self.mx,
                self.xn.as_ptr(),
                self.w.as_ptr(),
                &self.k,
                &s,
                &self.nest,
                &mut self.n,
                self.t.as_mut_ptr(),
                &nc,
                self.c.as_mut_ptr(),
                &mut fp,
                self.wrk.as_mut_ptr(),
                &lwrk,
                self.iwrk.as_mut_ptr(),
                &mut ierr,
            );
        }
        self.e_rms = Some((fp / self.m as f64).sqrt());
        ierr
    }

    /**
    Interpolating Spline


    */
    pub fn interpolating_spline(mut self) -> Result<SplineCurve<K, N>> {
        let ierr = self.clocur(0, Some(0.0), None);
        if ierr <= 0 {
            Ok(self.into())
        } else {
            Err(FitError(ierr).into())
        }
    }

    /**
     * Fit a smoothing spline
     *
     * nrguments:
     * - rms: root mean square error
     *
     */
    pub fn smoothing_spline(mut self, rms: f64) -> Result<SplineCurve<K, N>> {
        let ierr = self.clocur(0, Some(rms), None);
        if ierr > 0 {
            Err(FitError(ierr).into())
        } else {
            Ok(self.into())
        }
    }

    /**
     * Fit a best-fit smoothing spline by decreasing rms target
     *
     * nrguments:
     * - rms_start: root mean square error start value
     * - rms_scale_ratio: `rms *= rms_scale_ratio` for each step
     * - converged: boolean convergence function, with arguments
     *   - number of knots (usize),
     *   - number of added knots in the last iteration (usize)
     *   - root mean square error of the fit (f64),
     *   - root mean square error improvement in the last iteration
     * - rms_scale_ratio: amount to reduce the target rms value in each iteration, default = 0.8
     * - n_iter: number of iterations
     *
     */
    pub fn smoothing_spline_optimize(
        mut self,
        rms_start: f64,
        rms_scale_ratio: f64,
        converged: impl Fn(i32, i32, f64, f64) -> bool,
        n_iter: Option<usize>,
    ) -> Result<SplineCurve<K, N>> {
        let n_iter = n_iter.unwrap_or(40);
        let ierr = self.clocur(0, Some(rms_start), None);
        if ierr > 0 {
            return Err(FitError(ierr).into());
        }
        let mut rms = self.e_rms.unwrap();
        let mut n_prev;
        let mut rms_prev;
        for _ in 0..n_iter {
            n_prev = self.n;
            rms_prev = rms;
            let ierr = self.clocur(1, Some(rms * rms_scale_ratio), None);
            rms = self.e_rms.unwrap();
            if ierr > 0 {
                return Err(FitError(ierr).into());
            }
            if converged(self.n, self.n - n_prev, rms, rms_prev - rms) {
                let ierr = self.clocur(0, Some(rms_prev), None);
                if ierr > 0 {
                    return Err(FitError(ierr).into());
                } else {
                    return Ok(self.into());
                }
            };
        }
        Err(FitError(206).into())
    }
} // impl ClosedParametricCurveSplineFit

impl<const K: usize, const N: usize> From<ClosedParameterSplineCurveFit<K, N>>
    for SplineCurve<K, N>
{
    fn from(mut sp: ClosedParameterSplineCurveFit<K, N>) -> Self {
        sp.t.truncate(sp.n as usize);
        sp.t.shrink_to_fit();

        sp.c.truncate((sp.n * sp.idim) as usize);

        for dim in 0..sp.idim as usize {
            let ib = (dim + 1) * (sp.n - sp.k - 1) as usize;
            let ie = ib + sp.k as usize + 1;
            sp.c.drain(ib..ie);
        }
        sp.c.shrink_to_fit();

        Self::new(sp.t, sp.c)
    }
}
