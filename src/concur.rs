//! General Constrained Curve (K-Degree) Spline-Fit for Multi-Dimensional (N) Data
//! 
//! Rust wrapper of Dierckx' `concur` Fortran subroutine, with 28(!) input parameters:
//! 
//! ```fortran
//!  subroutine concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,nest,n,t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
//! ```
//! Wrapper `ConstrainedSpline<K,N>` encapsulates most of this data, and implements methods to build spline
//! representations of a set of parameterized input coordinates, and a set of end point constraints for both 
//! end points. These constraints are derivative constraints, such as for a cubic spline first and second derivatives,
//! for each of the dimension, e.g. (dx/du, dx<sup>2</sup>/du<sup>2</sup>, dy/du, dy<sup>2</sup>/du<sup>2</sup>).



use std::iter::repeat;
//use crate::dierckx::{concur_};
use dierckx_sys::{concur_};
use super::{FitError};
use crate::Result;
use spliny::SplineCurve;


#[derive(Clone)]
pub struct ParameterSplineCurveFit<const K:usize, const N:usize> {
    // input values
    xn: Vec<f64>, // data (x,y,..) coordinates
    u: Vec<f64>,
    w: Vec<f64>,    // weight factors, 
    xb: Vec<f64>,
    xe: Vec<f64>,

    t: Vec<f64>,
    c: Vec<f64>,
    e_rms: Option<f64>,
    n: i32,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries
    xx: Vec<f64>,
    cp: Vec<f64>,
    ib: i32,
    ie: i32,
    m: i32,
    mx: i32,
    nest: i32,
    k: i32,
    idim: i32,
}


/**
  
 Fit parametric B-Spline curve to a set of coordinates

 Wrapper for Dierckx' `concur` subroutine.

 */

impl<const K:usize, const N:usize> ParameterSplineCurveFit<K, N> {

    /// Constructor, taking as input cuve parameter, curve coordinates, and end point constraints, and setting up
    /// remaining datastructures and values for `concur`.
    /// The curve parameter is `u`.
    /// 
    /// Coordinates are represented by the vector `xn`, starting with the coordinates of the first point;
    /// for example, if N=3, a three dimensional space, with coordinates given as (x,y,z), the array can be
    /// constructed as [x0, y0, z0, x1, y1, z1, x2 ...]. Its the number of coordinates is m, its size is 
    //  N * m.
    
    pub fn new(
        u: Vec<f64>,
        xn: Vec<f64>,
    ) -> Result<Self> {

        let k = K as i32;
        if ![1,3,5].contains(&(k as i32)) { return Err(FitError(208).into()) };
        let idim =  if (1..=10).contains(&N) { N as i32 } else {
                return Err(FitError(200).into())
            };
        let m = u.len() as i32; // number of coordinates
        if m<2 {return Err(FitError(201).into())};
        let mx = m * idim;
        if xn.len() as i32!= mx { return Err(FitError::new(202).into())}
        let w_vec = vec![1.0;m as usize];

        let xb = Vec::new();
        let ib = 0;

        //let xe = xn[xn.len()-idim-1..xn.len()].to_vec();
        let xe = Vec::new();
        let ie = 0;

        let nest = m+k+1 + 2*(k-1); 
        let n = nest;  // length of tc
        let t_vec = vec![0.0; nest as usize];
        let c_vec = vec![0.0; (nest * idim) as usize];

        let iwrk_vec = vec![0i32; nest as usize];

        let wrk_vec = vec![0f64; (m*(k+1)+nest*(6+idim+3*k)) as usize];
        let xx_vec = vec![0.0; (idim*m) as usize];
        let cp_vec = vec![0.0; (2 * (k+1) * idim) as usize];

        Ok(Self { u, xn, w: w_vec, xb, xe, t: t_vec, c: c_vec, wrk: wrk_vec, iwrk: iwrk_vec, xx: xx_vec, cp: cp_vec, ib, ie, m, mx, nest, k, idim, n, e_rms: None})

    }

    pub fn u(&self) -> Vec<f64> {
        let mut v: Vec<f64> = Vec::with_capacity(self.n as usize);
        v.extend_from_slice(&self.u[0..self.n as usize]);
        v
    }

    pub fn xn(&self) -> Vec<f64> {
        let mut v: Vec<f64> = Vec::with_capacity((self.n*self.idim) as usize);
        v.extend_from_slice(&self.xn[0..(self.n*self.idim) as usize]);
        v
    }

    pub fn weights(mut self, weights: Vec<f64>) -> Result<Self> {
        if weights.len() == self.u.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err(FitError(203).into())
        }
    }

    pub fn begin_constraints<const D: usize>(mut self, ub: [[f64;N];D]) -> Result<Self> {
        if D<=(K+1)/2+1 {
            self.xb = ub.iter().flatten().cloned().collect();
            self.ib = D as i32 -1;
            Ok(self)
        } else {
            Err(FitError(204).into())
        }
    }

    pub fn end_constraints<const D: usize>(mut self, ub: [[f64;N];D]) -> Result<Self> {
        if D<=(K+1)/2+1 {
            self.xe = ub.iter().flatten().cloned().collect();
            self.ie = D as i32 -1;
            Ok(self)
        } else {
            Err(FitError(204).into())
        }
    }

    fn concur(&mut self, iopt:i32, e_rms:Option<f64>, knots: Option<Vec<f64>>) ->  i32 {
        let mut fp = 0.0;
        let s = if let Some(e) = e_rms {
            self.m as f64 * e.powi(2)
        } else {
            0.0
        };

        let nb = self.xb.len() as i32;
        let ne = self.xe.len() as i32;
        let np = self.cp.len() as i32;
        let nc = self.c.len() as i32;
        let lwrk = self.wrk.len() as i32;

        if let Some(knots) = knots {
            self.n = knots.len() as i32;
            self.t = knots;
        }
        let mut ierr = 0;
        unsafe {
            concur_(
                &iopt,
                &self.idim,
                &self.m,
                self.u.as_ptr(),
                &self.mx,
                self.xn.as_ptr(),
                self.xx.as_mut_ptr(), 
                self.w.as_ptr(),
                &self.ib,
                self.xb.as_ptr(),
                &nb,
                &self.ie,
                self.xe.as_ptr(),
                &ne,
                &self.k,
                &s,
                &self.nest,
                &mut self.n,
                self.t.as_mut_ptr(),
                &nc,
                self.c.as_mut_ptr(),
                &np,
                self.cp.as_mut_ptr(),
                &mut fp,
                self.wrk.as_mut_ptr(),
                &lwrk,
                self.iwrk.as_mut_ptr(),
                &mut ierr,
            );
        }
       // self.tc.t.truncate(n); //todo to in from
       // self.tc.c.truncate(n);
       self.e_rms = Some((fp/self.m as f64).sqrt());
       ierr
    }


    /**
     * Cardinal Spline: Weighted least squares spline with equidistant knots
     * 
     * Returns Spline, and rms error, with knots dt (input parameter) apart,
     * and aligned to integer multiples of it. Knots cover the range within
     * the bounds of x.
     */
    pub fn cardinal_spline(mut self, dt:f64) -> Result<SplineCurve<K,N>>{
        let m = self.u.len();
        let tb = ((self.u[0]+f64::EPSILON)/dt).ceil() * dt; // inner knots should bot be equal
        let te = ((self.u[m-1]-f64::EPSILON)/dt).floor() * dt;
        let n = ((te - tb)/dt).round() as usize;
        if n == 0 { return Err(FitError(205).into()) };
        let mut t: Vec<f64> = Vec::with_capacity(n + 2 * (K + 1) + 1);
        t.extend( repeat(tb).take(K+1) // begin padding, needed for spline evaluation
            .chain(
                repeat(dt).scan(tb, 
                    |s, dx|{
                        let t=*s; 
                        *s+=dx; 
                        if t<=te {
                            Some(t)
                        } else {
                            None
                        }
                    })
            )
            .chain(
                repeat(te).take(K+1) // end padding
            )
        );

        let ierr = self.concur(-1, Some(0.0),Some(t));
        if ierr <= 0 {
            Ok(self.into())
        } else {
            Err(FitError(ierr).into())
        }
    }

    /**
     Interpolating Spline
      
      
     */ 
    pub fn interpolating_spline(mut self) -> Result<SplineCurve<K,N>> {
        let ierr = self.concur(0, Some(0.0),None);
        if ierr<=0  {
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
    pub fn smoothing_spline(mut self, rms: f64) -> Result<SplineCurve<K,N>>{
        let ierr= self.concur(0, Some(rms), None);
        if ierr>0 {
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
    pub fn smoothing_spline_optimize(mut self, 
            rms_start: f64, 
            rms_scale_ratio: f64, 
            converged: impl Fn(i32, i32, f64, f64) -> bool, 
            n_iter: Option<usize>,
        ) -> Result<SplineCurve<K,N>>{
        let n_iter = n_iter.unwrap_or(40);
        let ierr= self.concur(0, Some(rms_start), None);
        if ierr>0 {
            return Err(FitError(ierr).into())
        }
        let mut rms = self.e_rms.unwrap();
        let mut n_prev;
        let mut rms_prev ;
        for _ in 0..n_iter {
            n_prev = self.n;
            rms_prev = rms;
            let ierr= self.concur(1, Some(rms * rms_scale_ratio), None);
            rms = self.e_rms.unwrap();
            if ierr>0 {
                return Err(FitError(ierr).into())
            }
            if converged(self.n, self.n-n_prev, rms, rms_prev - rms) { // finishing fit
                let ierr = self.concur(0, Some(rms_prev), None);
                if ierr>0 {
                    return Err(FitError(ierr).into())
                } else {
                    return Ok(self.into())
                }
            };
        }
        Err(FitError(206).into())
    }


    

} // impl ParametricCurveSplineFit



impl<const K:usize, const N:usize> From<ParameterSplineCurveFit<K,N>> for SplineCurve<K,N> {
    fn from(mut sp: ParameterSplineCurveFit<K,N>) -> Self {
        sp.t.truncate(sp.n as usize);
        sp.t.shrink_to_fit();

        sp.c.truncate((sp.n*sp.idim) as usize); // this is the size as returned, but this conains K+1 unused values at the end

        for dim in 0..sp.idim as usize {
            let ib = (dim+1) * (sp.n-sp.k-1) as usize;
            let ie = ib + sp.k as usize + 1;
            sp.c.drain(ib..ie);
        }
        sp.c.shrink_to_fit();

        Self::new(
            sp.t,
            sp.c,

        )
    }
}

