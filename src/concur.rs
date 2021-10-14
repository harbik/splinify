//! Spline Representation of Endpoint Constrainted Parameterized Curves
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
use crate::dierckx::{concur_};
use super::{Spline, DierckxError};

use crate::Result;

pub struct ParametricCurveSplineFit<const K:usize, const N:usize> {
    // input values
    xn: Vec<f64>, // data (x,y,..) coordinates
    u: Vec<f64>,
    w: Vec<f64>,    // weight factors, 
    xb: Vec<f64>,
    xe: Vec<f64>,

    t: Vec<f64>,
    c: Vec<f64>,
    n: usize,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries
    xx: Vec<f64>,
    cp: Vec<f64>,
    ib: usize,
    ie: usize,
    m: usize,
    mx: usize,
    nest: usize,
    k: usize,
    idim: usize,
}


/**
  
 A ConstrainedSplineCurve builder, using Dierckx' `concur` subroutine.

 */
impl<const K:usize, const N: usize> ParametricCurveSplineFit<K, N> {

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

        let k = K;
        if ![1,3,5].contains(&(k as i32)) { return Err("K should be 1, 3 or 5".into()) };
        let idim =  if (1..=10).contains(&N) { N } else {
                return Err(DierckxError(200).into())
            };
        let m = u.len(); // number of coordinates
        if m<2 {return Err(DierckxError(201).into())};
        let mx = m * N;
        if xn.len()!= mx { return Err(DierckxError::new(202).into())}
        let w = vec![1.0;m];

        let xb = Vec::new();
        let ib = 0;

        //let xe = xn[xn.len()-idim-1..xn.len()].to_vec();
        let xe = Vec::new();
        let ie = 0;

        let nest = m+K+1 + 2*(K-1); 
        let n = nest;  // length of tc
        let t = vec![0.0; nest];
        let c = vec![0.0; nest * N];

        let iwrk = vec![0i32; nest];

        let wrk = vec![0f64; m*(K+1)+nest*(6+N+3*K)];
        let xx = vec![0.0; N*m];
        let cp = vec![0.0; 2 * (K+1) * N];

        Ok(Self { u, xn, w, xb, xe, t, c, wrk, iwrk, xx, cp, ib, ie, m, mx, nest, k, idim, n})

    }

    pub fn weights(mut self, weights: Vec<f64>) -> Result<Self> {
        if weights.len() == self.u.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err(DierckxError(203).into())
        }
    }

    pub fn begin_constraints<const D: usize>(mut self, ub: [[f64;N];D]) -> Result<Self> {
        if D<=(K+1)/2+1 {
            self.xb = ub.iter().flatten().cloned().collect();
            self.ib = D-1;
            Ok(self)
        } else {
            Err(DierckxError(204).into())
        }
    }

    pub fn end_constraints<const D: usize>(mut self, ub: [[f64;N];D]) -> Result<Self> {
        if D<=(K+1)/2+1 {
            self.xe = ub.iter().flatten().cloned().collect();
            self.ie = D-1;
            Ok(self)
        } else {
            Err(DierckxError(204).into())
        }
    }

    fn concur(&mut self, iopt:i32, e_rms:Option<f64>, knots: Option<Vec<f64>>) ->  (i32, f64) {
        let mut fp = 0.0;
        let s = if let Some(e) = e_rms {
            self.m as f64 * e.powi(2)
        } else {
            0.0
        };

        let nb = self.xb.len();
        let ne = self.xe.len();
        let np = self.cp.len();
        let nc = self.c.len();
        let lwrk = self.wrk.len();

        if let Some(knots) = knots {
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
        (ierr, (fp/self.m as f64).sqrt()) // fit error, in space coordinates
    }


    /**
     * Cardinal Spline: Weighted least squares spline with equidistant knots
     * 
     * Returns Spline, and rms error, with knots dt (input parameter) apart,
     * and aligned to integer multiples of it. Knots cover the range within
     * the bounds of x.
     */
    pub fn cardinal_spline(mut self, dt:f64) -> Result<(Spline<K>,f64)>{
        let m = self.u.len();
        let tb = (self.u[0]/dt).ceil() * dt;
        let te = (self.u[m-1]/dt).floor() * dt;
        let n = ((te - tb)/dt).round() as usize;
        if n == 0 { return Err(DierckxError(205).into()) };
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
         //   .collect();

        let (ierr, fp) = self.concur(-1, Some(0.0),Some(t));
        if ierr<=0  {
            Ok((self.into(), fp))
        } else {
            Err(DierckxError(ierr).into())
        }
    }

    /**
     Interpolating Spline
      
      
     */ 
    pub fn interpolating_spline(mut self) -> Result<Spline<K>> {
        let (ierr, _fp) = self.concur(0, Some(0.0),None);
        if ierr<=0  {
            Ok(self.into())
        } else {
            Err(DierckxError(ierr).into())
        }
    }

    /**
     * Fit a smoothing spline
     * 
     * nrguments:
     * - rms_start: root mean square error start value
     * - converged: boolean convergence function, with arguments
     *   - number of knots (usize), 
     *   - number of added knots in the last iteration (usize)
     *   - root mean square error of the fit (f64),
     *   - root mean square error improvement in the last iteration
     * - rms_scale_ratio: amount to reduce the target rms value in each iteration, default = 0.8
     * - n_iter: number of iterations
     * 
     */
    pub fn smoothing_spline(mut self, 
            rms_start: f64, 
            converged: impl Fn(usize, usize, f64, f64) -> bool, 
            rms_scale_ratio: Option<f64>, 
            n_iter: Option<usize>,
        ) -> Result<Spline<K>>{
        let ratio = rms_scale_ratio.unwrap_or(0.85);
        let n_iter = n_iter.unwrap_or(20);
        let (ierr, fp)= self.concur(0, Some(rms_start), None);
        if ierr>0 {
            return Err(DierckxError(207).into())
        }
        let mut rms = fp;
        let mut n_prev;
        let mut rms_prev ;
        for _ in 0..n_iter {
            n_prev = self.n;
            rms_prev = rms;
            let (ierr, fp)= self.concur(1, Some(rms * ratio), None);
            rms = fp;
            if ierr>0 {
                return Err(DierckxError(ierr).into())
            }
            if converged(self.n, self.n-n_prev, rms, rms_prev - rms) { // finishing fit
                let (ierr, _fp)= self.concur(0, Some(rms), None);
                if ierr>0 {
                    return Err(DierckxError(ierr).into())
                } else {
                    return Ok(self.into())
                }
            };
        }
        Err(DierckxError(206).into())
    }

} // impl ParametricCurveSplineFit


impl<const K:usize, const N:usize> From<ParametricCurveSplineFit<K,N>> for Spline<K> {
    fn from(mut sp: ParametricCurveSplineFit<K,N>) -> Self {
        sp.t.truncate(sp.n);
        sp.t.shrink_to_fit();
        sp.c.truncate(sp.n*sp.idim);
        sp.c.shrink_to_fit();
        Spline::new(
            sp.t,
            sp.c,

        )
    }
}
