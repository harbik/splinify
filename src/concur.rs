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

pub struct ConstrainedSpline<const K:usize, const N:usize> {
    // input values
    xn: Vec<f64>, // data (x,y,..) coordinates
    u: Vec<f64>,
    w: Vec<f64>,    // weight factors, 
    xb: Vec<f64>,
    xe: Vec<f64>,

    pub tc: Spline<K>,

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
impl<const K:usize, const N: usize> ConstrainedSpline<K, N> {

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
        xb: Vec<f64>,
        xe: Vec<f64>,
    ) -> Result<Self> {

        let k = K;
        let idim = N;
        let m = u.len(); // number of coordinates
        if m<2 {return Err("need at least 2 parameter values".into())};
        let mx = m * N;
        if xn.len()!= mx { return Err("incorrect size of coordinate array xn".into())}
        let w = vec![1.0;m];

        let ib = xb.len()/N-1;
        let ie = xe.len()/N-1;

        let nest = m+K+1+(ib-1).max(0)+(ie-1).max(0);
        let tc = Spline::<K>::new(vec![0.0; nest], vec![0.0; nest * N]);
        let iwrk = vec![0i32; nest];

        let wrk = vec![0f64; m*(K+1)+nest*(6+N+3*K)];
        let xx = vec![0.0; N*m];
        let cp = vec![0.0; 2 * (K+1) * N];

        Ok(Self { u, xn, w, xb, xe, tc, wrk, iwrk, xx, cp, ib, ie, m, mx, nest, k, idim})

    }

    pub fn set_weights(&mut self, weights: Vec<f64>) -> Result<&mut Self> {
        if weights.len() == self.u.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err("Wrong size for weights array".into())
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
        let nc = self.tc.c.len();
        let lwrk = self.wrk.len();

        if let Some(knots) = knots {
            self.tc.t = knots;
        }
        let mut n = self.tc.t.len();
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
                &mut n,
                self.tc.t.as_mut_ptr(),
                &nc,
                self.tc.c.as_mut_ptr(),
                &np,
                self.cp.as_mut_ptr(),
                &mut fp,
                self.wrk.as_mut_ptr(),
                &lwrk,
                self.iwrk.as_mut_ptr(),
                &mut ierr,
            );
        }
        self.tc.t.truncate(n);
        self.tc.c.truncate(n);
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
        if n == 0 { return Err("Cardinal spline spacing too large: select smaller interval".into()) };
        let t: Vec<f64> = 
            repeat(tb).take(K-1)
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
                repeat(te).take(K-1)
            )
            .collect();

        let (ierr, fp) = self.concur(-1, Some(0.0),Some(t));
        if ierr<=0  {
            Ok((self.tc, fp))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     Interpolating Spline
      
      
     */ 
    pub fn interpolating_spline(mut self) -> Result<Spline<K>> {
        let (ierr, _fp) = self.concur(0, Some(0.0),None);
        if ierr<=0  {
            Ok(self.tc)
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     * Smoothing Spline
     * 
     * A spline with a minimal number of knots, with error less than the specifed rms value.
     * Repeat fit with smaller rms value using `smooth_more`.
     */
    pub fn smoothing_spline(&mut self, rms: f64) -> Result<(&Self,f64)>{
        let (ierr, fp) = self.concur(0, Some(rms), None);
        if ierr<=0  {
            Ok((self, fp))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     * Improved Smoothing Spline 
     * 
     * Repeat fit with smaller rms value after a first `smoothing_spline` attempt.
     */
    pub fn smooth_more(&mut self, rms: f64) -> Result<(&Self,f64)>{
        let (ierr, fp) = self.concur(1, Some(rms),None);
        if ierr<=0  {
            Ok((self, fp.sqrt()/self.u.len() as f64))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }
}

impl<const K:usize, const N:usize> AsRef<Spline<K>> for ConstrainedSpline<K,N> {
    fn as_ref(&self) -> &Spline<K> {
        &self.tc
    }
}
