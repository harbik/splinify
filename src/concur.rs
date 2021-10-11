use std::iter::repeat;
use crate::dierckx::{concur_};
use super::{Spline, DierckxError};
use nalgebra::{OMatrix, DVector, Dynamic, Const};

use crate::Result;

pub struct ConstrainedSpline<const K:usize, const N:usize> {
    // input values
    xn: OMatrix<f64, Const<N>, Dynamic>, // data (x,y,..) coordinates
    u: DVector<f64>,
    w: DVector<f64>,    // weight factors, 
    xb: OMatrix<f64, Dynamic, Const<N>>,
    xe: OMatrix<f64, Dynamic, Const<N>>,

    pub tc: Spline<K>,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries
    xx: Vec<f64>,
    cp: Vec<f64>,
}


/**
  
 A ConstrainedSplineCurve builder, using Dierckx' `concur` subroutine.

 */
impl<const K:usize, const N: usize> ConstrainedSpline<K, N> {

    /**
    Constructor, taking as input cuve parameter, curve coordinates, and end point constraints,
    and setting up remaining datastructures and values for `concur`.

    The matrices are represented by the `nalgebra::OMatrix` owned matrix type.  
    */
    pub fn new(
        u: DVector<f64>,
        xn: OMatrix<f64, Const<N>, Dynamic>, 
        xb: OMatrix<f64, Dynamic, Const<N>>, 
        xe: OMatrix<f64, Dynamic, Const<N>>,
    ) -> Self {

        let m = xn.ncols();
        let w = DVector::repeat(m, 1.0);

        let ib = xb.ncols();
        let ie = xe.ncols();

        let nest = m+K+1+(ib-1).max(0)+(ie-1).max(0);
        let tc = Spline::<K>::new(vec![0.0; nest], vec![0.0; nest * N]);
        let iwrk = vec![0i32; nest];

        let lwrk =  m*(K+1)+nest*(6+N+3*K);
        let wrk = vec![0f64; lwrk];
        let xx = vec![0.0; N*m];
        let cp = vec![0.0; 2 * (K+1) * N];

        Self { u, xn, w, xb, xe, tc, wrk, iwrk, xx, cp}

    }

    pub fn set_weights(&mut self, weights: DVector<f64>) -> Result<&mut Self> {
        if weights.len() == self.u.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err("Wrong size for weights array".into())
        }
    }

    fn concur(&mut self, iopt:i32, e_rms:Option<f64>, knots: Option<Vec<f64>>) ->  (i32, f64) {
        let k = K;
        let m = self.xn.ncols();
        let mx = N * m;
        let idim = N;

        let nest = m * K  + 1;
        let lwrk = self.wrk.len();
        let mut fp = 0.0;
        let mut ierr = 0;
        let s = if let Some(e) = e_rms {
            m as f64 * e.powi(2)
        } else {
            0.0
        };

        let ib = self.xb.ncols();
        let nb = self.xb.len();
        let ie = self.xe.ncols();
        let ne = self.xe.len();
        let np = self.cp.len();

        if let Some(knots) = knots {
            self.tc.t = knots;
        }
        let mut n = self.tc.t.len();
        unsafe {
            concur_(
                &iopt,
                &idim,
                &m,
                self.u.as_ptr(),
                &mx,
                self.xn.as_ptr(),
                self.xx.as_mut_ptr(), 
                self.w.as_ptr(),
                &ib,
                self.xb.as_ptr(),
                &nb,
                &ie,
                self.xe.as_ptr(),
                &ne,
                &k,
                &s,
                &nest,
                &mut n,
                self.tc.t.as_mut_ptr(),
                &self.tc.c.len(),
                self.tc.c.as_mut_ptr(),
                self.cp.as_mut_ptr(),
                &np,
                &mut fp,
                self.wrk.as_mut_ptr(),
                &lwrk,
                self.iwrk.as_mut_ptr(),
                &mut ierr,
            );
        }
        self.tc.t.truncate(n);
        self.tc.c.truncate(n);
        (ierr, (fp/m as f64).sqrt()) // fit error, in space coordinates
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
        let n = n + 2; 
        let t: Vec<f64> = 
            repeat(tb).take(K)
            .chain(
                repeat(dt).take(n).scan(tb, 
                    |s, dx|{
                        let t=*s; 
                        *s+=dx; 
                        Some(t)
                    })
            )
            .chain(
                repeat(te).take(K)
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
