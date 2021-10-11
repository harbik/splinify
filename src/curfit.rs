use std::iter::repeat;
use crate::dierckx::{curfit_};
use super::{Spline, DierckxError};
use crate::Result;


pub type CubicCurveFit = CurveFit::<3>;
pub type QuinticCurveFit = CurveFit::<5>;


pub struct CurveFit<const K:usize> {
    // input values
    x: Vec<f64>,    // data x coordinates
    y: Vec<f64>,    // data y coordinates
    w: Vec<f64>,    // weight factors, 

    pub tc: Spline<K>,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries
}


/**
 * 
 Fits a smooth B-Spline 1D curve from (x,y) data
 
 It is intended to fit a 1D Spline curve through a set of points,
 given as x and y vectors, and with the x values in strictly ascending order.
 As an example, the values could be temperature data, measured over a period of time,
 with the temperature sensor having a limited accuaracy.
 
 This is a Rust wrapper for Dierckx curfit Fortran subroutine, part
 of Dierckx FITPACK library.
 For the Foreign Function Interface to the Fortan subroutnine
 see [curfit][dierckx::curfit_];
 as you can see, it is not easy to use as it has 18(!) arguments, and requires quite a bit of reading and experimenting
 to grasp.
  
 This wrapper encapsulates the internal data, and breaks up the curfit call in multiple steps.

 */
impl<const K:usize> CurveFit<K> {

    /**
     Constructor, with inputs x and y vectors, and an optional weights vectors.

     The vectors should have equal length.
     */
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {

        let m = x.len();
        let w = vec![1.0; m];
        assert!(y.len()==m);
        assert!(w.len()==y.len());

        let nest = m * K  + 1;
        let tc = Spline::<K>::new(vec![0.0; nest], vec![0.0; nest]);
        let iwrk = vec![0i32; nest];

        let lwrk = m * (K + 1) + nest * (7 + 3 * K);
        let wrk = vec![0f64; lwrk];

        Self { x, y, w, tc, wrk, iwrk}

    }

    pub fn set_weights(&mut self, weights:Vec<f64>) -> Result<&mut Self> {
        if weights.len() == self.x.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err("Wrong size for weights array".into())
        }
    }

    fn curfit(&mut self, iopt:i32, noise_pct:Option<f64>, knots: Option<Vec<f64>>) ->  (i32, f64) {
        let k = K;
        let m = self.x.len();
        let nest = m * K  + 1;
        let lwrk = self.wrk.len();
        let mut fp = 0.0;
        let mut ierr = 0;
        let y_rms:f64;
        let s = if let Some(e) = noise_pct {
            y_rms = (self.y.iter().map(|y|y*y).sum::<f64>()/m as f64).sqrt();
            m as f64 * (e * y_rms / 100.0).powi(2)
        } else {
            y_rms = f64::MAX; // to force fp output to 0.0
            0.0
        };
        //let e_rms_pct = e_rms_pct.unwrap_or(0.0); // as percentage of y_rms
        //let y_rms = (self.y.iter().map(|y|y*y).sum::<f64>()/m as f64).sqrt();
        //let s = m as f64 * (e_rms_pct * y_rms / 100.0).powi(2);
        if let Some(knots) = knots {
            self.tc.t = knots;
        }
        let mut n = self.tc.t.len();
        unsafe {
            curfit_(&iopt, &m, 
                self.x.as_ptr(), self.y.as_ptr(), self.w.as_ptr(), 
                &self.x[0], &self.x[m-1], 
                &k, &s, &nest, &mut n, 
                self.tc.t.as_mut_ptr(), self.tc.c.as_mut_ptr(), 
                &mut fp, 
                self.wrk.as_mut_ptr(), &lwrk, self.iwrk.as_mut_ptr(), 
                &mut ierr
            );
        }
        self.tc.t.truncate(n);
        self.tc.c.truncate(n);
        (ierr, (fp/m as f64).sqrt() * 100.0/ y_rms) // fit error as percent rms
    }


    /**
     * Cardinal Spline: Weighted least squares spline with equidistant knots
     * 
     * Returns Spline, and rms error, with knots dt (input parameter) apart,
     * and aligned to integer multiples of it. Knots cover the range within
     * the bounds of x.
     */
    pub fn cardinal_spline(mut self, dt:f64) -> Result<(Spline<K>,f64)>{
        let m = self.x.len();
        let tb = (self.x[0]/dt).ceil() * dt;
        let te = (self.x[m-1]/dt).floor() * dt;
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

        let (ierr, fp) = self.curfit(-1, Some(0.0),Some(t));
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
        let (ierr, _fp) = self.curfit(0, Some(0.0),None);
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
        let (ierr, fp) = self.curfit(0, Some(rms), None);
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
        let (ierr, fp) = self.curfit(1, Some(rms),None);
        if ierr<=0  {
            Ok((self, fp.sqrt()/self.x.len() as f64))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }
}

impl<const K:usize> AsRef<Spline<K>> for CurveFit<K> {
    fn as_ref(&self) -> &Spline<K> {
        &self.tc
    }
}
