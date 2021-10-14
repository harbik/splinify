
#![doc = include_str!("../README.md")]


/// Foreign Function Interface definitions
mod dierckx;

pub mod curfit;
pub use curfit::*;

pub mod concur;
pub use concur::*;

use std::error;
use std::fmt;
use crate::dierckx::{splev_};


pub type Result<T> = std::result::Result<T, Box<dyn error::Error>>;

#[derive(Debug, Clone)]
pub struct DierckxError(i32);

impl DierckxError {
    fn new(ierr: i32) -> Self { Self(ierr) }
}

impl fmt::Display for DierckxError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.0 {
            // Dierckx
           -2 => write!(f, "normal return for weighted least squares spline, fp upper bound for smoothing factor"),
           -1 => write!(f, "normal return for interpolating spline"),
            0 => write!(f, "normal return"),
            1 => write!(f, "out of storage space; nest too small (m/2); or s too small"),
            2 => write!(f, "smoothing spline error, s too small"),
            3 => write!(f, "reached iteration limit (20) for finding smoothing spline; s too small"),
           10 => write!(f, "invalid input data; check if -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)"),
           // this library
           200 => write!(f, "N should be between 1 and 10"),
           201 => write!(f, "need at least 2 parameter values"),
           202 => write!(f, "incorrect size of coordinate array xn"),
           203 => write!(f, "wrong size for weights array"),
           204 => write!(f, "too many derivative contraints supplied"),
           205 => write!(f, "cardinal spline spacing too large: select smaller interval"),
           206 => write!(f, "smoothing_spline not converged"),
           207 => write!(f, "failed to initialize smoothing_spline"),
            _ => write!(f, "unknown error"),
        }
    }
}


impl error::Error for DierckxError {}


/**
 * General B-Spline Curve Knot/Coefficient Representation
 */
#[derive(Clone)]
pub struct Spline<const K:usize> {
    pub t: Vec<f64>,    // Knot values
    pub c: Vec<f64>,    // b-Spline coefficients
}

impl<const K:usize> Spline<K>  {
    pub fn new(t: Vec<f64>, c: Vec<f64>) -> Self {
        assert!(t.len()==c.len());
        Self {t, c}
    }
    
    pub fn evaluate(&self, x: &Vec<f64>) -> Result<Vec<f64>> {
        let k = K;
        let m = x.len();
        let mut y = vec![0.0; m];
        let n = self.t.len();
        let mut ierr = 0;
        unsafe {splev_(self.t.as_ptr(), &n, self.c.as_ptr(), &k, x.as_ptr(), y.as_mut_ptr(), &m, &mut ierr) }
        if ierr<=0 {
            Ok(y)
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }
}


/*

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
    pub fn new(x: Vec<f64>, y: Vec<f64>, weights: Option<Vec<f64>>) -> Self {

        let m = x.len();
        let w = weights.unwrap_or(vec![1.0; m]);
        assert!(y.len()==m);
        assert!(w.len()==y.len());

        let nest = m * K  + 1;
        let tc = Spline::<K>::new(vec![0.0; nest], vec![0.0; nest]);

        let lwrk = m * (K + 1) + nest * (7 + 3 * K);
        let wrk = vec![0f64; lwrk];
        let iwrk = vec![0i32; lwrk];

        Self { x, y, w, tc, wrk, iwrk}

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

 */