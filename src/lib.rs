mod dierckx;

use std::error;
use std::fmt;
use std::iter::repeat;
use crate::dierckx::{curfit_, splev_};


type Result<T> = std::result::Result<T, Box<dyn error::Error>>;

#[derive(Debug, Clone)]
pub struct DierckxError{
    pub ierr: i32,
}

impl DierckxError {
    fn new(ierr: i32) -> Self { Self { ierr } }
}

impl fmt::Display for DierckxError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.ierr {
           -2 => write!(f, "normal return for weighted least squares spline, fp upper bound for smoothing factor"),
           -1 => write!(f, "normal return for interpolating spline"),
            0 => write!(f, "normal return"),
            1 => write!(f, "error, out of storage space; nest too small (m/2); or s too small"),
            2 => write!(f, "error, smoothing spline error, s too small"),
            3 => write!(f, "error, reached iteration limit (20) for finding smoothing spline; s too small"),
           10 => write!(f, "error, invalid input data; check if -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)"),
            _ => write!(f, "unknown error"),
        }
    }
}


impl error::Error for DierckxError {}


/**
 * Spline (t,c) control points.
 */
pub struct Spline<const K:usize> {
    pub t: Vec<f64>,    // Knot values
    pub c: Vec<f64>,    // b-Spline coefficients
}

impl<const K:usize> Spline<K>  {
    pub fn new(t: Vec<f64>, c: Vec<f64>) -> Self {
        assert!(t.len()==c.len());
        Self {t, c}
    }
    
    pub fn values(&self, x: &Vec<f64>) -> Result<Vec<f64>> {
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


pub struct Dierckx<const K:usize> {
    // input values
    x: Vec<f64>,    // data x coordinates
    y: Vec<f64>,    // data y coordinates
    w: Vec<f64>,    // weight factors, 

    pub tc: Spline<K>,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries
}


impl<const K:usize> Dierckx<K> {
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

    fn curfit(&mut self, iopt:i32, e_rms_pct:Option<f64>, knots: Option<Vec<f64>>) ->  (i32, f64) {
        let k = K;
        let m = self.x.len();
        let nest = m * K  + 1;
        let lwrk = self.wrk.len();
        let mut fp = 0.0;
        let mut ierr = 0;
        let y_rms:f64;
        let s = if let Some(e) = e_rms_pct {
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
    pub fn cardinal_spline(&mut self, dt:f64) -> Result<(&Self,f64)>{
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
            Ok((self, fp))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     * Interpolating Spline
     * 
     * Knots at x values, no error: fp = s = 0.0;
     */ 
    pub fn interpolating_spline(&mut self) -> Result<&Self> {
        let (ierr, _fp) = self.curfit(0, Some(0.0),None);
        if ierr<=0  {
            Ok(self)
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

impl<const K:usize> AsRef<Spline<K>> for Dierckx<K> {
    fn as_ref(&self) -> &Spline<K> {
        &self.tc
    }
}

#[test]
fn test_smoothing() -> Result<()> {

    let xi = vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
    let yi = xi.iter().map(|x|x*x).collect();
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let mut d = Dierckx::<3>::new(xi, yi, None);
    let r = d.interpolating_spline()?;
    println!("{:.4?}", r.as_ref().values(&x));


    Ok(())
}