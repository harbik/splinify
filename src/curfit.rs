use std::iter::repeat;
//use crate::dierckx::{curfit_};
use dierckx_sys::{curfit_};
use super::{Spline, DierckxError};
use crate::FitResult;


pub type CubicCurveFit = CurveSplineFit::<3>;
pub type QuinticCurveFit = CurveSplineFit::<5>;


pub struct CurveSplineFit<const K:usize> {
    // input values
    x: Vec<f64>,    // data x coordinates
    y: Vec<f64>,    // data y coordinates
    w: Vec<f64>,    // weight factors, 

    pub tc: Spline<K,1>,
    e_rms: Option<f64>,

    // work space values
    wrk: Vec<f64>,  // used for successive tries
    iwrk: Vec<i32>, // used for successive tries
}


impl<const K:usize> CurveSplineFit<K> {

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
        let tc = Spline::<K,1>::new(vec![0.0; nest], vec![0.0; nest]);
        let iwrk = vec![0i32; nest];

        let lwrk = m * (K + 1) + nest * (7 + 3 * K);
        let wrk = vec![0f64; lwrk];

        Self { x, y, w, tc, wrk, iwrk, e_rms:None}

    }

    pub fn set_weights(mut self, weights:Vec<f64>) -> FitResult<Self> {
        if weights.len() == self.x.len() {
            self.w = weights;
            Ok(self)
        } else {
            Err(DierckxError(203).into())
        }
    }

    fn curfit(&mut self, iopt:i32, e_rms:Option<f64>, knots: Option<Vec<f64>>) ->  i32 {
        let k = K as i32;
        let m = self.x.len() as i32;
        let nest = m * k  + 1;
        let lwrk = self.wrk.len() as i32;
        let mut fp = 0.0;
        let s = if let Some(e) = e_rms {
            m as f64 * e.powi(2)
        } else {
            0.0
        };
        let mut ierr = 0;

        //let e_rms_pct = e_rms_pct.unwrap_or(0.0); // as percentage of y_rms
        //let y_rms = (self.y.iter().map(|y|y*y).sum::<f64>()/m as f64).sqrt();
        //let s = m as f64 * (e_rms_pct * y_rms / 100.0).powi(2);
        if let Some(knots) = knots {
            self.tc.t = knots;
        }
        let mut n = self.tc.t.len() as i32;
        unsafe {
            curfit_(&iopt, &m, 
                self.x.as_ptr(), self.y.as_ptr(), self.w.as_ptr(), 
                &self.x[0], &self.x[m as usize -1], 
                &k, &s, &nest, &mut n, 
                self.tc.t.as_mut_ptr(), self.tc.c.as_mut_ptr(), 
                &mut fp, 
                self.wrk.as_mut_ptr(), &lwrk, self.iwrk.as_mut_ptr(), 
                &mut ierr
            );
        }
        self.tc.t.truncate(n as usize);
        self.tc.c.truncate(n as usize);
        self.e_rms = Some((fp/m as f64).sqrt());
        ierr
    }


    /**
     * Cardinal Spline: Weighted least squares spline with equidistant knots
     * 
     * Returns Spline, and rms error, with knots dt (input parameter) apart,
     * and aligned to integer multiples of it. Knots cover the range within
     * the bounds of x.
     */
    pub fn cardinal_spline(mut self, dt:f64) -> FitResult<Spline<K,1>>{
        let m = self.x.len();
        let tb = (self.x[0]/dt).ceil() * dt;
        let te = (self.x[m-1]/dt).floor() * dt;
        let n = ((te - tb)/dt).round() as usize;
        if n == 0 { return Err(DierckxError(205).into())};
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

        let ierr = self.curfit(-1, Some(0.0),Some(t));
        if ierr<=0  {
            Ok(self.tc)
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     Interpolating Spline
     */ 
    pub fn interpolating_spline(mut self) -> FitResult<Spline<K,1>> {
        let ierr = self.curfit(0, Some(0.0),None);
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
    pub fn smoothing_spline(mut self, rms: f64) -> FitResult<Spline<K,1>>{
        let ierr = self.curfit(0, Some(rms), None);
        if ierr<=0  {
            Ok(self.tc)
        } else {
            Err(DierckxError(ierr).into())
        }
    }
}

