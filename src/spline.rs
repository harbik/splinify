
use super::{DierckxError, Result};
use dierckx_sys::{splev_, curev_};
use serde::{Deserialize, Serialize};

/**
 * General B-Spline Curve Knot/Coefficient Representation
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct Spline<const K:usize, const N:usize> {
    pub t: Vec<f64>,        // Knot values
    pub c: Vec<f64>,        // b-Spline coefficients
    k: usize,               // Spline degree
    n: usize,               // Spline dimension
    #[serde(skip)]
    pub e_rms: Option<f64>, // optional rms fit error
}

impl<const K:usize, const N:usize> Spline<K, N>  {
    pub fn new(t: Vec<f64>, c: Vec<f64>) -> Self {
        Self::with_e_rms(t, c, None)
    }

    pub fn with_e_rms(t: Vec<f64>, c: Vec<f64>, e_rms: Option<f64>) -> Self {
        assert!(N*t.len()==c.len());
        Self {t, c, e_rms,k: K, n:N}
    }
    
    pub fn evaluate(&self, x: &Vec<f64>) -> Result<Vec<f64>> {
        let (ierr, y)  = 
            match N {
                1  => self.splev(x),
                _ => self.curev(x),
            };
        if ierr<=0 {
            Ok(y)
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    fn splev(&self, x: &Vec<f64>) -> (i32, Vec<f64>) {
        let k = K as i32;
        let m = x.len() as i32;
        let mut y = vec![0.0; m as usize];
        let n = self.t.len() as i32;
        let mut ierr = 0;
        unsafe {
            splev_(
                self.t.as_ptr(), 
                &n, 
                self.c.as_ptr(), 
                &k, 
                x.as_ptr(), 
                y.as_mut_ptr(), 
                &m, 
                &mut ierr
            );
        }
        (ierr, y)
    }

    fn curev(&self, u: &Vec<f64>) -> (i32, Vec<f64>) {
        let k = K as i32;
        let idim = N as i32;
        let m = u.len() as i32;
        let mxy = m * idim;
        let mut xy = vec![0.0; mxy as usize];
        let n = self.t.len() as i32;
        let nc = self.c.len() as i32;
        let mut ierr = 0;
        unsafe {
            curev_(
                &idim,
                self.t.as_ptr(), 
                &n, 
                self.c.as_ptr(), 
                &nc,
                &k, 
                u.as_ptr(), 
                &m, 
                xy.as_mut_ptr(), 
                &mxy,
                &mut ierr
            );
        }
        (ierr, xy)
    }

}
