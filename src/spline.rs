
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
    pub e: Option<f64>, // optional rms fit error
}

impl<const K:usize, const N:usize> Spline<K, N>  {
    pub fn new(t: Vec<f64>, c: Vec<f64>) -> Self {
        Self::with_e_rms(t, c, None)
    }

    pub fn with_e_rms(t: Vec<f64>, c: Vec<f64>, e_rms: Option<f64>) -> Self {
        assert!(c.len() == N * (t.len() - K - 1));
        Self {t, c, e: e_rms,k: K, n:N}
    }


    pub fn eval(&self, xs: &[f64]) -> Result<Vec<f64>> {
        let n = self.t.len();
        let mut v:Vec<f64> = Vec::with_capacity(xs.len());

        let mut i = self. k;
        let mut x_prev = f64::NEG_INFINITY;

        for &x in xs {
            if x<=x_prev {
                return Err("x values should be sorted in strict increasing order".into());
            } else {
                x_prev = x;
            };
            if  x<=self.t[self.k] || x>=self.t[n-self.k+1] {
                v.push(x.clamp(self.t[self.k], self.t[n-K+1]))
            } else {
                while i<n-1 &&self.t[i+1]<=x {
                    i+=1
                };
                v.push(self.deboor(i, x))
            }
        };
        Ok(v)
    }

    pub fn deboor(&self, i: usize,  x: f64) -> f64 {
        let mut d = vec![0.0; self.k + 1];

        for j in 0..K+1 {
            d[j] = self.c[(j + i - self.k)];
        }

        for r in 1..self.k+1 {
            for j in (r..=self.k).into_iter().rev() {
                let alpha = (x - self.t[j + i - self.k]) / (self.t[j + 1 + i - r] - self.t[j + i - self.k]);
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]
            }
        }
        d[self.k]
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

    #[test]
    fn test_eval() -> Result<()> {
        let s: Spline<3,1> = Spline::new(
            vec![-2.0, -2.0, -2.0, -2.0, -1.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0],
            vec![0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, /*check extras*/0.0, 0.0, 0.0, 0.0]
        );
        println!("{:.3?}", s.eval(&[-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0])?);
        Ok(())
    }