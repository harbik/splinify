
use super::{FitError, Result};
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
        assert!(c.len() >= N * (t.len() - K - 1));
        Self {t, c, e: e_rms,k: K, n:N}
    }


    pub fn eval(&self, xs: &[f64]) -> Result<Vec<f64>> {
        let n = self.t.len();
        let mut v:Vec<f64> = Vec::with_capacity(xs.len());

        let mut i = self. k;
        let mut x_prev = f64::NEG_INFINITY;
        let mut d = [0.0;6]; // want to use K+1 here, but currently not allowed yet by the compiler

        for &x in xs {

            if x<=x_prev {
                return Err("x values should be sorted in strict increasing order".into());
            } else {
                x_prev = x;
            };

            // clamp x to interval tb..=te
            let arg = if  x<self.t[self.k] || x>self.t[n-self.k - 1] {
                x.clamp(self.t[self.k], self.t[n-K -1])
            } else {
                x
            };

            // find knot interval which contains x=arg
            while !(arg>=self.t[i] && arg<=self.t[i+1]) {
                i+=1
            };

            for (j, dm) in d.iter_mut().enumerate().take(K+1) {
                *dm = self.c[(j + i - self.k)];
            };
            v.push(self.deboor(i, arg, &mut d))

        };
        Ok(v)
    }

    pub fn eval_n(&self, us: &[f64]) -> Result<Vec<f64>> {
        let n = self.t.len();
        let nc = self.c.len()/N;
        let mut v:Vec<f64> = Vec::with_capacity(us.len() * N); // x,y,..x,y coordinates

        let mut i = self. k;
        let mut u_prev = f64::NEG_INFINITY;
        let mut d = [0.0;6]; // want to use K+1 here, but currently not allowed yet by the compiler

        for &u in us {

            if u<=u_prev {
                return Err("x values should be sorted in strict increasing order".into());
            } else {
                u_prev = u;
            };

            // clamp x to interval tb..=te
            let arg = if  u<self.t[self.k] || u>self.t[n-self.k - 1] {
                u.clamp(self.t[self.k], self.t[n-K -1])
            } else {
                u
            };

            // find knot interval which contains x=arg
            while !(arg>=self.t[i] && arg<=self.t[i+1]) {
                i+=1
            };

            // calculate spline values
            for dim in 0..N {
                // copy relevant c values into d
                for (j, dm) in d.iter_mut().enumerate().take(K+1) {
                    *dm = self.c[dim * nc + j + i - self.k];
                };

                v.push(self.deboor(i, arg, &mut d))
            }


        };
        Ok(v)
    }

    pub fn deboor(&self, i: usize,  x: f64, d: &mut [f64;6]) -> f64 {

        /*
        for j in 0..K+1 {
            d[j] = self.c[(j + i - self.k)];
        }
        */

        for r in 1..self.k+1 {
            for j in (r..=self.k).into_iter().rev() {
                let alpha = (x - self.t[j + i - self.k]) / (self.t[j + 1 + i - r] - self.t[j + i - self.k]);
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]
            }
        }
        d[self.k]
    }

    
    pub fn evaluate(&self, x: &[f64]) -> Result<Vec<f64>> {
        let (ierr, y)  = 
            match N {
                1  => self.splev(x),
                _ => self.curev(x),
            };
        if ierr<=0 {
            Ok(y)
        } else {
            Err(FitError::new(ierr).into())
        }
    }

    fn splev(&self, x: &[f64]) -> (i32, Vec<f64>) {
        let k = K as i32;
        let m = x.len() as i32;
        let mut y_v = vec![0.0; m as usize];
        let n = self.t.len() as i32;
        let mut ierr = 0;
        unsafe {
            splev_(
                self.t.as_ptr(), 
                &n, 
                self.c.as_ptr(), 
                &k, 
                x.as_ptr(), 
                y_v.as_mut_ptr(), 
                &m, 
                &mut ierr
            );
        }
        (ierr, y_v)
    }

    fn curev(&self, u: &[f64]) -> (i32, Vec<f64>) {
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

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use super::Spline;

    // spline test values from https://docs.rs/bspline/1.0.0/bspline/index.html crate 

    #[test]
    fn linear_bspline() {
        let x = vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
        let y = vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
        
        let s: Spline<1,1> = Spline::new(
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 1.0],
        );
        let yt = s.eval_n(&x).unwrap();
        //println!("{:?}", yt);
        y.iter().zip(yt.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-8));

        let ytx = s.evaluate(&x).unwrap();
        //println!("{:?}", ytx);
        y.iter().zip(ytx.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-8));
    }
    #[test]
    fn quadratic_bspline() {
        let x = [0.0, 0.5, 1.0, 1.4, 1.5, 1.6, 2.0, 2.5, 3.0];
        let y = [0.0, 0.125, 0.5, 0.74, 0.75, 0.74, 0.5, 0.125, 0.0];

        let s: Spline<2,1> = Spline::new(
        vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0],
        vec![0.0, 0.0, 1.0, 0.0, 0.0],
        );
        let yt = s.eval(&x).unwrap();
    //  println!("{:?}", yt);
        y.iter().zip(yt.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-8));

        let ytx = s.evaluate(&x).unwrap();
    //  println!("{:?}", ytx);
        y.iter().zip(ytx.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-8));
    }

    #[test]
    fn cubic_bspline() {
        let x = vec![-2.0, -1.5, -1.0, -0.6, 0.0, 0.5, 1.5, 2.0];
        let y = vec![ 0.0, 0.125, 1.0, 2.488, 4.0, 2.875, 0.12500001, 0.0];
        let s: Spline<3,1> = Spline::new(
        vec![-2.0, -2.0, -2.0, -2.0, -1.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0],
        vec![0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0],
        );
        let yt = s.eval(&x).unwrap();
        //println!("{:?}", yt);
        y.iter().zip(yt.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-7));

        let ytx = s.evaluate(&x).unwrap();
        //println!("{:?}", ytx);
        y.iter().zip(ytx.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-7));
    }
    #[test]
    fn quartic_bspline() {
        let x = vec![0.0, 0.4, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.1, 4.5, 5.0];
        let y = vec![ 0.0,  0.0010666668,  0.041666668,  0.19791667,  0.4583333,  0.5989583,  0.4583333,  0.35206667,  0.02733751,  0.002604167,  0.0];
        let s: Spline<4,1> = Spline::new(
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0],
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        );
        let yt = s.eval(&x).unwrap();
        //println!("{:?}", yt);
        y.iter().zip(yt.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-7));

        let ytx = s.evaluate(&x).unwrap();
        //println!("{:?}", ytx);
        y.iter().zip(ytx.iter()).for_each(|(&a,&b)| assert_abs_diff_eq!(a,b,epsilon=1E-7));
    }

}