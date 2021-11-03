
use super::{FitError, Result};
use spliny::SplineCurve;
use dierckx_sys::{splev_, curev_};


pub fn evaluate<const K: usize, const N: usize>(s: &SplineCurve<K,N>, x: &[f64]) -> Result<Vec<f64>> {
    let (ierr, y)  = 
        match N {
            1  => splev(s, x),
            _ => curev(s, x),
        };
    if ierr<=0 {
        Ok(y)
    } else {
        Err(FitError::new(ierr).into())
    }
}

fn splev<const K: usize, const N: usize>(s: &SplineCurve<K,N>, x: &[f64]) -> (i32, Vec<f64>) {
    let k = K as i32;
    let m = x.len() as i32;
    let mut y_v = vec![0.0; m as usize];
    let n = s.t.len() as i32;
    let mut ierr = 0;
    unsafe {
        splev_(
            s.t.as_ptr(), 
            &n, 
            s.c.as_ptr(), 
            &k, 
            x.as_ptr(), 
            y_v.as_mut_ptr(), 
            &m, 
            &mut ierr
        );
    }
    (ierr, y_v)
}

fn curev<const K: usize, const N: usize>(s: &SplineCurve<K,N>, u: &[f64]) -> (i32, Vec<f64>) {
    let k = K as i32;
    let idim = N as i32;
    let m = u.len() as i32;
    let mxy = m * idim;
    let mut xy = vec![0.0; mxy as usize];
    let n = s.t.len() as i32;
    let nc = s.c.len() as i32;
    let mut ierr = 0;
    unsafe {
        curev_(
            &idim,
            s.t.as_ptr(), 
            &n, 
            s.c.as_ptr(), 
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
