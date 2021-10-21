#![doc = include_str!("../README.md")]


/// Foreign Function Interface definitions
//mod dierckx;

pub mod curfit;
pub use curfit::*;

pub mod concur;
pub use concur::*;

pub mod spline;
pub use spline::*;

pub mod util;
pub use util::*;

use std::error;
use std::fmt;
//use crate::dierckx::{splev_};
//use dierckx_sys::{splev_};


pub type Result<T> = std::result::Result<T, Box<dyn error::Error>>;

// Single Output Spline Fit
pub type LinearSplineFit = CurveSplineFit<1>;
pub type CubicSplineFit = CurveSplineFit<3>;
pub type QuinticSplineFit = CurveSplineFit<5>;


// Multi-Output Parametrized Curve Fits
pub type LinearSplineFit1D = ParameterCurveSplineFit<1,1>;
pub type CubicSplineFit1D = ParameterCurveSplineFit<3,1>;
pub type QuinticSplineFit1D = ParameterCurveSplineFit<5,1>;
pub type LinearSplineFit2D = ParameterCurveSplineFit<1,2>;
pub type CubicSplineFit2D = ParameterCurveSplineFit<3,2>;
pub type QuinticSplineFit2D = ParameterCurveSplineFit<5,2>;
pub type LinearSplineFit3D = ParameterCurveSplineFit<1,3>;
pub type CubicSplineFit3D = ParameterCurveSplineFit<3,3>;
pub type QuinticSplineFit3D = ParameterCurveSplineFit<5,3>;

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
           208 => write!(f, "K should be 1, 3 or 5"),
            _ => write!(f, "unknown error"),
        }
    }
}


impl error::Error for DierckxError {}

