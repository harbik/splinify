use std::error;
use std::fmt;

extern {

    /**

    If iopt>=0 the number of knots of the spline s(x) and the position t(j),j=1,2,...,n is chosen automatically by the
    routine.  The smoothness of s(x) is then achieved by minimalizing the discontinuity jumps of the k<sup>th</sup>
    derivative of s(x) at the knots t(j),j=k+2,k+3,..., n-k-1. 

    The amount of smoothness is determined by the condition that f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a
    given non-negative constant, called the smoothing factor.

    By means of the parameter s, the user can control the tradeoff between closeness of fit and smoothness of fit of the
    approximation.  If s is too large, the spline will be too smooth and signal will be lost; if s is too small the
    spline will pick up too much noise.  In the extreme cases the program will return an interpolating spline if s=0 and
    the weighted least-squares polynomial of degree k if s is very large.  Between these extremes, a properly chosen s
    will result in a good compromise between closeness of fit and smoothness of fit. 

    To decide whether an approximation, corresponding to a certain s, is satisfactory, the user is highly recommended to
    inspect the fits graphically.  Recommended values for s depend on the weights w(i).  If these are taken as d(i) with
    d(i) an estimate of the standard deviation of y(i), a good s-value should be found in the range (m-sqrt(2m),m+
    sqrt(2m)).  If nothing is known about the statistical error in y(i) each w(i) can be set equal to one and s
    determined by trial and error, taking account of the comments above.  The best is then to start with a very large
    value of s ( to determine the least-squares polynomial and the corresponding upper bound fp0 for s) and then to
    progressively decrease the value of s (say by a factor 10 in the beginning, i.e. s=fp010, fp0100,...and more
    carefully as the approximation shows more detail) to obtain closer fits. 

    To economize the search for a good s-value the program provides with different modes of computation:
    -  iopt = 1: At the first call of the routine, or whenever he wants to restart with the initial set of knots the user must set
    iopt=0.  
    - If iopt=1 the program will continue with the set of knots found at the last call of the routine.  This will save a
    lot of computation time if curfit is called repeatedly for different values of s.  The number of knots of the spline
    returned and their location will depend on the value of s and on the complexity of the shape of the function
    underlying the data.  But, if the computation mode iopt=1 is used, the knots returned may also depend on the
    s-values at previous calls (if these were smaller). Therefore, if after a number of trials with different s-values
    and iopt=1, the user can finally accept a fit as satisfactory, it may be worthwhile for him to call curfit once more
    with the selected value for s but now with iopt=0. 

    */
    fn curfit(
        iopt: &i32,     // iopt -1: Least-squares spline fixed knots, 0,1: smoothing spline. iopt=0 and s=0: interpolating spline
        m: &usize,      // Number of data points supplied
        x: *const f64,  // Array of x coordinates (at least m values)
        y: *const f64,  // Array of y coordinates (at least m values)
        w: *const f64,  // Array weights (at least m values)
        xb: &f64,       // Bounderies of the approximation interval. xb<=x(1), xe>=x(m)
        xe: &f64, 
        k: &usize,      // Degree of the spline, Cubic = 3
        s: &f64,        // Smoothing factor to be used if iopt >= 0
        nest: &usize,   // nest = m + k + 1
        n: &mut usize,  // Number of knots returned. For iopt=-1 value needs to pe specified on entry
        t: *mut f64,    // Array of dimension of at least nest. For iopt=-1 array of knots to be used for lsq spline
        c: *mut f64,    // Double array of at least nest. Will contain the coefficients of the b-spline representation
        fp: &mut f64,   // Weighted sum of the squared residuals of the spline approximation.
        wrk: *mut f64,  // Double array of dimension at least (m(k+1)+nest(7+3k)).
        lwrk: &usize,   // Size of 'wrk'
        iwrk: *mut i32, // int Array of at least nest (m + k + 1)
        ier: &mut i32   // Error flag.
    );

    fn splev(
        t: *const f64,  // array,length n, which contains the position of the knots
        n: &usize,      // integer, giving the total number of knots of s(x). 
        c: *const f64,  // array,length n, which contains the b-spline coefficients
        k: &usize,      // integer, giving the degree of s(x)
        x: *const f64,  // array,length m, which contains the points where s(x) must be evaluated
        y: *mut f64,    // array,length m, giving the value of s(x) at the different points
        m: &usize,      // lenght of x and y
        ier: &mut i32,  // ier = 0 : normal return;  ier =10 : invalid input data : restrictions:  m >= 1, t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1
    ); 
}

pub type CubicSpline = Dierckx<3>;
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


pub struct Spline<const K:usize> {
    pub t: Vec<f64>,    // Knot values
    pub c: Vec<f64>,    // b-Spline coefficients
}




pub struct Dierckx<const K:usize> {
    // input values
    x: Vec<f64>,    // data x coordinates
    y: Vec<f64>,    // data y coordinates
    w: Vec<f64>,    // weight factors, 

    // output
    t: Vec<f64>,    // Knot values
    c: Vec<f64>,    // b-Spline coefficients
                    // (t,c) coordinates can 'control points'

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
        let t = vec![0.0; nest];
        let c = vec![0.0; nest];

        let lwrk = m * (K + 1) + nest * (7 + 3 * K);
        let wrk = vec![0f64; lwrk];
        let iwrk = vec![0i32; lwrk];

        Self { x, y, w, t, c, wrk, iwrk}

    }

    fn curfit(&mut self, iopt:i32, rms:Option<f64>, knots: Option<Vec<f64>>) ->  (i32, f64) {
        let k = K;
        let m = self.x.len();
        let nest = m * K  + 1;
        let lwrk = self.wrk.len();
        let mut fp = 0.0;
        let mut ierr = 0;
        let s = if iopt == 0 {
            let e_r = rms.unwrap_or({
                //  1% average deviation as starting point
                 let y_avg_sum = self.y.iter().sum::<f64>();
                 y_avg_sum/m as f64/100.0
                });
            m as f64 * e_r.powi(2) // curfit needs the integral square error
        } else {
            // not used for least squares spline
            0.0
        };
        if let Some(knots) = knots {
            self.t = knots;
        }
        let mut n = self.t.len();
        unsafe {
            curfit(&iopt, &m, 
                self.x.as_ptr(), self.y.as_ptr(), self.w.as_ptr(), 
                &self.x[0], &self.x[m-1], 
                &k, 
                &s, 
                &nest, 
                &mut n, 
                self.t.as_mut_ptr(),
                self.c.as_mut_ptr(), 
                &mut fp, 
                self.wrk.as_mut_ptr(), 
                &lwrk, 
                self.iwrk.as_mut_ptr(), 
                &mut ierr
            );
        }
        self.t.truncate(n);
        self.c.truncate(n);
        (ierr, fp.sqrt()/m as f64)
    }


    /**
     * Cardinal Spline: Weighted least squares spline with equidistant knots
     * 
     * Returns Spline, and rms error
     */
    pub fn cardinal_spline(&mut self, tb:f64, dt:f64, n:usize) -> Result<(Spline<K>,f64)>{
        let t: Vec<f64> = (0..n).map(|i|tb+(i as f64)*dt).collect();
        let (ierr, fp) = self.curfit(-1, Some(0.0),Some(t));
        if ierr<=0  {
            Ok((Spline::<K>{t: self.t.clone(), c: self.c.clone()}, fp.sqrt()/self.x.len() as f64))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     * Interpolating Spline
     * 
     * Knots at x,y values, no error: fp = s = 0.0;
     */ 
    pub fn interpolating_spline(&mut self) -> Result<Spline<K>> {
        let (ierr, _fp) = self.curfit(0, Some(0.0),None);
        if ierr<=0  {
            Ok(Spline::<K>{t: self.t.clone(), c: self.c.clone()})
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     * Smoothing Spline
     * 
     * A spline with a minimal number of knots, with error less than the specifed rms value.
     * If no rms value is given, an rms is guessed by setting it to 1% of the average y value.
     * Repeat fit with smaller rms value using `smooth_more`.
     */
    pub fn smoothing_spline(&mut self, rms: Option<f64>) -> Result<(Spline<K>,f64)>{
        let (ierr, fp) = self.curfit(0, rms,None);
        if ierr<=0  {
            Ok((Spline::<K>{t: self.t.clone(), c: self.c.clone()}, fp.sqrt()/self.x.len() as f64))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    /**
     * Improved Smoothing Spline 
     * 
     * Repeat fit with smaller rms value after a first `smoothing_spline` attempt.
     */
    pub fn smooth_more(&mut self, rms: f64) -> Result<(Spline<K>,f64)>{
        let (ierr, fp) = self.curfit(1, Some(rms),None);
        if ierr<=0  {
            Ok((Spline::<K>{t: self.t.clone(), c: self.c.clone()}, fp.sqrt()/self.x.len() as f64))
        } else {
            Err(DierckxError::new(ierr).into())
        }
    }

    fn spline_base(x: Vec<f64>, y: Vec<f64>, weights: Option<Vec<f64>>, rms_error: Option<f64>, knots: Option<Vec<f64>>) -> Self {
        let k = K;
        let iopt = if knots.is_some() {
            -1 // least squares spline
        } else {
            0 // smoothing splnie
        };

        let m = x.len();
        assert!(m == y.len());
        let w = weights.unwrap_or(vec![1.0; m]);
        assert!(w.len() == m);

        let xb = x[0];
        let xe = x[m-1];

        let s = if iopt == 0 {
            let e_r = rms_error.unwrap_or({
                //  1% average deviation as starting point
                 let y_avg_sum = y.iter().sum::<f64>();
                 y_avg_sum/m as f64/100.0
                });
            m as f64 * e_r.powi(2) // curfit needs the integral square error
        } else {
            // not used for least squares spline
            0.0
        };
        let nest = m * k  + 1;

        let mut t = knots.unwrap_or(vec![0.0; nest]);
        let mut c = vec![0.0; nest];
        let mut n = t.len();

        let lwrk = m * (k + 1) + nest * (7 + 3 * k);
        let mut wrk = vec![0f64; lwrk as usize];
        let mut iwrk = vec![0i32; lwrk as usize];
        let mut ierr = 0;

        let mut fp = 0.0;

        unsafe {
            curfit(&iopt, 
                &m, 
                x.as_ptr(), 
                y.as_ptr(), 
                w.as_ptr(), 
                &xb, 
                &xe, 
                &k, 
                &s, 
                &nest, 
                &mut n, 
                t.as_mut_ptr(),
                c.as_mut_ptr(), 
                &mut fp, 
                wrk.as_mut_ptr(), 
                &lwrk, 
                iwrk.as_mut_ptr(), 
                &mut ierr
            );
        }
        t.truncate(n); // opt 0, 1, no change for opt -1, as than t has fixed lenght
        c.truncate(n);
        Self {
            x,
            y,
            w,
            t,
            c,
            wrk,
            iwrk,
        }

    }

    /**
     * Create a smoothing spline, by allowing the spline to deviate from the data points by a root mean square
     * deviation `rms_error`. If no `rms_error` is given, an estimate is used, being a 1% RMS relative to the 
     * average of y. The number and location of knots is given by the algorithm, and determined by the value 
     * of `rms_error`: for an error of 0.0, a purely interpolating spline is obtained, with the knots at the given
     * x values.
     * 
     */
    pub fn new_smoothing_spline(x: Vec<f64>, y: Vec<f64>, weights: Option<Vec<f64>>, rms_error: Option<f64>) -> Self {
        Self::spline_base(x, y, weights, rms_error, None)
    }

    /**
     * Create a interpolating spline, fixed at the input points (x,y). 
     * 
     */
    pub fn new_interpolating_spline(x: Vec<f64>, y: Vec<f64>, weights: Option<Vec<f64>>) -> Self {
        Self::spline_base(x, y, weights, Some(0.0), None)
    }

    /**
     * Create a least-squares-spline, for a given number of knots. Useful for heavily oversampled data:
     * the number of knots to use is determined by the resolution.
     */
    pub fn new_lsq_spline(x: Vec<f64>, y: Vec<f64>, knots: Vec<f64>, weights: Option<Vec<f64>>) -> Self {
        Self::spline_base(x,y, weights, None, Some(knots))
    }

    pub fn values(&self, x: &Vec<f64>) -> Vec<f64> {
        let k = K;
        let m = x.len();
        let mut y = vec![0.0; m];
        let n = self.t.len();
        let mut ierr = 0;
        unsafe {
            splev(
                self.t.as_ptr(), 
                &n, 
                self.c.as_ptr(), 
                &k, 
                x.as_ptr(), 
                y.as_mut_ptr(), 
                &m, 
                &mut ierr)
        }
        y
    }
}



#[test]
fn test_smoothing() -> Result<()> {
    let xi = vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
    let yi = xi.iter().map(|x|x*x).collect();
    let d = CubicSpline::new_interpolating_spline(xi, yi, None);
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let y = d.values(&x);
    println!("{:.4?}", y);

    let xi = vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
    let yi = xi.iter().map(|x|x*x).collect();
    let mut d = Dierckx::<3>::new(xi, yi, None);
    let r = d.interpolating_spline()?;
    println!("{:.4?}", r.values());

    Ok(())
}