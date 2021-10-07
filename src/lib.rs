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
    fn curfit_(
        iopt: &i32, 
        m: &usize,  // nr of data points
        x: *const f64, 
        y: *const f64, 
        w: *const f64, 
        xb: &f64, 
        xe: &f64, 
        k: &usize, 
        s: &f64, 
        nest: &usize, 
        n: &mut usize, 
        t: *mut f64, 
        c: *mut f64, 
        fp: &mut f64, 
        wrk: *mut f64, 
        lwrk: &usize, 
        iwrk: *mut i32, 
        ier: &mut i32
    );

    fn splev_(
        t: *const f64, // array,length n, which contains the position of the knots
        n: &usize,       // integer, giving the total number of knots of s(x). 
        c: *const f64, // array,length n, which contains the b-spline coefficients
        k: &usize,       // integer, giving the degree of s(x)
        x: *const f64, // array,length m, which contains the points where s(x) must be evaluated
        y: *mut f64,   // array,length m, giving the value of s(x) at the different points
        m: &usize,       // lenght of x and y
        ier: &mut i32, // ier = 0 : normal return;  ier =10 : invalid input data 
        // restrictions:  m >= 1, t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1
    ); 
}

pub struct Spline {
    n: usize,       // total number fo knots
    t: Vec<f64>,    // Position of the knots
    c: Vec<f64>,    // b-Spline coefficients
    k: usize,       // Spline degree
    fp: f64,        // weighted sum of the square residuals
}


impl Spline {
    fn spline_base(x: Vec<f64>, y: Vec<f64>, weights: Option<Vec<f64>>, rms_error: Option<f64>, knots: Option<Vec<f64>>) -> Self {
        let iopt = if knots.is_some() {
            -1 // least squares spline
        } else {
            0 // smoothing splnie
        };
        let m = x.len();
        assert!(m == y.len());
        let w = weights.unwrap_or(vec![1.0; m]);
        assert!(w.len() == m);
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
        let k = 3;
        let nest = m * k  + 1;
        let mut this = Self {
            t: knots.unwrap_or(vec![0.0; nest]),
            c: vec![0.0; nest],
            n: 0,
            fp: 0.0,
            k,
        };
        let xb = x[0];
        let xe = x[(m-1) as usize];
        let lwrk = m * (k + 1) + nest * (7 + 3 * k);
        let mut wrk = vec![0f64; lwrk as usize];
        let mut iwrk = vec![0i32; lwrk as usize];
        let mut ierr = 0;
        unsafe {
            curfit_(&iopt, 
                &m, 
                x.as_ptr(), 
                y.as_ptr(), 
                w.as_ptr(), 
                &xb, 
                &xe, 
                &this.k, 
                &s, 
                &nest, 
                &mut this.n, 
                this.t.as_mut_ptr(),
                this.c.as_mut_ptr(), 
                &mut this.fp, 
                wrk.as_mut_ptr(), 
                &lwrk, 
                iwrk.as_mut_ptr(), 
                &mut ierr
            );
        }
        this

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
        let m = x.len();
        let mut y = vec![0.0; m];
        let mut ierr = 0;
        unsafe {
            splev_(
                self.t.as_ptr(), 
                &self.n, 
                self.c.as_ptr(), 
                &self.k, 
                x.as_ptr(), 
                y.as_mut_ptr(), 
                &m, 
                &mut ierr)
        }
        y
    }
}



#[test]
fn test_smoothing(){
    let xi = vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
    let yi = xi.iter().map(|x|x*x).collect();
    let d = Spline::new_interpolating_spline(xi, yi, None);
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let y = d.values(&x);
    println!("{:.4?}", y);
}