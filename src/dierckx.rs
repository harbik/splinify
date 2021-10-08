
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

    pub(super) fn curfit_(
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

    pub(super) fn splev_(
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

