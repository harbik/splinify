#pragma once

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32"

/**
 * 
 * @param iopt -1: Least-squares spline fixed knots, 0,1: smoothing spline. iopt=0 and s=0: interpolating spline
 * @param m Number of data points supplied
 * @param x Array of x coordinates (at least m values)
 * @param y Array of y coordinates (at least m values)
 * @param w Array weights (at least m values)
 * @param xb Bounderies of the approximation interval. xb<=x(1), xe>=x(m)
 * @param xe 
 * @param k Degree of the spline
 * @param s Smoothing factor to be used if iopt >= 0
 * @param nest nest = m + k + 1
 * @param n Number of knots returned. For iopt=-1 value needs to pe specified on entry
 * @param t Array of dimension of at least nest. For iopt=-1 array of knots to be used for lsq spline
 * @param c Double array of at least nest. Will contain the coefficients of the b-spline representation
 * @param fp Weighted sum of the squared residuals of the spline approximation.
 * @param wrk Double array of dimension at least (m(k+1)+nest(7+3k)).
 * @param lwrk Size of 'wrk'
 * @param iwrk int Array of at least nest (m + k + 1)
 * @param ier Error flag.
 */
int curfit_(int *iopt, int *m, double *x, double *y, double *w, double *xb, double *xe, int *k, 
	double *s, int *nest, int *n, double *t, double *c, double *fp, double *wrk, int *lwrk, int *iwrk, int *ier);

/**
 * 
 * @param t Array,length n, which contains the position of the knots.
 * @param n Integer, giving the total number of knots of s(x).
 * @param c Array,length n, which contains the b-spline coefficients.
 * @param k Integer, giving the degree of s(x).
 * @param x Array,length m, which contains the points where s(x) must be evaluated. 
 * @param y Array,length m, giving the value of s(x) at the different points.
 * @param m Integer, giving the number of points where s(x) must be evaluated.
 * @param ier Error value: ier = 0 normal return, ier = 10, invalid input data (see restrictions).
 */
int splev_(double *t, int *n, double *c, int *k, double *x, double *y, int *m, int *ier); 