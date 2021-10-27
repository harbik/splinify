
# Fit smoothing spline for LED spectrum with optimizations

This example demonstrates the use the `smoothing_spline_optimize` in `ParameterCurveSplinefit`.
The input data is an LED spectrum, read from a csv-file, named `leds4000.csv`.
This file is an export from a spreadsheet, with wavelength values in the first column, and spectral values in a second column.
Wavelengths range from 380 to 780 nanometer, and are irregularly spaced, with an average step size of about 0.4nm.
It contains 941 data-points.

Here we are going to fit a constrained spline, with first and second derivatives set to 0.0 at both ends.

<img src="https://www.harbik.com/img/dierckx/led26.png"/>

The original data is shown as a black line.
As you can see, the data is a bit noisy, which is not uncommon with light measurements.

If you can predict the noise level in the data, from physics principles, a smoothing spline fit can be obtained directly using `smoothing_spline` --- 
if you don't, you can start with a very high estimate, and try successive fits with improved fit precision.
This procedure is implemented with `smoothing_spline_optimize`, following Dierckx' instructions.

When running this program, this function will start with the `rms_start` error start value, and reduce it with successive approximations by multiplying this value with a scale ratio: `rms_scale_ratio`.
As third parameter a criterium for convergence has to be supplied, in form of a function, typically a Rust closure.
Here it is set to stop when the number of knots `k` is larger then 30 as an upper limit, based on an estimate of the spectral resolution in the graph.
This convergence function is also used to get a view of the number of knots in the fit, and the accuracy of the fit.
The initial target for `rms-start` is here set to `1.0`, and in each iteration it is reduced by a factor of 0.5.

```text
9/0.2434
12/0.1217
17/0.0608
18/0.0304
22/0.0152
26/0.0076
315/0.0038
```
As you can see, there is a big jump in the number of knots from 26 to 315 in the last iteration step: this is where the algorithm starts to fit the noise.
This method will give you the solution which is within the set convergence limits, so in this case that is the solution with the 26 knots, and refit with a smoothing spline with its `rms` target.

The fit results are here printed on the command line, and  plotted in an "fit.png" image file as shown above.

```json
{
  "t": [
    381.09, 381.09, 381.09, 381.09, 405.04, 429.15, 441.47, 447.64, 450.53, 
    453.41, 465.4, 477.42, 489.89, 501.99, 526.71, 551.58, 576.17, 601.33, 
    626.63, 652.07, 677.21, 728.3, 779.88, 779.88, 779.88, 779.88
  ],
  "c": [
    0.0022, 0.0022000000000000006, 0.011333202500664304, -0.007280989002300535, 
    0.3701504089598165, 1.0271690365100785, 0.9338346755132363, 0.5240967440473344, 
    0.2996392495021921, 0.2209199301647846, 0.5503228933662222, 0.643534450394433, 
    0.5184346552047464, 0.6488285856664529, 0.6323493134354996, 0.660478648460517,
    0.8052787571658832, 0.7380670192932405, 0.38635154293989304, 0.08489388281768777, 
    0.035, 0.035
  ],
  "k": 3,
  "n": 1,
}
```
This spline representation has 22 control points, using 26 knots.
The memory size of this representation --using type `f64`-- is 390 bytes, while the original data set is 15056 bytes, a reduction of almost 40x.
This is not that important for local (desktop) use, but for using large spectral datasets in a browser environment this is a substantial benefit.
