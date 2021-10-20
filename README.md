# Dierckx' FITPACK B-Spline Fortran Library Wrapper

Dierckx is a Fortran library for fitting curve and surface *B-Splines* to sets of data points.  It was written by Paul
Dierckx in the mid 1980s, and it is still the most advanced general B-spline fitting library today.
It is based on a solid mathematical foundation, implements many advanced algorithms, as described in detail in Paul Dierckx' book:
[Curve and Surface Fitting with Splines](https://www.google.com/books/edition/Curve_and_Surface_Fitting_with_Splines/-RIQ3SR0sZMC?hl=en "Paul Dierckx, Curve and Surface Fitting with Splines, Oxford University Press, 1993").

With all its features, and being programmed in a somewhat limited language (at least as concerned to API design), 
Dierckx' FITPACK is difficult to use:
for example, its `concur` subroutine has 28 numerical arguments!

This aim of this library is to make the Dierckx library more accessible from Rust,
by implementing a foreign function interface,
and by defining and implementing high level *Rustacean* application interface.

Although Dierckx' Fortran library covers fitting a spline evaluation methods for Univariate and Bivariate B-Splines,
currently this library only covers **univariate, or single input-parameter fitting methods**.


# Usage

This functionality of this library is split in two crates:

- the first crate is `dierckx-sys` containing the foreign function Rust interfaces, and Dierckx' Fortran source files,
- and the second is this crate, `dierckx`, which implements the higher level Rust object models, and associated methods.

As this crate is a Fortran wrapper, you need a Fortran compiler on your system to use it:
in specific, at the moment, it requires the Gnu Fortran **GFortran** to be installed.
See [Installing GFortran](https://fortran-lang.org/learn/os_setup/install_gfortran) how to do this on your machine.

Having to install a Fortran compiler is a big restriction, so I recommend not using it as a dependency in big projects.

To use this library add this to your `Cargo.toml` file:

```
[dependencies]
dierckx = "0.0.1"
```

# Examples


```

```

<div style="text-align:center;">
<img src="https://www.harbik.com/img/dierckx/led26.png" style="width:80%;"/>
</div>


# B-Splines

Some use for B-Spline representations are:
- **interpolation** in one or more dimensions, 
    for example when data is needed between a set of measured data points to make line plots, 
- calculate **derivatives** for measured data,
such as the momentary speed and acceleration of a car from a set of time and distance values,
- **smoothing**, to calculate a smooth line or shape through noisy data,
for example to smooth the hand written script on a tablet.

B-splines --or basis splines-- are polynomials, not fitted to to the whole set of data, but fitted to segments of the
data set, and connecting at the ends.
The collection of all the polynomials for all the segments is a curve (or surface) B-Spline representation.
Besides that the polynomials are all connected, to form a continuous line, they are also constructed to have the same
first, second, or higher (depending on the degrees polynomial) derivatives, to form smooth curves and surfaces.

Here is a quote from [Wikipedia](https://en.wikipedia.org/wiki/B-spline) with regard to B-Splines:
> A spline function of order n is a piecewise polynomial function of degree n−1 in a variable x. 
> The places where the pieces meet are known as knots. 

B-Splines can be used to:
- fit curves, in which case they are called *Univariate B-Splines*, to fit one or more outputs to a single input,
- or can be used to fit surfaces, called *Bivariate B-Splines*, used to fit one or more outputs to two inputs,
- or can be used to fit "volumes" using *Multivariate B-Splines*, by fitting one or more outputs to three or more inputs.

For example, the speed of a car as function of time could be described by a univariate spline, and weather-forecast
barometric air-pressure distribution can be described by a bivariate spline.
An example of a multivariate B-Spline representation would be a spline fit of 3D-temperature measurements in a room.




# Curve Fitting 

This wrapper provides various interface objects to Dierckx' library:

- [`CurveSplineFit<K>`][crate::curfit], which wraps `curfit`, used to fit a curve, with degree *K, to a set of
*(x<sub>i</sub>,y<sub>i</sub>)* data points, with the condition *x<sub>i+1</sub>&gt;x<sub>i+1</sub>*; an example of this
type of data set would be an array with time and temperatures, measured over a period of time, and,

- [`ParameterCurveSplineFit<K,N>`][crate::concur], wrapping `concur`, to fit an --optionally constrained-- parametrized
curve, with degree *K*, to data points consisting of a input-parameter value *u<sub>i</sub>*, also with the condition that
*u<sub>i+1</sub>&gt;u<sub>i+1</sub>*, and a set of *(x<sub>i</sub>,y<sub>i</sub>)* in two dimensions *(N=2)*, or a set of
*(x<sub>i</sub>,y<sub>i</sub>,z<sub>i</sub>)* in three dimensions *(N=3)*, or even more dimensions *(N<=10)*; an example of this
is a trajectory of a fly flying over a dinner table.

Both of these are generic spline implementations,
 using spline degree *`K`* and space dimension *`N`* as type parameters.
Dierckx strongly recommends to use only odd degree linear *(N=1)*, cubic *(N=3)*, and quintic *(N=5)* spline functions.
To avoid using type parameters, the following type aliases have been defined:

- **Type aliases for `CurveSplineFit<K>`**

|      Alias            | K |
|-----------------------|:-:|
|`LinearSplineFit`      | 1 |
|`CubicSplineFit`       | 3 |
|`QuinticSplineFit`     | 5 |



- **Type aliases for `ParameterCurveSplineFit<K,N>`**

|          Alias        | K | N |
|-----------------------|:-:|:-:|
| `LinearSplineFit2D`   | 1 | 2 |
| `CubicSplineFit2D`    | 3 | 2 |
| `QuinticSplineFit2D`  | 5 | 2 |
| `LinearSplineFit3D`   | 1 | 3 |
| `CubicSplineFit3D`    | 3 | 3 |
| `QuinticSplineFit3D`  | 5 | 3 |

To use a constrained spline-fit for one-dimensional curves,
 you can use a parametrized curve representation (use x as 'u', and y as 'xn') to use Dierckx' `concur` instead of
 the default `curfit` subroutine:

|      Alias            | K | N |
|-----------------------|:-:|:-:|
|`LinearSplineFit1D`    | 1 | 1 |
|`CubicSplineFit1D`     | 3 | 1 |
|`QuinticSplineFit1D`   | 5 | 1 |




# Spline-Fit Types

For a dataset, the following spline-types can be generated:

- **Interpolating Splines**   
These splines fit the given data exactly, with no error, with the knots at the location of the input parameter values.
Besides the input data, no further input is needed, unless you want to set specific boundary conditions.

- **Least-Squares Splines**  
These are obtained when you specify the location of the knots. 
The splines will be fitted to minimize the square deviation between the spline and data points.
The resulting spline values can deviate from the given output data values.
It is common to use equidistant locations of the knots, in which case the spline is referred to as a *cardinal spline**.
Additional input data here are the location of the knots, and optional boundary conditions.

- **Smoothing spline**   
Smoothing splines are generated by supplying a measure for the estimated 'noise' in the data, in form of an RMS error (Root Mean Square noise value).
The smoothing spline algorithm will add knots until the fit error is less than the given RMS noise estimate.
A weighting factor can be added to each datapoint too, to limit the effect of those points on the fit result.



# Spline` Spline Represenations
The result of the spline fits is an array of knot locations, 
 consist of curves as single input parameter values,
 and *N* ---for each dimension one--- coefficient values.
 
In this library these are represented by a `Spline<K,N>` object,
 with *K* the spline degree, 
 *N* the space dimension,
 and containing two vectors,
 a vector *t*, containing the knot locations,
 and a vector *c*, containing the spline coefficients.

For convenience, the following aliases have been defined:

|          Alias        | K | N |
|-----------------------|:-:|:-:|
| `LinearSpline`        | 1 | 1 |
| `CubicSpline`         | 3 | 1 |
| `QuinticSpline`       | 5 | 1 |
| `LinearSpline2D`      | 1 | 2 |
| `CubicSpline2D`       | 3 | 2 |
| `QuinticSpline2D`     | 5 | 2 |
| `LinearSpline3D`      | 1 | 3 |
| `CubicSpline3D`       | 3 | 3 |
| `QuinticSpline3D`     | 5 | 3 |


