
# Splinify: Fit Curve and Surface Splines
<center>
    <br>
    <img src="https://www.harbik.com/img/dierckx/lissajous.png"/>
    <br>
    <i>A Splinify generated b-Spline curve, with 15 control points</i>
    <br>
    <i>See <a href="https://github.com/harbik/splinify/tree/main/examples/lissajous">lissajous</a> example for how this was generated</i> 
</center>
  
# Introduction
    
*Warning: this library uses a Fortran library as fitting engine, and requires ---in addition to Rust--- a [Fortran](#fortran) compiler*

Splinify fits curve and surface *B-Splines* to sets of data points,
    using Dierckx' Fortran FITPACK library written by Paul Dierckx in the mid 1980s.
Dierckx' library is still one of the the most advanced general B-spline fitting libraries available today ---
    its mathematical foundations and algorithms are described in detail in Paul Dierckx' book:
[Curve and Surface Fitting with Splines](https://www.google.com/books/edition/Curve_and_Surface_Fitting_with_Splines/-RIQ3SR0sZMC?hl=en "Paul Dierckx, Curve and Surface Fitting with Splines, Oxford University Press, 1993").

With all its features --- and being programmed in a somewhat limited language, at least as concerned to API design ---
    Dierckx' FITPACK is difficult to use:
    for example, its `concur` subroutine has 28 numerical arguments!
This aim of this library is to make it more accessible from Rust by implementing a high level *Rustacean* application interface.

Its intended use is scientific analysis and modeling, for example to represent spectral distributions in colorimetry,
or to interpolate positional data and use spline derivatives to get speed and acceleration of moving objects in one or more dimensions.

Although Dierckx' Fortran library covers fitting a spline evaluation methods for Univariate and Bivariate B-Splines,
    currently only its **univariate, or single input-parameter fitting methods** are implemented.


# Fortran

As this crate is a Fortran wrapper, you need a Fortran compiler on your system to use it:
in specific it requires the Gnu Fortran **GFortran** to be installed.
See [Installing GFortran](https://fortran-lang.org/learn/os_setup/install_gfortran) for instructions.

Having to install a Fortran compiler is a big restriction, so I recommend not using it as a dependency in projects for generic use.


# Instructions

Splinify is part of a family of three crates:

- **splinify** fits (non-uniform) [B-Spline](b-splines) curves to input data,
and results in a fitted as a `spliny`-crate `CurveSpline`.
Data inputs are `x` and `y` vectors for 1 dimensional curves,
and `u` and `xyn` vectors in case of N-dimensional curves.

- Use **spliny** to to use the generated splines, for example to calculate curve coordinates, or a spline-curve's derivatives.
This package also implements basic tools for input and output of spline representations in form of JSON files, and spline plots.
It is completely written in Rust, and does **not** require a Fortran compiler. 

- **dierckx-sys** contains Fortran foreign function interfaces to Paul Dierckx' FITPACK library. 
It is used by `splinify`, but ---unless you want to explore Paul Dierckx library yourself--- can be ignored as concerned to using `splinify` and `spliny`.


To use this library add this to your `Cargo.toml` file:

```
[dependencies]
splinify= "0.1"
```

# Examples






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

# Spline Fit Types

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


# License

All content in this repository is &copy;2021 Harbers Bik LLC, and licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>?)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.


