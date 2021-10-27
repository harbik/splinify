
# Parametric 2-dimension Lissajous curve spline fit

![2D Lissajous SplineCurve](https://harbik.github.io/img/dierckx/lissajous.png)

Here we fit a two-dimensional cubic spline using `CubicSplineFit2D` `smoothing_spline` fit method.
The input data is here generated from the `lissajous` function, which produces `x` and `y`-values from an input angle as parameter, and two harmonic functions.

```rust
fn lissajous(t:f64, a: f64, kx: f64, b: f64, ky: f64) -> [f64;2] {
    [
        a * (kx * t).cos(),
        b * (ky * t).sin()
    ]
}
```

Input angle `t` is varied from 0 to 180º, with steps of 1º.
