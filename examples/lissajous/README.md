
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

A dataset is generated, with parameter `us`, by varying input angle `t` from 0 to 180º, with steps of 1º,
and an `xy` vector, containing corresponding  [x(t),y(t)] coordinates:
``` rust
    let us: Vec<f64> = (0..=180u32).into_iter().map(|v|(v as f64).to_radians()).collect();

    for t in &us {
        let [x,y] = lissajous(*t, 1.0, 3.0, 1.0, 5.0);
        xy.extend(&[x,y]);
        xx.push(x);
        yy.push(y);
    }
```
Also vectors `xx` and `yy` are generated here, which are used for plotting purposes, but are not needed otherwise.

The spline representation is calculated using:
```rust
    let s = CubicSplineFit2D::new(us.clone(), xy.clone())?
            .smoothing_spline(5e-3)?;
```
The inputs here are cloned, also to be able to reuse them for plotting; if you don't need a plot, clones are not needed.
The target accuracy ---as an rms error value--- of the spline fit is set to `5e-3`, or `0.005`;
with a smaller error the spline fit will be better, but will contain more control points.

With this precision, we the following spline representation, printed to the command line output as a JSON object:
```rust
    println!("{}", serde_json::to_string_pretty(&s)?);
```
Here is the JSON output:
```json
{
  "t": [
    0.0, 0.0, 0.0, 0.0, 0.4014257279586958, 0.7853981633974483, 0.9948376736367679,
    1.1868238913561442, 1.3788101090755203, 1.5707963267948966, 1.7802358370342162, 
    1.9722220547535925, 2.1642082724729685, 2.356194490192345, 2.7576202181510405, 
    3.141592653589793, 3.141592653589793, 3.141592653589793, 3.141592653589793
  ],
  "c": [
    0.9961805460172887, 1.01581609212485, 0.45785551737300106, -0.6743400479743561, 
    -1.04606086926188, -0.9655723250526262, -0.5757280827332156, 0.017487591989126007, 
    0.610401362313724, 0.9866605336566671, 1.0335232431543322, 0.6453772441164307, 
    -0.4846359310719992, -1.014195365951515, -0.9964502165043656, -0.027374797128379907, 
    0.770580186441133, 1.5844468932141083, -0.7504988105159386, -1.1533053154158592,
    -0.39940623356854926, 0.669953241001308, 1.1822672685309246, 0.6182210556297563, 
    -0.493261782484514, -1.1530136311265229, -0.6885569654105621, 1.6217843337199846, 
    0.7207977628265237, -0.02317389641793977,
  ],
  "k": 3,
  "n": 2,
}
```
The `t` array contains the knot values, in this case in form of parameter `t` angles, in which the curve is split up for this fit: it's size is 19, including three repeats at the start, and three repeats at the end, which are needed for spline evaluation.
The `c` array contains the coordinates of 15 control points with their `x` values in the first 15 values, and their `y` values in the second set of 15 values.
The values `k` and `n` are the spline degree, for a cubic spline 3, and the dimension of the spline, in this case 2.
