
# Fit interpolating spline for CIE daylight spectral data

This example demonstrates fitting a 2D parametric interpolating spline to the CIE standard
daylight spectral power distribution component S0, using `ParameterSplineCurveFit`.

The S0 spectral data covers wavelengths from 300 to 830 nm in 5 nm steps (107 data points).
Rather than reading from a CSV file, the standard tabulated values are embedded directly in
the source code, as they are part of the CIE standard and do not change.

The data is formatted as a 2D parametric curve: each sample is an *(x, y)* pair consisting
of the wavelength and the corresponding S0 spectral value, with wavelength used as the curve
parameter *u*.

Boundary constraints are set so that the first and second derivatives at both endpoints
are zero, producing a smooth, physically realistic fit at the spectral limits.

```text
cargo run --example cie-daylight --features plot
```

The fit results are printed to the command line as JSON, and plotted in a `fit.png` file.

```json
{
  "k": 3,
  "n": 2,
  "t": [...],
  "c": [...]
}
```

This spline representation can be used with the `spliny` crate to evaluate the interpolated
spectrum at any wavelength within the fitted range.
