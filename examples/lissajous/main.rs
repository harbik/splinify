
#![doc = include_str!("README.md")]

use splinify::{CubicSplineFit2D, Result, SplineCurveData};

fn lissajous(t:f64, a: f64, kx: f64, b: f64, ky: f64) -> [f64;2] {
    [
        a * (kx * t).cos(),
        b * (ky * t).sin()
    ]
}

fn main() -> Result<()> {


    let example_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(std::path::Path::new(file!()).parent().unwrap());

    // Generate Lissajous data points, with angle parameter ranging from 0 to 180º, with 1º-steps.
    let u: Vec<f64> = (0..=180u32).into_iter().map(|v|(v as f64).to_radians()).collect();
    let xy: Vec<f64> = u.iter().flat_map(|t|lissajous(*t,1.0, 3.0, 1.0, 5.0)).collect();

    // fit Cubic Spline with Splinify's CubicSplineFit
    let s = CubicSplineFit2D::new(u, xy)?.smoothing_spline(5e-3)?;

    // Output fit results as JSON file and plot
    println!("{}", serde_json::to_string_pretty(&SplineCurveData::from(&s))?);
    let plot_path = example_dir.join("lissajous.png");
    s.plot_with_control_points(plot_path.to_str().expect("non-UTF-8 path"), (800,800))?;

    Ok(())
}

