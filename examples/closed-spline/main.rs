use splinify::{ClosedCubicSplineFit2D, Result};

fn ellipse(t: f64, a: f64, b: f64) -> [f64; 2] {
    [a * t.cos(), b * t.sin()]
}

fn main() -> Result<()> {
    std::env::set_current_dir(std::path::Path::new(file!()).parent().unwrap())?;

    let n_points = 1000;
    let mut u: Vec<f64> = (0..n_points).map(|i| i as f64 / n_points as f64).collect();
    let mut xy: Vec<f64> = u
        .iter()
        .flat_map(|&t| ellipse(t * std::f64::consts::TAU, 2.0, 1.0))
        .collect();

    // Close the curve: append first point with u = 1.0
    xy.extend([xy[0], xy[1]]);
    u.push(1.0);

    let d = ClosedCubicSplineFit2D::new(u, xy)?.smoothing_spline(1E-6)?;

    let json = serde_json::to_string_pretty(&d)?;
    println!("{}", json);

    d.plot("fit.png", (2000, 1000))?;

    Ok(())
}
