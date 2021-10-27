
#![doc = include_str!("README.md")]

use splinify::{CubicSplineFit2D, plot2d, Plot2DFlags, Result};

fn lissajous(t:f64, a: f64, kx: f64, b: f64, ky: f64) -> [f64;2] {
    [
        a * (kx * t).cos(),
        b * (ky * t).sin()
    ]
}

fn main() -> Result<()> {

    let mut xy = Vec::new();
    let mut xx = Vec::new();
    let mut yy = Vec::new();

    // set path to directory of example script, used for saving plot
    std::env::set_current_dir(std::path::Path::new(file!()).parent().unwrap())?;   

    // angle parameter, from 0 to 180º, with 1º-steps
    let us: Vec<f64> = (0..=180u32).into_iter().map(|v|(v as f64).to_radians()).collect();

    for t in &us {
        let [x,y] = lissajous(*t, 1.0, 3.0, 1.0, 5.0);
        xy.extend(&[x,y]);
        xx.push(x);
        yy.push(y);
    }

    // fit Cubic 
    let d = CubicSplineFit2D::new(us.clone(), xy.clone())?
            .smoothing_spline(5E-3)?;

    // output spline fit results 
    println!("{}", serde_json::to_string_pretty(&d)?);

    // graph of spline with knots, produced with help of the `plotters` library
    plot2d("lissajous.png",d, us, xx, yy, Some(Plot2DFlags::SPLINE | Plot2DFlags::KNOTS_FILLED))?;

    Ok(())
}

