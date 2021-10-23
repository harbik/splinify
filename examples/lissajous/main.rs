
use dierckx::{ParameterCurveSplineFit, Result, read_csv_uxy, plot2d};


fn lissajous(t:f64, a: f64, kx: f64, b: f64, ky: f64) -> [f64;2] {
    [
        a * (kx * t).cos(),
        b * (ky * t).sin()
    ]

}

fn main() -> Result<()> {

    std::env::set_current_dir(std::path::Path::new(file!()).parent().unwrap())?;   

    let ts: Vec<f64> = (0..=180u32).into_iter().step_by(1).map(|v|(v as f64).to_radians()).collect();
    let mut xy: Vec<f64> = Vec::with_capacity(ts.len()*2);

    for t in &ts {
        xy.extend(lissajous(*t, 1.0, 3.0, 1.0, 2.0));
    }


    let d = 
        ParameterCurveSplineFit::<3,2>
            ::new(ts.clone(), xy.clone())?
            .smoothing_spline(1E-2)?;

    let json = serde_json::to_string_pretty(&d)?;
    println!("{}", json);

    let mut x: Vec<f64> = Vec::new();
    let mut y: Vec<f64> = Vec::new();
    xy.chunks(2).for_each(|xy|{x.push(xy[0]);y.push(xy[1]);});

   plot2d("fit.png",d, ts, x, y)?;


    Ok(())
}

