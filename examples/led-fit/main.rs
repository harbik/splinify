#![doc = include_str!("README.md")]

use splinify::{CubicSplineFit2D, Result, read_csv_xy, SplineCurveData};

fn main()-> Result<()> {

    let example_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(std::path::Path::new(file!()).parent().unwrap());

    let csv_path = example_dir.join("leds4000.csv");
    let (x,y) = read_csv_xy(csv_path.to_str().expect("non-UTF-8 path"))?;
    let xy: Vec<f64> = x.iter().zip(y.iter()).flat_map(|(&x,&y)|[x,y]).collect();

    let d = CubicSplineFit2D::new(x.clone(), xy.clone())?
            .begin_constraints([ [x[0],y[0]], [0.0,0.0], [0.0,0.0]])?
            .end_constraints([ [*x.last().unwrap(), *y.last().unwrap()], [0.0,0.0], [0.0,0.0] ])?
            .smoothing_spline_optimize(0.1, 0.5,  |k,_,rms,_| {println!("{}/{:.4}", k, rms); k>50}, None)?;

    // output spline fit results
    println!("{}", serde_json::to_string_pretty(&SplineCurveData::from(&d))?);

    let plot_path = example_dir.join("fit.png");
    d.plot_with_control_points_and_data(plot_path.to_str().expect("non-UTF-8 path"), (2000,1000), &xy)?;

    Ok(())
}

