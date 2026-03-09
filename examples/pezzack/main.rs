
use splinify::{ParameterSplineCurveFit, Result, read_csv_xy, SplineCurveData};



fn main() -> Result<()> {
    let example_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(std::path::Path::new(file!()).parent().unwrap());

    let csv_path = example_dir.join("pezzack_noisy.csv");
    let (x,y) = read_csv_xy(csv_path.to_str().expect("non-UTF-8 path"))?;
    let xy: Vec<f64> = x.iter().zip(y.iter()).flat_map(|(&x,&y)|[x,y]).collect();

    let d =
    ParameterSplineCurveFit::<5,2>::new(x.clone(), xy.clone())?
       // .begin_constraints([ [y[0]] ])?
       // .end_constraints([ [y[y.len()-1]] ])?
        .smoothing_spline_optimize(0.1, 0.8,  |k,_,_,_| k>40, None)?;

    let json = serde_json::to_string_pretty(&SplineCurveData::from(&d))?;
    println!("{}", json);

    let plot_path = example_dir.join("fit.png");
    d.plot_with_control_points_and_data(plot_path.to_str().expect("non-UTF-8 path"), (1500,1000), &xy)?;

    Ok(())
}

