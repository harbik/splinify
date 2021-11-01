
use splinify::{ParameterCurveSplineFit, Result, read_csv_xy, plot};



fn main() -> Result<()> {
    std::env::set_current_dir(std::path::Path::new(file!()).parent().unwrap())?;   

    let (x,y) =  read_csv_xy("pezzack_noisy.csv")? ;
    let xy = x.iter().zip(y.iter()).flat_map(|(x,y)|[*x,*y]).collect();

    let d = 
    ParameterCurveSplineFit::<5,1>::new(x.clone(), y.clone())?
        .begin_constraints([ [y[0]] ])?
        .end_constraints([ [y[y.len()-1]] ])?
        .smoothing_spline_optimize(0.1, 0.8,  |k,_,_,_| k>40, None)?;

    let json = serde_json::to_string_pretty(&d)?;
    println!("{}", json);

    d.plot_with_control_points_and_data("fit.png", (1500,1000), xy)?;

    Ok(())
}

