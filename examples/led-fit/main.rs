#![doc = include_str!("README.md")]

use splinify::{CubicSplineFit2D, Result, read_csv_xy};

fn main()-> Result<()> {

    std::env::set_current_dir(std::path::Path::new(file!()).parent().unwrap())?;   

    let (x,y) =  read_csv_xy("leds4000.csv")? ;
    let xy: Vec<f64> = x.iter().zip(y.iter()).flat_map(|(&x,&y)|[x,y]).collect();


    let d = CubicSplineFit2D::new(x.clone(), xy.clone())?
            .begin_constraints([ [x[0],y[0]], [0.0,0.0], [0.0,0.0]])?
            .end_constraints([ [*x.last().unwrap(), *y.last().unwrap()], [0.0,0.0], [0.0,0.0] ])?
            .smoothing_spline_optimize(0.1, 0.5,  |k,_,rms,_| {println!("{}/{:.4}", k, rms); k>50}, None)?;

    // output spline fit results 
    println!("{}", serde_json::to_string_pretty(&d)?);

    d.plot_with_control_points_and_data("fit.png", (2000,1000), &xy)?;

    Ok(())
}

