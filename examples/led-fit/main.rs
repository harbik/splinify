
use dierckx::{ParameterCurveSplineFit, Result, read_csv_xy, plot};



fn main()-> Result<()> {
    std::env::set_current_dir(std::path::Path::new(file!()).parent().unwrap())?;   
    let (x,y) =  read_csv_xy("leds4000.csv")? ;

    let d = 
    ParameterCurveSplineFit::<3,1>::new(x.clone(), y.clone())?
        .begin_constraints([ [y[0]], [0.0], [0.0]])?
        .end_constraints([ [y[y.len()-1]], [0.0], [0.0] ])?
        .smoothing_spline_optimize(1.0, 0.5,  |k,_,rms,_| {println!("{}/{:.4}", k, rms); k>30}, None)?;
    println!("knots {:?}", d.t);
    println!("number of knots: {}", d.t.len());

    plot("fit.png", x, y, d)?;

    Ok(())
}
