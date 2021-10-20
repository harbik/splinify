
use dierckx::{ParameterCurveSplineFit, FitResult, read_csv_xy, plot};



#[test]
fn test_smoothing() -> FitResult<()> {

    let (x,y) =  read_csv_xy("tests/data/leds4000.csv")? ;

    let d = 
    ParameterCurveSplineFit::<3,1>::new(x.clone(), y.clone())?
        .begin_constraints([ [y[0]], [0.0], [0.0]])?
        .end_constraints([ [y[y.len()-1]], [0.0], [0.0] ])?
        .smoothing_spline_optimize(1.0, 0.5,  |k,_,rms,_| {println!("{}/{:.4}", k, rms); k>30}, None)?;
    println!("knots {:?}", d.t);
    println!("number of knots: {}", d.t.len());

    plot("tests/img/fit.png", x, y, d)?;

    Ok(())
}


#[test]
fn test_interpolating_spline() -> FitResult<()> {

    let x: Vec<f64> = (0..=180).map(|i|(i as f64).to_radians()).collect();
    let y: Vec<f64> = x.iter().map(|x|x.sin().powi(8)).collect();
    let x_data: Vec<f64> = (0..=180).step_by(20).map(|i|(i as f64).to_radians()).collect();
    let y_data: Vec<f64> = x_data.iter().map(|x|x.sin().powi(8)).collect();

    let s =
        ParameterCurveSplineFit::<5,1>::new(x_data, y_data.clone())?
            .begin_constraints([ [y_data[0]], [0.0], [0.0]])?
            .end_constraints([ [y_data[y_data.len()-1]], [0.0], [0.0] ])?
            .interpolating_spline()?;

    plot("tests/img/fit.png", x, y, s) ?;
    Ok(())
}


#[test]
fn constrained_cardinal_spline() -> FitResult<()> {
    let (x,y) =  read_csv_xy("tests/data/leds4000.csv")? ;

    let cs = ParameterCurveSplineFit::<1,1>::new(x.clone(), y.clone())?;
    let tc = cs
       .begin_constraints([[y[0]],[0.0]])?
       .end_constraints([[y[y.len()-1]], [0.0]])?
       .cardinal_spline(5.0)?;
     
    println!("knots {:?}", tc.t);
    println!("number of knots: {}", tc.t.len());
    println!("fp: {:?}", tc.e_rms.unwrap());
    plot("tests/img/fit.png", x, y, tc) ?;


    Ok(())
}

