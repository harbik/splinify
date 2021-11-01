
use splinify::{CubicCurveFit, CurveSplineFit, ParameterCurveSplineFit, Result, read_csv_xy, plot};

#[test]
fn test_smoothing() -> Result<()> {
    let (x,y) =  read_csv_xy("tests/data/leds4000.csv")? ;

    let d = CurveSplineFit::<3>::new(x.clone(), y.clone());
    let d = d.smoothing_spline(0.01)?;
    println!("knots {:?}", d.t);
    println!("number of knots: {}", d.t.len());
    println!("fp: {:?}", d.e);

    let json = serde_json::to_string_pretty(&d)?;
    println!("{}", json);

//    plot("tests/img/curfit-smooth.png",x,y,d)?;


    Ok(())
}

#[test]
fn test_cardinal() -> Result<()> {

    let (x,y) =  read_csv_xy("tests/data/leds4000.csv")? ;

    let d = CubicCurveFit::new(x.clone(), y.clone());
    let tc = d.cardinal_spline(10.0)?;
    println!("knots {:?}", tc.t);
    println!("number of knots: {}", tc.t.len());
    println!("fp: {:?}", tc.e);

    Ok(())
}


#[test]
fn test_interpolating_spline() -> Result<()> {

    let x: Vec<f64> = (0..=180).map(|i|(i as f64).to_radians()).collect();
    let xy: Vec<f64> = x.iter().flat_map(|&x|[x,x.sin().powi(8)]).collect();

    let x_data: Vec<f64> = (0..=180).step_by(30).map(|i|(i as f64).to_radians()).collect();
    let xy_data: Vec<f64> = x_data.iter().flat_map(|&x|[x,x.sin().powi(8)]).collect();

    // ConstrainedSpline test
    let int_spline =
        ParameterCurveSplineFit::<3,2>::new(x_data, xy_data.clone())?
            .begin_constraints([ [xy_data[0], xy_data[1]], [0.0, 0.0], [0.0, 0.0]])?
            .end_constraints([ [xy_data[xy_data.len()-2],xy_data[xy_data.len()-1]],  [0.0, 0.0], [0.0, 0.0] ])?
            .interpolating_spline()?;
    
    int_spline.plot_with_control_points_and_data("fit.png", (2000,1000), &xy)?;
    Ok(())
}


#[test]
fn constrained_cardinal_spline() -> Result<()> {

    let (x,y) =  read_csv_xy("tests/data/leds4000.csv")? ;

    let mut xy: Vec<f64> = Vec::new();
    x.iter().zip(y.iter()).for_each(|(x,y)| xy.extend([x,y]));

    let cs = ParameterCurveSplineFit::<5,1>::new(x.clone(), y.clone())?;
    let tc = cs
       .begin_constraints([[y[0]],[0.0], [0.0]])?
       .end_constraints([[y[y.len()-1]], [0.0], [0.0]])?
       .cardinal_spline(10.0)?;
     
    println!("knots {:?}", tc.t);
    println!("number of knots: {}", tc.t.len());

    tc.plot("fit.png", (2000,1000))?;

    Ok(())
}
