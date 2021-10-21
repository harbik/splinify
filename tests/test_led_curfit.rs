
use dierckx::{CubicCurveFit, CurveSplineFit, ParameterCurveSplineFit, Result, read_csv_xy, plot};

#[test]
fn test_smoothing() -> Result<()> {
    let (x,y) =  read_csv_xy("tests/leds4000.csv")? ;

    let d = CurveSplineFit::<3>::new(x.clone(), y.clone());
    let d = d.smoothing_spline(1.25)?;
    println!("knots {:?}", d.t);
    println!("number of knots: {}", d.t.len());
    println!("fp: {:?}", d.e_rms);

    plot("tests/img/curfit-smooth.png",x,y,d)?;



    Ok(())
}

#[test]
fn test_cardinal() -> Result<()> {

    let (x,y) =  read_csv_xy("tests/leds4000.csv")? ;

    let d = CubicCurveFit::new(x.clone(), y.clone());
    let tc = d.cardinal_spline(10.0)?;
    println!("knots {:?}", tc.t);
    println!("number of knots: {}", tc.t.len());
    println!("fp: {:?}", tc.e_rms);

    Ok(())
}


#[test]
fn test_interpolating_spline() -> Result<()> {
    use plotters::prelude::*;


    let x: Vec<f64> = (0..=180).map(|i|(i as f64).to_radians()).collect();
    let y: Vec<f64> = x.iter().map(|x|x.sin().powi(8)).collect();
    let x_data: Vec<f64> = (0..=180).step_by(20).map(|i|(i as f64).to_radians()).collect();
    let y_data: Vec<f64> = x_data.iter().map(|x|x.sin().powi(8)).collect();

    let root = BitMapBackend::new("interpolating_spline.png", (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(50)
        .caption("Interpolating", ("sans-serif",24))
        .build_cartesian_2d(0f64..180f64.to_radians(), 0.0..1.0)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(x_data.iter().cloned().zip(y_data.iter().cloned()).map(|xy| Circle::new(xy, 5, RED.filled())))?;
    // CurveFit test
    //  let d = CurveFit::<3>::new(x_data, y_data);
    // let y_fit = d.interpolating_spline()?.evaluate(&x)?;

    // ConstrainedSpline test
    let y_fit =
        ParameterCurveSplineFit::<3,1>::new(x_data, y_data.clone())?
            .begin_constraints([ [y_data[0]], [0.0], [0.0]])?
            .end_constraints([ [y_data[y_data.len()-1]], [0.0], [0.0] ])?
            .interpolating_spline()?
            .evaluate(&x)?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y_fit.iter().cloned()), &HSLColor(0.5, 1.0, 0.5)))?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y.iter().cloned()), &BLACK))?;

   
    root.present()?;


    Ok(())
}


#[test]
fn constrained_cardinal_spline() -> Result<()> {

    let (x,y) =  read_csv_xy("tests/leds4000.csv")? ;

    let cs = ParameterCurveSplineFit::<5,1>::new(x.clone(), y.clone())?;
    let tc = cs
       .begin_constraints([[y[0]],[0.0], [0.0]])?
       .end_constraints([[y[y.len()-1]], [0.0], [0.0]])?
       .cardinal_spline(5.0)?;
     
    println!("knots {:?}", tc.t);
    println!("number of knots: {}", tc.t.len());
    println!("fp: {:?}", tc.e_rms.unwrap());

    Ok(())
}
