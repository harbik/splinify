
extern crate dierckx;
use dierckx::{ParameterCurveSplineFit, Result, read_csv_xy, plot};

/*

fn plot<const K:usize> (filepath: &str, x: Vec<f64>, y: Vec<f64>, s: Spline<K>) -> FitResult<()>{
    use plotters::prelude::*;

    let y_max = y.iter().cloned().fold(y[0], f64::max);
    let y_min = y.iter().cloned().fold(y[0], f64::min);
    let y_height = y_max-y_min;

    let root = BitMapBackend::new(filepath, (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(50)
        .build_cartesian_2d(x[0]..x[x.len()-1], y_min-y_height/10.0..y_max+y_height/10.0)?;
    chart.configure_mesh().draw()?;
   // chart.draw_series(LineSeries::new(x.iter().cloned().zip(y.iter().cloned()), BLACK.mix(0.5).stroke_width(8)))?;
    chart.draw_series(x.iter().cloned().zip(y.iter().cloned()).map(|xy| Circle::new(xy, 9, HSLColor(0.5,0.5,0.45).filled())))?;

    let y_s = s.evaluate(&x)?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y_s.iter().cloned()), BLACK.stroke_width(2)))?;
    chart.draw_series(s.t.iter().cloned().zip(s.c.iter().cloned()).map(|xy| Circle::new(xy, 10, HSLColor(0.1,0.0,0.25).stroke_width(3))))?;

    root.present()?;
    Ok(())

}
*/


#[test]
fn test_pezzack() -> Result<()> {

    let (x,y) =  read_csv_xy("tests/data/pezzack_noisy.csv")? ;

    let d = 
    ParameterCurveSplineFit::<3,1>::new(x.clone(), y.clone())?
        .begin_constraints([ [y[0]] ])?
        .end_constraints([ [y[y.len()-1]] ])?
        .smoothing_spline_optimize(0.1, 0.8,  |k,_,_,_| k>40, None)?;

    //let y_s = d.evaluate(&x)?;
    plot("tests/img/pezzack_noisy.png", x, y, d)?;

    Ok(())
}

