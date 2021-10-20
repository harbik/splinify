use super::{FitResult, Spline};
use plotters::prelude::*;

pub fn plot<const K: usize>(
    filepath: &str,
    x: Vec<f64>,
    y: Vec<f64>,
    s: Spline<K,1>,
) -> FitResult<()> {
    let y_max = y.iter().cloned().fold(y[0], f64::max);
    let y_min = y.iter().cloned().fold(y[0], f64::min);
    let y_height = y_max - y_min;

    let y_s = s.evaluate(&x)?;

    let spline_color = HSLColor(0.5, 0.5, 0.4);
    //let spline_coef_color = HSLColor(0.05, 0.5, 0.4);

    let root = BitMapBackend::new(filepath, (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root).margin(50).build_cartesian_2d(
        x[0]..x[x.len() - 1],
        y_min - y_height / 10.0..y_max + y_height / 10.0,
    )?;
    chart.draw_series(
        s.t.iter()
            .skip((K+1)/2)
            .cloned()
            .zip(s.c.iter().cloned())
            .map(|xy| Circle::new(xy, 8, spline_color.stroke_width(2))),
    )?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
        x.iter().cloned().zip(y_s.iter().cloned()),
        spline_color.mix(1.0).stroke_width(5),
    ))?;
    chart.draw_series(LineSeries::new(
        x.iter().cloned().zip(y.iter().cloned()),
        BLACK.mix(1.0).stroke_width(2),
    ))?;

    root.draw(
        &(EmptyElement::at((1550, 100))
            + Text::new(
                format!("{} knots  {:.4} rms", s.t.len(), s.e_rms.unwrap()),
                (0, 0),
                &"sans-serif".into_font().resize(40.0).color(&spline_color),
            )),
    )?;
    root.present()?;
    Ok(())
}

pub fn plot2D<const K: usize>(
    filepath: &str,
    x: Vec<f64>,
    y: Vec<f64>,
    s: Spline<K,2>,
) -> FitResult<()> {
    let y_max = y.iter().cloned().fold(y[0], f64::max);
    let y_min = y.iter().cloned().fold(y[0], f64::min);
    let y_height = y_max - y_min;

    let y_s = s.evaluate(&x)?;

    let spline_color = HSLColor(0.5, 0.5, 0.4);
    //let spline_coef_color = HSLColor(0.05, 0.5, 0.4);

    let root = BitMapBackend::new(filepath, (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root).margin(50).build_cartesian_2d(
        x[0]..x[x.len() - 1],
        y_min - y_height / 10.0..y_max + y_height / 10.0,
    )?;
    chart.draw_series(
        s.t.iter()
            .skip((K+1)/2)
            .cloned()
            .zip(s.c.iter().cloned())
            .map(|xy| Circle::new(xy, 8, spline_color.stroke_width(2))),
    )?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
        x.iter().cloned().zip(y_s.iter().cloned()),
        spline_color.mix(1.0).stroke_width(5),
    ))?;
    chart.draw_series(LineSeries::new(
        x.iter().cloned().zip(y.iter().cloned()),
        BLACK.mix(1.0).stroke_width(2),
    ))?;

    root.draw(
        &(EmptyElement::at((1550, 100))
            + Text::new(
                format!("{} knots  {:.4} rms", s.t.len(), s.e_rms.unwrap()),
                (0, 0),
                &"sans-serif".into_font().resize(40.0).color(&spline_color),
            )),
    )?;
    root.present()?;
    Ok(())
}