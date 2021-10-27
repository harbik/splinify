use super::{Result, Spline};
use csv::ReaderBuilder;
use plotters::prelude::*;
use bitflags::bitflags;

pub fn plot<const K: usize>(
    filepath: &str,
    x: Vec<f64>,
    y: Vec<f64>,
    s: Spline<K,1>,
) -> Result<()> {
    let y_max = y.iter().cloned().fold(y[0], f64::max);
    let y_min = y.iter().cloned().fold(y[0], f64::min);
    let y_height = y_max - y_min;

   // let y_s = s.evaluate(&x)?;
    let y_s = s.eval_n(&x)?;

    let spline_color = HSLColor(0.5, 0.5, 0.4);
    //let spline_coef_color = HSLColor(0.05, 0.5, 0.4);

    let mut root = BitMapBackend::new(filepath, (2000, 1000));
    root.draw_rect((0,0), (2000,1000), &WHITE, true)?;

    let chartarea= root
        .into_drawing_area()
        .margin(100, 100, 100, 100);
    chartarea.fill(&HSLColor(0.1, 0.5, 0.95))?;
    //root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&chartarea)
        .margin(50)
        .build_cartesian_2d(
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


    let leg_txt = if let Some(e) = s.e {
        format!("{} knots  {:.1e} rms", s.t.len(), e)
    } else {
        format!("{} knots", s.t.len())
    };

    chartarea.draw(
        &(EmptyElement::at((1450, 100))
            + Text::new(
                leg_txt,
                (0, 0),
                &"sans-serif".into_font().resize(40.0).color(&spline_color),
            )),
    )?;
    chartarea.present()?;
    Ok(())
}


bitflags! {
    pub struct Plot2DFlags: u32 {
        const SPLINE =          0b00000001;
        const LEGEND =          0b00000010;
        const DATA =            0b00000100;
        const PLOT_KNOTS =      0b00001000;
        const KNOTS_FILLED =    0b00010000;
        const ALL=              0b00011111;
    }
}



/// Plots a two-dimensional (xy) spline curve, its knots, and 
///  
pub fn plot2d<const K: usize, const N: usize>(
    filepath: &str,
    s: Spline<K,N>,
    u: Vec<f64>,
    x: Vec<f64>,
    y: Vec<f64>,
    flags: Option<Plot2DFlags>,
) -> Result<()> {
    let x_max = x.iter().cloned().reduce(f64::max).unwrap();
    let x_min = x.iter().cloned().reduce(f64::min).unwrap();
    let y_max = y.iter().cloned().reduce(f64::max).unwrap();
    let y_min = y.iter().cloned().reduce(f64::min).unwrap();
    let width = x_max - x_min;
    let height = y_max - y_min;
    let flags = flags.unwrap_or(Plot2DFlags::ALL);


    let spline_color = HSLColor(0.5, 0.5, 0.4);
    //let spline_coef_color = HSLColor(0.05, 0.5, 0.4);

    let mut chartarea = BitMapBackend::new(filepath, (2000, 1000)).into_drawing_area();
    chartarea.fill(&WHITE)?;
    chartarea = chartarea.margin(100, 100, 100, 100);
    chartarea.fill(&HSLColor(0.1, 0.5, 0.95))?;

    let mut chart = ChartBuilder::on(&chartarea).margin(50).build_cartesian_2d(
        x_min - width / 10.0..x_max + width / 10.0,
        y_min - height / 10.0..y_max + height / 10.0,
    )?;


    // draw the mesh
    chart.configure_mesh().draw()?;

    if flags.contains(Plot2DFlags::PLOT_KNOTS) || flags.contains(Plot2DFlags::KNOTS_FILLED) {
        // plot the knots, as represented in spline c
        let c_x = &s.c[0..s.c.len()/2-(K+1)];
        let c_y = &s.c[s.c.len()/2..s.c.len()-(K+1)];
        if flags.contains(Plot2DFlags::KNOTS_FILLED) {
            chart.draw_series(
                c_x.iter().zip(c_y.iter())
                .map(|(&x, &y)|Circle::new((x,y), 10, spline_color.filled())),
            )?;
        } else {
            chart.draw_series(
                c_x.iter().zip(c_y.iter())
                .map(|(&x, &y)|Circle::new((x,y), 10, spline_color.stroke_width(2))),
            )?;
        }
    }

    if flags.contains(Plot2DFlags::SPLINE) {
        // spline 
        let s_xy = s.eval_n(&u)?;
        chart.draw_series(LineSeries::new(
            s_xy.chunks(2).map(|xy|(xy[0],xy[1])),
            spline_color.mix(1.0).stroke_width(5)))?;
    }


    if flags.contains(Plot2DFlags::DATA) {
        // plot the original input data used for spline fitting
        chart.draw_series(LineSeries::new(
            x.iter().cloned().zip(y.iter().cloned()),
            BLACK.mix(1.0).stroke_width(2),
        ))?;
    };

    if flags.contains(Plot2DFlags::LEGEND) {
        let leg_txt = if let Some(e) = s.e {
            format!("{} knots  {:.1e} rms", s.t.len(), e)
        } else {
            format!("{} knots", s.t.len())
        };

        chartarea.draw(
            &(EmptyElement::at((1550, 100))
                + Text::new( leg_txt, (0, 0), &"sans-serif".into_font().resize(40.0).color(&spline_color),
                )),
        )?;
    }

    chartarea.present()?;
    Ok(())
}

pub fn read_csv_xy(csv_file: &str) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut rdr = csv::Reader::from_path(csv_file)?;
    //let mut r = csv::StringRecord::new();
    let mut x = Vec::<f64>::new();
    let mut y = Vec::<f64>::new();
    for r in rdr.records() {
        if let Ok(r) = r {
            x.push(r[0].parse::<f64>().unwrap());
            y.push(r[1].parse::<f64>().unwrap());
        } else {
            break
        }
    }
    Ok((x,y))
}

pub fn read_csv_uxy(csv_file: &str) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    //let mut rdr = csv::Reader::from_path(csv_file)?;
    let mut rdr = 
        ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#'))
        .trim(csv::Trim::All)
        .from_path(csv_file)?;
    //let mut r = csv::StringRecord::new();
    let mut u = Vec::<f64>::new();
    let mut x = Vec::<f64>::new();
    let mut y = Vec::<f64>::new();
    for r in rdr.records() {
        if let Ok(r) = r {
            u.push(r[0].parse::<f64>().unwrap());
            x.push(r[1].parse::<f64>().unwrap());
            y.push(r[2].parse::<f64>().unwrap());
        } else {
            break
        }
    }
    Ok((u, x,y))
}

pub fn write_csv_xy(csv_file: &str, x: &[f64], y: &[f64]) -> Result<()> {
    let mut wtr = csv::Writer::from_path(csv_file)?;
    wtr.write_record(&["wl[nm]", "spd[-]"])?;
    for (&x,&y) in x.iter().zip(y.iter()) {
        
        wtr.write_record(&[format!("{:.2}",x), format!("{:.4}",y)])?
    }
    wtr.flush()?;
    Ok(())
}