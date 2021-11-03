use super::{Result, };
use csv::ReaderBuilder;

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