
use splinify::{ParameterSplineCurveFit, Result, read_csv_uxy, SplineCurveData};

fn main() -> Result<()> {

    let example_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(std::path::Path::new(file!()).parent().unwrap());

    let csv_path = example_dir.join("bb-locus-31.csv");
    let (u, x,y) = read_csv_uxy(csv_path.to_str().expect("non-UTF-8 path"))?;
    let mut xy: Vec<f64> = Vec::with_capacity(x.len()*2);
    x.iter().zip(y.iter()).for_each(|(&x,&y)|{xy.push(x); xy.push(y);});

    let d =
        ParameterSplineCurveFit::<3,2>
            ::new(u.clone(), xy.clone())?
            .smoothing_spline(5E-4)?;

    let json = serde_json::to_string_pretty(&SplineCurveData::from(&d))?;
    println!("{}", json);

    let plot_path = example_dir.join("fit.png");
    d.plot(plot_path.to_str().expect("non-UTF-8 path"), (2000,1000))?;

    Ok(())
}

