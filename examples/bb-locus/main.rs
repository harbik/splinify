
use splinify::{ParameterSplineCurveFit, Result, read_csv_uxy};

fn main() -> Result<()> {

    std::env::set_current_dir(std::path::Path::new(file!()).parent().unwrap())?;   

    let (u, x,y) =  read_csv_uxy("bb-locus-31.csv")? ;
    let mut xy: Vec<f64> = Vec::with_capacity(x.len()*2);
    x.iter().zip(y.iter()).for_each(|(&x,&y)|{xy.push(x); xy.push(y);});

    let d = 
        ParameterSplineCurveFit::<3,2>
            ::new(u.clone(), xy.clone())?
            .smoothing_spline(1E-5)?;

    let json = serde_json::to_string_pretty(&d)?;
    println!("{}", json);

    d.plot("fit.png", (2000,1000))?;

    Ok(())
}

