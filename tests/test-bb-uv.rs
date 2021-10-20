
extern crate dierckx;
use dierckx::{ParameterCurveSplineFit, FitResult, read_csv_uxy, plot};


#[test]
fn test_bb() -> FitResult<()> {

    let (u, x,y) =  read_csv_uxy("tests/data/bb-locus-31.csv")? ;
    let mut xy: Vec<f64> = Vec::with_capacity(x.len()*2);
    
    x.iter().zip(y.iter()).for_each(|(&x,&y)|{xy.push(x); xy.push(y);});

    let d = 
    ParameterCurveSplineFit::<3,2>::new(u.clone(), xy.clone())?
       // .begin_constraints([ [u[0]] ])?
       // .end_constraints([ [u[u.len()-1]] ])?
      //  .smoothing_spline_optimize(1E-5, 0.8,  |k,_,rms,_| {println!("{}/{}", k, rms*1000.0); k>60}, None)?;
        .smoothing_spline(1E-5)?;
    println!("{}:{}", d.t.len(), d.e_rms.unwrap());

    //let y_s = d.evaluate(&x)?;
   // plot("tests/img/fit.png", x, y, d)?;

    Ok(())
}

