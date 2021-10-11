
use dierckx::{CubicCurveFit, CurveFit, ConstrainedSpline, Result};
use nalgebra::{Const, DMatrix, DVector, Dynamic, OMatrix};

#[test]
fn test_smoothing() -> Result<()> {
    use plotters::prelude::*;
    let (x,y) = xys();

    let mut d = CurveFit::<3>::new(x.clone(), y.clone());
    let (d, f) = d.smoothing_spline(1.25)?;
    println!("knots {:?}", d.tc.t);
    println!("number of knots: {}", d.tc.t.len());
    println!("fp: {:?}", f);


    let ymax = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    let root = BitMapBackend::new("led_smooth_spline.png", (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(50)
        .caption("Smoothing Spline", ("sans-serif",24))
        .build_cartesian_2d(380f64..780f64, -ymax/5.0..ymax)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y.iter().cloned()), &BLUE))?;

    let y_s = d.as_ref().evaluate(&x)?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y_s.iter().cloned()), &BLACK))?;
    chart.draw_series(d.tc.t.iter().cloned().zip(d.tc.c.iter().cloned()).map(|xy| Circle::new(xy, 3, RED.filled())))?;

    root.present()?;

    Ok(())
}

#[test]
fn test_cardinal() -> Result<()> {
    use plotters::prelude::*;
    let (x,y) = xys();

    let d = CubicCurveFit::new(x.clone(), y.clone());
    let (tc, f) = d.cardinal_spline(10.0)?;
    println!("knots {:?}", tc.t);
    println!("number of knots: {}", tc.t.len());
    println!("fp: {:?}", f);


    let ymax = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    let root = BitMapBackend::new("led_card_spline.png", (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(50)
        .caption("Cardinal Spline", ("sans-serif",24))
        .build_cartesian_2d(380f64..780f64, -ymax/5.0..ymax)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y.iter().cloned()), &BLUE))?;

    let y_s = tc.evaluate(&x)?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y_s.iter().cloned()), &BLACK))?;
    chart.draw_series(tc.t.iter().cloned().zip(tc.c.iter().cloned()).map(|xy| Circle::new(xy, 3, RED.filled())))?;

    root.present()?;

    Ok(())
}

#[test]
fn test_iter(){
    use std::iter::repeat;
    let v:Vec<f64> = repeat(5.0).take(81).scan(380.0, |s, dx| {let t=*s; *s+=dx; Some(t)}).collect();
    println!("{:?}", v);

    let w:Vec<f64> = repeat(5.0).scan(380.0, |s, dx| {let t=*s; *s+=dx; if t>780.0 {None } else {Some(t)}}).collect();
    println!("{:?}", w);
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
    let d = CurveFit::<3>::new(x_data, y_data);
    let y_fit = d.interpolating_spline()?.evaluate(&x)?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y_fit.iter().cloned()), &HSLColor(0.5, 1.0, 0.5)))?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y.iter().cloned()), &BLACK))?;

    /*
    let y_s = tc.values(&x)?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y_s.iter().cloned()), &BLACK))?;
    chart.draw_series(tc.t.iter().cloned().zip(tc.c.iter().cloned()).map(|xy| Circle::new(xy, 3, RED.filled())))?;
     */

    root.present()?;


    Ok(())
}


#[test]
fn constrained_cardinal_spline() -> Result<()> {
    use plotters::prelude::*;
    let (u,y) = xys();
    let u = DVector::<f64>::from_vec(u);
    let xn = OMatrix::<f64, Const<1>, Dynamic>::from_vec(y);
    let xb = OMatrix::<f64, Dynamic, Const<1>>::from_vec(vec![y[0], 0.0, 0.0]); // first and second derivatives 0.0
    let xe = OMatrix::<f64, Dynamic, Const<1>>::from_vec(vec![y[y.len()-1], 0.0, 0.0]);
    let cs = ConstrainedSpline::<3,1>::new(u, xn, xb, xe);

    let (tc, f) = cs.cardinal_spline(10.0)?;
    println!("knots {:?}", tc.t);
    println!("number of knots: {}", tc.t.len());
    println!("fp: {:?}", f);

    let ymax = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    let root = BitMapBackend::new("led_card_constrained_spline.png", (2000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(50)
        .caption("Constrained Cardinal Spline", ("sans-serif",24))
        .build_cartesian_2d(380f64..780f64, -ymax/5.0..ymax)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y.iter().cloned()), &BLUE))?;

    let y_s = tc.evaluate(&x)?;
    chart.draw_series(LineSeries::new(x.iter().cloned().zip(y_s.iter().cloned()), &BLACK))?;
    chart.draw_series(tc.t.iter().cloned().zip(tc.c.iter().cloned()).map(|xy| Circle::new(xy, 3, RED.filled())))?;

    root.present()?;

    Ok(())
}


fn xys() -> (Vec<f64>, Vec<f64>) {

    let x: Vec<f64> = vec![
        380.68753, 381.0921, 381.49673, 381.9014, 382.30609, 382.71085, 383.11566, 383.52051, 383.92538, 384.33032,
        384.73532, 385.14035, 385.54541, 385.95053, 386.35571, 386.76093, 387.16617, 387.57147, 387.97681, 388.3822,
        388.78766, 389.19312, 389.59866, 390.00424, 390.40985, 390.81552, 391.22122, 391.62698, 392.03278, 392.43863,
        392.84451, 393.25046, 393.65643, 394.06247, 394.46854, 394.87466, 395.28082, 395.68704, 396.09329, 396.4996,
        396.90594, 397.31232, 397.71878, 398.12524, 398.53177, 398.93835, 399.34497, 399.75162, 400.15833, 400.56509,
        400.97189, 401.37872, 401.78561, 402.19257, 402.59955, 403.00656, 403.41364, 403.82077, 404.22791, 404.63513,
        405.04239, 405.44968, 405.85703, 406.2644, 406.67184, 407.07932, 407.48685, 407.89441, 408.30203, 408.70969,
        409.11737, 409.52515, 409.93292, 410.34076, 410.74866, 411.15659, 411.56458, 411.9726, 412.38065, 412.78876,
        413.1969, 413.6051, 414.01337, 414.42166, 414.82999, 415.23837, 415.64679, 416.05527, 416.46378, 416.87234,
        417.28094, 417.68961, 418.0983, 418.50705, 418.91583, 419.32468, 419.73355, 420.14246, 420.55142, 420.96045,
        421.36951, 421.77859, 422.18774, 422.59695, 423.0062, 423.41547, 423.8248, 424.23416, 424.64359, 425.05304,
        425.46255, 425.8721, 426.28171, 426.69135, 427.10101, 427.51077, 427.92053, 428.33035, 428.74023, 429.15015,
        429.56009, 429.97009, 430.38013, 430.79022, 431.20038, 431.61053, 432.02078, 432.43103, 432.84137, 433.25171,
        433.66211, 434.07257, 434.48306, 434.89362, 435.3042, 435.71481, 436.12549, 436.53622, 436.94696, 437.35779,
        437.76865, 438.17953, 438.59048, 439.00146, 439.41251, 439.82358, 440.23468, 440.64584, 441.05707, 441.46832,
        441.87961, 442.29095, 442.70236, 443.1138, 443.52527, 443.9368, 444.34836, 444.75998, 445.17163, 445.58334,
        445.99509, 446.40686, 446.8187, 447.23059, 447.64252, 448.05447, 448.46649, 448.87854, 449.29065, 449.70279,
        450.11499, 450.52722, 450.93951, 451.35184, 451.76422, 452.17664, 452.58908, 453.00159, 453.41412, 453.82672,
        454.23938, 454.65204, 455.06476, 455.47754, 455.89035, 456.30322, 456.71613, 457.12906, 457.54205, 457.95511,
        458.36819, 458.78131, 459.19449, 459.6077, 460.02097, 460.43427, 460.8476, 461.26099, 461.67444, 462.08792,
        462.50143, 462.91501, 463.32861, 463.74228, 464.15598, 464.56973, 464.98352, 465.39737, 465.81125, 466.22516,
        466.63913, 467.05313, 467.46719, 467.88129, 468.29544, 468.70963, 469.12387, 469.53815, 469.95245, 470.36682,
        470.78125, 471.19571, 471.6102, 472.02475, 472.43933, 472.85394, 473.26862, 473.68335, 474.09811, 474.51291,
        474.92776, 475.34265, 475.7576, 476.17258, 476.58762, 477.00269, 477.41782, 477.83298, 478.24817, 478.66342,
        479.0787, 479.49405, 479.90942, 480.32486, 480.74033, 481.15582, 481.57138, 481.98697, 482.40262, 482.8183,
        483.23404, 483.64981, 484.06564, 484.48151, 484.8974, 485.31335, 485.72934, 486.14539, 486.56146, 486.9776,
        487.39377, 487.81, 488.22623, 488.64255, 489.0589, 489.47528, 489.89172, 490.3082, 490.7247, 491.14127, 491.55789,
        491.97452, 492.39124, 492.80795, 493.22476, 493.64157, 494.05844, 494.47537, 494.89233, 495.30933, 495.72638,
        496.14346, 496.56058, 496.97775, 497.39499, 497.81226, 498.22955, 498.64691, 499.0643, 499.48172, 499.8992,
        500.31674, 500.73431, 501.15192, 501.56958, 501.98727, 502.40503, 502.82281, 503.24063, 503.65851, 504.07642,
        504.49438, 504.91238, 505.33044, 505.74854, 506.16666, 506.58484, 507.00305, 507.42133, 507.83963, 508.25797,
        508.67636, 509.09482, 509.51328, 509.93179, 510.35037, 510.76898, 511.18765, 511.60632, 512.02509, 512.44385,
        512.86267, 513.28156, 513.7005, 514.11945, 514.53845, 514.95746, 515.37659, 515.79572, 516.21484, 516.63409,
        517.05334, 517.47266, 517.89197, 518.3114, 518.73083, 519.15027, 519.56982, 519.98938, 520.40894, 520.82861,
        521.24829, 521.66803, 522.08783, 522.50763, 522.92749, 523.34741, 523.76733, 524.18732, 524.60736, 525.02747,
        525.44757, 525.86774, 526.28796, 526.70819, 527.12848, 527.54883, 527.96918, 528.38959, 528.81006, 529.23059,
        529.65112, 530.07172, 530.49237, 530.91302, 531.33374, 531.75452, 532.17535, 532.59619, 533.01709, 533.43799,
        533.85901, 534.28003, 534.70105, 535.12219, 535.54333, 535.96454, 536.38574, 536.80701, 537.22833, 537.64972,
        538.07111, 538.49255, 538.91406, 539.33557, 539.75714, 540.17877, 540.60046, 541.02216, 541.44391, 541.86566,
        542.28748, 542.70935, 543.13129, 543.55322, 543.97522, 544.39728, 544.8194, 545.24152, 545.6637, 546.08588,
        546.50812, 546.93042, 547.35278, 547.77515, 548.19757, 548.62006, 549.04254, 549.46509, 549.8877, 550.31036,
        550.73303, 551.15576, 551.57849, 552.00128, 552.42413, 552.84705, 553.26996, 553.69293, 554.11597, 554.539,
        554.96216, 555.38525, 555.80847, 556.23169, 556.65497, 557.07825, 557.50165, 557.92499, 558.34845, 558.77191,
        559.19543, 559.61902, 560.0426, 560.46625, 560.88995, 561.31372, 561.73749, 562.16125, 562.58514, 563.00903,
        563.43298, 563.85699, 564.28101, 564.70508, 565.12915, 565.55334, 565.97754, 566.40173, 566.82605, 567.25037,
        567.67468, 568.09912, 568.52356, 568.94806, 569.37256, 569.79712, 570.22174, 570.64642, 571.07111, 571.49585,
        571.92059, 572.3454, 572.77026, 573.19519, 573.62012, 574.0451, 574.47015, 574.8952, 575.32031, 575.74548,
        576.17065, 576.59589, 577.02118, 577.44647, 577.87189, 578.29724, 578.72272, 579.14819, 579.57373, 579.99927,
        580.42487, 580.85052, 581.27625, 581.70197, 582.12775, 582.55353, 582.97943, 583.40533, 583.83124, 584.2572,
        584.68323, 585.10931, 585.5354, 585.96155, 586.38776, 586.81396, 587.24023, 587.66656, 588.0929, 588.51929,
        588.94574, 589.37219, 589.79871, 590.22528, 590.65186, 591.07849, 591.50519, 591.93195, 592.3587, 592.78546,
        593.21234, 593.63922, 594.06616, 594.4931, 594.9201, 595.34717, 595.77423, 596.20135, 596.62854, 597.05579,
        597.48303, 597.91034, 598.33765, 598.76501, 599.19244, 599.61987, 600.04742, 600.47491, 600.90253, 601.33014,
        601.75781, 602.18549, 602.61322, 603.04102, 603.46887, 603.89673, 604.32458, 604.75256, 605.18054, 605.60858,
        606.03662, 606.46472, 606.89288, 607.32111, 607.74933, 608.17761, 608.6059, 609.03424, 609.46265, 609.89105,
        610.31952, 610.74805, 611.17657, 611.60516, 612.03381, 612.46252, 612.89124, 613.31995, 613.74878, 614.17761,
        614.60645, 615.0354, 615.46436, 615.89331, 616.32239, 616.75146, 617.18054, 617.60974, 618.03894, 618.46814,
        618.8974, 619.32672, 619.7561, 620.18549, 620.61493, 621.04443, 621.47394, 621.9035, 622.33307, 622.7627, 623.19238,
        623.62213, 624.05188, 624.48163, 624.9115, 625.34137, 625.7713, 626.20123, 626.63123, 627.06128, 627.49133,
        627.92145, 628.35162, 628.7818, 629.21204, 629.64233, 630.07263, 630.50299, 630.93335, 631.36377, 631.79425,
        632.22479, 632.65533, 633.08594, 633.51654, 633.9472, 634.37793, 634.80865, 635.23944, 635.67029, 636.10114,
        636.53204, 636.96301, 637.39398, 637.82501, 638.25604, 638.68719, 639.11829, 639.5495, 639.98071, 640.41199,
        640.84326, 641.2746, 641.70599, 642.13739, 642.56885, 643.00037, 643.43188, 643.86346, 644.29504, 644.72668,
        645.15839, 645.59015, 646.02191, 646.45367, 646.88556, 647.31744, 647.74933, 648.18134, 648.61334, 649.04535,
        649.47742, 649.90955, 650.34174, 650.77393, 651.20612, 651.63843, 652.07074, 652.50305, 652.93549, 653.36786,
        653.80035, 654.23285, 654.66541, 655.09796, 655.53058, 655.96326, 656.39594, 656.82867, 657.26147, 657.69427,
        658.12714, 658.56, 658.99298, 659.4259, 659.85895, 660.29199, 660.72504, 661.1582, 661.59131, 662.02454, 662.45776,
        662.89105, 663.32434, 663.75769, 664.1911, 664.62451, 665.05798, 665.49152, 665.92505, 666.35858, 666.79224,
        667.22589, 667.65961, 668.09332, 668.5271, 668.96088, 669.39471, 669.82861, 670.26251, 670.69647, 671.13049,
        671.56451, 671.9986, 672.43274, 672.86688, 673.30103, 673.73529, 674.16956, 674.60382, 675.03821, 675.47253,
        675.90698, 676.34143, 676.77588, 677.21045, 677.64502, 678.07959, 678.51422, 678.94891, 679.38361, 679.81836,
        680.25317, 680.68799, 681.12286, 681.55774, 681.99268, 682.42767, 682.86267, 683.29773, 683.73285, 684.16797,
        684.60309, 685.03833, 685.47357, 685.90881, 686.34412, 686.77948, 687.2149, 687.65027, 688.08575, 688.52124,
        688.95679, 689.39233, 689.82794, 690.26361, 690.69928, 691.13501, 691.5708, 692.00659, 692.44238, 692.87823,
        693.31415, 693.75012, 694.1861, 694.62207, 695.05817, 695.49426, 695.93036, 696.36652, 696.80273, 697.23895,
        697.67523, 698.11151, 698.54785, 698.98425, 699.42065, 699.85712, 700.29364, 700.73016, 701.16669, 701.60333,
        702.03998, 702.47662, 702.91333, 703.3501, 703.78687, 704.22369, 704.66052, 705.09741, 705.53436, 705.97131,
        706.40833, 706.84534, 707.28241, 707.71954, 708.15668, 708.59387, 709.03107, 709.46832, 709.90558, 710.34296,
        710.78027, 711.21771, 711.65509, 712.09259, 712.53009, 712.96765, 713.40521, 713.84283, 714.28046, 714.71814,
        715.15588, 715.59363, 716.03137, 716.46924, 716.9071, 717.34497, 717.7829, 718.22089, 718.65887, 719.09692,
        719.53497, 719.97308, 720.41125, 720.84943, 721.28766, 721.72589, 722.16418, 722.60254, 723.04089, 723.47925,
        723.91766, 724.35614, 724.79462, 725.23315, 725.67175, 726.11035, 726.54901, 726.98767, 727.42639, 727.86511,
        728.30389, 728.74274, 729.18158, 729.62042, 730.05939, 730.49829, 730.93732, 731.37634, 731.81537, 732.25446,
        732.6936, 733.13275, 733.57196, 734.01117, 734.45044, 734.88977, 735.3291, 735.76849, 736.20789, 736.64734,
        737.08679, 737.52631, 737.96582, 738.4054, 738.84503, 739.28467, 739.72437, 740.16406, 740.60382, 741.04358,
        741.4834, 741.92328, 742.36316, 742.8031, 743.24304, 743.68304, 744.12305, 744.56311, 745.00317, 745.4433,
        745.88348, 746.32367, 746.76392, 747.20416, 747.64447, 748.08478, 748.52515, 748.96552, 749.40594, 749.84644,
        750.28693, 750.72748, 751.16803, 751.60864, 752.04926, 752.48993, 752.9306, 753.37134, 753.81213, 754.25293,
        754.69373, 755.13458, 755.5755, 756.01642, 756.4574, 756.89844, 757.33942, 757.78052, 758.22162, 758.66272,
        759.10388, 759.5451, 759.98633, 760.42761, 760.8689, 761.31024, 761.75159, 762.19299, 762.63446, 763.07587,
        763.5174, 763.95892, 764.40051, 764.8421, 765.28369, 765.72534, 766.16705, 766.60876, 767.05054, 767.49231,
        767.93414, 768.37604, 768.81793, 769.25983, 769.70178, 770.1438, 770.58582, 771.02783, 771.46991, 771.91205,
        772.35419, 772.79639, 773.23859, 773.68085, 774.12311, 774.56543, 775.00775, 775.45013, 775.89258, 776.33502,
        776.77747, 777.21997, 777.66248, 778.10504, 778.54767, 778.9903, 779.43298, 779.87567
    ];
 
    let y: Vec<f64> = vec![ 
        0.0, 1.88E-05, 1.93E-05, 1.53E-05, 2.78E-05, 2.53E-06, -6.27E-06, 1.28E-05, -1.61E-05, 1.61E-06,
        -3.04E-06, 9.98E-06, 7.33E-07, 4.86E-06, 0.0, 5.17E-06, 2.67E-06, 5.76E-06, 7.11E-06, -2.01E-06, -4.80E-06,
        3.90E-06, 6.19E-07, 1.94E-06, 3.97E-06, 9.62E-06, 6.71E-06, 1.18E-06, 1.39E-06, 6.54E-06, -2.89E-06, 1.90E-05,
        1.42E-05, -4.69E-05, -6.73E-05, 1.04E-05, 2.04E-05, 6.73E-05, 0.0, 2.00E-05, 9.22E-06, 7.56E-06, 9.56E-06,
        2.38E-05, 3.17E-05, 3.91E-05, 2.41E-05, 1.73E-05, 1.34E-05, 1.01E-05, 2.57E-05, 1.51E-05, 2.47E-05, 2.90E-05,
        4.12E-05, 7.98E-06, 3.06E-05, 2.24E-05, 2.16E-05, 7.18E-05, 3.66E-05, 2.23E-05, 4.58E-05, -2.81E-06, 5.66E-05,
        5.11E-05, 0.0, 5.58E-05, 4.25E-05, 7.51E-05, 6.77E-05, 4.37E-05, 5.00E-05, 6.34E-05, 9.15E-05, 6.69E-05, 3.88E-05,
        6.83E-05, 9.71E-05, 0.000111881, 0.000164116, 7.42E-05, 7.97E-05, 0.000130494, 0.000175154, 0.000164443,
        0.000153293, 0.000150988, 0.000178454, 0.000202917, 0.000214779, 0.000246923, 0.000250916, 0.000263674,
        0.000230919, 0.00025795, 0.000281972, 0.000339531, 0.000309252, 0.000370167, 0.000352638, 0.000411623,
        0.000478442, 0.000459483, 0.000496909, 0.000509134, 0.00052195, 0.000545041, 0.000609703, 0.000658039,
        0.000696053, 0.000778006, 0.000772076, 0.000864727, 0.00083969, 0.000947801, 0.000962223, 0.001004805,
        0.001151018, 0.001084543, 0.001174869, 0.00128682, 0.001316528, 0.001477899, 0.001572519, 0.001566959,
        0.001657925, 0.001698911, 0.001798549, 0.001929285, 0.001964873, 0.002119973, 0.002240692, 0.002453445,
        0.002565855, 0.002635181, 0.00276223, 0.002792229, 0.003029751, 0.003056762, 0.003355537, 0.003607889,
        0.003667268, 0.003976232, 0.004007224, 0.004285469, 0.004531726, 0.004663665, 0.005080742, 0.005251654,
        0.005612905, 0.005740715, 0.006057187, 0.006300707, 0.006652835, 0.006849675, 0.00695273, 0.007346191,
        0.007560199, 0.007746336, 0.007988869, 0.008090429, 0.008169225, 0.008410285, 0.008622333, 0.008566361,
        0.008709354, 0.00859989, 0.0084324, 0.008339626, 0.008310358, 0.008128845, 0.007872418, 0.007810225,
        0.007556766, 0.007320587, 0.006955465, 0.006547886, 0.006464886, 0.006179383, 0.00596529, 0.00579608,
        0.005488984, 0.005295227, 0.005081343, 0.004919291, 0.004735145, 0.004729892, 0.004566751, 0.004355087,
        0.004294773, 0.004086271, 0.00397839, 0.003894307, 0.003853339, 0.003798506, 0.003747407, 0.003602932,
        0.003485786, 0.003508372, 0.003419269, 0.003465368, 0.00339712, 0.003309464, 0.003204965, 0.003144598,
        0.003193706, 0.003047522, 0.003022552, 0.003065556, 0.002836976, 0.002771146, 0.002845049, 0.002759388,
        0.002609941, 0.002621447, 0.002548149, 0.002524689, 0.00239757, 0.002403131, 0.002374261, 0.00230472,
        0.00229082, 0.002301689, 0.00231146, 0.002382841, 0.002303139, 0.002275317, 0.002236171, 0.002303688,
        0.002312744, 0.002377662, 0.002318854, 0.002382924, 0.002418139, 0.002474635, 0.002457994, 0.002536467,
        0.002628063, 0.002577146, 0.002708557, 0.002798099, 0.002859973, 0.002993595, 0.002946404, 0.00300678,
        0.003015438, 0.003126139, 0.003254354, 0.003248937, 0.003354867, 0.003397296, 0.003488809, 0.003464737,
        0.003659145, 0.003656715, 0.003780805, 0.003713943, 0.003874292, 0.003932341, 0.004036786, 0.004053155,
        0.004093671, 0.004247804, 0.004231356, 0.004300484, 0.004267295, 0.004340806, 0.004414844, 0.004591048,
        0.004542998, 0.004668816, 0.004696992, 0.004667009, 0.00478497, 0.004921008, 0.005068902, 0.005017227,
        0.004960955, 0.004991987, 0.004999028, 0.005050673, 0.005113024, 0.005250935, 0.005161613, 0.005261939,
        0.005256636, 0.005202534, 0.005359141, 0.005311876, 0.005383568, 0.005290797, 0.005296554, 0.005357758,
        0.005115685, 0.005166963, 0.005241775, 0.005254473, 0.005292706, 0.005355761, 0.005325135, 0.005305052,
        0.005336053, 0.005300144, 0.005338578, 0.00535152, 0.005239748, 0.005230251, 0.005314321, 0.005198107,
        0.005311927, 0.005281714, 0.00528349, 0.005182703, 0.005162988, 0.005250734, 0.005249176, 0.005267455,
        0.005183063, 0.00514709, 0.005118011, 0.005075561, 0.005164674, 0.005145983, 0.005020297, 0.005050018,
        0.005112761, 0.00504338, 0.004983208, 0.005063999, 0.004991397, 0.005102702, 0.00507015, 0.005093383,
        0.005033624, 0.004994579, 0.005063294, 0.005034323, 0.004962576, 0.004937514, 0.004965405, 0.005029959,
        0.004947572, 0.004899276, 0.005058983, 0.004981111, 0.004955026, 0.004966856, 0.004968359, 0.005019622,
        0.004876757, 0.004961764, 0.004882302, 0.004952333, 0.004928261, 0.00493845, 0.004861353, 0.004936358,
        0.004883653, 0.004963166, 0.004985452, 0.004916546, 0.004951123, 0.00493857, 0.004940946, 0.005049864,
        0.005037436, 0.004902182, 0.005039033, 0.005107977, 0.004963665, 0.004981274, 0.004963765, 0.005044605,
        0.005055717, 0.005117237, 0.005035179, 0.005066915, 0.005186214, 0.005074454, 0.005054481, 0.005031101,
        0.005070759, 0.005100157, 0.00514338, 0.00516031, 0.005222086, 0.005177075, 0.005305445, 0.005110662,
        0.005140294, 0.005122001, 0.005128784, 0.005278984, 0.005257196, 0.005306066, 0.005218274, 0.005194386,
        0.005324285, 0.005298422, 0.005288776, 0.005310475, 0.005320677, 0.005356461, 0.005349436, 0.005337677,
        0.005308257, 0.005346688, 0.005360222, 0.005334019, 0.005402685, 0.005461058, 0.005412766, 0.005454604,
        0.005490141, 0.005451533, 0.005395066, 0.005507046, 0.005414781, 0.005467052, 0.005499635, 0.005414355,
        0.005515161, 0.005446905, 0.00549317, 0.005574624, 0.005450216, 0.00554141, 0.005567138, 0.005500748,
        0.005499675, 0.005603284, 0.005534298, 0.005468897, 0.005541022, 0.005500267, 0.005663436, 0.005593217,
        0.005471381, 0.005562747, 0.005611538, 0.00552141, 0.005473552, 0.005421688, 0.005458955, 0.00551314,
        0.005595853, 0.005622043, 0.005613429, 0.005562971, 0.005693684, 0.005675371, 0.005745274, 0.005671814,
        0.005603876, 0.005513085, 0.005630184, 0.005518861, 0.005494108, 0.005457217, 0.005429008, 0.005516675,
        0.005556013, 0.00557938, 0.005647647, 0.005522619, 0.005645283, 0.005527002, 0.005542572, 0.005692406,
        0.005748662, 0.005575783, 0.005692891, 0.005545016, 0.005633183, 0.005577462, 0.005562638, 0.005545202,
        0.005624392, 0.005557599, 0.005543399, 0.005556915, 0.005639771, 0.005574035, 0.005651176, 0.005492705,
        0.005533426, 0.005720685, 0.005511465, 0.005522836, 0.005523132, 0.005544914, 0.005716985, 0.005600734,
        0.005707, 0.005515064, 0.005651278, 0.005720394, 0.00566799, 0.005691845, 0.005551718, 0.005606448, 0.005651456,
        0.005573341, 0.00568448, 0.005694585, 0.005695297, 0.005714192, 0.005737385, 0.005740552, 0.00576112,
        0.005628755, 0.005728273, 0.005710599, 0.005804245, 0.005775043, 0.005812277, 0.005839625, 0.005776786,
        0.005704135, 0.005824538, 0.005809952, 0.005922952, 0.00585627, 0.005770726, 0.005780453, 0.005850903,
        0.005894681, 0.005904093, 0.006083122, 0.005942714, 0.005953457, 0.006005428, 0.005863444, 0.006040091,
        0.006096551, 0.005973846, 0.005933246, 0.005874732, 0.005897629, 0.006080405, 0.005945555, 0.006072457,
        0.006162772, 0.006182151, 0.006192696, 0.006183777, 0.006192586, 0.006173669, 0.006278231, 0.006279832,
        0.006246925, 0.006291591, 0.006217598, 0.006289922, 0.006313355, 0.006255723, 0.006408779, 0.006295424,
        0.006276653, 0.006387811, 0.0063439, 0.006241645, 0.006414533, 0.006326277, 0.006439147, 0.00652289,
        0.006414831, 0.006397001, 0.006542845, 0.006555151, 0.00662147, 0.006499148, 0.006559197, 0.006551219,
        0.006651951, 0.006721432, 0.006660607, 0.006611188, 0.006562711, 0.006575733, 0.006744721, 0.006569618,
        0.006678002, 0.006761916, 0.006802026, 0.006734482, 0.006674493, 0.006623236, 0.006635374, 0.006645309,
        0.006676288, 0.006715584, 0.00681364, 0.006794264, 0.006778293, 0.006680547, 0.006722165, 0.006701445,
        0.006815904, 0.006664699, 0.006840439, 0.006860296, 0.006649026, 0.006699006, 0.00662937, 0.006685491,
        0.006727609, 0.006681109, 0.006761219, 0.006614255, 0.006723322, 0.006588429, 0.006816092, 0.006761434,
        0.006759277, 0.00679049, 0.006767593, 0.006749422, 0.00667812, 0.006623796, 0.006674083, 0.006500121,
        0.00660879, 0.006499457, 0.006487289, 0.00654175, 0.006584444, 0.006541966, 0.006476061, 0.006540419,
        0.00642384, 0.006467193, 0.006334045, 0.006385906, 0.006343367, 0.006422316, 0.006291378, 0.006397335,
        0.006332044, 0.00627681, 0.006309, 0.006223083, 0.006082194, 0.006181721, 0.006219845, 0.006159267, 0.006096018,
        0.006107702, 0.005936448, 0.005945856, 0.005966925, 0.005909676, 0.006054096, 0.006022497, 0.006005811,
        0.005919962, 0.005973577, 0.005892728, 0.005872133, 0.005788168, 0.005699977, 0.005721101, 0.005713616,
        0.005681178, 0.005723145, 0.005623939, 0.005627688, 0.005495157, 0.005619369, 0.005456762, 0.005513314,
        0.005523157, 0.005423628, 0.005388115, 0.005331789, 0.00528674, 0.005329008, 0.005223737, 0.005267281,
        0.005062176, 0.005125968, 0.005096705, 0.005017403, 0.005104591, 0.005045589, 0.004975312, 0.004882117,
        0.00487338, 0.004890399, 0.004893463, 0.004834899, 0.004805757, 0.004654402, 0.004808625, 0.004720221,
        0.004573786, 0.004596691, 0.004495814, 0.004516444, 0.004453201, 0.004453518, 0.004385377, 0.004329897,
        0.004290876, 0.004301486, 0.004230069, 0.004209026, 0.004241599, 0.004138981, 0.004127437, 0.004059171,
        0.00407716, 0.003982399, 0.004028507, 0.003912763, 0.003995082, 0.003840026, 0.003858029, 0.003781074,
        0.003845232, 0.003649395, 0.003667857, 0.003642053, 0.003506956, 0.00358287, 0.003546082, 0.003506902,
        0.003584901, 0.00350646, 0.003392768, 0.003348357, 0.003335634, 0.003168942, 0.003341862, 0.003227355,
        0.003249683, 0.003207419, 0.003254827, 0.00312274, 0.003158623, 0.003159637, 0.003017083, 0.003086396,
        0.003060954, 0.003021917, 0.002893632, 0.002813183, 0.002762099, 0.002768613, 0.00283683, 0.002792178,
        0.002787752, 0.002778684, 0.002689691, 0.002650511, 0.002640715, 0.002639378, 0.002479359, 0.002588093,
        0.002522736, 0.002452144, 0.002404659, 0.002371471, 0.002429273, 0.002429867, 0.002368371, 0.002325418,
        0.002314183, 0.002324532, 0.002191607, 0.002222679, 0.002209439, 0.002172881, 0.002225934, 0.002150524,
        0.002123775, 0.002071722, 0.002013805, 0.001980454, 0.00194621, 0.001983818, 0.001903308, 0.001950574,
        0.001956442, 0.001833421, 0.00185855, 0.001821971, 0.001794478, 0.001786102, 0.001835611, 0.001757625,
        0.001697285, 0.001657263, 0.001682565, 0.00163971, 0.001597143, 0.001596157, 0.001623898, 0.00163077,
        0.001615173, 0.001565703, 0.001562672, 0.001473475, 0.001498031, 0.001485372, 0.0013924, 0.00142923,
        0.001392181, 0.001393245, 0.001370315, 0.001386469, 0.001364985, 0.001273426, 0.001266168, 0.001346913,
        0.001293797, 0.001271513, 0.001179394, 0.001223479, 0.001198157, 0.00120214, 0.001196996, 0.001223884,
        0.001111913, 0.001132051, 0.001102065, 0.001074705, 0.00111025, 0.001068429, 0.001068823, 0.001013254,
        0.001072231, 0.000990055, 0.000997406, 0.000981928, 0.000944022, 0.000978616, 0.000959139, 0.000931135,
        0.000986568, 0.000928374, 0.000903904, 0.000911729, 0.000864837, 0.000895501, 0.000873878, 0.000947796,
        0.000825282, 0.000830474, 0.000820672, 0.000838655, 0.000803606, 0.000780193, 0.00081895, 0.000746315,
        0.000764061, 0.000725589, 0.000713775, 0.000699769, 0.000702114, 0.000689355, 0.000702181, 0.000664502,
        0.000711214, 0.000640277, 0.000642958, 0.000651476, 0.000674373, 0.000584814, 0.000614535, 0.000629801,
        0.000623424, 0.000569846, 0.000610521, 0.000607921, 0.000532349, 0.000589896, 0.000555316, 0.000545977,
        0.000521261, 0.000510669, 0.000506295, 0.000527331, 0.000514903, 0.000523479, 0.000485247, 0.000447268,
        0.000466381, 0.000480821, 0.000483258, 0.000464119, 0.00044007, 0.000458776, 0.00043584, 0.000442365,
        0.000425599, 0.000438773, 0.000444373, 0.000427679, 0.000411666, 0.000405452, 0.000363392, 0.000374099,
        0.000381374, 0.00038405, 0.000381639, 0.000346745, 0.000389552, 0.00035469, 0.000337453, 0.00033749,
        0.000388477, 0.00035062, 0.000325969, 0.000366554, 0.00034222, 0.000287808, 0.000275384, 0.00028275,
        0.000317003, 0.000326354, 0.000271796, 0.00029947, 0.000327808, 0.000270022, 0.00033277, 0.000243539,
        0.000278997, 0.000286905, 0.000262354, 0.000272509, 0.000245036, 0.000265087, 0.000253833, 0.000263894,
        0.000242344, 0.000236343, 0.000239197, 0.000228471, 0.00030493 
    ];
    (x,y)

}