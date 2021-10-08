fn main(){
    cc::Build::new()
       // .file("src/dierckx.c")
        .file("src/dierckx-src/curfit.f")
        .file("src/dierckx-src/splev.f")
        .file("src/dierckx-src/fpback.f")
        .file("src/dierckx-src/fpbspl.f")
        .file("src/dierckx-src/fpchec.f")
        .file("src/dierckx-src/fpcurf.f")
        .file("src/dierckx-src/fpdisc.f")
        .file("src/dierckx-src/fpgivs.f")
        .file("src/dierckx-src/fpknot.f")
        .file("src/dierckx-src/fprati.f")
        .file("src/dierckx-src/fprota.f")
        .compiler("gfortran")
        .flag("-w") // no warnings
        .flag("-O3") // opitmize level 3
        .flag("-fdefault-real-8") // use f64
        .compile("dierckx")
}

