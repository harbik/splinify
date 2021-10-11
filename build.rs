fn main(){
    cc::Build::new()
       // .file("src/dierckx.c")
        .file("src/dierckx-src/concur.f")
        .file("src/dierckx-src/curev.f")
        .file("src/dierckx-src/curfit.f")
        .file("src/dierckx-src/fpadpo.f")
        .file("src/dierckx-src/fpback.f")
        .file("src/dierckx-src/fpbspl.f")
        .file("src/dierckx-src/fpchec.f")
        .file("src/dierckx-src/fpched.f")
        .file("src/dierckx-src/fpcons.f")
        .file("src/dierckx-src/fpcurf.f")
        .file("src/dierckx-src/fpdisc.f")
        .file("src/dierckx-src/fpgivs.f")
        .file("src/dierckx-src/fpinst.f")
        .file("src/dierckx-src/fpknot.f")
        .file("src/dierckx-src/fppocu.f")
        .file("src/dierckx-src/fprati.f")
        .file("src/dierckx-src/fprota.f")
        .file("src/dierckx-src/splev.f")
        .compiler("gfortran")
        .flag("-w") // no warnings
        .flag("-O3") // opitmize level 3
        .flag("-fdefault-real-8") // use f64
        .compile("dierckx")
}

