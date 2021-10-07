fn main(){
    cc::Build::new()
        .file("src/dierckx.c")
        .compile("dierckx")
}