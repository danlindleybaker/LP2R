extern crate cc;

fn main() {
    cc::Build::new()
        .file("kww.c")
        .compile("kww.so");
}
