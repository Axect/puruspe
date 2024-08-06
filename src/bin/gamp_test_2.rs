extern crate puruspe;
use puruspe::*;

fn main() {
    let a = 100.;
    for x in [1., 10., 50., 75., 99., 100., 101., 110., 1000.] {
        println!("{}\t{}", x, gammp(a, x));
    }
}
