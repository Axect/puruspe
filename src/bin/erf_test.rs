extern crate puruspe;
use puruspe::*;

fn main() {
    for x in 0..10 {
        let erf1 = erf(x as f64);
        let erfc1 = erfc(x as f64);
        println!(
            "X: {}, ERF: {}, ERFC: {}, INVERF:{}, INVERFC: {}",
            x,
            erf1,
            erfc1,
            inverf(erf1),
            inverfc(erfc1)
        );
    }
}
