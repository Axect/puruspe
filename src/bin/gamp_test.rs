extern crate puruspe;
use puruspe::gammp;

fn main() {
    for a in 1 .. 11 {
        for x in 0 .. 14 {
            let a_f64 = a as f64;
            let x_f64 = x as f64;
            println!("a:{}, x:{}, P:{}", a, x, gammp(a_f64, x_f64));
        }
    }
}
