extern crate puruspe;
use puruspe::*;

fn main() {
    for a in 1 .. 11 {
        for x in 0 .. 14 {
            let a_f64 = a as f64;
            let x_f64 = x as f64;
            let p = gammp(a_f64, x_f64);
            println!("a:{}, x:{}, P:{}", a, x, p);
            println!("a:{}, P:{}, x:{}", a, p, invgammp(p, a_f64));
        }
    }
}
