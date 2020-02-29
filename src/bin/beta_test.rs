extern crate puruspe;
use puruspe::*;

fn main() {
    for a in 1 .. 6 {
        for b in 1 .. 6 {
            for x in 0 .. 11 {
                let a_f64 = a as f64;
                let b_f64 = b as f64;
                let x_f64 = x as f64 / 10f64;
                let i = betai(a_f64, b_f64, x_f64);
                println!("a:{}, b: {}, x:{}, I:{}", a, b, x_f64, i);
                println!("a:{}, b: {}, I:{}, x:{}", a, b, i, invbetai(i, a_f64, b_f64));
            }
        }
    }
}
