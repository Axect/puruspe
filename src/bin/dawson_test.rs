extern crate puruspe;
use puruspe::dawson;

fn main() {
    let expected = 0.13818492867352312; // From Wolfram|Alpha
    let res = dawson(0.14);
    dbg!(expected, res);
    assert!((expected - res).abs() < 1e-7);
}
