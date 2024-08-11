use puruspe::{gammp, invgammp};

#[test]
fn gammp_test() {
    let expected_1 = 0.05265301734; // ɣ(5,2)/Γ(5) according to Wolfram|Alpha.
    let a = 5_f64;
    let x = 2_f64;
    let res_1 = gammp(a, x);
    assert!((res_1 - expected_1).abs() < 1e-11);

    let expected_2 = 0.96739453751; // Ditto.
    let a = 0.8;
    let x = 3_f64;
    let res_2 = gammp(a, x);
    assert!((res_2 - expected_2).abs() < 1e-11);
}

#[test]
fn invgammp_test() {
    // Invert the first test above.
    let expected_1 = 2_f64;
    let a = 5_f64;
    let p = 0.05265301734371115;
    let res_1 = invgammp(p, a);
    assert!((res_1 - expected_1).abs() < 1e-10);

    // Invert the second test.
    let expected_2 = 3_f64;
    let a = 0.8;
    let p = 0.96739453751512363;
    let res_2 = invgammp(p, a);
    assert!((res_2 - expected_2).abs() < 1e-10);
}
