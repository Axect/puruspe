use approx::assert_abs_diff_eq;
use puruspe::{erf, erfc, inverf, inverfc};

#[test]
fn test_erf() {
    for &(x, expected) in ERF_TABLE.iter() {
        let result = erf(x);
        let epsilon = f64::EPSILON + 1e-15 * result.abs();
        assert_abs_diff_eq!(result, expected, epsilon = epsilon);
    }
}

#[test]
fn test_erfc() {
    for &(x, expected) in ERFC_TABLE.iter() {
        let result = erfc(x);
        let epsilon = f64::EPSILON + 1e-14 * result.abs();
        assert_abs_diff_eq!(result, expected, epsilon = epsilon);
    }
}

#[test]
fn test_inverf() {
    for &(p, expected) in INVERF_TABLE.iter() {
        let result = inverf(p);
        let epsilon = f64::EPSILON + 1e-15 * result.abs();
        assert_abs_diff_eq!(result, expected, epsilon = epsilon);
    }
}

#[test]
fn test_inverfc() {
    for &(p, expected) in INVERFC_TABLE.iter() {
        let result = inverfc(p);
        let epsilon = f64::EPSILON + 1e-15 * result.abs();
        assert_abs_diff_eq!(result, expected, epsilon = epsilon);
    }
}

// ERF_TABLE, ERFC_TABLE, INVERF_TABLE, INVERFC_TABLE은 여기에 추가될 예정입니다.

const ERF_TABLE: [(f64, f64); 11] = [
    (-3.00000000000000e+00, -9.99977909503001e-01),
    (-1.50000000000000e+00, -9.66105146475311e-01),
    (-5.00000000000000e-01, -5.20499877813047e-01),
    (0.00000000000000e+00, 0.00000000000000e+00),
    (5.00000000000000e-01, 5.20499877813047e-01),
    (1.50000000000000e+00, 9.66105146475311e-01),
    (3.00000000000000e+00, 9.99977909503001e-01),
    (-1.00000000000000e-10, -1.12837916709551e-10),
    (1.00000000000000e-10, 1.12837916709551e-10),
    (-1.00000000000000e+01, -1.00000000000000e+00),
    (1.00000000000000e+01, 1.00000000000000e+00),
];

const ERFC_TABLE: [(f64, f64); 11] = [
    (-3.00000000000000e+00, 1.99997790950300e+00),
    (-1.50000000000000e+00, 1.96610514647531e+00),
    (-5.00000000000000e-01, 1.52049987781305e+00),
    (0.00000000000000e+00, 1.00000000000000e+00),
    (5.00000000000000e-01, 4.79500122186953e-01),
    (1.50000000000000e+00, 3.38948535246893e-02),
    (3.00000000000000e+00, 2.20904969985854e-05),
    (-1.00000000000000e-10, 1.00000000011284e+00),
    (1.00000000000000e-10, 9.99999999887162e-01),
    (-1.00000000000000e+01, 2.00000000000000e+00),
    (1.00000000000000e+01, 2.08848758376254e-45),
];

const INVERF_TABLE: [(f64, f64); 6] = [
    (-9.90000000000000e-01, -1.82138636771845e+00),
    (-5.00000000000000e-01, -4.76936276204470e-01),
    (0.00000000000000e+00, 0.00000000000000e+00),
    (5.00000000000000e-01, 4.76936276204470e-01),
    (9.90000000000000e-01, 1.82138636771845e+00),
    (1.00000000000000e-10, 8.86226925452758e-11),
];

const INVERFC_TABLE: [(f64, f64); 5] = [
    (5.00000000000000e-01, 4.76936276204470e-01),
    (9.90000000000000e-01, 8.86250128095061e-03),
    (1.00000000000000e-10, 4.57282496738949e+00),
    (1.00000000000000e+00, -0.00000000000000e+00),
    (1.99000000000000e+00, -1.82138636771845e+00),
];
