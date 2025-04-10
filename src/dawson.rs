//! This module provides an implementation of Dawson's integral.
//!
//! It includes the following main function:
//! - `dawson`: Calculates Dawson's integral for real values.

use crate::utils::sign;
use core::f64::consts::PI;

/// Dawson's integral for real values.
/// Accurate to about 2e-7 for small inputs.
pub fn dawson(x: f64) -> f64 {
    const NMAX: usize = 6;
    let h = 0.4;
    let a1 = 2. / 3.;
    let a2 = 0.4;
    let a3 = 2. / 7.;
    let n0: u64; // i is declared in the loops.
    let mut d1: f64;
    let mut d2: f64;
    let mut e1: f64;
    let e2: f64;
    let mut sum: f64;
    let x2: f64;
    let xp: f64;
    let xx: f64;
    let ans: f64;

    if x.abs() < 0.2 {
        x2 = x * x;
        ans = x * (1.0 - a1 * x2 * (1.0 - a2 * x2 * (1.0 - a3 * x2)));
    } else {
        xx = x.abs();
        n0 = 2 * (0.5 * xx / h + 0.5) as u64; // Simply truncates the decimal places, as in NR.
        xp = xx - n0 as f64 * h;
        e1 = f64::exp(2.0 * xp * h);
        e2 = e1 * e1;
        d1 = n0 as f64 + 1.;
        d2 = d1 - 2.0;
        sum = 0.0;

        for i in 0..NMAX {
            let cval = f64::exp(-((2. * i as f64 + 1.) * h).powi(2));
            d1 += 2.;
            d2 -= 2.;
            e1 *= e2;

            sum += cval * (e1 / d1 + 1.0 / (d2 * e1));
        }

        ans = (1. / PI.sqrt()) * sign(f64::exp(-xp * xp), x) * sum;
    }
    ans
}
