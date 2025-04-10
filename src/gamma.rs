//! This module provides implementations of the gamma function and related functions.
//!
//! It includes the following main functions:
//! - `ln_gamma`: Calculates the natural logarithm of the gamma function.
//! - `gamma`: Calculates the gamma function.
//! - `gammp`: Calculates the regularized lower incomplete gamma function.
//! - `gammq`: Calculates the regularized upper incomplete gamma function.
//! - `invgammp`: Calculates the inverse of the regularized lower incomplete gamma function.

use crate::utils::factorial;
use crate::{EPS, FPMIN, W, Y};
use core::f64::consts::PI;
const ASWITCH: usize = 100;
const NGAU: usize = 18;
const G: f64 = 5f64;
const N: usize = 7;

/// Calculates the natural logarithm of the gamma function.
///
/// The gamma function is defined as:
///
/// $$ \Gamma(z) = \int_0^\infty t^{z-1} e^{-t} dt $$
///
/// This function computes $\ln(\Gamma(z))$.
///
/// # Arguments
///
/// * `z` - The input value
///
/// # Returns
///
/// The natural logarithm of the gamma function at `z`
pub fn ln_gamma(z: f64) -> f64 {
    let z = z - 1f64;
    let base = z + G + 0.5;
    let mut s = 0f64;
    for (i, &coeff) in LG5N7.iter().enumerate().skip(1) {
        s += coeff / (z + i as f64);
    }
    s += LG5N7[0];
    (2f64 * PI).sqrt().ln() + s.ln() - base + base.ln() * (z + 0.5)
}

/// Calculates the gamma function.
///
/// The gamma function is defined as:
///
/// $$ \Gamma(z) = \int_0^\infty t^{z-1} e^{-t} dt $$
///
/// # Arguments
///
/// * `z` - The input value
///
/// # Returns
///
/// The value of the gamma function at `z`
pub fn gamma(z: f64) -> f64 {
    if z >= 1f64 {
        let z_int = z as usize;
        if (z - (z_int as f64)).abs() < EPS {
            return factorial(z_int - 1);
        }
    }

    if z < 0.5 {
        PI / ((PI * z).sin() * gamma(1f64 - z))
    } else {
        ln_gamma(z).exp()
    }
}

// =============================================================================
// Incomplete Gamma function (regularized)
// =============================================================================
/// Calculates the regularized lower incomplete gamma function P(a,x).
///
/// The regularized lower incomplete gamma function is defined as:
///
/// $$ P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)} = \frac{1}{\Gamma(a)} \int_0^x t^{a-1} e^{-t} dt $$
///
/// # Arguments
///
/// * `a` - The shape parameter
/// * `x` - The upper limit of integration
///
/// # Returns
///
/// The value of P(a,x)
///
/// # Panics
///
/// Panics if `x` < 0 or if `a` ≤ 0.
pub fn gammp(a: f64, x: f64) -> f64 {
    assert!(x >= 0f64 && a > 0f64, "Bad args in gammp");
    if x == 0f64 {
        0f64
    } else if (a as usize) >= ASWITCH {
        // Quadrature
        gammpapprox(a, x, IncGamma::P)
    } else if x < a + 1f64 {
        // Series representation
        gser(a, x)
    } else {
        // Continued fraction representation
        1f64 - gcf(a, x)
    }
}

/// Calculates the regularized upper incomplete gamma function Q(a,x).
///
/// The regularized upper incomplete gamma function is defined as:
///
/// $$ Q(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)} = \frac{1}{\Gamma(a)} \int_x^\infty t^{a-1} e^{-t} dt $$
///
/// Note that $Q(a,x) = 1 - P(a,x)$.
///
/// # Arguments
///
/// * `a` - The shape parameter
/// * `x` - The lower limit of integration
///
/// # Returns
///
/// The value of Q(a,x)
///
/// # Panics
///
/// Panics if `x` < 0 or if `a` ≤ 0
pub fn gammq(a: f64, x: f64) -> f64 {
    assert!(x >= 0f64 && a > 0f64, "Bad args in gammp");
    if x == 0f64 {
        1f64
    } else if (a as usize) >= ASWITCH {
        // Quadrature
        gammpapprox(a, x, IncGamma::Q)
    } else if x < a + 1f64 {
        // Series representation
        1f64 - gser(a, x)
    } else {
        // Continued fraction representation
        gcf(a, x)
    }
}

/// Series expansion
fn gser(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut del = 1f64 / a;
    let mut sum = 1f64 / a;
    loop {
        ap += 1f64;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * EPS {
            return sum * (-x + a * x.ln() - gln).exp();
        }
    }
}

/// Continued Fraction
fn gcf(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut b = x + 1f64 - a;
    let mut c = 1f64 / FPMIN;
    let mut d = 1f64 / b;
    let mut h = d;
    let mut an: f64;
    for i in 1.. {
        an = -i as f64 * (i as f64 - a);
        b += 2f64;
        d = an * d + b;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = b + an / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1f64 / d;
        let del = d * c;
        h *= del;
        if (del - 1f64).abs() < EPS {
            break;
        }
    }
    (-x + a * x.ln() - gln).exp() * h
}

/// Kinds of Incomplete Gamma function
#[derive(Debug, Copy, Clone)]
enum IncGamma {
    P,
    Q,
}

/// Gauss Legendre Quadrature (order of 18)
fn gammpapprox(a: f64, x: f64, psig: IncGamma) -> f64 {
    let a1 = a - 1f64;
    let lna1 = a1.ln();
    let sqrta1 = a1.sqrt();
    let gln = ln_gamma(a);
    let xu = if x > a1 {
        (a1 + 11.5 * sqrta1).max(x + 6f64 * sqrta1)
    } else {
        0f64.max((a1 - 7.5 * sqrta1).min(x - 5f64 * sqrta1))
    };
    let mut sum = 0f64;
    let mut t: f64;
    for j in 0..NGAU {
        t = x + (xu - x) * Y[j];
        sum += W[j] * (-(t - a1) + a1 * (t.ln() - lna1)).exp();
    }
    let ans = sum * (xu - x) * (a1 * (lna1 - 1f64) - gln).exp();
    match psig {
        IncGamma::P => {
            if ans > 0f64 {
                1f64 - ans
            } else {
                -ans
            }
        }
        IncGamma::Q => {
            if ans >= 0f64 {
                ans
            } else {
                1f64 + ans
            }
        }
    }
}

/// Calculates the inverse of the regularized lower incomplete gamma function.
///
/// This function finds x such that:
///
/// $$ p = P(a,x) = \frac{1}{\Gamma(a)} \int_0^x t^{a-1} e^{-t} dt $$
///
/// # Arguments
///
/// * `p` - The probability value (between 0 and 1)
/// * `a` - The shape parameter
///
/// # Returns
///
/// The value of x for which P(a,x) = p
///
/// # Panics
///
/// Panics if `a` ≤ 0.
pub fn invgammp(p: f64, a: f64) -> f64 {
    let gln = ln_gamma(a);
    let a1 = a - 1f64;
    let lna1 = a1.ln();
    let mut afac = 0f64;
    let pp: f64;
    let mut t: f64;

    assert!(a > 0f64, "a must be positive in invgammp");
    if p >= 1f64 {
        return 100f64.max(a + 100f64 * a.sqrt());
    } else if p <= 0f64 {
        return 0f64;
    }

    // Initial guess
    let mut x = if a > 1f64 {
        afac = (a1 * (lna1 - 1f64) - gln).exp();
        pp = if p < 0.5 { p } else { 1f64 - p };
        t = (-2f64 * pp.ln()).sqrt();
        let mut x = (2.30753 + t * 0.27061) / (1f64 + t * (0.99229 + t * 0.04481)) - t;
        if p < 0.5 {
            x = -x;
        }
        1e-3_f64.max(a * (1f64 - 1f64 / (9f64 * a) - x / (3f64 * a.sqrt())).powi(3))
    } else {
        t = 1f64 - a * (0.253 + a * 0.12);
        if p < t {
            (p / t).powf(1f64 / a)
        } else {
            1f64 - (1f64 - (p - t) / (1f64 - t)).ln()
        }
    };

    for _j in 0..12 {
        // x is too small to compute accurately
        if x <= 0f64 {
            return 0f64;
        }
        let err = gammp(a, x) - p;
        t = if a > 1f64 {
            afac * (-(x - a1) + a1 * (x.ln() - lna1)).exp()
        } else {
            (-x + a1 * x.ln() - gln).exp()
        };
        let u = err / t;
        // Halley's method
        t = u / (1f64 - 0.5 * 1f64.min(u * (a1 / x - 1f64)));
        x -= t;
        if x <= 0f64 {
            x = 0.5 * (x + t);
        }
        if t.abs() < (x * EPS).max(EPS) {
            break;
        }
    }
    x
}

// =============================================================================
// Lanczos approximation of Gamma
// =============================================================================
/// Lanczos g=5, n=7
const LG5N7: [f64; N] = [
    1.000000000189712,
    76.18009172948503,
    -86.50532032927205,
    24.01409824118972,
    -1.2317395783752254,
    0.0012086577526594748,
    -0.00000539702438713199,
];
