//! This module provides implementations of the gamma function and related functions.
//!
//! It includes the following main functions:
//! - `ln_gamma`: Calculates the natural logarithm of the gamma function.
//! - `gamma`: Calculates the gamma function.
//! - `gammp`: Calculates the regularized lower incomplete gamma function.
//! - `gammq`: Calculates the regularized upper incomplete gamma function.
//! - `invgammp`: Calculates the inverse of the regularized lower incomplete gamma function.

use crate::utils::polynomial;
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
// Based on the Fortran implementation by Toshio Fukushima: <https://www.researchgate.net/publication/336578125_xgamtxt_Fortran_90_test_program_package_of_qgam_dgam_and_sgam_fast_computation_of_Gamma_function_with_quadruple_double_and_single_precision_accuracy_respectively>
pub fn gamma(mut x: f64) -> f64 {
    /// This input gives an output of [`f64::MAX`].
    /// Computed with WolframAlpha: https://www.wolframalpha.com/input?i=solve+gamma%28x%29+%3D+%281+%E2%88%92+2%5E%28%E2%88%9253%29%29*2%5E1024+for+x+%3E+1.
    const MAX_INPUT: f64 = 171.624_376_956_302_7;

    // Special cases.
    if x > MAX_INPUT {
        // Output is too large to represent.
        return f64::INFINITY;
    } else if x == 0.0 {
        // The gamma function diverges for an input of zero.
        // It diverges to positive or negative infinity depending on which direction zero is approached from.
        // According to entry F.9.5.4 of the ISO C 99 standard (https://www.open-std.org/jtc1/sc22/wg14/www/docs/n1256.pdf)
        // gamma(+/- 0) = +/- infinity.
        // This is the standard followed by `scipy.special.gamma` as well as the `tgamma` function in `libc`.
        // Since the Gamma function in (nightly) Rust's standard library currently corresponds to `tgamma`
        // and thus follows the standard we do the same here.
        return f64::INFINITY.copysign(x);
    } else if x < 0.0 && x.fract() == 0.0 {
        // The Gamma function diverges for non-positive integers.
        // There is however no clear way for the caller to signal whether the input approached the integer from above or below.
        // According to the same standard as above the gamma function should return NaN for these inputs.
        return f64::NAN;
    }

    let f = if x > 3.5 {
        let mut f = 1.0;
        // Decrement x by 1 until less than 3.5.
        // This will not be a massively long loop since we checked for too large input above.
        while x >= 3.5 {
            x -= 1.0;
            f *= x;
        }
        f
    } else if x < 2.5 {
        let mut f = 1.0;
        // Increment x by 1 until larger than 2.5.
        while x <= 2.5 {
            f *= x;
            if f.is_infinite() {
                // There is no point in continuing the calculation,
                // as `f` will remain infinite when multiplied by any finite number.
                // (except zero, but that can only happen if x is an integer, which we have handled above.)
                // The final result will thus be `g` / inf,
                // which underflows to 0 since `g` (defined below) is always finite.
                // This early return stops the loop from running too long for very large negative inputs.
                return 0.0;
            }
            x += 1.0;
        }
        f.recip()
    } else {
        1.0
    };

    let g = if x > 3.0 {
        polynomial(
            x - 3.25,
            &[
                2.549_256_966_718_529,
                2.592_571_165_129_980_8,
                1.776_919_804_709_978_3,
                0.858_853_822_880_646_4,
                0.346_268_553_360_568_7,
                0.115_260_482_799_219_55,
                0.034_350_198_197_417_38,
                0.008_882_603_286_354_344,
                0.002_149_849_572_114_427,
                0.000_455_766_746_785_133_16,
                0.000_095_940_184_093_793_52,
                0.000_016_738_398_919_923_317,
            ],
        )
    } else if x < 3.0 {
        polynomial(
            x - 2.75,
            &[
                1.608_359_421_985_545_5,
                1.317_087_179_192_858_7,
                0.891_167_948_026_476_6,
                0.384_779_120_148_185_27,
                0.155_955_856_853_634_36,
                0.046_735_160_424_814_17,
                0.014_576_803_693_379_53,
                0.003_228_437_222_320_849,
                0.000_914_252_863_201_387_9,
                0.000_132_763_577_755_739_43,
                0.000_048_420_884_483_658_204,
                7.377_657_774_398_435e-7,
            ],
        )
    } else {
        2.0
    };

    g * f
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

    // Handle edge cases
    if x == 0f64 {
        return 0f64;
    }

    // For very small x relative to a, use series (always more accurate)
    if x < a * EPS {
        return 0f64;
    }

    // For very large x relative to a, result approaches 1
    if x > a + 20.0 * (a.sqrt() + 1.0) {
        return 1f64;
    }

    if (a as usize) >= ASWITCH {
        // For very small x relative to a, use series instead of quadrature
        if x < 0.2 * a {
            gser(a, x)
        } else {
            // Quadrature for large a
            gammpapprox(a, x, IncGamma::P)
        }
    } else if x < a + 1f64 {
        // Series representation - more accurate for x < a+1
        gser(a, x)
    } else {
        // Continued fraction representation
        // Use complementary function to avoid catastrophic cancellation
        let result = 1f64 - gcf(a, x);
        result.clamp(0.0, 1.0)
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
    assert!(x >= 0f64 && a > 0f64, "Bad args in gammq");

    // Handle edge cases
    if x == 0f64 {
        return 1f64;
    }

    // For very small x relative to a, result approaches 1
    if x < a * EPS {
        return 1f64;
    }

    // For very large x relative to a, result approaches 0
    if x > a + 20.0 * (a.sqrt() + 1.0) {
        return 0f64;
    }

    if (a as usize) >= ASWITCH {
        // For very small x relative to a, use series instead of quadrature
        if x < 0.2 * a {
            let result = 1f64 - gser(a, x);
            result.clamp(0.0, 1.0)
        } else {
            // Quadrature for large a
            gammpapprox(a, x, IncGamma::Q)
        }
    } else if x < a + 1f64 {
        // Series representation
        // Use complementary function to avoid catastrophic cancellation
        let result = 1f64 - gser(a, x);
        result.clamp(0.0, 1.0)
    } else {
        // Continued fraction representation - more accurate for x >= a+1
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
            // Compute in log space for better numerical stability
            let log_result = -x + a * x.ln() - gln;

            // Guard against overflow/underflow
            if log_result > 700.0 {
                return sum * f64::INFINITY;
            } else if log_result < -700.0 {
                return 0.0;
            }

            return sum * log_result.exp();
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

    // Compute in log space for better numerical stability
    let log_result = -x + a * x.ln() - gln;

    // Guard against overflow/underflow
    if log_result > 700.0 {
        return h * f64::INFINITY;
    } else if log_result < -700.0 {
        return 0.0;
    }

    log_result.exp() * h
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

    // Use log-sum-exp trick for better numerical stability
    let mut log_values = Vec::with_capacity(NGAU);
    let mut t: f64;
    for j in 0..NGAU {
        t = x + (xu - x) * Y[j];
        let log_term = W[j].ln() + (-(t - a1) + a1 * (t.ln() - lna1));
        log_values.push(log_term);
    }

    // Find maximum for log-sum-exp
    let log_max = log_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Compute sum using log-sum-exp trick
    let mut sum = 0f64;
    for log_val in log_values {
        if log_val - log_max > -700.0 {
            sum += (log_val - log_max).exp();
        }
    }

    // Compute final answer with overflow/underflow protection
    let log_scale = a1 * (lna1 - 1f64) - gln;
    let log_ans = log_max + sum.ln() + (xu - x).ln() + log_scale;

    let ans = if log_ans > 700.0 {
        f64::INFINITY
    } else if log_ans < -700.0 {
        0.0
    } else {
        log_ans.exp()
    };

    // Clamp results to [0, 1] and handle the P vs Q case
    match psig {
        IncGamma::P => {
            let result = 1f64 - ans;
            result.clamp(0.0, 1.0)
        }
        IncGamma::Q => ans.clamp(0.0, 1.0),
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
