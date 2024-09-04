//! This module provides implementations of the error function and related functions.
//! 
//! It includes the following main functions:
//! - `erf`: Calculates the error function.
//! - `erfc`: Calculates the complementary error function.
//! - `inverf`: Calculates the inverse of the error function.
//! - `inverfc`: Calculates the inverse of the complementary error function.

/// Calculates the error function.
///
/// The error function is defined as:
///
/// $$ \text{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt $$
///
/// # Arguments
///
/// * `x` - The input value
///
/// # Returns
///
/// The value of the error function at `x`
pub fn erf(x: f64) -> f64 {
    if x >= 0f64 {
        1.0 - erfccheb(x)
    } else {
        erfccheb(-x) - 1f64
    }
}

/// Calculates the complementary error function.
///
/// The complementary error function is defined as:
///
/// $$ \text{erfc}(x) = 1 - \text{erf}(x) = \frac{2}{\sqrt{\pi}} \int_x^\infty e^{-t^2} dt $$
///
/// # Arguments
///
/// * `x` - The input value
///
/// # Returns
///
/// The value of the complementary error function at `x`
pub fn erfc(x: f64) -> f64 {
    if x >= 0f64 {
        erfccheb(x)
    } else {
        2f64 - erfccheb(-x)
    }
}

/// Calculates the inverse of the error function.
///
/// This function finds x such that:
///
/// $$ p = \text{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt $$
///
/// # Arguments
///
/// * `p` - The probability value (between -1 and 1)
///
/// # Returns
///
/// The value of $x$ for which $\text{erf}(x) = p$
pub fn inverfc(p: f64) -> f64 {
    // Return arbitrary large pos or neg value
    if p >= 2f64 {
        return -100f64;
    } else if p <= 0f64 {
        return 100f64;
    }

    let pp = if p < 1f64 { p } else { 2f64 - p };
    let t = (-2f64 * (pp / 2f64).ln()).sqrt();
    let mut x = -std::f64::consts::FRAC_1_SQRT_2
        * ((2.30753 + t * 0.27061) / (1f64 + t * (0.99229 + t * 0.04481)) - t);
    for _j in 0..2 {
        let err = erfc(x) - pp;
        x += err / (std::f64::consts::FRAC_2_SQRT_PI * (-x.powi(2)).exp() - x * err);
    }
    if p < 1f64 {
        x
    } else {
        -x
    }
}

/// Calculates the inverse of the complementary error function.
///
/// This function finds x such that:
///
/// $$ p = \text{erfc}(x) = 1 - \text{erf}(x) = \frac{2}{\sqrt{\pi}} \int_x^\infty e^{-t^2} dt $$
///
/// # Arguments
///
/// * `p` - The probability value (between 0 and 2)
///
/// # Returns
///
/// The value of $x$ for which $\text{erfc}(x) = p$
pub fn inverf(p: f64) -> f64 {
    inverfc(1f64 - p)
}

// =============================================================================
// Chebyshev coefficients
// =============================================================================
const NCOEF: usize = 28;
const COF: [f64; 28] = [
    -1.3026537197817094,
    6.4196979235649026e-1,
    1.9476473204185836e-2,
    -9.561514786808631e-3,
    -9.46595344482036e-4,
    3.66839497852761e-4,
    4.2523324806907e-5,
    -2.0278578112534e-5,
    -1.624290004647e-6,
    1.303655835580e-6,
    1.5626441722e-8,
    -8.5238095915e-8,
    6.529054439e-9,
    5.059343495e-9,
    -9.91364156e-10,
    -2.27365122e-10,
    9.6467911e-11,
    2.394038e-12,
    -6.886027e-12,
    8.94487e-13,
    3.13092e-13,
    -1.12708e-13,
    3.81e-16,
    7.106e-15,
    -1.523e-15,
    -9.4e-17,
    1.21e-16,
    -2.8e-17,
];

/// Chebyshev coefficients
fn erfccheb(z: f64) -> f64 {
    let mut d = 0f64;
    let mut dd = 0f64;

    assert!(z >= 0f64, "erfccheb requires nonnegative argument");
    let t = 2f64 / (2f64 + z);
    let ty = 4f64 * t - 2f64;
    for j in (1..NCOEF - 1).rev() {
        let tmp = d;
        d = ty * d - dd + COF[j];
        dd = tmp;
    }
    t * (-z.powi(2) + 0.5 * (COF[0] + ty * d) - dd).exp()
}

