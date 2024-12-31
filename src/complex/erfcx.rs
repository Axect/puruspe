use super::erfcx_chebyshev_coeffs::*;
use crate::utils::frexp;

/// Helper function for erfcx(x).
fn cheb_interpolant(x: f64) -> f64 {
    // Application-specific constants:
    const M: i32 = 6; // 2^M subranges
    let j0 = -2; // first octave runs from 2^(j0-1) to 2^j0
    let l0 = 0; // index of x_min in full first octave
    let loff = (j0 + 1) * (1 << M) + l0; // precomputed offset

    // For given x, obtain mantissa xm and exponent je:
    let (xm, je) = frexp(x); // sets xm and je

    // Integer arithmetics to obtain reduced coordinate t:
    let ip = ((1 << (M + 1)) as f64 * xm) as i32; // index in octave + 2^M
    let lij = (je * (1 << M) + ip - loff) as i32; // index in lookup table
    let t = ((1 << (M + 2)) as f64 * xm) - (1 + 2 * ip) as f64;

    let p_idx = lij * 8;
    let q_idx = lij * 2;

    assert!(
        p_idx >= 0 && q_idx >= 0,
        "Negative indices not permitted in Rust."
    );

    // hard-coded for N1cheb=10
    ((((((((CHEBYSHEV_COEFFS_1[(p_idx + 7) as usize] * t
        + CHEBYSHEV_COEFFS_1[(p_idx + 6) as usize])
        * t
        + CHEBYSHEV_COEFFS_1[(p_idx + 5) as usize])
        * t
        + CHEBYSHEV_COEFFS_1[(p_idx + 4) as usize])
        * t
        + CHEBYSHEV_COEFFS_1[(p_idx + 3) as usize])
        * t
        + CHEBYSHEV_COEFFS_1[(p_idx + 2) as usize])
        * t
        + CHEBYSHEV_COEFFS_1[(p_idx + 1) as usize])
        * t
        + CHEBYSHEV_COEFFS_1[p_idx as usize])
        * t
        + CHEBYSHEV_COEFFS_0[(q_idx + 1) as usize])
        * t
        + CHEBYSHEV_COEFFS_0[q_idx as usize]
}

/// The scaled complementary error function, or erfcx(x).
///
/// ### Definition:
/// erfcx(x) = erfc(x) / exp(-x^2) = (1 - erf(x)) / exp(-x^2)
///
/// ### Input:
/// - x: A real number
///
/// ### Output:
/// The result of the computation, guaranteed to match `libcerf` to at least 1e-6.
///
/// ### Notes on Implementation:
/// This is a Rust port of the erfcx routine found in [`libcerf`](https://jugit.fz-juelich.de/mlz/libcerf),
/// which is under the MIT License.
///
/// It uses the following methods:
/// - Asymptotic expansion for large positive x,
/// - Chebyshev polynomials for medium positive x,
/// - Maclaurin series for small |x|,
/// - Asymptote 2exp(x^2) for large negative x,
/// - 2*exp(x^2)-erfcx(-x) for medium negative x.
pub fn erfcx(x: f64) -> f64 {
    let ax = x.abs();

    if ax < 0.125 {
        /* **************************************************************** */
        /* Use the Maclaurin series.                                        */
        /* **************************************************************** */
        return ((((((((((((((1.9841269841269841e-04) * x - 5.3440090793734269e-04) * x
            + 1.3888888888888889e-03)
            * x
            - 3.4736059015927274e-03)
            * x
            + 8.3333333333333332e-03)
            * x
            - 1.9104832458760001e-02)
            * x
            + 4.1666666666666664e-02)
            * x
            - 8.5971746064419999e-02)
            * x
            + 1.6666666666666666e-01)
            * x
            - 3.0090111122547003e-01)
            * x
            + 5.0000000000000000e-01)
            * x
            - 7.5225277806367508e-01)
            * x
            + 1.0000000000000000e+00)
            * x
            - 1.1283791670955126e+00)
            * x
            + 1.0000000000000000e+00;
    }

    if x < 0. {
        if x < -26.7 {
            return f64::INFINITY;
        }
        if x < -6.1 {
            return 2. * (x * x).exp();
        }
        return 2. * (x * x).exp() - cheb_interpolant(-x);
    }

    if x <= 12. {
        return cheb_interpolant(x);
    }

    /* **************************************************************** */
    /*     Use asymptotic expansion.                                    */
    /* Coefficients are a_0 = 1 / sqrt(PI),                             */
    /*                  a_N = (2N-1)!! / 2^N / sqrt(PI)                 */
    /* **************************************************************** */
    let r = 1. / x;

    if x < 150. {
        if x < 23.2 {
            return (((((((((((3.6073371500083758e+05) * (r * r) - 3.7971970000088164e+04)
                * (r * r)
                + 4.4672905882456671e+03)
                * (r * r)
                - 5.9563874509942218e+02)
                * (r * r)
                + 9.1636730015295726e+01)
                * (r * r)
                - 1.6661223639144676e+01)
                * (r * r)
                + 3.7024941420321507e+00)
                * (r * r)
                - 1.0578554691520430e+00)
                * (r * r)
                + 4.2314218766081724e-01)
                * (r * r)
                - 2.8209479177387814e-01)
                * (r * r)
                + 5.6418958354775628e-01)
                * r;
        }
        return (((((((9.1636730015295726e+01) * (r * r) - 1.6661223639144676e+01) * (r * r)
            + 3.7024941420321507e+00)
            * (r * r)
            - 1.0578554691520430e+00)
            * (r * r)
            + 4.2314218766081724e-01)
            * (r * r)
            - 2.8209479177387814e-01)
            * (r * r)
            + 5.6418958354775628e-01)
            * r;
    }
    if x < 6.9e7 {
        return ((((-1.0578554691520430e+00) * (r * r) + 4.2314218766081724e-01) * (r * r)
            - 2.8209479177387814e-01)
            * (r * r)
            + 5.6418958354775628e-01)
            * r;
    }
    // 1-term expansion, important to avoid overflow
    return 0.56418958354775629 / x;
}
