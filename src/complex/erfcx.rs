use super::erfcx_chebyshev_coeffs::*;
use crate::utils::frexp;

fn cheb_interpolant(x: f64) -> f64 {
    // Application-specific constants:
    const M: usize = 6; // 2^M subranges
    let j0 = -2; // first octave runs from 2^(j0-1) to 2^j0
    let l0 = 0; // index of x_min in full first octave
    let loff = (j0 + 1) * (1 << M) + l0; // precomputed offset

    // For given x, obtain mantissa xm and exponent je:
    let (xm, je) = frexp(x); // sets xm and je

    // Integer arithmetics to obtain reduced coordinate t:
    let ip = (1 << (M + 1)) * xm as usize; // index in octave + 2^M
    let lij = (je * (1 << M) + ip as i32 - loff) as usize; // index in lookup table
    let t = ((1 << (M + 2)) * xm as usize - (1 + 2 * ip)) as usize;

    let p_idx = lij * 8;
    let q_idx = lij * 2;

    // hard-coded for N1cheb=10
    ((((((((CHEBYSHEV_COEFFS_0[p_idx + 7] * t as f64 + CHEBYSHEV_COEFFS_0[p_idx + 6])
        * t as f64
        + CHEBYSHEV_COEFFS_0[p_idx + 5])
        * t as f64
        + CHEBYSHEV_COEFFS_0[p_idx + 4])
        * t as f64
        + CHEBYSHEV_COEFFS_0[p_idx + 3])
        * t as f64
        + CHEBYSHEV_COEFFS_0[p_idx + 2])
        * t as f64
        + CHEBYSHEV_COEFFS_0[p_idx + 1])
        * t as f64
        + CHEBYSHEV_COEFFS_0[p_idx])
        * t as f64
        + CHEBYSHEV_COEFFS_1[q_idx + 1])
        * t as f64
        + CHEBYSHEV_COEFFS_1[q_idx]
}

/******************************************************************************/
/*  Library function erfcx                                                    */
/******************************************************************************/

pub fn erfcx(x: f64) -> f64 {
    // Steven G. Johnson, October 2012.
    // Rewritten for better accuracy by Joachim Wuttke, Sept 2024.

    // Uses the following methods:
    // - Asymptotic expansion for large positive x,
    // - Chebyshev polynomials for medium positive x,
    // - Taylor (Maclaurin) series for small |x|,
    // - Asymptote 2exp(x^2) for large negative x,
    // - 2*exp(x^2)-erfcx(-x) for medium negative x.

    let ax = x.abs();

    if ax < 0.125 {
        // Use Taylor expansion
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

    /* else */
    {
        // Use asymptotic expansion
        //
        // Coefficient are a_0 = 1/sqrt(pi), a_N = (2N-1)!!/2^N/sqrt(pi).

        let r = 1. / x;

        if x < 150. {
            if x < 23.2 {
                return (((((((((((3.6073371500083758e+05) * (r * r)
                    - 3.7971970000088164e+04)
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
            return (((((((9.1636730015295726e+01) * (r * r) - 1.6661223639144676e+01)
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
} // erfcx
