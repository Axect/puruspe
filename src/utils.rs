//! This module provides utility functions used throughout the library.
//!
//! It includes the following main functions:
//! - `factorial`: Calculates the factorial of a non-negative integer.
//! - `rising_factorial`: Computes the rising factorial (also known as Pochhammer function).
//! - `frexp`: Decomposes a floating-point number into its mantissa and exponent.
//! - `ldexp`: Multiplies a floating-point number by a power of two.
//! - `sign`: Returns a value with the magnitude of one argument and the sign of another.
//!
//! These utility functions are primarily used internally by other functions in the library,
//! but they may also be useful for users who need to perform similar calculations.

/// Precise implementation of factorial
pub fn factorial(n: usize) -> f64 {
    rising_factorial(1f64, n)
}

/// Rising factorial via rectangular splitting
///
/// $$ z(z+1) \cdots (z + n - 1) $$
fn rising_factorial(z: f64, n: usize) -> f64 {
    let m = (n as f64).sqrt().floor().min(60f64) as usize;

    let mut k = 0usize;
    let mut r = 1f64;

    while k < n {
        let l = m.min(n - k);
        let t = (k..k + l).fold(1f64, |acc, i| acc * (z + i as f64));
        r *= t;
        k += m;
    }
    r
}

/// Pure Rust implementation of frexp
pub fn frexp(x: f64) -> (f64, i32) {
    // If the input is zero, return (0.0, 0)
    if x == 0.0 {
        return (0.0, 0);
    }

    // Convert the input to its IEEE 754 binary representation
    let bits = x.to_bits();

    // Determine the sign of the input
    // If the most significant bit (bit 63) is 1, the input is negative
    let sign = if (bits >> 63) != 0 { -1.0 } else { 1.0 };

    // Extract the exponent from the binary representation
    // Bits 52 to 62 represent the exponent in IEEE 754 format
    // Subtract 1023 to obtain the actual exponent value
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023;

    // Extract the mantissa (significand) from the binary representation
    // Bits 0 to 51 represent the mantissa
    // Set the implicit leading bit (bit 52) to 1 to normalize the mantissa
    // Multiply the mantissa by the sign to handle negative inputs correctly
    let mantissa = sign * f64::from_bits((bits & 0xfffffffffffff) | 0x3fe0000000000000);

    // Return the normalized mantissa and the exponent incremented by 1
    (mantissa, exponent + 1)
}

/// Pure Rust implementation of ldexp
pub fn ldexp(x: f64, exp: i32) -> f64 {
    // If the input is zero or the exponent is zero, return the input unchanged
    if x == 0.0 || exp == 0 {
        return x;
    }

    // Convert the input to its IEEE 754 binary representation
    let bits = x.to_bits();

    // Extract the exponent from the binary representation
    // Bits 52 to 62 represent the exponent in IEEE 754 format
    let exponent = ((bits >> 52) & 0x7ff) as i32;

    // Calculate the new exponent by adding the input exponent to the existing exponent
    let new_exponent = exponent + exp;

    // Check if the new exponent is within the valid range for IEEE 754 format
    if !(0..=0x7ff).contains(&new_exponent) {
        // If the exponent is out of range, return infinity or zero depending on the input sign
        return if (bits >> 63) != 0 {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }

    // Combine the new exponent with the mantissa and sign bits to create the result
    let result_bits = (bits & 0x800fffffffffffff) | ((new_exponent as u64) << 52);
    f64::from_bits(result_bits)
}

/// The sign function is a helper that returns
/// a value with the magnitude of a and sign of b.
pub fn sign(a: f64, b: f64) -> f64 {
    if b.is_sign_positive() {
        a.abs()
    } else {
        -(a.abs())
    }
}
