// =============================================================================
// Utils
// =============================================================================
/// Just factorial
pub fn factorial(n: usize) -> usize {
    let mut p = 1usize;
    for i in 1..(n + 1) {
        p *= i;
    }
    p
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
