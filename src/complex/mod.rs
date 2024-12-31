//! This module is gated behind the `cmplx` feature flag.
//! Enable this feature to use the functions provided here.
//!
//! The Faddeeva function is a complex error function which
//! can be used to compute complex values of the more commonly
//! used error functions, such as `erfc(z)`. Consequently, this
//! module also provides complex-valued implementations of these
//! error functions.

mod erfcx;
mod erfcx_chebyshev_coeffs;
mod faddeeva;
mod faddeeva_taylor_coeffs;

const INV_SQRT_PI: f64 = 0.5641895835_4775628694; // Inverse square root of pi from WolframAlpha.

pub use erfcx::erfcx;
pub use faddeeva::w;
