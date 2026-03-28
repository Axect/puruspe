// =============================================================================
// Import from other crates
// =============================================================================
// Lambert W functions
pub use lambert_w::{lambert_w0, lambert_wm1, sp_lambert_w0, sp_lambert_wm1};
use lambert_w::lambert_w as complex_lambert_w;

// Re-export the complex Lambert W function as the version of `lambert_w` from version 2 of the crate takes an error tolerance,
// which we set to floating point epsilon to match the API and performance of version 1.
/// Branch k of the complex valued Lambert W function computed on 64-bit floats with Halley's method.
/// The return value is a tuple where the first element is the real part and the second element is the imaginary part.
///
/// The function iterates until the current and previous iterations are within floating point epsilon, or it has iterated a maximum number of times.
///
/// This function may be slightly less accurate close to the branch cut at -1/e, as well as close to zero on branches other than k=0.
///
/// If you know you want the principal or secondary branches where they are real-valued, take a look at the [`lambert_w0`] or [`lambert_wm1`] functions instead. They can be up to two orders of magnitude faster.
// We inline this function since it's just a trivial wrapper.
#[inline]
pub fn lambert_w(k: i32, z_re: f64, z_im: f64) -> (f64, f64) {
    complex_lambert_w(k, z_re, z_im, f64::EPSILON)
}

// =============================================================================
// Export funcionts
// =============================================================================
pub mod utils;

pub mod gamma;
pub use gamma::{gamma, gammp, gammq, invgammp, ln_gamma};

pub mod beta;
pub use beta::{beta, betai, invbetai};

pub mod error;
pub use error::{erf, erfc, erfcx, inverf, inverfc};

pub mod dawson;
pub use dawson::dawson;

pub mod faddeeva;
pub use faddeeva::faddeeva;

pub mod bessel;
pub use bessel::{
    besselik, besseljy, CachedBesselIK, CachedBesselJY, CachedInuKnu, CachedJnuYnu, In, Inu_Knu,
    Jn, Jnu_Ynu, Kn, Yn,
};

// =============================================================================
// Constants
// =============================================================================
const EPS: f64 = f64::EPSILON;
const FPMIN: f64 = f64::MIN_POSITIVE / EPS;
const INV_SQRT_PI: f64 = 0.564_189_583_547_756_3;
// For gamma functions & beta functions
const Y: [f64; 18] = [
    0.0021695375159141994,
    0.011413521097787704,
    0.027972308950302116,
    0.051_727_015_600_492_42,
    0.082_502_225_484_340_94,
    0.12007019910960293,
    0.164_152_833_007_524_7,
    0.21442376986779355,
    0.27051082840644336,
    0.33199876341447887,
    0.39843234186401943,
    0.46931971407375483,
    0.544_136_055_566_579_7,
    0.622_327_452_880_310_8,
    0.703_315_004_655_971_7,
    0.786_499_107_683_134_5,
    0.871_263_896_190_615_2,
    0.956_981_801_526_291_4,
];
const W: [f64; 18] = [
    0.005_565_719_664_244_557,
    0.012_915_947_284_065_42,
    0.020181515297735382,
    0.027298621498568734,
    0.034_213_810_770_299_54,
    0.040_875_750_923_643_26,
    0.047_235_083_490_265_58,
    0.053_244_713_977_759_69,
    0.058_860_144_245_324_8,
    0.064_039_797_355_015_48,
    0.068_745_323_835_736_41,
    0.072_941_885_005_653_09,
    0.076_598_410_645_870_64,
    0.079_687_828_912_071_67,
    0.082_187_266_704_339_7,
    0.084_078_218_979_661_95,
    0.085_346_685_739_338_72,
    0.085_983_275_670_394_82,
];
