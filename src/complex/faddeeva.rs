use super::{faddeeva_taylor_coeffs::*, INV_SQRT_PI};

use crate::{complex::erfcx::erfcx, dawson};

/// The complex Faddeeva function, or w(z). This is Abramowitz
/// & Stegun function (7.1.3), also known as the plasma dispersion
/// function or Kramp function.
///
/// ### Definition:
/// w(z) = exp(-z^2) erfc(-iz)
///
/// ## Inputs:
/// - re: ℜ(z)
/// - im: ℑ(z)
///
/// ## Outputs:
/// A tuple of ([`f64`], [`f64`]) where the first
/// is ℜ(w(z)) and the second is ℑ(w(z)).
///
/// ### Notes on Implementation:
/// For sufficiently large |z|, we use a continued-fraction expansion,
/// which is an independent reimplementation of W. Gautschi's algorithm
/// in [this paper](http://dx.doi.org/10.1137/0707012). It is roughly
/// similar to [Algorithm 680](http://dx.doi.org/10.1145/77626.77629).
///
/// For ℑ(z) near the real axis or small |z|, we use an independent reimplementation
/// of M. Zaghloul's [Algorithm 916](http://dx.doi.org/10.1145/2049673.2049679).
///
/// This is effectively a Rust port of `libcerf` and the
/// [Faddeeva package](http://ab-initio.mit.edu/faddeeva/)
/// written by Steven G. Johnson at MIT. (MIT License)
pub fn w(re: f64, im: f64) -> (f64, f64) {
    println!("\n----------------------------------------------------");
    use num_complex::Complex64 as c64;
    let mut res = c64::default();

    let z = c64::new(re, im);
    let x = re;
    let xabs = re.abs();
    let y = im;
    let yabs = im.abs();
    let z_sqared = xabs * xabs + y * y;

    /* **************************************************************** */
    /*     Case when |ℜ(z)| >> |ℑ(z)|: In this scenario we may have    */
    /* precision issues because the complex norm can be very accurate,  */
    /* but either of the components may be off by orders of magnitude.  */
    /* **************************************************************** */

    if yabs < 1e-8 * xabs {
        let wi = im_w_of_x(x);
        let e2 = if xabs > 27. {
            0.
        } else {
            f64::exp(-xabs * xabs)
        };
        if yabs == 0. {
            return (e2, wi);
        }
        return (e2 + y * 2. * (x * wi - INV_SQRT_PI), wi - 2. * x * y * e2);
    }

    if xabs < 1e-8 * yabs {
        println!("DEBUG: CALLED ERFCX(Y).");
        let wr = erfcx(y);
        if xabs == 0. {
            return (wr, 0.);
        }
        return (wr, x * (2. * (INV_SQRT_PI - y * wr)));
    }

    // When z tends to zero, we use the MacLaurin series...

    if z_sqared < 0.053 {
        if z_sqared < 0.00689 {
            if z_sqared < 4e-7 {
                res = ((((
                              c64::new(5.0000000000000000e-01, 0.) ) * z // z^4
                          + c64::new(0., -7.5225277806367508e-01) ) * z // z^3
                         + c64::new(-1.0000000000000000e+00, 0.) ) * z // z^2
                        + c64::new(0., 1.1283791670955126e+00) ) * z // z^1
                    + 1.;

                return (res.re, res.im);
            }

            res = (((((((((((((
                                   c64::new(0., 5.3440090793734269e-04) ) * z // z^13
                               + c64::new(1.3888888888888889e-03, 0.) ) * z // z^12
                              + c64::new(0., -3.4736059015927274e-03) ) * z // z^11
                             + c64::new(-8.3333333333333332e-03, 0.) ) * z // z^10
                            + c64::new(0., 1.9104832458760001e-02) ) * z // z^9
                           + c64::new(4.1666666666666664e-02, 0.) ) * z // z^8
                          + c64::new(0., -8.5971746064419999e-02) ) * z // z^7
                         + c64::new(-1.6666666666666666e-01, 0.) ) * z // z^6
                        + c64::new(0., 3.0090111122547003e-01) ) * z // z^5
                       + c64::new(5.0000000000000000e-01, 0.) ) * z // z^4
                      + c64::new(0., -7.5225277806367508e-01) ) * z // z^3
                     + c64::new(-1.0000000000000000e+00, 0.) ) * z // z^2
                    + c64::new(0., 1.1283791670955126e+00) ) * z // z^1
                + 1.;

            return (res.re, res.im);
        }

        res = (((((((((((((((((((
				     c64::new(0., -8.8239572002038009e-07) ) * z // z^19
				 + c64::new(-2.7557319223985893e-06, 0.) ) * z // z^18
				+ c64::new(0., 8.3827593401936105e-06) ) * z // z^17
			       + c64::new(2.4801587301587302e-05, 0.) ) * z // z^16
			      + c64::new(0., -7.1253454391645692e-05) ) * z // z^15
			     + c64::new(-1.9841269841269841e-04, 0.) ) * z // z^14
			    + c64::new(0., 5.3440090793734269e-04) ) * z // z^13
			   + c64::new(1.3888888888888889e-03, 0.) ) * z // z^12
			  + c64::new(0., -3.4736059015927274e-03) ) * z // z^11
			 + c64::new(-8.3333333333333332e-03, 0.) ) * z // z^10
			+ c64::new(0., 1.9104832458760001e-02) ) * z // z^9
		       + c64::new(4.1666666666666664e-02, 0.) ) * z // z^8
		      + c64::new(0., -8.5971746064419999e-02) ) * z // z^7
		     + c64::new(-1.6666666666666666e-01, 0.) ) * z // z^6
		    + c64::new(0., 3.0090111122547003e-01) ) * z // z^5
		   + c64::new(5.0000000000000000e-01, 0.) ) * z // z^4
		  + c64::new(0., -7.5225277806367508e-01) ) * z // z^3
		 + c64::new(-1.0000000000000000e+00, 0.) ) * z // z^2
		+ c64::new(0., 1.1283791670955126e+00) ) * z // z^1
	    + 1.;

        return (res.re, res.im);
    }

    if z_sqared >= 49. {
        res = c64::new(0., 0.);
        let xs = if y < 0. { -z.re } else { z.re }; // compute for -z if y < 0

        if z_sqared > 4.8e15 {
            // Scale to prevent overflow.
            if xabs > yabs {
                let yax = yabs / xs;
                let denom = INV_SQRT_PI / (xs + yax * yabs);
                res = c64::new(denom * yax, denom);
            } else if yabs.is_infinite() {
                res = if xabs.is_nan() || y < 0. {
                    c64::new(f64::NAN, f64::NAN)
                } else {
                    c64::new(0., 0.)
                };
                return (res.re, res.im);
            } else {
                let xya = xs / yabs;
                let denom = INV_SQRT_PI / (xya * xs + yabs);
                res = c64::new(denom, denom * xya);
            }
        } else {
            let zm2 = 1. / z_sqared; // 1/|z|^2
            let r = c64::new(yabs * zm2, xs * zm2); // i/z
            let zm4 = zm2 * zm2; // 1/|z|^4
            let r2 = c64::new(zm4 * (xs + yabs) * (xs - yabs), -2. * zm4 * xs * yabs); // 1/z^2

            if z_sqared > 540. {
                if z_sqared > 22500. {
                    res = ((((
				    1.0578554691520430e+00) * r2 // n=3
			     + 4.2314218766081724e-01) * r2 // n=2
			    + 2.8209479177387814e-01) * r2 // n=1
			   + 5.6418958354775628e-01)
                        * r; // n=0
                } else {
                    res = ((((((((((((
					    3.7877040075087948e+06) * r2 // n=11
				     + 3.6073371500083758e+05) * r2 // n=10
				    + 3.7971970000088164e+04) * r2 // n=9
				   + 4.4672905882456671e+03) * r2 // n=8
				  + 5.9563874509942218e+02) * r2 // n=7
				 + 9.1636730015295726e+01) * r2 // n=6
				+ 1.6661223639144676e+01) * r2 // n=5
			       + 3.7024941420321507e+00) * r2 // n=4
			      + 1.0578554691520430e+00) * r2 // n=3
			     + 4.2314218766081724e-01) * r2 // n=2
			    + 2.8209479177387814e-01) * r2 // n=1
			   + 5.6418958354775628e-01)
                        * r; // n=0
                }
            } else {
                res = ((((((((((((((((((((
                                            8.8249260943025370e+15) * r2 // n=19
                                         + 4.7702303212446150e+14) * r2 // n=18
                                        + 2.7258458978540656e+13) * r2 // n=17
                                       + 1.6520278168812520e+12) * r2 // n=16
                                      + 1.0658243979879044e+11) * r2 // n=15
                                     + 7.3505130895717545e+09) * r2 // n=14
                                    + 5.4448245107938921e+08) * r2 // n=13
                                   + 4.3558596086351141e+07) * r2 // n=12
                                  + 3.7877040075087948e+06) * r2 // n=11
                                 + 3.6073371500083758e+05) * r2 // n=10
                                + 3.7971970000088164e+04) * r2 // n=9
                               + 4.4672905882456671e+03) * r2 // n=8
                              + 5.9563874509942218e+02) * r2 // n=7
                             + 9.1636730015295726e+01) * r2 // n=6
                            + 1.6661223639144676e+01) * r2 // n=5
                           + 3.7024941420321507e+00) * r2 // n=4
                          + 1.0578554691520430e+00) * r2 // n=3
                         + 4.2314218766081724e-01) * r2 // n=2
                        + 2.8209479177387814e-01) * r2 // n=1
                       + 5.6418958354775628e-01)
                    * r; // n=0
            }
        }
        if y < 0. {
            // Use w(z) = 2.0*exp(-z*z) - w(-z),
            // but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
            res = 2.0 * c64::new((yabs - xs) * (xs + yabs), 2. * xs * y).exp() - res;
            return (res.re, res.im);
        } else {
            return (res.re, res.im);
        }
    }

    if xabs.is_nan() {
        return (xabs, xabs);
    }
    if y.is_nan() {
        return (y, y);
    }

    println!("DEBUG: USING TAYLOR EXPANSION.");
    let i_tile =
        2 * ((INV_A * xabs) as isize * N_X_COVER as isize + (INV_A * yabs) as isize) as usize;
    let kp = TILES[i_tile];
    assert!(kp >= 0, "kP is always >= 0 in the libcerf code.");
    let idx = 2 * (NTAY + 1) * kp as usize;
    let t = TAYLOR_COEFFS[idx];
    let t_next = TAYLOR_COEFFS[idx + 1];
    dbg!(t, t_next);
    let dz = c64::new(xabs - t, yabs - t_next);

    // Something is wrong here...
    let nk = TILES[i_tile + 1];
    assert!(nk >= 0); // Not sure how a negative value would affect the loop...
    res = c64::new(
        TAYLOR_COEFFS[(idx as isize + 2 * nk as isize) as usize],
        TAYLOR_COEFFS[(idx as isize + 2 * nk as isize) as usize + 1],
    );
    // Equivalent to `for (int k = Nk-1; k >= 1; --k)` in libcerf.
    for k in (1..nk).rev() {
        res = res * dz
            + c64::new(
                TAYLOR_COEFFS[(idx as isize + (2 * k as isize)) as usize],
                TAYLOR_COEFFS[(idx as isize + (2 * k as isize)) as usize + 1],
            );
    }

    if y < 0. {
        if x < 0. {
            res = 2.0 * c64::new((y - x) * (x + y), -2. * x * y).exp() - res;
            return (res.re, res.im);
        }
        res = 2. * (c64::new((y - x) * (x + y), -2. * x * y)).exp() - c64::new(res.re, -(res.im));
        return (res.re, res.im);
    }
    if x < 0. {
        res = c64::new(res.re, -res.im);
    }
    return (res.re, res.im);
}

/// Computes the scaled Dawson integral of x,
/// which is 2 * Dawson(x) / sqrt(pi).
/// This equals the imaginary part of the Faddeeva
/// function of x, for real x.
fn im_w_of_x(x: f64) -> f64 {
    println!("DEBUG: CALLED DAWSON(X).");
    2. * dawson(x) * INV_SQRT_PI
}
