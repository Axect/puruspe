//! This module provides implementations of Bessel functions and related functions.
//!
//! It includes the following main functions:
//! - `Jn`: Calculates the Bessel function of the first kind of integer order.
//! - `Yn`: Calculates the Bessel function of the second kind of integer order.
//! - `In`: Calculates the modified Bessel function of the first kind of integer order.
//! - `Kn`: Calculates the modified Bessel function of the second kind of integer order.
//! - `Jnu_Ynu`: Calculates Bessel functions of the first and second kind for non-integer order.
//! - `Inu_Knu`: Calculates modified Bessel functions of the first and second kind for non-integer order.

use crate::utils::{frexp, ldexp};
use std::collections::HashMap;

/// Calculates the Bessel function of the first kind of order n.
///
/// The Bessel function of the first kind is defined as the solution to the differential equation:
///
/// $$ x^2 \frac{d^2y}{dx^2} + x \frac{dy}{dx} + (x^2 - n^2)y = 0 $$
///
/// # Arguments
///
/// * `n` - The order of the Bessel function (non-negative integer)
/// * `x` - The input value (non-negative real number)
///
/// # Returns
///
/// The value of $J_n(x)$
#[allow(non_snake_case)]
pub fn Jn(n: u32, x: f64) -> f64 {
    let n_f64 = n as f64;
    let acc = 160f64;
    let iexp = f64::MAX_EXP / 2;

    if n == 0 {
        j0(x)
    } else if n == 1 {
        j1(x)
    } else if x.powi(2) <= 8.0 * f64::MIN_POSITIVE {
        0f64
    } else {
        let mut ans = 0f64;
        let ax = x.abs();
        if ax > n_f64 {
            let tox = 2.0 / ax;
            let mut bjm = j0(ax);
            let mut bj = j1(ax);
            for j in 1..n {
                let bjp = j as f64 * tox * bj - bjm;
                bjm = bj;
                bj = bjp;
            }
            ans = bj;
        } else {
            let tox = 2.0 / ax;
            let m = 2 * ((n + (acc * n_f64).sqrt() as u32) / 2);
            let mut jsum = false;
            let mut bjp = 0f64;
            let mut sum = 0f64;
            let mut bj = 1f64;
            for j in (1..=m).rev() {
                let bjm = j as f64 * tox * bj - bjp;
                bjp = bj;
                bj = bjm;
                let (_, k) = frexp(bj);
                if k > iexp {
                    bj = ldexp(bj, -iexp);
                    bjp = ldexp(bjp, -iexp);
                    ans = ldexp(ans, -iexp);
                    sum = ldexp(sum, -iexp);
                }
                if jsum {
                    sum += bj;
                }
                jsum = !jsum;
                if j == n {
                    ans = bjp;
                }
            }
            sum = 2f64 * sum - bj;
            ans /= sum;
        }
        if x < 0f64 && n % 2 == 1 {
            -ans
        } else {
            ans
        }
    }
}

/// Calculates the Bessel function of the second kind of order n.
///
/// The Bessel function of the second kind is another solution to the Bessel differential equation:
///
/// $$ x^2 \frac{d^2y}{dx^2} + x \frac{dy}{dx} + (x^2 - n^2)y = 0 $$
///
/// It is also known as the Neumann function.
///
/// # Arguments
///
/// * `n` - The order of the Bessel function (non-negative integer)
/// * `x` - The input value (positive real number)
///
/// # Returns
///
/// The value of $Y_n(x)$
#[allow(non_snake_case)]
pub fn Yn(n: u32, x: f64) -> f64 {
    if n == 0 {
        y0(x)
    } else if n == 1 {
        y1(x)
    } else {
        let tox = 2.0 / x;
        let mut by = y1(x);
        let mut bym = y0(x);
        for j in 1..n {
            let byp = j as f64 * tox * by - bym;
            bym = by;
            by = byp;
        }
        by
    }
}

/// Calculates the modified Bessel function of the first kind of order n.
///
/// The modified Bessel function of the first kind is defined as:
///
/// $$ I_n(x) = i^{-n} J_n(ix) $$
///
/// where J_n is the Bessel function of the first kind.
///
/// # Arguments
///
/// * `n` - The order of the Bessel function (non-negative integer)
/// * `x` - The input value (non-negative real number)
///
/// # Returns
///
/// The value of I_n(x)
#[allow(non_snake_case)]
pub fn In(n: u32, x: f64) -> f64 {
    let acc = 200f64;
    let iexp = f64::MAX_EXP / 2;
    let n_f64 = n as f64;
    if n == 0 {
        i0(x)
    } else if n == 1 {
        i1(x)
    } else if x.powi(2) <= 8.0 * f64::MIN_POSITIVE {
        0f64
    } else {
        let tox = 2f64 / x.abs();
        let mut bip = 0f64;
        let mut ans = 0f64;
        let mut bi = 1f64;
        for j in (1..2 * (n + (acc * n_f64).sqrt() as u32)).rev() {
            let bim = j as f64 * tox * bi + bip;
            bip = bi;
            bi = bim;
            let (_, k) = frexp(bi);
            if k > iexp {
                bi = ldexp(bi, -iexp);
                bip = ldexp(bip, -iexp);
                ans = ldexp(ans, -iexp);
            }
            if j == n {
                ans = bip;
            }
        }
        ans *= i0(x) / bi;
        if x < 0f64 && n % 2 == 1 {
            -ans
        } else {
            ans
        }
    }
}

/// Calculates the modified Bessel function of the second kind of order n.
///
/// The modified Bessel function of the second kind is defined as:
///
/// $$ K_n(x) = \frac{\pi}{2} \frac{I_{-n}(x) - I_n(x)}{\sin(n\pi)} $$
///
/// where I_n is the modified Bessel function of the first kind.
///
/// # Arguments
///
/// * `n` - The order of the Bessel function (non-negative integer)
/// * `x` - The input value (positive real number)
///
/// # Returns
///
/// The value of $K_n(x)$
#[allow(non_snake_case)]
pub fn Kn(n: u32, x: f64) -> f64 {
    if n == 0 {
        k0(x)
    } else if n == 1 {
        k1(x)
    } else {
        let tox = 2.0 / x;
        let mut bkm = k0(x);
        let mut bk = k1(x);
        for j in 1..n {
            let bkp = j as f64 * tox * bk + bkm;
            bkm = bk;
            bk = bkp;
        }
        bk
    }
}

/// Calculates the Bessel functions of the first and second kind for non-integer order
///
/// # Arguments
///
/// - `nu` - The order of the Bessel function (non-negative real number)
/// - `x` - The input value (positive real number)
///
/// # Returns
///
/// - `J_nu(x)` - Bessel function of the first kind
/// - `Y_nu(x)` - Bessel function of the second kind
/// - `J_nu'(x)` - Derivative of the Bessel function of the first kind
/// - `Y_nu'(x)` - Derivative of the Bessel function of the second kind
///
/// # Panics
///
/// Panics if `x` ≤ 0 or if `nu` < 0.
/// Also panics if `x` is too large or if the function fails to converge.
pub fn besseljy(nu: f64, x: f64) -> (f64, f64, f64, f64) {
    const MAXIT: usize = 10000;
    const EPS: f64 = f64::EPSILON;
    const FPMIN: f64 = f64::MIN_POSITIVE / EPS;
    const XMIN: f64 = 2.0;
    const PI: f64 = std::f64::consts::PI;

    if x <= 0f64 || nu < 0f64 {
        panic!("bad arguments in besseljy");
    }
    let nl = if x < XMIN {
        (nu + 0.5) as usize
    } else {
        0usize.max((nu - x + 1.5) as usize)
    };
    let xmu = nu - nl as f64;
    let xmu2 = xmu.powi(2);
    let xi = 1.0 / x;
    let xi2 = 2.0 * xi;
    let w = xi2 / PI;
    let mut isign = 1i32;
    let mut h = nu * xi;
    if h < FPMIN {
        h = FPMIN;
    }
    let mut b = xi2 * nu;
    let mut d = 0f64;
    let mut c = h;
    let mut i = 0usize;
    while i < MAXIT {
        b += xi2;
        d = b - d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = b - 1.0 / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = c * d;
        h *= del;
        if d < 0f64 {
            isign = -isign;
        }
        if (del - 1f64).abs() <= EPS {
            break;
        }
        i += 1;
    }
    if i >= MAXIT {
        panic!("x too large in besseljy; try asymptotic expansion");
    }
    let mut rjl = isign as f64 * FPMIN;
    let mut rjpl = h * rjl;
    let rjl1 = rjl;
    let rjp1 = rjpl;

    let rjmu: f64;
    let mut rymu: f64;
    let rymup: f64;
    let mut ry1: f64;

    let mut fact = nu * xi;
    for _ in (0..nl).rev() {
        let rjtemp = fact * rjl + rjpl;
        fact -= xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
    }
    if rjl == 0f64 {
        rjl = EPS;
    }
    let f = rjpl / rjl;
    if x < XMIN {
        let x2 = 0.5 * x;
        let pimu = PI * xmu;
        fact = if pimu.abs() < EPS {
            1f64
        } else {
            pimu / pimu.sin()
        };
        d = -x2.ln();
        let e = xmu * d;
        let fact2 = if e.abs() < EPS { 1f64 } else { e.sinh() / e };
        let xx = 8.0 * xmu.powi(2) - 1f64;
        let gam1 = chebev(xx, &C1, 7);
        let gam2 = chebev(xx, &C2, 8);
        let gampl = gam2 - xmu * gam1;
        let gammi = gam2 + xmu * gam1;
        let mut ff = TWOOPI * fact * (gam1 * e.cosh() + gam2 * fact2 * d);
        let e = e.exp();
        let mut p = e / (gampl * PI);
        let mut q = 1f64 / (e * PI * gammi);
        let pimu2 = 0.5 * pimu;
        let fact3 = if pimu2.abs() < EPS {
            1f64
        } else {
            pimu2.sin() / pimu2
        };
        let r = PI * pimu2 * fact3 * fact3;
        let mut c = 1f64;
        let d = -x2.powi(2);
        let mut sum = ff + r * q;
        let mut sum1 = p;
        let mut i = 1usize;
        while i <= MAXIT {
            let i_f64 = i as f64;
            ff = (i_f64 * ff + p + q) / (i_f64 * i_f64 - xmu2);
            c *= d / i_f64;
            p /= i_f64 - xmu;
            q /= i_f64 + xmu;
            let del = c * (ff + r * q);
            sum += del;
            let del1 = c * p - i_f64 * del;
            sum1 += del1;
            if del.abs() < (1f64 + sum.abs()) * EPS {
                break;
            }
            i += 1;
        }
        if i > MAXIT {
            panic!("besseljy: failure to converge in cf1");
        }
        rymu = -sum;
        ry1 = -sum1 * xi2;
        rymup = xmu * xi * rymu - ry1;
        rjmu = w / (rymup - f * rymu);
    } else {
        let mut a = 0.25 - xmu2;
        let mut p = -0.5 * xi;
        let mut q = 1f64;
        let br = 2f64 * x;
        let mut bi = 2f64;
        let mut fact = a * xi / (p * p + q * q);
        let mut cr = br + q * fact;
        let mut ci = bi + p * fact;
        let mut den = br * br + bi * bi;
        let mut dr = br / den;
        let mut di = -bi / den;
        let dlr = cr * dr - ci * di;
        let dli = cr * di + ci * dr;
        let temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        let mut i = 1usize;
        while i < MAXIT {
            let i_f64 = i as f64;
            a += 2f64 * i_f64;
            bi += 2f64;
            dr = a * dr + br;
            di = a * di + bi;
            if dr.abs() + di.abs() < FPMIN {
                dr = FPMIN;
            }
            fact = a / (cr * cr + ci * ci);
            cr = br + cr * fact;
            ci = bi - ci * fact;
            if cr.abs() + ci.abs() < FPMIN {
                cr = FPMIN;
            }
            den = dr * dr + di * di;
            dr /= den;
            di /= -den;
            let dlr = cr * dr - ci * di;
            let dli = cr * di + ci * dr;
            let temp = p * dlr - q * dli;
            q = p * dli + q * dlr;
            p = temp;
            if (dlr - 1f64).abs() + dli.abs() < EPS {
                break;
            }
            i += 1;
        }
        if i >= MAXIT {
            panic!("besseljy: failure to converge in cf2");
        }
        let gam = (p - f) / q;
        let rjmu_temp = (w / ((p - f) * gam + q)).sqrt();
        rjmu = if rjl < 0f64 { -rjmu_temp } else { rjmu_temp };
        rymu = gam * rjmu;
        rymup = rymu * (p + q / gam);
        ry1 = xmu * xi * rymu - rymup;
    }
    let fact = rjmu / rjl;
    let jo = rjl1 * fact;
    let jpo = rjp1 * fact;
    for i in 1..=nl {
        let rytemp = (xmu + i as f64) * xi2 * ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
    }
    let yo = rymu;
    let ypo = nu * xi * rymu - ry1;

    (jo, yo, jpo, ypo)
}

/// This macro implements the convenience functions that are all the same among the cached Bessel functions.
/// It takes in the name of the cache struct as well as the type of that cache's values.
macro_rules! impl_cached_bessel_convenience_functions {
    ($name:ty, $val:ty) => {
        impl $name {
            /// Create a new cache.
            #[inline]
            pub fn new() -> Self {
                Self(HashMap::new())
            }

            /// Creates an empty cache, like [`new`](Self::new), but with at least the specified capacity.
            ///
            /// The cache will be able to hold at least `capacity` elements without reallocating.
            /// This method is allowed to allocate for more elements than `capacity`.
            /// If `capacity` is 0, the cache will not allocate.
            #[inline]
            pub fn with_capacity(capacity: usize) -> Self {
                Self(HashMap::with_capacity(capacity))
            }

            /// Returns true if the cache contains function values for the arguments `nu` and `x`.
            #[inline]
            pub fn contains(&self, nu: f64, x: f64) -> bool {
                self.0.contains_key(&(nu.to_bits(), x.to_bits()))
            }

            /// Clears the cache. Keeps the allocated memory for reuse.
            #[inline]
            pub fn clear(&mut self) {
                self.0.clear();
            }

            /// If the given `nu` and `x` values have function values associated with them in the cache,
            /// this function returns a reference to them.
            #[inline]
            pub fn get(&self, nu: f64, x: f64) -> Option<&$val> {
                self.0.get(&(nu.to_bits(), x.to_bits()))
            }

            /// Returns true if the cache is empty.
            #[inline]
            pub fn is_empty(&self) -> bool {
                self.0.is_empty()
            }

            /// Returns the number of elements in the cache.
            #[inline]
            pub fn len(&self) -> usize {
                self.0.len()
            }

            /// Shrink the cache as much as possible without removing elements.
            #[inline]
            pub fn shrink_to_fit(&mut self) {
                self.0.shrink_to_fit();
            }

            /// Reserves capacity for at least `additional` more elements to be inserted in the cache.
            /// The collection may reserve more space to speculatively avoid frequent reallocations.
            /// `capacity` will be greater than or equal to `self.len() + additional`.
            /// Does nothing if capacity is already sufficient.
            ///
            /// # Panics
            ///
            /// Panics if the new allocation size overflows `usize`.
            #[inline]
            pub fn reserve(&mut self, additional: usize) {
                self.0.reserve(additional)
            }

            /// Retains only the argument-function value pairs specified by the predicate.
            ///
            /// This removes all elements from the cache for which `f` returns false.
            /// The elements are visited in unspecified order.
            #[inline]
            pub fn retain<F>(&mut self, mut f: F)
            where
                F: FnMut((f64, f64), $val) -> bool,
            {
                self.0
                    .retain(|k, v| f((f64::from_bits(k.0), f64::from_bits(k.1)), *v))
            }
        }
    };
}

/// A cache of the values and derivaties of the Bessel functions
/// of the first and second kind for non-integer order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CachedBesselJY(HashMap<(u64, u64), (f64, f64, f64, f64)>);

impl CachedBesselJY {
    /// Get the values and derivatives of the Bessel functions of the first and second kind for non-integer order.
    ///
    /// If the values corresponding to the *exact* inputs are in the cache they are just returned, otherwise they are calculated and
    /// inserted into the cache.
    ///
    ///  # Arguments
    ///
    /// - `nu` - The order of the Bessel function (non-negative real number)
    /// - `x` - The input value (positive real number)
    ///
    /// # Returns
    ///
    /// `(J_nu(x), Y_nu(x), J_nu'(x), Y_nu'(x))` where:
    ///
    /// - `J_nu(x)` - Bessel function of the first kind
    /// - `Y_nu(x)` - Bessel function of the second kind
    /// - `J_nu'(x)` - Derivative of the Bessel function of the first kind
    /// - `Y_nu'(x)` - Derivative of the Bessel function of the second kind
    ///
    /// # Panics
    ///
    /// Panics if `x` is less than or equal to 0 or if `nu` is less than 0.
    /// Also panics if `x` is too large or the implementation fails to converge.
    pub fn besseljy(&mut self, nu: f64, x: f64) -> (f64, f64, f64, f64) {
        if let Some(&res) = self.0.get(&(nu.to_bits(), x.to_bits())) {
            res
        } else {
            let res = besseljy(nu, x);
            self.0.insert((nu.to_bits(), x.to_bits()), res);
            res
        }
    }
}

impl_cached_bessel_convenience_functions!(CachedBesselJY, (f64, f64, f64, f64));

/// Calculate the modified Bessel functions of the first and second kind for non-integer order
///
/// # Arguments
///
/// * `nu` - The order of the Bessel function (non-negative real number)
/// * `x` - The input value (positive real number)
///
/// # Returns
///
/// * `I_nu(x)` - Modified Bessel function of the first kind
/// * `K_nu(x)` - Modified Bessel function of the second kind
/// * `I_nu'(x)` - Derivative of the modified Bessel function of the first kind
/// * `K_nu'(x)` - Derivative of the modified Bessel function of the second kind
///
/// # Panics
///
/// Panics if `x` ≤ 0 or if `nu` < 0.
/// Also panics if `x` is too large or if the function fails to converge.
pub fn besselik(nu: f64, x: f64) -> (f64, f64, f64, f64) {
    const MAXIT: usize = 10000;
    const EPS: f64 = f64::EPSILON;
    const FPMIN: f64 = f64::MIN_POSITIVE / EPS;
    const XMIN: f64 = 2.0;
    const PI: f64 = std::f64::consts::PI;

    if x <= 0f64 || nu < 0f64 {
        panic!("bad arguments in besselik");
    }
    let nl = (nu + 0.5) as usize;
    let xmu = nu - nl as f64;
    let xmu2 = xmu.powi(2);
    let xi = 1.0 / x;
    let xi2 = 2.0 * xi;
    let mut h = nu * xi;
    if h < FPMIN {
        h = FPMIN;
    }
    let mut b = xi2 * nu;
    let mut d = 0f64;
    let mut c = h;
    let mut i = 0usize;
    while i < MAXIT {
        b += xi2;
        d = 1.0 / (b + d);
        c = b + 1.0 / c;
        let del = c * d;
        h *= del;
        if (del - 1f64).abs() <= EPS {
            break;
        }
        i += 1;
    }
    if i >= MAXIT {
        panic!("x too large in besselik; try asymptotic expansion");
    }
    let mut ril = FPMIN;
    let mut ripl = h * ril;
    let ril1 = ril;
    let rip1 = ripl;
    let mut fact = nu * xi;
    for _ in (0..nl).rev() {
        let ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
    }
    let f = ripl / ril;

    let mut rkmu: f64;
    let mut rk1: f64;
    if x < XMIN {
        let x2 = 0.5 * x;
        let pimu = PI * xmu;
        fact = if pimu.abs() < EPS {
            1f64
        } else {
            pimu / pimu.sin()
        };
        d = -x2.ln();
        let e = xmu * d;
        let fact2 = if e.abs() < EPS { 1f64 } else { e.sinh() / e };
        let xx = 8.0 * xmu.powi(2) - 1f64;
        let gam1 = chebev(xx, &C1, 7);
        let gam2 = chebev(xx, &C2, 8);
        let gampl = gam2 - xmu * gam1;
        let gammi = gam2 + xmu * gam1;
        let mut ff = fact * (gam1 * e.cosh() + gam2 * fact2 * d);
        let mut sum = ff;
        let e = e.exp();
        let mut p = 0.5 * e / gampl;
        let mut q = 0.5 / (e * gammi);
        let mut c = 1f64;
        let d = x2.powi(2);
        let mut sum1 = p;
        let mut i = 1usize;
        while i <= MAXIT {
            let i_f64 = i as f64;
            ff = (i_f64 * ff + p + q) / (i_f64 * i_f64 - xmu2);
            c *= d / i_f64;
            p /= i_f64 - xmu;
            q /= i_f64 + xmu;
            let del = c * ff;
            sum += del;
            let del1 = c * (p - i_f64 * ff);
            sum1 += del1;
            if del.abs() < sum.abs() * EPS {
                break;
            }
            i += 1;
        }
        if i > MAXIT {
            panic!("besselik: failure to converge in cf1");
        }
        rkmu = sum;
        rk1 = sum1 * xi2;
    } else {
        let mut b = 2.0 * (1.0 + x);
        let mut d = 1f64 / b;
        let mut h = d;
        let mut delh = d;
        let mut q1 = 0f64;
        let mut q2 = 1f64;
        let a1 = 0.25 - xmu2;
        let mut q = a1;
        let mut c = a1;
        let mut a = -a1;
        let mut s = 1.0 + q * delh;
        let mut i = 1usize;
        while i < MAXIT {
            let i_f64 = i as f64;
            a -= 2.0 * i_f64;
            c = -a * c / (i_f64 + 1.0);
            let q_new = (q1 - b * q2) / a;
            q1 = q2;
            q2 = q_new;
            q += c * q_new;
            b += 2.0;
            d = 1.0 / (b + a * d);
            delh *= b * d - 1.0;
            h += delh;
            let dels = q * delh;
            s += dels;
            if (dels / s).abs() <= EPS {
                break;
            }
            i += 1;
        }
        if i >= MAXIT {
            panic!("besselik: failure to converge in cf2");
        }
        h *= a1;
        rkmu = (PI / (2.0 * x)).sqrt() * (-x).exp() / s;
        rk1 = rkmu * (xmu + x + 0.5 - h) * xi;
    }
    let rkmup = xmu * xi * rkmu - rk1;
    let rimu = xi / (f * rkmu - rkmup);
    let io = (rimu * ril1) / ril;
    let ipo = (rimu * rip1) / ril;
    for i in 1..=nl {
        let i_f64 = i as f64;
        let rktemp = (xmu + i_f64) * xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    let ko = rkmu;
    let kpo = nu * xi * rkmu - rk1;

    (io, ko, ipo, kpo)
}

/// A cache of the values and derivaties of the modified Bessel functions
/// of the first and second kind for non-integer order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CachedBesselIK(HashMap<(u64, u64), (f64, f64, f64, f64)>);

impl CachedBesselIK {
    /// Calculates the values and derivatives of the modified Bessel functions of the first and second kind
    /// for non-integer order with cached results.
    ///
    /// If the values corresponding to the *exact* inputs are in the cache they are just returned, otherwise they are calculated and
    /// inserted into the cache.
    ///
    /// # Arguments
    ///
    /// - `nu` - The order of the Bessel function (non-negative real number)
    /// - `x` - The input value (positive real number)
    ///
    /// # Returns
    ///
    /// `(I_nu(x), K_nu(x), I_nu'(x), K_nu'(x))` where:
    ///
    /// - `I_nu(x)` - Modified Bessel function of the first kind
    /// - `K_nu(x)` - Modified Bessel function of the second kind
    /// - `I_nu'(x)` - Derivative of the modified Bessel function of the first kind
    /// - `K_nu'(x)` - Derivative of the modified Bessel function of the second kind
    ///
    /// # Panics
    ///
    /// Panics if `x` is less than or equal to 0, or if `nu` is less than zero.
    /// Also panics `x` is too large or the implementation fails to converge.
    pub fn besselik(&mut self, nu: f64, x: f64) -> (f64, f64, f64, f64) {
        if let Some(&res) = self.0.get(&(nu.to_bits(), x.to_bits())) {
            res
        } else {
            let res = besselik(nu, x);
            self.0.insert((nu.to_bits(), x.to_bits()), res);
            res
        }
    }
}

impl_cached_bessel_convenience_functions!(CachedBesselIK, (f64, f64, f64, f64));

/// Calculates the Bessel functions of the first and second kind for non-integer order.
///
/// # Arguments
///
/// - `nu` - The order of the Bessel function (non-negative real number)
/// - `x` - The input value (positive real number)
///
/// # Returns
///
/// A tuple (J_nu(x), Y_nu(x)) containing the values of the Bessel functions of the first and second kind
#[allow(non_snake_case)]
pub fn Jnu_Ynu(nu: f64, x: f64) -> (f64, f64) {
    let (jo, yo, _, _) = besseljy(nu, x);
    (jo, yo)
}

/// A cache of the values of the Bessel functions of the first and second kind for non-integer order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CachedJnuYnu(HashMap<(u64, u64), (f64, f64)>);

impl CachedJnuYnu {
    /// Calculates the values of the Bessel functions of the first and second kind for non-integer order.
    ///
    /// If the values corresponding to the *exact* inputs are in the cache they are just returned, otherwise they are calculated and
    /// inserted into the cache.
    ///
    /// # Arguments
    ///
    /// - `nu` - Order of the Bessel function (nu >= 0)
    /// - `x` - Argument of the Bessel function (x > 0)
    ///
    /// # Returns
    ///
    /// `(J_nu(x), Y_nu(x))`
    ///
    /// - `J_nu(x)` - Bessel function of the first kind
    /// - `Y_nu(x)` - Bessel function of the second kind
    ///
    /// # Panics
    ///
    /// Panics if `x` is less than or equal to zero or if `nu` is less than zero.
    /// Also panics if `x` is too large or the implementation fails to converge.
    #[allow(non_snake_case)]
    pub fn Jnu_Ynu(&mut self, nu: f64, x: f64) -> (f64, f64) {
        if let Some(&res) = self.0.get(&(nu.to_bits(), x.to_bits())) {
            res
        } else {
            let res = Jnu_Ynu(nu, x);
            self.0.insert((nu.to_bits(), x.to_bits()), res);
            res
        }
    }
}

impl_cached_bessel_convenience_functions!(CachedJnuYnu, (f64, f64));

/// Calculates the modified Bessel functions of the first and second kind for non-integer order.
///
/// # Arguments
///
/// - `nu` - The order of the Bessel function (non-negative real number)
/// - `x` - The input value (positive real number)
///
/// # Returns
///
/// A tuple (I_nu(x), K_nu(x)) containing the values of the modified Bessel functions of the first and second kind
#[allow(non_snake_case)]
pub fn Inu_Knu(nu: f64, x: f64) -> (f64, f64) {
    let (io, ko, _, _) = besselik(nu, x);
    (io, ko)
}

/// A cache of the values of the modified Bessel functions
/// of the first and second kind for non-integer order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CachedInuKnu(HashMap<(u64, u64), (f64, f64)>);

impl CachedInuKnu {
    /// Cached modified Bessel functions of the first and second kind for non-integer order
    ///
    /// If the values corresponding to the *exact* inputs are in the cache they are just returned, otherwise they are calculated and
    /// inserted into the cache.
    ///
    /// # Arguments
    ///
    /// - `nu` - Order of the Bessel function (nu >= 0)
    /// - `x` - Argument of the Bessel function (x > 0)
    ///
    /// # Returns
    ///
    /// `(I_nu(x), K_nu(x))`
    ///
    /// - `I_nu(x)` - Modified Bessel function of the first kind
    /// - `K_nu(x)` - Modified Bessel function of the second kind
    ///
    /// # Panics
    ///
    /// Panics if `x` is smaller than or equal to zero, or if `nu` is smaller than 0.
    /// Also panics if `x` is too large or if the implementation fails to converge.
    #[allow(non_snake_case)]
    pub fn Inu_Knu(&mut self, nu: f64, x: f64) -> (f64, f64) {
        if let Some(&res) = self.0.get(&(nu.to_bits(), x.to_bits())) {
            res
        } else {
            let res = Inu_Knu(nu, x);
            self.0.insert((nu.to_bits(), x.to_bits()), res);
            res
        }
    }
}

impl_cached_bessel_convenience_functions!(CachedInuKnu, (f64, f64));

// =============================================================================
// Building Blocks
// =============================================================================
/// Rational approximation
fn rat(x: f64, r: &[f64], s: &[f64], n: usize) -> (f64, f64, f64) {
    let y = x.powi(2);
    let z = 64f64 - y;
    let mut nump = r[n];
    let mut denp = s[n];
    for i in (0..n).rev() {
        nump = r[i] + z * nump;
        denp = s[i] + y * denp;
    }
    (y, nump, denp)
}

/// Asymptotic approximation
fn asp(
    ax: f64,
    pn: &[f64],
    pd: &[f64],
    qn: &[f64],
    qd: &[f64],
    fac: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let z = 8f64 / ax;
    let y = z.powi(2);
    let xx = ax - fac * PIO4;
    let mut nump = pn[4];
    let mut denp = pd[4];
    let mut numq = qn[4];
    let mut denq = qd[4];
    for i in (0..4).rev() {
        nump = pn[i] + y * nump;
        denp = pd[i] + y * denp;
        numq = qn[i] + y * numq;
        denq = qd[i] + y * denq;
    }
    (z, xx, nump, denp, numq, denq)
}

/// Bessel function of the first kind (n=0)
fn j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8f64 {
        let (y, nump, denp) = rat(x, &J0R, &J0S, 6);
        nump * (y - XJ00) * (y - XJ10) / denp
    } else {
        let (z, xx, nump, denp, numq, denq) = asp(ax, &J0PN, &J0PD, &J0QN, &J0QD, 1f64);
        (TWOOPI / ax).sqrt() * (xx.cos() * nump / denp - z * xx.sin() * numq / denq)
    }
}

/// Bessel function of the first kind (n=1)
fn j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8f64 {
        let (y, nump, denp) = rat(x, &J1R, &J1S, 6);
        x * nump * (y - XJ01) * (y - XJ11) / denp
    } else {
        let (z, xx, nump, denp, numq, denq) = asp(ax, &J1PN, &J1PD, &J1QN, &J1QD, 3f64);
        let ans = (TWOOPI / ax).sqrt() * (xx.cos() * nump / denp - z * xx.sin() * numq / denq);
        if x < 0f64 {
            -ans
        } else {
            ans
        }
    }
}

/// Bessel function of the second kind (n=0)
fn y0(x: f64) -> f64 {
    if x < 8f64 {
        let j0x = j0(x);
        let (_, nump, denp) = rat(x, &Y0R, &Y0S, 8);
        nump / denp + TWOOPI * j0x * x.ln()
    } else {
        let ax = x;
        let (z, xx, nump, denp, numq, denq) = asp(ax, &Y0PN, &Y0PD, &Y0QN, &Y0QD, 1f64);
        (TWOOPI / ax).sqrt() * (xx.sin() * nump / denp + z * xx.cos() * numq / denq)
    }
}

/// Bessel function of the second kind (n=1)
fn y1(x: f64) -> f64 {
    if x < 8f64 {
        let j1x = j1(x);
        let (_, nump, denp) = rat(x, &Y1R, &Y1S, 7);
        x * nump / denp + TWOOPI * (j1x * x.ln() - 1f64 / x)
    } else {
        let ax = x;
        let (z, xx, nump, denp, numq, denq) = asp(ax, &Y1PN, &Y1PD, &Y1QN, &Y1QD, 3f64);
        (TWOOPI / ax).sqrt() * (xx.sin() * nump / denp + z * xx.cos() * numq / denq)
    }
}

fn poly(x: f64, coef: &[f64], n: usize) -> f64 {
    let mut ans = coef[n];
    for i in (0..n).rev() {
        ans = ans * x + coef[i];
    }
    ans
}

/// Modified Bessel function of the first kind (n=0)
fn i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 15f64 {
        let y = x.powi(2);
        poly(y, &I0P, 13) / poly(225f64 - y, &I0Q, 4)
    } else {
        let z = 1f64 - 15f64 / ax;
        ax.exp() * poly(z, &I0PP, 4) / (ax.sqrt() * poly(z, &I0QQ, 5))
    }
}

/// Modified Bessel function of the first kind (n=1)
fn i1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 15f64 {
        let y = x.powi(2);
        x * poly(y, &I1P, 13) / poly(225f64 - y, &I1Q, 4)
    } else {
        let z = 1f64 - 15f64 / ax;
        let ans = x.exp() * poly(z, &I1PP, 4) / (ax.sqrt() * poly(z, &I1QQ, 5));
        if x < 0f64 {
            -ans
        } else {
            ans
        }
    }
}

/// Modified Bessel function of the second kind (n=0)
fn k0(x: f64) -> f64 {
    if x <= 1f64 {
        let z = x.powi(2);
        let term = poly(z, &K0PI, 4) * x.ln() / poly(1f64 - z, &K0QI, 2);
        poly(z, &K0P, 4) / poly(1f64 - z, &K0Q, 2) - term
    } else {
        let z = 1f64 / x;
        (-x).exp() * poly(z, &K0PP, 7) / (x.sqrt() * poly(z, &K0QQ, 7))
    }
}

/// Modified Bessel function of the second kind (n=1)
fn k1(x: f64) -> f64 {
    if x <= 1f64 {
        let z = x.powi(2);
        let term = poly(z, &K1PI, 4) * x.ln() / poly(1f64 - z, &K1QI, 2);
        x * (poly(z, &K1P, 4) / poly(1f64 - z, &K1Q, 2) + term) + 1f64 / x
    } else {
        let z = 1f64 / x;
        (-x).exp() * poly(z, &K1PP, 7) / (x.sqrt() * poly(z, &K1QQ, 7))
    }
}

fn chebev(x: f64, c: &[f64], m: usize) -> f64 {
    let mut d = 0f64;
    let mut dd = 0f64;
    for j in (1..m).rev() {
        let sv = d;
        d = 2f64 * x * d - dd + c[j];
        dd = sv;
    }
    x * d - dd + 0.5 * c[0]
}

const C1: [f64; 7] = [
    -1.142022680371168e0,
    6.5165112670737e-3,
    3.087090173086e-4,
    -3.4706269649e-6,
    6.9437664e-9,
    3.67795e-11,
    -1.356e-13,
];

const C2: [f64; 8] = [
    1.843740587300905e0,
    -7.68528408447867e-2,
    1.2719271366546e-3,
    -4.9717367042e-6,
    -3.31261198e-8,
    2.423096e-10,
    -1.702e-13,
    -1.49e-15,
];

const XJ00: f64 = 5.783185962946785;
const XJ10: f64 = 3.047126234366209e1;
const XJ01: f64 = 1.468197064212389e1;
const XJ11: f64 = 4.921_845_632_169_46e1;
const TWOOPI: f64 = core::f64::consts::FRAC_2_PI;
const PIO4: f64 = core::f64::consts::FRAC_PI_4;

const J0R: [f64; 7] = [
    1.682397144220462e-4,
    2.058861258868952e-5,
    5.288_947_320_067_75e-7,
    5.557173907680151e-9,
    2.865540042042604e-11,
    7.398972674152181e-14,
    7.925088479679688e-17,
];
const J0S: [f64; 7] = [
    1.0,
    1.019685405805929e-2,
    5.130296867064666e-5,
    1.659702063950243e-7,
    3.728997574317067e-10,
    5.709292619977798e-13,
    4.932979170744996e-16,
];
const J0PN: [f64; 5] = [
    9.999999999999999e-1,
    1.039698629715637,
    2.576910172633398e-1,
    1.504152485749669e-2,
    1.052_598_413_585_27e-4,
];
const J0PD: [f64; 5] = [
    1.0,
    1.040797262528109,
    2.588070904043728e-1,
    1.529954477721284e-2,
    1.168931211650012e-4,
];
const J0QN: [f64; 5] = [
    -1.562499999999992e-2,
    -1.920039317065641e-2,
    -5.827951791963418e-3,
    -4.372674978482726e-4,
    -3.895839560412374e-6,
];
const J0QD: [f64; 5] = [
    1.0,
    1.237_980_436_358_39,
    3.838793938147116e-1,
    3.100323481550864e-2,
    4.165515825072393e-4,
];
const J1R: [f64; 7] = [
    7.309637831891357e-5,
    3.551248884746503e-6,
    5.820673901730427e-8,
    4.500650342170622e-10,
    1.831596352149641e-12,
    3.891583573305035e-15,
    3.524978592527982e-18,
];
const J1S: [f64; 7] = [
    1.0,
    9.398354768446072e-3,
    4.328_946_737_100_23e-5,
    1.271526296341915e-7,
    2.566305357932989e-10,
    3.477378203574266e-13,
    2.593535427519985e-16,
];
const J1PN: [f64; 5] = [
    1.0,
    1.014039111045313,
    2.426762348629863e-1,
    1.350_308_200_342e-2,
    9.516522033988099e-5,
];
const J1PD: [f64; 5] = [
    1.0,
    1.012208056357845,
    2.408580305488938e-1,
    1.309511056184273e-2,
    7.746422941504713e-5,
];
const J1QN: [f64; 5] = [
    4.687499999999991e-2,
    5.652407388406023e-2,
    1.676531273460512e-2,
    1.231216817715814e-3,
    1.178364381441801e-5,
];
const J1QD: [f64; 5] = [
    1.0,
    1.210119370463693,
    3.626494789275638e-1,
    2.761695824829316e-2,
    3.240517192670181e-4,
];
const Y0R: [f64; 9] = [
    -7.653778457189104e-3,
    -5.854760129990403e-2,
    3.720671300654721e-4,
    3.313722284628089e-5,
    4.247761237036536e-8,
    -4.134562661019613e-9,
    -3.382190331837473e-11,
    -1.017764126587862e-13,
    -1.107646382675456e-16,
];
const Y0S: [f64; 9] = [
    1.0,
    1.125494540257841e-2,
    6.427_210_537_081_4e-5,
    2.462520624294959e-7,
    7.029372432344291e-10,
    1.560784108184928e-12,
    2.702374957564761e-15,
    3.468496737915257e-18,
    2.716600180811817e-21,
];
const Y0PN: [f64; 5] = [
    9.999999999999999e-1,
    1.039698629715637,
    2.576910172633398e-1,
    1.504152485749669e-2,
    1.052_598_413_585_27e-4,
];
const Y0PD: [f64; 5] = [
    1.0,
    1.040797262528109,
    2.588070904043728e-1,
    1.529954477721284e-2,
    1.168931211650012e-4,
];
const Y0QN: [f64; 5] = [
    -1.562499999999992e-2,
    -1.920039317065641e-2,
    -5.827951791963418e-3,
    -4.372674978482726e-4,
    -3.895839560412374e-6,
];
const Y0QD: [f64; 5] = [
    1.0,
    1.237_980_436_358_39,
    3.838793938147116e-1,
    3.100323481550864e-2,
    4.165515825072393e-4,
];
const Y1R: [f64; 8] = [
    -1.041835425863234e-1,
    -1.135093963908952e-5,
    2.212118520638132e-4,
    1.270981874287763e-6,
    -3.982892100836748e-8,
    -4.820712110115943e-10,
    -1.929392690596969e-12,
    -2.725259514545605e-15,
];
const Y1S: [f64; 8] = [
    1.0,
    1.186694184425838e-2,
    7.121205411175519e-5,
    2.847142454085055e-7,
    8.364240962784899e-10,
    1.858128283833724e-12,
    3.018846060781846e-15,
    3.015_798_735_815_98e-18,
];
const Y1PN: [f64; 5] = [
    1.0,
    1.014039111045313,
    2.426762348629863e-1,
    1.350_308_200_342e-2,
    9.516522033988099e-5,
];
const Y1PD: [f64; 5] = [
    1.0,
    1.012208056357845,
    2.408580305488938e-1,
    1.309511056184273e-2,
    7.746422941504713e-5,
];
const Y1QN: [f64; 5] = [
    4.687499999999991e-2,
    5.652407388406023e-2,
    1.676531273460512e-2,
    1.231216817715814e-3,
    1.178364381441801e-5,
];
const Y1QD: [f64; 5] = [
    1.0,
    1.210119370463693,
    3.626494789275638e-1,
    2.761695824829316e-2,
    3.240517192670181e-4,
];

const I0P: [f64; 14] = [
    9.999999999999997e-1,
    2.466405579426905e-1,
    1.478980363444585e-2,
    3.826_993_559_940_36e-4,
    5.395676869878828e-6,
    4.700912200921704e-8,
    2.733894920915608e-10,
    1.115830108455192e-12,
    3.301093025084127e-15,
    7.209167098020555e-18,
    1.166898488777214e-20,
    1.378948246502109e-23,
    1.124884061857506e-26,
    5.498556929587117e-30,
];
const I0Q: [f64; 5] = [
    4.463598170691436e-1,
    1.702205745042606e-3,
    2.792125684538934e-6,
    2.369902034785866e-9,
    8.965900179621208e-13,
];
const I0PP: [f64; 5] = [
    1.192_273_748_120_67e-1,
    1.947452015979746e-1,
    7.629241821600588e-2,
    8.474_903_580_801_55e-3,
    2.023821945835647e-4,
];
const I0QQ: [f64; 6] = [
    2.962898424533095e-1,
    4.866115913196384e-1,
    1.938352806477617e-1,
    2.261671093400046e-2,
    6.450448095075585e-4,
    1.529_835_782_400_45e-6,
];
const I1P: [f64; 14] = [
    5e-1,
    6.090824836578078e-2,
    2.407_288_574_545_34e-3,
    4.622311145544158e-5,
    5.161743818147913e-7,
    3.712362374847555e-9,
    1.833983433811517e-11,
    6.493125133990706e-14,
    1.693074927497696e-16,
    3.299609473102338e-19,
    4.813071975603122e-22,
    5.164_275_442_089_09e-25,
    3.846870021788629e-28,
    1.712948291408736e-31,
];
const I1Q: [f64; 5] = [
    4.665973211630446e-1,
    1.677754477613006e-3,
    2.583049634689725e-6,
    2.045930934253556e-9,
    7.166133240195285e-13,
];
const I1PP: [f64; 5] = [
    1.286515211317124e-1,
    1.930915272916783e-1,
    6.965689298161343e-2,
    7.345978783504595e-3,
    1.963602129240502e-4,
];
const I1QQ: [f64; 6] = [
    3.309385098860755e-1,
    4.878218424097628e-1,
    1.663088501568696e-1,
    1.473541892809522e-2,
    1.964131438571051e-4,
    -1.034524660214173e-6,
];
const K0PI: [f64; 5] = [
    1.0,
    2.346487949187396e-1,
    1.187082088663404e-2,
    2.150707366040937e-4,
    1.425433617130587e-6,
];
const K0QI: [f64; 3] = [
    9.847324170755358e-1,
    1.518_396_076_767_77e-2,
    8.362215678646257e-5,
];
const K0P: [f64; 5] = [
    1.159315156584126e-1,
    2.770731240515333e-1,
    2.066458134619875e-2,
    4.574734709978264e-4,
    3.454715527986737e-6,
];
const K0Q: [f64; 3] = [
    9.836249671709183e-1,
    1.627693622304549e-2,
    9.809660603621949e-5,
];
const K0PP: [f64; 8] = [
    1.253314137315499,
    1.475_731_032_429_9e1,
    6.123767403223466e1,
    1.121012633939949e2,
    9.285288485892228e1,
    3.198_289_277_679_66e1,
    3.595376024148513,
    6.160228690102976e-2,
];
const K0QQ: [f64; 8] = [
    1.0,
    1.189963006673403e1,
    5.027773590829784e1,
    9.496513373427093e1,
    8.318077493230258e1,
    3.181399777449301e1,
    4.443672926432041,
    1.408_295_601_966_6e-1,
];
const K1PI: [f64; 5] = [
    0.5,
    5.598072040178741e-2,
    1.818666382168295e-3,
    2.397509908859959e-5,
    1.239567816344855e-7,
];
const K1QI: [f64; 3] = [
    9.870_202_601_341_15e-1,
    1.292092053534579e-2,
    5.881933053917096e-5,
];
const K1P: [f64; 5] = [
    -3.079657578292062e-1,
    -8.109417631822442e-2,
    -3.477550948593604e-3,
    -5.385594871975406e-5,
    -3.110372465429008e-7,
];
const K1Q: [f64; 3] = [
    9.861813171751389e-1,
    1.375_094_061_153_16e-2,
    6.774221332947002e-5,
];
const K1PP: [f64; 8] = [
    1.253314137315502,
    1.457171340220454e1,
    6.063161173098803e1,
    1.147386690867892e2,
    1.040442011439181e2,
    4.356596656837691e1,
    7.265_230_396_353_69,
    3.144418558991021e-1,
];
const K1QQ: [f64; 8] = [
    1.0,
    1.125154514806458e1,
    4.427_488_496_597_63e1,
    7.616113213117645e1,
    5.863377227890893e1,
    1.850303673841586e1,
    1.857244676566022,
    2.538540887654872e-2,
];
