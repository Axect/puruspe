use std::collections::HashMap;
use std::f64::consts::PI;

// =============================================================================
// Import from other crates
// =============================================================================
mod local_lambert_w;
pub use local_lambert_w::{lambert_w0, lambert_wm1, sp_lambert_w0, sp_lambert_wm1};
// =============================================================================
// For testing
// =============================================================================
#[cfg(test)]
mod tests;

// =============================================================================
// Constants
// =============================================================================
const EPS: f64 = f64::EPSILON;
const FPMIN: f64 = f64::MIN_POSITIVE / EPS;
const G: f64 = 5f64;
const N: usize = 7;
const ASWITCH: usize = 100;
const NGAU: usize = 18;
const Y: [f64; 18] = [
    0.0021695375159141994,
    0.011413521097787704,
    0.027972308950302116,
    0.051727015600492421,
    0.082502225484340941,
    0.12007019910960293,
    0.16415283300752470,
    0.21442376986779355,
    0.27051082840644336,
    0.33199876341447887,
    0.39843234186401943,
    0.46931971407375483,
    0.54413605556657973,
    0.62232745288031077,
    0.70331500465597174,
    0.78649910768313447,
    0.87126389619061517,
    0.95698180152629142,
];
const W: [f64; 18] = [
    0.0055657196642445571,
    0.012915947284065419,
    0.020181515297735382,
    0.027298621498568734,
    0.034213810770299537,
    0.040875750923643261,
    0.047235083490265582,
    0.053244713977759692,
    0.058860144245324798,
    0.064039797355015485,
    0.068745323835736408,
    0.072941885005653087,
    0.076598410645870640,
    0.079687828912071670,
    0.082187266704339706,
    0.084078218979661945,
    0.085346685739338721,
    0.085983275670394821,
];
// Error functions
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
// Incomplete beta function
const SWITCH: usize = 3000;

// =============================================================================
// Incomplete Gamma function (regularized)
// =============================================================================
/// Incomplete Gamma function P(a,x)
///
/// $$P(a,x) = \frac{1}{\Gamma(a)} \int_0^x t^{a-1} e^{-t} dt$$
pub fn gammp(a: f64, x: f64) -> f64 {
    assert!(x >= 0f64 && a > 0f64, "Bad args in gammp");
    if x == 0f64 {
        0f64
    } else if (a as usize) >= ASWITCH {
        // Quadrature
        gammpapprox(a, x, IncGamma::P)
    } else if x < a + 1f64 {
        // Series representation
        gser(a, x)
    } else {
        // Continued fraction representation
        1f64 - gcf(a, x)
    }
}

/// Incomplete Gamma function Q(a,x)
///
/// $$Q(a,x) = 1 - P(a,x)$$
pub fn gammq(a: f64, x: f64) -> f64 {
    assert!(x >= 0f64 && a > 0f64, "Bad args in gammp");
    if x == 0f64 {
        1f64
    } else if (a as usize) >= ASWITCH {
        // Quadrature
        gammpapprox(a, x, IncGamma::Q)
    } else if x < a + 1f64 {
        // Series representation
        1f64 - gser(a, x)
    } else {
        // Continued fraction representation
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
            return sum * (-x + a * x.ln() - gln).exp();
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
    (-x + a * x.ln() - gln).exp() * h
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
    let mut sum = 0f64;
    let mut t: f64;
    for j in 0..NGAU {
        t = x + (xu - x) * Y[j];
        sum += W[j] * (-(t - a1) + a1 * (t.ln() - lna1)).exp();
    }
    let ans = sum * (xu - x) * (a1 * (lna1 - 1f64) - gln).exp();
    match psig {
        IncGamma::P => {
            if ans > 0f64 {
                1f64 - ans
            } else {
                -ans
            }
        }
        IncGamma::Q => {
            if ans >= 0f64 {
                ans
            } else {
                1f64 + ans
            }
        }
    }
}

/// Iunverse Incomplete Gamma function
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
const LG5N7: [f64; 7] = [
    1.000000000189712,
    76.18009172948503,
    -86.50532032927205,
    24.01409824118972,
    -1.2317395783752254,
    0.0012086577526594748,
    -0.00000539702438713199,
];

/// Logarithm Gamma
pub fn ln_gamma(z: f64) -> f64 {
    let z = z - 1f64;
    let base = z + G + 0.5;
    let mut s = 0f64;
    for i in 1..N {
        s += LG5N7[i] / (z + i as f64);
    }
    s += LG5N7[0];
    (2f64 * PI).sqrt().ln() + s.ln() - base + base.ln() * (z + 0.5)
}

/// Gamma function
pub fn gamma(z: f64) -> f64 {
    if z > 1f64 {
        let z_int = z as usize;
        if z - (z_int as f64) == 0f64 {
            return factorial(z_int - 1) as f64;
        }
    }

    if z < 0.5 {
        PI / ((PI * z).sin() * gamma(1f64 - z))
    } else {
        ln_gamma(z).exp()
    }
}

// =============================================================================
// Beta function
// =============================================================================
/// Beta function
pub fn beta(z: f64, w: f64) -> f64 {
    (ln_gamma(z) + ln_gamma(w) - ln_gamma(z + w)).exp()
}

// =============================================================================
// Error functions
// =============================================================================
/// Error function
pub fn erf(x: f64) -> f64 {
    if x >= 0f64 {
        1.0 - erfccheb(x)
    } else {
        erfccheb(-x) - 1f64
    }
}

/// Complementary error function
pub fn erfc(x: f64) -> f64 {
    if x >= 0f64 {
        erfccheb(x)
    } else {
        2f64 - erfccheb(-x)
    }
}

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

/// Inverse of complementary error function
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

pub fn inverf(p: f64) -> f64 {
    inverfc(1f64 - p)
}

// =============================================================================
// Incomplete Beta function
// =============================================================================
pub fn betai(a: f64, b: f64, x: f64) -> f64 {
    assert!(a > 0f64 && b > 0f64, "Bad a or b in routine betai");
    assert!(x >= 0f64 && x <= 1f64, "Bad x in routine betai");
    if x == 0f64 || x == 1f64 {
        return x;
    }
    let switch = SWITCH as f64;
    if a > switch && b > switch {
        return betaiapprox(a, b, x);
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1f64 - x).ln()).exp();
    if x < (a + 1f64) / (a + b * 2f64) {
        bt * betacf(a, b, x) / a
    } else {
        1f64 - bt * betacf(b, a, 1f64 - x) / b
    }
}

/// Continued fraction beta
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1f64;
    let qam = a - 1f64;
    let mut c = 1f64;
    let mut d = 1f64 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1f64 / d;
    let mut h = d;
    for m in 1..10000 {
        let m = m as f64;
        let m2 = 2f64 * m;
        let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1f64 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1f64 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1f64 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1f64 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1f64 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1f64 / d;
        let del = d * c;
        h *= del;
        if (del - 1f64).abs() <= EPS {
            break;
        }
    }
    h
}

/// Incomplete beta by Gauss Legendre quadrature
fn betaiapprox(a: f64, b: f64, x: f64) -> f64 {
    let a1 = a - 1f64;
    let b1 = b - 1f64;
    let mu = a / (a + b);
    let lnmu = mu.ln();
    let lnmuc = (1f64 - mu).ln();
    let mut t = (a * b / ((a + b).powi(2) * (a + b + 1f64))).sqrt();
    let xu = if x > a / (a + b) {
        if x >= 1f64 {
            return 1f64;
        }
        1f64.min((mu + 10f64 * t).max(x + 5f64 * t))
    } else {
        if x <= 0f64 {
            return 0f64;
        }
        0f64.max((mu - 10f64 * t).min(x - 5f64 * t))
    };
    let mut sum = 0f64;
    for j in 0..18 {
        t = x + (xu - x) * Y[j];
        sum += W[j] * (a1 * (t.ln() - lnmu) + b1 * (1f64 - t).ln() - lnmuc).exp();
    }
    let ans = sum
        * (xu - x)
        * (a1 * lnmu - ln_gamma(a) + b1 * lnmuc - ln_gamma(b) + ln_gamma(a + b)).exp();
    if ans > 0f64 {
        1f64 - ans
    } else {
        -ans
    }
}

pub fn invbetai(p: f64, a: f64, b: f64) -> f64 {
    let a1 = a - 1f64;
    let b1 = b - 1f64;
    let mut t: f64;
    let mut x: f64;
    let mut u: f64;
    if p <= 0f64 {
        return 0f64;
    } else if p >= 1f64 {
        return 1f64;
    } else if a >= 1f64 && b >= 1f64 {
        let pp = if p < 0.5 { p } else { 1f64 - p };
        t = (-2f64 * pp.ln()).sqrt();
        x = (2.30753 + t * 0.27061) / (1f64 + t * (0.99229 + t * 0.04481)) - t;
        if p < 0.5 {
            x = -x;
        }
        let al = (x.powi(2) - 3f64) / 6f64;
        let h = 2f64 / (1f64 / (2f64 * a - 1f64) + 1f64 / (2f64 * b - 1f64));
        let w = (x * (al + h).sqrt() / h)
            - (1f64 / (2f64 * b - 1f64) - 1f64 / (2f64 * a - 1f64))
                * (al + 5f64 / 6f64 - 2f64 / (3f64 * h));
        x = a / (a + b * (2f64 * w).exp());
    } else {
        let lna = (a / (a + b)).ln();
        let lnb = (b / (a + b)).ln();
        t = (a * lna).exp() / a;
        u = (b * lnb).exp() / b;
        let w = t + u;
        x = if p < t / w {
            (a * w * p).powf(1f64 / a)
        } else {
            1f64 - (b * w * (1f64 - p)).powf(1f64 / b)
        };
    }
    let afac = -ln_gamma(a) - ln_gamma(b) + ln_gamma(a + b);
    for j in 0..10 {
        if x == 0f64 || x == 1f64 {
            return x;
        }
        let err = betai(a, b, x) - p;
        t = (a1 * x.ln() + b1 * (1f64 - x).ln() + afac).exp();
        u = err / t;
        t = u / (1f64 - 0.5 * 1f64.min(u * (a1 / x - b1 / (1f64 - x))));
        x -= t;
        if x <= 0f64 {
            x = 0.5 * (x + t);
        }
        if x >= 1f64 {
            x = 0.5 * (x + t + 1f64);
        }
        if t.abs() < EPS * x && j > 0 {
            break;
        }
    }
    x
}

// =============================================================================
// Dawson's integral
// =============================================================================
/// Dawson's integral for real values.
/// Accurate to about 2e-7 for small inputs.
pub fn dawson(x: f64) -> f64 {
    const NMAX: usize = 6;
    let mut c = [0_f64; NMAX];
    let h = 0.4;
    let a1 = 2. / 3.;
    let a2 = 0.4;
    let a3 = 2. / 7.;
    let n0: u64; // i is declared in the loops.
    let mut d1: f64;
    let mut d2: f64;
    let mut e1: f64;
    let e2: f64;
    let mut sum: f64;
    let x2: f64;
    let xp: f64;
    let xx: f64;
    let ans: f64;

    for i in 0..NMAX {
        c[i] = f64::exp(-((2. * i as f64 + 1.) * h).powf(2.));
    }

    if x.abs() < 0.2 {
        x2 = x * x;
        ans = x * (1.0 - a1 * x2 * (1.0 - a2 * x2 * (1.0 - a3 * x2)));
    } else {
        xx = x.abs();
        n0 = 2 * (0.5 * xx / h + 0.5) as u64; // Simply truncates the decimal places, as in NR.
        xp = xx - n0 as f64 * h;
        e1 = f64::exp(2.0 * xp * h);
        e2 = e1 * e1;
        d1 = n0 as f64 + 1.;
        d2 = d1 - 2.0;
        sum = 0.0;

        for i in 0..NMAX {
            d1 += 2.;
            d2 -= 2.;
            e1 *= e2;

            sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
        }

        ans = (1. / PI.sqrt()) * sign(f64::exp(-xp * xp), x) * sum;
    }
    return ans;
}

// =============================================================================
// Bessel functions
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

/// Bessel function of the second kind (n>=0)
///
/// # Arguments
/// * `n` - Order of Bessel function (non-negative integer)
/// * `x` - Argument of Bessel function
#[allow(non_snake_case)]
pub fn Yn(n: usize, x: f64) -> f64 {
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

/// Bessel function of the first kind (n>=0)
///
/// # Arguments
/// * `n` - Order of Bessel function (non-negative integer)
/// * `x` - Argument of Bessel function
#[allow(non_snake_case)]
pub fn Jn(n: usize, x: f64) -> f64 {
    let n_f64 = n as f64;
    let acc = 160f64;
    let iexp = f64::MAX_EXP / 2;

    if n == 0 {
        return j0(x);
    } else if n == 1 {
        return j1(x);
    } else if x.powi(2) <= 8.0 * f64::MIN_POSITIVE {
        return 0f64;
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
            let m = 2 * ((n + (acc * n_f64).sqrt() as usize) / 2);
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

/// Modified Bessel function of the second kind (n>=0)
///
/// # Arguments
/// * `n` - Order of the Bessel function (non-negative integer)
/// * `x` - Argument of the Bessel function
#[allow(non_snake_case)]
pub fn Kn(n: usize, x: f64) -> f64 {
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

/// Modified Bessel function of the first kind (n>=0)
///
/// # Arguments
/// * `n` - Order of the Bessel function (non-negative integer)
/// * `x` - Argument of the Bessel function
#[allow(non_snake_case)]
pub fn In(n: usize, x: f64) -> f64 {
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
        for j in (1..2 * (n + (acc * n_f64).sqrt() as usize)).rev() {
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

/// Bessel functions for fractional order
///
/// # Arguments
/// * `nu` - Order of the Bessel function (nu >= 0)
/// * `x` - Argument of the Bessel function (x > 0)
///
/// # Returns
/// * `J_nu(x)` - Bessel function of the first kind
/// * `Y_nu(x)` - Bessel function of the second kind
/// * `J_nu'(x)` - Derivative of the Bessel function of the first kind
/// * `Y_nu'(x)` - Derivative of the Bessel function of the second kind
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

/// Cached Bessel functions for fractional order
///
/// # Arguments
/// * `nu` - Order of the Bessel function (nu >= 0)
/// * `x` - Argument of the Bessel function (x > 0)
/// * `cache` - Cache for the Bessel functions (HashMap)
///
/// # Returns
/// * `J_nu(x)` - Bessel function of the first kind
/// * `Y_nu(x)` - Bessel function of the second kind
/// * `J_nu'(x)` - Derivative of the Bessel function of the first kind
/// * `Y_nu'(x)` - Derivative of the Bessel function of the second kind
pub fn cached_besseljy(
    nu: f64,
    x: f64,
    cache: &mut HashMap<(u64, u64), (f64, f64, f64, f64)>,
) -> (f64, f64, f64, f64) {
    if let Some(&res) = cache.get(&(nu.to_bits(), x.to_bits())) {
        res
    } else {
        let res = besseljy(nu, x);
        cache.insert((nu.to_bits(), x.to_bits()), res);
        res
    }
}

/// Modified Bessel functions for fractional order
///
/// # Arguments
/// * `nu` - Order of the Bessel function (nu >= 0)
/// * `x` - Argument of the Bessel function (x > 0)
///
/// # Returns
/// * `I_nu(x)` - Modified Bessel function of the first kind
/// * `K_nu(x)` - Modified Bessel function of the second kind
/// * `I_nu'(x)` - Derivative of the modified Bessel function of the first kind
/// * `K_nu'(x)` - Derivative of the modified Bessel function of the second kind
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
            delh = (b * d - 1.0) * delh;
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
        h = a1 * h;
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

/// Cached modified Bessel functions for fractional order
pub fn cached_besselik(
    nu: f64,
    x: f64,
    cache: &mut HashMap<(u64, u64), (f64, f64, f64, f64)>,
) -> (f64, f64, f64, f64) {
    if let Some(&res) = cache.get(&(nu.to_bits(), x.to_bits())) {
        res
    } else {
        let res = besselik(nu, x);
        cache.insert((nu.to_bits(), x.to_bits()), res);
        res
    }
}

/// Bessel functions for fractional order
///
/// # Arguments
/// * `nu` - Order of the Bessel function (nu >= 0)
/// * `x` - Argument of the Bessel function (x > 0)
///
/// # Returns
/// * `J_nu(x)` - Bessel function of the first kind
/// * `Y_nu(x)` - Bessel function of the second kind
#[allow(non_snake_case)]
pub fn Jnu_Ynu(nu: f64, x: f64) -> (f64, f64) {
    let (jo, yo, _, _) = besseljy(nu, x);
    (jo, yo)
}

/// Cached Bessel functions for fractional order
///
/// # Arguments
/// * `nu` - Order of the Bessel function (nu >= 0)
/// * `x` - Argument of the Bessel function (x > 0)
/// * `cache` - Cache for the Bessel functions (HashMap)
///
/// # Returns
/// * `J_nu(x)` - Bessel function of the first kind
/// * `Y_nu(x)` - Bessel function of the second kind
#[allow(non_snake_case)]
pub fn cached_Jnu_Ynu(nu: f64, x: f64, cache: &mut HashMap<(u64, u64), (f64, f64)>) -> (f64, f64) {
    if let Some(&res) = cache.get(&(nu.to_bits(), x.to_bits())) {
        res
    } else {
        let res = Jnu_Ynu(nu, x);
        cache.insert((nu.to_bits(), x.to_bits()), res);
        res
    }
}

/// Modified Bessel functions for fractional order
///
/// # Arguments
/// * `nu` - Order of the Bessel function (nu >= 0)
/// * `x` - Argument of the Bessel function (x > 0)
///
/// # Returns
/// * `I_nu(x)` - Modified Bessel function of the first kind
/// * `K_nu(x)` - Modified Bessel function of the second kind
#[allow(non_snake_case)]
pub fn Inu_Knu(nu: f64, x: f64) -> (f64, f64) {
    let (io, ko, _, _) = besselik(nu, x);
    (io, ko)
}

/// Cached modified Bessel functions for fractional order
///
/// # Arguments
/// * `nu` - Order of the Bessel function (nu >= 0)
/// * `x` - Argument of the Bessel function (x > 0)
/// * `cache` - Cache for the Bessel functions (HashMap)
///
/// # Returns
/// * `I_nu(x)` - Modified Bessel function of the first kind
/// * `K_nu(x)` - Modified Bessel function of the second kind
#[allow(non_snake_case)]
pub fn cached_Inu_Knu(nu: f64, x: f64, cache: &mut HashMap<(u64, u64), (f64, f64)>) -> (f64, f64) {
    if let Some(&res) = cache.get(&(nu.to_bits(), x.to_bits())) {
        res
    } else {
        let res = Inu_Knu(nu, x);
        cache.insert((nu.to_bits(), x.to_bits()), res);
        res
    }
}

const XJ00: f64 = 5.783185962946785;
const XJ10: f64 = 3.047126234366209e1;
const XJ01: f64 = 1.468197064212389e1;
const XJ11: f64 = 4.921845632169460e1;
const TWOOPI: f64 = std::f64::consts::FRAC_2_PI;
const PIO4: f64 = std::f64::consts::FRAC_PI_4;

const J0R: [f64; 7] = [
    1.682397144220462e-4,
    2.058861258868952e-5,
    5.288947320067750e-7,
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
    1.052598413585270e-4,
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
    1.237980436358390,
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
    4.328946737100230e-5,
    1.271526296341915e-7,
    2.566305357932989e-10,
    3.477378203574266e-13,
    2.593535427519985e-16,
];
const J1PN: [f64; 5] = [
    1.0,
    1.014039111045313,
    2.426762348629863e-1,
    1.350308200342000e-2,
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
    6.427210537081400e-5,
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
    1.052598413585270e-4,
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
    1.237980436358390,
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
    3.015798735815980e-18,
];
const Y1PN: [f64; 5] = [
    1.0,
    1.014039111045313,
    2.426762348629863e-1,
    1.350308200342000e-2,
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
    3.826993559940360e-4,
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
    1.192273748120670e-1,
    1.947452015979746e-1,
    7.629241821600588e-2,
    8.474903580801549e-3,
    2.023821945835647e-4,
];
const I0QQ: [f64; 6] = [
    2.962898424533095e-1,
    4.866115913196384e-1,
    1.938352806477617e-1,
    2.261671093400046e-2,
    6.450448095075585e-4,
    1.529835782400450e-6,
];
const I1P: [f64; 14] = [
    5.000000000000000e-1,
    6.090824836578078e-2,
    2.407288574545340e-3,
    4.622311145544158e-5,
    5.161743818147913e-7,
    3.712362374847555e-9,
    1.833983433811517e-11,
    6.493125133990706e-14,
    1.693074927497696e-16,
    3.299609473102338e-19,
    4.813071975603122e-22,
    5.164275442089090e-25,
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
    1.518396076767770e-2,
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
    1.475731032429900e1,
    6.123767403223466e1,
    1.121012633939949e2,
    9.285288485892228e1,
    3.198289277679660e1,
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
    1.408295601966600e-1,
];
const K1PI: [f64; 5] = [
    0.5,
    5.598072040178741e-2,
    1.818666382168295e-3,
    2.397509908859959e-5,
    1.239567816344855e-7,
];
const K1QI: [f64; 3] = [
    9.870202601341150e-1,
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
    1.375094061153160e-2,
    6.774221332947002e-5,
];
const K1PP: [f64; 8] = [
    1.253314137315502,
    1.457171340220454e1,
    6.063161173098803e1,
    1.147386690867892e2,
    1.040442011439181e2,
    4.356596656837691e1,
    7.265230396353690,
    3.144418558991021e-1,
];
const K1QQ: [f64; 8] = [
    1.0,
    1.125154514806458e1,
    4.427488496597630e1,
    7.616113213117645e1,
    5.863377227890893e1,
    1.850303673841586e1,
    1.857244676566022,
    2.538540887654872e-2,
];

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
fn sign(a: f64, b: f64) -> f64 {
    if b.is_sign_positive() {
        a.abs()
    } else {
        -(a.abs())
    }
}
