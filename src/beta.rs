use crate::{EPS, FPMIN, ln_gamma, Y, W};
const SWITCH: usize = 3000;

// =============================================================================
// Beta function
// =============================================================================
/// Beta function
pub fn beta(z: f64, w: f64) -> f64 {
    (ln_gamma(z) + ln_gamma(w) - ln_gamma(z + w)).exp()
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
