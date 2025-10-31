# PURUSPE

[![On crates.io](https://img.shields.io/crates/v/puruspe.svg)](https://crates.io/crates/puruspe)
[![On docs.rs](https://docs.rs/puruspe/badge.svg)](https://docs.rs/puruspe)

**PUR**e **RUS**t **SPE**cial function library.

## Features

- Gamma functions
- Beta functions
- Error functions
- Bessel functions
- Lambert W functions
- Dawson function

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
puruspe = "0.4.1"
```

## Example

```rust
use puruspe::gamma;

fn main() {
    let x = 5.0;
    let result = gamma(x);
    println!("Gamma({}) = {}", x, result);
}
```

## Available Functions

### Gamma Functions
- `gamma(x)`: Gamma function
- `ln_gamma(x)`: Natural logarithm of the gamma function
- `gammp(a, x)`: Regularized lower incomplete gamma function P(a,x)
- `gammq(a, x)`: Regularized upper incomplete gamma function Q(a,x)
- `invgammp(p, a)`: Inverse of the regularized lower incomplete gamma function

### Beta Functions
- `beta(z, w)`: Beta function
- `betai(a, b, x)`: Regularized incomplete beta function I_x(a,b)
- `invbetai(p, a, b)`: Inverse of the regularized incomplete beta function

### Error Functions
- `erf(x)`: Error function
- `erfc(x)`: Complementary error function
- `erfcx(x)`: Scaled complementary error function
- `inverf(p)`: Inverse error function
- `inverfc(p)`: Inverse complementary error function

### Faddeeva Function
- `faddeeva(x)`: Faddeeva function

### Bessel Functions
- `Jn(n, x)`: Bessel function of the first kind of integer order
- `Yn(n, x)`: Bessel function of the second kind of integer order
- `In(n, x)`: Modified Bessel function of the first kind of integer order
- `Kn(n, x)`: Modified Bessel function of the second kind of integer order
- `Jnu_Ynu(nu, x)`: Bessel functions of the first and second kind of fractional order
- `Inu_Knu(nu, x)`: Modified Bessel functions of the first and second kind of fractional order
- `besseljy(nu, x)`: Bessel functions of the first and second kind with derivatives
- `besselik(nu, x)`: Modified Bessel functions of the first and second kind with derivatives

### Lambert W Functions
- `lambert_w(k, z_re, z_im)`: Branch `k` of the complex valued Lambert W function computed at point `z`
 in the complex plane.
- `lambert_w0(x)`: The principal branch of the Lambert W function computed to 50 bits of accuracy.
- `lambert_wm1(x)`: The secondary branch of the Lambert W function computed to 50 bits of accuracy.
- `sp_lambert_w0(x)`: The principal branch of the Lambert W function computed to 24 bits of accuracy on `f64`s.
- `sp_lambert_wm1(x)`: The secondary branch of the Lambert W function computed to 24 bits of accuracy on `f64`s.

### Dawson Function
- `dawson(x)`: Dawson's integral

## Precision

The precision of each function can vary depending on the input values and the complexity of the calculation.
For detailed information about the precision of specific functions, please refer to the test files in the `tests/` directory.

- Gamma function precision: see `tests/gamma_test.rs`
- Beta function precision: see `tests/beta_test.rs`
- Bessel function precision: see `tests/bessel_test.rs`
- Error function precision: see `tests/erf_test.rs`
- Faddeeva function precision: see `tests/faddeeva_test.rs`
- Lambert W function precision: see `tests/lambert_w_test.rs`
- Dawson function precision: see `tests/dawson_test.rs`

These test files contain comparisons between the results of our implementations and the corresponding functions in SciPy, a widely-used scientific computing library in Python.
This comparison provides insights into the precision of each function across various input ranges.

The test tables used for these comparisons are generated using SciPy, and the scripts for generating these tables can be found in the `scripts/` directory.
This allows for transparent verification and updating of our test cases.

Note that while we strive for high accuracy, the actual precision in your use case may differ slightly from the test cases.
If you require guaranteed precision for a specific input range, we recommend additional testing for your particular use case.

Also, please be aware that there might be small discrepancies between our implementation and SciPy's results due to differences in algorithms or internal precision.
These discrepancies are generally within acceptable margins for most applications, but if you need exact agreement with SciPy or any other specific implementation, you should perform detailed comparisons.

## Contributing

Contributions are welcome! Here are some ways you can contribute to this project:

1. Report bugs and request features by opening issues.
2. Submit pull requests to fix bugs or add new features.
3. Improve documentation or add examples.
4. Add new special functions or optimize existing ones.

When contributing code, please ensure that:

1. Your code follows the existing style of the project.
2. You add appropriate tests for your changes.
3. All tests pass when you run `cargo test`.

If you're adding new functions or making significant changes, you may need to update or create new test tables. You can use the Python scripts in the `scripts/` directory to generate these tables using SciPy.

## License

This project is dual-licensed under either of

* Apache License, Version 2.0, [LICENSE-APACHE](LICENSE-Apache2.0)
* MIT license [LICENSE-MIT](LICENSE-MIT)

at your option.

## Reference

- Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007.
