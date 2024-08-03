# PURUSPE

[![On crates.io](https://img.shields.io/crates/v/puruspe.svg)](https://crates.io/crates/puruspe)
[![On docs.rs](https://docs.rs/puruspe/badge.svg)](https://docs.rs/puruspe)

**PUR**e **RUS**t **SPE**cial function library.
There are no dependencies.

## Implemented Functions

### Gamma functions

* `ln_gamma` : Logarithmic gamma function
* `gamma` : Gamma function
* `gammp` : Regularized lower gamma function
* `gammq` : Regularized upper gamma function
* `invgammp` : Inverse regularized lower gamma function

### Beta function

* `beta` : Beta function
* `betai` : Regularized Incomplete beta function
* `invbetai` : Inverse regularized incomplete beta function

### Error functions

* `erf` : Error function
* `erfc` : Complementary Error function
* `inverf` : Inverse error function
* `inverfc` : Inverse complementary error function

### Bessel functions

* `Jn` : Bessel function of the first kind (integer order)
* `Yn` : Bessel function of the second kind (integer order)
* `In` : Modified Bessel function of the first kind (integer order)
* `Kn` : Modified Bessel function of the second kind (integer order)
* `Jnu_Ynu` : Bessel function of the first kind and second kind (fractional order)
* `Inu_Knu` : Modified Bessel function of the first kind and second kind (fractional order)
* `besseljy` : Bessel function of the first and second kinds (include derivatives)
* `besselik` : Modified Bessel function of the second kind (include derivatives)

### Dawson's integral

* `dawson`: Dawson's integral

### Note on precision of Bessel functions

Bessel functions exhibit poor precision for extremely large arguments.
For instance, Jn(1,1e10) yields only 6 correct digits, in contrast to the typical 14 digits of double precision.
This behavior is comparable to Python's scipy.special.j1, which is an optimized variant of the Bessel J function.

## Reference

*  Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007. 
