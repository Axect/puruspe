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

### Note on precision

The precision is poor for huge arguments. E.g. Jn(1,1e10) gives only 6 correct digits, compared to 14 of typical double precision. Thus, it works similar to Python's `scipy.special.j1`, a fast variant of Bessel J function.

## Reference

*  Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007. 
