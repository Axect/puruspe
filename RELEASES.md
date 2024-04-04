# Ver 0.2.4 (2024-04-04)

* Implement Bessel & Modified Bessel functions
  * Integer order
    * `Jn(n: usize, x: f64) -> f64`: Bessel function of the first kind
    * `Yn(n: usize, x: f64) -> f64`: Bessel function of the second kind
    * `In(n: usize, x: f64) -> f64`: Modified Bessel function of the first kind
    * `Kn(n: usize, x: f64) -> f64`: Modified Bessel function of the second kind
  * Fractional order
    * `besseljy(nu: f64, x: f64) -> (f64, f64, f64, f64)`: Bessel functions (include derivatives)
    * `besselik(nu: f64, x: f64) -> (f64, f64, f64, f64)`: Modified Bessel functions (include derivatives)
    * `Jnu_Ynu(nu: f64, x: f64) -> (f64, f64)`: Bessel functions
    * `Inu_Knu(nu: f64, x: f64) -> (f64, f64)`: Modified Bessel functions

# Ver 0.2.3 (2024-02-18)

* Add `docs.rs` badge to README.md
* Add latex support to `gammp` and `gammq` in docs

# Ver 0.2.1 (2024-02-17)

* Update README.md
  * Fix typo (`gammap` -> `gammp`, `gammaq` -> `gammq`)

# Ver 0.2.0 (2022-05-02)

* Fix bugs (https://github.com/Axect/puruspe/issues/4)

# Ver 0.1.5 (2021-01-24)

* LICENSE Update : BSD-3-Clause -> MIT or Apache-2.0

# Ver 0.1.4 (2020-08-06)

* Make public `ln_gamma`. (Thanks to [rw](https://github.com/rw))

