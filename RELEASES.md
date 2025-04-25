# Ver 0.4.1 (2025-04-25)

- Add all branches of the complex Lambert W function
  - Contributors: [**JSorngard**](https://github.com/JSorngard)
  - Related issue: [#40](https://github.com/Axect/puruspe/issues/40)
  - Related PR: [#45](https://github.com/Axect/puruspe/pull/45)
- Enhance CI / Fix proptest & lint
  - Contributors: [**JSorngard**](https://github.com/JSorngard)
  - Related issue: [#42](https://github.com/Axect/puruspe/issues/42)
  - Related PRs:
    - [#41](https://github.com/Axect/puruspe/pull/41)
    - [#43](https://github.com/Axect/puruspe/pull/43)
    - [#44](https://github.com/Axect/puruspe/pull/44)

# Ver 0.4.0 (2025-02-24)

## Implement Faddeeva function

- Contributors: [**ethanbarry**](https://github.com/ethanbarry), [inkydragon](https://github.com/inkydragon), [Axect](https://github.com/Axect), [JSorngard](https://github.com/JSorngard)
- Related PRs
  - [#30](https://github.com/Axect/puruspe/pull/30)
- Newly added dependencies
  - [num-complex](https://crates.io/crates/num-complex)

## Fix bugs

- [#32](https://github.com/Axect/puruspe/issues/32) : Fix swapped doc for `inverf` and `inverfc` (Thanks to [ethanbarry](https://github.com/ethanbarry))
- [#36](https://github.com/Axect/puruspe/issues/36) : Fix `gamma` function for `1f64` (Thanks to [arihant2math](https://github.com/arihant2math))

# Ver 0.3.0 (2024-10-17)

## Implement comprehensive Tests & Integrate `lambert_w`

- Contributors: [JSorngard](https://github.com/JSorngard), [Axect](https://github.com/Axect), [ethanbarry](https://github.com/ethanbarry)
- Related issue: [#12](https://github.com/Axect/puruspe/issues/12)
- Related PRs
  - [#10](https://github.com/Axect/puruspe/pull/10)
  - [#13](https://github.com/Axect/puruspe/pull/13)
  - [#15](https://github.com/Axect/puruspe/pull/15)
  - [#16](https://github.com/Axect/puruspe/pull/16)
  - [#17](https://github.com/Axect/puruspe/pull/17)
  - [#18](https://github.com/Axect/puruspe/pull/18)
  - [#19](https://github.com/Axect/puruspe/pull/19)
  - [#21](https://github.com/Axect/puruspe/pull/21)
  - [#22](https://github.com/Axect/puruspe/pull/22)
  - [#23](https://github.com/Axect/puruspe/pull/23)
  - [#24](https://github.com/Axect/puruspe/pull/24)
  - [#25](https://github.com/Axect/puruspe/pull/25)
  - [#26](https://github.com/Axect/puruspe/pull/26)
- Newly added dependencies
  - [lambert_w](https://crates.io/crates/lambert_w)

# Ver 0.2.5 (2024-08-03)

* Implement Dawson's integral ([#9](https://github.com/Axect/puruspe/pull/9)) (Thanks to [ethanbarry](https://github.com/ethanbarry))
  * `dawson(x: f64) -> f64`: Dawson's integral

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
