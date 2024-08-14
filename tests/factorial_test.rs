use puruspe::utils::factorial;
use approx::assert_ulps_eq;

#[test]
fn test_factorial() {
    for i in 0..30 {
        assert_ulps_eq!(
            factorial(i),
            (2..=i).map(|j| j as f64).product()
        );
    }
}