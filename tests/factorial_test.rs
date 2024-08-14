use puruspe::utils::factorial;
use approx::assert_ulps_eq;

#[test]
fn test_factorial() {
    for i in 0..30 {
        let mut fact = 1.0;
        for j in 2..=i {
            fact *= j as f64;
        }
        
        assert_ulps_eq!(factorial(i), fact);
    }
}