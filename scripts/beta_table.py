from scipy.special import beta, betainc, betaincinv

# a and b values to test for beta, betai, and invbetai
ab_values = [
    (0.1, 0.1),  # Both parameters very small, test edge case
    (0.5, 0.5),  # Both parameters 0.5, important in statistical applications
    (1.0, 1.0),  # Uniform distribution case
    (2.0, 2.0),  # Symmetric case
    (0.5, 5.0),  # Asymmetric case, a < 1, b > 1
    (5.0, 0.5),  # Asymmetric case, a > 1, b < 1
    (10.0, 10.0) # Larger values, approaches normal distribution
]

# Generate table for beta function
beta_table = [(a, b, beta(a, b)) for a, b in ab_values]

# x values to test for betai
x_values = [
    0.1,  # Small x, tests left tail
    0.25, # Lower quartile
    0.5,  # Median
    0.75, # Upper quartile
    0.9   # Large x, tests right tail
]

# Generate table for betai function
betai_table = [(a, b, x, betainc(a, b, x)) for a, b in ab_values for x in x_values]

# p values to test for invbetai
p_values = [
    0.01,  # Very small probability, tests left tail
    0.1,   # Small probability
    0.25,  # Lower quartile
    0.5,   # Median
    0.75,  # Upper quartile
    0.9,   # Large probability
    0.99   # Very large probability, tests right tail
]

# Generate table for invbetai function
invbetai_table = [(a, b, p, betaincinv(a, b, p)) for a, b in ab_values for p in p_values]

# Print tables in Rust code format
print("const BETA_TABLE: [(f64, f64, f64); {}] = [".format(len(beta_table)))
for a, b, y in beta_table:
    print("    ({:.14e}, {:.14e}, {:.14e}),".format(a, b, y))
print("];")

print("\nconst BETAI_TABLE: [(f64, f64, f64, f64); {}] = [".format(len(betai_table)))
for a, b, x, y in betai_table:
    print("    ({:.14e}, {:.14e}, {:.14e}, {:.14e}),".format(a, b, x, y))
print("];")

print("\nconst INVBETAI_TABLE: [(f64, f64, f64, f64); {}] = [".format(len(invbetai_table)))
for a, b, p, x in invbetai_table:
    print("    ({:.14e}, {:.14e}, {:.14e}, {:.14e}),".format(a, b, p, x))
print("];")
