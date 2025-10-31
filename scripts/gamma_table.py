from scipy.special import gammaln, gamma, gammainc, gammaincc, gammaincinv
import numpy as np

# ┌──────────────────────────────────────────────────────────┐
#  ln_gamma
# └──────────────────────────────────────────────────────────┘
# x values to test
ln_gamma_x_values = [
    0.1, 0.2, 0.5,  # Values between 0 and 1
    1.0, 1.5, 2.0, 2.5,  # Values between 1 and 3
    3.0, 4.0, 5.0,  # Small integers and values between them
    10.0, 20.0, 50.0,  # Larger values
    0.25, 0.75,  # 1/4 and 3/4 points
    1e-5, 1e-10,  # Values very close to 0
    1e5, 1e10  # Very large values
]

# Generate table
table = [(x, gammaln(x)) for x in ln_gamma_x_values]

# Print in Rust code format
print("const LN_GAMMA_TABLE: [(f64, f64); {}] = [".format(len(table)))
for x, y in table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

# ┌──────────────────────────────────────────────────────────┐
#  gamma
# └──────────────────────────────────────────────────────────┘
# x values to test (excluding inputs that would result in complex outputs)
gamma_x_values = [
    -1.5, -0.5, -100.5,  # Negative non-integer values that result in real outputs
    -4503599627370495.5,  # Large negative non-integer value that underflows to zero.
    0.1, 0.2, 0.5,  # Values between 0 and 1
    1.0, 1.5, 2.0, 2.5,  # Values between 1 and 3
    3.0, 4.0, 5.0,  # Small integers and values between them
    10.0, 20.0, 50.0,  # Larger values
    0.25, 0.75,  # 1/4 and 3/4 points
    1e-5, 1e-10,  # Values very close to 0
    100.0, 170.0  # Larger values (but not too large to avoid overflow)
]

# Generate table
table = [(x, gamma(x)) for x in gamma_x_values if not np.iscomplex(gamma(x))]

# Print in Rust code format
print("const GAMMA_TABLE: [(f64, f64); {}] = [".format(len(table)))
for x, y in table:
    print(f"    ({x:.16e}, {y:.16e}),")
print("];")

# a values to test for gammp, gammq, and invgammp
a_values = [
    0.1,  # Very small shape parameter, tests behavior near zero
    0.5,  # Half-integer shape parameter, important in statistical applications
    1.0,  # Exponential distribution case
    2.0,  # Chi-squared distribution with 4 degrees of freedom
    5.0,  # Medium-sized shape parameter
    10.0  # Larger shape parameter, approaches normal distribution
]

# x values to test for gammp and gammq
x_values = [
    0.1,  # Very small x, tests behavior near zero
    0.5,  # Tests behavior for x < 1
    1.0,  # Important reference point, x = a for some cases
    2.0,  # Tests behavior for x > 1
    5.0,  # Larger x value
    10.0, # Much larger x value
    20.0  # Very large x value, tests tail behavior
]

# Generate tables for gammp and gammq
gammp_table = [(a, x, gammainc(a, x)) for a in a_values for x in x_values]
gammq_table = [(a, x, gammaincc(a, x)) for a in a_values for x in x_values]

# p values to test for invgammp
p_values = [
    0.01,  # Very small probability, tests left tail
    0.1,   # Small probability
    0.25,  # Lower quartile
    0.5,   # Median
    0.75,  # Upper quartile
    0.9,   # Large probability
    0.99   # Very large probability, tests right tail
]

# Generate table for invgammp
invgammp_table = [(a, p, gammaincinv(a, p)) for a in a_values for p in p_values]

# Print tables in Rust code format
print("const GAMMP_TABLE: [(f64, f64, f64); {}] = [".format(len(gammp_table)))
for a, x, y in gammp_table:
    print("    ({:.14e}, {:.14e}, {:.14e}),".format(a, x, y))
print("];")

print("\nconst GAMMQ_TABLE: [(f64, f64, f64); {}] = [".format(len(gammq_table)))
for a, x, y in gammq_table:
    print("    ({:.14e}, {:.14e}, {:.14e}),".format(a, x, y))
print("];")

print("\nconst INVGAMMP_TABLE: [(f64, f64, f64); {}] = [".format(len(invgammp_table)))
for a, p, x in invgammp_table:
    print("    ({:.14e}, {:.14e}, {:.14e}),".format(a, p, x))
print("];")
