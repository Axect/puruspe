from scipy.special import gammaln, gamma
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
    print("    ({:.16e}, {:.16e}),".format(x, y))
print("];")

# ┌──────────────────────────────────────────────────────────┐
#  gamma
# └──────────────────────────────────────────────────────────┘
# x values to test (excluding inputs that would result in complex outputs)
gamma_x_values = [
    -1.5, -0.5,  # Negative non-integer values that result in real outputs
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
