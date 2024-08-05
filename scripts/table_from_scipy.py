import numpy as np
from scipy.special import gammaln, lambertw

# x values to test
x_values = [
    0.1, 0.2, 0.5,  # Values between 0 and 1
    1.0, 1.5, 2.0, 2.5,  # Values between 1 and 3
    3.0, 4.0, 5.0,  # Small integers and values between them
    10.0, 20.0, 50.0,  # Larger values
    0.25, 0.75,  # 1/4 and 3/4 points
    1e-5, 1e-10,  # Values very close to 0
    1e5, 1e10  # Very large values
]

# Generate table
table = [(x, gammaln(x)) for x in x_values]

# Print in Rust code format
print("const LN_GAMMA_TABLE: [(f64, f64); {}] = [".format(len(table)))
for x, y in table:
    print("    ({:.16e}, {:.16e}),".format(x, y))
print("];")

#                  Test close to the edge of the domain                                                 Also test massive input
lambert_w0_table = [(-0.367_879_441_171_442_33, -1.0)] + [(x, np.real(lambertw(x))) for x in x_values] + [(1e308, np.real(lambertw(1e308)))]

print("const LAMBERT_W0_TABLE: [(f64, f64); {}] = [".format(len(lambert_w0_table)))
for x, y in lambert_w0_table:
    print("    ({:.16e}, {:.16e}),".format(x, y))
print("];")