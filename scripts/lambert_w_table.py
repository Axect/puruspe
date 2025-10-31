import numpy as np
from scipy.special import lambertw

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

#                                                                   Also test massive input
lambert_w0_table = [(x, np.real(lambertw(x))) for x in x_values] + [(1e308, np.real(lambertw(1e308)))]

print("const LAMBERT_W0_TABLE: [(f64, f64); {}] = [".format(len(lambert_w0_table)))
for x, y in lambert_w0_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

lambert_wm1_x_values = [
    # Random numbers between -1/e and 0
    -0.1623304668493971, -0.0141318890794931, -0.1781317247582059,
    -0.1038435920310578, -0.364876111085733 , -0.0546821782149597,
    -0.3632700388723338, -0.010217109406601 , -0.2376990995256563,
    -0.3037676546798408,
    # Small values
    -1e-3, -3.1e-5, -1e-100,
]

lambert_wm1_table = [(x, np.real(lambertw(x, -1))) for x in lambert_wm1_x_values]

print("const LAMBERT_WM1_TABLE: [(f64, f64); {}] = [".format(len(lambert_wm1_table)))
for x, y in lambert_wm1_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

# We hard code the result here because scipy.lambertw gives nan for the branch point,
# while we define the value of the principal and secondary branch at the branch point
# to be -1. This aligns with Wolfram Mathematica as well as Wikipedia.
print("const BRANCH_POINT: (f64, f64) = (-0.36787944117144232, -1.0);")

z_values = [
    # Some complex numbers
    0.1 + 0.2j, -1.5 - 0.5j, 1.0 + 1.0j, -2.0 + 2.0j,
    # Small complex numbers
    1e-3 + 1e-3j, -1e-5 - 1e-5j, 1e-10 + 1e-10j,
    # Large complex numbers
    1e5 + 1e5j, -1e10 - 1e10j
]

k_values = [
    0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5,
    10, -10, 20, -20, 50, -50
]

print("const COMPLEX_LAMBERT_W_TABLE: [(i32, (f64, f64), (f64, f64)); {}] = [".format(len(k_values) * len(z_values)))
for k in k_values:
    for z in z_values:
        w = lambertw(z, k)
        print("    ({}, ({:.14e}, {:.14e}), ({:.14e}, {:.14e})),".format(k, np.real(z), np.imag(z), np.real(w), np.imag(w)))
print("];")
