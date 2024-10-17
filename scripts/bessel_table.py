from scipy.special import yn, kn, jv, iv, yv, kv

# x values to test
x_values = [
    0.2, 0.5,  # Values between 0 and 1
    1.0, 1.5, 2.0, 2.5,  # Values between 1 and 3
    3.0, 4.0, 5.0,  # Small integers and values between them
    0.25, 0.75,  # 1/4 and 3/4 points
    10.0, 20.0 # Larger values
]

j0_table = [(x, jv(0, x)) for x in x_values]

print("const J0_TABLE: [(f64, f64); {}] = [".format(len(j0_table)))
for x, y in j0_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

j5_table = [(x, jv(5, x)) for x in x_values]

print("const J5_TABLE: [(f64, f64); {}] = [".format(len(j5_table)))
for x, y in j5_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

j20_table = [(x, jv(20, x)) for x in x_values]

print("const J20_TABLE: [(f64, f64); {}] = [".format(len(j20_table)))
for x, y in j20_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

y0_table = [(x, yn(0, x)) for x in x_values]

print("const Y0_TABLE: [(f64, f64); {}] = [".format(len(y0_table)))
for x, y in y0_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

y5_table = [(x, yn(5, x)) for x in x_values]

print("const Y5_TABLE: [(f64, f64); {}] = [".format(len(y5_table)))
for x, y in y5_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

k0_table = [(x, kn(0, x)) for x in x_values]

print("const K0_TABLE: [(f64, f64); {}] = [".format(len(k0_table)))
for x, y in k0_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

k5_table = [(x, kn(5, x)) for x in x_values]

print("const K5_TABLE: [(f64, f64); {}] = [".format(len(k5_table)))
for x, y in k5_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

i0_table = [(x, iv(0, x)) for x in x_values]

print("const I0_TABLE: [(f64, f64); {}] = [".format(len(i0_table)))
for x, y in i0_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

i5_table = [(x, iv(5, x)) for x in x_values]

print("const I5_TABLE: [(f64, f64); {}] = [".format(len(i5_table)))
for x, y in i5_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")


# nu values to test
nu_values = [0.5, 1.5, 2.5, 5.5, 10.5]
# Explanation: We choose these nu values to test:
# - Fractional orders (0.5, 1.5, 2.5) to check behavior between integer orders
# - A medium-large value (5.5) to test performance for larger orders
# - A large value (10.5) to test asymptotic behavior for high orders

# x values to test
x_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
# Explanation: We choose these x values to test:
# - Small values (0.1, 0.5) to check behavior near x = 0
# - Medium values (1.0, 2.0, 5.0) to test typical use cases
# - Large values (10.0, 20.0, 50.0) to test asymptotic behavior for large x
# This range covers most practical applications and helps identify any issues 
# in different regimes of the functions.

print("\nconst JNU_YNU_TABLE: [(f64, f64, f64, f64); {}] = [".format(len(nu_values) * len(x_values)))
for nu in nu_values:
    for x in x_values:
        jnu = jv(nu, x)
        ynu = yv(nu, x)
        print(f"    ({nu:.1f}, {x:.1f}, {jnu:.14e}, {ynu:.14e}),")
print("];")

print("\nconst INU_KNU_TABLE: [(f64, f64, f64, f64); {}] = [".format(len(nu_values) * len(x_values)))
for nu in nu_values:
    for x in x_values:
        inu = iv(nu, x)
        knu = kv(nu, x)
        print(f"    ({nu:.1f}, {x:.1f}, {inu:.14e}, {knu:.14e}),")
print("];")

print("\nconst BESSELJY_TABLE: [(f64, f64, f64, f64, f64, f64); {}] = [".format(len(nu_values) * len(x_values)))
for nu in nu_values:
    for x in x_values:
        j = jv(nu, x)
        y = yv(nu, x)
        jp = jv(nu-1, x) - nu/x * jv(nu, x)  # Derivative of J_nu(x)
        yp = yv(nu-1, x) - nu/x * yv(nu, x)  # Derivative of Y_nu(x)
        print(f"    ({nu:.1f}, {x:.1f}, {j:.14e}, {y:.14e}, {jp:.14e}, {yp:.14e}),")
print("];")

print("\nconst BESSELIK_TABLE: [(f64, f64, f64, f64, f64, f64); {}] = [".format(len(nu_values) * len(x_values)))
for nu in nu_values:
    for x in x_values:
        i = iv(nu, x)
        k = kv(nu, x)
        ip = iv(nu-1, x) - nu/x * iv(nu, x)  # Derivative of I_nu(x)
        kp = -kv(nu-1, x) - nu/x * kv(nu, x)  # Derivative of K_nu(x)
        print(f"    ({nu:.1f}, {x:.1f}, {i:.14e}, {k:.14e}, {ip:.14e}, {kp:.14e}),")
print("];")
