from scipy.special import jv, yn, kn, iv

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