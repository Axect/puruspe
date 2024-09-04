from scipy.special import erf, erfc, erfinv, erfcinv

# x values to test for erf and erfc
x_values = [
    -3.0,  # Large negative value
    -1.5,  # Moderate negative value
    -0.5,  # Small negative value
    0.0,   # Zero
    0.5,   # Small positive value
    1.5,   # Moderate positive value
    3.0,   # Large positive value
    -1e-10, # Very small negative value
    1e-10,  # Very small positive value
    -10.0,  # Very large negative value
    10.0    # Very large positive value
]

# Generate tables for erf and erfc functions
erf_table = [(x, erf(x)) for x in x_values]
erfc_table = [(x, erfc(x)) for x in x_values]

# p values to test for inverf and inverfc
p_values = [
    -0.99,  # Near -1 (for inverf only)
    -0.5,   # Negative value (for inverf only)
    0.0,    # Zero
    0.5,    # Moderate value
    0.99,   # Near 1 (for inverf) / Near 0 (for inverfc)
    1e-10,  # Very small value (for inverfc only)
    1.0,    # 1 (for inverfc only)
    1.99    # Near 2 (for inverfc only)
]

# Generate tables for inverf and inverfc functions
inverf_table = [(p, erfinv(p)) for p in p_values if -1 < p < 1]
inverfc_table = [(p, erfcinv(p)) for p in p_values if 0 < p < 2]

# Print tables in Rust code format
print("const ERF_TABLE: [(f64, f64); {}] = [".format(len(erf_table)))
for x, y in erf_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

print("\nconst ERFC_TABLE: [(f64, f64); {}] = [".format(len(erfc_table)))
for x, y in erfc_table:
    print("    ({:.14e}, {:.14e}),".format(x, y))
print("];")

print("\nconst INVERF_TABLE: [(f64, f64); {}] = [".format(len(inverf_table)))
for p, x in inverf_table:
    print("    ({:.14e}, {:.14e}),".format(p, x))
print("];")

print("\nconst INVERFC_TABLE: [(f64, f64); {}] = [".format(len(inverfc_table)))
for p, x in inverfc_table:
    print("    ({:.14e}, {:.14e}),".format(p, x))
print("];")
