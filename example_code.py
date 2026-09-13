import math
from operators import expected_val, uncertainty_x, uncertainty_p

# User define their wavefunctions as callable Python-functions
def psi_gauss(x):
    # Normalized Gaussian wavefunction
    return (1.0 / math.pi) ** 0.25 * math.exp(-x ** 2 / 2.0)


# 1) Calculate Expectation Values
print("1) Calculating Expectation Values:")

# Position operators
mean_x = expected_val(psi_gauss, 'x', '-inf', 'inf')
mean_x2 = expected_val(psi_gauss, 'x2', '-inf', 'inf')

# Momentum operators (using numerical central differences internally)
mean_p = expected_val(psi_gauss, 'p', '-inf', 'inf')
mean_p2 = expected_val(psi_gauss, 'p2', '-inf', 'inf')

print(f"<x> = {mean_x.real:.6f}")
print(f"<x^2> = {mean_x2.real:.6f}")
print(f"<p> = {mean_p.real:.6f}")
print(f"<p^2> = {mean_p2.real:.6f}\n")

# 2) Evaluate Heisenberg's Uncertainty Principle
print("2) Evaluating Heisenberg's Uncertainty Principle:")

# Calculate individual standard deviations (Uncertainties)
delta_x = uncertainty_x(psi_gauss, '-inf', 'inf')
delta_p = uncertainty_p(psi_gauss, '-inf', 'inf')

uncertainty_product = delta_x * delta_p
theoretical_minimum = 0.5  # h_bar / 2 where h_bar = 1

print(f"Delta x (Position Uncertainty) : {delta_x:.6f}")
print(f"Delta p (Momentum Uncertainty) : {delta_p:.6f}")
print(f"Product (Delta x * Delta p)     : {uncertainty_product:.6f}")
print(f"Theoretical Minimum (hbar/2)    : {theoretical_minimum:.1f}")

if uncertainty_product >= (theoretical_minimum - 1e-5):
    print("\n[SUCCESS] Heisenberg's Uncertainty Principle holds true! (dx * dp >= 0.5)")
else:
    print("\n[WARNING] Uncertainty principle violated. Check operator precision settings.")
