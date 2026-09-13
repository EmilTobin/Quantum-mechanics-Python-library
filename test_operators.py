# Grudat26, Emil Risan Tobin, övning 7
"""
python3 -m unittest test_operators.py
"""

import unittest
import math
from operators import expected_val, uncertainty_x, uncertainty_p

class TestQuantumOperators(unittest.TestCase):

    def setUp(self):
        """
        Defines two wavefunctions with known expectations values and uncertainties.
        """
        # 1) Gaussian wavefunction:
        # Tests: Real-valued case, simple integration
        # <x> = 0, <x^2> = 0.5, dx = sqrt(0.5)
        # <p> = 0, <p^2> = 0.5, dp = sqrt(0.5)
        self.psi_gauss = lambda x: (1.0 / math.pi)**0.25 * math.exp(-x**2 / 2.0)

        # 2) Shifted Gaussian curve with momentum
        # Tests: Complex-valued, translated, rapidly oscillating
        # <x> = 2.0, <x2> = 4.5, dx = sqrt(0.5)
        # <p> = 3.0, <p2> = 9.5, dp = sqrt(0.5
        x0 = 2.0
        p0 = 3.0
        self.psi_shifted = lambda x: self.psi_gauss(x - x0) * math.e**(1j * p0 * x)

    def test_expected_val_inf_bounds(self):
        # Function 1)
        # Operator x
        res = expected_val(self.psi_gauss, 'x', '-inf', 'inf')
        self.assertAlmostEqual(res, 0.0, delta=1e-7)

        # Operator x^2
        res = expected_val(self.psi_gauss, 'x2', '-inf', 'inf')
        self.assertAlmostEqual(res, 0.5, delta=1e-7)

        # Operator p
        res_p = expected_val(self.psi_gauss, 'p', '-inf', 'inf')
        self.assertAlmostEqual(res_p, 0.0, delta=1e-7)

        # Operator p^2
        res_p2 = expected_val(self.psi_gauss, 'p2', '-inf', 'inf')
        self.assertAlmostEqual(res_p2, 0.5, delta=1e-7)


        # Function 2)
        # Operator x
        res = expected_val(self.psi_shifted, 'x', '-inf', 'inf')
        self.assertAlmostEqual(res, 2.0, delta=1e-7)

        # Operator x^2
        res = expected_val(self.psi_shifted, 'x2', '-inf', 'inf')
        self.assertAlmostEqual(res, 4.5, delta=1e-7)

        # Operator p
        res_p = expected_val(self.psi_shifted, 'p', '-inf', 'inf')
        self.assertAlmostEqual(res_p, 3.0, delta=1e-7)

        # Operator p^2
        res_p2 = expected_val(self.psi_shifted, 'p2', '-inf', 'inf')
        self.assertAlmostEqual(res_p2, 9.5, delta=1e-7)

    def test_uncertainty_inf_bounds(self):
        # Function 1)
        # Operator x
        expected_dx = math.sqrt(0.5)
        res = uncertainty_x(self.psi_gauss, '-inf', 'inf')
        self.assertAlmostEqual(res, expected_dx, delta=1e-7)

        # Operator p
        expected_dp = math.sqrt(0.5)
        res_dp = uncertainty_p(self.psi_gauss, '-inf', 'inf')
        self.assertAlmostEqual(res_dp, expected_dp, delta=1e-7)


        # Function 2)
        # Operator x
        expected_dx = math.sqrt(0.5)
        res = uncertainty_x(self.psi_shifted, '-inf', 'inf')
        self.assertAlmostEqual(res, expected_dx, delta=1e-7)

        # Operator p
        expected_dp = math.sqrt(0.5)
        res_dp = uncertainty_p(self.psi_shifted, '-inf', 'inf')
        self.assertAlmostEqual(res_dp, expected_dp, delta=1e-7)

    def test_finite_bounds(self):

        # Function 1)
        # Operator x
        res_x = expected_val(self.psi_gauss, 'x', -2, 2, tol=1e-8)
        self.assertAlmostEqual(res_x, 0.0, delta=1e-7)

        # Operator x^2
        res_x2 = expected_val(self.psi_gauss, 'x2', -2, 2, tol=1e-8)
        self.assertAlmostEqual(res_x2, 0.47699, delta=1e-5)

        # Operator p
        res_p = expected_val(self.psi_gauss, 'p', -2, 2, tol=1e-8)
        self.assertAlmostEqual(res_p, 0.0, delta=1e-7)

        # Operator p^2
        res_p2 = expected_val(self.psi_gauss, 'p2', -2, 2, tol=1e-8)
        self.assertAlmostEqual(res_p2, 0.51833, delta=1e-5)

    def test_uncertainty_finite_bounds(self):
        """Test that uncertainty calculations work correctly within finite intervals"""
        # Function 1)
        # Operator x
        res_dx = uncertainty_x(self.psi_gauss, -2, 2, tol=1e-9)
        expected_dx = 0.69064
        self.assertAlmostEqual(res_dx, expected_dx, delta=1e-3)

        # Operator p
        res_dp = uncertainty_p(self.psi_gauss, -2, 2, tol=1e-8)
        expected_dp = 0.71995
        self.assertAlmostEqual(res_dp, expected_dp, delta=1e-5)


if __name__ == '__main__':
    unittest.main()

