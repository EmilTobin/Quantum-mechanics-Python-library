# Grudat26, Emil Risan Tobin, övning 7

import math
from integrate import integrate_adaptive, _integrate_simpson

def expected_val(psi, operator, a, b, tol=1e-7):
    """
    Calculates the expectation value of a given operator for a wavefunction.
    Used internally by uncertainty_x and uncertainty_p.

    Input arguments:
        psi (callable): A Python function representing the wavefunction.
        operator (str): The type of operator to apply. Supported values: "x", "x2", "p", "p2".
        a, b: The lower and upper bound of the integration interval, float/str.
        tol: The user can choose to increase this variable to increase the answer's accuracy.
    Note:
        To integrate over an infinite range, you can pass float('inf')
        or "-inf". The library will approximate this using a large
        finite boundary.
    """
    # Explain this

    # Hantera undre gränsen a
    if a == '-inf' or a == float('-inf'):
        a, _ = _prepare_limits(psi)

    # Hantera övre gränsen b
    if b == 'inf' or b == float('inf'):
        _, b = _prepare_limits(psi)

    # Explain this
    def integrand(x):
        return _apply_operator(x, psi, operator)

    val = integrate_adaptive(integrand, a, b, tol)
    return val.real

def uncertainty_x(psi, a, b, tol=1e-7):
    """
    Calculates the position uncertainty (Delta x) for a given wavefunction.

    Input arguments:
        psi (callable): A Python function representing the wavefunction.
            Must accept a single float (x) and return a float or complex number.
        a: The lower bound of the integration interval, float/str.
        b: The upper bound of the integration interval, float/str.
        tol: The user can choose to increase this variable to increase the answer's accuracy.
    Note:
        To integrate over an infinite range, you can pass float('inf')
        or "-inf". The library will approximate this using a large
        finite boundary.
    """
    val_1 = expected_val(psi, 'x', a, b, tol)
    val_2 = expected_val(psi, 'x2', a, b, tol)

    return math.sqrt(val_2 - val_1**2)

def uncertainty_p(psi, a, b, tol=1e-7):
    """
    Calculates the momentum uncertainty (Delta p) for a given wavefunction.

    Input arguments:
        psi (callable): A Python function representing the wavefunction.
        a: The lower bound of the integration interval, float/str.
        b: The upper bound of the integration interval, float/str.
        tol: The user can choose to increase this variable to increase the answer's accuracy.
    Note:
        To integrate over an infinite range, you can pass float('inf')
        or "-inf". The library will approximate this using a large
        finite boundary.
    """

    val_1 = expected_val(psi, 'p', a, b, tol)
    val_2 = expected_val(psi, 'p2', a, b, tol)

    return math.sqrt(val_2 - val_1**2)


def _prepare_limits(psi, tol=1e-7):
    """
    Finds the integration boundary [a, b] for a given wavefunction.

    Locates an active region of the probability density |psi(x)|^2 and expands outwards until
    its integral is close to 1.0. Stops expanding when the error is within the specified tolerance.
    """

    def psi2(x):
        # Function that returns the probability density: |psi(x)|^2
        return abs(psi(x)) ** 2

    def find_bounds_recursive(f, a, b, current_area, target=1.0):
        # Base case: The target area has been reached
        if current_area >= target - tol:
            return a, b

        # If boundaries exceed one million, stop to prevent infinite loops.
        if abs(a) > 1e6:
            return a, b

        width = b - a
        new_a = a - width
        new_b = b + width

        # Divide the area and integrate each one
        area_in_new_wings = _integrate_simpson(f, new_a, a, n=10) + \
                            _integrate_simpson(f, b, new_b, n=10)

        # Stops the recursion if the current integral gives a value close to 0.0
        if area_in_new_wings < 1e-15:
            if current_area > 1e-4:
                return a, b

        return find_bounds_recursive(f, new_a, new_b, current_area + area_in_new_wings, target)

    # 1. Finds a starting point where the probability function is active.
    start_x = 0
    for scale in [1, 10, 100, 1000, 10000]:
        test_points = [-scale, scale]
        if any(psi2(p) > 1e-12 for p in test_points):
            start_x = test_points[0] if psi2(test_points[0]) > psi2(test_points[1]) else test_points[1]
            break

    # 2. Start the recursion from a small area around the starting point
    return find_bounds_recursive(psi2, start_x - 0.5, start_x + 0.5,
                                 _integrate_simpson(psi2, start_x - 0.5, start_x + 0.5, n=10))


def _apply_operator(x, psi, operator):
    """
    Helper function to expected_val.
    Applies an operator to a wavefunction and returns the whole integrand:
    psi(x) * [operator*psi(x)].
    For position ('x'), it performs simple multiplication.
    For momentum ('p'), it calculates numerical derivatives using the
    central difference method.

    Args:
        x (float): Coordinate where the operator is applied.
        psi (callable): Wavefunction.
        operator (str): The type of operator, either "x", "x2", "p" or "p2".
    """
    # Plancks reduced constant (J*s)
    h_bar = 1  # 1.054571817e-34

    # l, l2: step sizes for the numerical central difference approximations.
    l = 1e-5
    l2 = 1e-4
    psi_val = psi(x)
    psi_conj = psi_val.conjugate()

    if operator == 'x':
        return psi_conj * x * psi_val

    elif operator == 'x2':
        return psi_conj * (x**2) * psi_val

    elif operator == 'p':
        # Calculate the wavefunction's (psi) derivative using central difference
        derivative1 = (psi(x + l) - psi(x - l)) / (2 * l)

        # Apply momentum-operator: -i * h_bar * d/dx
        p1_psi = -1j * h_bar * derivative1

        return psi_conj * p1_psi

    elif operator == 'p2':
        # Calculate the wavefunction's (psi) 2nd derivative using central difference
        derivative2 = (psi(x + l2) - 2 * psi(x) + psi(x - l2)) / (l2 ** 2)

        # Apply momentum-operator twice: (-i * h_bar)^2 * d^2/dx^2 = -h_bar^2 * d^2/dx^2
        p2_psi = -(h_bar ** 2) * derivative2

        return psi_conj * p2_psi

