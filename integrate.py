# Grudat26, Emil Risan Tobin, övning 7

def integrate_adaptive(f, a, b, tol, total_val=None):
    """
    Recursively integrates a function f(x) over the interval [a, b] using Simpson's method
    with the function _integrate_simpson.

    Input arguments:
        f (callable): A Python function representing the integrand.
        a: The lower bound of the integration interval, float/str.
        b: The upper bound of the integration interval, float/str.
        tol: The user can choose to increase this variable to increase the answer's accuracy.
    """

    mid = (a + b) / 2

    left_area = _integrate_simpson(f, a, mid)
    right_area = _integrate_simpson(f, mid, b)

    if total_val is None:
        total_val = _integrate_simpson(f, a, b)

    if abs(left_area + right_area - total_val) < 15 * tol:
        return left_area + right_area + (left_area + right_area - total_val) / 15

    return (integrate_adaptive(f, a, mid, tol / 2, left_area) +
            integrate_adaptive(f, mid, b, tol / 2, right_area))


def _integrate_simpson(f, a, b, n=1000):
    """
    Simpson integrator to numerically integrate function f on the interval [a, b].

    Helper function to the functions integrate_adaptive and _prepare_limits.
    """
    # Make sure n is an even number
    if n % 2 > 0:
        n += 1

    # Step length
    h = (b-a)/n

    res = f(a) + f(b)
    for i in range(1, n):
        xi = a + h*i
        if i % 2 == 1:
            res += 4*f(xi)
        else:
            res += 2*f(xi)
    res = (h/3)*res
    return res
