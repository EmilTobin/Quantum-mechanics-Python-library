# Numerical quantum mechanics 

A Python library for numerical integration and computation of quantum mechanical expectation values and uncertainties.

## Description

This library numerically integrates expectation values and uncertainties using an adaptive Simpsons integrator and can handle the most common wavefunctions defined by the user. 

## Functions
* expected_val: Numerically calculates the expectation value given a wavefunction. Works for x, x^2, p and p^2 operators. 
* uncertainty_x: Calculates the uncertainty of the position-operator (x) given a wavefunction.  
* uncertainty_p: Calculates the uncertainty of the position-operator (p) given a wavefunction.  

## Example code
Example code showcasing the librarie's functions is shown in the file example_code.py

## Future development
This library is currently considered feature-complete. No active development or future releases are planned. The only accepted reason to modify the API of this package is to handle issues that can't be resolved in any other reasonable way.
