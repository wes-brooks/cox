# Task list for cox package:

- Correct calculations for derivative of log determinant of Cholesky factor
- Include the entropy term in variational approximation and its score functions (will require LogDetDerChol working)
- Use the corrected LogDetDerChol in calculation of the Laplace approximation (write a gradient function for Laplace optim)
- Complete the documentation, write vignette
- Create unit tests
- Clean out cruft code
