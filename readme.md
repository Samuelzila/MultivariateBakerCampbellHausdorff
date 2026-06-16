# Multivariate Baker–Campbell–Hausdorff

A SageMath implementation of the multivariate Baker–Campbell–Hausdorff (BCH) expansion.

This package was developed primarily for research and experimentation with multivariate BCH formulas in free Lie algebras.

## Requirements

* SageMath (tested with Sage 10.9)
* Python 3 (tested with Python 3.14)


## Quick Start

You can start a jupyter notebook environment with the following command:

```bash
sage -n jupyter
```

Then, you can import the `MultivariateBakerCampbellHausdorff` class and compute the BCH terms as follows:

```python
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from bch import MultivariateBakerCampbellHausdorff

g = LieAlgebra(SR, "H", 2)

bch = MultivariateBakerCampbellHausdorff(g, g.gens())

bch(3) # Returns the BCH expansion up to degree 3

```
Alternatively, you can compute a specific term of the BCH expansion directly:

```python
term = bch.nth(3) # Returns the third order term of the BCH expansion
print(term)
```

## Acknowledgements

Thank you to [Maxime Dion](https://github.com/oneminimax) for providing the mathematical formulas on which this implementation is based.
