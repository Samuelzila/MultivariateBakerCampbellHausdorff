from sage.all import SR, MatrixSpace, FreeAlgebra, I, factorial
from sage.categories.lie_algebras import LieAlgebras
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra

from bch.combinatorics import gamma_coefficient, compositions

from functools import cache

class MultivariateBakerCampbellHausdorff:
    def __init__(self, lie_algebra, generators):
        """
        Initializes the BCH generator for a given Lie algebra and its generators.
            - `lie_algebra`: The Lie algebra for which we want to compute the BCH formula
            - `generators`: A tuple of lie algebra elements that will be used as the input to the BCH formula.
            The order of the generators in the tuple determines the order of the inputs to the BCH formula, from left to right.

        Usage:

            ```
            g = LieAlgebra(SR, "H", n)
            H = g.gens()
            bch_gen = MultivariateBakerCampbellHausdorff(g, H)
            result = bch_gen(k)  # Compute the BCH formula up to order k
            ```
        """
        self.lie_algebra = lie_algebra
        self.H = tuple(generators)
        self.n = len(generators)


    def __call__(self, k):
        """
        Computes the BCH formula up to order `k` by summing the nth terms from 1 to `k`.
        """
        return sum(self.omega(i) for i in range(1, k+1))

    @cache
    def F(self, k):
        g = self.lie_algebra
        H = self.H
        n = self.n

        p = [1] * (k+1)
        total = 0
        while p[0] <= self.n:

            # Compute the commutator term
            commutator = g.bracket(H[p[1]-1], H[p[0]-1]) if len(p) > 1 else H[p[0]-1]
            for j in p[2:]:
                commutator = g.bracket(H[j-1], commutator)

            # Compute gamma
            γ = gamma_coefficient(n, p[1:])

            # Add the term to the total
            total += γ * commutator

            # Ensure we stay within the simplex
            p[-1] += 1
            for j in reversed(range(1, len(p))):
                if p[j] > p[j-1]:
                    p[j] = 1
                    p[j-1] += 1
                else:
                    break

        return total


    @cache
    def omega(self, k):
        """
        Computes the `k`-th order term of the multivariate Baker-Campbell-Hausdorff formula.
        """
        g = self.lie_algebra
        F = self.F
        omega = self.omega

        term1 = (SR(1) / k) * F(k-1)

        term2 = (SR(1) / k)

        term3 = g.zero()
        for N in range(1, k):
            commutator_total = g.zero()
            for r in compositions(N+1, k):
                commutator = r[0]*omega(r[0])
                for j in reversed(range(1, len(r))):
                    commutator = g.bracket(omega(r[j]), commutator)

                commutator_total += commutator

            term3 += (SR(1)/factorial(SR(N+1)))*commutator_total

        return term1 - (term2*term3)

    def nth(self, n):
        """
        Computes the `n`-th order term of the BCH formula.
        Alias for `omega(n)`.
        """
        return self.omega(n)
