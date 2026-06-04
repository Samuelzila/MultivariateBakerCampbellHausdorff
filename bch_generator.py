from sage.all import SR, MatrixSpace, FreeAlgebra, I, factorial
from sage.categories.lie_algebras import LieAlgebras
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra

n = 2

g = LieAlgebra(SR, "H", n)
g = LieAlgebra(SR, 'X,Y')
H = g.gens()


def magic_iterator(n,k):
    if n <= 0 or k < n:
        return

    if n == 1:
        yield (k,)
        return

    # Leave at least 1 for each remaining position
    for x in range(1, k - n + 2):
        for rest in magic_iterator(n - 1, k - x):
            yield (x,) + rest


def gamma(k, indices):
    product = SR(1)
    for q in range(1, n+1):
        product *= factorial(sum(SR(1) for j in indices if j == q))
    return SR(1) / product


def F(k):
    p = [1] * (k+1)
    total = 0
    while p[0] <= n:

        # Compute the commutator term
        commutator = g.bracket(H[p[1]-1], H[p[0]-1]) if k > 0 else H[p[0]-1]
        for j in p[2:]:
            commutator = g.bracket(H[j-1], commutator)

        # Compute gamma
        γ = gamma(k, p[1:])

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

def omega(k):
    term1 = (SR(1) / k) * F(k-1)

    commutator_total = g.zero()
    for r in magic_iterator(n, k):
        commutator = g.bracket(omega(r[1]), omega(r[0])) if len(r) > 1 else omega(r[0])
        print(f"Initial commutator for r={r}: {commutator}")
        for j in r[2:]:
            commutator = g.bracket(omega(j), commutator)

        commutator_total += commutator
        print(f"Updated commutator_total after r={r}: {commutator_total}")
    term2 = - (SR(1) / k)
    term3 = g.zero()
    for j in range(1, n):
        term3 += (SR(1)/factorial(SR(n+1)))*commutator_total

    print(f"term1: {term1}, term2: {term2}, term3: {term3}")

    return term1 + (term2*term3)


print(omega(3))

