from sage.all import factorial, SR


def compositions(size,total):
    """
    Generate all compositions of `total` into `size` parts, where each part is a positive integer.
    """
    if size <= 0 or total < size:
        return

    if size == 1:
        yield (total,)
        return

    # Leave at least 1 for each remaining position
    for x in range(1, total - size + 2):
        for rest in compositions(size - 1, total - x):
            yield (x,) + rest


def gamma_coefficient(n, indices):
    """
    Computes the gamma coefficients for the commutator terms.
    `n` is the total number of generators, and `indices` is a list of the indices
    of the generators involved in the commutator.
    """
    product = SR(1)
    for q in range(1, n+1):
        product *= factorial(sum(SR(1) for j in indices if j == q))
    return SR(1) / product
