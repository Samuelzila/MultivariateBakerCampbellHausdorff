# from sage.all import *
#
# # Number of exponentials and truncation order
# n = 2
# N = 2
#
# # Free associative algebra
# names = [f'X{i+1}' for i in range(n)]
# A = FreeAlgebra(SR, names)
# X = A.gens()
#
# def truncate_degree(expr, N):
#     A = expr.parent()
#     out = A.zero()
#
#     for m, c in expr.monomial_coefficients().items():
#         if len(m) <= N:
#             out += c * A.monomial(m)
#
#     return out
#
# def exp_truncated(x, N):
#     out = A.one()
#     p = A.one()
#
#     for k in range(1, N+1):
#         p *= x
#         out += p/factorial(k)
#
#     return truncate_degree(out, N)
#
# def log_truncated(u, N):
#     """
#     log(u) assuming u = 1 + higher order terms.
#     """
#     a = truncate_degree(u - 1, N)
#
#     out = A.zero()
#     p = a
#
#     for k in range(1, N+1):
#         out += (-1)**(k+1) * p / k
#         p = truncate_degree(p*a, N)
#
#     return truncate_degree(out, N)
#
# # Product of exponentials
# U = A.one()
#
# for Xi in X:
#     U = truncate_degree(U * exp_truncated(Xi, N), N)
#
# Omega = log_truncated(U, N)
#
# print(Omega)

from sage.all import *

##############################################################################
# PARAMETERS
##############################################################################

n = 3
N = 3

##############################################################################
# FREE LIE ALGEBRA
##############################################################################

L = LieAlgebra(SR, ["H_T", "H_S"])

HT, HS = L.gens()
gamma, beta, epsilon = var('gamma beta epsilon')
X = [gamma*HT, beta*HS, epsilon*HT]
t = var('t')
X = [I*Xi for Xi in X]  # Make them imaginary for better readability

##############################################################################
# COMMUTATOR
##############################################################################

def ad(a, b):
    return a.bracket(b)

##############################################################################
# BASIC MAGNUS TERMS (CORRECT STRUCTURE)
##############################################################################

def omega_1():
    return sum(X)

def omega_2():
    out = L.zero()
    for i in range(n):
        for j in range(i):
            out += ad(X[i], X[j])
    return SR(1)/2 * out

def omega_3():
    out = L.zero()

    for i in range(n):
        for j in range(i):
            for k in range(j):

                xi, xj, xk = X[i], X[j], X[k]

                out += ad(xi, ad(xj, xk))
                out += ad(ad(xi, xj), xk)

    return SR(1)/12 * out

def omega_4():
    out = L.zero()

    for i in range(n):
        for j in range(i):
            for k in range(j):
                for l in range(k):

                    xi, xj, xk, xl = X[i], X[j], X[k], X[l]

                    out += ad(xi, ad(xj, ad(xk, xl)))
                    out += ad(ad(ad(xi, xj), xk), xl)

    return SR(1)/24 * out

##############################################################################
# FULL EXPANSION
##############################################################################

Omega = L.zero()
Omega += omega_1()

if N >= 2:
    Omega += omega_2()
if N >= 3:
    Omega += omega_3()
if N >= 4:
    Omega += omega_4()

print(Omega)

