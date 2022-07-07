from scipy import sparse
import scipy as sp
import numpy as np
import cvxpy
import ecos

min
    q @ x
s.t.
    F @ x = 1
    A1 @ x <= 0
    x @ A2 @ x + A3 @ x + A4 - h^2 <= 0

# -----------------------------------------
# Version 1 (Fails because h not right dim)
# -----------------------------------------
q   = np.array([1, 1, 1, 1, 0, 0, 0, 0])
F   = np.array([0, 0, 0, 0, 1, 1, 1, 1])
A1  = np.array([
     [-1,   0,   0,   0,   1,   0,   0,   0],
     [ 0,  -1,   0,   0,   1,   1,   0,   0],
     [ 0,   0,  -1,   0,   1,   1,   1,   0],
     [ 0,   0,   0,  -1,   1,   1,   1,   1],
     [-1,   0,   0,   0,  -1,   0,   0,   0],
     [ 0,  -1,   0,   0,  -1,  -1,   0,   0],
     [ 0,   0,  -1,   0,  -1,  -1,  -1,   0],
     [ 0,   0,   0,  -1,  -1,  -1,  -1,  -1]
])
A2 = np.array([
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,  .0005316404, -.0000287359, -.0001288079, -8.12950e-06],
    [0, 0, 0, 0, -.0000287359,   .000368052, -.0001560354,  .0000317829],
    [0, 0, 0, 0, -.0001288079, -.0001560354,  .0004084173, -.0000971909],
    [0, 0, 0, 0, -8.12950e-06,  .0000317829, -.0000971909,  .0003089595]
])
A3 = np.array([0, 0, 0, 0, -.0001051178, .0000722476, .0001285712, .0002395732])
hh = .02599883
A4 = .0003599704
postPeriodIndices = [4, 5, 6, 7]
threshold_sumweights = 1

E0   = A2[postPeriodIndices,:][:,postPeriodIndices]
L, C = np.linalg.eigh(E0)
E1   = (C * np.sqrt(L)) @ C.T
E2   = (C / np.sqrt(L)) @ C.T
E3   = (C / L) @ C.T
D1   = A2.copy()
D2   = A2.copy()
D3   = A2.copy()
D1[postPeriodIndices,:][:,postPeriodIndices] = E1
D2[postPeriodIndices,:][:,postPeriodIndices] = E2
D3[postPeriodIndices,:][:,postPeriodIndices] = E3
assert np.abs((C @ np.diag(L) @ C.T) - E0).max() < 1e-12
assert np.abs(E1 @ E1 - E0).max() < 1e-12
assert np.abs(E1 @ E2 - np.eye(E0.shape[0])).max() < 1e-12
assert np.abs(E3 @ E0 - np.eye(E0.shape[0])).max() < 1e-12

n = 8
x = [1.86114e-08, .03646372, .3908837, 1, 1.82187e-08, .0364637, .35442, .6091163]
x = np.concatenate((x, A3))
A = np.vstack([
    np.concatenate((F, np.zeros(n)))[None, :],
    np.hstack((np.zeros((n, n)), np.eye(n)))
])
G1 = np.hstack([A2, np.zeros((n, n))])
G2 = np.hstack([D1, 0.5 * D2])
G  = sparse.csr_matrix(np.vstack([G1, G2]))
A  = sparse.csr_matrix(A)
h  = np.array([np.sqrt(0.25 * (A3 @ D3 @ A3) + hh ** 2 - A4)])
c  = np.concatenate([q.astype(float), np.zeros(n)])
b  = np.concatenate([[threshold_sumweights], A3])

solution = ecos.solve(c, G, h, {'l': 8, 'q': 8, 'e': []}, A, b[:,None])

# ------------------------------------------------------------------
# Version 2 (xx you are here; treat h (paper not) as b (wiki not)...
# ------------------------------------------------------------------

# Wait, wait, so the requirement is that h - Gx has to be part of
# the second-order cone. A vector u = (u_1, u_2) in R^1 x R^{n-1}
# is part of the second orer cone if ||u_2||_2 <= u_1. Therefore in
# wikipedia notation, G^T = [c A^T] and h^T = [d b^T] for the second
# order cone. For the first-order cone, it's simpler, since you just
# need Gx <= h.

q   = np.array([1, 1, 1, 1, 0, 0, 0, 0])
F   = np.array([0, 0, 0, 0, 1, 1, 1, 1])
A1  = np.array([
     [-1,   0,   0,   0,   1,   0,   0,   0],
     [ 0,  -1,   0,   0,   1,   1,   0,   0],
     [ 0,   0,  -1,   0,   1,   1,   1,   0],
     [ 0,   0,   0,  -1,   1,   1,   1,   1],
     [-1,   0,   0,   0,  -1,   0,   0,   0],
     [ 0,  -1,   0,   0,  -1,  -1,   0,   0],
     [ 0,   0,  -1,   0,  -1,  -1,  -1,   0],
     [ 0,   0,   0,  -1,  -1,  -1,  -1,  -1]
])
A2 = np.array([
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,            0,            0,            0,            0],
    [0, 0, 0, 0,  .0005316404, -.0000287359, -.0001288079, -8.12950e-06],
    [0, 0, 0, 0, -.0000287359,   .000368052, -.0001560354,  .0000317829],
    [0, 0, 0, 0, -.0001288079, -.0001560354,  .0004084173, -.0000971909],
    [0, 0, 0, 0, -8.12950e-06,  .0000317829, -.0000971909,  .0003089595]
])
A3 = np.array([0, 0, 0, 0, -.0001051178, .0000722476, .0001285712, .0002395732])
hh = .02599883
hh = 0.02927835
hh = .0215963422
A4 = .0003599704
postPeriodIndices = [4, 5, 6, 7]
threshold_sumweights = 1

E0   = A2[np.ix_(postPeriodIndices, postPeriodIndices)]
L, C = np.linalg.eigh(E0)
E1   = (C * np.sqrt(L)) @ C.T
E2   = (C / np.sqrt(L)) @ C.T
E3   = (C / L) @ C.T
D1   = A2.copy()
D2   = A2.copy()
D3   = A2.copy()
D1[np.ix_(postPeriodIndices, postPeriodIndices)] = E1
D2[np.ix_(postPeriodIndices, postPeriodIndices)] = E2
D3[np.ix_(postPeriodIndices, postPeriodIndices)] = E3
assert np.abs((C @ np.diag(L) @ C.T) - E0).max() < 1e-12
assert np.abs(E1 @ E1 - E0).max() < 1e-12
assert np.abs(E1 @ E2 - np.eye(E0.shape[0])).max() < 1e-12
assert np.abs(E3 @ E0 - np.eye(E0.shape[0])).max() < 1e-12

n  = 8
x  = [1.86114e-08, .03646372, .3908837, 1, 1.82187e-08, .0364637, .35442, .6091163]
A  = sparse.csc_matrix(F[None, :].astype(float))
b  = np.array([threshold_sumweights]).astype(float)
G1 = sparse.csc_matrix(A1.astype(float))
h1 = np.zeros(n)
G2 = sparse.csc_matrix(np.vstack([np.zeros((1, 8)), -D1]))
h2 = np.concatenate([[np.sqrt(0.25 * (A3 @ D3 @ A3) - (A4 - hh**2))], (0.5 * D2 @ A3)])
G  = sparse.vstack([G1, G2])
h  = np.concatenate([h1, h2])
c  = q.astype(float)

solution = ecos.solve(c, G, h, {'l': 8, 'q': [9]}, A, b)
x_r = np.array([
     0.370467255,
     0.692228716,
     1.005513232,
     1.000000018,
     0.370467239,
     0.321761461,
     0.313284517,
    -0.005513214
])
x_p = solution['x']
x_p = solution['x']

c @ x_p
F @ x_p # = 1
A1 @ x_p # <= 0
x_p @ A2 @ x_p + A3 @ x_p + A4 - hh**2 # <= 0

# F @ x_r = 1
# A1 @ x_r <= 0
# x @ A2 @ x + A3 @ x + A4 - hh**2 <= 0
#
# F @ x_p = 1
# A1 @ x_p <= 0
# x @ A2 @ x + A3 @ x + A4 - hh**2 <= 0



































































































# # cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
# soc_constraints = [
#       cp.SOC(c[i].T @ x + d[i], A[i] @ x + b[i])
# ]
#
# obj = cvxpy.Minimize(q @ x)
# constraints = [
#     F @ x == 1,
#     @ X <= ..
# ]
# prob = cvxpy.Problem(obj, constraints)
#
# # Solve with ECOS.
# prob.solve(solver=ECOS)
# print "optimal value with ECOS:", prob.value








import cvxpy as cp
import numpy as np

# Generate a random feasible SOCP.
m = 3
n = 10
p = 5
n_i = 5
np.random.seed(2)
f = np.random.randn(n)
A = []
b = []
c = []
d = []
x0 = np.random.randn(n)
for i in range(m):
    A.append(np.random.randn(n_i, n))
    b.append(np.random.randn(n_i))
    c.append(np.random.randn(n))
    d.append(np.linalg.norm(A[i] @ x0 + b, 2) - c[i].T @ x0)

F = np.random.randn(p, n)
g = F @ x0

# Define and solve the CVXPY problem.
x = cp.Variable(n)
# We use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
soc_constraints = [
      cp.SOC(c[i].T @ x + d[i], A[i] @ x + b[i]) for i in range(m)
]
prob = cp.Problem(cp.Minimize(f.T@x),
                  soc_constraints + [F @ x == g])
prob.solve()

# Print result.
print("The optimal value is", prob.value)
print("A solution x is")
print(x.value)
for i in range(m):
    print("SOC constraint %i dual variable solution" % i)
    print(soc_constraints[i].dual_value)
