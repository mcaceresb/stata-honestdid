from scipy.optimize import minimize
import numpy as np


sigma = np.array([
    [ 8.428358e-04,  4.768687e-04, 2.618051e-04, 0.0002354220, 0.0001676371, 0.0001128708, 1.992816e-05, -1.368265e-04],
    [ 4.768687e-04,  6.425420e-04, 3.987425e-04, 0.0002435515, 0.0002201960, 0.0001804591, 3.843765e-05, -2.960422e-05],
    [ 2.618051e-04,  3.987425e-04, 5.229950e-04, 0.0002117686, 0.0001840722, 0.0001458528, 7.005197e-05,  5.952995e-05],
    [ 2.354220e-04,  2.435515e-04, 2.117686e-04, 0.0003089595, 0.0001197866, 0.0001334081, 1.016335e-04,  1.079052e-04],
    [ 1.676371e-04,  2.201960e-04, 1.840722e-04, 0.0001197866, 0.0003599704, 0.0002478819, 1.749579e-04,  1.654257e-04],
    [ 1.128708e-04,  1.804591e-04, 1.458528e-04, 0.0001334081, 0.0002478819, 0.0004263950, 2.171438e-04,  2.892748e-04],
    [ 1.992816e-05,  3.843765e-05, 7.005197e-05, 0.0001016335, 0.0001749579, 0.0002171438, 4.886698e-04,  3.805322e-04],
    [-1.368265e-04, -2.960422e-05, 5.952995e-05, 0.0001079052, 0.0001654257, 0.0002892748, 3.805322e-04,  7.617394e-04]
])

A_quadratic_sd = np.array([
    [0, 0, 0, 0,             0,             0,             0,             0],
    [0, 0, 0, 0,             0,             0,             0,             0],
    [0, 0, 0, 0,             0,             0,             0,             0],
    [0, 0, 0, 0,             0,             0,             0,             0],
    [0, 0, 0, 0,   .0005316404,  -.0000287359,  -.0001288079,  -8.12950e-06],
    [0, 0, 0, 0,  -.0000287359,    .000368052,  -.0001560354,   .0000317829],
    [0, 0, 0, 0,  -.0001288079,  -.0001560354,   .0004084173,  -.0000971909],
    [0, 0, 0, 0,  -8.12950e-06,   .0000317829,  -.0000971909,   .0003089595]
])

A_lienar_sd = np.array([
     0.000000e+00,
     0.000000e+00,
     0.000000e+00,
     0.000000e+00,
    -1.051178e-04,
     7.224761e-05,
     1.285712e-04,
     2.395732e-04
])

A_constant_sd = 0.0003599704

K = 4
lowerTriMat = np.tri(4)
A_absolutevalue = np.row_stack((
    np.column_stack((-np.eye(K),  lowerTriMat)),
    np.column_stack((-np.eye(K), -lowerTriMat))
))
A_sumweights = np.concatenate((np.zeros(K), np.ones(K)))


def func(x, *args):
    return x.T.dot(A_quadratic_sd.dot(x)) + A_lienar_sd.T.dot(x) + A_constant_sd


def func_deriv(x, *args):
    return A_quadratic_sd.dot(x) + A_lienar_sd


cons = ({'type': 'eq',
         'fun' : lambda x: A_sumweights.dot(x) - 1,
         'jac' : lambda x: A_sumweights},
        {'type': 'ineq',
         'fun' : lambda x: -A_absolutevalue.dot(x),
         'jac' : lambda x: -A_absolutevalue})

cons = ({'type': 'eq',
         'fun' : lambda x: A_sumweights.dot(x) - 1},
        {'type': 'ineq',
         'fun' : lambda x: -A_absolutevalue.dot(x)})

x0  = np.ones(x.shape)
res = minimize(func, x0, -A_absolutevalue.dot(x), constraints=cons)
cons[0]['fun'](res.x)
cons[1]['fun'](res.x)
res.fun, func(r)
for xx, rr in zip(res.x, r):
    print(f"{xx:.5f}\t{rr:.5f}")

r = np.array([
     0.767721234,
     1.500062829,
     2.307930431,
     2.555812619,
     0.370475145,
     0.321767143,
     0.313283593,
    -0.005525881
])

func(r)
cons[0]['fun'](r)
cons[1]['fun'](r)

c = np.array([
    0.76772654,
    1.50007612,
    2.30795471,
    2.55584621,
    0.370465624,
    0.321752658,
    0.313218465,
    -0.00566719151
])

func(c)
cons[0]['fun'](c)
cons[1]['fun'](c)










Hello!

Here's a good news, bad news sandwhich:

- I was able to replicate the createSensitivityResults from the vignette from Stata, including the plot.

- Stata's constrained optimization library does not allow for inequality costraints, which are used in `.findLowestH` and `.findWorstCaseBiasGivenH`.

- The underlying libraries used by CVXR are available for C and C++, which I used. Stata has a plugin interface that allows it to interact with C and C++ code without having the user install anything other than the Stata package.

I've used plugins in other Stata projects and they work reasonably well. The only drawback is that there's a possibility they fail to load on some systems. This is not very common and is fixable, but the fix, which is for the user to compile the plugin themselves, can be very easy or very hard depending on the user.

The alternative to plugins would be a wrapper, be it a Python wrapper (Stata allows calling Python from within since version 16) or an R wrapper (e.g. rcall). Neither option seems as appealing since it would require the user to install the requisite packages in another program, but LMK.





import osqp
import numpy as np
import scipy as sp
from scipy import sparse


K  = 4
Ab = 2 * A_quadratic_sd[K:, K:]
Ab[np.tril_indices_from(Ab, -1)] = 0
P = sparse.block_diag([sparse.csc_matrix((K, K)), sparse.triu(Ab)], format='csc')
# P  = sparse.block_diag([sparse.csc_matrix((K, K)), Ab], format='csc')
q  = A_lienar_sd
A  = sparse.csc_matrix(np.vstack([A_sumweights, A_absolutevalue]))
u  = np.array([1] + [0] * 8)
l  = np.array([1] + [-np.inf] * 8)

# Create an OSQP object
prob = osqp.OSQP()

# Setup workspace
prob.setup(P, q, A, l, u)

# Solve problem
res = prob.solve()

func(r), func(res.x), func(c)
cons[0]['fun'](res.x)
cons[1]['fun'](res.x)
