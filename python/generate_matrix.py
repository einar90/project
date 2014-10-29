import numpy as np
import matplotlib.pyplot as pl

M = 10
A = np.zeros([(M+1)**2, (M+1)**2], dtype=float)
delta = np.zeros((M+1)**2, dtype=float)

# Create schema ===============================================================
eq = []

for i in range(M+1):
    for j in range(M+1):
        eq.append([
            [i-1, j], [i+1, j], [i, j-1], [i, j+1], [i, j]
        ])

# Apply boundaries
for line in eq:
    for entry in range(5):
        if line[entry][0] == -1:
            line[entry][0] = 1
        elif line[entry][0] == M+1:
            line[entry][0] = M-1

        if line[entry][1] == -1:
            line[entry][1] = 1
        elif line[entry][1] == M+1:
            line[entry][1] = M-1

# Populate A matrix ===========================================================
for i in range((M+1)**2):
    for j in range((M+1)**2):
        col_i = j / 3
        col_j = j % 3
        if i == j:
            A[i][j] = -4
        else:
            A[i][j] = eq[i].count([col_i, col_j])

# Populate delta vector =======================================================
delta[0] = 1
delta[-1] = -1

# Solve system ================================================================
P = np.linalg.solve(A, delta)

# Calculating pressure difference from producer
P_delta = []
for p in P:
    P_delta.append(p - P[0])

# Creating contour plot
P_delta = np.reshape(P_delta, [M+1, M+1])
P = np.reshape(P, [M+1, M+1])
print P_delta

# fig = pl.figure()
# ax = pl.gca()
# ax.scatter(r, P_delta)
# ax.set_xscale('log')
# pl.xlim([0,6])
# pl.ylim([0,.6])
pl.imshow(P)
pl.colorbar(orientation='vertical')
pl.show()