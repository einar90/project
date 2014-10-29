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

print A


# Populate delta vector =======================================================
delta[0] = 1
delta[-1] = -1
print delta


# Solve system ================================================================
P = np.linalg.solve(A, delta)
P = np.reshape(P, (M+1, M+1))
print P

# Compute pressure drops ======================================================
P_diff = np.zeros([M+1, M+1], dtype=float)
for i in range(M+1):
    for j in range(M+1):
        P_diff[i, j] = P[i, j] - P[0, 0]

print P_diff

# Compute radius matrix =======================================================
r = np.zeros([M+1, M+1], dtype=float)
for i in range(M+1):
    for j in range(M+1):
        r[i, j] = np.sqrt(i**2 + j**2)

print r

# Plotting ====================================================================
r_plot = []
P_plot = []
for i in range(1, 3):
    r_plot.append(r[i, i])
    P_plot.append(P_diff[i, i])

print 'r to be plotted:'
print r_plot

print 'p to be plotted:'
print P_plot

fig = pl.figure()
ax = fig.add_subplot(2, 1, 1)
ax.scatter(r_plot, P_plot)
ax.set_xscale('log')

pl.show()
