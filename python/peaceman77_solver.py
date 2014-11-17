import numpy as np
import matplotlib.pyplot as plt

M = 10
A = np.zeros([(M+1)**2, (M+1)**2], dtype=float)
delta = np.zeros((M+1)**2, dtype=float)

print '============================================================'
print 'Solving for M = ', M
print '============================================================'

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
        col_i = j / (M+1)
        col_j = j % (M+1)
        if i == j:
            A[i][j] = -4
        else:
            A[i][j] = eq[i].count([col_i, col_j])


# Populate delta vector =======================================================
delta[0] = 1
delta[-1] = -1

# Solve system  and save results ==============================================
P = np.linalg.solve(A, delta)
np.savetxt((str(M) + 'x' + str(M) + '-pressure.dat'), np.reshape(P, [M+1, M+1]))
