import numpy as np
import matplotlib.pyplot as plt

M = 20
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

# Solve system ================================================================
P = np.linalg.solve(A, delta)

# Calculating pressure difference from producer ===============================
P_delta = []
for p in P:
    P_delta.append(p - P[0])

np.savetxt('pressure_difference.dat', P_delta)

# Reshaping pressure matrices to fit on the modelled MxM grid =================
P_delta = np.reshape(P_delta, [M+1, M+1])
P = np.reshape(P, [M+1, M+1])

# Creating radius matrix from producer in upper-left corner ===================
r = np.zeros([M+1, M+1], dtype=float)
for i in range(M+1):
    for j in range(M+1):
        r[i, j] = np.sqrt(float(i)**2.0 + float(j)**2.0)

# Exact calculation of r_0 ====================================================
print 'Pressure drop: ', P[M, M] - P[0, 0]
print 'Exact solution: ', np.sqrt(2) * M * np.exp(-0.6190 - np.pi * (P[M, M]-P[0, 0]))

# Creating plots ==============================================================

# Cutting matrices
r = r[0:(M+1)/2, 0:(M+1)/2]
P_delta = P_delta[0:(M+1)/2, 0:(M+1)/2]

# Relinearizing matrices
r = np.reshape(r, [((M+1)/2)**2, ])
P_delta = np.reshape(P_delta, [((M+1)/2)**2, ])

# Removing first element (represents the well-block)
r = r[1:]
P_delta = P_delta[1:]

# Create the regression line and polynomial
reg = np.polyfit(np.log(r), P_delta, deg=1)
reg_poly = np.poly1d(reg)
print 'Regression polynomial: ', reg_poly

# r vector for extended plotting
reg_r = [i/10.0 for i in range(1, 61)]
print 'Graphical solution: r_0 = ', np.exp(np.roots(reg_poly)[0])

# Creating plots ==============================================================
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.imshow(P)
ax1.set_title('Pressure')

ax2.scatter(r, P_delta)
ax2.semilogx(reg_r, reg_poly(np.log(reg_r)))
ax2.set_title('Pressure drop')
ax2.set_xscale('log')
ax2.set_xlim(1e-1, .6e1)
ax2.set_ylim(0, .6)

plt.draw()
plt.show()

# Writing data to files
save = raw_input('Save arrays to file? (y/n)  ')
if save == 'y':
    np.savetxt('scatter.dat', (r, P_delta), delimiter=',')
    np.savetxt('regression.dat', (reg_r, reg_poly(np.log(reg_r))), delimiter=',')
