import numpy as np
import matplotlib.pyplot as plt

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

# Reshaping pressure matrices to fit on the modelled MxM grid =================
P_delta = np.reshape(P_delta, [M+1, M+1])
P = np.reshape(P, [M+1, M+1])

# Creating radius matrix from producer in upper-left corner ===================
r = np.zeros([M+1, M+1], dtype=float)
for i in range(M+1):
    for j in range(M+1):
        r[i, j] = np.sqrt(float(i)**2.0 + float(j)**2.0)

r_vector = np.reshape(r, [(M+1)**2, ])  # To be used for regression line

# Creating regression line ====================================================
r_reg = np.reshape(r, [(M+1)**2])[1:(M+1)**2-1]
r_reg = np.log10(r_reg)

P_reg = np.reshape(P_delta, [(M+1)**2])[1:(M+1)**2-1]

reg = np.polyfit(r_reg, P_reg, 1)  # Create regression line
reg_poly = np.poly1d(reg)  # Polynomial for regression
print (reg_poly(0), reg_poly(5))


# Creating plots ==============================================================
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.imshow(P)
ax1.set_title('Pressure')
# fig.colorbar(ax1)

ax2.scatter(r, P_delta)
# ax2.plot(r_reg, reg_poly(r_reg))
ax2.set_title('Pressure drop')
ax2.set_xscale('log')
ax2.set_xlim(1e-1, .6e1)
ax2.set_ylim(0, 0.6)

plt.draw()
plt.show()
