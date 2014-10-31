import numpy as np
import matplotlib.pyplot as plt

N = 10

# Parameters
k = 300.0  # mD
h = 30.0   # m
q = 150.0  # m^3/day
mu = .5    # cP

total_conversion_factor = 8.5298e-3
dimensionless_pressure_factor = k*h/q/mu*total_conversion_factor

# Loading data from files
P_ecl = np.loadtxt('eclipse/10x10-pressure.dat')*dimensionless_pressure_factor

# Calculating pressure difference matrix
P_ecl_delta = np.zeros([N, N])
for i in range(N):
    for j in range(N):
        P_ecl_delta[i, j] = P_ecl[i, j] - P_ecl[0, 0]

# Cutting pressure difference matrix
P_ecl_delta = P_ecl_delta[0:N/2, 0:N/2]

# Linearizing pressure difference matrix
P_ecl_delta = np.reshape(P_ecl_delta, [(N/2)**2, ])

# Removing first element (represents well block)
P_ecl_delta = P_ecl_delta[1:]

# Calculating radius matrix
r = np.zeros([N/2, N/2])
for i in range(N/2):
    for j in range(N/2):
        r[i, j] = np.sqrt(float(i**2 + j**2))


# Linearizing radius matrix
r = np.reshape(r, [(N/2)**2, ])

# Removing first element in radius matrix (represents well-block)
r = r[1:]

# Creating regression line
reg_ecl = np.polyfit(np.log(r), P_ecl_delta, deg=1)
reg_ecl_poly = np.poly1d(reg_ecl)

print 'Delta P = 0 at r = ', np.exp(np.roots(reg_ecl_poly)[0])

# r vector for extended plotting
reg_r = [i/10.0 for i in range(1, 61)]

# Plotting pressure difference

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.scatter(r, P_ecl_delta)
ax1.semilogx(reg_r, reg_ecl_poly(np.log(reg_r)))
ax1.set_xscale('log')
ax1.set_title('Pressure drop (ECL)')
# ax1.set_xlim(1e-1, .6e1)
ax1.set_ylim(0, 1.5)

# ax2.scatter(r, P_delta)
# ax2.semilogx(reg_r, reg_poly(np.log(reg_r)))
# ax2.set_title('Pressure drop (MRST)')
# ax2.set_xscale('log')
# ax2.set_xlim(1e-1, .6e1)
# ax2.set_ylim(0, .6)

plt.draw()
plt.show()
