import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt('peaceman/32x32-pressure.dat')  # Load pressure from file
M = int(np.sqrt(p.size))-1  # Calculate M


def r_eq_exact(p, M):
    ''' Calculate equivalent wellbore radius using the 'exact' method. '''
    return np.sqrt(2.0) * float(M) * np.exp(
        - np.pi * (p[-1, -1] - p[0, 0]) - 0.6190
        )


def r_eq_regression(p, M):
    ''' Calculate pressure using the regression method '''
    p_diff = (p - p[0, 0])[0:M/2, 0:M/2]  # Calculate pressure difference

    # Calculate radius matrix
    r = np.zeros([M/2, M/2])
    for i in range(M/2):
        for j in range(M/2):
            r[i, j] = np.sqrt(i**2 + j**2)

    # linearize matrices and remove wellblock
    r = r.reshape([(M/2)**2, ])[1:]
    p_diff = p_diff.reshape([(M/2)**2, ])[1:]

    # create regression line and -polynomial
    reg = np.polyfit(np.log(r), p_diff, deg=1)
    reg_poly = np.poly1d(reg)
    intersection = np.exp(np.roots(reg_poly)[0])
    return (intersection, reg_poly, r, p_diff)


def make_plots(p, M):
    # Get regression results
    (intersection, reg_poly, r_reg, p_diff) = r_eq_regression(p, M)

    # Create extended radius vector for plotting the regression line
    r_reg_extended = range(1, 60, 1)
    r_reg_extended = [r/10.0 for r in r_reg_extended]

    # Calculate radius matrix
    r = np.zeros([M, M])
    for i in range(M):
        for j in range(M):
            r[i, j] = np.sqrt(i**2 + j**2)

    # Scatter plot pressure on semilogx axis
    ax1 = plt.subplot(1, 2, 1)
    ax1.scatter(r_reg, p_diff)
    ax1.set_xscale('log')

    #  Plot regression line
    ax1.semilogx(r_reg_extended, reg_poly(np.log(r_reg_extended)))

    # Set axis labels and title
    ax1.set_title('Pressure drop')
    ax1.set_xlabel('$r/\Delta x$')
    ax1.set_ylabel('$\Delta p$')

    # Set axis constraints
    ax1.set_xlim(1e-1, .6e1)
    ax1.set_ylim(0, .6)

    # Create contour plot of the pressure
    ax2 = plt.subplot(1, 2, 2)
    cntplot = ax2.contour(p, 20, colors='k')
    plt.clabel(cntplot, fontsize=9, inline=1)
    ax2.set_title('Pressure distribution')

    # Display plots
    plt.draw()
    plt.show()

# Print solution values
print 'Solution for M = ', str(M)
print 'r_eq from exact method: ', r_eq_exact(p, M)
print 'r_eq from regression method: ', r_eq_regression(p, M)[0]

make_plots(p, M)
