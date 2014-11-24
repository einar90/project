import numpy as np
import matplotlib.pyplot as plt

# Filenames
ecl10_file = 'eclipse/10x10-pressure-corner.dat'
ecl20_file = 'eclipse/20x20-pressure-corner.dat'
ecl50_file = 'eclipse/50x50-pressure-corner.dat'
mrst10_file = 'mrst/10x10-pressure-corner.dat'
mrst20_file = 'mrst/20x20-pressure-corner.dat'
mrst50_file = 'mrst/50x50-pressure-corner.dat'

# Load pressures from file
p = {
    'ecl10': np.loadtxt(ecl10_file)*0.986923267,    # bar -> atm
    'ecl20': np.loadtxt(ecl20_file)*0.986923267,    # bar -> atm
    'ecl50': np.loadtxt(ecl50_file)*0.986923267,    # bar -> atm
    'mrst10': np.loadtxt(mrst10_file)*0.986923267,  # bar -> atm
    'mrst20': np.loadtxt(mrst20_file)*0.986923267,  # bar -> atm
    'mrst50': np.loadtxt(mrst50_file)*0.986923267,  # bar -> atm
}

# Set common properties
props = {
    'k':  1000.0/1000.0,       # mD -> D
    'h':  30.0*100.0,         # m -> cm
    'q':  150.0*11.5740741,   # m3/day -> cc/sec
    'mu': 0.5,                # cP
    'dx': 30.0*100.0,         # m -> cm
}


def process_grid_pressures(p):
    ''' Modify the pressure grid to better fit Peaceman's model '''
    # Set corner points to average of the four corner blocks
    p[1, 1] = np.average([p[0, 0], p[1, 0], p[0, 1], p[1, 1]])
    p[-2, -2] = np.average([p[-1, -1], p[-2, -1], p[-1, -2], p[-2, -2]])

    # Drop outermost rows
    p = np.delete(p, 0, 0)
    p = np.delete(p, 0, 1)
    p = np.delete(p, -1, 0)
    p = np.delete(p, -1, 1)

    return p


def r_eq_exact(p):
    ''' Calculate equivalent wellbore radius using the 'exact' method. '''
    M = p.shape[0]-1
    return np.sqrt(2.0) * M * np.exp(
        - np.pi * props.get('k') * props.get('h')
        / (props.get('q') * props.get('mu')) * (p[-1, -1] - p[0, 0])
        - 0.6190
    )


def r_eq_regression(p):
    ''' Calculate pressure using the regression method '''
    # Calculate dimensionless pressure and pressure difference
    p_D = p * (props['k'] * props['h'] / (props['q'] * props['mu']))
    p_diff = (p_D - p_D[0, 0])

    # Cutting matrix
    half = p_diff.shape[0]/2
    p_diff = p_diff[1:half, 1:half]

    # Calculate radius matrix
    r = np.zeros([(half-1), (half-1)])
    for i in range(1, half):
        for j in range(1, half):
            r[i-1, j-1] = np.sqrt(i**2 + j**2)

    # linearize matrices and remove well-block
    r = r.reshape([(half-1)**2, ])
    p_diff = p_diff.reshape([(half-1)**2, ])

    # create regression line and -polynomial
    reg = np.polyfit(np.log(r), p_diff, deg=1)
    reg_poly = np.poly1d(reg)
    intersection = np.exp(np.roots(reg_poly)[0])
    return (intersection, reg_poly, r, p_diff)


def make_plots(p, title):
    # Get regression results
    (intersection, reg_poly, r_reg, p_diff) = r_eq_regression(p)
    r_reg_extended = range(1, 60, 1)
    r_reg_extended = [r/10.0 for r in r_reg_extended]

    # plot regression
    plt.figure(title)
    ax1 = plt.subplot(1, 2, 1)
    ax1.scatter(r_reg, p_diff)
    ax1.semilogx(r_reg_extended, reg_poly(np.log(r_reg_extended)))
    ax1.set_title('Pressure drop')
    ax1.set_xlabel('$r$')
    ax1.set_ylabel('$\Delta p$')
    ax1.set_xscale('log')
    ax1.set_xlim(1e-1, .6e1)
    ax1.set_ylim(0, 2)

    # plot contour
    ax2 = plt.subplot(1, 2, 2)
    ax2.set_title('Pressure distribution')
    cntplot = ax2.contour(p, 20, colors='k')
    plt.clabel(cntplot, fontsize=9, inline=1)

# Calculate and print results; make plots
for key, value in p.iteritems():
    p[key] = process_grid_pressures(value)
    print 'Regression on ', key, ': r_eq = ', r_eq_regression(p[key])[0]
    print 'Exact method on ', key, ': r_eq = ', r_eq_exact(p[key])
    make_plots(p[key], key)

# Show plots
plt.draw()
plt.show()
