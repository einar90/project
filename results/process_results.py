import numpy as np
import matplotlib.pyplot as plt
import json  # Only used for prettyprinting

ecl_file = 'eclipse/11x11-pressure.dat'         # ECL100, metric units
ecll_file = 'eclipse/11x11-pressure-lab.dat'    # ECL100, lab units
mrst_file = 'mrst/11x11-pressure.dat'           # MRST, metric units
pcm_file = 'peaceman/10x10-pressure.dat'        # Peaceman, dimensionless

p = {
    'ecl': np.loadtxt(ecl_file)*0.986923267,  # bar -> atm
    'ecll': np.loadtxt(ecll_file),            # atm
    'mrst': np.loadtxt(mrst_file),            # atm
    'pcm': np.loadtxt(pcm_file)               # Dimensionless
}

props = {
    'ecl': {
        'k':  300.0/1000.0,       # mD -> D
        'h':  30.0*100.0,         # m -> cm
        'q':  150.0*11.5740741,   # m3/day -> cc/sec
        'mu': 0.5,                # cP
        'dx': 30.0*100.0,         # m -> cm
        'M': int(np.sqrt(p['ecl'].size))-1,
        'p_prod': p['ecl'][0, 0],
        'p_inj':  p['ecl'][-1, -1],
    },
    'ecll': {
        'k':  300.0/1000.0,       # mD -> D
        'h':  30.0,               # cm
        'q':  50000.0/60.0/60.0,  # cc/hr -> cc/sec
        'mu': 0.5,                # cP
        'dx': 30.0,               # cm
        'M': int(np.sqrt(p['ecll'].size))-1,
        'p_prod': p['ecll'][0, 0],
        'p_inj':  p['ecll'][-1, -1],
    },
    'mrst': {
        'k':  .3,                 # D
        'h':  30.0*100.0,         # m -> cm
        'q':  150.0*11.5740741,   # m3/day -> cc/sec
        'mu': 0.5,                # cP
        'dx': 30.0*100,           # m -> cm
        'M': int(np.sqrt(p['mrst'].size))-1,
        'p_prod': p['mrst'][0, 0],
        'p_inj':  p['mrst'][-1, -1],
    },
    'pcm': {                      # Dimensionless
        'k':  1.0,
        'h':  1.0,
        'q':  1.0,
        'mu': 1.0,
        'dx': 1.0,
        'M': int(np.sqrt(p['pcm'].size))-1,
        'p_prod': p['pcm'][0, 0],
        'p_inj':  p['pcm'][-1, -1],
    }
}


def r_eq_exact(props):
    return np.sqrt(2.0) * float(props['M']) * np.exp(
        - np.pi * props['k'] * props['h']
        / (props['q'] * props['mu'])
        * (props['p_inj'] - props['p_prod'])
        - 0.6190
    )


def r_eq_regression(props, p):
    # calculate dimensionless pressure and pressure difference
    p_D = p * props['k'] * props['h'] \
        / (props['q'] * props['mu'])
    p_diff = (p_D - p_D[0, 0])[0:props['M']/2, 0:props['M']/2]

    # compute radius matrix
    r = np.zeros([props['M']/2, props['M']/2])
    for i in range(props['M']/2):
        for j in range(props['M']/2):
            r[i, j] = np.sqrt(i**2 + j**2)

    # linearize matrices and remove well-block
    r = r.reshape([(props['M']/2)**2, ])[1:]
    p_diff = p_diff.reshape([(props['M']/2)**2, ])[1:]

    # create regression line and -polynomial
    reg = np.polyfit(np.log(r), p_diff, deg=1)
    reg_poly = np.poly1d(reg)
    intersection = np.exp(np.roots(reg_poly)[0])
    return (intersection, reg_poly, r, p_diff)


def make_plots(props, p, r_eq, title):
    # get regression results
    (intersection, reg_poly, r_reg, p_diff) = r_eq_regression(props, p)
    r_reg_extended = range(1, 60, 1)
    r_reg_extended = [r/10.0 for r in r_reg_extended]

    # compute full radius matrix
    r = np.zeros([props['M'], props['M']])
    for i in range(props['M']):
        for j in range(props['M']):
            r[i, j] = np.sqrt(i**2 + j**2)

    # plot regression
    titlestring = title + ', M = ' + str(props['M']) \
                        + '; exact: ' + str(r_eq['exact']) \
                        + '; regression: ' + str(r_eq['regression'])
    plt.figure(titlestring)
    ax1 = plt.subplot(1, 2, 1)
    ax1.scatter(r_reg, p_diff)
    ax1.semilogx(r_reg_extended, reg_poly(np.log(r_reg_extended)))
    ax1.set_title('Pressure drop')
    ax1.set_xlabel('$r$')
    ax1.set_ylabel('$\Delta p$')
    ax1.set_xscale('log')
    ax1.set_xlim(1e-1, .6e1)
    ax1.set_ylim(0, .6)

    # plot contour
    ax2 = plt.subplot(1, 2, 2)
    ax2.set_title('Pressure distribution')
    cntplot = ax2.contour(p, 20, colors='k')
    plt.clabel(cntplot, fontsize=9, inline=1)

r_eq = {
    'ecl': {
        'exact': r_eq_exact(props['ecl']),
        'regression': r_eq_regression(
            props['ecl'], p['ecl']
        )[0]
    },
    'ecll': {
        'exact': r_eq_exact(props['ecll']),
        'regression': r_eq_regression(
            props['ecll'], p['ecll']
        )[0]
    },
    'mrst': {
        'exact': r_eq_exact(props['mrst']),
        'regression': r_eq_regression(
            props['mrst'], p['mrst']
        )[0]
    },
    'pcm': {
        'exact': r_eq_exact(props['pcm']),
        'regression': r_eq_regression(
            props['pcm'], p['pcm']
        )[0]
    }
}


print json.dumps(props, sort_keys=True, indent=2)
print 'Calculated r_eq: '
print json.dumps(r_eq, sort_keys=True, indent=2)

make_plots(props['ecl'], p['ecl'], r_eq['ecl'], 'ECL100 Metric')
make_plots(props['ecll'], p['ecll'], r_eq['ecll'], 'ECL100 Lab')
make_plots(props['mrst'], p['mrst'], r_eq['mrst'], 'MRST')
make_plots(props['pcm'], p['pcm'], r_eq['pcm'], 'Peaceman')

plt.draw()
plt.show()
