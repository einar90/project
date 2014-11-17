import numpy as np
import matplotlib.pyplot as plt
import json  # Only used for prettyprinting

ecl_file = 'eclipse/11x11-pressure.dat'         # ECL100, metric units
ecll_file = 'eclipse/11x11-pressure-lab.dat'    # ECL100, lab units
mrst_file = 'mrst/11x11-pressure.dat'           # MRST, metric units
pcm_file = 'peaceman/10x10-pressure.dat'        # Peaceman solver, dimensionless

p = {
    'ecl': np.loadtxt(ecl_file)*0.986923267,  # bar -> atm
    'ecll': np.loadtxt(ecll_file),            # atm
    'mrst': np.loadtxt(mrst_file),            # atm
    'pcm': np.loadtxt(pcm_file)               # Dimensionless
}

M = {
    'ecl': int(np.sqrt(p.get('ecl').size))-1,
    'ecll': int(np.sqrt(p.get('ecll').size))-1,
    'mrst': int(np.sqrt(p.get('mrst').size))-1,
    'pcm': int(np.sqrt(p.get('pcm').size))-1
}

props = {
    'ecl': {
        'k':  300.0/1000.0,       # mD -> D
        'h':  30.0*100.0,         # m -> cm
        'q':  150.0*11.5740741,   # m3/day -> cc/sec
        'mu': 0.5,                # cP
        'dx': 30.0*100.0,         # m -> cm
        'p_prod': p.get('ecl')[0, 0],
        'p_inj':  p.get('ecl')[M.get('ecl'), M.get('ecl')],
    },
    'ecll': {
        'k':  300.0/1000.0,       # mD -> D
        'h':  30.0,               # cm
        'q':  50000.0/60.0/60.0,  # cc/hr -> cc/sec
        'mu': 0.5,                # cP
        'dx': 30.0,               # cm
        'p_prod': p.get('ecll')[0, 0],
        'p_inj':  p.get('ecll')[M.get('ecll'), M.get('ecll')],
    },
    'mrst': {
        'k':  .3,                 # D
        'h':  30.0*100.0,         # m -> cm
        'q':  150.0*11.5740741,   # m3/day -> cc/sec
        'mu': 0.5,                # cP
        'dx': 30.0*100,           # m -> cm
        'p_prod': p.get('mrst')[0, 0],
        'p_inj':  p.get('mrst')[M.get('mrst'), M.get('mrst')],
    },
    'pcm': {                      # Dimensionless
        'k':  1.0,
        'h':  1.0,
        'q':  1.0,
        'mu': 1.0,
        'dx': 1.0,
        'p_prod': p.get('pcm')[0, 0],
        'p_inj':  p.get('pcm')[M.get('pcm'), M.get('pcm')],
    }
}


def r_eq_exact(M, props):
    return np.sqrt(2.0) * float(M) * np.exp(
        - np.pi * props.get('k') * props.get('h')
        / (props.get('q') * props.get('mu'))
        * (props.get('p_inj') - props.get('p_prod'))
        - 0.6190
    )


def r_eq_regression(M, props, p):
    # calculate dimensionless pressure and pressure difference
    p_D = p * props.get('k') * props.get('h') \
        / (props.get('q') * props.get('mu'))
    p_diff = (p_D - p_D[0, 0])[0:M/2, 0:M/2]

    # compute radius matrix
    r = np.zeros([M/2, M/2])
    for i in range(M/2):
        for j in range(M/2):
            r[i, j] = np.sqrt(i**2 + j**2)

    # linearize matrices and remove well-block
    r = r.reshape([(M/2)**2, ])[1:]
    p_diff = p_diff.reshape([(M/2)**2, ])[1:]

    # create regression line and -polynomial
    reg = np.polyfit(np.log(r), p_diff, deg=1)
    reg_poly = np.poly1d(reg)
    intersection = np.exp(np.roots(reg_poly)[0])
    return (intersection, reg_poly, r, p_diff)


def make_plots(M, props, p, r_eq, title):
    # get regression results
    (intersection, reg_poly, r_reg, p_diff) = r_eq_regression(M, props, p)
    r_reg_extended = range(1, 60, 1)
    r_reg_extended = [r/10.0 for r in r_reg_extended]

    # compute full radius matrix
    r = np.zeros([M, M])
    for i in range(M):
        for j in range(M):
            r[i, j] = np.sqrt(i**2 + j**2)

    # plot regression
    print r_eq
    titlestring = title + ', M = ' + str(M) \
                        + '; exact: ' + str(r_eq.get('exact')) \
                        + '; regression: ' + str(r_eq.get('regression'))
    fig = plt.figure(titlestring)
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
        'exact': r_eq_exact(M.get('ecl'), props.get('ecl')),
        'regression': r_eq_regression(
            M.get('ecl'), props.get('ecl'), p.get('ecl')
        )[0]
    },
    'ecll': {
        'exact': r_eq_exact(M.get('ecll'), props.get('ecll')),
        'regression': r_eq_regression(
            M.get('ecll'), props.get('ecll'), p.get('ecll')
        )[0]
    },
    'mrst': {
        'exact': r_eq_exact(M.get('mrst'), props.get('mrst')),
        'regression': r_eq_regression(
            M.get('mrst'), props.get('mrst'), p.get('mrst')
        )[0]
    },
    'pcm': {
        'exact': r_eq_exact(M.get('pcm'), props.get('pcm')),
        'regression': r_eq_regression(
            M.get('pcm'), props.get('pcm'), p.get('pcm')
        )[0]
    }
}


print 'M: ', M
print json.dumps(props, sort_keys=True, indent=2)
print 'Calculated r_eq: '
print json.dumps(r_eq, sort_keys=True, indent=2)

make_plots(M.get('ecl'), props.get('ecl'), p.get('ecl'), r_eq.get('ecl'),
           'ECL100 Metric')
make_plots(M.get('ecll'), props.get('ecll'), p.get('ecll'), r_eq.get('ecll'),
           'ECL100 Lab')
make_plots(M.get('mrst'), props.get('mrst'), p.get('mrst'), r_eq.get('mrst'),
           'MRST')
make_plots(M.get('pcm'), props.get('pcm'), p.get('pcm'), r_eq.get('pcm'),
           'Peaceman')

plt.draw()
plt.show()
