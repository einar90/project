import numpy as np
import matplotlib.pyplot as plt
import json  # Only used for prettyprinting

ecl_file = 'eclipse/20x20-pressure-corner.dat'

p = {
    'ecl': np.loadtxt(ecl_file)*0.986923267,  # bar -> atm
}

# Set corner points to average of the four corner blocks
p['ecl'][1, 1] = np.average([p['ecl'][0, 0],
                             p['ecl'][1, 0],
                             p['ecl'][0, 1],
                             p['ecl'][1, 1]])

p['ecl'][-2, -2] = np.average([p['ecl'][-1, -1],
                               p['ecl'][-2, -1],
                               p['ecl'][-1, -2],
                               p['ecl'][-2, -2]])

# Drop outermost rows
p['ecl'] = np.delete(p['ecl'], 0, 0)
p['ecl'] = np.delete(p['ecl'], 0, 1)
p['ecl'] = np.delete(p['ecl'], -1, 0)
p['ecl'] = np.delete(p['ecl'], -1, 1)

props = {
    'ecl': {
        'k':  300.0/1000.0,       # mD -> D
        'h':  30.0*100.0,         # m -> cm
        'q':  150.0*11.5740741,   # m3/day -> cc/sec
        'mu': 0.5,                # cP
        'dx': 30.0*100.0,         # m -> cm
        'M': int(np.sqrt(p.get('ecl').size))-1,
        'p_prod': p.get('ecl')[0, 0],
        'p_inj':  p.get('ecl')[-1, -1],
    }
}


def r_eq_exact(props):
    return np.sqrt(2.0) * 10 * np.exp(
        - np.pi * props.get('k') * props.get('h')
        / (props.get('q') * props.get('mu'))
        * (props.get('p_inj') - props.get('p_prod'))
        - 0.6190
    )


def r_eq_regression(props, p):
    # calculate dimensionless pressure and pressure difference
    p_D = p * (props.get('k') * props.get('h')
               / (props.get('q') * props.get('mu')))
    p_diff = (p_D - p_D[0, 0])
    print p_diff
    p_diff = p_diff[1:10, 1:10]

    print p_diff

    # compute radius matrix
    r = np.zeros([9, 9])
    for i in range(1, 10):
        for j in range(1, 10):
            r[i-1, j-1] = np.sqrt(i**2 + j**2)
    print r

    # linearize matrices and remove well-block
    r = r.reshape([9**2, ])
    p_diff = p_diff.reshape([9**2, ])

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
    r = np.zeros([props.get('M'), props.get('M')])
    for i in range(props.get('M')):
        for j in range(props.get('M')):
            r[i, j] = np.sqrt(i**2 + j**2)

    # plot regression
    titlestring = title + ', M = ' + str(props.get('M')) \
                        + '; exact: ' + str(r_eq.get('exact')) \
                        + '; regression: ' + str(r_eq.get('regression'))
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
        'exact': r_eq_exact(props.get('ecl')),
        'regression': r_eq_regression(
            props.get('ecl'), p.get('ecl')
        )[0]
    }
}


print 'M: ', props.get('M')
print json.dumps(props, sort_keys=True, indent=2)
print 'Calculated r_eq: '
print json.dumps(r_eq, sort_keys=True, indent=2)

make_plots(props.get('ecl'), p.get('ecl'), r_eq.get('ecl'), 'ECL100 Metric')



plt.draw()
plt.show()
