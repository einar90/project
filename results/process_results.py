import numpy as np
import matplotlib.pyplot as plt
import json  # Only used for prettyprinting

pcm_file = 'peaceman/10x10-pressure.dat'        # Peaceman, dimensionless

p = {
    'pcm': np.loadtxt(pcm_file)               # Dimensionless
}

props = {
    'pcm': {                      # Dimensionless
        'M': int(np.sqrt(p['pcm'].size))-1,
        'p_prod': p['pcm'][0, 0],
        'p_inj':  p['pcm'][-1, -1],
    }
}


def r_eq_exact(props):
    return np.sqrt(2.0) * float(props['M']) * np.exp(
        - np.pi * (props['p_inj'] - props['p_prod']) - 0.6190)


def r_eq_regression(props, p):
    # calculate dimensionless pressure and pressure difference
    p_D = p
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

make_plots(props['pcm'], p['pcm'], r_eq['pcm'], 'Peaceman')

plt.draw()
plt.show()
