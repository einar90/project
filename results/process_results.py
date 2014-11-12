import numpy as np
import matplotlib.pyplot as plt

props = {
    'ecl': {
        'k':  300.0/1000.0,      # D
        'h':  30.0*100.0,        # cm
        'q':  150.0*11.5740741,  # cc/sec
        'mu': 0.5,               # cP
        'dx': 30.0*100.0,        # cm
    }
}


def r_eq_exact(M, p_prod, p_inj, props):
    print props
    print 'p_00: ', p_prod
    print 'p_MM: ', p_inj
    return np.sqrt(2.0) * float(M) * np.exp(-np.pi * props.get('k') * props.get('h')
                                            / (props.get('q') * props.get('mu'))
                                            * (p_inj - p_prod)
                                            - 0.6190
                                            )

p = {
    'ecl': np.loadtxt('eclipse/11x11-pressure.dat')*0.986923267  # atm
}

M = {
    'ecl': int(np.sqrt(p.get('ecl').size))-1
}


r_eq = {
    'ecl': {
        'exact': r_eq_exact(
            M.get('ecl'),
            p.get('ecl')[0, 0],
            p.get('ecl')[M.get('ecl'), M.get('ecl')],
            props.get('ecl')
        )
    }
}

print M
#print p
print r_eq
