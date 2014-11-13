import numpy as np
import matplotlib.pyplot as plt
import json  # Only used for prettyprinting

pcm_file = 'peaceman/32x32-pressure.dat'   # Peaceman solver results
ecl_file = 'eclipse/11x11-pressure.dat'    # ECL100 simulation results

p = {
    'ecl': np.loadtxt(ecl_file)*0.986923267,  # bar -> atm
    'pcm': np.loadtxt(pcm_file)               # Dimensionless
}

M = {
    'ecl': int(np.sqrt(p.get('ecl').size))-1,
    'pcm': int(np.sqrt(p.get('pcm').size))-1
}

props = {
    'ecl': {
        'k':  300.0/1000.0,      # mD -> D
        'h':  30.0*100.0,        # m -> cm
        'q':  150.0*11.5740741,  # stb/day -> cc/sec
        'mu': 0.5,               # cP
        'dx': 30.0*100.0,        # m -> cm
        'p_prod': p.get('ecl')[0, 0],
        'p_inj':  p.get('ecl')[M.get('ecl'), M.get('ecl')],
    },
    'pcm': {                     # Dimensionless
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


r_eq = {
    'ecl': {
        'exact': r_eq_exact(M.get('ecl'), props.get('ecl'))
    },
    'pcm': {
        'exact': r_eq_exact(M.get('pcm'), props.get('pcm'))
    }
}

print 'M: ', M
print json.dumps(props, sort_keys=True, indent=2)
print 'Calculated r_eq: ', r_eq
